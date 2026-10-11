
"""GROMACS Adapter.

Ehemals unter ``electrofit.external.gromacs`` – Pfad vollständig entfernt.
"""
from __future__ import annotations

import logging
import os
import re
import shutil
import subprocess
import matplotlib.pyplot as plt
import pandas as pd
from electrofit.cli.run_commands import run_command
from electrofit.io.ff import (
    include_ff,
    include_ions,
    include_tip3p,
    remove_defaults_section_lines,
    replace_posres_in_file,
)
from electrofit.io.files import strip_extension
from electrofit.infra.logging import setup_logging
from electrofit.cli.safe_run import ensure_finalized
from electrofit.infra.scratch_manager import scratch_dir_for, setup_scratch_directory
from electrofit.simulation.ions import DerivedSimulation

__all__ = ["mdrun_flags", "set_up_production", "resume_production"]

_OFFLOAD_CHOICES = ("auto", "gpu", "cpu")

# Files mdrun needs to append to an interrupted production run. The checkpoint is
# last, so a directory holding md.cpt after a move also holds the other files.
_RESUME_FILES = ("md.log", "md.xtc", "md.trr", "md.edr", "md_prev.cpt", "md.tpr", "md.cpt")


def plot_svg(svg: str) -> None:
    logging.disable(logging.CRITICAL)
    name = strip_extension(svg)
    df = pd.read_csv(svg, sep=r"\s+", header=None, names=["time", name])
    plt.figure(figsize=(8, 6))
    plt.plot(df["time"], df[name], color="darkblue", linestyle="-", label=f"{name}")
    plt.xlabel("time (s)", fontsize=14)
    plt.ylabel(f"{name}", fontsize=14)
    plt.axhline(0, color="black", linewidth=0.8, linestyle="--")
    plt.tight_layout()
    plt.savefig(f"{name}.pdf", format="pdf")
    plt.close()
    logging.disable(logging.NOTSET)


def _format_lengths(lengths: tuple[float, float, float]) -> str:
    return " ".join(f"{value:.5f}" for value in lengths)


def _build_editconf_command(
    m_gro: str,
    output: str,
    box_type: str,
    lengths: tuple[float, float, float] | None,
    angles: tuple[float, float, float] | None,
    edge_nm: float | None,
) -> str:
    base = f"echo 0 | gmx editconf -f {m_gro} -o {output} -bt {box_type}"
    if lengths is not None:
        base += f" -box {_format_lengths(lengths)}"
        if angles and any(abs(a - 90.0) > 1e-6 for a in angles):
            base += " -angles " + " ".join(f"{a:.2f}" for a in angles)
    else:
        if edge_nm is None:
            raise ValueError("Box edge distance is required when no explicit lengths are provided.")
        base += f" -d {edge_nm}"
    base += " -c -princ -nobackup"
    return base


def _build_genion_command(
    m_top: str,
    derived: DerivedSimulation,
    cation: str,
    anion: str,
    salt_concentration: float | None,
) -> str:
    if derived.mode == "target":
        pos = derived.ions.positive_count
        neg = derived.ions.negative_count
        if pos is None or neg is None:
            raise ValueError("Target ion mode requires explicit positive and negative counts.")
        return (
            f'echo "SOL" | gmx genion -s ions.tpr -o {m_top}_tip3p_ions.gro '
            f"-p {m_top}.top -pname {cation} -nname {anion} -np {pos} -nn {neg}"
        )
    conc = salt_concentration if salt_concentration is not None else 0.0
    return (
        f'echo "SOL" | gmx genion -s ions.tpr -o {m_top}_tip3p_ions.gro '
        f"-conc {conc} -p {m_top}.top -pname {cation} -nname {anion} -neutral"
    )


def _clamp_threads(threads: int | None) -> int | None:
    try:
        import psutil

        logical_cpus = psutil.cpu_count(logical=True) or (os.cpu_count() or 1)
    except Exception:
        logical_cpus = os.cpu_count() or 1
    if threads is not None and threads > logical_cpus:
        logging.warning(
            "Requested threads=%d exceeds logical CPUs=%d -> clamping",
            threads,
            logical_cpus,
        )
        return logical_cpus
    return threads


def mdrun_flags(
    *,
    threads: int | None,
    pin: bool | None,
    gpu: bool | None,
    gpu_id: str | int | None = None,
    pme: str = "auto",
    bonded: str = "auto",
    update: str = "auto",
    checkpoint_minutes: float | None = None,
    minimization: bool = False,
) -> str:
    """Return the resource flags for one ``gmx mdrun`` call.

    CPU runs use thread-MPI with ``-nt``. GPU runs use one rank (``-ntmpi 1``) with
    ``threads`` OpenMP threads and nonbonded interactions on the GPU. PME, bonded and
    update placement default to GROMACS' automatic choice, which offloads each task
    where the stage supports it: energy minimisation keeps PME, bonded and update on
    the CPU, and GROMACS 2025 runs the update on the GPU for the md integrator but
    not for sd. Explicit 'gpu' or 'cpu' placements apply to dynamics stages only.
    ``checkpoint_minutes`` sets the checkpoint interval of dynamics stages (``-cpt``;
    GROMACS writes one every 15 minutes by default).
    """
    for option, value in (("pme", pme), ("bonded", bonded), ("update", update)):
        if value not in _OFFLOAD_CHOICES:
            raise ValueError(f"gpu_{option} must be one of {_OFFLOAD_CHOICES}, got {value!r}")
    flags: list[str] = []
    if gpu:
        flags += ["-ntmpi", "1"]
        if threads is not None:
            flags += ["-ntomp", str(threads)]
    elif threads is not None:
        flags += ["-nt", str(threads)]
        # For tiny systems, forcing single rank avoids DD issues on some builds
        if threads == 1:
            flags += ["-ntmpi", "1"]
    if pin is not None:
        flags += ["-pin", "on" if pin else "off"]
    if gpu:
        flags += ["-nb", "gpu"]
        if not minimization:
            for option, value in (("-pme", pme), ("-bonded", bonded), ("-update", update)):
                if value != "auto":
                    flags += [option, value]
        if gpu_id is not None and str(gpu_id) != "":
            flags += ["-gpu_id", str(gpu_id)]
    if checkpoint_minutes is not None and not minimization:
        flags += ["-cpt", f"{float(checkpoint_minutes):g}"]
    flags.append("-nobackup")
    return " ".join(flags)


def _postprocess_production(scratch_dir: str, postprocess_on_scratch: bool) -> None:
    if not postprocess_on_scratch:
        logging.info(
            "[postprocess] Skipping md_nojump/md_center/mindist/gyrate on scratch "
            "(gromacs.runtime.postprocess_on_scratch=false)."
        )
        return
    run_command(
        'echo 0 | gmx trjconv -s md.tpr -f md.xtc -o md_nojump.xtc -pbc nojump',
        cwd=scratch_dir,
    )
    run_command(
        'echo "1\n0\n" | gmx trjconv -s md.tpr -f md_nojump.xtc -o md_center.xtc -center -pbc mol',
        cwd=scratch_dir,
    )
    run_command(
        'echo "1\n" | gmx mindist -s md.tpr -f md_center.xtc -pi -od mindist.xvg',
        cwd=scratch_dir,
    )
    run_command(
        'echo "1" | gmx gyrate -f md_center.xtc -s md.tpr -o gyrate.xvg -xvg none',
        cwd=scratch_dir,
    )
    plot_svg("gyrate.xvg")


def _production_complete(state_dir: str) -> bool:
    """True if md.cpt in state_dir is at the nsteps recorded in md.log."""
    try:
        with open(os.path.join(state_dir, "md.log"), errors="ignore") as handle:
            header = handle.read(200000)
        nsteps = int(re.search(r"^\s+nsteps\s+=\s+(-?\d+)", header, re.M).group(1))
        dump = subprocess.run(
            ["gmx", "dump", "-cp", os.path.join(state_dir, "md.cpt")],
            capture_output=True, text=True, check=True,
        )
        step = int(re.search(r"^\s*step\s+=\s+(\d+)", dump.stdout, re.M).group(1))
    except Exception:
        return False
    logging.info("[resume] Checkpoint at step %d of nsteps %d in %s", step, nsteps, state_dir)
    return nsteps >= 0 and step >= nsteps


def _move_verified(src: str, dst: str) -> None:
    """Copy src to dst, check the size, then remove src."""
    shutil.copy2(src, dst)
    if os.path.getsize(src) != os.path.getsize(dst):
        raise OSError(f"Size mismatch after copying {src} -> {dst}; source kept")
    os.remove(src)
    logging.info("[resume] Moved %s -> %s", src, dst)


def set_up_production(
    m_gro,
    MDP_dir,
    base_scratch_dir,
    molecule_name,
    simulation,
    *,
    derived: DerivedSimulation | None = None,
    ff="amber14sb.ff",
    threads: int | None = None,
    pin: bool | None = None,
    gpu: bool | None = None,
    gpu_id: str | int | None = None,
    gpu_pme: str = "auto",
    gpu_bonded: str = "auto",
    gpu_update: str = "auto",
    checkpoint_minutes: float | None = None,
    postprocess_on_scratch: bool = True,
    prepare_only: bool = False,
):
    """Build, equilibrate and run production MD with GROMACS in a scratch directory.

    A force-field folder named ``ff`` inside the run directory is staged into
    scratch with the topology, so grompp uses that copy.
    """
    fullpath = os.getcwd()
    log_file_path = os.path.join(fullpath, "process.log")
    suppress = any(isinstance(h, logging.FileHandler) for h in logging.getLogger().handlers)
    setup_logging(log_file_path, suppress_initial_message=suppress)

    previous_scratch = scratch_dir_for(fullpath, base_scratch_dir)
    if os.path.isfile(os.path.join(previous_scratch, "md.cpt")):
        raise RuntimeError(
            f"{previous_scratch} holds a production checkpoint from an earlier run; continue it "
            "with 'electrofit step8 --resume' or move it away before starting a new run"
        )

    if derived is None:
        derived = getattr(simulation, "derived", None)
    if derived is None:
        raise ValueError("Simulation derived settings unavailable; ensure load_config was used.")

    m_gro_name = strip_extension(m_gro)
    m_itp = f"{m_gro_name}.itp"
    m_top = f"{m_gro_name}.top"
    m_posres = f"posre_{m_gro_name.replace('_GMX', '')}.itp"
    mdp_dirname = os.path.basename(os.path.normpath(MDP_dir))

    for path_check in (m_itp, m_top, m_posres):
        if not os.path.isfile(path_check):
            raise FileNotFoundError(f"Missing required file: {path_check}")

    replace_posres_in_file(m_top)
    include_ff(m_top, ff)
    include_tip3p(m_top, ff)
    include_ions(m_top, ff)
    remove_defaults_section_lines(m_top)

    input_files = [m_gro, m_itp, m_top, m_posres, MDP_dir]
    if os.path.isdir(ff):
        input_files.append(ff)
    scratch_dir, original_dir = setup_scratch_directory(input_files, base_scratch_dir)

    with ensure_finalized(original_dir=original_dir, scratch_dir=scratch_dir, input_files=input_files):
        os.chdir(scratch_dir)

        if derived.mode == "target":
            logging.info(
                "[ions] Target insertion: %s=%d, %s=%d (neutrality=%s)",
                derived.ions.cation,
                derived.ions.positive_count,
                derived.ions.anion,
                derived.ions.negative_count,
                derived.ions.enforce_neutrality,
            )
        else:
            logging.info(
                "[ions] Salt concentration mode: %s/%s at %.3f mol/L",
                simulation.ions.cation,
                simulation.ions.anion,
                simulation.ions.salt_concentration,
            )

        editconf_cmd = _build_editconf_command(
            m_gro,
            f"{m_gro_name}_box.gro",
            simulation.box.type,
            derived.box.lengths_nm,
            derived.box.angles_deg,
            simulation.box.edge_nm,
        )
        run_command(editconf_cmd, cwd=scratch_dir)

        run_command(
            f"gmx solvate -cp {m_gro_name}_box.gro -cs spc216 -o {m_gro_name}_tip3p.gro -p {m_gro_name}.top -nobackup",
            cwd=scratch_dir,
        )
        run_command(f"tail {m_gro_name}.top", cwd=scratch_dir)
        run_command("touch ions.mdp", cwd=scratch_dir)
        run_command(
            f"gmx grompp -f ions.mdp -c {m_gro_name}_tip3p.gro -p {m_gro_name}.top -o ions.tpr",
            cwd=scratch_dir,
        )

        genion_cmd = _build_genion_command(
            m_gro_name,
            derived,
            simulation.ions.cation,
            simulation.ions.anion,
            simulation.ions.salt_concentration,
        )
        run_command(genion_cmd, cwd=scratch_dir)
        run_command(f"tail {m_gro_name}.top", cwd=scratch_dir)

        if prepare_only:
            logging.info(
                "[prepare-only] Solvated and ionized GROMACS inputs are complete; "
                "skipping EM, NVT, NPT, production, and postprocessing."
            )
            return

        threads = _clamp_threads(threads)
        offload = dict(
            gpu_id=gpu_id, pme=gpu_pme, bonded=gpu_bonded, update=gpu_update, checkpoint_minutes=checkpoint_minutes
        )
        em_flags = mdrun_flags(threads=threads, pin=pin, gpu=gpu, minimization=True, **offload)
        md_flags = mdrun_flags(threads=threads, pin=pin, gpu=gpu, **offload)
        logging.info("[mdrun] minimization flags: %s | dynamics flags: %s", em_flags, md_flags)

        run_command(
            f"gmx grompp -f {mdp_dirname}/em_steep.mdp -c {m_gro_name}_tip3p_ions.gro -p {m_gro_name}.top -o em_steep.tpr",
            cwd=scratch_dir,
        )
        run_command(f"gmx mdrun -deffnm em_steep {em_flags}", cwd=scratch_dir)
        run_command(
            'echo "Potential\n0\n" | gmx energy -f em_steep.edr -o potential.xvg -xvg none',
            cwd=scratch_dir,
        )
        plot_svg("potential.xvg")

        run_command(
            f"gmx grompp -f {mdp_dirname}/NVT.mdp -c em_steep.gro -r em_steep.gro -p {m_gro_name}.top -o nvt.tpr -nobackup",
            cwd=scratch_dir,
        )
        run_command(f"gmx mdrun -deffnm nvt {md_flags}", cwd=scratch_dir)
        run_command(
            'echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg -xvg none -b 20',
            cwd=scratch_dir,
        )
        plot_svg("temperature.xvg")
        run_command(
            f"gmx grompp -f {mdp_dirname}/NPT.mdp -c nvt.gro -r nvt.gro -p {m_gro_name}.top -o npt.tpr -nobackup",
            cwd=scratch_dir,
        )
        run_command(f"gmx mdrun -deffnm npt {md_flags}", cwd=scratch_dir)
        run_command(
            'echo "Pressure" | gmx energy -f npt.edr -o pressure.xvg -xvg none',
            cwd=scratch_dir,
        )
        plot_svg("pressure.xvg")

        # Pick latest equilibrated state: prefer NPT if available, else NVT
        prod_conf = "npt.gro" if os.path.isfile(os.path.join(scratch_dir, "npt.gro")) else "nvt.gro"
        prod_cpt = (
            "npt.cpt"
            if prod_conf.startswith("npt") and os.path.isfile(os.path.join(scratch_dir, "npt.cpt"))
            else "nvt.cpt"
        )
        if prod_conf.startswith("npt"):
            logging.info("[prod] Using NPT output (%s, %s) for production start", prod_conf, prod_cpt)
        else:
            logging.warning(
                "[prod] NPT outputs missing -> falling back to NVT (%s, %s)",
                prod_conf,
                prod_cpt,
            )

        run_command(
            f"gmx grompp -f {mdp_dirname}/Production.mdp -c {prod_conf} -t {prod_cpt} -p {m_gro_name}.top -o md.tpr",
            cwd=scratch_dir,
        )
        run_command(f"gmx mdrun -deffnm md {md_flags}", cwd=scratch_dir)
        _postprocess_production(scratch_dir, postprocess_on_scratch)


def resume_production(
    m_gro,
    MDP_dir,
    base_scratch_dir,
    *,
    ff="amber14sb.ff",
    threads: int | None = None,
    pin: bool | None = None,
    gpu: bool | None = None,
    gpu_id: str | int | None = None,
    gpu_pme: str = "auto",
    gpu_bonded: str = "auto",
    gpu_update: str = "auto",
    checkpoint_minutes: float | None = None,
    postprocess_on_scratch: bool = True,
) -> bool:
    """Continue an interrupted production run from its checkpoint (``mdrun -cpi``).

    The state is taken from the run's scratch directory when it holds md.tpr and
    md.cpt (for example after a reboot or a killed session). Otherwise it is taken
    from the run directory, where an interrupted run copies its scratch contents;
    those production files are moved into a new scratch directory so that the run
    keeps a single copy of its trajectory. mdrun appends to the existing output
    files, and the results are copied back as in a normal run. Returns False
    without changes if the checkpoint is already at the target step.
    """
    run_dir = os.getcwd()
    log_file_path = os.path.join(run_dir, "process.log")
    suppress = any(isinstance(h, logging.FileHandler) for h in logging.getLogger().handlers)
    setup_logging(log_file_path, suppress_initial_message=suppress)

    m_gro_name = strip_extension(m_gro)
    input_files = [
        m_gro,
        f"{m_gro_name}.itp",
        f"{m_gro_name}.top",
        f"posre_{m_gro_name.replace('_GMX', '')}.itp",
        MDP_dir,
    ]
    if os.path.isdir(ff):
        input_files.append(ff)

    scratch_dir = scratch_dir_for(run_dir, base_scratch_dir)

    def holds(directory, *names):
        return all(os.path.isfile(os.path.join(directory, name)) for name in names)

    if holds(scratch_dir, "md.cpt") and holds(run_dir, "md.cpt"):
        raise RuntimeError(
            f"Both {scratch_dir} and {run_dir} hold md.cpt; keep the state to continue "
            "and move the other away"
        )
    state_dir = scratch_dir if holds(scratch_dir, "md.cpt") else run_dir
    if holds(state_dir, "md.cpt") and _production_complete(state_dir):
        logging.info("[resume] Production in %s is already complete; nothing to continue", state_dir)
        return False
    if holds(scratch_dir, "md.tpr", "md.cpt"):
        logging.info("[resume] Continuing the checkpoint in scratch %s", scratch_dir)
    elif holds(run_dir, "md.cpt") and (holds(run_dir, "md.tpr") or holds(scratch_dir, "md.tpr")):
        logging.info("[resume] Continuing the checkpoint in run directory %s", run_dir)
        scratch_dir, _ = setup_scratch_directory(input_files, base_scratch_dir)
        for name in _RESUME_FILES:
            if holds(run_dir, name):
                _move_verified(os.path.join(run_dir, name), os.path.join(scratch_dir, name))
    else:
        raise FileNotFoundError(
            f"No production checkpoint (md.tpr and md.cpt) in {scratch_dir} or {run_dir}"
        )

    threads = _clamp_threads(threads)
    flags = mdrun_flags(
        threads=threads, pin=pin, gpu=gpu, gpu_id=gpu_id, pme=gpu_pme, bonded=gpu_bonded,
        update=gpu_update, checkpoint_minutes=checkpoint_minutes,
    )
    logging.info("[mdrun] dynamics flags: %s", flags)
    with ensure_finalized(original_dir=run_dir, scratch_dir=scratch_dir, input_files=input_files):
        os.chdir(scratch_dir)
        run_command(f"gmx mdrun -deffnm md -cpi md.cpt {flags}", cwd=scratch_dir)
        _postprocess_production(scratch_dir, postprocess_on_scratch)
    return True
