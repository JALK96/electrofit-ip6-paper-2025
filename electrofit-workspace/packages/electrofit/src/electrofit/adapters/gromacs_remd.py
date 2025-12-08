"""GROMACS Adapter for REMD preparation (equilibrate + per-replica .tpr)."""
from __future__ import annotations

import logging
import os
from typing import List

from electrofit.cli.run_commands import run_command
from electrofit.adapters.gromacs import _build_editconf_command, _build_genion_command, plot_svg
from electrofit.io.ff import (
    include_ff,
    include_ions,
    include_tip3p,
    remove_defaults_section_lines,
    replace_posres_in_file,
)
from electrofit.io.files import strip_extension
from electrofit.infra.scratch_manager import setup_scratch_directory
from electrofit.cli.safe_run import ensure_finalized


def _mdrun_flags(threads: int | None, pin: bool | None, gpu: bool | None = None) -> str:
    flags: list[str] = []
    # On small boxes, prefer a single rank to avoid DD issues; clamp threads.
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
        threads = logical_cpus

    if threads is not None:
        flags += ["-nt", str(threads)]
        if threads == 1:
            flags += ["-ntmpi", "1"]
    if pin is not None:
        flags += ["-pin", "on" if pin else "off"]
    if gpu:
        flags += ["-nb", "gpu"]
    flags.append("-nobackup")
    return " ".join(flags)


def _read_lines(path: str) -> List[str]:
    with open(path, "r") as f:
        return f.read().splitlines()


def _write_lines(path: str, lines: List[str]) -> None:
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")


def _norm_key(s: str) -> str:
    return s.strip().lower().replace("_", "-")


def _set_key(lines: List[str], key: str, value: str) -> List[str]:
    """Replace or append a single key using its exact name."""
    out: List[str] = []
    found = False
    k = key.strip().lower()
    for line in lines:
        parts = line.split("=", 1)
        if len(parts) == 2 and parts[0].strip().lower() == k:
            out.append(f"{key} = {value}")
            found = True
        else:
            out.append(line)
    if not found:
        out.append(f"{key} = {value}")
    return out


def _set_key_alias(lines: List[str], primary: str, aliases: List[str], value: str) -> List[str]:
    """Replace any of the alias keys; if none found, append the primary key.

    Aliases are compared case-insensitively and with '-' vs '_' normalized.
    """
    out: List[str] = []
    found = False
    keys_norm = {_norm_key(primary)} | {_norm_key(a) for a in aliases}
    for line in lines:
        parts = line.split("=", 1)
        if len(parts) == 2 and _norm_key(parts[0]) in keys_norm:
            out.append(f"{primary} = {value}")
            found = True
        else:
            out.append(line)
    if not found:
        out.append(f"{primary} = {value}")
    return out


def _make_replica_mdp(template: str, out_path: str, T: float, seed: int | None = None) -> None:
    lines = _read_lines(template)
    # Set temperature
    lines = _set_key_alias(lines, "ref_t", ["ref-t"], f"{T:.3f}")
    # Force NVT
    lines = _set_key_alias(lines, "pcoupl", [], "no")
    # Continuation and velocity generation
    lines = _set_key_alias(lines, "continuation", [], "yes")
    lines = _set_key_alias(lines, "gen-vel", ["gen_vel"], "no")
    # Langevin seed (if applicable)
    if seed is not None:
        lines = _set_key_alias(lines, "ld-seed", ["ld_seed"], str(seed))
    _write_lines(out_path, lines)


def _geometric_ladder(nrep: int, tmin: float, tmax: float) -> list[float]:
    if nrep < 2:
        return [tmin]
    ratio = tmax / tmin
    return [tmin * (ratio ** (i / (nrep - 1))) for i in range(nrep)]


def set_up_remd(
    *,
    m_gro: str,
    MDP_dir: str,
    base_scratch_dir: str,
    molecule_name: str,
    simulation,
    remd_cfg,
    ff: str = "amber14sb.ff",
    threads: int | None = None,
    pin: bool | None = None,
    gpu: bool | None = None,
):
    """Equilibrate and generate REMD inputs in a scratch dir; copy back to run dir."""
    m_gro_name = strip_extension(m_gro)
    m_itp = f"{m_gro_name}.itp"
    m_top = f"{m_gro_name}.top"
    m_posres = f"posre_{m_gro_name.replace('_GMX', '')}.itp"
    mdp_dirname = os.path.basename(os.path.normpath(MDP_dir))

    for path_check in (m_itp, m_top, m_posres, MDP_dir):
        if not os.path.exists(path_check):
            raise FileNotFoundError(f"Missing required input: {path_check}")

    # Prepare topology includes
    replace_posres_in_file(m_top)
    include_ff(m_top, ff)
    include_tip3p(m_top, ff)
    include_ions(m_top, ff)
    remove_defaults_section_lines(m_top)

    input_items = [m_gro, m_itp, m_top, m_posres, MDP_dir]
    scratch_dir, original_dir = setup_scratch_directory(input_items, base_scratch_dir)

    with ensure_finalized(original_dir=original_dir, scratch_dir=scratch_dir, input_files=input_items):
        os.chdir(scratch_dir)

        # Box → solvate → ions
        # Prefer derived box lengths from ion targets; fall back to edge padding.
        derived = getattr(simulation, "derived", None)
        lengths = getattr(getattr(derived, "box", None), "lengths_nm", None) if derived else None
        angles = getattr(getattr(derived, "box", None), "angles_deg", None) if derived else None
        editconf_cmd = _build_editconf_command(
            m_gro,
            f"{m_gro_name}_box.gro",
            simulation.box.type,
            lengths,
            angles,
            simulation.box.edge_nm,
        )
        run_command(editconf_cmd, cwd=scratch_dir)
        run_command(
            f"gmx solvate -cp {m_gro_name}_box.gro -cs spc216 -o {m_gro_name}_tip3p.gro -p {m_gro_name}.top -nobackup",
            cwd=scratch_dir,
        )
        run_command(f"tail {m_gro_name}.top", cwd=scratch_dir)
        # ions.mdp only needs to exist
        run_command("touch ions.mdp", cwd=scratch_dir)
        run_command(
            f"gmx grompp -f ions.mdp -c {m_gro_name}_tip3p.gro -p {m_gro_name}.top -o ions.tpr",
            cwd=scratch_dir,
        )

        # Decide ion mode
        derived = getattr(simulation, "derived", None)
        if derived is None:
            raise ValueError("simulation.derived missing; ensure load_config() was used.")
        genion_cmd = _build_genion_command(
            m_gro_name,
            derived,
            simulation.ions.cation,
            simulation.ions.anion,
            simulation.ions.salt_concentration,
        )
        run_command(genion_cmd, cwd=scratch_dir)
        run_command(f"tail {m_gro_name}.top", cwd=scratch_dir)

        # EM → NVT → NPT
        run_command(
            f"gmx grompp -f {mdp_dirname}/em_steep.mdp -c {m_gro_name}_tip3p_ions.gro -p {m_gro_name}.top -o em_steep.tpr",
            cwd=scratch_dir,
        )
        run_command(f"gmx mdrun -deffnm em_steep {_mdrun_flags(threads, pin, gpu)}", cwd=scratch_dir)
        run_command(
            'echo "Potential\n0\n" | gmx energy -f em_steep.edr -o potential.xvg -xvg none',
            cwd=scratch_dir,
        )
        plot_svg("potential.xvg")

        run_command(
            f"gmx grompp -f {mdp_dirname}/NVT.mdp -c em_steep.gro -r em_steep.gro -p {m_gro_name}.top -o nvt.tpr -nobackup",
            cwd=scratch_dir,
        )
        run_command(f"gmx mdrun -deffnm nvt {_mdrun_flags(threads, pin, gpu)}", cwd=scratch_dir)
        run_command(
            'echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg -xvg none -b 20',
            cwd=scratch_dir,
        )
        plot_svg("temperature.xvg")

        run_command(
            f"gmx grompp -f {mdp_dirname}/NPT.mdp -c nvt.gro -r nvt.gro -p {m_gro_name}.top -o npt.tpr -nobackup",
            cwd=scratch_dir,
        )
        run_command(f"gmx mdrun -deffnm npt {_mdrun_flags(threads, pin, gpu)}", cwd=scratch_dir)
        run_command(
            'echo "Pressure" | gmx energy -f npt.edr -o pressure.xvg -xvg none',
            cwd=scratch_dir,
        )
        plot_svg("pressure.xvg")

        # Use equilibrated structure for REMD
        conf = "npt.gro" if os.path.isfile(os.path.join(scratch_dir, "npt.gro")) else "nvt.gro"

        # Temperatures
        if getattr(remd_cfg, "temperatures", None) and remd_cfg.spacing == "list":
            temps = list(remd_cfg.temperatures)
            if len(temps) != remd_cfg.nreplicas:
                raise ValueError("Length of remd.temperatures must equal remd.nreplicas")
        else:
            temps = _geometric_ladder(remd_cfg.nreplicas, remd_cfg.tmin, remd_cfg.tmax)

        # Template
        template = os.path.join(mdp_dirname, getattr(remd_cfg, "mdp_template", "REMD.mdp"))
        if not os.path.isfile(template):
            raise FileNotFoundError(f"REMD template MDP not found: {template}")

        # Write per-replica mdp + tpr (one subdirectory per replica for -multidir)
        with open("temperatures.txt", "w") as f:
            for t in temps:
                f.write(f"{t:.3f}\n")
        for i, T in enumerate(temps):
            rep_dir = f"rep_{i}"
            os.makedirs(rep_dir, exist_ok=True)
            # Store each replica's .mdp inside its own directory for clarity
            mdp_i = os.path.join(rep_dir, f"remd{i}.mdp")
            _make_replica_mdp(template, mdp_i, T, seed=1000 + i)
            run_command(
                f"gmx grompp -f {mdp_i} -c {conf} -p {m_gro_name}.top -o {rep_dir}/remd.tpr -maxwarn 1",
                cwd=scratch_dir,
            )

        # Write a run helper script using mpirun + -multidir with MPI/OpenMP knobs
        replex = getattr(remd_cfg, "replex", 100)
        nrep = remd_cfg.nreplicas
        with open("run_remd.sh", "w") as f:
            f.write(
                "\n".join(
                    [
                        "#!/usr/bin/env bash",
                        "set -euo pipefail",
                        'RUN_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
                        'cd "$RUN_DIR"',
                        "",
                        'GMX="${GMX:-gmx_mpi}"',
                        f'NREP="${{NREP:-{nrep}}}"',
                        f'REPLEX="${{REPLEX:-{replex}}}"',
                        'NTOMP="${NTOMP:-1}"',
                        'NP="${NP:-$NREP}"',
                        'EXTRA_MDRUN_FLAGS="${EXTRA_MDRUN_FLAGS:-}"',
                        'REP_PREFIX="${REP_PREFIX:-rep_}"',
                        "",
                        "rep_dirs=()",
                        'for i in $(seq 0 $((NREP-1))); do',
                        '    rep_dirs+=("${REP_PREFIX}${i}")',
                        "done",
                        "",
                        'echo "Using $GMX with:"',
                        'echo "  replicas (NREP)      = $NREP"',
                        'echo "  exchange every       = $REPLEX steps"',
                        'echo "  MPI ranks (NP)       = $NP"',
                        'echo "  OpenMP threads/rank  = $NTOMP"',
                        'echo "  replica dirs         = ${rep_dirs[*]}"',
                        'echo "  extra mdrun flags    = $EXTRA_MDRUN_FLAGS"',
                        "",
                        'exec mpirun -np "$NP" "$GMX" mdrun \\',
                        '    -v \\',
                        '    -s remd.tpr \\',
                        '    -deffnm remd \\',
                        '    -multidir "${rep_dirs[@]}" \\',
                        '    -replex "$REPLEX" \\',
                        '    -ntomp "$NTOMP" \\',
                        '    $EXTRA_MDRUN_FLAGS',
                        "",
                    ]
                )
            )
        os.chmod("run_remd.sh", 0o755)
