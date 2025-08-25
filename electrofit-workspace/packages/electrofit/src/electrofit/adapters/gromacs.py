"""GROMACS Adapter.

Ehemals unter ``electrofit.external.gromacs`` – Pfad vollständig entfernt.
"""
from __future__ import annotations

import logging
import os
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
from electrofit.infra.scratch_manager import setup_scratch_directory

__all__ = ["set_up_production"]


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


def set_up_production(
	m_gro,
	MDP_dir,
	base_scratch_dir,
	molecule_name,
	box_type="dodecahedron",
	cation="NA",
	anion="CL",
	d="1.2",
	conc="0.15",
	ff="amber14sb.ff",
	*,
	threads: int | None = None,
	pin: bool | None = None,
):
	"""Set up und führe eine Produktions‑MD mit GROMACS aus (migrierte Version)."""
	fullpath = os.getcwd()
	log_file_path = os.path.join(fullpath, "process.log")
	suppress = any(isinstance(h, logging.FileHandler) for h in logging.getLogger().handlers)
	setup_logging(log_file_path, suppress_initial_message=suppress)

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
	scratch_dir, original_dir = setup_scratch_directory(input_files, base_scratch_dir)

	with ensure_finalized(original_dir=original_dir, scratch_dir=scratch_dir, input_files=input_files):
		os.chdir(scratch_dir)

		run_command(
			f"echo 0 | gmx editconf -f {m_gro} -o {m_gro_name}_box.gro -bt {box_type} -d {d} -c -princ -nobackup",
			cwd=scratch_dir,
		)
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
		run_command(
			f'echo "SOL" | gmx genion -s ions.tpr -o {m_gro_name}_tip3p_ions.gro -conc {conc} -p {m_gro_name}.top -pname {cation} -nname {anion} -neutral',
			cwd=scratch_dir,
		)
		run_command(f"tail {m_gro_name}.top", cwd=scratch_dir)

		# -------- Decide thread flags ONCE and reuse for all mdrun calls --------
		try:
			import psutil
			logical_cpus = psutil.cpu_count(logical=True) or (os.cpu_count() or 1)
		except Exception:
			logical_cpus = os.cpu_count() or 1

		if threads is not None and threads > logical_cpus:
			logging.warning("Requested threads=%d exceeds logical CPUs=%d -> clamping",
							threads, logical_cpus)
			threads = logical_cpus

		def _mdrun_flags() -> str:
			flags = []
			if threads is not None:
				flags += ["-nt", str(threads)]
				# For tiny systems, forcing single rank avoids DD issues on some builds
				if threads == 1:
					flags += ["-ntmpi", "1"]
			if pin is not None:
				flags += ["-pin", "on" if pin else "off"]
			flags.append("-nobackup")
			return " ".join(flags)

		run_command(
			f"gmx grompp -f {mdp_dirname}/em_steep.mdp -c {m_gro_name}_tip3p_ions.gro -p {m_gro_name}.top -o em_steep.tpr",
			cwd=scratch_dir,
		)
		run_command(f"gmx mdrun -deffnm em_steep {_mdrun_flags()}", cwd=scratch_dir)
		run_command(
			'echo "Potential\n0\n" | gmx energy -f em_steep.edr -o potential.xvg -xvg none',
			cwd=scratch_dir,
		)
		plot_svg("potential.xvg")

		run_command(
			f"gmx grompp -f {mdp_dirname}/NVT.mdp -c em_steep.gro -r em_steep.gro -p {m_gro_name}.top -o nvt.tpr -nobackup",
			cwd=scratch_dir,
		)
		run_command(f"gmx mdrun -deffnm nvt {_mdrun_flags()}", cwd=scratch_dir)
		run_command(
			'echo "Temperature" | gmx energy -f nvt.edr -o temperature.xvg -xvg none -b 20',
			cwd=scratch_dir,
		)
		plot_svg("temperature.xvg")
		run_command(
			f"gmx grompp -f {mdp_dirname}/NPT.mdp -c nvt.gro -r nvt.gro -p {m_gro_name}.top -o npt.tpr -nobackup",
			cwd=scratch_dir,
		)
		run_command(f"gmx mdrun -deffnm npt {_mdrun_flags()}", cwd=scratch_dir)
		run_command(
			'echo "Pressure" | gmx energy -f npt.edr -o pressure.xvg -xvg none',
			cwd=scratch_dir,
		)
		plot_svg("pressure.xvg")

		# Pick latest equilibrated state: prefer NPT if available, else NVT
		prod_conf = "npt.gro" if os.path.isfile(os.path.join(scratch_dir, "npt.gro")) else "nvt.gro"
		prod_cpt  = "npt.cpt" if prod_conf.startswith("npt") and os.path.isfile(os.path.join(scratch_dir, "npt.cpt")) else "nvt.cpt"
		if prod_conf.startswith("npt"):
			logging.info("[prod] Using NPT output (%s, %s) for production start", prod_conf, prod_cpt)
		else:
			logging.warning("[prod] NPT outputs missing -> falling back to NVT (%s, %s)", prod_conf, prod_cpt)
		
		run_command(
			f"gmx grompp -f {mdp_dirname}/Production.mdp -c {prod_conf} -t {prod_cpt} -p {m_gro_name}.top -o md.tpr",
			cwd=scratch_dir,
		)
		run_command(f"gmx mdrun -deffnm md {_mdrun_flags()}", cwd=scratch_dir)
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