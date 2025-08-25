"""RESP adapter (skeleton).

Encapsulates generation & execution of RESP charge fitting stages.

Minimal surface (initial):
- prepare_ac(mol2_file, name, net_charge, cwd) -> ac_file
- generate_inputs(name, cwd) -> (resp1, resp2)
- apply_symmetry(resp1, scratch_dir, adjust: bool, ignore: bool) -> resp1_effective
- run_two_stage(name, esp_file, resp1, resp2, cwd) -> mol2_resp

Later: unify with prep vs conformer pipelines, expose structured results / timings.
"""
from __future__ import annotations
import os
import time
from .results import RespStageResult

__all__ = [
	"prepare_ac",
	"generate_inputs",
	"apply_symmetry",
	"run_two_stage",
	"run_two_stage_with_result",
]


def prepare_ac(mol2_file: str, name: str, net_charge: int, cwd: str) -> str:
	from electrofit.cli.run_commands import run_command
	from electrofit.io.files import replace_charge_in_ac_file
	run_command(f"bondtype -i {mol2_file} -o {name}.ac -f mol2 -j part", cwd=cwd)
	replace_charge_in_ac_file(file_path=f"{name}.ac", new_charge_float=net_charge, cwd=cwd)
	return f"{name}.ac"


def generate_inputs(name: str, cwd: str) -> tuple[str, str]:
	from electrofit.cli.run_commands import run_command
	run_command(f"respgen -i {name}.ac -o ANTECHAMBER_RESP1.IN -f resp1", cwd=cwd)
	run_command(f"respgen -i {name}.ac -o ANTECHAMBER_RESP2.IN -f resp2", cwd=cwd)
	return ("ANTECHAMBER_RESP1.IN", "ANTECHAMBER_RESP2.IN")


def apply_symmetry(scratch_dir: str, adjust_sym: bool, ignore_sym: bool) -> str:
	if not adjust_sym:
		return os.path.join(scratch_dir, "ANTECHAMBER_RESP1.IN")
	from electrofit.domain.symmetry.resp_constraints import apply_and_optionally_modify, RespSymmetryConfig
	from electrofit.cli.run_commands import run_python
	from electrofit.io.files import find_file_with_extension
	from electrofit.io.resp_edit import edit_resp_input
	from electrofit.domain.symmetry.write_symmetry import write_symmetry
	return apply_and_optionally_modify(
		run_python,
		find_file_with_extension,
		edit_resp_input,
		write_symmetry,
		RespSymmetryConfig(scratch_dir=scratch_dir, adjust_sym=adjust_sym, ignore_sym=ignore_sym),
	)


def run_two_stage(name: str, esp_file: str, resp1: str, resp2: str, cwd: str) -> str:
	from electrofit.cli.run_commands import run_resp, run_python
	from electrofit.io.mol2_ops import update_mol2_charges
	resp_output1 = f"{name}-resp1.out"
	resp_pch1 = f"{name}-resp1.pch"
	resp_chg1 = f"{name}-resp1.chg"
	run_resp(resp1, esp_file, resp_output1, resp_pch1, resp_chg1, cwd)
	resp_output2 = f"{name}-resp2.out"
	resp_pch2 = f"{name}-resp2.pch"
	resp_chg2 = f"{name}-resp2.chg"
	run_resp(resp2, esp_file, resp_output2, resp_pch2, resp_chg2, cwd, resp_chg1)
	mol2_resp = f"{name}_resp.mol2"
	run_python(update_mol2_charges, f"{name}.mol2", resp_chg2, mol2_resp, cwd=cwd)
	return mol2_resp


def run_two_stage_with_result(name: str, esp_file: str, resp1: str, resp2: str, cwd: str, symmetry_applied: bool) -> RespStageResult:
	"""Timing + strukturierte Metadaten (nicht von Domain genutzt)."""
	t0 = time.perf_counter()
	mol2_resp = run_two_stage(name, esp_file, resp1, resp2, cwd)
	t1 = time.perf_counter()
	return RespStageResult(
		name=name,
		resp1_in=resp1,
		resp2_in=resp2,
		resp1_out=f"{name}-resp1.out",
		resp2_out=f"{name}-resp2.out",
		final_mol2=mol2_resp,
		charges_file=f"{name}-resp2.chg",
		symmetry_applied=symmetry_applied,
		wall_time_s=t1 - t0,
	)
