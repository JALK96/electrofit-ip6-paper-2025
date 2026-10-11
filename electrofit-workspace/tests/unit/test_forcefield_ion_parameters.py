"""Ion Lennard-Jones parameters of the shipped Amber force-field folders."""
import math
import pathlib

DATA = pathlib.Path(__file__).resolve().parents[2] / "data"
CORRECTED = DATA / "amber14sb_corrected-Na.ff"
HISTORICAL = DATA / "amber14sb.ff"

# Joung & Cheatham, J. Phys. Chem. B 112, 9020 (2008), TIP3P set as in
# AmberTools frcmod.ionsjc_tip3p: R_min/2 in Angstrom, epsilon in kcal/mol.
JOUNG_CHEATHAM_TIP3P = {"Na": (1.369, 0.0874393), "K": (1.705, 0.1936829), "Cl": (2.513, 0.0355910)}


def _atomtypes(ff_dir):
    """Return {name: (sigma_nm, epsilon_kj_mol)} from the [ atomtypes ] block."""
    types, section = {}, None
    for raw in (ff_dir / "ffnonbonded.itp").read_text().splitlines():
        line = raw.split(";", 1)[0].strip()
        if line.startswith("["):
            section = line.strip("[] ").lower()
        elif line and section == "atomtypes":
            fields = line.split()
            types[fields[0]] = (float(fields[-2]), float(fields[-1]))
    return types


def _gromacs_units(rmin_half_angstrom, epsilon_kcal):
    sigma_nm = 2 * rmin_half_angstrom / 2 ** (1 / 6) / 10
    return sigma_nm, epsilon_kcal * 4.184


def test_corrected_folder_has_joung_cheatham_monovalent_ions():
    types = _atomtypes(CORRECTED)
    for name, (rmin_half, eps_kcal) in JOUNG_CHEATHAM_TIP3P.items():
        sigma, epsilon = _gromacs_units(rmin_half, eps_kcal)
        assert math.isclose(types[name][0], sigma, rel_tol=1e-5), name
        assert math.isclose(types[name][1], epsilon, rel_tol=1e-5), name


def test_corrected_folder_differs_from_historical_only_in_na_epsilon():
    for historical_file in sorted(p for p in HISTORICAL.rglob("*") if p.is_file()):
        relative = historical_file.relative_to(HISTORICAL)
        old = historical_file.read_bytes().splitlines()
        new = (CORRECTED / relative).read_bytes().splitlines()
        changed = [(a, b) for a, b in zip(old, new) if a != b]
        assert len(old) == len(new), relative
        if relative.as_posix() == "ffnonbonded.itp":
            assert len(changed) == 1
            assert changed[0][0].split()[0] == changed[0][1].split()[0] == b"Na"
            assert _atomtypes(HISTORICAL)["Na"] == (2.43928e-01, 3.65846e-02)
        else:
            assert not changed, relative
