# amber14sb_corrected-Na.ff

Amber ff14SB GROMACS port with Joung–Cheatham TIP3P monovalent ions. Source:
`amber14sb.ff.tar.gz` from the [GROMACS user contributions](https://www.gromacs.org/user_contributions.html)
(dated 23 Dec 2015; archive SHA-256
`e77087de23d35b78ba5247e368613482289ca4b7ca5b9beb94cbc241bed6c0b6`). The
neighbouring `amber14sb.ff` folder is byte-identical to that source and records
the model used by simulations that reference it.

This folder differs from the source in one value: the `Na` Lennard-Jones epsilon
in `ffnonbonded.itp` is 3.65846e-01 kJ/mol. This is the Joung–Cheatham TIP3P
value (R_min/2 = 1.369 Å, ε = 0.0874393 kcal/mol; Joung & Cheatham,
J. Phys. Chem. B 112, 9020 (2008); `frcmod.ionsjc_tip3p` in AmberTools). The
source contains 3.65846e-02 kJ/mol, ten times too small. The same value is
distributed in `amber14sb_OL15.ff_corrected-Na-cation-params` (P. Jurečka, 2019).

All other parameters keep the source values: Joung–Cheatham K⁺ and Cl⁻, parm99
(Åqvist) Mg²⁺ and parm99 Ca²⁺ (`C0`), and TIP3P water. `forcefield.doc` is the
source file and names AMBER99SB-ILDN. The source's HID/HIP improper dihedral
assignment is unchanged; it affects protein systems only.

Topologies include this folder by name (`#include "amber14sb_corrected-Na.ff/..."`).
`tests/unit/test_forcefield_ion_parameters.py` checks the ion values and the
single difference from `amber14sb.ff`.
