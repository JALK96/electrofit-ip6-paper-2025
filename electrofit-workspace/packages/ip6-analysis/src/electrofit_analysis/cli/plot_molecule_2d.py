import os
import argparse
import logging

from electrofit.io.files import find_file_with_extension
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

logger = logging.getLogger(__name__)


def main(project_path: str, subdir: str = "process") -> None:
    """Create 2D depictions (SVG) for molecules found under project_path/subdir.

    The script searches each microstate directory for a MOL2 file in
    'run_gau_create_gmx_in', then renders a labeled 2D structure into
    'analyze_final_sim/molecule.svg'.
    """
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

    process_dir = os.path.join(project_path, subdir)
    if not os.path.isdir(process_dir):
        logger.error("Process directory not found: %s", process_dir)
        return

    for folder_name in os.listdir(process_dir):
        folder_path = os.path.join(process_dir, folder_name)
        if not os.path.isdir(folder_path):
            continue

        run_final_sim_dir = os.path.join(folder_path, "run_final_gmx_simulation")
        run_gau_create_gmx_in_dir = os.path.join(folder_path, "run_gau_create_gmx_in")

        if not os.path.isdir(run_final_sim_dir):
            continue

        dest_dir = os.path.join(folder_path, "analyze_final_sim")
        os.makedirs(dest_dir, exist_ok=True)

        prev_cwd = os.getcwd()
        try:
            # Locate MOL2 file
            os.chdir(run_gau_create_gmx_in_dir)
            mol2_file_name = find_file_with_extension("mol2")
            mol2_file = os.path.join(run_gau_create_gmx_in_dir, mol2_file_name)

            # Load the molecule
            molecule = Chem.MolFromMol2File(mol2_file, removeHs=False)
            if molecule is None:
                logger.warning("Failed to load molecule from %s", mol2_file)
                continue
            logger.info("Loaded molecule for %s", folder_name)

            # Coordinates and depiction
            rdDepictor.SetPreferCoordGen(True)
            rdDepictor.Compute2DCoords(molecule)

            # Custom atom labels
            atom_counters: dict[str, int] = {}
            atom_labels: dict[int, str] = {}
            for atom in molecule.GetAtoms():
                symbol = atom.GetSymbol()
                idx = atom.GetIdx()
                atom_counters[symbol] = atom_counters.get(symbol, 0) + 1
                atom_labels[idx] = f"{symbol}{atom_counters[symbol]}"

            # Draw
            svg_size = 500
            drawer = rdMolDraw2D.MolDraw2DSVG(svg_size, svg_size)
            opts = drawer.drawOptions()
            opts.addAtomIndices = False
            opts.addBondIndices = False
            opts.baseFontSize = 0.3
            for idx, label in atom_labels.items():
                opts.atomLabels[idx] = label
            rdMolDraw2D.PrepareMolForDrawing(molecule)
            drawer.DrawMolecule(molecule)
            drawer.FinishDrawing()
            svg = drawer.GetDrawingText()

            out_svg = os.path.join(dest_dir, "2d-structure.svg")
            with open(out_svg, "w") as f:
                f.write(svg)
            logger.info("Wrote %s", out_svg)
        except Exception as e:
            logger.exception("Failed processing %s: %s", folder_name, e)
        finally:
            os.chdir(prev_cwd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Render 2D molecule depictions for IP6 microstates.")
    parser.add_argument(
        "--project",
        "-p",
        required=True,
        help="Path to the project root directory (contains data/todo).",
    )
    parser.add_argument(
        "--subdir",
        default="process",
        help="Relative subdirectory under project root to traverse (default: process).",
    )
    args = parser.parse_args()
    main(os.path.abspath(args.project), args.subdir)
