import os
import sys
import logging

from electrofit.io.mol2_ops import update_mol2_charges, Mol2ChargeError

def main():
    # Check if the correct number of arguments is provided
    if len(sys.argv) != 4:
        print(
            "Usage: python update_mol2_charges.py <input_MOL2_file> <chg_file> <output_MOL2_file>"
        )
        print(
            "Example: python update_mol2_charges.py mol2_input.mol2 charges.chg fixed_out.mol2"
        )
        sys.exit(1)

    # Assign command-line arguments to variables
    mol2_file = sys.argv[1]
    chg_file = sys.argv[2]
    mol2_output = sys.argv[3]

    # Validate input file paths
    if not os.path.isfile(mol2_file):
        raise Mol2ChargeError(f"Input MOL2 file '{mol2_file}' does not exist")
    if not os.path.isfile(chg_file):
        raise Mol2ChargeError(f"Charges file '{chg_file}' does not exist")

    # Ensure the output directory exists
    output_dir = os.path.dirname(os.path.abspath(mol2_output))
    if output_dir and not os.path.exists(output_dir):
        raise Mol2ChargeError(f"Output directory '{output_dir}' does not exist")

    # Update the MOL2 charges
    try:
        update_mol2_charges(mol2_file, chg_file, mol2_output)
    except Mol2ChargeError as e:
        logging.error("mol2 charge update failed: %s", e)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
