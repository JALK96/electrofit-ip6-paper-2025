import argparse

from electrofit.io.resp_edit import edit_resp_input


def main():
    # Initialize the argument parser
    parser = argparse.ArgumentParser(
        description="Edit RESP input file with equivalence groups."
    )

    # Define required positional arguments
    parser.add_argument("input_RESP_file", help="Path to the input RESP file")
    parser.add_argument(
        "equiv_groups_file", help="Path to the equivalence groups JSON file"
    )
    parser.add_argument("output_RESP_file", help="Path to save the output RESP file")

    # Define optional arguments
    parser.add_argument(
        "--ignore_sym",
        action="store_true",
        help="Ignore symmetry groups when editing RESP file",
    )

    # Parse the arguments
    args = parser.parse_args()

    # Call the edit_resp_input function with parsed arguments
    edit_resp_input(
        input_file=args.input_RESP_file,
        equiv_groups_file=args.equiv_groups_file,
        output_file=args.output_RESP_file,
        ignore_sym=args.ignore_sym,
    )


if __name__ == "__main__":
    main()
