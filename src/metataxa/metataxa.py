import argparse
from taxdistinct import *


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="A tool for computing taxonomic distinctness for inhomogeneous taxa information"
    )

    parser.add_argument(
        "--tax-file",
        type=str,
        required=True,
        help="File with taxonomy information for ORF ids",
    )
    parser.add_argument(
        "--pwy-file", type=str, required=True, help="File with pathways and ORF ids"
    )

    args = parser.parse_args()

    compute_pwy_tax_distance(args.tax_file, args.pwy_file)
