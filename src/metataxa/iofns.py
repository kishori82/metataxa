from nodes import NODE
from typing import Pattern, Tuple, List, Dict, Type, Set
from constants import *

import re
import sys
import gzip


def read_tax_file(tax_file_name: str) -> NODE:
    """
    Reads the taxa file

    :param tax_file_name: Input file name for the  taxa in the sample
    :return: The root NODE for the taxa in the file
    """

    print(f"TODO: Make taxa names case insensitive")
    root_regex_patt: Pattern[str] = re.compile(r"root")
    taxid_regex_patt: Pattern[str] = re.compile(r"\([0-9]+\)$")

    tree: NODE = NODE()
    taxstring: List[str] = []

    with open(tax_file_name, "rb") if newfile.endswith(".gz") else open(
        tax_file_name, "r"
    ) as newfile:
        # open a file to perform read operation using file object

        print("TODO: make the taxstring columns and count flexible")
        k: int = 0
        # read data from file object and put it into string.
        for tp in newfile.readline():
            # print the data of the string
            # print(tp)
            taxstring: List[str] = tp.strip().split("\t")
            if len(taxstring) >= 10:
                strmatch_result = root_regex_patt.search(taxstring[8])
                if strmatch_result:
                    taxon: [str] = re.sub(taxid_regex_patt, "", taxstring[8])
                    """
                   for(auto x: split(trim(taxon, " "), ';') ) {
                   std::cout << '\t' << "<" <<  x << ">" << std::endl;
                   std::cout << taxon << std::endl;
                  """
                    k += 1
                    insert_node(
                        tree,
                        [x for x in taxon.trim(" ").split(";") if x.strip()],
                        0,
                        false,
                    )
    print(f"No of nodes inserted {k}")
    return tree


def read_pwy_file(pwy_filename: str) -> Dict[str, List[str]]:
    """
    Creates a map from pwy name to a vector of ORF names.

    :param pwy_filename: File name with the pathway and corresponding ORFs in each pathway in the sample.
    """

    pwy_orf: Dict[str, Set[str]] = dict()

    # open a file to perform read operation using file object
    with gzip.open(pwy_filename, "rt") if pwy_filename.endswith(".gz") else open(
        pwy_filename, "r"
    ) as infile:
        for tp in infile:
            # read data from file object and put it into string.

            # print the data of the string
            # print(tp)

            orfstring: List[str] = tp.strip().split("\t")
            if len(orfstring) >= 2:
                pwy_orf[orfstring[0]] = {orf for orf in orfstring[1:]}

    print(f"No of pwys inserted {len(pwy_orf)}")
    return pwy_orf


def delete_pwy_orfs(pwy_orf: Dict[str, List[str]]) -> None:
    """
    Deletes the pathway to ORF mapping.

    :param pwy_orf: Map of pathways to the list of ORFs.
    """

    # delete the vectors allocated
    for pwy in pwy_orf:
        pwy_orf[pwy].clear()
    pwy_orf.clear()
