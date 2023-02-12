from nodes import NODE
from typing import Pattern, Tuple, List, Dict, Type, Set
from constants import *
from treeops import *

import re
import sys
import gzip
import numpy as np
from scipy import stats


def printf(fmt, *args):
    sys.stdout.write(fmt % args)
    sys.stdout.flush()


def _print_tree_(node: NODE, taxons) -> None:
    """
    Prints the sub-tree rooted at the node "node".

    :param node: The NODE for the node of the tree.
    :param taxons:  List for implementing stack operations of string used for printing the taxa in the taxa-tree.
    """

    taxons.append(node.name)

    if node.count >= 0:
        print(f"\t<{';'.join(taxons)}>\t{node.count}")
        pass

    for child_node in node.children.values():
        _print_tree_(child_node, taxons)

    taxons.pop()


def print_tree(node: NODE) -> None:
    """
    Prints the taxons and the counts in the tree.

    :param node: The root NODE for the tree.
    """
    _print_tree_(node, [])

def delete_pwy_orfs(pwy_orf: Dict[str, List[str]]) -> None:
    """
    Deletes the pathway to ORF mapping.

    :param pwy_orf: Map of pathways to the list of ORFs.
    """

    # delete the vectors allocated
    for pwy in pwy_orf:
        pwy_orf[pwy].clear()
    pwy_orf.clear()


def num_species_in_pwy(
    orf_taxon: Dict[str, str],
    pwy_orf: Dict[str, List[str]],
    pwy: str,
) -> int:
    """
    Computes the number of species in the pwy.

    :param orf_taxon: ORF-to-taxon map.
    :param pwy_orf: A map keyed by pathways to a  list of ORFs as  Dict[str, List[str]].
    :return: The root node for the tree created.
    """

    root_regex_patt: Pattern[str] = re.compile(r"root")
    taxid_regex_patt: Pattern[str] = re.compile(r"\([0-9]+\)$")

    species = set()

    if pwy in pwy_orf:
        for orfid in pwy_orf[pwy]:
            if orfid in orf_taxon:
                match_results = root_regex_patt.search(orf_taxon[orfid])
                if match_results:
                    taxon: str = re.sub(taxid_regex_patt, "", orf_taxon[orfid])
                    species.add(
                        ";".join([x for x in taxon.strip(" ").split(";") if x.strip()])
                    )
    return len(species)


def create_pwy_tree(
    orf_taxon: Dict[str, str],
    pwy_orf: Dict[str, List[str]],
    pwy: str,
    present_absent: bool,
) -> NODE:
    """
    Creates the tree for just one pathway.

    :param orf_taxon: ORF-to-taxon map.
    :param pwy_orf: A map keyed by pathways to a  list of ORFs as  Dict[str, List[str]].
    :return: The root node for the tree created.
    """

    root_regex_patt: Pattern[str] = re.compile(r"root")
    taxid_regex_patt: Pattern[str] = re.compile(r"\([0-9]+\)$")

    tree: NODE = NODE()

    if pwy in pwy_orf:
        for orfid in pwy_orf[pwy]:
            if orfid in orf_taxon:
                match_results = root_regex_patt.search(orf_taxon[orfid])
                if match_results:
                    taxon: str = re.sub(taxid_regex_patt, "", orf_taxon[orfid])
                    insert_node(
                        tree,
                        [x for x in taxon.strip(" ").split(";") if x.strip()],
                        0,
                        present_absent,
                    )
    return tree


def create_orf_tax_map(tax_file_name: str) -> Dict[str, str]:
    """
    Creates an ORF-to-taxon map from the functional and taxonomic file.

    :param tax_file_name: Input file name for taxonomy.
    :returns: A Dict[str, str] as a ORF to taxon map.
    """

    orf_tax_map: Dict[str, str] = dict()
    root_regex_patt: Pattern[str] = re.compile(r"root")
    taxid_regex_patt: Pattern[str] = re.compile(r"\([0-9]+\)$")
    # open a file to perform read operation using file object
    with gzip.open(tax_file_name, "rt") if tax_file_name.endswith(".gz") else open(
        tax_file_name, "r"
    ) as infile:

        for tp in infile:
            # read data from file object and put it into string.
            # print the data of the string

            taxstring: List[str] = tp.strip().split("\t")
            if len(taxstring) < 9:
                continue

            match_results = root_regex_patt.search(taxstring[8])
            if match_results:
                taxon: str = re.sub(taxid_regex_patt, "", taxstring[8])
                orf_tax_map[taxstring[0]] = taxon

    return orf_tax_map
