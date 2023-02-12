from nodes import NODE
from typing import Pattern, Tuple, List, Dict, Type, Set
from constants import *
from utils import *
from iofns import *
from treeops import *
from estimates import *


import re
import sys
import gzip
import numpy as np
from scipy import stats


def subsample_tree_count(node: NODE, p: float) -> None:
    """
    Updates the count variables for each node by generating a new count from a binomial(actual_count, p)
    sampling, where p is the percentage of a samples to retain during the subsampling.

    :param node: A node in the tree taxonomic tree.
    :param p: The percentage of counts to  subsample.
    """
    if node.actual_count:
        node.count = np.random.binomial(node.actual_count, p, 1)[0]
        # node.count = node.true_count
        # print(len(node.children), node.count, node.true_count)
    else:
        node.count = 0

    for child_node in node.children.values():
        # print(child_node.name)
        subsample_tree_count(child_node, p)


def compute_tax_distance(tax_file_name: str) -> None:
    """
    Computes the taxonomic distance for the data in the taxonomy file.

    :param tax_file_name: The input file with the list of taxonomic annotations in the sample.
    """

    # read the file
    tree = read_tax_file(tax_file_name)

    # count the subtree and update in each  node
    sum_tree_count(tree)

    print(f"Subtree count at root: {subtree_count(tree)}")

    print(f"Delta*  {compute_delta(tree, DELTA_STAR, false)}")


def estimate_delta_confidence_interval(
    orf_taxon, pwy_orfs, pwy, delta_type, use_wtd, n=100, p=0.9
):

    switch = {DELTA: False, DELTA_STAR: False, DELTA_PLUS: True}
    present_absent: bool = switch.get(delta_type, True)

    root_node = create_pwy_tree(orf_taxon, pwy_orfs, pwy, present_absent=present_absent)
    # present_absent =  false
    sum_tree_count(root_node)
    subtree_count(root_node)
    delta = compute_delta(root_node, delta_type, use_wtd=use_wtd)

    # save the origin count
    save_tree_count(root_node)
    np.random.seed(1)
    deltas = []
    for i in range(0, n):
        subsample_tree_count(root_node, p)
        sum_tree_count(root_node)
        subtree_count(root_node)

        deltas.append(compute_delta(root_node, delta_type, use_wtd=use_wtd))
        # print_tree(root_node)
    #    delete_tree(root_node)
    confidence_interval = stats.t.interval(
        alpha=0.95, df=len(deltas) - 1, loc=np.mean(deltas), scale=stats.sem(deltas)
    )

    return delta, confidence_interval


def compute_pwy_tax_distance(tax_file_name: str, pwy_orf_filename: str) -> None:
    """
    Computes the taxonomic distance for the data in the taxonomy file for the give pwy-to-ORF map file.

    :param tax_filename: The input file with the list of taxonomic annotations in the sample.
    :param pwy_orf_filename: The input file with the list of ORFs for individual pathways in the sample.
    """

    # read the file
    orf_taxon: Dict[str, str] = create_orf_tax_map(tax_file_name)
    pwy_orfs: Dict[str, List[str]] = read_pwy_file(pwy_orf_filename)

    headers = [
        [
            "PWY",
            "num_orfs",
            "-",
            "Delta",
            "-",
            "-",
            "Delta*",
            "-",
            "-",
            "Delta+",
            "-",
            "-",
            "Delta (WTD)",
            "-",
            "-",
            "Delta* (WTD)",
            "-",
            "-",
            "Delta+ (WTD)",
            "-",
        ],
        [
            "-",
            "-",
            "mu",
            "CI",
            "-",
            "mu",
            "CI",
            "-",
            "mu",
            "CI",
            "-",
            "mu",
            "CI",
            "-",
            "mu",
            "CI",
            "-",
            "mu",
            "CI",
            "-",
        ],
        [
            "-",
            "-",
            "-",
            "low",
            "hi",
            "-",
            "low",
            "hi",
            "-",
            "low",
            "hi",
            "-",
            "low",
            "hi",
            "-",
            "low",
            "hi",
            "-",
            "low",
            "hi",
        ],
    ]

    for header in headers:
        print("\t".join([f"{x:^8}" for x in header]))

    np.random.seed(1)
    outputs = []
    for pwy, orflist in pwy_orfs.items():
        # if pwy!='PWY-15':
        #   continue

        num_species = num_species_in_pwy(orf_taxon, pwy_orfs, pwy)
        pwy_output = [pwy]
        pwy_output.append(len(orflist))

        # Delta
        results, confidence_interval = estimate_delta_confidence_interval(
            orf_taxon, pwy_orfs, pwy, DELTA, False, n=100, p=0.9
        )
        pwy_output.append(results)
        pwy_output = pwy_output + list(confidence_interval)

        # Delta*
        results, confidence_interval = estimate_delta_confidence_interval(
            orf_taxon, pwy_orfs, pwy, DELTA_STAR, False, n=100, p=0.9
        )
        pwy_output.append(results)
        pwy_output = pwy_output + list(confidence_interval)

        # Delta+
        results, confidence_interval = estimate_delta_confidence_interval(
            orf_taxon, pwy_orfs, pwy, DELTA_PLUS, False, n=100, p=0.9
        )
        pwy_output.append(results)
        pwy_output = pwy_output + list(confidence_interval)

        # Delta (WTD)
        results, confidence_interval = estimate_delta_confidence_interval(
            orf_taxon, pwy_orfs, pwy, DELTA, True, n=100, p=0.9
        )
        pwy_output.append(results)
        pwy_output = pwy_output + list(confidence_interval)

        # Delta* (WTD)
        results, confidence_interval = estimate_delta_confidence_interval(
            orf_taxon, pwy_orfs, pwy, DELTA_STAR, True, n=100, p=0.9
        )
        pwy_output.append(results)
        pwy_output = pwy_output + list(confidence_interval)

        # Delta+ (WTD)
        results, confidence_interval = estimate_delta_confidence_interval(
            orf_taxon, pwy_orfs, pwy, DELTA_PLUS, True, n=100, p=0.9
        )
        pwy_output.append(results)
        pwy_output = pwy_output + list(confidence_interval)

        #print(pwy, len(orflist), results, confidence_interval)
        outputs.append(pwy_output)

    for output in outputs:
        fields = [] 
        for field in output:
            if isinstance(field, float):
                fields.append(f"{field:<.6f}") 
            else: 
                fields.append(f"{field:<8}") 

        print("\t".join(fields))
