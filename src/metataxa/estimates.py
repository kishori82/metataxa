from typing import Pattern, Tuple, List, Dict, Type, Set
from nodes import *
from constants import *
from treeops import *


def _sum_tree_count_(node: NODE, _sum: int) -> None:
    """
    Computes the sum of all the taxon nodes in the sub-tree rooted at node.

    :param node: A node in the tree taxonomic tree.
    :param _sum: The _sum variable retains the value/result after the fuction call..
    """
    _sum = _sum + node.count
    for taxon, child_node in node.children.items():
        _sum_tree_count_(child_node, _sum)


def sum_tree_count(node: NODE) -> int:
    """
    Computes the sum of all the counts of taxon nodes in the tree rooted at node

    :param node:  A node in the tree
    """

    _sum: int = 0
    #  put the value of the sum in the sum variable
    _sum_tree_count_(node, _sum)
    return _sum


def subtree_count(node: NODE) -> int:
    """
    Updates the subtree_count variable for each NODE in the sub-tree rooted at "node".

    :param node: Any NODE in the tree.
    """

    _sum: int = node.count

    for child_node in node.children.values():
        _sum = _sum + subtree_count(child_node)

    node.subtree_count = _sum

    return _sum


def _count_edge_crossings_(node: NODE, tot_count: int, use_wtd: bool) -> float:
    """
    Counts the total number of edge crossings from the various species. This can also compute with
    WTD consideration as weights on the edges.

    :param node: The root node of the tree/sub-tree for which the edge-corossing is counted.
    :param total_count: The total count of the tree at that node.
    :param use_wtd: Use WTD edge weighting if it is True.
    :return: Number of edge-crossings in the tree.
    """
    _sum: float = 0
    for child_node in node.children.values():
        if use_wtd:
            _sum = _sum + (
                tot_count - child_node.subtree_count
            ) * child_node.subtree_count * pow(0.5, child_node.depth)
        else:
            _sum = (
                _sum + (tot_count - child_node.subtree_count) * child_node.subtree_count
            )
        _sum = _sum + _count_edge_crossings_(child_node, tot_count, use_wtd)

        # print(f"use_wtd: {use_wtd}   _sum: {_sum}   child_node: {child_node.name}\n")

    return _sum


def __get_count_vector(node: NODE, count_vector: List[int]) -> None:
    """
    Returns the vector of counts for all nodes in the tree.

    :param node: The root node of the tree.
    :param count_vector: The vector of the taxonomic counts for the nodes in the sub-tree rooted in the "node".
    """
    if node.count:
        count_vector.append(node.count)

    for child_node in node.children.values():
        __get_count_vector(child_node, count_vector)


def compute_delta(node: NODE, delta_type: Type[DeltaType], use_wtd: bool) -> float:
    """
    Computes the delta star for the tree rooted in node.

    :param node: Node for the subtree to compute delta star.
    :param delta_type: Type of the delta.
    :param use_wtd: Flag to indicate whether WTD should be used in the calculation of the edge corssings.
    :return: A double representing the computed delta.
    """

    tot_count: int = node.subtree_count
    prod_count: float = _count_edge_crossings_(node, tot_count, use_wtd)

    # print(f"\nprod_count {prod_count}")
    count_vector: List[int] = []
    __get_count_vector(node, count_vector)

    # sum-product for all distinct counts for the taxons
    delta: float = 0
    if len(count_vector) < 2:
        return 0.0

    if delta_type == DELTA_STAR:
        sum_product: int = 0
        for i in range(0, len(count_vector)):
            for j in range(i + 1, len(count_vector)):
                sum_product += count_vector[i] * count_vector[j]

        delta = float(prod_count) / float(sum_product)
    else:
        N: int = sum(count_vector)
        nums: float = float(N) * (float(N) - 1.0) / 2.0
        delta = float(prod_count) / nums

    # print(f"prod_count {prod_count}  {sum_product}")
    return delta

def Delta(taxa: List[str], delta_type: Type[DeltaType], use_wtd: bool) -> float:
    """
    Computes the  Delta for the taxa with or without WTD consideration.

    :param taxa: List of the taxa name.
    :param delta_type: Type of the delta.
    :param use_wtd: Flag to indicate whether adjustments of the edge weights are adjusted to WTD.
    :return: A double representing the computed delta.
    """

    switch = {DELTA: False, DELTA_STAR: False, DELTA_PLUS: True}
    present_absent: bool = switch.get(delta_type, True)

    root_node: NODE = create_tree(taxa, present_absent)

    sum_tree_count(root_node)
    subtree_count(root_node)

    delta: float = compute_delta(
        root_node, delta_type, use_wtd=use_wtd
    )  # use_wtd = use_wtd
    delete_tree(root_node)

    return delta


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



