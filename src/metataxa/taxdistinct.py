from collections.abc import Mapping
from nodes import NODE
from typing import Pattern, Tuple, List, Dict, Type
from constants import DeltaType
import re


def insert_node(
    node: NODE, tax_hierarchy: List[str], depth: int = 0, present_absent: bool = False
) -> None:
    """
    Inserts a node in to the tree based on the name of the taxon
    to create the taxonomy tree.

    :param node: A NODE to start creating the taxonomy tree because it is implemented in a recursive manner.
    :param tax_hierarchy: The vector of strings that encodes the taxonomic hierachy, e.g., "<taxA, taxB, taxC, taxD> for taxA;taxB;taxC;taxD
    :param depth: The depth of the taxa hierarchy we can insert into the tree.

    """
    n: NODE = None

    if node == None:
        node = NODE()
        node.name = tax_hierarchy[depth]
        node.depth = depth

    if len(tax_hierarchy) == 1:
        if present_absent:
            node.count = 0
        node.count += 1

    if len(tax_hierarchy) > depth + 1:
        depth = depth + 1

        if tax_hierarchy[depth] in node.children:
            n = NODE()
            n.name = tax_hierarchy[depth]
            n.depth = depth
            node.children.append((tax_hierarchy[depth], n))
            print(f"\t{n.name} : {depth}")
        else:
            n = node.children[tax_hierarchy[depth]]

        if len(tax_hierarchy) == depth + 1:
            if present_absent:
                n.count = 1
            else:
                n.count += 1

        # print(f"\t{depth}")
        insert_node(n, tax_hierarchy, depth, present_absent)


def _sum_tree_count_(node: NODE, _sum: int) -> None:
    """
    Computes the sum of all the taxon nodes in the sub-tree rooted at node.

    :param node: A node in the tree taxonomic tree.
    :param _sum: The _sum variable retains the value/result after the fuction call..
    """
    sum = sum + node.count
    for taxon, child_node in node.children.items():
        _sum_tree_count_(child_node, sum)


def sum_tree_count(node: NODE) -> int:
    """
    Computes the sum of all the counts of taxon nodes in the tree rooted at node

    :param node:  A node in the tree
    """

    sum: int = 0
    #  put the value of the sum in the sum variable
    _sum_tree_count_(node, sum)
    return sum


def subtree_count(node: NODE) -> int:
    """
    Updates the subtree_count variable for each NODE in the sub-tree rooted at "node".

    :param node: Any NODE in the tree.
    """

    _sum: int = node.count

    for taxon, child_node in node.children.items():
        _sum = _sum + subtree_count(child_node)

    node.subtree_count = _sum

    return _sum


def ___delete_tree(node: NODE) -> None:
    """
    Recursively deletes the nodes in the sub-tree rooted in node ``node".

    :param node: Root NODE of the sub-tree to start deleting from.
    """

    for child_node in node.children:
        ___delete_tree(child_node)

    node.children.clear()


def delete_tree(node: NODE) -> None:
    """
    Deletes the nodes in the sub-tree rooted in node.

    :param node: Root NODE of the tree .
    """

    __delete_tree(node)
    node = None


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
    for child_node in node.children:
        if use_wtd:
            _sum = _sum + (tot_count - child_node[1].subtree_count) * child_node[
                1
            ].subtree_count * pow(0.5, child_node[1].depth)
        else:
            _sum = (
                _sum
                + (tot_count - child_node[1].subtree_count)
                * child_node[1].subtree_count
            )
        _sum = _sum + _count_edge_crossings_(child_node[1], tot_count, use_wtd)

    return sum


def __get_count_vector(node: NODE, count_vector: List[int]) -> None:
    """
    Returns the vector of counts for all nodes in the tree.

    :param node: The root node of the tree.
    :param count_vector: The vector of the taxonomic counts for the nodes in the sub-tree rooted in the "node".
    """

    if node.count:
        count_vector.append(node.count)

    for child_node in node.children:
        __get_count_vector(child_node[1], count_vector)


def compute_delta(node: NODE, delta_type: Type[DeltaType], use_wtd: bool) -> float:
    """
    Computes the delta star for the tree rooted in node.

    :param node: Node for the subtree to compute delta star.
    :param delta_type: Type of the delta.
    :param use_wtd: Flag to indicate whether WTD should be used in the calculation of the edge corssings.
    :return: A double representing the computed delta.
    """

    tot_count: int = node.subtree_count
    prod_count: double = _count_edge_crossings_(node, tot_count, use_wtd)

    count_vector: List[int] = []

    __get_count_vector(node, count_vector)

    # sum-product for all distinct counts for the taxons
    delta: float = 0
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


def _print_tree_(node: NODE, taxons: List[str]) -> None:
    """
    Prints the sub-tree rooted at the node "node".

    :param node: The NODE for the node of the tree.
    :param taxons:  List for implementing stack operations of string used for printing the taxa in the taxa-tree.
    """

    taxons.append(node.name)

    if node.count >= 0:
        # print(f"\t<{taxon_concat(taxons)}>\t{node.count}\t{taxon}")
        pass

    for child_node in node.children.values():
        _print_tree_(child_node, taxons)

    taxons.pop()


def print_tree(node: NODE) -> None:
    """
    Prints the taxons and the counts in the tree.

    :param node: The root NODE for the tree.
    """
    taxons: List[str] = []
    _print_tree_(node, taxons)


def read_tax_file(tax_file_name: str) -> NODE:
    """
    Reads the taxa file

    :param tax_file_name: Input file name for the  taxa in the sample
    :return: The root NODE for the taxa in the file
    """

    print(f"TODO: Make taxa names case insensitive")
    root_regex_patt: Pattern[str] = re.compile("root")
    count_regex_patt: Pattern[str] = re.compile("(\\([0-9]+\\)$)")

    tree: NODE = None
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
            taxstring: List[str] = tp.split("\t")
            if len(taxstring) >= 10:
                strmatch_result = re.search(taxstring[8], root_regex_patt)
                if strmatch_result:
                    taxon: [str] = re.sub(taxstring[8], count_regex_patt, "")
                    """
                   for(auto x: split(trim(taxon, " "), ';') ) {
                   std::cout << '\t' << "<" <<  x << ">" << std::endl;
                   std::cout << taxon << std::endl;
                  """
                    k += 1
                    insert_node(tree, taxon.trim(" ").split(";"), 0, false)
    print(f"No of nodes inserted {k}")
    return tree


def create_orf_tax_map(tax_file_name: str) -> Dict[str, str]:
    """
    Creates an ORF-to-taxon map from the functional and taxonomic file.

    :param tax_file_name: Input file name for taxonomy.
    :returns: A Dict[str, str] as a ORF to taxon map.
    """

    orf_tax_map: Dict[str, str] = dict()

    # open a file to perform read operation using file object
    with gzip.open(tax_file_name, "rb") if tax_file_name.endswith(".gz") else open(
        tax_file_name, "rt"
    ) as infile:
        root_regex_patt: Pattern[str] = re.compile("root")
        count_regex_patt: Pattern[str] = "(\\([0-9]+\\)$)"

        for tp in infile.readline():
            # read data from file object and put it into string.

            # print the data of the string
            # print(tp)

            taxstring: List[str] = tp.split("\t")
            if len(taxstring) < 9:
                continue

            match_results = re.search(taxstring[8], root_regex_patt)
            if match_results:
                taxon: str = re.sub(taxstring[8], count_regex_patt, "")
                orf_tax_map.append((taxstring, taxon))

    return orf_tax_map


def read_pwy_file(pwy_filename: str) -> Dict[str, List[str]]:
    """
    Creates a map from pwy name to a vector of ORF names.

    :param pwy_filename: File name with the pathway and corresponding ORFs in each pathway in the sample.
    """

    pwy_orf: Dict[str, List[str]] = dict()

    k: int = 0

    # open a file to perform read operation using file object
    with gzip.open(pwy_filename, "rb") if pwy_filename.endswith(".gz") else open(
        pwy_filename, "rt"
    ) as infile:
        for tp in infile.readline():
            # read data from file object and put it into string.

            # print the data of the string
            # print(tp)

            orfstring: List[str] = tp.split("\t")
            if len(taxstring) >= 2:
                pwy: str = taxstring[0]
                pwy_orflist: Tuple[str, List[str]] = (
                    pwy,
                    [orf for orf in orfstring[1:]],
                )
                pwy_orf.append(pwy_orflist)
    print(f"No of pwys inserted {k}")
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


def create_tree(taxa: List[str], present_absent: bool) -> NODE:
    """
    Creates the tree with the taxons.

    :param taxa: List of subtaxa that defines a lineage when read from left-to-right.
    :param present_absent: Present or absent boolean flag, i.e., the value of the taxa are either 0 or 1.
    :return: The root node for the tree created.
    """

    root_regex_patt: Pattern[str] = re.compile("root")
    count_regex_patt: Pattern[str] = re.compile("(\\([0-9]+\\)$)")

    tree: NODE = None

    for i in range(0, len(taxa)):
        match_results = re.search(taxa[i], root_regex_patt)
        if match_results:
            taxon: str = re.sub(taxa[i], count_regex_patt, "")
            insert_node(tree, taxon.trim(" ").split(";"), 0, present_absent)

    return tree


def Delta(taxa: List[str], delta_type: DeltaType, wtd_based: bool) -> float:
    """
    Computes the  Delta for the taxa with or without WTD consideration.

    :param taxa: List of the taxa name.
    :param delta_type: Type of the delta.
    :param use_wtd: Flag to indicate whether adjustments of the edge weights are adjusted to WTD.
    :return: A double representing the computed delta.
    """

    present_absent: bool = True

    switch = {DELTA: False, DELTA_STAR: False, DELTA_PLUT: True}

    present_absent = switch.get(delta_type, "None")
    root_node: NODE = create_tree(taxa, present_absent)

    sum_tree_count(root_node)
    subtree_count(root_node)

    delta: float = compute_delta(
        root_node, delta_type, wtd_based
    )  # use_wtd = wtd_based
    delete_tree(root_node)

    return delta


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

    root_regex_patt: Pattern[str] = re.compile("root")
    count_regex_patt: Pattern[str] = re.compile("(\\([0-9]+\\)$)")

    i: int = 0
    tree: NODE = None

    if pwy in pwy_orf:
        for orfid in pwy_orf[pwy]:
            if orfid in orf_taxon:
                match_results = re.search(orf_taxon[orfid], root_regex_patt)
                if match_results > 0:
                    taxon: str = re.sub(orf_taxon[ordid], count_regex_patt, "")
                    insert_node(tree, taxon.trim(" ").split(";"), 0, present_absent)
                    i += 1
    return tree


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


def compute_pwy_tax_distance(tax_file_name: str, pwy_orf_filename: str) -> None:
    """
    Computes the taxonomic distance for the data in the taxonomy file for the give pwy-to-ORF map file.

    :param tax_filename: The input file with the list of taxonomic annotations in the sample.
    :param pwy_orf_filename: The input file with the list of ORFs for individual pathways in the sample.
    """

    # read the file
    orf_taxon: Dict[str, str] = create_orf_tax_map(tax_file_name)

    pwy_orfs: Dict[str, List[str]] = read_pwy_file(pwy_orf_filename)

    taxa: List[str] = []

    print("PWY\tN\tDelta\tDelta*\tDelta+\tDelta (WTD)\tDelta* (WTD)\tDelta+ (WTD)")
    for pwy, orflist in pwy_orfs:
        print(f"{pwy}\t{len(orflist)}")
        root_node: NODE = None

        # delta
        root_node = create_pwy_tree(
            orf_taxon, pwy_orfs, pwy, false
        )  # present_absent =  false)
        sum_tree_count(root_node)
        subtree_count(root_node)

        print(
            f"\t {compute_delta(root_node, DELTA, false)}"
        )  # compute_delta(root_node, DELTA, /* use_wtd = */ false)

        taxa.clear()
        for orfid in orflist:
            if orfid in orf_taxon:
                taxa.append(orf_taxon[orfid])

        print(
            f" ({Delta(taxa, DELTA, false)})"
        )  # Delta(taxa, DELTA, /*use_wtd = */false);
        delete_tree(root_node)

        # delta*
        root_node = create_pwy_tree(orf_taxon, pwy_orfs, pwy, false)
        #  present_absent =  false
        sum_tree_count(root_node)
        subtree_count(root_node)
        print(f"\t{compute_delta(root_node, DELTA_STAR, false)}")
        print(
            f"({Delta(taxa, DELTA_STAR, false)})"
        )  # Delta(taxa, DELTA_STAR, /*use_wtd = */false)
        delete_tree(root_node)

        # delta+
        root_node = create_pwy_tree(orf_taxon, pwy_orfs, pwy, true)
        # present_absent = true
        sum_tree_count(root_node)
        subtree_count(root_node)
        print(f"\t{compute_delta(root_node, DELTA_PLUS, false)}")
        print(
            f"({Delta(taxa, DELTA_PLUS, false)})"
        )  # Delta(taxa, DELTA_PLUS, /*use_wtd = */false)
        delete_tree(root_node)
