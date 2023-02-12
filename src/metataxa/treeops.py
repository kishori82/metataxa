from nodes import *
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
    if depth == 0:
        node.name = tax_hierarchy[depth]
        node.depth = depth

    if len(tax_hierarchy) == 1:
        if present_absent:
            node.count = 0
        node.count += 1

    if len(tax_hierarchy) > depth + 1:
        depth = depth + 1

        if tax_hierarchy[depth] not in node.children:
            n: NODE = NODE()
            n.name = tax_hierarchy[depth]
            n.depth = depth
            node.children[tax_hierarchy[depth]] = n
        else:
            n = node.children[tax_hierarchy[depth]]

        if len(tax_hierarchy) == depth + 1:
            if present_absent:
                n.count = 1
            else:
                n.count += 1

        insert_node(n, tax_hierarchy, depth, present_absent)

def create_tree(taxa: List[str], present_absent: bool) -> NODE:
    """
    Creates the tree with the taxons.

    :param taxa: List of subtaxa that defines a lineage when read from left-to-right.
    :param present_absent: Present or absent boolean flag, i.e., the value of the taxa are either 0 or 1.
    :return: The root node for the tree created.
    """

    root_regex_patt: Pattern[str] = re.compile(r"root")
    taxid_regex_patt: Pattern[str] = re.compile(r"([0-9]+)$")

    tree: NODE = NODE()

    for i in range(0, len(taxa)):
        match_results = root_regex_patt.search(taxa[i])
        if match_results:
            taxon: str = re.sub(taxid_regex_patt, "", taxa[i])
            insert_node(
                tree,
                [x for x in taxon.strip(" ").split(";") if x.strip()],
                0,
                present_absent,
            )

    return tree

def save_tree_count(node: NODE) -> None:
    """
    Saves the values of the count variable in the NODEs of the subtree rooted in
    node to the actual_count variable.

    :param node: A node in the tree taxonomic tree.
    """
    node.actual_count = node.count
    for child_node in node.children.values():
        save_tree_count(child_node)


def __delete_tree(node: NODE) -> None:
    """
    Recursively deletes the nodes in the sub-tree rooted in node ``node".

    :param node: Root NODE of the sub-tree to start deleting from.
    """

    for child_node in node.children.values():
        __delete_tree(child_node)

    node.children.clear()


def delete_tree(node: NODE) -> None:
    """
    Deletes the nodes in the sub-tree rooted in node.

    :param node: Root NODE of the tree .
    """

    __delete_tree(node)
    node = None
