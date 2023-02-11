from collections.abc import Mapping
from nodes import NODE
from typing import Pattern, Tuple, List
import re

def insert_node( node: NODE , tax_hierarchy list[str], depth: int = 0, present_absent: bool = False)->

    """
    Inserts a node in to the tree based on the name of the taxon 
    to create the taxonomy tree

    Parameters:
    node (NODE): A NODE to start creating the taxonomy tree because it is implemented in a recursive manner   
    tax_hierarchy (str): The vetor of strings that encodes the taxonomic hierachy, e.g., "<taxA, taxB, taxC, taxD> for taxA;taxB;taxC;taxD
    depth (int): The depth of the taxa hierarchy we can insert into the tree.

    Returns: 
    None
    """
  n: NODE = None

  if node == None: 
    node  = NODE();
    node.name = tax_hierarchy[depth]
    node.depth = depth

  if len((tax_hierarchy) == 1:  
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
       n = node.children[tax_hierarchy[depth]];

     if len(tax_hierarchy) == depth + 1:
       if present_absent:
         n.count = 1
       else:
         n->count += 1

     # print(f"\t{depth}")
     insert_node(n, tax_hierarchy, depth, present_absent)

def _sum_tree_count_(node: NODE, sum:int)-> None:
   """
   Computes the sum of all the taxon nodes in the 
  
   Parameters:
   node [NODE]: A node in the tree
   sum: The sum variable

   Returns:
   None
   """
   sum = sum + node.count;
   for taxon, child_node in node.children.items():
      _sum_tree_count_(child_node, sum);


def sum_tree_count(node: NODE) -> int:
   """
   Computes the sum of all the taxon nodes in the tree

   Parameters:
   param [NODE]:  node in the tree

   Returns:
     the sum total of the node counts
   """

   sum: int = 0
   #  put the value of the sum in the sum variable   
   _sum_tree_count_(node, sum)
   return sum

def subtree_count(node: NODE)->int:
   """
   Updates the subtree_count variable for each NODE in the tree

   Parameters:
   node [NODE]: Any NODE in the tree

   Returns:
   Count for subtree rooted at NODE node
   """

   _sum: int = node.count

   for taxon, child_node in node.children.items():
     _sum = _sum + subtree_count(child_node)

   node.subtree_count = _sum

   return _sum

def ___delete_tree(node: NODE)->None:
   """
   Recursively deletes the nodes in the sub-tree rooted in node

   Parameters:
   node [NODE]: Root NODE of the sub-tree 

   Returns:
   None
   """

   for child_node in node.children:
      ___delete_tree(child_node);

   node.children.clear();

def delete_tree(node: NODE)->None: 
   """
   Deletes the nodes in the sub-tree rooted in node

   Parameters:
   node [NODE]: Root NODE of the tree 

   Returns:
   None
   """

   __delete_tree(node)
   node = None;


/**
 * @brief _count_edge_crossings The sum of the number of pairs of taxons in the tree 
 *   that crosses the edge across all edges
 *
 * @param node         The root node
 * @param total_count  The total count of the tree at that node
 * @return total sum-product across all edges
 */
double _count_edge_crossings_(NODE *node, const int tot_count, bool use_wtd) {


   double sum = 0;
   for (auto it = node->children.begin(); it != node->children.end(); it++) {
     if (use_wtd) {
       sum = sum + (tot_count - it->second->subtree_count) * it->second->subtree_count * 
         pow(0.5, it->second->depth);
     } else {
       sum = sum + (tot_count - it->second->subtree_count) * it->second->subtree_count;
     }
     sum = sum + _count_edge_crossings_(it->second, tot_count, use_wtd);
   }
   return sum;
}

/**
 * @brief _get_count_vector  This function returns the vector of counts for all 
 *   nodes in the tree
 *
 * @param node          The pointer to the node, which is the root of the tree
 * @param count_vector  The vector of the taxonomic counts for the nodes
 */
void _get_count_vector(NODE *node, vector<int> &count_vector) {

   if (node->count > 0) {
      count_vector.push_back(node->count);
   }

   for (auto it = node->children.begin(); it != node->children.end(); it++) {
      _get_count_vector(it->second, count_vector);
   }
}

def compute_delta(node: NODE, deltatype: DeltaType, use_wtd: bool):
   """
   Compute delta star

   Parameters:
   node [NODE]: Node for the subtree to compute delta star
   deltatype [DeltaType]: Type of the delta 
   use_wtd [bool]: Flag to indicate whether it is WTD

   Returns:
    A double representing the computed delta
   """

   tot_count: int = node.subtree_count
   prod_count: double = _count_edge_crossings_(node, tot_count, use_wtd)

   count_vector: list[int]=[]

   _get_count_vector(node, count_vector)

   # sum-product for all distinct counts for the taxons
   delta: float = 0; 
   if deltatype == DELTA_STAR:
      sum_product: int = 0;
      for i in range(0, len(count_vector)): 
         for j in range(i+1,  len(count_vector)):
            sum_product +=  count_vector[i] * count_vector[j]

      delta = float(prod_count)/float(sum_product)
   else:
     N: int = sum(count_vector)
     nums: float = float(N)*(float(N) - 1.0)/2.0
     delta = float(prod_count)/nums

   # print(f"prod_count {prod_count}  {sum_product}")
   return delta;

def _print_tree_(node: NODE, taxons: list[str]): 
   """* 
   The node that prints taxon and node count for the tree 
  
   Parameters:
   node [NODE]: The NODE for the node of the tree
   taxons list[str]:  Reference to the vector of the taxons 

   Returns:
   None
   """

   taxons.append(node.name)

   if node.count >= 0:
     # print(f"\t<{taxon_concat(taxons)}>\t{node.count}\t{taxon}")
     pass

   for child_node in node.children.values():
      _print_tree_(child_node, taxons)

   taxons.pop()


def print_tree(node: NODE)-> None:
   """
   Prints the taxons and the counts in the tree

   Parameters:
   node [NODE]:  The root NODE for the tree
   """
   taxons: list[str] = []
   _print_tree_(node, taxons)
 

def read_tax_file(tax_file_name: str)-> NODE:
  """
  Reads the taxa file 

  Parameters:
  tax_file_name [str]: Taxa file name

  Returns:
  The root NODE for the taxa in the file 
  """

  print(f"TODO: Make taxa names case insensitive")
  root_regex_patt: Pattern[str] = re.compile("root")
  count_regex_patt: Pattern[str]  = re.compile("(\\([0-9]+\\)$)")

  tree: NODE = None
  taxstring: list[str] = []

  with open(tax_file_name, 'rb') if newfile.endswith('.gz') else open(tax_file_name, 'r') as newfile: 
  # open a file to perform read operation using file object

      print("TODO: make the taxstring columns and count flexible")
      int k = 0;
      # read data from file object and put it into string.
      for tp: str in newfile.readline():
          # print the data of the string
          # print(tp)
          taxstring: list[str]  = tp.split('\t')
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
                  insert_node(tree, taxon.trim(" ").split(';'), 0, false)
  print(f"No of nodes inserted {k}")
  return tree;


def create_orf_tax_map(tax_file_name: str) -> dict[str, str]:
   """
   Create an ORF to taxon map from the functional and taxonomic file
 
   Parameters:
   tax_file_name [str]: Taxonomy file name

   Returns:
   dict[str, str], an ORF to taxon map
   """

   orf_tax_map: dict[str, str]  = dict[str, str]
 
   # open a file to perform read operation using file object
   with gzip.open(tax_file_name, 'rb') if tax_file_name.endswith('.gz') else open(tax_file_name, 'rt') as file: 
       root_regex_patt: Pattern[str] = re.compile("root")
       count_regex_patt: Pattern[str] = ("(\\([0-9]+\\)$)")
 
       for tp: str in file.readline():
           # read data from file object and put it into string.
 
           # print the data of the string
           # print(tp)
 
         taxstring: list[str] = tp.split('\t')
         if len(taxstring) < 9:
            continue
 
         match_results = re.search(taxstring[8], root_regex_patt)
         if match_results: 
            std::string taxon: [str] = re.sub(taxstring[8], count_regex_patt, "")
            orf_tax_map.append((taxstring, taxon))

  return orf_tax_map;

def read_pwy_file(const string &pwy_filename)-> dict[str, list[str]]:
    """
    Creates a map from pwy name to a vector of ORF names
    
    Parameters:
    pwy_filename [str]: Filename with the pathway and corresponding ORFs in each pathway in the sample
    """

    pwy_orf: dict[str, list[str]= dict()
  
    k: int = 0
  
    # open a file to perform read operation using file object
    with gzip.open(pwy_filename, 'rb') if pwy_filename.endswith('.gz') else open(pwy_filename, 'rt') as infile: 
        for tp: str in file.readline():
           # read data from file object and put it into string.
   
           # print the data of the string
           # print(tp)
   
           orfstring: list[str] = tp.split('\t')
           if len(taxstring) >= 2:
              pwy: str = taxstring[0]
              pwy_orflist: Tuple[str, list[str]] = (pwy, [orf for orf in orfstring[1:]]:List[str][])
              pwy_orf.append(pwy_orflist)
    print(f"No of pwys inserted {k}")
    return pwy_orf;
             

def delete_pwy_orfs(pwy_orf: dict[str, List[str]]) {
    """
    Deletes the pathway orfs
    
    Parameters:
    pwy_orf dict[str, List[str]]: Map of pathways to the list of ORFs.
    """

    # delete the vectors allocated
    for pwy in pwy_orf:
        pwy_orf[pwy].clear()
    pwy_orf.clear()


def create_tree(taxa: list[str], present_absent: bool): 
   """
   Create the tree with the taxons @copydoc create_pwy_tree */

   Parameters:
   taxa list[str]: List of subtaxa that defines a lineage 
   present_absent [bool]: Present or absent boolean flag

   Returns:
   [NODE], The root node for the tree created
   """

   root_regex_patt: Pattern[str] = re.compile("root")
   count_regex_patt: Pattern[str] = re.compile("(\\([0-9]+\\)$)")
 
   tree: NODE = None
 
   for i in range(0, len(taxa)):
     match_results = re.search(taxa[i], root_regex_patt)
     if match_results: 
       taxon: str = re.sub(taxa[i], count_regex_patt, "")
       insert_node(tree, taxon.trim(" ").split(';'), 0, present_absent)
 
   return tree;

def Delta(taxa: list[str], DeltaType deltatype, wtd_based: bool):
  """ 
  Computes the  Delta 

  Parameters:
  taxa list[str]: List of the taxa name
  deltatype [DeltaType]: Type of the delta 
  use_wtd [bool]: Flag to indicate whether it is WTD

  Returns:
  A double representing the computed delta
  """
 
  bool present_absent = True

  switch = {DELTA: False, DELTA_STAR: False, DELTA_PLUT: True} 

  present_absent = switch.get(deltatype, 'None') 
  root_node: NODE = create_tree(taxa, present_absent)

  sum_tree_count(root_node)
  subtree_count(root_node) 

  double delta: float  = compute_delta(root_node, deltatype,  wtd_based) # use_wtd = wtd_based
  delete_tree(root_node)

  return delta

def create_pwy_tree(orf_taxon: dict[str, str], pwy_orf: dict[str, List[str]],  pwy: str, present_absent: bool)-> NODE:
   """
   Create the tree for just one pathway 

   Parameters:
   orf_taxon dict[str, str]: Map of ORF to taxon
   pwy_orf dict[str, List[str]]: List of ORFs for the pathways

   Returns:
   [NODE], The root node for the tree created
   """

   root_regex_patt: Pattern[str] = re.compile("root")
   count_regex_patt: Pattern[str] = re.compile("(\\([0-9]+\\)$)")
 
  i:int  = 0
  tree: NODE = None

  if pwy in pwy_orf:
      for orfid in pwy_orf[pwy]:
         if orfid in orf_taxon:
             match_results = re.search(orf_taxon[orfid], root_regex_patt);
             if match_results > 0:  
                 taxon: str = re.sub(orf_taxon[ordid], count_regex_patt, "")
                 insert_node(tree, taxon.trim(" ").split(';'), 0, present_absent)
                 i += 1
  return tree


def compute_tax_distance(tax_file_name: str) -> None :
  """ 
  Compute the taxonomic distance for the data in the taxonomy file
  
  Parameters:
  tax_file_name [str]: The input file with the list of taxonomic annotations in the sample 

  Returns:
  None
  """

  # read the file
  tree = read_tax_file(tax_file_name) 

  # count the subtree and update in each  node
  sum_tree_count(tree)

  print(f"Subtree count at root: {subtree_count(tree)}") 

  print(f"Delta*  {compute_delta(tree, DELTA_STAR, false)}") 

def compute_pwy_tax_distance(tax_file_name: str, pwy_orf_filename: str)->None:
  """ 
  Compute the taxonomic distance for the data in the taxonomy file for the give pwy-to-ORF map file
  
  Parameters:
  tax_filename [str]: The input file with the list of taxonomic annotations in the sample 
  pwy_orf_filename [str]: The input file with the list of ORFs for individual pathways in the sample 

  Returns:
  None
  """

  # read the file
  orf_taxon: dict[str, str] = create_orf_tax_map(tax_file_name)

  pwy_orfs: dict[str, list[str]] = read_pwy_file(pwy_orf_filename)

  taxa: list[str] = []
 
  print("PWY\tN\tDelta\tDelta*\tDelta+\tDelta (WTD)\tDelta* (WTD)\tDelta+ (WTD)") 
  for pwy, orflist in pwy_orfs:
    print(f"{pwy}\t{len(orflist)}") 
    root_node: NODE = None

    # delta
    root_node = create_pwy_tree(orf_taxon, pwy_orfs, it->first,  false) # present_absent =  false)
    sum_tree_count(root_node)
    subtree_count(root_node) 

    print(f"\t {compute_delta(root_node, DELTA, /* use_wtd = */ false)}")  # compute_delta(root_node, DELTA, /* use_wtd = */ false)  

    taxa.clear()
    for orfid in orflist:
       if orfid in orf_taxon:
           taxa.append(orf_taxon[orfid])

    print(f" ({Delta(taxa, DELTA, false)})"  # Delta(taxa, DELTA, /*use_wtd = */false);  
    delete_tree(root_node);

    # delta*
    root_node = create_pwy_tree(orf_taxon, pwy_orfs, it->first, /* present_absent = */ false);
    sum_tree_count(root_node)
    subtree_count(root_node) 
    print(f"\t{compute_delta(root_node, DELTA_STAR, false)}")  
    print(f"({Delta(taxa, DELTA_STAR, /*use_wtd = */false)})")  # Delta(taxa, DELTA_STAR, /*use_wtd = */false)   
    delete_tree(root_node)

    # delta+
    root_node = create_pwy_tree(orf_taxon, pwy_orfs, it->first, /* present_absent = */ true);
    sum_tree_count(root_node)
    subtree_count(root_node) 
    print(f"\t{compute_delta(root_node, DELTA_PLUS, false)}")
    print(f"({Delta(taxa, DELTA_PLUS, false)})")  # Delta(taxa, DELTA_PLUS, /*use_wtd = */false)   
    delete_tree(root_node)

