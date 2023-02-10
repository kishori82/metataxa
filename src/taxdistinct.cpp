#include "utilities.h"
#include <cinttypes>
#include "taxdistinct.h"


using namespace std;

/**
 * @brief insert_node   Inserts a node in to the tree based on the name of the taxon 
 *    to create the taxonomy tree
 * @param node a NODE ** to start creating the taxonomy tree because it is implemented
 *    in a recursive manner   
 * @param tax_hierarchy The vetor of strings that encodes the taxonomic hierachy, e.g., 
 *   "<taxA, taxB, taxC, taxD> for taxA;taxB;taxC;taxD
 * @param depth         The depth of the taxa hierarchy we can insert into the tree.
 */
void insert_node(
    NODE **node, 
    const vector<std::string> &tax_hierarchy, 
    unsigned int depth,
    bool present_absent
) {
  NODE *n;

  if (*node == nullptr) {
    *node  = new NODE;
    (*node)->name = tax_hierarchy[depth];
    (*node)->depth = depth;
  }
  if (tax_hierarchy.size() == 1 ) {  
    if (present_absent) {
      (*node)->count = 1;
    } else {
      (*node)->count++;
    }
  }

  if (tax_hierarchy.size() > depth + 1 ) {
     depth = depth + 1;
     if ((*node)->children.find(tax_hierarchy[depth]) == (*node)->children.end() ) {
       n = new NODE;
       n->name = tax_hierarchy[depth];
       n->depth = depth;
       (*node)->children.insert( {tax_hierarchy[depth], n }) ;
         //   std::cout <<  "\t" << n->name << " : " << depth << std::endl;
     } else {
       n = (*node)->children.at(tax_hierarchy[depth]);
     }

     if (tax_hierarchy.size() == depth + 1) {
       if (present_absent) {
         n->count = 1;
       } else {
         n->count++;
       }
     }

     //std::cout << "\t" << depth << std::endl;
     insert_node(&n, tax_hierarchy, depth, present_absent) ;
  }
}

/**
 * @brief _sum_tree_count computes the sum of all the taxon nodes in the 
 *
 * @param node pointer to a NODE
 * @param sum pointer to the sum variable
 */
void _sum_tree_count_(NODE *node, int &sum) {
   sum = sum + node->count;
   for (auto it = node->children.begin(); it != node->children.end(); it++) {
      _sum_tree_count_(it->second, sum);
   }
}

/**
 * @brief sum_tree_count computes the sum of all the taxon nodes in the 
 *   tree
 *
 * @param node pointer to a tree node
 * @return the sum total of the node counts
 */
int sum_tree_count(NODE *node) {
   int sum = 0;
   // put the value of the sum in the sum variable   
   _sum_tree_count_(node, sum);
   return sum;
}

/**
 * @brief subtree_count updates the subtree_count variable for each 
 *   node in the tree
 *
 * @param node any NODE * in the tree
 * @return sub-tree count at NODE * node
 */
int subtree_count(NODE *node) {
   int sum = node->count;

   for (auto it = node->children.begin(); it != node->children.end(); it++) {
     sum = sum + subtree_count(it->second);
   }
   node->subtree_count = sum;

   return sum;
}

/**
 * @brief _delete_tree Deletes the nodes of the tree recursively
 *
 * @param node pointer to a NODE
 */
void _delete_tree(NODE *node) {
   for (auto it = node->children.begin(); it != node->children.end(); it++) {
      _delete_tree(it->second);
   }

   for (auto it = node->children.begin(); it != node->children.end(); it++) {
      delete it->second;
   }
   node->children.clear();
}

/**
 * @brief delete_tree deletes the tree of nodes 
 *
 * @param node pointer to a NODE
 */
void delete_tree(NODE *node) {
   _delete_tree(node);
   delete node;
}


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

/** @copydoc compute_delta_star */
double compute_delta(NODE *node, DeltaType deltatype, bool use_wtd) {
   int tot_count = node->subtree_count;
   double prod_count = _count_edge_crossings_(node, tot_count, use_wtd);

   vector<int> count_vector;

   _get_count_vector(node, count_vector);

   // sum-product for all distinct counts for the taxons
   double delta = 0; 
   if (deltatype == DELTA_STAR ) {
     int sum_product = 0;
     for (uint32_t i = 0;  i < count_vector.size(); i++) {
       for (uint32_t j = i+1;  j < count_vector.size(); j++) {
         sum_product +=  count_vector[i] * count_vector[j];
       }
     }
     delta = static_cast<double>(prod_count)/static_cast<double>(sum_product);
   } else {
     int N = std::accumulate(count_vector.begin(), count_vector.end(), 0);
     double nums =static_cast<double>(N* (N - 1)/2);
     delta = static_cast<double>(prod_count)/nums;
   }

   //std::cout << "prod_count " << prod_count << "  " << sum_product << std::endl;
   return delta;
} 

/**
 * @brief _print_tree_ The node that prints taxon and node count for the tree 
 *
 * @param node    The NODE* for the node of the tree
 * @param taxons  Reference to the vector of the taxons 
*/ 
void _print_tree_(NODE *node, std::vector<std::string> &taxons) {
   taxons.push_back(node->name);
   if (node->count >= 0) {
   //   std::cout << "\t<" << taxon_concat(taxons) << ">\t" << node->count << std::endl;
     // std::cout <<  taxon <<  std::endl;
   }

   for (auto it = node->children.begin(); it != node->children.end(); it++) {
      _print_tree_(it->second, taxons);
   }
   taxons.pop_back();
}

/**
 * @brief print_tree Prints the taxons and the counts in the tree
 * 
 * @param node The root NODE* for the tree
*/
void print_tree(NODE *node) {
   std::vector<std::string> taxons;

   _print_tree_(node, taxons);
}

/** @copydoc read_tax_file */
NODE *read_tax_file(const string &tax_file_name) {
  fstream newfile;

  std::cout << "TODO: Make taxa names case insensitive" << std::endl;
  regex root_regex_patt("root");
  regex count_regex_patt("(\\([0-9]+\\)$)");
  smatch sm;
  NODE *tree = nullptr;

  vector<std::string> taxstring;
  newfile.open(tax_file_name, ios::in); //open a file to perform read operation using file object
  int k = 0;
  if (newfile.is_open()){   //checking whether the file is open
    string tp;
    while (getline(newfile, tp)){ //read data from file object and put it into string.
    //  cout << tp << "\n"; //print the data of the string
      taxstring = split(tp, '\t');
      if (taxstring.size() >= 10) {
        regex_search(taxstring[8], sm,  root_regex_patt);
        if (sm.size() > 0) { 
          std::string taxon = regex_replace(taxstring[8], count_regex_patt, "");
          if (sm.size() >= 0) {
            /*
            for(auto x: split(trim(taxon, " "), ';') ) {
            std::cout << '\t' << "<" <<  x << ">" << std::endl;
            }
            */
            // std::cout << taxon << std::endl;
            k++;
            insert_node(&tree, split(trim(taxon, " "),';'), 0, false);
          }
      }
     }
    }
    newfile.close(); //close the file object.
  }
  std::cout << "No of nodes inserted " << k << std::endl;
  return tree;
} 

/**
 * @brief read_orf_tax_map  Create a orf to taxon map from the functional 
 *   and taxonomic file
 *
 * @return Pointer to the string to (pointer to) a vector of strings 
 *
*/
STRING_STRING_MAP *create_orf_tax_map(const std::string &tax_file_name) {
  STRING_STRING_MAP *orf_tax_map = new STRING_STRING_MAP;

  fstream file;
  file.open(tax_file_name, ios::in); //open a file to perform read operation using file object

  if (!file.is_open()) {   //checking whether the file is open
     std::cerr << "Failed to open file " << tax_file_name << std::endl; 
     abort();
  }

  regex root_regex_patt("root");
  regex count_regex_patt("(\\([0-9]+\\)$)");
  smatch sm;

  vector<std::string> taxstring;
  string tp;
  while (getline(file, tp)) { //read data from file object and put it into string.
    //  cout << tp << "\n"; //print the data of the string
    taxstring = split(tp, '\t');
    if (taxstring.size() < 9) continue;
    regex_search(taxstring[8], sm,  root_regex_patt);
    if (sm.size() > 0) { 
       std::string taxon = regex_replace(taxstring[8], count_regex_patt, "");
       std::pair<std::string, std::string> tuple(taxstring[0], taxon);
       orf_tax_map->insert(tuple);
    }
  }

  return orf_tax_map;
};

/**
 * @brief read_pwy_file Creates a map from pwy name to a vector of orf names
 *
 * @return Pointer to an unordered map, keys are pathway name and values are 
     pointer to the vector of strings 
*/
//TODO: deallocate the nodes
STRING_STRING_VECTOR_MAP *read_pwy_file(const string &pwy_filename) {
  fstream newfile;
  std::vector<std::string> taxstring;
  newfile.open(pwy_filename, ios::in); //open a file to perform read operation using file object

  STRING_STRING_VECTOR_MAP *pwy_orf = new STRING_STRING_VECTOR_MAP;

  int k = 0;
  if (newfile.is_open()) {   //checking whether the file is open
    string tp;
    while (getline(newfile, tp)) { //read data from file object and put it into string.
    //  cout << tp << "\n"; //print the data of the string
      taxstring = split(tp, '\t');
      if (taxstring.size() >= 2) {
         std::string pwy = taxstring[0];
         std::pair<std::string, std::vector<std::string> *> tuple(pwy, new std::vector<std::string>);
         pwy_orf->insert(tuple);

         auto it = taxstring.begin();
         std::advance(it, 1);
         for( ; it != taxstring.end(); it++) {
           pwy_orf->at(pwy)->push_back(*it);
         }
      }
    }
    newfile.close(); //close the file object.
  }
  std::cout << "No of pwys inserted " << k << std::endl;
  return pwy_orf;
} 

/**
 * @brief delete_pwy_orfs delete the orf to pointer to vector map
 *
 * @param pwy_orf Pointer to STRING_STRING_VECTOR_MAP 
*/
void delete_pwy_orfs(STRING_STRING_VECTOR_MAP *pwy_orf) {
   // delete the vectors allocated
   for (auto it = pwy_orf->begin(); it != pwy_orf->end(); it++) {
      delete it->second;
   }
   delete pwy_orf;
}

/** @copydoc create_pwy_tree */
NODE *create_tree(
    const STRING_VECTOR &taxa,
    bool present_absent
) 
{
  regex root_regex_patt("root");
  regex count_regex_patt("(\\([0-9]+\\)$)");
  smatch sm;

  NODE *tree = nullptr;

  for (unsigned int i = 0; i < taxa.size(); i++) {
    regex_search(taxa[i], sm,  root_regex_patt);
    if (sm.size() > 0) { 
      std::string taxon = regex_replace(taxa[i], count_regex_patt, "");
      insert_node(&tree, split(trim(taxon, " "),';'), 0, present_absent);
    }
  }
  return tree;
}
/** @copydoc Delta */
double Delta(const STRING_VECTOR &taxa,
             DeltaType deltatype,
             bool wtd_based) {

  bool present_absent = true; 

  switch(deltatype) {
    case DELTA:
    case DELTA_STAR:
      present_absent = false; 
      break;
    case DELTA_PLUS:
      present_absent = true; 
      break;
    default:
      break;
  }  

  NODE *root_node = create_tree(taxa, present_absent);
  sum_tree_count(root_node);
  subtree_count(root_node); 
  double delta  = compute_delta(root_node, deltatype, /* use_wtd = */ wtd_based);
  delete_tree(root_node);
  return delta;
}

/** @copydoc create_pwy_tree */
NODE *create_pwy_tree(
    const STRING_STRING_MAP *orf_taxon, 
    const STRING_STRING_VECTOR_MAP *pwy_orf, 
    const std::string pwy,
    bool present_absent
) {
  regex root_regex_patt("root");
  regex count_regex_patt("(\\([0-9]+\\)$)");
  smatch sm;

  int i = 0;
  NODE *tree = nullptr;
  if (pwy_orf->find(pwy) != pwy_orf->end()) {
    for (auto it = pwy_orf->at(pwy)->begin(); it != pwy_orf->at(pwy)->end(); it++) {
      if (orf_taxon->find(*it) != orf_taxon->end()) {
        regex_search(orf_taxon->at(*it), sm,  root_regex_patt);
        if (sm.size() > 0) { 
          std::string taxon = regex_replace(orf_taxon->at(*it), count_regex_patt, "");
          insert_node(&tree, split(trim(taxon, " "),';'), 0, present_absent);
          i++;
        }
      }
    }
  }
  return tree;
}

/** @copydoc compute_pwy_tax_distance */
void compute_tax_distance(const std::string &tax_file_name) {

  // read the file
  NODE *tree = read_tax_file(tax_file_name) ;

  // count the subtree and update in each  node
  sum_tree_count(tree);
  std::cout << "Subtree count at root: " << subtree_count(tree) << endl; 

  std::cout << "Delta* " << compute_delta(tree, DELTA_STAR, false) << endl; 
}

/** @copydoc compute_pwy_tax_distance */
void compute_pwy_tax_distance(const std::string &tax_file_name, \
  const std::string &pwy_orf_filename) {

  // read the file
  STRING_STRING_MAP *orf_taxon = create_orf_tax_map(tax_file_name);

  STRING_STRING_VECTOR_MAP *pwy_orfs = read_pwy_file(pwy_orf_filename);

  STRING_VECTOR taxa;
 
  std::cout << std::setw(10);
  std::cout << "PWY\tN\tDelta\tDelta*\tDelta+\tDelta (WTD)\tDelta* (WTD)\tDelta+ (WTD)" << std::endl; 
  for (auto it = pwy_orfs->begin(); it != pwy_orfs->end(); it++) {
    std::cout << it->first; 
    std::cout << "\t" << it->second->size(); 
    NODE *root_node = nullptr;

    // delta
    root_node = create_pwy_tree(orf_taxon, pwy_orfs, it->first, /* present_absent = */ false);
    sum_tree_count(root_node);
    subtree_count(root_node); 
    std::cout << "\t" << compute_delta(root_node, DELTA, /* use_wtd = */ false);  

    taxa.clear();
    for (auto orfid_it = it->second->begin(); orfid_it != it->second->end(); orfid_it++) {
       if (orf_taxon->find(*orfid_it) != orf_taxon->end()) {
           taxa.push_back(orf_taxon->at(*orfid_it));
       }
    }

    std::cout << " (" << Delta(taxa, DELTA, /*use_wtd = */false) << ")";  
    delete_tree(root_node);

    // delta*
    root_node = create_pwy_tree(orf_taxon, pwy_orfs, it->first, /* present_absent = */ false);
    sum_tree_count(root_node);
    subtree_count(root_node); 
    std::cout << "\t" << compute_delta(root_node, DELTA_STAR, false);  
    std::cout << " (" << Delta(taxa, DELTA_STAR, /*use_wtd = */false) << ")";  
    delete_tree(root_node);

    // delta+
    root_node = create_pwy_tree(orf_taxon, pwy_orfs, it->first, /* present_absent = */ true);
    sum_tree_count(root_node);
    subtree_count(root_node); 
    std::cout << "\t" << compute_delta(root_node, DELTA_PLUS, false);
    std::cout << " (" << Delta(taxa, DELTA_PLUS, /*use_wtd = */false) << ")";  
    delete_tree(root_node);

    // delta WTD
    root_node = create_pwy_tree(orf_taxon, pwy_orfs, it->first, /* present_absent = */ false);
    sum_tree_count(root_node);
    subtree_count(root_node); 
    std::cout << "\t" << compute_delta(root_node, DELTA, true);  
    delete_tree(root_node);

    // delta* WTD
    root_node = create_pwy_tree(orf_taxon, pwy_orfs, it->first, /* present_absent = */ false);
    sum_tree_count(root_node);
    subtree_count(root_node); 
    std::cout << "\t" << compute_delta(root_node, DELTA_STAR, true);  
    delete_tree(root_node);

    // delta+ WTD
    root_node = create_pwy_tree(orf_taxon, pwy_orfs, it->first, /* present_absent = */ true);
    sum_tree_count(root_node);
    subtree_count(root_node); 
    std::cout << "\t" << compute_delta(root_node, DELTA_PLUS, /* use_wtd = */ true);
    delete_tree(root_node);
    std::cout << std::endl; 
  }

  // clear orf to taxon map 
  orf_taxon->clear();
  delete orf_taxon;

  // clear pwy to orfs maps
  delete_pwy_orfs(pwy_orfs);

}

