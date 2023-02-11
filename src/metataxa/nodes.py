class NODE:
    # map of taxon name to NODE  for the tree
    children: dict[str, : NODE] = None;   

    # number of counts for the taxa
    count: int = 0

    # sum total of the counts in the subtree at the node
    int subtree_count: int = 0

    # depth in the tree
    depth: int = 0
    # initialize the counts to 0
   
    # name of the taxon
    name: str = None
