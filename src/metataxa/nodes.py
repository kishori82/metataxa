from typing import List


class NODE:
    def __init__(self):
        # map of taxon name to NODE  for the tree
        self.children = dict()

        # number of counts for the taxa, this value may change during simulation
        # confidence interval estimation
        self.count: int = 0

        # actual counts for the taxa
        self.actual_count: int = 0

        # sum total of the counts in the subtree at the node
        self.subtree_count: int = 0

        # depth in the tree
        self.depth: int = 0
        # initialize the counts to 0

        # name of the taxon
        self.name: str = None
