"""
A binary search tree for peak classification.
"""

class SearchNode:
    def __init__(self, value):
        self.left = None
        self.right = None
        self.value = value

    def insert_node(self, value):
        """
        Add a node into the search tree.

        Parameters
        ----------
        value: float

        Returns
        -------
        None
        """

        if self.value == None:
            self.value = value
        else:
            if self.value > value:
                if self.left == None:
                    self.left = SearchNode(value)
                else:
                    self.left.insert_node(value)
            else:
                if self.right == None:
                    self.right = SearchNode(value)
                else:
                    self.right.insert_node(value)

