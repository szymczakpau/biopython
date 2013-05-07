# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.


from functools import total_ordering

class DiGraph:
    """
    Base class for directed graph representation.

    Nodes' labels can be any hashable objects.

    Examples
    --------

    >>> g = DiGraph()

    """
    
    def __init__(self, edges=None):
        """
        Initialize graph with edges.
    
        Parameters
        ----------
        data - list of edges in a graph.
        """
        self.nodes = {}
        if edges != None:
            for (u, v) in edges:
                self.add_edge(u, v)

    def add_edge(self, u, v):
        """
        Adds an edge u->v to the graph

        Parameters
        ----------
        u,v - nodes connected by the edge
        """
        if u not in self.nodes:
            self.add_node(u)
        if v not in self.nodes:
            self.add_node(v)
        
        u_node = self.nodes[u]
        v_node = self.nodes[v]
        v_node.pred.add(u_node)
        u_node.succ.add(v_node)
        
    def add_node(self, u, data = None):
        """
        Adds node to the graph

        Parameters
        ----------
        u - node to add
        """
        if u not in self.nodes:
            self.nodes[u] = DiNode(u, data)
        else:
            raise ValueError("Node already exists")

    def get_node(self, u):
        """
        Gets node from the graph

        Parameters
        ----------
        u - node id
        """
        return self.nodes[u]

@total_ordering
class DiNode:
    """
    Class containing information about graph structure. Only used
    internally in DiGraph.

    Nodes with the same label are not distinguishable.
    
    >>> a = DiNode(1)
    >>> b = DiNode(1)
    >>> a >= b
    True
    >>> a
    1
    """

    def __init__(self, label, data = None):
        """
        Initialize node with data.

        Parameters
        ----------
        label - node label
        data - internal node data
        """
        self.label = label
        self.data = data
        self.pred = set()
        self.succ = set()
        self.attr = {}

    def __lt__(self, other):
        return self.label< other.label

    def __eq__(self, other):
        return self.label == other.label

    def __hash__(self):
        return hash(self.label)

    def __str__(self):
        return str(self.label)

    def __repr__(self):
        return repr(self.label)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
