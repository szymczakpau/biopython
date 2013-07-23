# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

import unittest, string
from Bio.Ontology.IO.GraphIO import GmlWriter
from StringIO import StringIO
from Bio.Ontology.Graph import DiGraph

class GraphIOTest(unittest.TestCase):

    
    def test_data_to_gml(self):
        correct_output = """  l1 [
    l2a 2
    l2b [
      l3b [
        l4 4
      ]
      l3c "[1, 2, 3]"
      l3a 3
    ]
  ]"""
        writer = GmlWriter(None)
        lines = writer.data_to_gml("l1", {"l2a" : 2, "l2b" : { "l3a" : 3, "l3b" : {"l4" : 4}, "l3c" : [1, 2, 3]}}, 1)
        output = string.join(lines, "\n")
        self.assertEqual(correct_output, output)
    
    def test_write_file(self):
        correct_output = """graph [
  directed 1
  node [
    id 0
    label "1"
    a 1
    abb "[1, 2]"
    bb 2
  ]
  node [
    id 1
    label "2"
  ]
  node [
    id 2
    label "3"
  ]
  node [
    id 3
    label "4"
  ]
  node [
    id 4
    label "5"
  ]
  edge [
    source 0
    target 1
  ]
  edge [
    source 0
    target 3
    x "x"
  ]
  edge [
    source 1
    target 2
  ]
  edge [
    source 2
    target 3
  ]
  edge [
    source 2
    target 4
  ]
  edge [
    source 4
    target 1
    label "zzzz"
  ]
]"""
        out = StringIO()
        writer = GmlWriter(out)
        graph = DiGraph([(1,2), (2,3), (3,4), (3,5)])
        graph.update_node(1, {'a' : 1, 'bb' : 2, 'abb' : [1, 2] })
        graph.add_edge(1, 4, {'x' : 'x'})
        graph.add_edge(5, 2, "zzzz")
        writer.write_file(graph)
        self.assertEqual(correct_output, out.getvalue())

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)