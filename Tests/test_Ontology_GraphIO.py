# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

from Bio._py3k import StringIO

import unittest
from Bio.Ontology.IO.GraphIO import GmlWriter
from Bio.Ontology.Graph import DiGraph

class GraphIOTest(unittest.TestCase):

    
    def test_data_to_gml(self):
        correct_output = """  l1 [
    l2b [
      l3c "[1, 2, 3]"
    ]
  ]"""
        writer = GmlWriter(None)
        lines = writer.data_to_gml("l1", {"l2b" : { "l3c" : [1, 2, 3]}}, 1)
        output = "\n".join(lines)
        self.assertEqual(correct_output, output)
    
    def test_write(self):
        correct_output_a = """graph [
  directed 1
  node [
    id 0
    label "1"
    a 1
  ]
  node [
    id 1
    label "2"
  ]
  edge [
    source 0
    target 1
    x "x"
  ]
  edge [
    source 1
    target 0
    label "zzzz"
  ]
]"""
        correct_output_b =  """graph [
  directed 1
  node [
    id 0
    label "2"
  ]
  node [
    id 1
    label "1"
    a 1
  ]
  edge [
    source 1
    target 0
    x "x"
  ]
  edge [
    source 0
    target 1
    label "zzzz"
  ]
]"""
        out = StringIO()
        writer = GmlWriter(out)
        graph = DiGraph()
        graph.add_node(1, {'a' : 1 })
        graph.add_edge(1, 2, {'x' : 'x'})
        graph.add_edge(2, 1, "zzzz")
        writer.write(graph)
        self.assertIn(out.getvalue(), set([correct_output_a, correct_output_b]))

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
