# Copyright 2014 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

import unittest
from Bio.Ontology.Graph import DiGraph

class GraphTest(unittest.TestCase):


    def test__get_reachable(self):
        
        g = DiGraph([(1,2), (2,3), (3,4), (3,5), (5,2), (5,6), (6,8), (6,7), (2,9), (9,2)])
        reachable_from_2 = g._get_reachable(g.get_node(2))
        expected = ([], set([2, 3, 4, 5, 6, 7, 8, 9]))
        self.assertEqual(expected, reachable_from_2)
        reachable_from_6 = g._get_reachable(g.get_node(6))
        expected = ([], set([8, 7]))
        self.assertEqual(expected, reachable_from_6)
        # Found cycles
        expected =  [[2, 9], [2, 5, 3]]
        self.assertEqual(expected, g.cycles)
        

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)