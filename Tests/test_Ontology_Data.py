# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

import unittest

from Bio.Ontology.IO.OboIO import terms_to_graph

class OntologyGraphTest(unittest.TestCase):

    def setUp(self):
        
        self.terms_pro = [("Term", {"id" : ["GO:0030534"],
                                 "name" : ["adult behavior"],
                                 "is_a" : ["GO:0050896"]}),
                          ("Term", {"id" : ["GO:0045471"],
                                 "name" : ["response to ethanol"],
                                 "is_a" : ["GO:0050896"]}),
                          ("Term", {"id" : ["GO:0050896"],
                                 "name" : ["response to stimulus"]})]
        
        self.terms_diam = self.terms_pro + [("Term", {"id" : ["GO:0048149"],
                                 "name" : ["behavioral response to ethanol"],
                                 "is_a" : ["GO:0030534", "GO:0045471"]})]
        
        self.terms_bad_rel = self.terms_diam + [("Term",
                            {"id" : ["GO:0005374"],
                             "name" : ["ethanol_intoxication"],
                             "relationship": ["regulates"]})]
        
        self.terms_cycle_no = self.terms_pro + [("Term", {"id" : ["GO:0048149"],
                                 "name" : ["behavioral response to ethanol"],
                                 "relationship": ["regulates GO:0005374"],
                                 "is_a" : ["GO:0030534", "GO:0045471"]}),
                                             ("Term", {"id" : ["GO:0005374"],
                                 "name" : ["ethanol intoxication"],
                                 "relationship": ["regulates GO:0048149"]})]
        
        self.terms_cycle = self.terms_cycle_no + [("Typedef",
                                                  {"id" : ["regulates"],
                                                   "name" : ["regulates"]})]

    def test_get_term(self):
        g = terms_to_graph(self.terms_diam)
        self.assertEqual("GO:0045471", g.get_term("GO:0045471").id)
        
    def test_get_parents(self):
        g = terms_to_graph(self.terms_diam)
        self.assertEqual(set(), g.get_parents("GO:0050896"))
        self.assertEqual(set(["GO:0030534", "GO:0045471"]),
                         g.get_parents("GO:0048149"))

    def test_get_ancestors(self):
        g = terms_to_graph(self.terms_diam)
        self.assertEqual(set(["GO:0030534", "GO:0045471", "GO:0050896"]),
                         g.get_ancestors("GO:0048149"))
        
        self.assertEqual(set(["GO:0050896"]), g.get_ancestors("GO:0030534"))
    
    
    def test_to_networkx(self):
        g = terms_to_graph(self.terms_diam)
        na = g.to_networkx({"GO:0050896" : ["X1", "X2"]})
        self.assertEqual(["GO:0030534", "GO:0045471"], na.successors("GO:0048149"))
        self.assertEqual(["GO:0050896"], na.successors("GO:0030534"))
        self.assertEqual(["GO:0050896"], na.successors("GO:0045471"))
        self.assertEqual(["X1", "X2"], na.node["GO:0050896"]["annotated_genes"])
        
    def test_bad_rel(self):
        with self.assertRaises(ValueError):
            terms_to_graph(self.terms_bad_rel)
            
    def test_no_typedef(self):
        with self.assertRaises(ValueError):
            terms_to_graph(self.terms_cycle_no)
    
    def test_get_parents_cycle(self):
        g = terms_to_graph(self.terms_cycle)
        self.assertEqual(set(['GO:0030534', 'GO:0050896', 'GO:0005374',
                              'GO:0048149', 'GO:0045471']),
                         g.get_ancestors("GO:0048149"))
        
        self.assertEqual([['GO:0048149', 'GO:0005374']], g.cycles)
        
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)