# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

import unittest
import Bio.Ontology.IO as OntoIO
from Bio.Ontology.Data import OntologyGraph
from Bio.Ontology import EnrichmentFinder, RankedEnrichmentFinder


class EnrichmentFinderTest(unittest.TestCase):

    def setUp(self):
        gonodes_iter = OntoIO.parse("Ontology/go_test.obo", "obo")
        assocs_iter = OntoIO.parse("Ontology/ga_test.fb", "gaf")
        go_graph = OntologyGraph(gonodes_iter)
        self.ef = EnrichmentFinder(assocs_iter, go_graph)
    
    def test_find_terms_associations(self):
        expected = {'GO:0032502': set(['FBgn0010340']),'GO:0005737': set(['FBgn0026615']), 'GO:0045177': set(['FBgn0026615']), 'GO:0032501': set(['FBgn0010340']), 'GO:0044707': set(['FBgn0010340']), 'GO:0044464': set(['FBgn0026615']), 'GO:0005623': set(['FBgn0026615']), 'GO:0005622': set(['FBgn0026615']), 'GO:0005575': set(['FBgn0026615']), 'GO:0008150': set(['FBgn0010340', 'FBgn0026615']), 'GO:0044424': set(['FBgn0026615']), 'GO:0003674': set(['FBgn0026615']), 'GO:0007275': set(['FBgn0010340']), 'GO:0044699': set(['FBgn0010340'])}
        self.assertEqual(expected, self.ef._find_terms_associations(['FBgn0010340', 'FBgn0026615']))
        #TODO hand check and rest of tests
        
    def test_parent_child_union(self):
        pass
    
    def test_parent_child_intersection(self):
        pass
    
    def test_term_by_term(self):
        pass

class RankedEnrichmentFinderTest(unittest.TestCase):
    
    def setUp(self):
        gonodes_iter = OntoIO.parse("Ontology/go_test.obo", "obo")
        assocs_iter = OntoIO.parse("Ontology/ga_test.fb", "gaf")
        go_graph = OntologyGraph(gonodes_iter)
        self.ef = RankedEnrichmentFinder(assocs_iter, go_graph)
    
    def test_find_enriched_terms(self):
        expected = set(['GO:1901700', 'GO:0030534', 'GO:0003674', 'GO:0032501',
                        'GO:0044707', 'GO:0048149', 'GO:0044708', 'GO:0050896',
                        'GO:0007610', 'GO:0005575', 'GO:0045471', 'GO:0008150',
                        'GO:0042221', 'GO:0097305', 'GO:0010033', 'GO:0044699'])
        self.assertEqual(expected, self.ef._find_enriched_terms(['FBgn0043467', 'FBgn0004367']))
    
    
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)