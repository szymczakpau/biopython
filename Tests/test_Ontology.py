# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

"""
Basic tests for EnrichmentFinder and Ranked EnrichmentFinder.
"""

import unittest
import Bio.Ontology.IO as OntoIO
from Bio.Ontology import EnrichmentFinder, RankedEnrichmentFinder

class EnrichmentFinderTest(unittest.TestCase):

    def setUp(self):
        go_graph = OntoIO.read("Ontology/go_test.obo", "obo")
        assocs = OntoIO.read("Ontology/ga_test.fb", "gaf")
        self.ef = EnrichmentFinder(assocs, go_graph)
    
    def test_find_terms_associations(self):
        expected = {'GO:0032502': set(['FBgn0010340']),'GO:0005737': set(['FBgn0026615']), 'GO:0045177': set(['FBgn0026615']), 'GO:0032501': set(['FBgn0010340']), 'GO:0044707': set(['FBgn0010340']), 'GO:0044464': set(['FBgn0026615']), 'GO:0005623': set(['FBgn0026615']), 'GO:0005622': set(['FBgn0026615']), 'GO:0005575': set(['FBgn0026615']), 'GO:0008150': set(['FBgn0010340', 'FBgn0026615']), 'GO:0044424': set(['FBgn0026615']), 'GO:0003674': set(['FBgn0026615']), 'GO:0007275': set(['FBgn0010340']), 'GO:0044699': set(['FBgn0010340'])}
        self.assertEqual(expected, self.ef._find_terms_associations(['FBgn0010340', 'FBgn0026615']))
        
    def test_parent_child_union(self):
        genes = ['FBgn0043467', 'FBgn0010339', 'FBgn0070057', 'FBgn0070052']
        en = self.ef.find_enrichment(genes, ['bh_fdr'], 'parent_child_union')
        self.assertEqual(7, len(en.filter_p_val(0.9).entries))
    
    def test_parent_child_intersection(self):
        genes = ['FBgn0043467', 'FBgn0010339', 'FBgn0070057', 'FBgn0070052']
        en = self.ef.find_enrichment(genes, ['bonferroni'], 'parent_child_intersection')
        self.assertEqual(7, len(en.filter_p_val(0.9).entries))
        
    def test_term_for_term(self):
        genes = ['FBgn0043467', 'FBgn0010339', 'FBgn0070057', 'FBgn0070052']
        en = self.ef.find_enrichment(genes)
        self.assertEqual(24, len(en.filter_p_val(0.9).entries))

class RankedEnrichmentFinderTest(unittest.TestCase):
    
    def setUp(self):
        go_graph = OntoIO.read("Ontology/go_test.obo", "obo")
        assocs_iter = OntoIO.read("Ontology/ga_test.fb", "gaf")
        self.ef = RankedEnrichmentFinder(assocs_iter, go_graph)

    def test_ranked_parent_child(self):
        gene_rank = [('FBgn0043467', 0.1), ('FBgn0010339', 0.7), ('FBgn0070057', 0.4), ('FBgn0070052', 0.9)]
        en = self.ef.find_enrichment_parent_child(gene_rank, "+", ["bh_fdr"], False, "parent_child_intersection")
        self.assertEqual(8, len(en.filter_p_val(0.9).entries))

    def test_gsea(self):
        gene_rank = [('FBgn0043467', 0.1), ('FBgn0010339', 0.7), ('FBgn0070057', 0.4), ('FBgn0070052', 0.9)]
        en = self.ef.find_enrichment_from_rank(gene_rank, 10, 2)
        self.assertEqual(6, len(en.entries))
        
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)