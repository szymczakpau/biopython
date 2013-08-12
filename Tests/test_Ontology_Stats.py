# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.


import unittest
from math import *
from Bio.Ontology.Stats import *

class StatisticalFunctionsTest(unittest.TestCase):
    
    def test_lnfactorial(self):
        n = 1
        for i in xrange(1, 1000):
            n *= i
            self.assertAlmostEqual(log(n), lnfactorial(i), 11)
    
    def test_lncombination(self):
        for n in xrange(45):
            s = 0
            for k in xrange(n + 1):
                s += exp(lncombination(n, k))
            self.assertEqual(pow(2, n), round(s))
    
    def test_bonferroni_correction(self):
        expected = [0.5, 0.01, 1.0, 0.2, 0.001]
        computed = bonferroni_correction([0.1, 0.002, 0.3, 0.04, 0.0002])
        for i in xrange(len(expected)):
            self.assertEqual(expected[i], computed[i])
            
    def test_bh_fdr_correction(self):
        expected = [0.16, 0.05, 0.0625, 0.0625, 0.0625]
        computed = bh_fdr_correction([0.16, 0.01, 0.027, 0.025, 0.026])
        for i in xrange(len(expected)):
            self.assertEqual(expected[i], computed[i])
    
    def test_kolmogorov_smirnov_rank_test(self):
        genes_list = ['YIL152W', 'YDL077C', 'YIL149C', 'YNL110C',
                      'YBR043C', 'YLR033W', 'YDL096C', 'YNL039W',
                      'YNL136W', 'YIL160C']
        gene_corr = [-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4]
        genes_set = set([genes_list[i] for i in [1, 2, 4]] + ['YBB0101'])
        
        val, plot =  kolmogorov_smirnov_rank_test(genes_set, genes_list, gene_corr, 1)
        
        self.assertAlmostEqual(val, 0.7321428571428572, 15)
        self._almostAssertLists(plot, [-0.14285714285714285, 0.35714285714285726,
            0.7321428571428572, 0.5892857142857144, 0.7142857142857144,
            0.5714285714285716, 0.42857142857142877, 0.2857142857142859,
            0.14285714285714307, 0], 15)
        
    def _almostAssertLists(self, la, lb, places):
        self.assertEqual(len(la), len(lb), msg = "List not equal.")
        for i in xrange(len(la)):
            self.assertAlmostEqual(la[i], lb[i], places)
        
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)