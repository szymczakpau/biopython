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
    
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)