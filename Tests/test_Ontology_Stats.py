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
            
    
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)