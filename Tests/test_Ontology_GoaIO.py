# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

import unittest
from Bio.Ontology.IO.GoaIO import GafIterator
from Bio.Ontology.Data import TermAssociation, GeneAnnotation

class GoaIOTest(unittest.TestCase):


    def test_no_version(self):
        with open('Ontology/GoaIO/no_ver.fb', 'r') as f:
            it = GafIterator(f)
            with self.assertRaises(ValueError):
                list(it)
                
    def test_incorrect_line_len(self):
        with open('Ontology/GoaIO/bad_line.fb', 'r') as f:
            it = GafIterator(f)
            with self.assertRaises(ValueError):
                list(it)

    def test_read_file(self):
        to = [GeneAnnotation('FB', 'FBgn0026615', '10-4', '10-4',
                       [], 'gene_product', ['taxon:7227'], [], '',
                       [TermAssociation([], 'GO:0005737', ['FB:FBrf0106275'],
                                      'IDA', [], 'C', '20060803', 'FlyBase' ),
                        TermAssociation([], 'GO:0045177', ['FB:FBrf0106275'],
                                      'IDA', [], 'C', '20060803', 'FlyBase' )]),
              GeneAnnotation('FB', 'FBgn0043467', '064Ya', '064Ya',
                       [], 'gene_product', ['taxon:7227'], [], '',
                       [TermAssociation([], 'GO:0048149',
                                      ['FB:FBrf0131396','PMID:11086999'],
                                      'IMP', [], 'P', '20060803', 'FlyBase')])]
        with open('Ontology/GoaIO/correct20.fb', 'r') as f:
            objs = list(GafIterator(f))
            objs.sort(key = lambda x : x.oid)
            self.assertEqual(to, objs)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)