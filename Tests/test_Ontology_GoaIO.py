# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

import unittest
from Bio.Ontology.IO.GoaIO import GafReader, InSqlAssoc, GAF20FIELDS, TsvIterator
from Bio.Ontology.Data import TermAssociation, GeneAnnotation


class InSqlAssocTest(unittest.TestCase):

    def test_insert_get_len(self):
        assoc = InSqlAssoc(GAF20FIELDS, [1, 4], lambda x: x)
        self.assertEquals(0, len(assoc))
        
        with open('Ontology/GoaIO/no_ver.fb', 'r') as f:
            for r in TsvIterator(f):
                if len(r) == 17:
                    assoc.add_row(r)
                    
        self.assertEquals(3, len(assoc))
        self.assertEquals([(u'FB', u'FBgn0043467', u'064Ya', u'', u'GO:0005575', u'FB:FBrf0159398', u'ND', u'', u'C', u'064Ya', u'', u'gene_product', u'taxon:7227', u'20060803', u'FlyBase', u'', u''), (u'FB', u'FBgn0043467', u'064Ya', u'', u'GO:0048149', u'FB:FBrf0131396|PMID:11086999', u'IMP', u'', u'P', u'064Ya', u'', u'gene_product', u'taxon:7227', u'20060803', u'FlyBase', u'', u'')],
                          sorted(assoc['FBgn0043467']))
        
        self.assertTrue('FBgn0043467' in assoc)
    
    def test_iter(self):
        assoc = InSqlAssoc(GAF20FIELDS, [1, 4], lambda x: x)
        
        with open('Ontology/GoaIO/no_ver.fb', 'r') as f:
            for r in TsvIterator(f):
                if len(r) == 17:
                    assoc.add_row(r)
        
        self.assertEqual([(u'FBgn0026619', [(u'FB', u'FBgn0026615', u'10-4', u'', u'GO:0005737', u'FB:FBrf0106275', u'IDA', u'', u'C', u'10-4', u'', u'gene_product', u'taxon:7227', u'20060803', u'FlyBase', u'', u''), (u'FB', u'FBgn0026615', u'10-4', u'', u'GO:0045177', u'FB:FBrf0106275', u'IDA', u'', u'C', u'10-4', u'', u'gene_product', u'taxon:7227', u'20060803', u'FlyBase', u'', u'')]), (u'FBgn0043467', [(u'FB', u'FBgn0026619', u'10-4', u'', u'GO:0045177', u'FB:FBrf0106275', u'IDA', u'', u'C', u'10-4', u'', u'gene_product', u'taxon:7227', u'20060803', u'FlyBase', u'', u'')]), (u'FBgn0043467', [(u'FB', u'FBgn0043467', u'064Ya', u'', u'GO:0048149', u'FB:FBrf0131396|PMID:11086999', u'IMP', u'', u'P', u'064Ya', u'', u'gene_product', u'taxon:7227', u'20060803', u'FlyBase', u'', u''), (u'FB', u'FBgn0043467', u'064Ya', u'', u'GO:0005575', u'FB:FBrf0159398', u'ND', u'', u'C', u'064Ya', u'', u'gene_product', u'taxon:7227', u'20060803', u'FlyBase', u'', u'')])],
                         list(assoc))
            
            
class GoaIOTest(unittest.TestCase):


    def test_no_version(self):
        with open('Ontology/GoaIO/no_ver.fb', 'r') as f:
            it = GafReader(f, assoc_format = "in_mem_sql")
            with self.assertRaises(ValueError):
                it.read()
                
    def test_incorrect_line_len(self):
        with open('Ontology/GoaIO/bad_line.fb', 'r') as f:
            it = GafReader(f)
            with self.assertRaises(ValueError):
                it.read()

    def test_read_file(self):
        to = {'FBgn0026615' : GeneAnnotation('FBgn0026615',
                             [TermAssociation('GO:0005737',
                                           {GAF20FIELDS[3] : [],
                                            GAF20FIELDS[5] : ['FB:FBrf0106275'],
                                            GAF20FIELDS[6] : 'IDA',
                                            GAF20FIELDS[7] : [],
                                            GAF20FIELDS[8] : 'C',
                                            GAF20FIELDS[13] : '20060803',
                                            GAF20FIELDS[14] : 'FlyBase'}),
                              TermAssociation('GO:0045177',
                                           {GAF20FIELDS[3] : [],
                                            GAF20FIELDS[5] : ['FB:FBrf0106275'],
                                            GAF20FIELDS[6] : 'IDA',
                                            GAF20FIELDS[7] : [],
                                            GAF20FIELDS[8] : 'C',
                                            GAF20FIELDS[13] : '20060803',
                                            GAF20FIELDS[14] : 'FlyBase'})],
                             {GAF20FIELDS[0] : 'FB',
                              GAF20FIELDS[2] : '10-4',
                              GAF20FIELDS[9] : '10-4',
                              GAF20FIELDS[10] : [],
                              GAF20FIELDS[11] : 'gene_product',
                              GAF20FIELDS[12] : ['taxon:7227'],
                              GAF20FIELDS[15] : [],
                              GAF20FIELDS[16] : ''}),
                'FBgn0043467' : GeneAnnotation('FBgn0043467',
                             [TermAssociation('GO:0048149',
                                           {GAF20FIELDS[3] : [],
                                            GAF20FIELDS[5] : ['FB:FBrf0131396','PMID:11086999'],
                                            GAF20FIELDS[6] : 'IMP',
                                            GAF20FIELDS[7] : [],
                                            GAF20FIELDS[8] : 'P',
                                            GAF20FIELDS[13] : '20060803',
                                            GAF20FIELDS[14] : 'FlyBase'})],
                             {GAF20FIELDS[0] : 'FB',
                              GAF20FIELDS[2] : '064Ya',
                              GAF20FIELDS[9] : '064Ya',
                              GAF20FIELDS[10] : [],
                              GAF20FIELDS[11] : 'gene_product',
                              GAF20FIELDS[12] : ['taxon:7227'],
                              GAF20FIELDS[15] : [],
                              GAF20FIELDS[16] : ''})}
        with open('Ontology/GoaIO/correct20.fb', 'r') as f:
            objs = GafReader(f).read()
            self.assertEqual(to, objs)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)