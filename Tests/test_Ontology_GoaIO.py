# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

import unittest
from Bio.Ontology.IO.GoaIO import GafReader
from Bio.Ontology.Data import TermAssociation, GeneAnnotation

class GoaIOTest(unittest.TestCase):


    def test_no_version(self):
        with open('Ontology/GoaIO/no_ver.fb', 'r') as f:
            it = GafReader(f)
            with self.assertRaises(ValueError):
                it.read()
                
    def test_incorrect_line_len(self):
        with open('Ontology/GoaIO/bad_line.fb', 'r') as f:
            it = GafReader(f)
            with self.assertRaises(ValueError):
                it.read()

    def test_read_file(self): #TODO
        to = [GeneAnnotation('FBgn0026615',
                             [TermAssociation('GO:0005737',
                                              {'Qualifier' : [],
                                            'DB:Reference' : ['FB:FBrf0106275'],
                                            'Evidence Code' : 'IDA',
                                            'With (or) From' : [],
                                            'Aspect' : 'C',
                                            'Date' : '20060803',
                                            'Assigned By' : 'FlyBase'}),
                              TermAssociation('GO:0045177',
                                              {'Qualifier' : [],
                                            'DB:Reference' : ['FB:FBrf0106275'],
                                            'Evidence Code' : 'IDA',
                                            'With (or) From' : [],
                                            'Aspect' : 'C',
                                            'Date' : '20060803',
                                            'Assigned By' : 'FlyBase'})],
                             {'DB' : 'FB', 'DB Object Symbol' : '10-4',
                              'DB Object Name' : '10-4',
                              'DB Object Synonym' : [],
                              'DB Object Type' : 'gene_product',
                              'Taxon' : ['taxon:7227'],
                              'Annotation Extension' : [],
                              'Gene Product Form ID' : ''}),
              GeneAnnotation('FBgn0043467',
                             [TermAssociation('GO:0048149',
                                              {'Qualifier' : [],
                                            'DB:Reference' : ['FB:FBrf0131396','PMID:11086999'],
                                            'Evidence Code' : 'IMP',
                                            'With (or) From' : [],
                                            'Aspect' : 'P',
                                            'Date' : '20060803',
                                            'Assigned By' : 'FlyBase'})],
                             {'DB' : 'FB',
                              'DB Object Symbol' : '064Ya',
                              'DB Object Name' : '064Ya',
                              'DB Object Synonym' : [],
                              'DB Object Type' : 'gene_product',
                              'Taxon' :['taxon:7227'],
                              'Annotation Extension' : [],
                              'Gene Product Form ID' : ''})]
        with open('Ontology/GoaIO/correct20.fb', 'r') as f:
            objs = GafReader(f).read()
            objs.sort(key = lambda x : x.id)
            self.assertEqual(to, objs)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)