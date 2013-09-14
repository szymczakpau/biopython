# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

import unittest
from StringIO import StringIO
from Bio.Ontology.IO.EnrichmentIO import EnrichmentWriter, EnrichmentReader
from Bio.Ontology import EnrichmentEntry, Enrichment

class EnrichmentWriterTest(unittest.TestCase):

    def test_write(self):
        result = StringIO()

        e1 = EnrichmentEntry("9951", "structure-specific DNA binding", 0.032301032301)
        e1.corrections = {'bh_fdr': 1.0, 'bonferroni': 1.0}
        e1.attrs = {'plot' : [0.1, 0.2, 1.0, 0.1], 'score' : 1.0}
        
        e2 = EnrichmentEntry("9916", "polysomal ribosome", 0.025)
        e2.corrections = {'bh_fdr': 1.0, 'bonferroni': 1.0}
        e2.attrs = {'plot' : [0.1, 0.2, 1.0, 0.1], 'score' : 1.0}
        
        en = Enrichment("ranked parent-child", [e1, e2], ["Cycles found..."], ['bh_fdr', 'bonferroni'])
        
        writer = EnrichmentWriter(result)
        writer.write_file(en)
        expected = ("# ranked parent-child\r\n"
                    "# 2 1\r\n"
                    "id\tname\tp-value\tbh_fdr|bonferroni\tattributes\r\n"
                    "9951\tstructure-specific DNA binding\t0.032301032301\t1.0|1.0\t{'plot': [0.1, 0.2, 1.0, 0.1], 'score': 1.0}\r\n"
                    "9916\tpolysomal ribosome\t0.025\t1.0|1.0\t{'plot': [0.1, 0.2, 1.0, 0.1], 'score': 1.0}\r\n"
                    "!\tCycles found...\r\n")
        self.assertEqual(expected, result.getvalue())


class EnrichmentReaderTest(unittest.TestCase):

    def test_read_correct(self):
        inputstr = ("# ranked parent-child\r\n"
                    "# 2 1\r\n"
                    "id\tname\tp-value\tbh_fdr|bonferroni\tattributes\r\n"
                    "9951\tstructure-specific DNA binding\t0.032301032301\t1.0|1.0\t{'plot': [0.1, 0.2, 1.0, 0.1], 'score': 1.0}\r\n"
                    "9916\tpolysomal ribosome\t0.025\t1.0|1.0\t{'plot': [0.1, 0.2, 1.0, 0.1], 'score': 1.0}\r\n"
                    "!\tCycles found...\r\n")
        test_input = StringIO(inputstr)

        e1 = EnrichmentEntry("9951", "structure-specific DNA binding", 0.032301032301)
        e1.corrections = {'bh_fdr': 1.0, 'bonferroni': 1.0}
        e1.attrs = {'plot' : [0.1, 0.2, 1.0, 0.1], 'score' : 1.0}
        
        e2 = EnrichmentEntry("9916", "polysomal ribosome", 0.025)
        e2.corrections = {'bh_fdr': 1.0, 'bonferroni': 1.0}
        e2.attrs = {'plot' : [0.1, 0.2, 1.0, 0.1], 'score' : 1.0}
        
        expected_en = Enrichment("ranked parent-child", [e1, e2], ["Cycles found..."], ['bh_fdr', 'bonferroni'])
        
        reader = EnrichmentReader(test_input)
        en = reader.read()
        
        self.assertEqual(expected_en, en)
        
if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)