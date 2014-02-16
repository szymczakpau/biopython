# Copyright 2013 by Kamil Koziara. All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.  

from Bio._py3k import StringIO

import unittest

from Bio.Ontology.IO.OboIO import OboIterator, OboWriter
from Bio.Ontology.Data import OntologyTerm

class OboWriterTest(unittest.TestCase):
    
    def test_write(self):
        terms_to_write = [OntologyTerm("GO:0009628", "response to abiotic stimulus",
                                 {"is_a" : ["GO:0050896"]}),
                          OntologyTerm("GO:0022627", "cytosolic small ribosomal subunit",
                                 {"is_a" : ["GO:0015935", "GO:0044445"]})]
        f = StringIO()
        writer = OboWriter(f, version = "1.2")
        writer.write(terms_to_write)
        
        expected_output = """format-version:1.2

[Term]
id: GO:0009628
name: response to abiotic stimulus
is_a: GO:0050896

[Term]
id: GO:0022627
name: cytosolic small ribosomal subunit
is_a: GO:0015935
is_a: GO:0044445
"""
        self.assertEqual(expected_output, f.getvalue())     

class OboIteratorTest(unittest.TestCase):

    def test_strip_comment(self):
        a = "plastikowe fiolki "
        b = a + "! zakwitly ! "
        x = OboIterator(None)
        self.assertEqual(a, x._strip_comment(b))


    def test_line_reading(self):
        
        a = "Plastikowe fiolki zakwitly na zime."
        b = "Poczytaj Schopenhauera"
        line = "Plastikowe fiolki \\\nzakwitly na zime.! serio \nPoczytaj Schopenhauera"
        x = OboIterator(StringIO(line))
        self.assertEqual(a, x._read_line())
        self.assertEqual(b, x._read_line())
        self.assertEqual(None, x._read_line())
    
    def test_wrong_tag(self):
        a = """
[Term]
id : a
is_a
"""
        with self.assertRaises(ValueError):
            list(OboIterator(StringIO(a)))
    
    def test_simple_terms(self):
        a = """
format-version: 1.2
data-version: 2013-05-04
saved-by: kamil

[Term]
id : a
is_a : b

[Typedef]
id: y
name: x

[Term]
id : b \\
aha
is_a : c
is_a : d

[Typedef]
id: x
name: y

"""
        x = list(OboIterator(StringIO(a)))
        expected = [("Term", {"id" : ["a"], "is_a" : ["b"] }),
                    ("Typedef", {"id" : ["y"], "name" : ["x"]}),
                    ("Term", {"id" : ["b aha"], "is_a" : ["c", "d"] }),
                    ("Typedef", {"id" : ["x"], "name" : ["y"]})]
        self.assertEqual(expected, x)

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
