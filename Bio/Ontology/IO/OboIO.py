# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

_START = 0
_STANZA = 1
_READ_STANZA = 2
_EOF = 3


import collections
import re

from Bio.Ontology.Data import OntologyTerm

class OboWriter(object):
    """
    Writes obo files.
    
    Writes OntologyTerms to obo files.
    """
    
    def __init__(self, file_handle):
        self.handle = file_handle
    
    def write_file(self, terms_list, version):
        # now only terms are valid for writing
        if version != None:
            self.handle.write("format-version:" + version + "\n")
        for term in terms_list:
            if isinstance(term, OntologyTerm):
                self.handle.write("\n[Term]\n")
                self.handle.write("id: " + term.id + "\n")
                self.handle.write("name: " + term.name + "\n")
                for k, v in term.attrs.iteritems():
                    for vi in v:
                        self.handle.write(k + ": " + vi + "\n")
                
        
class OboIterator(object):
    """
    Parses obo files.
    """
    
    def __init__(self, file_handle):
        self._handle = file_handle
        self._reg = re.compile("^\[([\w]+)\]$")
        self._state = _START

    def __iter__(self):
        return self
    
    def _strip_comment(self, line):
        pos = line.find('!')
        if pos > 0:
            return line[0:pos]
        else:
            return line
    
    def _read_line(self):
        iterate = True
        result = ''
        while iterate:
            line = self._handle.readline()
            if line:
                line = self._strip_comment(line)
                pos = line.find("\\\n")
                if pos > 0:
                    line = line[:pos]
                else:
                    iterate = False
                result += line
            else:
                return None
        return result
    
    def _split_tag(self, line):
        pos = line.find(':')
        if pos < 0:
            raise ValueError(("Invalid obo file: Incorrect tag: ':'"
                              " expected in line '{0}'.").format(line))
        return (line[0:pos].strip(), line[(pos+1):].strip())
    
    def _read_stanza(self):
        while self._state != _STANZA and self._state != _EOF:
            line = self._read_line()
            if line is None:
                self._state = _EOF
            else:
                line = line.strip()
                match = self._reg.match(line)
                if match != None:
                    self._state = _STANZA
                    self._found_stanza_type, = match.groups()
                elif len(line) > 0 and self._state == _READ_STANZA:
                    k, v = self._split_tag(line)
                    self._dict[k].append(v)


    def next(self):
        self._read_stanza()
        if self._state == _EOF:
            raise StopIteration
        self._dict = collections.defaultdict(list)
        self._state = _READ_STANZA
        stanza_type = self._found_stanza_type
        self._read_stanza()
        return (stanza_type, self._dict)
