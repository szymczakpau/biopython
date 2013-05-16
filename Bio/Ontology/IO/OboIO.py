# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

START = 0
STANZA = 1
READ_STANZA = 2
EOF = 3


import collections
import re

class OboIterator:
    """
    Parses obo files.
    """
    
    def __init__(self, file_handle):
        self._handle = file_handle
        self._reg = re.compile("^\[([\w]+)\]$")
        self._state = START

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
        pos = line.find(':') # TODO add exception here if : not found
        return (line[0:pos].strip(), line[(pos+1):].strip())
    
    def _read_stanza(self):
        while self._state != STANZA and self._state != EOF:
            line = self._read_line()
            if line is None:
                self._state = EOF
            else:
                line = line.strip()
                match = self._reg.match(line)
                if match != None:
                    self._state = STANZA
                    self._stanza_type, = match.groups()
                elif len(line) > 0 and self._state == READ_STANZA:
                    k, v = self._split_tag(line)
                    self._dict[k].append(v)


    def next(self):
        self._read_stanza()
        if self._state == EOF:
            raise StopIteration
        self._dict = collections.defaultdict(list)
        self._state = READ_STANZA
        self._read_stanza()
        return (self._stanza_type, self._dict)
