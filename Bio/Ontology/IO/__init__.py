# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

from Bio.File import as_handle

import OboIO
import GoaIO

_FormatToIterator = { "obo" : OboIO.OboIterator,
                      "tsv" : GoaIO.TsvIterator,
                      "gaf" : GoaIO.GafIterator }

_FormatToWriter = {"obo" : OboIO.OboWriter }

def write(data, handle, file_format, version = None):
    """
    Write an ontology data to file.

    Arguments:
     - data - data to write to a file,
     - handle - File handle object to write to, or filename as string
                   (note older versions of Biopython only took a handle),
     - file_format - lower case string describing the file format to write,
     - version - file format version .
     
    You should close the handle after calling this function.

    """
    
    if not isinstance(file_format, basestring):
        raise TypeError("Need a string for the file format (lower case)")
    if not file_format:
        raise ValueError("Format required (lower case string)")

    with as_handle(handle, 'w') as fp:
        #Map the file format to a writer class
        if file_format in _FormatToWriter:
            writer_class = _FormatToWriter[file_format]
            writer_class(fp).write_file(data, version)
        else:
            raise ValueError("Unknown format '%s'" % file_format)

def parse(handle, file_format):
    """
    Iterate over a gene ontology file.
    """

    if not isinstance(file_format, basestring):
        raise TypeError("Need a string for the file format (lower case)")
    if not file_format:
        raise ValueError("Format required (lower case string)")          
    if file_format != file_format.lower():
        raise ValueError("Format string '%s' should be lower case" % format)
    with as_handle(handle, 'rU') as fp:
        if file_format in _FormatToIterator:
            iterator_generator = _FormatToIterator[file_format]
            it = iterator_generator(fp)

            for el in it:
                yield el
        else:
            raise ValueError("Unknown format '%s'" % file_format)
