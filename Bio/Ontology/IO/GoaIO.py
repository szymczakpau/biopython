# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

"""
I/O operations for gene annotation files.
"""

import csv
import collections
from Bio.Ontology.Data import GeneAnnotation, TermAssociation
from Interfaces import OntoIterator, OntoReader

class TsvIterator(OntoIterator):
    """
    Parses TSV files
    """

    def __init__(self, file_handle):
        self._reader = csv.reader(file_handle, delimiter='\t')

    def __iter__(self):
        return self._reader


class GafReader(OntoReader):
    """
    Reads GAF files into list of GeneAnnotation.
    
    GAF file is list of tab separated values in the following order:
        'DB', 'DB Object ID', 'DB Object Symbol', 'Qualifier', 'GO ID',
        'DB:Reference', 'Evidence Code', 'With (or) From', 'Aspect',
        'DB Object Name', 'DB Object Synonym', 'DB Object Type',
        'Taxon', 'Date', 'Assigned By', 'Annotation Extension',
        'Gene Product Form ID'
    """

    _ID_IDX = 1
    
    def __init__(self, file_handle):
        self.handle = file_handle
        self.version = None
    
    def _split_multi(self, value):
        if len(value) > 0:
            return value.split('|')
        else:
            return []
    
    def _to_goa(self, obj_rows):
        row = obj_rows[0]
        
        obj_id = row[1]
        obj_attrs = {'DB' : row[0], 'DB Object Symbol' : row[2],
                      'DB Object Name' : row[9],
                      'DB Object Synonym' : self._split_multi(row[10]),
                      'DB Object Type' : row[11],
                      'Taxon' : self._split_multi(row[12])}
        
        if self.version == "1.0":
            row_len = 15
        else:
            row_len = 17
            obj_attrs['Annotation Extension'] = self._split_multi(row[15])
            obj_attrs['Gene Product Form ID'] = row[16]
            
        assocs = []
        for row in obj_rows:
            if len(row) == row_len:
                assocs.append(TermAssociation(row[4],
                                           {'Qualifier' : self._split_multi(row[3]),
                                            'DB:Reference' : self._split_multi(row[5]),
                                            'Evidence Code' : row[6],
                                            'With (or) From' :self._split_multi(row[7]),
                                            'Aspect' : row[8],
                                            'Date' : row[13],
                                            'Assigned By' : row[14]}
                                              ))
            else:
                raise ValueError("Invalid gaf file: Incorrect row length.")
        
        return GeneAnnotation(obj_id, assocs, obj_attrs)
    
    def read(self):
        tsv_iter = TsvIterator(self.handle)
        raw_records = collections.defaultdict(list)
        for row in tsv_iter:
            first = row[0]
            if not first.startswith('!'):
                raw_records[row[self._ID_IDX]].append(row)
            elif first.startswith('!gaf-version:'):
                self.version = first[(first.find(':') + 1):].strip()
        if self.version is None:
            raise ValueError("Invalid gaf file: No version specified.")
        return [self._to_goa(v) for v in raw_records.values()]