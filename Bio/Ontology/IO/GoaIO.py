# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

import csv
import collections

from Bio.Ontology.GOData import GOAObject, GOAssociation



class TsvIterator(object):
    """
    Parses TSV files
    """

    def __init__(self, file_handle):
        self._reader = csv.reader(file_handle, delimiter='\t')

    def __iter__(self):
        return self._reader


class GafIterator(object):
    """
    Parses GAF files into list of GOAObject.
    
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
        self.records = None
    
    def _split_multi(self, value):
        if len(value) > 0:
            return value.split('|')
        else:
            return []
    
    def _to_goa(self, obj_rows):
        row = obj_rows[0]
        obj_params = [row[0], row[1], row[2], row[9], self._split_multi(row[10]),
                      row[11], self._split_multi(row[12])]
        
        if self.version == "1.0":
            row_len = 15
            obj_params += [None, None]
        else:
            row_len = 17
            obj_params += [self._split_multi(row[15]), row[16]]
            
        assocs = []
        for row in obj_rows:
            if len(row) == row_len:
                assocs.append(GOAssociation(self._split_multi(row[3]), row[4],
                                            self._split_multi(row[5]), row[6],
                                            self._split_multi(row[7]), row[8],
                                            row[13], row[14]))
            else:
                raise ValueError("Invalid gaf file: Incorrect row length.")
        
        return GOAObject(*obj_params, associations = assocs)
    
    def __iter__(self):
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
        self.records = [self._to_goa(v) for v in raw_records.values()]
        return iter(self.records)