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
    Parses GAF files
    """
    
    
    _GAF_LABELS = ['DB', 'DB Object ID', 'DB Object Symbol', 'Qualifier', 'GO ID',
                   'DB:Reference', 'Evidence Code', 'With (or) From', 'Aspect',
                   'DB Object Name', 'DB Object Synonym', 'DB Object Type',
                   'Taxon', 'Date', 'Assigned By', 'Annotation Extension',
                   'Gene Product Form ID']

    _LABEL_MAP = { "1.0" : _GAF_LABELS[:15], "2.0" : _GAF_LABELS }
    

    _ID_IDX = 1
    _AID_IDX = 10

    _OBJ_IDXS = [0, 1, 2, 9, 11, 12, 15, 16]
    _AS_IDXS = [3, 4, 5, 6, 7, 8, 13, 14]
    
    def __init__(self, file_handle):
        tsv_iter = TsvIterator(file_handle)
        raw_records = collections.defaultdict(list)
        for row in tsv_iter:
            first = row[0]
            if not first.startswith('!'):
                raw_records[row[self._ID_IDX]].append(row)
            elif first.startswith('!gaf-version:'):
                    self.version = first[(first.find(':') + 1):].strip()
        self.records = [self._to_goa(v) for v in raw_records.values()]
        
    def _to_goa(self, obj_rows):
        assocs = [GOAssociation(*[row[i] for i in self._AS_IDXS]) for row in obj_rows]
        if len(row[self._AID_IDX]) > 0:
            syns = row[self._AID_IDX].split('|')
        else:
            syns = []
        
        row = obj_rows[0] # TODO sprawdzic czy wszedzie gra
        
        if self.version == "1.0":
            idxs = self._OBJ_IDXS[:-2]
        else:
            idxs = self._OBJ_IDXS
            
        return GOAObject(*[row[i] for i in idxs], synonyms = syns, associations = assocs)
    
    def __iter__(self):
        return iter(self.records)