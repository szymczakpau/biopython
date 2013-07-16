# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.
"""
Module containing resolvers for ambigious entities ids.
"""


import collections

class Resolver(object):
    """
    Resolver which simply returns key given for disambiguation.
    """
    
    def __init__(self, assocs):
        pass
    
    def resolve(self, oid):
        return oid


class FirstOneResolver(object):
    """
    Resolver which picks first possible key for ambiguous entry.
    """
    
    def __init__(self, assocs):
        self.base_keys = set()
        alter = collections.defaultdict(set)
        for obj in assocs:
            self.base_keys.add(obj.oid)
            for aid in obj.synonyms:
                alter[aid].add(obj.oid)
        self.alter_keys = dict([(k,list(v)) for (k, v) in alter.iteritems()])
        
    
    def resolve(self, oid):
        if oid in self.base_keys:
            return oid
        elif oid in self.alter_keys:
            return self.alter_keys[oid][0]

class ManualResolver(object):
    """
    Manual resolver. Asks user for each ambiguous item.
    """
    pass

if __name__ == '__main__':
    pass