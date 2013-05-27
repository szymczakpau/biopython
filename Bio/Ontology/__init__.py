# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

import collections

class Resolver(object):
    """
    Default resolver. Simply returns key given for disambiguation.
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
            self.alter_keys[oid][0]

class ManualResolver(object):
    """
    Manual resolver. Asks user for each ambiguous item.
    """
    pass


class OntologyMachine(object):
    
    def __init__(self, assocs, go_graph, resolver_generator = Resolver):
        self.go_graph = go_graph
        self.assocs = dict([(v.oid, v) for v in assocs])
        self.resolver = resolver_generator(self.assocs.itervalues())
        
        
    def get_go_terms(self, gene_list):
        resolved_list = [self.resolver.resolve(x) for x in gene_list]
        result = []
        for gene in resolved_list:
            obj = self.assocs[gene]
            ascs = []
            for asc in obj.associations:
                gnode = self.go_graph.get_node(asc.go_id)
                ascs.append((asc, gnode))
            result.append((obj, ascs))
        return result
    

    
