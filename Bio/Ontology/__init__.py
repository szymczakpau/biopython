# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

class OntologyMachine(object):
    
    def __init__(self, assocs, go_graph, resolver_generator = None):
        self.go_graph = go_graph
        self.assocs = dict([(v.oid, v) for v in assocs])
        if resolver_generator != None:
            self.resolver = resolver_generator(self.assocs.itervalues())
        
        
    def get_go_terms(self, gene_list, depth = 0):
        pass # TODO resolved_genes = self.resolver.resolve(gene_list)
        