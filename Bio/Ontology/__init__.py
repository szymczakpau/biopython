# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

import collections
import Stats

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

class EnrichmentEntry(object):

    def __init__(self, oid, p_value, study_count, study_n, population_count, population_n):
        self.oid = oid
        self.p_value = p_value
        self.study_count = study_count
        self.study_n = study_n
        self.population_count = population_count
        self.population_n = population_n

    def __repr__(self):
        return 'EnrichmentEntry([("ID" : {0}), ("p-value" : {1}), ("hits in study" : {2}), ("elements in study" : {3}), ("hits in population" : {4}), ("population count" : {5})])'.format(self.oid, self.p_value, self.study_count, self.study_n, self.population_count, self.population_n)

    def __str__(self):
        return """ID : {0}
p-value : {1}
hits in study : {2}
elements in study : {3}
hits in population : {4}
population count : {5}
""".format(self.oid, self.p_value, self.study_count, self.study_n, self.population_count, self.population_n)

class EnrichmentFinder(object):
    
    def __init__(self, assocs, o_graph, population = None, resolver_generator = FirstOneResolver):
        self.o_graph = o_graph
        self.assocs = dict([(v.oid, v) for v in assocs])
        self.resolver = resolver_generator(self.assocs.itervalues())
        if population != None:
            self.population = population
        else:
            self.population = self.assocs.keys()
        self.population_counts = self._count_terms(self.population)
        
    def _count_terms(self, gene_list):
        terms_counts = collections.defaultdict(int)
        for gene in gene_list:
            if gene in self.assocs:
                enriched_terms = set()
                for go_terms in self.assocs[gene].associations: # qualifier co z nim?
                    enriched_terms.add(go_terms.go_id) # czy mozliwy jest mapping do ancestora taki podwojny
                    enriched_terms |= self.o_graph.get_ancestors(go_terms.go_id)
                for t in enriched_terms:
                    terms_counts[t] += 1
        return terms_counts
    
    def find_enrichment(self, gene_list):
        resolved_list = [self.resolver.resolve(x) for x in gene_list]
        study_counts = self._count_terms(resolved_list)
        
        population_n = len(self.population)
        study_n = len(resolved_list)
        
        result = []
        
        for term, study_count in study_counts.items():
            population_count = self.population_counts[term]
            pval = Stats.hypergeometric_test(study_count, study_n, population_count, population_n)

            entry = EnrichmentEntry(term, pval,
                                    study_count, study_n, population_count,
                                    population_n)

            result.append(entry)
        return result
    
    def get_go_terms(self, gene_list):
        resolved_list = [self.resolver.resolve(x) for x in gene_list]
        result = []
        for gene in resolved_list:
            obj = self.assocs[gene]
            ascs = []
            for asc in obj.associations:
                gnode = self.o_graph.get_node(asc.go_id)
                ascs.append((asc, gnode))
            result.append((obj, ascs))
        return result
    

    
