# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

import collections
import Stats, IdResolver


class EnrichmentEntry(object):
    """
    Represents result returned by EnrichmentFinder.
    """
    def __init__(self, oid, name, p_value, study_count, study_n, population_count, population_n):
        self.oid = oid
        self.name = name
        self.p_value = p_value
        self.study_count = study_count
        self.study_n = study_n
        self.population_count = population_count
        self.population_n = population_n
        self.corrections = []
        
    def __repr__(self):
        return 'EnrichmentEntry([("ID" : {0}), ("name" : {2}), ("p-value" : {1})'.format(self.oid, self.p_value, self.name)

    def __str__(self):
        return """ID : {0}
name : {6}
p-value : {1}
corrected p-values: {7}
hits in study : {2}
elements in study : {3}
hits in population : {4}
population count : {5}""".format(self.oid, self.p_value, self.study_count, self.study_n,
           self.population_count, self.population_n, self.name, self.corrections)


class EnrichmentFinder(object):
    """
    Utility for finding enrichment of group of entities given population,
    ontology graph and associations to this graph.
    
    EnrichmentFinder could be used for finding enrichment in many cases of
    entities but in this example let's focus on finding gene enrichment in
    gene ontology graph.

    To initialize the finder you need at least ontology graph and associations
    from genes to it's nodes. Typically you would read them from files
    using readers from Bio.Ontology.IO module.
    
    >>> import Bio.Ontology.IO as OntoIO
    >>> gonodes_iter = OntoIO.parse("Ontology/go_test.obo", "obo")
    >>> assocs_iter = OntoIO.parse("Ontology/ga_test.fb", "gaf")
    
    Now you've got iterators to files containing both gene ontology graph and
    associations. You need to init the graph. 
    
    >>> from Bio.Ontology.GOData import GOGraph
    >>> go_graph = GOGraph(gonodes_iter)
    
    Now is the time to create EnrichmentFinder. Besides the arguments
    mentioned before you could specify an id resolver, which basically tries
    to find synonyms for gene ids that finder can understand, and population
    of genes as reference. If none of this is specified default resolver
    is used and all genes from association are used as the population.
    
    >>> from Bio.Ontology import EnrichmentFinder
    >>> ef = EnrichmentFinder(assocs_iter, go_graph)
    
    To run finder you just need to call find_enrichment method with list
    of genes as an argument:
    
    >>> genes_to_study = ['FBgn0070057', 'FBgn0004364', 'FBgn0043467']
    
    Additionally you can add a list of corrections which should be used
    (by default it's empty):
    
    >>> corrections = ['bonferroni', 'bh_fdr']
    >>> result = ef.find_enrichment(genes_to_study, corrections)
    
    The result is a list of EnrichmentEntry instances. Each for a node in graph
    which is enriched.
    
    >>> print result[0]
    ID : GO:0016020
    name : membrane
    p-value : 0.333333333333
    corrected p-values: [('bonferroni', 1.0), ('bh_fdr', 1.0)]
    hits in study : 1
    elements in study : 3
    hits in population : 1
    population count : 9
    
    """
    
    def __init__(self, assocs, o_graph, population = None,
                 resolver_generator = IdResolver.FirstOneResolver):
        """
        Initialize Enrichment Finder.
        
        Parameters
        ----------
        assocs - iterable containing associations
        o_graph - graph with ontology
        population - population used as reference to study sample
        resolver_generator - constructor of resolver used to disambiguate ids
        """
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
    
    def find_enrichment(self, gene_list, corrections = []):
        """
        Finds enrichment of specified group of genes.
        
        Parameters
        ----------
        gene_list - list of genes to study
        corrections - list of corrections that should be applied to result.
            Possible values are:
                o "bonferroni" - Bonferroni correction,
                o "fh_fdr" - Benjamin-Hochberg FDR correction.
        """
        
        resolved_list = [self.resolver.resolve(x) for x in gene_list]
        study_counts = self._count_terms(resolved_list)
        
        population_n = len(self.population)
        study_n = len(resolved_list)
        
        result = []
        
        # Calculate enrichment for every term given the counts
        for term, study_count in study_counts.items():
            population_count = self.population_counts[term]
            pval = Stats.hypergeometric_test(study_count, study_n,
                                             population_count, population_n)

            entry = EnrichmentEntry(term, self.o_graph.get_term(term).name,
                                    pval, study_count, study_n,
                                    population_count, population_n)

            result.append(entry)
        
        # Calculate chosen corrections
        if len(corrections) > 0:
            pvals = [x.p_value for x in result]
            corr_pvals = []
            for c_id in corrections:
                cfun = Stats.corrections[c_id]
                corr_pvals.append((c_id, cfun(pvals)))
            print corr_pvals
            for i in xrange(len(result)):
                result[i].corrections = [(c_id, pv[i]) for c_id, pv in corr_pvals] #TODO: tests
        return result
    
if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()