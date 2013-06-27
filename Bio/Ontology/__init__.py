# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

import collections
import Stats, IdResolver


class EnrichmentEntry(object):
    """
    Represents one result returned by EnrichmentFinder.
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

class Enrichment(object):
    """
    Contains all results found by EnrichmentFinder
    """
    
    def __init__(self, method, entries, warnings):
        self.entries = entries
        self.warnings = warnings
        self.method = method
    def __repr__(self):
        return "Enrichment(method = {2}, entries_num = {0} , warnings_num = {1})".format(len(self.entries), len(self.warnings), self.method)

    def __str__(self):
        return "Enrichment using {2} method: {0} entries, {1} warnings.".format(len(self.entries),
                                                               len(self.warnings), self.method)

_METHOD_PARENT_CHILD = "parent-child"

_METHOD_TERM_BY_TERM = "term_by_term"

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
    
    The result contains a list of EnrichmentEntry instances - each for a node in graph
    which is enriched - and a list of warnings.
    
    >>> print result
    Enrichment using term_by_term method: 64 entries, 0 warnings.
    >>> print result.entries[0]
    ID : GO:0044707
    name : single-multicellular organism process
    p-value : 1.0
    corrected p-values: [('bonferroni', 1.0), ('bh_fdr', 1.0)]
    hits in study : 2
    elements in study : 3
    hits in population : 5
    population count : 9
    
    You can also specify a method of computing p-values. Default method is
    term by term, but additionaly you can choose method that takes parent-child
    relationship into account:
    
    >>> result = ef.find_enrichment(genes_to_study, corrections, "parent-child")
    
    >>> print result
    Enrichment using parent-child method: 64 entries, 0 warnings.
    >>> print result.entries[0]
    ID : GO:0044707
    name : single-multicellular organism process
    p-value : 1.0
    corrected p-values: [('bonferroni', 1.0), ('bh_fdr', 1.0)]
    hits in study : 2
    elements in study : 2
    hits in population : 5
    population count : 5
    
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
        terms_counts = collections.defaultdict(set)
        for gene in gene_list:
            if gene in self.assocs:
                enriched_terms = set()
                for go_terms in self.assocs[gene].associations: # qualifier co z nim?
                    enriched_terms.add(go_terms.go_id)
                    enriched_terms |= self.o_graph.get_ancestors(go_terms.go_id)
                for t in enriched_terms:
                    terms_counts[t].add(gene)
        return terms_counts
    
    def find_enrichment(self, gene_list, corrections = [], method = _METHOD_TERM_BY_TERM):
        """
        Finds enrichment of specified group of genes.
        
        Parameters
        ----------
        gene_list - list of genes to study
        corrections - list of corrections that should be applied to result.
            Possible values are:
                o "bonferroni" - Bonferroni correction,
                o "fh_fdr" - Benjamin-Hochberg FDR correction.
        method - method of computing the p-values
            Possible values are:
                o "term_by_term" - standard method in which each term is
                    treated independently,
                o "parent_child" - method in which we take the parent-child
                    relation into account when computing p-value.
                    Reference: http://bioinformatics.oxfordjournals.org/content/23/22/3024.long
        """
        
        resolved_list = [self.resolver.resolve(x) for x in gene_list]
        study_counts = self._count_terms(resolved_list)
        
        population_n = len(self.population)
        study_n = len(resolved_list)
        
        result = []
        
        if method == _METHOD_PARENT_CHILD:
            # Calculate enrichment for every term taking parent-child
            # relation into account
            for term, study_set in study_counts.items():
                study_hits = len(study_set)
                population_hits = len(self.population_counts[term])
                fst = True
                parents_study = set()
                parents_population = set()
                # calculate sets of genes annotated to parents
                for parent in self.o_graph.get_parents(term):
                    s_set = study_counts[parent] if parent in study_counts else set()
                    p_set = self.population_counts[parent]
                    if fst:
                        parents_study |= s_set
                        parents_population |= p_set
                        fst = False
                    else:
                        parents_study &= s_set
                        parents_population &= p_set
                        
                pval = Stats.hypergeometric_test(study_hits, len(parents_study),
                                                 population_hits, len(parents_population))
    
                entry = EnrichmentEntry(term, self.o_graph.get_term(term).name,
                                        pval, study_hits, len(parents_study),
                                        population_hits, len(parents_population))
    
                result.append(entry)
                
        else:
            # Calculate enrichment for every term given the counts (term by term)
            for term, study_set in study_counts.items():
                study_hits = len(study_set)
                population_hits = len(self.population_counts[term])
                pval = Stats.hypergeometric_test(study_hits, study_n,
                                                 population_hits, population_n)
    
                entry = EnrichmentEntry(term, self.o_graph.get_term(term).name,
                                        pval, study_hits, study_n,
                                        population_hits, population_n)
    
                result.append(entry)
        
        # Calculate chosen corrections
        if len(corrections) > 0:
            pvals = [x.p_value for x in result]
            corr_pvals = []
            for c_id in corrections:
                cfun = Stats.corrections[c_id]
                corr_pvals.append((c_id, cfun(pvals)))
            for i in xrange(len(result)):
                result[i].corrections = [(c_id, pv[i]) for c_id, pv in corr_pvals] #TODO: tests
        
        warnings = []
        # check for warnings
        if len(self.o_graph.cycles) > 0:
            warnings.append("Graph contains cycles: " + str(self.o_graph.cycles))
        return Enrichment(method, result, warnings)
    
if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()