# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

import collections
import Stats, IdResolver
import random, math

class EnrichmentEntry(object):
    """
    Represents one result returned by EnrichmentFinder.
    """
    def __init__(self, oid, name, p_value):
        self.oid = oid
        self.name = name
        self.p_value = p_value
        self.corrections = []
        self.attrs = {}
    
    def __repr__(self):
        return 'EnrichmentEntry([("ID" : {0}), ("name" : {2}), ("p-value" : {1})'.format(self.oid, self.p_value, self.name)

    def __str__(self):
        return """ID : {0}
name : {1}
p-value : {2}
corrected p-values: {3}""".format(self.oid, self.name, self.p_value, self.corrections)

class Enrichment(object):
    """
    Contains all results found by EnrichmentFinder
    """
    
    def __init__(self, method, entries, warnings, corrections = []):
        self.entries = entries
        self.warnings = warnings
        self.method = method
        self.corrections = corrections
    
    def filter(self, filter_fun):
        return Enrichment(self.method, filter(filter_fun, self.entries),
                          list(self.warnings), list(self.corrections))
        
    def get_significant(self, p_val):
        """
        Returns enrichment entries with specified siginificance level.
        """
        return self.filter(lambda x: x.p_value <= p_val)
    
    def __repr__(self):
        return "Enrichment(method = {2}, entries_num = {0} , warnings_num = {1})".format(len(self.entries), len(self.warnings), self.method)

    def __str__(self):
        return "Enrichment found using {2} method: {0} entries, {1} warnings.".format(len(self.entries),
                                                               len(self.warnings), self.method)
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
    
    >>> from Bio.Ontology.Data import OntologyGraph
    >>> go_graph = OntologyGraph(gonodes_iter)
    
    Now is the time to create EnrichmentFinder. Besides the arguments
    mentioned before you could specify an id resolver, which basically tries
    to find synonyms for gene ids that finder can understand, and population
    of genes as reference. If none of this is specified default resolver
    is used and all genes from association are used as the population.
    
    >>> from Bio.Ontology import EnrichmentFinder
    >>> ef = EnrichmentFinder(assocs_iter, go_graph)
    
    To run finder you just need to call find_enrichment method with list
    of genes as an argument:
    
    >>> genes_to_study = ['FBgn0070057', '18-wheeler', 'FBgn0043467']
    
    Additionally you can add a list of corrections which should be used
    (by default it's empty):
    
    >>> corrections = ['bonferroni', 'bh_fdr']
    >>> result = ef.find_enrichment(genes_to_study, corrections)
    
    The result contains a list of EnrichmentEntry instances - each for a node in graph
    which is enriched - and a list of warnings.
    
    >>> print result
    Enrichment found using term_by_term method: 64 entries, 1 warnings.
    
    Notice there is one warning. Let's print the list:
    >>> print result.warnings
    ["Unknown id: '18-wheeler' was resolved to: 'FBgn0004364'"]
    
    Earlier when specifying list of genes we used non-standardized gene id:
    '18-wheeler'. Thanks to default resolver (FirstOneResolver)
    of EnrichmentFinder standard id was inferred.
    
    >>> print result.entries[0]
    ID : GO:0044707
    name : single-multicellular organism process
    p-value : 1.0
    corrected p-values: [('bonferroni', 1.0), ('bh_fdr', 1.0)]
    
    You can also specify a method of computing p-values. Default method is
    term by term, but additionaly you can choose method that takes parent-child
    relationship into account:
    
    >>> result = ef.find_enrichment(genes_to_study, corrections, "parent_child_intersection")
    
    >>> print result
    Enrichment found using parent_child_intersection method: 64 entries, 1 warnings.
    >>> print result.entries[0]
    ID : GO:0044707
    name : single-multicellular organism process
    p-value : 1.0
    corrected p-values: [('bonferroni', 1.0), ('bh_fdr', 1.0)]
    
    """
    
    def __init__(self, annotations, o_graph, population = None,
                 resolver_generator = IdResolver.FirstOneResolver):
        """
        Initialize Enrichment Finder.
        
        Parameters
        ----------
        annotations - iterable containing annotations
        o_graph - graph with ontology
        population - population used as reference to study sample
        resolver_generator - constructor of resolver used to disambiguate ids
        """
        self.o_graph = o_graph
        self.annotations = dict([(v.oid, v) for v in annotations])
        
        self.resolver = resolver_generator(self.annotations.itervalues())
        if population != None:
            self.population = population
        else:
            self.population = self.annotations.keys()
        self.terms_to_population_genes = self._find_terms_associations(self.population)
        
    def _find_terms_associations(self, gene_list):
        terms_assocs = collections.defaultdict(set)
        for gene in gene_list:
            if gene in self.annotations:
                enriched_terms = set()
                for term in self.annotations[gene].associations: # qualifier co z nim?
                    node = self.o_graph.get_node(term.go_id)
                    if node != None:
                        nid = node.data.id # because enrichments may use synonyms instead of ids
                        enriched_terms.add(nid)
                        enriched_terms |= self.o_graph.get_ancestors(nid)
                for t in enriched_terms:
                    terms_assocs[t].add(gene)
        return terms_assocs
    
    def _find_term_by_term_enrichment(self, terms_to_study_genes, study_size):
        result = []
        population_size = len(self.population)
        # Calculate enrichment for every term given the counts (term by term)
        for term, study_set in terms_to_study_genes.items():
            study_hits = len(study_set)
            population_hits = len(self.terms_to_population_genes[term])
            pval = Stats.hypergeometric_test(study_hits, study_size,
                                             population_hits, population_size)

            entry = EnrichmentEntry(term, self.o_graph.get_term(term).name, pval)

            result.append(entry)
        
        return result
    
    def _count_op_items(self, set_list, set_op):
        return len(set_op(*set_list)) if len(set_list) > 0 else 0
    
    def _find_parent_child_enrichment(self, terms_to_study_genes, method):
        # Calculate enrichment for every term taking parent-child
        # relation into account
        result = []
        for term, study_set in terms_to_study_genes.items():
            study_hits = len(study_set)
            population_hits = len(self.terms_to_population_genes[term])
            study_set_list = []
            pop_set_list = []
            # calculate sets of genes annotated to parents
            for parent in self.o_graph.get_parents(term):
                study_set_list.append(terms_to_study_genes[parent])
                pop_set_list.append(self.terms_to_population_genes[parent])
            
            if method == "parent_child_union":
                set_op = set.union
            elif method == "parent_child_intersection":
                set_op = set.intersection
            else:
                raise ValueError("{0} is not correct method type.".format(method))
            
            parents_study_size = self._count_op_items(study_set_list, set_op)
            population_parents_size = self._count_op_items(pop_set_list, set_op)
            
            if study_hits <= parents_study_size and population_hits <= population_parents_size:
                pval = Stats.hypergeometric_test(study_hits, parents_study_size,
                                             population_hits, population_parents_size)

                entry = EnrichmentEntry(term, self.o_graph.get_term(term).name, pval)
                entry.attrs = {"study_hits" : study_hits, "parents_study_size": parents_study_size,
                            "population_hits" : population_hits, "population_parents_size" : population_parents_size}
                result.append(entry)
        
        return result
    
    def _resolve_ids(self, gene_list, warnings):
        resolved_list = []
        for x in gene_list:
            rx = self.resolver.resolve(x)
            if x != rx:
                warnings.append("Unknown id: '{0}' was resolved to: '{1}'".format(x, rx))
            resolved_list.append(rx)
        return resolved_list
        
    def find_enrichment(self, gene_list, corrections = [], method = "term_by_term"):
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
                o "parent_child_union", "parent_child_intersection" - methods
                    in which we take the parent-child relationship into account
                    when computing p-value.
                    Reference: http://bioinformatics.oxfordjournals.org/content/23/22/3024.long
        """
        
        result = []
        warnings = []
        
        resolved_list = self._resolve_ids(gene_list, warnings)
            
        terms_to_study_genes = self._find_terms_associations(resolved_list)
        
        
        study_size = len(resolved_list)

        if method == "term_by_term":
            result = self._find_term_by_term_enrichment(terms_to_study_genes,
                                                        study_size)
        else:
            result = self._find_parent_child_enrichment(terms_to_study_genes,
                                                        method)
        
        # Calculate chosen corrections
        if len(corrections) > 0:
            pvals = [x.p_value for x in result]
            corr_pvals = []
            for c_id in corrections:
                cfun = Stats.corrections[c_id]
                corr_pvals.append((c_id, cfun(pvals)))
            for i in xrange(len(result)):
                result[i].corrections = [(c_id, pv[i]) for c_id, pv in corr_pvals]
        
        # check for warnings
        if len(self.o_graph.cycles) > 0:
            warnings.append("Graph contains cycles: " + str(self.o_graph.cycles))
        return Enrichment(method, result, warnings, corrections)
    
    def _find_enriched_terms(self, gene_list):
        enriched_terms = set()
        for gene in gene_list:
            if gene in self.annotations:
                for assoc in self.annotations[gene].associations:
                    enriched_terms.add(assoc.go_id)
        return enriched_terms
    
    def _get_perms(self, gene_list, perms_no):
        perms = []
        permutation = list(gene_list)
        for i in xrange(perms_no):
            random.shuffle(permutation)
            perms.append(list(permutation))
        return perms
    
    def _kolmogorov_smirnov_rank_test(self, gene_set, gene_list, gene_corr, p):
        cval = 0
        Dn = 0
        Nr = 0
        adj_corr = [math.pow(x, p) for x in gene_corr]
        for i in xrange(len(gene_list)):
            if gene_list[i] in gene_set:
                Nr += adj_corr[i] 
        
        miss_pen = float(1) / (len(gene_list) - len(gene_set))
         
        for i in xrange(len(gene_list)):
            if gene_list[i] in gene_set:
                cval += adj_corr[i] / Nr
            else:
                cval -= miss_pen
            if abs(cval) > abs(Dn):
                Dn = cval
        return Dn
    
    def find_enrichment_from_rank(self, gene_rank, perms_no = 100): #TODO test!!!
        sorted_gene_rank = sorted(gene_rank, key = lambda x: x[1])
        gene_corr = []
        gene_list = []
        
        for pos in sorted_gene_rank:
            gene_list.append(pos[0])
            gene_corr.append(pos[1])
        
        warnings = []
        result = []
        resolved_list = self._resolve_ids(gene_list, warnings)
        enriched_terms = self._find_enriched_terms(resolved_list)
        perms = self._get_perms(resolved_list, perms_no)
        
        for term in enriched_terms:
            gene_set = self.terms_to_population_genes[term]
            orig_dn = self._kolmogorov_smirnov_rank_test(self, gene_set, resolved_list, gene_corr, 1)
            
            pcount = 0
            for perm in perms:
                perm_dn = self._kolmogorov_smirnov_rank_test(self, gene_set, perm, gene_corr, 1)
                if orig_dn < 0:
                    pcount += int(orig_dn > perm_dn)
                else:
                    pcount += int(orig_dn < perm_dn)
            pval = pcount / float(perms_no)
            result.append(EnrichmentEntry(term, self.o_graph.get_term(term).name,
                                    pval))
        
        if len(self.o_graph.cycles) > 0:
            warnings.append("Graph contains cycles: " + str(self.o_graph.cycles))
        return Enrichment("GSEA", result, warnings)

class RankedEnrichmentFinder(object):
    
    def __init__(self, annotations, o_graph,
                 resolver_generator = IdResolver.FirstOneResolver):
        """
        Initialize Enrichment Finder.
        
        Parameters
        ----------
        annotations - iterable containing annotations
        o_graph - graph with ontology
        resolver_generator - constructor of resolver used to disambiguate ids
        """
        
        self.o_graph = o_graph
        self.annotations = dict([(v.oid, v) for v in annotations])
        
        self.resolver = resolver_generator(self.annotations.itervalues())
    
    def _resolve_ids(self, gene_list, warnings):
        resolved_list = []
        for x in gene_list:
            rx = self.resolver.resolve(x)
            if x != rx:
                warnings.append("Unknown id: '{0}' was resolved to: '{1}'".format(x, rx))
            resolved_list.append(rx)
        return resolved_list
    
    
    def _get_half_results(self, resolved_list, ef, warnings):
        list_slice = []
        results = collections.defaultdict(list)

        for gene in resolved_list:
            list_slice.append(gene)
            slice_res = ef.find_enrichment(list_slice, method="parent_child_intersection")
            warnings += slice_res.warnings
            for e in slice_res.entries:
                results[e.oid].append(e)
        return results
    
    def find_enrichment_parent_child(self, gene_rank):        
        """
        Finds enrichment by applying parent-child analysis to list slices.
        """
        
        warnings = []
        
        sorted_gene_rank = sorted(gene_rank, key = lambda x: x[1])
        gene_list = [g for g, _ in sorted_gene_rank]
        resolved_list = self._resolve_ids(gene_list, warnings)
        
        
        ef = EnrichmentFinder(self.o_graph, self.annotations.itervalues(),
                                  resolved_list, IdResolver.Resolver)
        
        plus_results = self._get_half_results(resolved_list, ef, warnings)
        minus_results = self._get_half_results(resolved_list[::-1], ef, warnings)
        
        return (plus_results, minus_results)        
        
if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()