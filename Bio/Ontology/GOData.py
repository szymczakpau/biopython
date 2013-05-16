# Copyright 2013 by Kamil Koziara. All rights reserved
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


from Bio.Ontology.Graph import DiGraph

class GOGraph(DiGraph):
    """
    Represents Gene Ontology graph.
    """
    
    def __init__(self, terms):
        DiGraph.__init__(self)
        for (term_type, data) in terms:
            if term_type == "Term": # Add only terms for now
                nid = data["id"][0]
                self.set_node(nid, data)
                for edge in data["is_a"]:
                    self.add_edge(nid, edge)

    def get_term(self, go_id):
        return self.get_node(go_id).data


class GOAObject(object):
    """
    Represents one gene ontology association object
    """
    
    def __init__(self, db, oid, symbol, name, otype, taxon, ext = None, gp_id = None, synonyms = None, associations = None):
        self.db = db
        self.oid = oid
        self.symbol = symbol
        self.name = name
        self.otype = otype
        self.taxon = taxon
        
        self.ext = ext
        self.gp_id = gp_id
        
        self.synonyms = synonyms
        self.associations = associations
    
    
class GOAssociation(object):
    """
    Represents one gene ontology association
    """
    
    def __init__(self, qualifier, go_id, db_ref, evidence, wf, aspect, date, assigned_by):
        self.qualifier = qualifier
        self.go_id = go_id
        self.db_ref = db_ref
        self.evidence = evidence
        self.wf = wf
        self.aspect = aspect
        self.date = date
        self.assigned_by = assigned_by

        