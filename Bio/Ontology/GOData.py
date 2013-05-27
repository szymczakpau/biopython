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
    
    def __repr__(self):
        return "GOAObject(db_object_id = {0})".format(self.oid)
    
    def __str__(self):
        b1 = """DB: {0}
DB Object ID: {1}
DB Object Symbol: {2}
DB Object Name: {3}
DB Object Synonyms: {4}
DB Object Type: {5}
Taxon: {6}
""".format(self.db, self.oid, self.symbol, self.name, self.synonyms, self.otype, self.taxon)
        if self.ext != None: # only in gaf 2.0
            b1 += """Annotation Extension: {0}
Gene Product Form ID: {1}
""".format(self.ext, self.gp_id)
        if self.associations != None:
            b1 += "Associations:"
            for a in self.associations:
                b1 += str(a)
        return b1
    
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

    def __repr__(self):
        return "GOAssociation(go_id = {0})".format(self.go_id)
    
    def __str__(self):
        return """
    Qualifier: {0}
    GO ID: {1}
    DB:Reference: {2}
    Evidence Code: {3}
    With (or) From: {4}
    Aspect: {5}
    Date: {6}
    Assigned by: {7}
    """.format(self.qualifier, self.go_id, self.db_ref,\
self.evidence, self.wf, self.aspect, self.date, self.assigned_by)
