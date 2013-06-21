# Copyright 2013 by Kamil Koziara. All rights reserved
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


from Bio.Ontology.Graph import DiGraph

class GOGraph(DiGraph):
    """
    Represents Gene Ontology graph.
    """
    
    _ANCESTORS = "ancestors"
    
    def __init__(self, terms):
        DiGraph.__init__(self)
        for (term_type, data) in terms:
            if term_type == "Term": # Add only terms for now
                nid = data.pop("id")[0]
                name = data.pop("name")[0]
                term = GOTerm(nid, name, data)
                if self.node_exists(nid):
                    self.update_node(nid, term)
                else:
                    self.add_node(nid, term)
                for edge in data["is_a"]:
                    self.add_edge(nid, edge, "is_a")
                for edge in data["relationship"]: #TODO: parse relationship better 
                    p = edge.find("part_of")
                    if p >= 0:
                        r_edge = edge[p + 7:].strip()
                        self.add_edge(nid, r_edge, "part_of")
                        
    def get_term(self, go_id):
        return self.get_node(go_id).data
    
    def get_ancestors(self, go_id):
        node = self.get_node(go_id)
        return self._get_ancestors(node)
    
    def _get_ancestors(self, node):
        if GOGraph._ANCESTORS in node.attr:
            return node.attr[GOGraph._ANCESTORS]
        else:
            anc_set = set()
            for edge in node.succ:
                anc_set |= self._get_ancestors(edge.to_node)
                anc_set.add(edge.to_node.label)
            node.attr[GOGraph._ANCESTORS] = anc_set
            return anc_set    


class GOTerm(object):
    
    def __init__(self, nid, name, attrs):
        self.id = nid
        self.name = name
        self.attrs = attrs
        
    def __str__(self):
        s = self.name + "\n" + "id: " + self.id + "\n"
        for k, v in self.attrs.items():
            s += "{0} : {1}\n".format(k, v)
            
    def __repr__(self):
        return "GOTerm(id = " + self.id + ", name = " + self.name + ")" 

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
