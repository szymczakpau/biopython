# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

from Bio.Ontology.Data import OntologyGraph, OntologyTerm, TermAssociation, GeneAnnotation
import xml.sax
import collections
import ast

_SKIP = 0
_NODE = 1
_EDGE = 2

class NexoContentHandler(xml.sax.ContentHandler):
    def __init__(self):
        xml.sax.ContentHandler.__init__(self)
        
        self.state = _SKIP
        self.old_state = _SKIP
        
        self.nodes = {}
        self.annotations = collections.defaultdict(list)
        self.edges = []
        
        self.current_term = None
        self.current_edge = None
    
    def _split_list(self, val):
        if val.startswith('['):
            return ast.literal_eval(val)
        else:
            return [val]
        
    def startElement(self, name, attrs):
        if name == "node":
            self.state = _NODE
            term_id = attrs["label"]
            self.current_term = OntologyTerm(term_id, term_id, {})
            self.nodes[attrs["id"]] = self.current_term

        elif name == "edge":
            self.state = _EDGE
            self.current_edge = [attrs["source"], attrs["target"]]
            self.edges.append(self.current_edge)
            
        elif name == "graphics":
            self.old_state = self.state
            self.state = _SKIP
            
        elif name == "att": #TODO add or no rest of the values ?
            if self.state == _NODE:
                if attrs["name"] == "Term":
                    val = attrs.get("value")
                    if val != None:
                        self.current_term.name = val
                elif attrs["name"] == "Assigned Genes":
                    val = attrs.get("value")
                    if val != None:
                        for gene in self._split_list(val):
                            self.annotations[gene].append(self.current_term.id)
            elif self.state == _EDGE:
                if attrs["name"] == "NeXO relation type":
                    self.current_edge.append(attrs.get("value"))
                    
    def endElement(self, name):
        if name == "node" or name == "edge":
            self.state = _SKIP
        elif name == "graphics":
            self.state = self.old_state
 
    def characters(self, content):
        pass
    
class NexoReader(object):
    """
    Class for reading Nexo xgmml network.
    """


    def __init__(self, file_handle):
        self.handle = file_handle
        
    
    def read(self):
        """
        Returns gene annotation list and ontology graph read from nexo file.
        """
        
        content_handler = NexoContentHandler()
        xml.sax.parse(self.handle, content_handler)
        
        annotations = []
        for obj, assocs in content_handler.annotations.iteritems():
            annotations.append(GeneAnnotation(obj,
                        associations = [TermAssociation(x) for x in assocs]))
        
        graph = OntologyGraph()

        for _, node in content_handler.nodes.iteritems():
            graph.add_node(node.id, node)
        
        for edge in content_handler.edges:
            source = content_handler.nodes[edge[0]].id
            target = content_handler.nodes[edge[1]].id
            graph.add_edge(source, target, edge[2])
        
        return (annotations, graph)
        