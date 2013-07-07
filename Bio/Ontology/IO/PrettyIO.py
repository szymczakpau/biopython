# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

from Bio.Ontology.IO.GraphIO import GmlWriter
from Bio.Ontology.Graph import DiGraph

def rgb_to_triple(rgb):
    """
    Returns triple of ints from rgb color in format #xxxxxx.
    
    >>> rgb_to_triple("#00bd28")
    (0, 189, 40)
    
    """
    if rgb[0] == '#':
        return (int(rgb[1:3], 16), int(rgb[3:5], 16), int(rgb[5:8], 16))
    else:
        raise ValueError("Not an rgb value.")

def triple_to_rgb(triple):
    """
    Returns rgb color in format #xxxxxx from triple of ints
    
    >>> triple_to_rgb((0, 189, 40))
    '#0bd28'
    
    """
    r, g, b = triple
    return "#{0}{1}{2}".format(hex(r)[2:], hex(g)[2:], hex(b)[2:])

def get_gradient(color_a, color_b, k):
    """
    Returns gradient of colors from a to b with k steps
    
    >>> get_gradient("#aabbcc", "#001122", 10)
    ['#aabbcc', '#99aabb', '#8899aa', '#778899', '#667788', '#556677', '#445566', '#334455', '#223344', '#112233']
    
    """
    r1, g1, b1 = rgb_to_triple(color_a)
    r2, g2, b2 = rgb_to_triple(color_b)
    d = float(k)
    r2 = (r2 - r1) / d
    g2 = (g2 - g1) / d
    b2 = (b2 - b1) / d
    grad = []
    for i in xrange(k):
        grad.append(triple_to_rgb((int(r1), int(g1), int(b1))))
        r1 += r2
        g1 += g2
        b1 += b2
    return grad

class GmlPrinter(object):
    """
    Stores GOGraph as graph visualization in gml format.
    """
    
    def __init__(self, file_handle, params = {"step" : 10,
                                              "color_a" : "#00bd28",
                                              "color_b" : "#f7ffc8"}):
        self.handle = file_handle
        self.params = params
        self.gradient_step = params["step"]
        self.gradient = get_gradient(params["color_a"], params["color_b"],
                                     self.gradient_step)
        
    def to_printable_data(self, e_entry):
        return {"name" : e_entry.name,
                "pvalue" : e_entry.p_value,
                "graphics" : { "fill" : self.gradient[int(e_entry.p_value * (self.gradient_step - 1))]
                              }
                }
    def entry_to_label(self, entry):
        return str(entry.oid)
    
    def to_printable_graph(self, enrichment, graph):
        viz_graph = DiGraph()
        viz_graph.attrs["defaultnodesize"] = "labelsize"
        viz_graph.attrs["label"] = str(enrichment)
        
        entry_labels = {}
        
        for entry in enrichment.entries:
            new_label = self.entry_to_label(entry)
            viz_graph.add_node(new_label, self.to_printable_data(entry))
            entry_labels[entry.oid] = new_label
            
        for entry in enrichment.entries:
            u = graph.get_node(entry.oid)
            for edge in u.succ:
                viz_graph.add_edge(entry_labels[entry.oid], entry_labels[edge.to_node.label])
        
        return viz_graph
    
    def pretty_print(self, enrichment, graph):
        nodes_list = [x.oid for x in enrichment.entries]
        g = graph.get_induced_subgraph(nodes_list)
        
        vg = self.to_printable_graph(enrichment, g)
        GmlWriter(self.handle).write_file(vg)

class GraphVizPrinter(object):
    """
    Stores GOGraph as visualization in png format using graphviz library.
    """
    
    def __init__(self, file_handle, params = {"dpi" : 96,
                                              "step" : 10,
                                              "color_a" : "#00bd28",
                                              "color_b" : "#f7ffc8"}):
        self.handle = file_handle
        self.params = params
        self.gradient_step = params["step"]
        self.gradient = get_gradient(params["color_a"], params["color_b"],
                                     self.gradient_step)
    
    def entry_to_label(self, entry):
        return "{0}\n{1}\np:{2:.2f}".format(entry.oid, entry.name, entry.p_value)

    def to_printable_graph(self, enrichment, graph):
        import pygraphviz
        
        viz_graph = pygraphviz.AGraph()
        viz_graph.graph_attr.update(dpi = str(self.params["dpi"]))
        viz_graph.node_attr.update(shape="box", style="rounded,filled")
        viz_graph.edge_attr.update(shape="normal", color="black", dir="back", label="is_a")
        
        entry_labels = {}
        
        for entry in enrichment.entries:
            new_label = self.entry_to_label(entry)
            col = self.gradient[int(entry.p_value * (self.gradient_step - 1))]
            viz_graph.add_node(new_label, fillcolor = col)
            entry_labels[entry.oid] = new_label
            
        for entry in enrichment.entries:
            u = graph.get_node(entry.oid)
            for edge in u.succ:
                viz_graph.add_edge(entry_labels[edge.to_node.label], entry_labels[entry.oid])
                
        return viz_graph
            
            
        
    def pretty_print(self, enrichment, graph):
        nodes_list = [x.oid for x in enrichment.entries]
        g = graph.get_induced_subgraph(nodes_list)
        
        vg = self.to_printable_graph(enrichment, g)
        vg.draw(self.handle, prog="dot")
    
if __name__ == '__main__':
    from Bio._utils import run_doctest
    run_doctest()