# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

from Bio.Ontology.IO.GraphIO import GmlWriter
from Bio.Ontology.Graph import DiGraph
from Bio.Ontology.Stats import corrections_labels

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
    Stores found enrichments as a graph in gml format.
    """
    
    def __init__(self, file_handle, params = {"step" : 10,
                                              "color_a" : "#00bd28",
                                              "color_b" : "#f7ffc8",
                                              "color_none" : "#c3c3c3"}):
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
    def term_to_printable(self, term):
        return {"name" : term.name,
                "graphics" : { "fill" : self.params["color_none"]
                              }
                }
    def entry_to_label(self, entry):
        return str(entry.oid)
    
    def term_to_label(self, term):
        return str(term.id)
    
    def to_printable_graph(self, enrichment, graph):
        viz_graph = DiGraph()
        viz_graph.attrs["defaultnodesize"] = "labelsize"
        viz_graph.attrs["label"] = str(enrichment)
        
        entry_labels = {}
        
        for entry in enrichment.entries:
            new_label = self.entry_to_label(entry)
            viz_graph.add_node(new_label, self.to_printable_data(entry))
            entry_labels[entry.oid] = new_label
        
        for label, node in graph.nodes.items():
            if label not in entry_labels:
                new_label = self.term_to_label(node.data)
                viz_graph.add_node(new_label, self.term_to_printable(node.data))
                entry_labels[label] = new_label
        
        for label, u in graph.nodes.items():
            for edge in u.succ:
                viz_graph.add_edge(entry_labels[label], entry_labels[edge.to_node.label])
                
        
#         for entry in enrichment.entries:
#             new_label = self.entry_to_label(entry)
#             viz_graph.add_node(new_label, self.to_printable_data(entry))
#             entry_labels[entry.oid] = new_label
#             
#         for entry in enrichment.entries:
#             u = graph.get_node(entry.oid)
#             for edge in u.succ:
#                 viz_graph.add_edge(entry_labels[entry.oid], entry_labels[edge.to_node.label])
        
        return viz_graph
    
    def pretty_print(self, enrichment, graph):
        nodes_ids = set()
        for x in enrichment.entries:
            nodes_ids.add(x.oid)
            nodes_ids = nodes_ids.union(graph.get_ancestors(x.oid))
        
        g = graph.get_induced_subgraph(nodes_ids)
        
        vg = self.to_printable_graph(enrichment, g)
        GmlWriter(self.handle).write_file(vg)

class GraphVizPrinter(object):
    """
    Stores found enrichments as visualization in png format using graphviz library.
    """
    
    def __init__(self, file_handle, params = {"dpi" : 96,
                                              "step" : 10,
                                              "color_a" : "#00bd28",
                                              "color_b" : "#f7ffc8",
                                              "color_none" : "#c3c3c3"}):
        self.handle = file_handle
        self.params = params
        self.gradient_step = params["step"]
        self.gradient = get_gradient(params["color_a"], params["color_b"],
                                     self.gradient_step)
    
    def entry_to_label(self, entry):
        return "{0}\n{1}\np:{2}".format(entry.oid, entry.name, entry.p_value)

    def term_to_label(self, term):
        return "{0}\n{1}".format(term.id, term.name)
    
    def to_printable_graph(self, enrichment, graph):
        import pygraphviz #TODO exception + warning
        
        viz_graph = pygraphviz.AGraph()
        viz_graph.graph_attr.update(dpi = str(self.params["dpi"]))
        viz_graph.node_attr.update(shape="box", style="rounded,filled")
        viz_graph.edge_attr.update(shape="normal", color="black", dir="back")
        
        entry_labels = {}
        
        for entry in enrichment.entries:
            new_label = self.entry_to_label(entry)
            col = self.gradient[int(entry.p_value * (self.gradient_step - 1))]
            viz_graph.add_node(new_label, fillcolor = col)
            entry_labels[entry.oid] = new_label
        
        col = self.params["color_none"]
        for label, node in graph.nodes.items():
            if label not in entry_labels:
                new_label = self.term_to_label(node.data)
                viz_graph.add_node(new_label, fillcolor = col)
                entry_labels[label] = new_label
        
        for label, u in graph.nodes.items():
            for edge in u.succ:
                viz_graph.add_edge(entry_labels[edge.to_node.label], entry_labels[label],  label=edge.data)
                
        return viz_graph
            
            
        
    def pretty_print(self, enrichment, graph):
        nodes_ids = set()
        for x in enrichment.entries:
            nodes_ids.add(x.oid)
            nodes_ids = nodes_ids.union(graph.get_ancestors(x.oid))
        g = graph.get_induced_subgraph(nodes_ids)
        
        vg = self.to_printable_graph(enrichment, g)
        vg.draw(self.handle, prog="dot")

class TxtPrinter(object):
    """
    Prints found enrichments to txt file.
    """
    
    def __init__(self, file_handle, params = None):
        self.handle = file_handle
        self.params = params
        
    def pretty_print(self, enrichment, graph):
        self.handle.write("Enrichments found using {0} method.\n\nEnrichments:\n\n"
                          .format(enrichment.method))
        sorted_entries = sorted(enrichment.entries, key = lambda x: x.p_value)
        for x in sorted_entries:
            self.handle.write(str(x))
            self.handle.write("\n\n")
        if (len(enrichment.warnings) > 0):
            self.handle.write("Warnings:\n")
            for x in enrichment.warnings:
                self.handle.write(str(x))


class HtmlPrinter(object):
    """
    Prints found enrichments to html file.
    """
    
    def __init__(self, file_handle, params = {"go_to_url" : "http://amigo.geneontology.org/cgi-bin/amigo/term_details?term="}):
        self.handle = file_handle
        self.params = params
        self.style = """<style type="text/css">
.warning
{
    color:#D13F31;
}
body
{
    margin:45px;
}
h1
{
    color:#669;
}
table
{
    font-family: Sans-Serif;
    font-size: 14px;
    background: #fff;
    margin: 45px 0px;
    border-collapse: collapse;
    text-align: left;
}
th
{
    font-size: 16px;
    font-weight: normal;
    color: #009;
    padding: 12px 10px;
    border-bottom: 2px solid #6678b0;
}
td
{
    color: #669;
    padding: 8px 10px;
    border-bottom: 1px solid #ccc;
}
tbody tr:hover td
{
    color: #009;
}
</style>
"""
    
    
    def write_tag(self, tag, text, attrs = None):
        self.open_tag(tag, attrs)
        self.handle.write(text)
        self.close_tag(tag)
    
    def open_tag(self, tag, attrs = None):
        ot = "<" + tag
        if attrs != None:
            for k, v in attrs.items():
                ot += ' {0}="{1}"'.format(k,v)
        ot += ">\n"
        self.handle.write(ot)
        
    def close_tag(self, tag):
        self.handle.write("</" + tag + ">\n")
    
    def pretty_print(self, enrichment, graph):
        self.handle.write("<!DOCTYPE html>")
        
        self.open_tag("html")
        self.open_tag("head")
        self.handle.write(self.style)
        self.close_tag("head")
        self.open_tag("body")
        self.write_tag("h1", "Enrichments found using {0} method."
                          .format(enrichment.method))
        
        sorted_entries = sorted(enrichment.entries, key = lambda x: x.p_value)
        
        self.open_tag("table")
        
        self.open_tag("tr")        
        for header in ["ID", "name", "p-value"] + [corrections_labels[x] for x in enrichment.corrections]:
            self.write_tag("th", header)
        self.close_tag("tr")
        
        for x in sorted_entries:
            self.open_tag("tr")
            self.open_tag("td")
            self.write_tag("a", str(x.oid), {"href" : self.params["go_to_url"] + str(x.oid)})
            self.close_tag("td")
            self.write_tag("td", str(x.name))
            self.write_tag("td", str(x.p_value))
            for _, corr in x.corrections:
                self.write_tag("td", str(corr))
            self.close_tag("tr")
            
        self.close_tag("table")
        
        if (len(enrichment.warnings) > 0):
            self.write_tag("h1", "Warnings:", {"class" : "warning"})
            self.open_tag("ul", {"class" : "warning"})
            for x in enrichment.warnings:
                self.write_tag("li", str(x))
            self.close_tag("ul")
        
        self.close_tag("body")
        self.close_tag("html")

if __name__ == '__main__':
    from Bio._utils import run_doctest
    run_doctest()