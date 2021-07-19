"""
Graph utilities used in the python implementation of the sequential
delta-stepping algorithm.

These are used mainly to build and draw a graph for reference.
"""
from matplotlib import pyplot as plt
import networkx as nx
import pandas as pd
from numpy.random import rand

class graphBuilder:

    def __init__(self, filePath):
        self.G = None

        # This dataframe should have two cols: source and sink 
        # node ids.
        self.inOutMatrix = pd.read_csv(filePath, header=None).values
        
        self.nodes_list = self.initNodes()
        self.edges_list = self.initEdges()
    
    def initNodes(self):
        return set(list(self.inOutMatrix[:, 0]) + 
                   list(self.inOutMatrix[:, 1]))

    def initEdges(self):
        edgeWeight = 0x1000
        return list(zip(self.inOutMatrix[:, 0], 
                        self.inOutMatrix[:, 1],
                        [edgeWeight * v for v in
                         range(len(self.nodes_list))]))

    def buildGraphFromFile(self):
        self.G = nx.Graph()

        self.G.add_nodes_from(self.nodes_list)
        self.G.add_weighted_edges_from(self.edges_list)

        return self.G

    def drawGraph(self):
        pos = nx.drawing.nx_agraph.graphviz_layout(self.G, 
                                                   prog="dot")

        # nodes
        nx.draw_networkx_nodes(self.G, pos, node_color="k")

        # edges
        nx.draw_networkx_edges(self.G, pos, width=1.0, alpha=0.5)
        
        # labels
        nx.draw_networkx_labels(self.G, pos, font_color="white")
        edge_attributes = nx.get_edge_attributes(self.G, "weight")
        edge_labels = {e: round(w, 2) for (e, w) in
                       zip(edge_attributes.keys(), 
                           edge_attributes.values())}

        nx.draw_networkx_edge_labels(self.G, pos, 
                                     edge_labels=edge_labels)
