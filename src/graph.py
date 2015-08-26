__author__ = 'snow'
"""
 utils is a package for graph handing in coldnet: a combinatorial disease network.
"""
import networkx as nx


def build_graph(G, data):
    nodes = set(row[0] for row in data)
    edges = [(row[0], row[1]) for row in data]
    G.add_nodes_from(nodes)
    G.add_edges_from(edges)


class InteractomeGraph:
    """
    Interactome graph implementation
    """

    def __init__(self):
        self.ppi_graph = nx.Graph()
        self.reg_graph = nx.DiGraph()

    def build(self, ppi_edges=None, reg_edges=None):
        """
        Build the graph from edges
        :param ppi_edges: tuple of edges for PPI network
        :param reg_edges: tuple of edges for regulatory network
        :return: None
        """
        if ppi_edges is not None:
            build_graph(self.ppi_graph, ppi_edges)
        if reg_edges is not None:
            build_graph(self.reg_graph, reg_edges)

    def propagate(self, v):
        """
        Finding the personalized pagerank of the random walk
        :param v: set of genes with mutations
        :return: a propagated pagerank of all genes
        """
        personalized = dict.fromkeys(self.ppi_graph.V, 0.0)
        pv = 1.0 / len(v)
        for n in v:
            personalized[n] = pv

        return nx.algorithms.pagerank(G=self.ppi_graph, personalization=personalized, max_iter=100, tol=1.0e-6)

