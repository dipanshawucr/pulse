# class SplicingGraphDB()


class SplicingGraph:
    """
    Represents a splicing graph for a particular gene.
    """
    def __init__(self, gene_id):
        """ (str) -> SplicingGraph

        Initializes a SplicingGraph object for a particular gene, and an empty list of transcripts contained.

        :param gene_id:
        :return:
        """
        self.__splicing_graph_dict = {}
        self.gene_id = gene_id
        self.transcripts_list = []

    def nodes(self):
        """
        Returns all of the nodes in the graph

        :return:
        """
        return list(self.__splicing_graph_dict.keys())

    def edges(self):
        """
        Returns all of the edges in the graph
        :return:
        """
        return self.__generate_exon_edges()

    def add_node(self, node):
        """
        If the node is not in self.__graph_dict, then the node with an empty list as its value
        is added to the dictionary. Else, nothing done.
        """
        if node not in self.__splicing_graph_dict:
            self.__splicing_graph_dict[node] = []

    def add_edge(self, edge):
        """
        Adds an edge to self.__graph_dict. There can be multiple edges between two nodes.
        """
        (node1, node2) = (edge.node_from, edge.node_to)
        if node1 in self.__splicing_graph_dict:
            self.__splicing_graph_dict[node1].append(node2)
        else:
            self.__splicing_graph_dict[node1] = [node2]

    def __generate_exon_edges(self):
        """
        Generates all of the edges in the splicing graph.
        """
        edges = []
        for node1 in self.__splicing_graph_dict:
            for node2 in self.__splicing_graph_dict[node1]:
                if {node2, node1} not in edges:
                    edges.append(Edge(True, node1, node2))
        return edges

    # def __str__(self):
    #     """ (self) -> str
    #
    #     Returns a string representation of the Splicing Graph object.
    #
    #     :return:
    #     """
    #     return "gene_id: {}, transcripts_list: {}".format(self.gene_id, self.transcripts_list)

    def __str__(self):
        res = "nodes: " + "\n"
        for k in self.__splicing_graph_dict:
            res += str(k) + " " + "\n"
        res += "\nexon edges: " + "\n"
        for edge in self.__generate_exon_edges():
            res += str(edge) + " " + "\n"
        return res

    # def plot_graph(self):
    #     """
    #     Draws a graph.
    #
    #     :return:
    #     """
    #     G = nx.Graph()
    #     for node in self.nodes():
    #         G.add_node(node)
    #     for edge in self.edges():
    #         G.add_edge(edge)
    #     nx.draw(G)
    #     plt.show()


class SplicingSite:
    """
    Represents a splicing site at a given gene.
    """
    # TODO: Need to specify whether the strand is forward or reverse (+ or -)
    def __init__(self, position, transcript_id, site_type):
        """ (int, str, str) -> SplicingSite

        Initializes a Splicing Site node for the Splicing Graph.

        :param position:
        :param transcript_id:
        :param site_type:
        :return:
        """
        # self.number = number
        self.position = position
        self.transcript_ids = [transcript_id]
        self.site_type = site_type

    def add_transcripts(self, transcript_id):
        """
        Allows to keep track of all of the trancripts that have passed through that node.

        :param transcript_id:
        :return:
        """

        self.transcript_ids.append(transcript_id)

    def __repr__(self):
        return "position: {}, transcript_id: {}, site_type: {}".format(self.position, self.transcript_ids, self.site_type)


class Edge:
    """
    Represents an intron or an exon at a point in a given gene.
    """
    def __init__(self, is_exon, from_node, to_node):
        """ (Edge, bool, SplicingSite, SplicingSite) => NoneType

        Initializes an edge in the splicing graph.

        :param is_exon:
        :param node1:
        :param node2:
        :return:
        """
        self.is_exon = is_exon
        self.node_from = from_node
        self.node_to = to_node

    def __repr__(self):
        return "is_exon: {}, node_from: {}, to_node: {}".format(self.is_exon, self.node_from, self.node_to)