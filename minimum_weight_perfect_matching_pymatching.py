import os
from layout import Hexagonal_layout
import bidict as bd
import itertools as it
import networkx as nx
import itertools
import copy
from create_matching_graph import Matching_graph
import numpy as np


class MWPM:
    """
    MWPM for hexagonal colour code
    """

    def __init__(self, distance, layout, error_model):

        #        self.initialize_blossom()
        self.layout = layout
        self.corr_to_ind = bd.bidict(
            zip(
                sorted(self.layout.ancilla_qubits),
                range(len(self.layout.ancilla_qubits)),
            )
        )
        g5 = Matching_graph(distance, n_layers=1, error_model=error_model)

        self.coords_matching_graph_X = g5.graphX
        self.matching_graph_X = copy.copy(g5.graphX)

        self.matching_graph_X.green_red = g5.green_red_X
        self.matching_graph_X.blue_green = g5.blue_green_X
        self.matching_graph_X.red_blue = g5.red_blue_X

        self.coords_matching_graph_Z = g5.graphZ
        self.matching_graph_Z = copy.copy(g5.graphZ)

        self.matching_graph_Z.green_red = g5.green_red_Z
        self.matching_graph_Z.blue_green = g5.blue_green_Z
        self.matching_graph_Z.red_blue = g5.red_blue_Z

    def initialize_blossom(self):
        dir_path = os.path.dirname(os.path.realpath(__file__))
        blossom_path = dir_path + "/../blossom5/libblossom.so"

        cdef_str = """
            typedef struct {
                int uid;
                int vid;
                int weight;
            } Edge;
            int Init();
            int Process(int node_num, int edge_num, Edge *edges);
            int PrintMatching();
            int GetMatching(int * matching);
            int Clean();
            """

        self.ffi = FFI()
        self.blossom = self.ffi.dlopen(blossom_path)
        self.ffi.cdef(cdef_str)

    def single_run(self, matching_graph, coords_matching_graph, error_cor=None):

        if error_cor is None:
            error_cor = self.layout.create_error()

        red_green_matching, red_blue_matching, green_blue_matching = self.get_matching(
            error_cor, matching_graph, coords_matching_graph
        )
        correction = set()
        graph = nx.compose_all(
            [red_green_matching, red_blue_matching, green_blue_matching]
        )
        (
            sigma_vertex_groups,
            sigma_graphs,
            rho_vertices,
        ) = self.find_sigma_and_rho_vertices(copy.deepcopy(graph))

        for index, sigma_vertices in enumerate(sigma_vertex_groups):
            for sigma_vertex in sigma_vertices:
                correction = self.find_correction_for_sigma_vertix(
                    sigma_vertex, sigma_graphs[index], correction
                )
        rho_graph = nx.Graph()
        rho_graph.add_edges_from(red_green_matching.edges())
        rho_graph.add_edges_from(red_blue_matching.edges())
        for s_graph in sigma_graphs:
            rho_graph.remove_edges_from(s_graph.edges())

        for node in rho_vertices:
            correction = self.find_correction_for_sigma_vertix(
                node, rho_graph, correction
            )
        return correction

    def find_correction_for_sigma_vertix(self, node, graph, correction):
        """
        finds correction for a vertex representing a connected component
        """
        n = graph.neighbors(node)
        all_neighbors = [neighbor for neighbor in n]
        boundary_neighbors = [
            neighbor for neighbor in all_neighbors if type(neighbor) == str
        ]
        neighbors = [neighbor for neighbor in all_neighbors if type(neighbor) != str]
        neighbors = self.translate_boundary_neighbors(
            boundary_neighbors, node, neighbors
        )
        new_correction = self.lift(node, neighbors)
        correction ^= new_correction
        return correction

    def translate_boundary_neighbors(self, boundary_neighbor, node, neighbors):
        # colour_dict, neighbors):
        """
        Translate boundary nodes to coordinates. At these coordinates there
        is not actually an ancilla.
        """
        if boundary_neighbor == []:
            neighbor = boundary_neighbor
        for neighbor in boundary_neighbor:
            colour = neighbor[1]
            if colour == "R":
                if node in self.layout.ancilla_qubits_blue:

                    if (node[0] - 3, node[1] + 1) in neighbors:
                        neighbors += [(node[0], node[1] + 2)]
                    else:
                        neighbors += [(node[0] + 3, node[1] - 1)]

                else:  # the node must be in green
                    neighbors += [(node[0] + 3, node[1] + 1)]
            elif colour == "G":
                if node in self.layout.ancilla_qubits_red:
                    if (node[0] - 3, node[1] + 1) in neighbors:
                        neighbors += [(node[0] - 3, node[1] - 1)]
                    else:
                        neighbors += [(node[0] + 3, node[1] - 1)]
                else:  # node must be in blue
                    neighbors += [(node[0], node[1] - 2)]
            else:
                if node in self.layout.ancilla_qubits_green:
                    if (node[0] + 3, node[1] + 1) in neighbors:
                        neighbors += [(node[0], node[1] + 2)]

                    else:
                        neighbors += [(node[0] - 3, node[1] - 1)]
                else:
                    neighbors += [(node[0] - 3, node[1] + 1)]
                # pass  # TODO: green boundary
        return neighbors

    def find_sigma_and_rho_vertices(self, graph):
        """
        see equation 4 in "triangular color codes on trivalent graphs with
        flag qubits

        graph are the highlighted vertices!
        """
        boundary_nodes = [node for node in graph.nodes() if type(node) == str]

        red_boundaries = [node for node in boundary_nodes if node[1] == "R"]

        sigma_graphs = []
        sigma_vertices_of_one_component = []
        rho_vertices = set()
        colours = {"R", "B", "G"}
        try:
            red_boundary_neighbors = [n for n in graph.neighbors("vR")]
        except:
            red_boundary_neighbors = []

        while red_boundary_neighbors != []:

            start_node = red_boundary_neighbors[0]
            graph.remove_edge(start_node, "vR")
            path = [(("vR", start_node))]
            neighbors = self.find_neighbors(start_node, graph.edges())
            all_paths = []
            if neighbors != []:
                all_paths = self.find_path(
                    start_node, neighbors[0], list(graph.edges()), path, all_paths
                )
            else:
                all_paths = [path]
            path = max(all_paths, key=len)

            path_graph = nx.Graph()
            path_graph.add_edges_from(path)

            sigma_graphs.append(path_graph)

            red_boundary_neighbors = red_boundary_neighbors[1:]
            colours_boundary = {"R"}

            colours_boundary.add(path[-1][-1][1])
            colour_difference = colours.difference(colours_boundary)
            sigma_vertices = self.add_to_sigma_vertices(colour_difference, path_graph)

            sigma_vertices_of_one_component.append(sigma_vertices)

            if path[-1][-1] == "vR":
                red_boundary_neighbors.remove(path[-1][0])
            graph.remove_edges_from(path[1:])

        rho_vertices = self.add_to_rho_vertices(graph, rho_vertices)

        return (sigma_vertices_of_one_component, sigma_graphs, rho_vertices)

    def add_to_sigma_vertices(self, colour_difference, component):
        sigma_vertices = []
        if "G" not in colour_difference:
            for node in component:
                if node in self.layout.ancilla_qubits_blue:
                    sigma_vertices.append(node)
        else:  # if 'r'&'b' or 'r'&'r' we choose green
            for node in component:
                if node in self.layout.ancilla_qubits_green:
                    sigma_vertices.append(node)
        return sigma_vertices

    def find_path(self, start_node, next_node, graph_edges, path, all_paths):

        new_path = path + [(start_node, next_node)]
        new_graph_edges = copy.copy(graph_edges)

        try:
            new_graph_edges.remove((start_node, next_node))
        except:
            new_graph_edges.remove((next_node, start_node))

        path_neighbors = self.find_neighbors(next_node, new_graph_edges)
        if type(next_node) == str:
            all_paths.append(new_path)
            return all_paths
        elif path_neighbors == []:
            all_paths.append([])
            return all_paths
        else:

            for next_next_node in path_neighbors:
                all_paths = self.find_path(
                    next_node, next_next_node, new_graph_edges, new_path, all_paths
                )

        return all_paths

    def build_sigma_vertex_graph(self, component, graph):
        boundary_nodes = [node for node in component if type(node) == str]
        component.remove(boundary_nodes[0])
        component.remove(boundary_nodes[1])
        split_b_node_0 = boundary_nodes[0][1:-2].split(",")
        start_node = (int(split_b_node_0[0]), int(split_b_node_0[1]))

        split_b_node_1 = boundary_nodes[1][1:-2].split(",")
        end_node = (int(split_b_node_1[0]), int(split_b_node_1[1]))
        path = []
        path_nodes = {start_node}
        component_graph_edges = list(graph.subgraph(component).edges())
        path_neighbors = self.find_neighbors(start_node, component_graph_edges)
        if path_neighbors != []:
            for next_node in path_neighbors:
                new_path = self.find_path_edge(
                    start_node,
                    next_node,
                    end_node,
                    component,
                    copy.copy(component_graph_edges),
                    path,
                    path_nodes,
                )
                if new_path:
                    longest_path = new_path
                    break

            path_graph = nx.Graph()
            path_graph.add_edges_from(longest_path)
            path_graph.add_edge(start_node, boundary_nodes[0])
            path_graph.add_edge(end_node, boundary_nodes[1])
        return path_graph

    def find_neighbors(self, start_node, graph_edges):
        path_neighbors = []
        for pair in graph_edges:
            if pair[0] == start_node:
                path_neighbors.append(pair[1])
            elif pair[1] == start_node:
                path_neighbors.append(pair[0])
        return path_neighbors

    def add_to_rho_vertices(self, component, rho_vertices):
        for node in component:
            if node in self.layout.ancilla_qubits_red:
                if len(list(component.neighbors(node))) > 1:
                    rho_vertices.add(node)
        return rho_vertices

    def get_matching(self, error_cor, matching_graph, coords_matching_graph):
        red_green_error_index = [
            self.layout.ancilla_coords_to_index[coords]
            for coords in error_cor["red"].union(error_cor["green"])
        ]
        red_green_matching = self.create_two_colour_matching(
            red_green_error_index, matching_graph.green_red
        )

        red_green_matching = self.translate_pymatching(
            red_green_matching,
            coords_matching_graph.green_red,
            [len(self.layout.ancilla_qubits), len(self.layout.ancilla_qubits) + 1],
        )

        blue_green_error_index = [
            self.layout.ancilla_coords_to_index[coords]
            for coords in error_cor["blue"].union(error_cor["green"])
        ]
        blue_green_matching = self.create_two_colour_matching(
            blue_green_error_index, matching_graph.blue_green
        )
        blue_green_matching = self.translate_pymatching(
            blue_green_matching,
            coords_matching_graph.blue_green,
            [len(self.layout.ancilla_qubits), len(self.layout.ancilla_qubits) + 2],
        )

        red_blue_error_index = [
            self.layout.ancilla_coords_to_index[coords]
            for coords in error_cor["red"].union(error_cor["blue"])
        ]

        red_blue_matching = self.create_two_colour_matching(
            red_blue_error_index, matching_graph.red_blue
        )

        red_blue_matching = self.translate_pymatching(
            red_blue_matching,
            coords_matching_graph.red_blue,
            [len(self.layout.ancilla_qubits) + 1, len(self.layout.ancilla_qubits) + 2],
        )

        return (red_green_matching, red_blue_matching, blue_green_matching)

    def lift(self, vertex, neighbors):
        correction = self.find_correct_faces(vertex, neighbors)
        return correction  # TODO fix

    def translate_face_to_data_qubit(self, face):
        correction = set()
        for i in range(0, len(face), 3):
            correction.add(
                (
                    round((face[i][0] + face[i + 1][0] + face[i + 2][0]) / 3),
                    round((face[i][1] + face[i + 1][1] + face[i + 2][1]) / 3),
                )
            )
        return correction

    def find_correct_faces(self, vertex, neighbors):

        faces = self.find_faces(vertex)
        correction = set()

        i = 1
        neighbors_not_in_rho = self.find_neighbors_not_in_rho(vertex, neighbors)

        while i < len(faces):
            for face_subset_tuples in itertools.combinations(faces, i):
                face_subset = set()
                for face in face_subset_tuples:
                    face_subset ^= set(face)
                correct_faces = tuple(j for i in face_subset_tuples for j in i)
                should_break = False
                for node in neighbors:
                    if node not in face_subset:
                        should_break = True
                        break
                    # check if node appears in face subset twice
                else:
                    for node in neighbors_not_in_rho:
                        if node in face_subset:
                            should_break = True
                            break

                if should_break == False:
                    correction = self.translate_face_to_data_qubit(correct_faces)
                    for ancilla_qubit in correction:
                        if ancilla_qubit not in self.layout.data_qubits:
                            break
                    else:
                        i += len(faces)  # just to exit while loop
                        break
            i += 1
        return correction

    def find_neighbors_not_in_rho(self, vertex, neighbors_in_rho):
        neighbors_not_in_rho = set(
            {
                (vertex[0], vertex[1] + 2),
                (vertex[0] + 3, vertex[1] + 1),
                (vertex[0] + 3, vertex[1] - 1),
                (vertex[0], vertex[1] - 2),
                (vertex[0] - 3, vertex[1] - 1),
                (vertex[0] - 3, vertex[1] + 1),
            }
        )

        neighbors_not_in_rho ^= set(neighbors_in_rho)
        neighbors_not_in_rho &= self.layout.ancilla_qubits
        return neighbors_not_in_rho

    def find_faces(self, vertex):
        face_1 = (
            vertex,
            (vertex[0] + 3, vertex[1] + 1),
            (vertex[0] + 3, vertex[1] - 1),
        )
        face_2 = (vertex, (vertex[0] + 3, vertex[1] - 1), (vertex[0], vertex[1] - 2))
        face_3 = (vertex, (vertex[0], vertex[1] - 2), (vertex[0] - 3, vertex[1] - 1))
        face_4 = (
            vertex,
            (vertex[0] - 3, vertex[1] - 1),
            (vertex[0] - 3, vertex[1] + 1),
        )
        face_5 = (vertex, (vertex[0] - 3, vertex[1] + 1), (vertex[0], vertex[1] + 2))
        face_6 = (vertex, (vertex[0], vertex[1] + 2), (vertex[0] + 3, vertex[1] + 1))
        faces = [face_1, face_2, face_3, face_4, face_5, face_6]
        return faces

    def create_two_colour_matching(self, error_index, matching_graph):
        syndrome_vector = np.zeros(
            (len(self.layout.ancilla_coords_to_index.keys()) - 1)
        )
        for index in error_index:
            syndrome_vector[index] = 1
        correction = matching_graph.decode_to_matched_dets_dict(syndrome_vector)
        return correction  # matching

    def translate_pymatching(self, matching, graph, boundary_nodes):
        matching_graph = nx.Graph()
        pairs = []
        keys = list(matching.keys())
        #        print(keys, "keys")
        #        print(matching)
        for u in keys:
            v = matching[u]

            if u == None and v == None:
                pass

            elif u == None:
                try:
                    path_0 = nx.shortest_path(
                        graph, v, boundary_nodes[0], weight="weight"
                    )
                except:
                    path_0 = []
                try:
                    path_1 = nx.shortest_path(
                        graph, v, boundary_nodes[1], weight="weight"
                    )
                except:
                    path_1 = []

                if len(path_0) > len(path_1):
                    path_edges = [
                        (path_1[i], path_1[i + 1]) for i in range(len(path_1) - 1)
                    ]
                elif len(path_1) > len(path_0):
                    path_edges = [
                        (path_0[i], path_0[i + 1]) for i in range(len(path_0) - 1)
                    ]

                matching_graph.add_edges_from(path_edges)

            elif v == None:
                try:
                    path_0 = nx.shortest_path(
                        graph, u, boundary_nodes[0], weight="weight"
                    )

                except:
                    path_0 = []
                try:
                    path_1 = nx.shortest_path(
                        graph, u, boundary_nodes[1], weight="weight"
                    )
                except:
                    path_1 = []
                if len(path_0) > len(path_1):
                    path_edges = [
                        (path_1[i], path_1[i + 1]) for i in range(len(path_1) - 1)
                    ]

                elif len(path_1) > len(path_0):
                    path_edges = [
                        (path_0[i], path_0[i + 1]) for i in range(len(path_0) - 1)
                    ]
                matching_graph.add_edges_from(path_edges)
            else:
                path = nx.shortest_path(graph, u, v, weight="weight")
                path_edges = [(path[i], path[i + 1]) for i in range(len(path) - 1)]

                for node in path[:-1]:
                    matching_graph.add_edges_from(path_edges)
            if v != None:
                keys.remove(v)
        matching_graph = nx.relabel_nodes(
            matching_graph, self.layout.ancilla_index_to_coords
        )
        return matching_graph

    def translate_cmatching(self, matching, cmatching, node2id, graph, boundary_nodes):
        matching_graph = nx.Graph()
        pairs = []
        id2node = {v: k for k, v in node2id.items()}

        for i in range(0, matching, 2):
            u, v = id2node[cmatching[i]], id2node[cmatching[i + 1]]
            if type(u) == str and type(v) == str:
                pass

            elif type(u) == str:
                try:
                    path_0 = nx.shortest_path(graph, v, boundary_nodes[0])
                except:
                    path_0 = []
                try:
                    path_1 = nx.shortest_path(graph, v, boundary_nodes[1])
                except:
                    path_1 = []

                if len(path_0) > len(path_1):
                    path_edges = [
                        (path_1[i], path_1[i + 1]) for i in range(len(path_1) - 1)
                    ]
                elif len(path_1) > len(path_0):
                    path_edges = [
                        (path_0[i], path_0[i + 1]) for i in range(len(path_0) - 1)
                    ]

                matching_graph.add_edges_from(path_edges)

            elif type(v) == str:

                try:
                    path_0 = nx.shortest_path(graph, u, boundary_nodes[0])

                except:
                    path_0 = []
                try:
                    path_1 = nx.shortest_path(graph, u, boundary_nodes[1])
                except:
                    path_1 = []

                if len(path_0) > len(path_1):
                    path_edges = [
                        (path_1[i], path_1[i + 1]) for i in range(len(path_1) - 1)
                    ]

                elif len(path_1) > len(path_0):
                    path_edges = [
                        (path_0[i], path_0[i + 1]) for i in range(len(path_0) - 1)
                    ]
                matching_graph.add_edges_from(path_edges)
            else:
                path = nx.shortest_path(graph, u, v)
                path_edges = [(path[i], path[i + 1]) for i in range(len(path) - 1)]

                for node in path[:-1]:
                    matching_graph.add_edges_from(path_edges)
        return matching_graph

    def create_graph(self, error_cor, graph, boundary_vertices):
        num_edges = len(error_cor) ** 2
        num_nodes = len(error_cor) * 2
        error_ind = []
        nodes = []

        for er in error_cor:
            error_ind.append(er)
            nodes.append(er)
            nodes.append(str(er) + ", b")

        edges, node2id = self.create_edges(
            error_ind, nodes, num_edges, graph, boundary_vertices
        )
        return (num_nodes, num_edges, edges, node2id)

    def create_edges(self, error_ind, nodes, edge_num, graph, boundary_vertices):
        """
        graph is red-blue, green-blue, red-green etc.
        """
        edges = self.ffi.new("Edge[%d]" % (edge_num))
        e_ind = 0
        node2id = {val: index for index, val in enumerate(nodes)}
        for v1, v2 in it.combinations(error_ind, 2):
            edges[e_ind].uid = node2id[v1]
            edges[e_ind].vid = node2id[v2]

            # this runs dijkstra, need to add weight
            weight = nx.shortest_path(graph, source=v1, target=v2, weight="weight")
            edges[e_ind].weight = weight
            e_ind += 1

            v1_boundary = str(v1) + ", b"
            v2_boundary = str(v2) + ", b"
            edges[e_ind].uid = node2id[v1_boundary]
            edges[e_ind].vid = node2id[v2_boundary]
            edges[e_ind].weight = 0
            e_ind += 1

        for vertex in error_ind:
            boundary_vertex = str(vertex) + ", b"
            edges[e_ind].uid = node2id[vertex]
            edges[e_ind].vid = node2id[boundary_vertex]
            weight_b0 = nx.shortest_path_length(
                graph, source=vertex, target=boundary_vertices[0]
            )
            weight_b1 = nx.shortest_path_length(
                graph, source=vertex, target=boundary_vertices[1]
            )

            edges[e_ind].weight = int(min(weight_b0, weight_b1))
            e_ind += 1

        return (edges, node2id)

    def run_blossom(self, node_num, edge_num, edges):
        cmatching = self.ffi.new("int[%d]" % (node_num))  # when edges added * 2

        retVal = self.blossom.Init()
        retVal = self.blossom.Process(node_num, edge_num, edges)
        nMatching = self.blossom.GetMatching(cmatching)
        retVal = self.blossom.Clean()
        return (nMatching, cmatching)


if __name__ == "__main__":
    test = MWPM()
    test.single_run()
