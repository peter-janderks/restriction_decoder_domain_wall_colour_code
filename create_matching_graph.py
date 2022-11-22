from layout import Hexagonal_layout
import networkx as nx
import matplotlib.pyplot as plt
import itertools as it
import numpy as np
from pymatching import Matching


class X_or_Z_graph(object):
    def __init__(self):
        self.red_blue = nx.Graph()
        self.green_red = nx.Graph()
        self.blue_green = nx.Graph()


class Matching_graph:
    def __init__(self, distance, error_model=None, n_layers=None):
        self.n_layers = 1
        self.distance = distance
        self.layout = Hexagonal_layout(distance)

        self.initialize_nodes()
        self.initialize_edges()
        if error_model != None:
            self.error_model = error_model
            self.translate_error_model()
            self.update_weights_of_edges(self.graphX)
            self.update_weights_of_edges(self.graphZ)
            self.init_pymatching_graph()

    def initialize_edges(self):
        self.graphX = X_or_Z_graph()

        for i in range(self.n_layers):
            self.create_body_edges_one_layer()
            self.create_boundary_edges_one_layer()

        self.graphZ = X_or_Z_graph()
        self.graphZ.red_blue = self.graphX.red_blue.copy()
        self.graphZ.green_red = self.graphX.green_red.copy()
        self.graphZ.blue_green = self.graphX.blue_green.copy()

    def initialize_nodes(self):
        self.blue_nodes = set(
            self.layout.ancilla_coords_to_index[coords]
            for coords in self.layout.ancilla_qubits_blue
        )
        self.blue_nodes.add("vB")
        self.green_nodes = set(
            self.layout.ancilla_coords_to_index[coords]
            for coords in self.layout.ancilla_qubits_green
        )
        self.green_nodes.add("vG")

        self.red_nodes = set(
            self.layout.ancilla_coords_to_index[coords]
            for coords in self.layout.ancilla_qubits_red
        )

        self.red_nodes.add("vR")

        self.green_nodes_coords = self.layout.ancilla_qubits_green.copy()
        #        self.green_nodes.add("vG")

        self.red_nodes_coords = self.layout.ancilla_qubits_red.copy()
        #        self.red_nodes.add("vR")
        self.blue_nodes_coords = self.layout.ancilla_qubits_blue.copy()

    def create_body_edges_one_layer(self):
        # make it depend on error model!

        for node in self.blue_nodes_coords:
            if type(node) == tuple:
                potential_neighbors = self.ancilla_neighbors(node)
                for pot_neigh in potential_neighbors:
                    if pot_neigh in self.red_nodes_coords:
                        self.graphX.red_blue.add_edge(
                            self.layout.ancilla_coords_to_index[(node[0], node[1])],
                            self.layout.ancilla_coords_to_index[
                                (pot_neigh[0], pot_neigh[1])
                            ],
                            error_probability=0,
                        )

                    elif pot_neigh in self.green_nodes_coords:
                        self.graphX.blue_green.add_edge(
                            self.layout.ancilla_coords_to_index[(node[0], node[1])],
                            self.layout.ancilla_coords_to_index[
                                (pot_neigh[0], pot_neigh[1])
                            ],
                            error_probability=0,
                        )

        for node in self.red_nodes_coords:
            if type(node) == tuple:
                potential_neighbors = self.ancilla_neighbors(node)
                for pot_neigh in potential_neighbors:
                    if pot_neigh in self.green_nodes_coords:
                        self.graphX.green_red.add_edge(
                            self.layout.ancilla_coords_to_index[(node[0], node[1])],
                            self.layout.ancilla_coords_to_index[
                                (pot_neigh[0], pot_neigh[1])
                            ],
                            error_probability=0,
                        )

    def create_boundary_edges_one_layer(self):
        for node in self.layout.ancilla_blue_boundary:
            node = self.layout.ancilla_coords_to_index[node]
            if node in self.red_nodes:
                self.graphX.red_blue.add_edge(
                    node,
                    self.layout.ancilla_coords_to_index["vB"],
                    error_probability=0,
                )
            elif node in self.green_nodes:
                self.graphX.blue_green.add_edge(
                    node,
                    self.layout.ancilla_coords_to_index["vB"],
                    error_probability=0,
                )

        for node in self.layout.ancilla_green_boundary:
            node = self.layout.ancilla_coords_to_index[node]
            if node in self.red_nodes:
                self.graphX.green_red.add_edge(
                    node,
                    self.layout.ancilla_coords_to_index["vG"],
                    error_probability=0,
                )

            elif node in self.blue_nodes:
                self.graphX.blue_green.add_edge(
                    node,
                    self.layout.ancilla_coords_to_index["vG"],
                    error_probability=0,
                )

        for node in self.layout.ancilla_red_boundary:
            node = self.layout.ancilla_coords_to_index[node]
            if node in self.blue_nodes:
                self.graphX.red_blue.add_edge(
                    node,
                    self.layout.ancilla_coords_to_index["vR"],
                    error_probability=0,
                )

            elif node in self.green_nodes:
                self.graphX.green_red.add_edge(
                    node,
                    self.layout.ancilla_coords_to_index["vR"],
                    error_probability=0,
                )

        self.graphX.blue_green.add_edge(
            self.layout.ancilla_coords_to_index["vB"],
            self.layout.ancilla_coords_to_index["vG"],
            error_probability=0,
            weight=0.0,
        )
        self.graphX.red_blue.add_edge(
            self.layout.ancilla_coords_to_index["vB"],
            self.layout.ancilla_coords_to_index["vR"],
            error_probability=0,
            weight=0.0,
        )
        self.graphX.green_red.add_edge(
            self.layout.ancilla_coords_to_index["vG"],
            self.layout.ancilla_coords_to_index["vR"],
            error_probability=0,
            weight=0.0,
        )

    def update_weights_of_edges(self, graph):
        graph.flipped_red_blue_edges = set()
        graph.flipped_blue_green_edges = set()
        graph.flipped_green_red_edges = set()
        for node_1, node_2 in graph.red_blue.edges():
            old_weight = graph.red_blue[node_1][node_2]["error_probability"]
            if old_weight != 0:
                graph.red_blue[node_1][node_2]["weight"] = np.log(
                    (1 - old_weight) / old_weight
                )
        for node_1, node_2 in graph.blue_green.edges():
            old_weight = graph.blue_green[node_1][node_2]["error_probability"]
            if old_weight != 0:
                graph.blue_green[node_1][node_2]["weight"] = np.log(
                    (1 - old_weight) / old_weight
                )

        for node_1, node_2 in graph.green_red.edges():
            old_weight = graph.green_red[node_1][node_2]["error_probability"]
            if old_weight != 0:
                graph.green_red[node_1][node_2]["weight"] = np.log(
                    (1 - old_weight) / old_weight
                )

    def init_pymatching_graph(self):
        self.red_blue_X = Matching(self.graphX.red_blue)
        self.blue_green_X = Matching(self.graphX.blue_green)
        self.green_red_X = Matching(self.graphX.green_red)
        self.red_blue_X.set_boundary_nodes(
            {
                len(self.layout.ancilla_qubits) + 2,
                len(self.layout.ancilla_qubits) + 1,
            }
        )
        self.blue_green_X.set_boundary_nodes(
            {
                len(self.layout.ancilla_qubits) + 2,
                len(self.layout.ancilla_qubits),
            }
        )
        self.green_red_X.set_boundary_nodes(
            {
                len(self.layout.ancilla_qubits),
                len(self.layout.ancilla_qubits) + 1,
            }
        )

        self.red_blue_Z = Matching(self.graphZ.red_blue)
        self.blue_green_Z = Matching(self.graphZ.blue_green)
        self.green_red_Z = Matching(self.graphZ.green_red)
        self.red_blue_Z.set_boundary_nodes(
            {
                len(self.layout.ancilla_qubits) + 2,
                len(self.layout.ancilla_qubits) + 1,
            }
        )
        self.blue_green_Z.set_boundary_nodes(
            {
                len(self.layout.ancilla_qubits) + 2,
                len(self.layout.ancilla_qubits),
            }
        )
        self.green_red_Z.set_boundary_nodes(
            {
                len(self.layout.ancilla_qubits),
                len(self.layout.ancilla_qubits) + 1,
            }
        )

    def translate_error_model(self):
        # translate single error needs to be given a graph...

        for (
            data_qubit,
            error_probability,
        ) in self.error_model.low_error_probability_dict_X.items():
            self.translate_single_error(data_qubit, error_probability, self.graphX)

        for (
            data_qubit,
            error_probability,
        ) in self.error_model.low_error_probability_dict_Z.items():
            self.translate_single_error(data_qubit, error_probability, self.graphZ)

    def translate_single_error(self, data_qubit, probability, graph):
        ancilla_qubits = [
            self.layout.ancilla_coords_to_index[ancilla_qubit]
            for ancilla_qubit in self.layout.data_qubits_to_ancilla_qubits([data_qubit])
        ]
        if len(ancilla_qubits) > 2:
            for ancilla_pair in it.combinations(ancilla_qubits, 2):
                if (
                    ancilla_pair[0] in graph.red_blue.nodes()
                    and ancilla_pair[1] in graph.red_blue.nodes()
                ):
                    graph.red_blue[ancilla_pair[0]][ancilla_pair[1]][
                        "error_probability"
                    ] += probability

                elif (
                    ancilla_pair[0] in graph.green_red.nodes()
                    and ancilla_pair[1] in graph.green_red.nodes()
                ):
                    graph.green_red[ancilla_pair[0]][ancilla_pair[1]][
                        "error_probability"
                    ] += probability

                elif (
                    ancilla_pair[0] in graph.blue_green.nodes()
                    and ancilla_pair[1] in graph.blue_green.nodes()
                ):
                    graph.blue_green[ancilla_pair[0]][ancilla_pair[1]][
                        "error_probability"
                    ] += probability

        # for two need to add to boundaries:
        elif len(ancilla_qubits) == 2:
            if (
                ancilla_qubits[0] in graph.red_blue.nodes()
                and ancilla_qubits[1] in graph.red_blue.nodes()
            ):
                graph.red_blue[ancilla_qubits[0]][ancilla_qubits[1]][
                    "error_probability"
                ] += probability
                # is on green boundary
                for qubit in ancilla_qubits:
                    if qubit in graph.green_red:
                        graph.green_red[qubit][
                            self.layout.ancilla_coords_to_index["vG"]
                        ]["error_probability"] += probability
                    else:
                        graph.blue_green[qubit][
                            self.layout.ancilla_coords_to_index["vG"]
                        ]["error_probability"] += probability
            elif (
                ancilla_qubits[0] in graph.green_red.nodes()
                and ancilla_qubits[1] in graph.green_red.nodes()
            ):
                graph.green_red[ancilla_qubits[0]][ancilla_qubits[1]][
                    "error_probability"
                ] += probability
                # is on blue boundary
                for qubit in ancilla_qubits:
                    if qubit in graph.blue_green:
                        graph.blue_green[qubit][
                            self.layout.ancilla_coords_to_index["vB"]
                        ]["error_probability"] += probability
                    else:
                        graph.red_blue[qubit][
                            self.layout.ancilla_coords_to_index["vB"]
                        ]["error_probability"] += probability
            elif (
                ancilla_qubits[0] in graph.blue_green.nodes()
                and ancilla_qubits[1] in graph.blue_green.nodes()
            ):
                graph.blue_green[ancilla_qubits[0]][ancilla_qubits[1]][
                    "error_probability"
                ] += probability
                # is on red boundary
                for qubit in ancilla_qubits:
                    if qubit in graph.green_red:
                        graph.green_red[qubit][
                            self.layout.ancilla_coords_to_index["vR"]
                        ]["error_probability"] += probability
                    else:
                        graph.red_blue[qubit][
                            self.layout.ancilla_coords_to_index["vR"]
                        ]["error_probability"] += probability

        else:
            # this handles corner ancillas
            # oh this checking is weird
            if ancilla_qubits[0] in self.red_nodes:
                graph.green_red[ancilla_qubits[0]][
                    self.layout.ancilla_coords_to_index["vG"]
                ]["error_probability"] += probability

                graph.red_blue[ancilla_qubits[0]][
                    self.layout.ancilla_coords_to_index["vB"]
                ]["error_probability"] += probability

            elif ancilla_qubits[0] in self.green_nodes:
                graph.green_red[ancilla_qubits[0]][
                    self.layout.ancilla_coords_to_index["vR"]
                ]["error_probability"] += probability

                graph.blue_green[ancilla_qubits[0]][
                    self.layout.ancilla_coords_to_index["vB"]
                ]["error_probability"] += probability

            else:
                graph.blue_green[ancilla_qubits[0]][
                    self.layout.ancilla_coords_to_index["vG"]
                ]["error_probability"] += probability

                graph.red_blue[ancilla_qubits[0]][
                    self.layout.ancilla_coords_to_index["vR"]
                ]["error_probability"] += probability

    def ancilla_neighbors(self, n):
        neighbors = (
            (n[0], n[1] + 2),
            (n[0] + 3, n[1] + 1),
            (n[0] + 3, n[1] - 1),
            (n[0], n[1] - 2),
            (n[0] - 3, n[1] - 1),
            (n[0] - 3, n[1] + 1),
        )
        return neighbors


#        for data_qubit in self.error_model = ColourCodeXErrorModel(per, self.layout.data_qubits)
