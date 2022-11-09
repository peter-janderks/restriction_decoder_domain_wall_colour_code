import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import itertools
from layout import Hexagonal_layout


class visualize():
    def __init__(self, n_layers, layout):
        """
        points: list with three entries: dict containing cords of red ancillas
                                    dict containing cord of blue anicallas
                                    dict containing cord of green ancillas
        """
        self.fig = plt.figure()
        self.distance = layout.distance
        self.n_layers = n_layers
        self.layout = layout

    def plot_matching_graph(self, matching_graph, n_graphs=6):
        self.matching_graph = matching_graph
        # = matching_graph
        if n_graphs == 1:
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            self.single_plot_matching_graph(
                ax, ('red', 'blue'), self.matching_graph.graphX.red_blue.edges)

        else:
            self.single_plot_matching_graph(self.fig.add_subplot(2, 3, 1, projection='3d'),
                                            ('red', 'green'), self.matching_graph.graphX.green_red.edges)
            self.single_plot_matching_graph(self.fig.add_subplot(2, 3, 2, projection='3d'),
                                            ('red', 'blue'), self.matching_graph.graphX.red_blue.edges)

            self.single_plot_matching_graph(self.fig.add_subplot(2, 3, 3, projection='3d'),
                                            ('green', 'blue'), self.matching_graph.graphX.blue_green.edges)

            self.single_plot_matching_graph(self.fig.add_subplot(2, 3, 4, projection='3d'),
                                            ('red', 'green'), self.matching_graph.graphZ.green_red.edges)

            self.single_plot_matching_graph(self.fig.add_subplot(2, 3, 5, projection='3d'),
                                            ('red', 'blue'), self.matching_graph.graphZ.red_blue.edges)

            self.single_plot_matching_graph(self.fig.add_subplot(2, 3, 6, projection='3d'),
                                            ('green', 'blue'), self.matching_graph.graphZ.blue_green.edges)

        plt.show()

    def single_plot_matching_graph(self, ax, colours, edges):
        self.plot_body_nodes(ax, colours)
        self.plot_boundary_nodes(ax, colours)
        self.plot_edges(ax, edges, colours)

    def plot_2d_layer(self, colours, plot_edges=False, matching_graph=None):
        ax = plt.axes()
        self.n_layers = 1
        self.plot_body_nodes(ax, colours)
        self.plot_boundary_nodes(ax, colours)
        if plot_edges == True:
            self.matching_graph = matching_graph
            if 'green' in colours and 'blue' in colours:
                edges = self.matching_graph.graphZ.blue_green.edges
                self.plot_edges_2d(ax, edges, ('green', 'blue'))

            if 'red' in colours and 'blue' in colours:
                edges = self.matching_graph.graphZ.red_blue.edges
                self.plot_edges_2d(ax, edges, ('red', 'blue'))

            if 'green' in colours and 'red' in colours:

                edges = self.matching_graph.graphZ.green_red.edges
                self.plot_edges_2d(ax, edges, ('green', 'red'))

        plt.show()

    def plot_syndrome(self, syndrome_x, syndrome_z, error_history=None):
        z_syndrome_ax = self.fig.add_subplot(1, 2, 1, projection='3d')
        x_syndrome_ax = self.fig.add_subplot(1, 2, 2, projection='3d')

        self.single_plot_syndrome(
            z_syndrome_ax, syndrome_x, error_history, 'Z')
        self.single_plot_syndrome(
            x_syndrome_ax, syndrome_z, error_history, 'X')

        if error_history:
            for layer, error_layer in enumerate(error_history):
                if error_layer.x_set:
                    for error in error_layer.x_set:
                        z_syndrome_ax.scatter(error[0], error[1],
                                              layer, color='purple', marker='s')

                if error_layer.z_set:
                    for error in error_layer.x_set:
                        x_syndrome_ax.scatter(error[0], error[1],
                                              layer, color='purple', marker='s')

        plt.show()

    def single_plot_syndrome(self, ax, syndrome, error_history, syndrome_type='X'):
        """ should add highlighting of data errors that fire"""
        self.plot_body_nodes(ax, ('red', 'green', 'blue'))
        self.plot_boundary_nodes(ax, ('red', 'green', 'blue'))
        self.plot_firing_ancillas(ax, syndrome)

    def plot_firing_ancillas(self, ax, syndrome_hist):
        """does not plot firing flags!"""
        """ seems to be some bug in syndrome hist"""
        for layer, syndrome_one_round in enumerate(syndrome_hist):
            for syndrome in syndrome_one_round:
                if type(syndrome) != str:
                    ax.scatter(syndrome[0], syndrome[1],
                               layer, color='grey', marker='s')

    def plot_body_nodes(self, ax, colours):
        if 'blue' in colours:

            for node in self.layout.ancilla_qubits_blue:
                print(node)
                self.plot_point_multiple_layers(ax, node[0], node[1], 'blue')

        if 'red' in colours:
            for node in self.layout.ancilla_qubits_red:
                self.plot_point_multiple_layers(ax, node[0], node[1], 'red')

        if 'green' in colours:
            for node in self.layout.ancilla_qubits_green:
                self.plot_point_multiple_layers(ax, node[0], node[1], 'green')

        for node in self.layout.data_qubits:
            self.plot_point_multiple_layers(
                ax, node[0], node[1], 'black')

    def plot_point_multiple_layers(self, ax, x_cor, y_cor, p_color):

        if self.n_layers > 1:
            i = 0
            while i < self.n_layers:
                ax.scatter(x_cor, y_cor, i, color=p_color)
                i += 1
        else:
            ax.scatter(x_cor, y_cor, color=p_color)

    def plot_boundary_nodes(self, ax, colours):
        for layer in range(self.distance):
            if 'green' in colours:
                ax.scatter((self.distance * 3 - 3)/2, -
                           (self.distance/2), layer, color='green')

            if 'red' in colours:
                ax.scatter((self.distance * 3 - 3) + (self.distance/2),
                           ((self.distance+1)//3 * 3)/2, layer, color='red')

            if 'blue' in colours:
                ax.scatter(-self.distance/2,
                           ((self.distance+1)//3 * 3)/2, layer, color='blue')

    def plot_edges_2d(self, ax, edges, colours):
        for edge in edges:
            if type(edge[0]) == tuple and type(edge[1]) == tuple:
                ax.plot([edge[0][0], edge[1][0]],
                        [edge[0][1], edge[1][1]],
                        color='black')
            elif type(edge[0]) == tuple or type(edge[1]) == tuple:
                self.plot_boundary_edge(
                    ax, edge, colours, two_dimensional=True)

            else:
                print(edge, 'edge')
                print(colours, 'colours')
                self.plot_corner_edge(ax, edge, colours)

    def plot_edges(self, ax, edges, colours):
        for edge in edges:
            if type(edge[0]) == tuple and type(edge[1]) == tuple:
                ax.plot([edge[0][0], edge[1][0]],
                        [edge[0][1], edge[1][1]],
                        [edge[0][2], edge[1][2]], color='black')
            elif type(edge[0]) == tuple or type(edge[1]) == tuple:
                self.plot_boundary_edge(ax, edge, colours)

            else:
                self.plot_corner_edge(ax, edge, colours)

    def plot_boundary_edge(self, ax, edge, colours, two_dimensional=False):

        if 'vR' in edge[1]:
            x_1 = (self.distance * 3 - 3) + (self.distance/2)
            y_1 = ((self.distance+1)//3 * 3)/2

        elif 'vB' in edge[1]:
            x_1 = (-self.distance/2)
            y_1 = ((self.distance+1)//3 * 3)/2

        else:
            x_1 = (self.distance * 3 - 3)/2
            y_1 = -(self.distance/2)
        if two_dimensional == False:

            ax.plot([edge[0][0], x_1],
                    [edge[0][1], y_1], [int(edge[1][-1]), int(edge[1][-1])], color='black')
        else:
            ax.plot([edge[0][0], x_1],
                    [edge[0][1], y_1], color='black')

    def plot_corner_edge(self, ax, edge, colours):

        if 'blue' in colours and 'green' in colours:
            ax.plot([(-self.distance/2), -2],
                    [((self.distance+1)//3 * 3)/2, -2], color='black')
            ax.plot([(self.distance * 3 - 3)/2, -2],
                    [-(self.distance/2), -2], color='black')

        elif 'blue' in colours and 'red' in colours:
            ax.plot([(-self.distance/2), self.distance+2],
                    [((self.distance+1)//3 * 3)/2, self.distance+2], color='black')
            ax.plot([(self.distance * 3 - 3) + (self.distance/2), self.distance+2],
                    [((self.distance+1)//3 * 3)/2, self.distance+2], color='black')

        else:
            ax.plot([(self.distance * 3 - 3)/2, self.distance*3-1],
                    [-(self.distance/2), -2], color='black')
            ax.plot([(self.distance * 3 - 3) + (self.distance/2), self.distance*3-1],
                    [((self.distance+1)//3 * 3)/2, -2], color='black')
