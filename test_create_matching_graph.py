from create_matching_graph import Matching_graph
from error_model import BiasedNoiseModel

# from sim3d import Sim3D
from layout import Hexagonal_layout
import itertools as it
import numpy as np

g5 = Matching_graph(5, n_layers=1)
# d5_setup = Hexagonal_layout(5)


def test_init():
    assert (g5.graphX.blue_green.has_edge((5, 5), "vB")) == True
    assert (g5.graphX.green_red.has_edge((5, 5), "vR")) == True
    # assert (g5.graphX.red_blue.has_edge((5, 5), 'vB')) == True
    assert (8, 2) in g5.green_nodes


d5layout = Hexagonal_layout(5)

noise_model = BiasedNoiseModel(0.1, 10, d5layout)
g5biased = Matching_graph(5, noise_model)


def test_translate_error_model():
    g5biased.translate_error_model()


# write a test to check that only edges between two boundary nodes are weight 0
def test_update_weights_of_edges():
    g5biased.update_weights_of_edges()


def test_biased_noise_model():
    biased_noise_model = BiasedNoiseModel(0.1, 10, d5layout)
    g5biased = Matching_graph(5, biased_noise_model)
    print(g5biased.red_blue.edges(), "red blue")
    print(g5biased.blue_green.edges(), "blue green")
    print(g5biased.blue_green.boundary)
    print(g5biased.green_red.edges(), "green red")

    unbiased_noise_model = BiasedNoiseModel(0.1, 1, d5layout)
    g5unbiased = Matching_graph(5, unbiased_noise_model)
    print(g5unbiased.red_blue.edges(), "red blue")
    

test_biased_noise_model()
