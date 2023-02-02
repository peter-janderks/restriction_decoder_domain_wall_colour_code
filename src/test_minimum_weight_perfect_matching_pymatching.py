import pytest
from create_matching_graph import Matching_graph
from minimum_weight_perfect_matching_pymatching import MWPM
import networkx as nx
from error_model import BiasedNoiseModel
from layout import Hexagonal_layout


# layout = Hexagonal_layout(7)
layout5 = Hexagonal_layout(5)
# decoder = MWPM(
#    7,
#    layout,
#    BiasedNoiseModel(0.4, 100, layout),
# )
g5biased = Matching_graph(5, BiasedNoiseModel(0.3, 1000, layout5))

decoder5 = MWPM(
    5,
    layout5,
    BiasedNoiseModel(0.3, 1000, layout5),
)

# print(decoder.matching_graph_X.green_red.edges())
# unbiased_decoder = MWPM(
#    7,
#    layout,
#    BiasedNoiseModel(0.1, 0.5, layout),
# )
# decoder_5 = MWPM(5)
# decoder_3 = MWPM(3)


def data_qubit_errors_to_ancilla(data_qubit_error_X, data_qubit_error_Z, decoder):
    ancilla_error_cords_X = decoder.layout.data_qubits_to_ancilla_qubits(
        data_qubit_error_X
    )
    ancilla_error_X = {
        "blue": ancilla_error_cords_X.intersection(decoder.layout.ancilla_qubits_blue),
        "red": ancilla_error_cords_X.intersection(decoder.layout.ancilla_qubits_red),
        "green": ancilla_error_cords_X.intersection(
            decoder.layout.ancilla_qubits_green
        ),
    }

    ancilla_error_cords_Z = decoder.layout.data_qubits_to_ancilla_qubits(
        data_qubit_error_Z
    )
    ancilla_error_Z = {
        "blue": ancilla_error_cords_Z.intersection(decoder.layout.ancilla_qubits_blue),
        "red": ancilla_error_cords_Z.intersection(decoder.layout.ancilla_qubits_red),
        "green": ancilla_error_cords_Z.intersection(
            decoder.layout.ancilla_qubits_green
        ),
    }

    return (ancilla_error_X, ancilla_error_Z)


def test_breaking_error():
    data_qubit_error_X = {(6, 6), (10, 4)}
    data_qubit_error_Z = {(6, 6)}
    (
        ancilla_error_X,
        ancilla_error_Z,
    ) = data_qubit_errors_to_ancilla(data_qubit_error_X, data_qubit_error_Z)
    correction_Z = decoder.single_run(
        decoder.matching_graph_Z, decoder.coords_matching_graph_Z, ancilla_error_Z
    )
    assert correction_Z == set()

    correction_X = decoder.single_run(
        decoder.matching_graph_X, decoder.coords_matching_graph_X, ancilla_error_X
    )
    resulting_operator = correction_X ^ data_qubit_error_X
    assert len(resulting_operator) % 2 == 0

def test_error_infinite_bias():

    layout3 = Hexagonal_layout(3)
    g3_inft_biased = Matching_graph(3, BiasedNoiseModel(0.4,'infty', layout3))

    decoder3 = MWPM(
        3,
        layout3,
        BiasedNoiseModel(0.4,'infty', layout3)
    )
    data_qubit_error_X = set()
    data_qubit_error_Z = {(3,3),(4,2)}
    (
        ancilla_error_X,
        ancilla_error_Z,
    ) = data_qubit_errors_to_ancilla(data_qubit_error_X, data_qubit_error_Z,decoder3)
    print(decoder3.coords_matching_graph_Z.blue_green.edges.data(),'Z')
    correction_Z = decoder3.single_run(
        decoder3.matching_graph_Z, decoder3.coords_matching_graph_Z, ancilla_error_Z
    )
#    print(
    resulting_operator_Z = data_qubit_error_Z ^ correction_Z

    correction_X = decoder3.single_run(
        decoder3.matching_graph_X, decoder3.coords_matching_graph_X, ancilla_error_X
    )
    resulting_operator_X = correction_X ^ data_qubit_error_X
    assert len(resulting_operator_X) % 2 == 0
    assert len(resulting_operator_Z) % 2 == 0

def test_error_bias():
    layout5 = Hexagonal_layout(5)
    error_model = BiasedNoiseModel(0.2,30, layout5)
    g5_inft_biased = Matching_graph(5, error_model)

    decoder5 = MWPM(
        5,
        layout5,
        error_model
    )
    print(error_model.error_probability_dict_X)
    print(error_model.error_probability_dict_Z,'Z')
    data_qubit_error_X =  {(3,1),(4,2)}
    data_qubit_error_Z = set()
    (
        ancilla_error_X,
        ancilla_error_Z,
    ) = data_qubit_errors_to_ancilla(data_qubit_error_X, data_qubit_error_Z,decoder5)
    print(decoder5.coords_matching_graph_Z.blue_green.edges.data(),'Z')
    correction_Z = decoder5.single_run(
        decoder5.matching_graph_Z, decoder5.coords_matching_graph_Z, ancilla_error_Z
    )
    print(correction_Z)
#    print(
    resulting_operator_Z = data_qubit_error_Z ^ correction_Z

    correction_X = decoder5.single_run(
        decoder5.matching_graph_X, decoder5.coords_matching_graph_X, ancilla_error_X
    )
    print(correction_X)
    resulting_operator_X = correction_X ^ data_qubit_error_X
    assert len(resulting_operator_X) % 2 == 0
    assert len(resulting_operator_Z) % 2 == 0

test_error_bias()
def test_weight_two_error():

    data_qubit_error_X = {(15, 1), (16, 0)}
    data_qubit_error_Z = {(12, 0), (10, 0)}
    (
        ancilla_error_X,
        ancilla_error_Z,
    ) = data_qubit_errors_to_ancilla(data_qubit_error_X, data_qubit_error_Z, decoder)
    correction_Z = decoder.single_run(
        decoder.matching_graph_Z, decoder.coords_matching_graph_Z, ancilla_error_Z
    )   
    resulting_operator_Z = correction_Z ^ data_qubit_error_Z
    assert resulting_operator_Z == set()

    correction_X = decoder.single_run(
        decoder.matching_graph_X, decoder.coords_matching_graph_X, ancilla_error_X
    )
    resulting_operator_X = correction_X ^ data_qubit_error_X
    assert resulting_operator_X == set()


def test_weird_two_qubit_error():
    e_model = BiasedNoiseModel(0.3, 1000, layout5)
    data_qubit_error_X = set()
    data_qubit_error_Z = {(12, 0), (7, 1)}
    (
        ancilla_error_X,
        ancilla_error_Z,
    ) = data_qubit_errors_to_ancilla(data_qubit_error_X, data_qubit_error_Z, decoder5)
    correction_Z = decoder5.single_run(
        decoder5.matching_graph_Z, decoder5.coords_matching_graph_Z, ancilla_error_Z
    )

    resulting_operator_Z = correction_Z ^ data_qubit_error_Z
    print(resulting_operator_Z)
    assert resulting_operator_Z == set()

    correction_X = decoder5.single_run(
        decoder5.matching_graph_X, decoder5.coords_matching_graph_X, ancilla_error_X
    )
    resulting_operator_X = correction_X ^ data_qubit_error_X
    assert resulting_operator_X == set()


def test_weight_three_error():
    data_qubit_error = {(16, 2), (10, 0), (7, 1)}
    ancilla_error_cords = decoder.layout.data_qubits_to_ancilla_qubits(data_qubit_error)

    ancilla_error = {
        "blue": ancilla_error_cords.intersection(decoder.layout.ancilla_qubits_blue),
        "red": ancilla_error_cords.intersection(decoder.layout.ancilla_qubits_red),
        "green": ancilla_error_cords.intersection(decoder.layout.ancilla_qubits_green),
    }
    print(ancilla_error, "a error")
    correction = decoder.single_run(ancilla_error)
    print(correction, "correction")
    resulting_operator = correction ^ data_qubit_error
    print(resulting_operator)
    #    assert resulting_operator == {(9, 1), (10, 0), (7, 1), (6, 0)}

    data_qubit_error = {(6, 2), (7, 7), (9, 3)}
    ancilla_error_cords = decoder.layout.data_qubits_to_ancilla_qubits(data_qubit_error)
    ancilla_error = {
        "blue": ancilla_error_cords.intersection(decoder.layout.ancilla_qubits_blue),
        "red": ancilla_error_cords.intersection(decoder.layout.ancilla_qubits_red),
        "green": ancilla_error_cords.intersection(decoder.layout.ancilla_qubits_green),
    }
    correction = decoder.single_run(ancilla_error)
    resulting_operator = correction ^ data_qubit_error
    assert resulting_operator == set()

    data_qubit_error = {(6, 6), (6, 4), (4, 2)}
    ancilla_error_cords = decoder.layout.data_qubits_to_ancilla_qubits(data_qubit_error)

    ancilla_error = {
        "blue": ancilla_error_cords.intersection(decoder.layout.ancilla_qubits_blue),
        "red": ancilla_error_cords.intersection(decoder.layout.ancilla_qubits_red),
        "green": ancilla_error_cords.intersection(decoder.layout.ancilla_qubits_green),
    }
    correction = decoder.single_run(ancilla_error)

    resulting_operator = correction ^ data_qubit_error
    print(resulting_operator, "resulting")
    assert resulting_operator == set()


# test_weight_three_error()


def test_error_breaking():
    errors = {"blue": {(17, 1), (11, 1)}, "red": set(), "green": set()}
    correction = decoder.single_run(error_cor=errors)
    assert correction == {(12, 0), (16, 0)}


def test_correction_without_boundaries():
    errors_blue = set({(11, 7)})
    errors_red = set({(8, 6)})
    errors_green = set({(11, 5)})

    errors = {"blue": errors_blue, "red": errors_red, "green": errors_green}
    correction = decoder.single_run(error_cor=errors)
    assert correction == {(10, 6)}

    errors = {"blue": {}, "red": {(5, 3)}, "green": {(5, 5)}}
    correction = decoder.single_run(error_cor=errors)
    assert correction == {(4, 4)}

    #    errors = {'blue': {(11, 1)}, 'red': {(8, 6)}, 'green': {(8, 8), (8, 2)}}
    #    correction = decoder.single_run(error_cor=errors)
    #    assert correction == {(15, 3), (13, 3), (7, 7),
    #                     (10, 2)}  # TODO: check this

    errors = {"blue": {}, "green": {(5, 5)}, "red": {(8, 6)}}
    correction = decoder.single_run(error_cor=errors)
    assert correction == {(6, 6)}

    errors = {"blue": {(11, 1)}, "red": {(8, 0)}, "green": set()}
    correction = decoder.single_run(error_cor=errors)
    assert correction == {(10, 0)}

    errors = {"blue": {(14, 4)}, "red": set(), "green": {(14, 2)}}
    correction = decoder.single_run(error_cor=errors)
    assert correction == {(15, 3)}

    errors = {"blue": {(14, 4)}, "red": set(), "green": {(11, 5)}}
    correction = decoder.single_run(error_cor=errors)
    assert correction == {(13, 5)}

    errors = {"blue": {(11, 7)}, "red": set(), "green": {(11, 5)}}
    correction = decoder.single_run(error_cor=errors)
    assert correction == {(12, 6)}

    errors = {"blue": {(5, 1)}, "red": {(8, 0)}, "green": set()}
    correction = decoder.single_run(error_cor=errors)
    assert correction == {(6, 0)}

    errors = {"blue": {}, "red": {(2, 0)}, "green": {(2, 2)}}
    correction = decoder.single_run(error_cor=errors)
    assert correction == {(1, 1)}

    errors = {"blue": {(11, 1)}, "red": {(8, 0)}, "green": set()}
    correction = decoder.single_run(error_cor=errors)
    assert correction == {(10, 0)}

    errors = {"blue": {(11, 1)}, "red": {(14, 0)}, "green": set()}
    correction = decoder.single_run(error_cor=errors)
    assert correction == {(12, 0)}

    errors = {"blue": {(5, 1)}, "red": {(8, 0)}, "green": {(8, 2)}}
    correction = decoder.single_run(error_cor=errors)
    assert correction == {(7, 1)}

    data_qubit_error = {(12, 0), (7, 1)}
    ancilla_error_cords = decoder.layout.data_qubits_to_ancilla_qubits(data_qubit_error)
    errors = {"blue": {(11, 1), (5, 1)}, "red": {(8, 0), (14, 0)}, "green": {(8, 2)}}
    correction = decoder.single_run(error_cor=errors)
    assert correction == {(12, 0), (7, 1)}


def test_corner():
    errors = {"blue": {(11, 7)}, "red": set(), "green": set()}
    correction = decoder.single_run(error_cor=errors)
    assert correction == {(7, 7), (9, 7)}


def test_single_run():
    correction = decoder.single_run()
    assert correction == {(10, 6), (16, 0), (18, 0)}


def test_create_red_green_matching():
    matching = decoder.create_red_green_matching({(2, 0), (5, 3)})
    assert set(matching.nodes()) == {(2, 0), (2, 2), (5, 3)}


def test_create_green_blue_matching():
    matching = decoder.create_blue_green_matching({(14, 4), (11, 1)})
    assert set(matching.nodes()) == {(11, 1), (14, 2), (14, 4)}


def test_create_graph():
    errors_red = [(8, 6), (5, 3), (2, 0), (14, 0)]
    errors_green = [(8, 2), (11, 5)]

    error_cor = errors_red + errors_green
    node_num, edge_num, edges, node2id = decoder.create_graph(
        error_cor, decoder.matching_graph.green_red, ["vR", "vG"]
    )
    assert node_num == 12
    assert edge_num == 36


def test_errors_figure_1():
    """figure 1 from chamberlands paper"""
    errors_blue = [(11, 7)]
    errors_red = [(8, 6), (5, 3), (2, 0), (14, 0)]
    errors_green = [(8, 2), (11, 5)]

    error_cor = errors_red + errors_green

    node_num, edge_num, edges, node2id = decoder.create_graph(
        error_cor, decoder.matching_graph.green_red, ["vR", "vG"]
    )
    matching, cmatching = decoder.run_blossom(node_num, edge_num, edges)
    graph = decoder.translate_cmatching(
        matching, cmatching, node2id, decoder.matching_graph.green_red, ["vR", "vG"]
    )

    assert graph.has_edge((8, 6), (11, 5)) is True
    assert graph.has_edge((5, 3), (8, 2)) is True
    assert graph.has_edge((2, 0), "vG") is True
    assert graph.has_edge((14, 0), "vG") is True

    error_cor = errors_red + errors_blue
    node_num, edge_num, edges, node2id = decoder.create_graph(
        error_cor, decoder.matching_graph.red_blue, ["vR", "vB"]
    )
    matching, cmatching = decoder.run_blossom(node_num, edge_num, edges)
    graph = decoder.translate_cmatching(
        matching, cmatching, node2id, decoder.matching_graph.red_blue, ["vR", "vB"]
    )
    assert graph.has_edge((8, 6), (11, 7)) is True
    assert graph.has_edge((5, 3), ("vB")) is True
    assert graph.has_edge((2, 0), ("vB")) is True
    assert graph.has_edge((14, 0), (17, 1)) is True
    assert graph.has_edge(("vR"), (17, 1)) is True

    error_cor = errors_green + errors_blue
    node_num, edge_num, edges, node2id = decoder.create_graph(
        error_cor, decoder.matching_graph.blue_green, ["vB", "vG"]
    )
    matching, cmatching = decoder.run_blossom(node_num, edge_num, edges)
    graph = decoder.translate_cmatching(
        matching, cmatching, node2id, decoder.matching_graph.blue_green, ["vB", "vG"]
    )

    assert graph.has_edge((8, 2), (11, 1)) is True
    assert graph.has_edge(("vG"), (11, 1)) is True
    assert graph.has_edge((11, 7), (11, 5)) is True


def test_neighbors_not_in_rho():
    vertex = (11, 5)
    neighbors = []
    n_not_in_rho = decoder.find_neighbors_not_in_rho(vertex, neighbors)
    assert n_not_in_rho == {(11, 7), (14, 4), (11, 3), (8, 4), (8, 6)}


def test_lift_non_boundary():
    nodes = (11, 5)
    neighbors = [(8, 6), (11, 7)]
    correction = decoder.lift(nodes, neighbors)
    assert correction == {(10, 6)}

    nodes = (11, 3)
    neighbors = [(11, 5), (8, 4), (8, 2), (11, 1)]
    correction = decoder.lift(nodes, neighbors)
    assert correction == {(10, 4), (10, 2)}

    nodes = (17, 1)
    neighbors = [(14, 0), (20, 0)]
    correction = decoder.lift(nodes, neighbors)
    assert correction == {(16, 0), (18, 0)}


def test_translate_boundary_neighbors():
    node = (17, 1)
    neighbor = ["vR"]
    G = nx.Graph()
    G.add_node("vR")

    boundary_neighbors = decoder.translate_boundary_neighbors(neighbor, node, [])
    assert boundary_neighbors == [(20, 0)]

    node = (17, 1)
    neighbor = ["vG"]
    G = nx.Graph()
    G.add_node("vG")

    boundary_neighbors = decoder.translate_boundary_neighbors(neighbor, node, [])
    assert boundary_neighbors == [(17, -1)]

    G2 = nx.Graph()
    node = (2, 2)
    neighbor = ["vB"]
    G2.add_node("vB")

    boundary_neighbors = decoder.translate_boundary_neighbors(neighbor, node, [])
    assert boundary_neighbors == [(-1, 1)]

    G3 = nx.Graph()
    node = (8, 0)
    neighbor = ["vG"]
    G3.add_node("vG")

    boundary_neighbors = decoder.translate_boundary_neighbors(neighbor, node, [])

    assert boundary_neighbors == [(11, -1)]


def test_single_distance_5():
    data_qubit_error = {(9, 1), (6, 0), (7, 3)}
    ancilla_error_cords = decoder.layout.data_qubits_to_ancilla_qubits(data_qubit_error)
    ancilla_error = {
        "blue": ancilla_error_cords.intersection(decoder.layout.ancilla_qubits_blue),
        "red": ancilla_error_cords.intersection(decoder.layout.ancilla_qubits_red),
        "green": ancilla_error_cords.intersection(decoder.layout.ancilla_qubits_green),
    }
    correction = decoder_5.single_run(ancilla_error)
    ancilla_correction_cords = decoder.layout.data_qubits_to_ancilla_qubits(
        data_qubit_error
    )
    ancilla_correction = {
        "blue": ancilla_error_cords.intersection(decoder.layout.ancilla_qubits_blue),
        "red": ancilla_error_cords.intersection(decoder.layout.ancilla_qubits_red),
        "green": ancilla_error_cords.intersection(decoder.layout.ancilla_qubits_green),
    }

    assert ancilla_correction == ancilla_error


def test_corner_correction():
    ancilla_error = {"blue": {(11, 1)}, "red": set(), "green": set()}
    correction = decoder_5.single_run(ancilla_error)
    assert correction == {(12, 0)}


def test_run_3():
    ancilla_error = {
        "blue": {(11, 1), (8, 4), (5, 1)},
        "red": {(2, 0), (8, 0)},
        "green": {(8, 2)},
    }
    correction = decoder_5.single_run(ancilla_error)
    print(correction, "cor")


def test_run_4():
    data_qubit_error = {(6, 2), (9, 3), (12, 0)}
    ancilla_error_cords = decoder_5.layout.data_qubits_to_ancilla_qubits(
        data_qubit_error
    )
    print(ancilla_error_cords, "ancilla error coordinates")

    ancilla_error = {"blue": {(11, 1), (8, 4), (5, 1)}, "red": {(5, 3)}, "green": set()}
    correction = decoder_5.single_run(ancilla_error)
    print(correction, "cor")
    ancilla_error_cords = decoder_5.layout.data_qubits_to_ancilla_qubits(correction)
    print(ancilla_error_cords, "ancilla error coordinates")


# test_run_4()


def test_run_5():
    ancilla_error = {
        "blue": {(11, 1), (8, 4), (5, 1)},
        "red": {(5, 3)},
        "green": {(2, 2)},
    }
    correction = decoder_5.single_run(ancilla_error)
    ancilla_error_cords = decoder_5.layout.data_qubits_to_ancilla_qubits(correction)
    print(ancilla_error_cords, "ancilla error cords")

    print("\n \n distance 7 \n \n \n")
    correction = decoder.single_run(ancilla_error)
    ancilla_error_cords = decoder.layout.data_qubits_to_ancilla_qubits(correction)
    print(ancilla_error_cords, "ancilla error cords")


# test_run_5()


def test_run_6():
    ancilla_error = {"blue": {(5, 1)}, "red": set(), "green": {(8, 2), (5, 5), (2, 2)}}
    correction = decoder_5.single_run(ancilla_error)
    ancilla_error_cords = decoder_5.layout.data_qubits_to_ancilla_qubits(correction)
    print(correction, "cor")
    print(ancilla_error_cords)


def test_run_7():
    ancilla_error = {"blue": {(11, 7), (17, 1)}, "red": set(), "green": {(5, 5)}}
    correction = decoder.single_run(ancilla_error)
    # correction.remove((9, 7))
    ancilla_error_cords = decoder.layout.data_qubits_to_ancilla_qubits(correction)
    print(correction, "cor")

    print(ancilla_error_cords, "ancilla")

    # remove duplicates from correction
