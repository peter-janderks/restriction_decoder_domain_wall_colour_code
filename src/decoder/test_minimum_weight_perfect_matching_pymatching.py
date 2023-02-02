import pytest
from create_matching_graph import Matching_graph
from minimum_weight_perfect_matching_pymatching import MWPM
import networkx as nx
from utils.error_model import BiasedNoiseModel
from utils.layout import Hexagonal_layout


layout5 = Hexagonal_layout(5)
layout7 = Hexagonal_layout(7)

g5biased = Matching_graph(5, BiasedNoiseModel(0.3, 1000, layout5))

decoder5 = MWPM(
    5,
    layout5,
    BiasedNoiseModel(0.1, 0.5, layout5),
)

decoder7 = MWPM(
    7,
    layout7,
    BiasedNoiseModel(0.1, 0.5, layout7),
)

decoder5_biased = MWPM(
    5,
    layout5,
    BiasedNoiseModel(0.3, 1000, layout5),
)




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


def test_weight_2_error_d7():
    data_qubit_error_X = {(6, 6), (10, 4)}
    data_qubit_error_Z = {(6, 6)}
    (
        ancilla_error_X,
        ancilla_error_Z,
    ) = data_qubit_errors_to_ancilla(data_qubit_error_X, data_qubit_error_Z, decoder7)
    correction_Z = decoder7.single_run(
        decoder7.matching_graph_Z, decoder7.coords_matching_graph_Z, ancilla_error_Z
    )
    assert correction_Z == {(6, 6)}

    correction_X = decoder7.single_run(
        decoder7.matching_graph_X, decoder7.coords_matching_graph_X, ancilla_error_X
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
    error_model = BiasedNoiseModel(0.2, 2, layout5)

    decoder5 = MWPM(
        5,
        layout5,
        error_model
    )
    data_qubit_error_X =  {(3,1),(4,2)}
    data_qubit_error_Z = set()
    (
        ancilla_error_X,
        ancilla_error_Z,
    ) = data_qubit_errors_to_ancilla(data_qubit_error_X, data_qubit_error_Z,decoder5)

    correction_Z = decoder5.single_run(
        decoder5.matching_graph_Z, decoder5.coords_matching_graph_Z, ancilla_error_Z
    )
    resulting_operator_Z = data_qubit_error_Z ^ correction_Z

    correction_X = decoder5.single_run(
        decoder5.matching_graph_X, decoder5.coords_matching_graph_X, ancilla_error_X
    )

    resulting_operator_X = correction_X ^ data_qubit_error_X
    assert len(resulting_operator_X) % 2 == 0
    assert len(resulting_operator_Z) % 2 == 0


def test_weight_two_error():

    data_qubit_error_X = {(15, 1), (16, 0)}
    data_qubit_error_Z = {(12, 0), (10, 0)}
    (
        ancilla_error_X,
        ancilla_error_Z,
    ) = data_qubit_errors_to_ancilla(data_qubit_error_X, data_qubit_error_Z, decoder7)
    correction_Z = decoder7.single_run(
        decoder7.matching_graph_Z, decoder7.coords_matching_graph_Z, ancilla_error_Z
    )   
    resulting_operator_Z = correction_Z ^ data_qubit_error_Z
    assert len(resulting_operator_Z)%2 == 0

    correction_X = decoder7.single_run(
        decoder7.matching_graph_X, decoder7.coords_matching_graph_X, ancilla_error_X
    )
    resulting_operator_X = correction_X ^ data_qubit_error_X
    assert len(resulting_operator_X)%2 == 0

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
    assert len(resulting_operator_Z)%2 == 0

    correction_X = decoder5.single_run(
        decoder5.matching_graph_X, decoder5.coords_matching_graph_X, ancilla_error_X
    )
    resulting_operator_X = correction_X ^ data_qubit_error_X

    assert len(resulting_operator_X)%2 == 0

def test_weight_three_error():
    data_qubit_error = {(16, 2), (10, 0), (7, 1)}
    ancilla_error_cords = decoder7.layout.data_qubits_to_ancilla_qubits(data_qubit_error)

    ancilla_error = {
        "blue": ancilla_error_cords.intersection(decoder7.layout.ancilla_qubits_blue),
        "red": ancilla_error_cords.intersection(decoder7.layout.ancilla_qubits_red),
        "green": ancilla_error_cords.intersection(decoder7.layout.ancilla_qubits_green),
    }
    correction = decoder7.single_run(decoder7.matching_graph_Z, decoder7.coords_matching_graph_Z,ancilla_error)
    resulting_operator = correction ^ data_qubit_error

    assert len(resulting_operator) == 0

    data_qubit_error = {(6, 2), (7, 7), (9, 3)}
    ancilla_error_cords = decoder7.layout.data_qubits_to_ancilla_qubits(data_qubit_error)
    ancilla_error = {
        "blue": ancilla_error_cords.intersection(decoder7.layout.ancilla_qubits_blue),
        "red": ancilla_error_cords.intersection(decoder7.layout.ancilla_qubits_red),
        "green": ancilla_error_cords.intersection(decoder7.layout.ancilla_qubits_green),
    }
    correction = decoder7.single_run(decoder7.matching_graph_Z, decoder7.coords_matching_graph_Z,ancilla_error)

    resulting_operator = correction ^ data_qubit_error
    assert len(resulting_operator) == 0

    data_qubit_error = {(6, 6), (6, 4), (4, 2)}
    ancilla_error_cords = decoder7.layout.data_qubits_to_ancilla_qubits(data_qubit_error)

    ancilla_error = {
        "blue": ancilla_error_cords.intersection(decoder7.layout.ancilla_qubits_blue),
        "red": ancilla_error_cords.intersection(decoder7.layout.ancilla_qubits_red),
        "green": ancilla_error_cords.intersection(decoder7.layout.ancilla_qubits_green),
    }
    correction = decoder7.single_run(decoder7.matching_graph_Z, decoder7.coords_matching_graph_Z,ancilla_error)
    resulting_operator = correction ^ data_qubit_error
    assert len(resulting_operator) % 2== 0


    errors_blue = set({(11, 7)})
    errors_red = set({(8, 6)})
    errors_green = set({(11, 5)})

    errors = {"blue": errors_blue, "red": errors_red, "green": errors_green}
    correction = decoder7.single_run(decoder7.matching_graph_Z, decoder7.coords_matching_graph_Z,errors)

    assert correction == {(10, 6)}

    errors = {"blue": set(), "red": {(5, 3)}, "green": {(5, 5)}}
    correction = decoder7.single_run(decoder7.matching_graph_Z, decoder7.coords_matching_graph_Z,errors)
 #   correction = decoder7.single_run(error_cor=errors)
    assert correction == {(4, 4)}

    errors = {"blue": set(), "green": {(5, 5)}, "red": {(8, 6)}}
    correction = decoder7.single_run(decoder7.matching_graph_Z, decoder7.coords_matching_graph_Z,errors)
#    correction = decoder7.single_run(error_cor=errors)
    assert correction == {(6, 6)}

    errors = {"blue": {(11, 1)}, "red": {(8, 0)}, "green": set()}
    correction = decoder7.single_run(decoder7.matching_graph_Z, decoder7.coords_matching_graph_Z,errors)
#    correction = decoder7.single_run(error_cor=errors)
    assert correction == {(10, 0)}

    errors = {"blue": {(14, 4)}, "red": set(), "green": {(14, 2)}}
    correction = decoder7.single_run(decoder7.matching_graph_Z, decoder7.coords_matching_graph_Z,errors)
#    correction = decoder7.single_run(error_cor=errors)
    assert correction == {(15, 3)}

    errors = {"blue": {(14, 4)}, "red": set(), "green": {(11, 5)}}
    correction = decoder7.single_run(decoder7.matching_graph_Z, decoder7.coords_matching_graph_Z,errors)
    assert correction == {(13, 5)}

    errors = {"blue": {(11, 7)}, "red": set(), "green": {(11, 5)}}
    correction = decoder7.single_run(decoder7.matching_graph_Z, decoder7.coords_matching_graph_Z,errors)
    assert correction == {(12, 6)}

    errors = {"blue": {(5, 1)}, "red": {(8, 0)}, "green": set()}
    correction = decoder7.single_run(decoder7.matching_graph_Z, decoder7.coords_matching_graph_Z,errors)
    assert correction == {(6, 0)}

    errors = {"blue": set(), "red": {(2, 0)}, "green": {(2, 2)}}
    correction = decoder7.single_run(decoder7.matching_graph_Z, decoder7.coords_matching_graph_Z,errors)
    assert correction == {(1, 1)}

    errors = {"blue": {(11, 1)}, "red": {(8, 0)}, "green": set()}
    correction = decoder7.single_run(decoder7.matching_graph_Z, decoder7.coords_matching_graph_Z,errors)
    assert correction == {(10, 0)}

    errors = {"blue": {(11, 1)}, "red": {(14, 0)}, "green": set()}
    correction = decoder7.single_run(decoder7.matching_graph_Z, decoder7.coords_matching_graph_Z,errors)
    assert correction == {(12, 0)}

    errors = {"blue": {(5, 1)}, "red": {(8, 0)}, "green": {(8, 2)}}
    correction = decoder7.single_run(decoder7.matching_graph_Z, decoder7.coords_matching_graph_Z,errors)
    assert correction == {(7, 1)}

    data_qubit_error = {(12, 0), (7, 1)}
    ancilla_error_cords = decoder7.layout.data_qubits_to_ancilla_qubits(data_qubit_error)
    errors = {"blue": {(11, 1), (5, 1)}, "red": {(8, 0), (14, 0)}, "green": {(8, 2)}}
    correction = decoder7.single_run(decoder7.matching_graph_Z, decoder7.coords_matching_graph_Z,errors)
    assert correction == {(12, 0), (7, 1)}

def test_corner():
    errors = {"blue": {(11, 7)}, "red": set(), "green": set()}
    correction = decoder7.single_run(decoder7.matching_graph_Z, decoder7.coords_matching_graph_Z,errors)
    assert correction == {(7, 7), (9, 7)}

def test_neighbors_not_in_rho():
    vertex = (11, 5)
    neighbors = []
    n_not_in_rho = decoder7.find_neighbors_not_in_rho(vertex, neighbors)
    assert n_not_in_rho == {(11, 7), (14, 4), (11, 3), (8, 4), (8, 6)}


def test_lift_non_boundary():
    nodes = (11, 5)
    neighbors = [(8, 6), (11, 7)]
    correction = decoder7.lift(nodes, neighbors)
    assert correction == {(10, 6)}

    nodes = (11, 3)
    neighbors = [(11, 5), (8, 4), (8, 2), (11, 1)]
    correction = decoder7.lift(nodes, neighbors)
    assert correction == {(10, 4), (10, 2)}

    nodes = (17, 1)
    neighbors = [(14, 0), (20, 0)]
    correction = decoder7.lift(nodes, neighbors)
    assert correction == {(16, 0), (18, 0)}


def test_translate_boundary_neighbors():
    node = (17, 1)
    neighbor = ["vR"]
    G = nx.Graph()
    G.add_node("vR")

    boundary_neighbors = decoder7.translate_boundary_neighbors(neighbor, node, [])
    assert boundary_neighbors == [(20, 0)]

    node = (17, 1)
    neighbor = ["vG"]
    G = nx.Graph()
    G.add_node("vG")

    boundary_neighbors = decoder7.translate_boundary_neighbors(neighbor, node, [])
    assert boundary_neighbors == [(17, -1)]

    G2 = nx.Graph()
    node = (2, 2)
    neighbor = ["vB"]
    G2.add_node("vB")

    boundary_neighbors = decoder7.translate_boundary_neighbors(neighbor, node, [])
    assert boundary_neighbors == [(-1, 1)]

    G3 = nx.Graph()
    node = (8, 0)
    neighbor = ["vG"]
    G3.add_node("vG")

    boundary_neighbors = decoder7.translate_boundary_neighbors(neighbor, node, [])

    assert boundary_neighbors == [(11, -1)]