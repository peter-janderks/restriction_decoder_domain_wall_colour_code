from error_model import BiasedNoiseModel
from layout import Hexagonal_layout

d7layout = Hexagonal_layout(7)
d5layout = Hexagonal_layout(5)
d3layout = Hexagonal_layout(3)
d9layout = Hexagonal_layout(9)
d11layout = Hexagonal_layout(11)
d13layout = Hexagonal_layout(13)
d15layout = Hexagonal_layout(15)
d17layout = Hexagonal_layout(17)
d19layout = Hexagonal_layout(19)


class BiasedNoiseModelMock(BiasedNoiseModel):
    def __init__(self, p_code_capicity: float, bias: float, layout: Hexagonal_layout):
        self.error_probability = p_code_capicity
        self.layout = layout
        self.bias = bias


def test_create_random_error():
    BiasedNoiseModel(0.1, 10, d5layout)
    pass


# test_create_random_error()


def calc_ratio():
    for i in range(101, 110, 2):
        layout = Hexagonal_layout(i)
        noise_model = BiasedNoiseModelMock(0.1, 1, layout)
        flip_data_qubits = noise_model.get_data_qubits_to_flip()
        print(len(flip_data_qubits) / len(layout.data_qubits), "layout", i)


# calc_ratio()


def get_data_qubits_to_flip():
    noise_model = BiasedNoiseModelMock(0.1, 1, d5layout)
    flip_data_qubits = noise_model.get_data_qubits_to_flip()
    assert flip_data_qubits == {
        (1, 1),
        (3, 1),
        (4, 0),
        (6, 0),
        (4, 4),
        (6, 4),
        (7, 3),
        (9, 3),
        (10, 2),
    }

    noise_model = BiasedNoiseModelMock(0.1, 1, d3layout)
    flip_data_qubits = noise_model.get_data_qubits_to_flip()
    assert flip_data_qubits == {
        (1, 1),
        (3, 1),
        (4, 0),
        (6, 0),
    }

    noise_model = BiasedNoiseModelMock(0.1, 1, d7layout)
    flip_data_qubits = noise_model.get_data_qubits_to_flip()
    assert flip_data_qubits == {
        (1, 1),
        (3, 1),
        (4, 0),
        (6, 0),
        (4, 4),
        (6, 4),
        (7, 3),
        (9, 3),
        (10, 2),
        (7, 7),
        (9, 7),
        (10, 6),
        (12, 6),
        (13, 5),
        (12, 2),
        (13, 1),
        (15, 1),
        (16, 0),
        (18, 0),
    }

get_data_qubits_to_flip()

def test_create_error_probabilities():
    noise_model = BiasedNoiseModelMock(0.1, 10, d3layout)
    (
        error_probability_dict_X,
        error_probability_dict_Y,
        error_probability_dict_Z,
        _,
        _
    ) = noise_model.create_error_probabilities()
    pz = 10 / (10 + 1) *0.1
    px = 1 / (2 * 10 + 2) * 0.1
    py = 1 / (2 * 10 + 2) * 0.1
    comparison_probability_dict_X = {
        (0, 0): px,
        (4, 0): pz,
        (6, 0): pz,
        (1, 1): pz,
        (3, 1): pz,
        (3, 3): px,
        (4, 2): px,
    }

    comparison_probability_dict_Y = {
        (0, 0): py,
        (4, 0): py,
        (6, 0): py,
        (1, 1): py,
        (3, 1): py,
        (3, 3): py,
        (4, 2): py,
    }

    comparison_probability_dict_Z = {
        (0, 0): pz,
        (4, 0): px,
        (6, 0): px,
        (1, 1): px,
        (3, 1): px,
        (3, 3): pz,
        (4, 2): pz,
    }
    assert error_probability_dict_X == comparison_probability_dict_X
    assert error_probability_dict_Y == comparison_probability_dict_Y
    assert error_probability_dict_Z == comparison_probability_dict_Z
   
test_create_error_probabilities()

def test_create_random_error():
    for _ in range(10):
        noise_model = BiasedNoiseModel(0.1, 1000, d3layout)
        error_X,error_Z = noise_model.create_random_error()
        print(error_X, error_Z)

