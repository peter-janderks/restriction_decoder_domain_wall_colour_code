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
        print(len(flip_data_qubits) / len(layout.data_qubits), "layout",i)


calc_ratio()


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

    noise_model = BiasedNoiseModelMock(0.1, 1, d7layout)
    flip_data_qubits = noise_model.get_data_qubits_to_flip()
    print(len(flip_data_qubits) / len(d7layout.data_qubits))

    noise_model = BiasedNoiseModelMock(0.1, 1, d9layout)
    flip_data_qubits = noise_model.get_data_qubits_to_flip()
    print(len(flip_data_qubits) / len(d9layout.data_qubits))

    noise_model = BiasedNoiseModelMock(0.1, 1, d11layout)
    flip_data_qubits = noise_model.get_data_qubits_to_flip()
    print(len(flip_data_qubits) / len(d11layout.data_qubits))

    noise_model = BiasedNoiseModelMock(0.1, 1, d13layout)
    flip_data_qubits = noise_model.get_data_qubits_to_flip()
    print(len(flip_data_qubits) / len(d13layout.data_qubits))

    noise_model = BiasedNoiseModelMock(0.1, 1, d15layout)
    flip_data_qubits = noise_model.get_data_qubits_to_flip()
    print(len(flip_data_qubits) / len(d15layout.data_qubits))

    noise_model = BiasedNoiseModelMock(0.1, 1, d17layout)
    flip_data_qubits = noise_model.get_data_qubits_to_flip()
    print(len(flip_data_qubits) / len(d17layout.data_qubits))

    noise_model = BiasedNoiseModelMock(0.1, 1, d19layout)
    flip_data_qubits = noise_model.get_data_qubits_to_flip()
    print(len(flip_data_qubits) / len(d19layout.data_qubits))


# get_data_qubits_to_flip()


def test_create_error_probabilities():
    noise_model = BiasedNoiseModelMock(0.1, 10, d3layout)
    error_probability_dict = noise_model.create_error_probabilities()
    pzy = (10 / (10 + 1) + 1 / (2 * 10 + 2)) * 0.1
    pxy = 2 / (2 * 10 + 2) * 0.1

    comparison_probability_dict = {
        (0, 0): pxy,
        (4, 0): pzy,
        (6, 0): pzy,
        (1, 1): pzy,
        (3, 1): pzy,
        (3, 3): pxy,
        (4, 2): pxy,
    }
    assert error_probability_dict == comparison_probability_dict


def test_create_random_error():
    noise_model = BiasedNoiseModel(0.1, 10, d3layout)
    error = noise_model.create_random_error()
