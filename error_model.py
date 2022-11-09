from typing import Set, Tuple
import numpy as np
from layout import Hexagonal_layout


class BiasedNoiseModel:
    """
    Noise model DW(2/3, -pi/4)
    """

    def __init__(
        self, error_probability: float, bias: int, layout: Hexagonal_layout
    ):  # data_qubit_set: Set[Tuple[int, int]]):
        """

            Probality for X,Y,Z is p/3, so probability for error causing an
        ancilla qubit to fire 2p/3
        """
        self.error_probability = error_probability
        self.bias = bias
        self.layout = layout
        self.error_probability_dict = self.create_error_probabilities()

    def create_random_error(self) -> set:
        # do this!

        error_cords = set(
            qubit_coordinates
            for qubit_coordinates, probability in self.error_probability_dict.items()
            if np.random.rand() < probability
        )
        # yeah this is a data qubit error of cours
        #      error_cords = set(i for i in self.qubits if int_sample(self.p_array) == 1)
        return error_cords

    def get_data_qubits_to_flip(self):
        data_qubits_to_flip = set()
        flip_parity = 0
        for index, i in enumerate(range(0, 2 * self.layout.distance + 5, 2)):
            if index > 0 and (index - 1) // 3 % 2 == 0:
                flip_parity = 1
            else:
                flip_parity = 0
            for j in range(0, self.layout.distance * 2):
                if (i + j, j) in self.layout.data_qubits:
                    if flip_parity == 1:
                        data_qubits_to_flip.add((i + j, j))
                        flip_parity = 0
                    else:
                        flip_parity = 1
        return data_qubits_to_flip

    def create_error_probabilities(self):
        error_probability_dict = dict()
        if self.bias == 1:
            for qubit in self.layout.data_qubits:
                error_probability_dict[qubit] = 2 / 3 * self.error_probability
        else:
            pzy = (
                self.bias / (self.bias + 1) + 1 / (2 * self.bias + 2)
            ) * self.error_probability  # pz+py
            pxy = 2 / (2 * self.bias + 2) * self.error_probability  # px+py

            data_qubits_to_flip = self.get_data_qubits_to_flip()
            for qubit in self.layout.data_qubits:
                if qubit in data_qubits_to_flip:
                    error_probability_dict[qubit] = pzy
                else:
                    error_probability_dict[qubit] = pxy
#        print(error_probability_dict)
        return error_probability_dict


class ColourCodeXErrorModel(object):
    def __init__(self, p_code_capacity: int, data_qubit_set: Set[Tuple[int, int]]):
        """
        Probality for X,Y,Z is p/3, so probability for error causing an
        ancilla qubit to fire 2p/3
        """
        error_prob = 2 * p_code_capacity / 3
        # error_prob = 2*p_code_capacity/3

        # p_array[0] is idling + Z prob, p_array[1] is X+Y prob
        self.p_array = [1 - error_prob, error_prob]
        self.qubits = data_qubit_set

    def create_random_error(self) -> set:
        error_cords = set(i for i in self.qubits if int_sample(self.p_array) == 1)
        return error_cords


# def create_p_array():

"""
"""


def int_sample(probs):
    """
    Code from Ben Criger:
    (TODO: check in numpy has something fast to do this)

    This is a little fancy, so I'll explain it.
    If we are given an array of probabilities [p_0, p_1, ..., p_n],
    we can sample, yielding k if a unit uniform deviate is
    within the interval (sum_{k<j}p_k, sum_{k<j}p_k + p_k). To
    accomplish this, we first take such a sample, then subtract off
    p_k's as we proceed. In theory, distributions which are sorted
    descending will be most efficiently sampled, but in practice it
    doesn't make a difference. We explicitly don't sort.
    """
    value = np.random.rand()
    # value = np.random.seed(random.randint(10000))

    for idx, p in enumerate(probs):
        if value < p:
            return idx
        else:
            value -= p

    raise ValueError("Enumeration Failure (probability > 1?)")


# """
