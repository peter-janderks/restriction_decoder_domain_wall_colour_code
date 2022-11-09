from error_model import BiasedNoiseModel
from layout import Hexagonal_layout
from minimum_weight_perfect_matching_pymatching import MWPM


class Calculate_LER(object):
    """
    Simulating the triangular color code on a trivalent graph with flag qubits
    Uses the code capacity noise model and the adapted restriction decoder for
    2D decoding.
    """

    def __init__(self, distance, per, bias):
        """
        distance = distance of the qec code
        per = physical error rate
        """
        self.layout = Hexagonal_layout(distance)
        self.error_model = BiasedNoiseModel(per, bias, self.layout)
        self.decoder = MWPM(distance, self.layout, self.error_model)

    def run(self, trials, max_logicals):
        logical_errors = 0
        # error = self.create_error()
        for n_runs in range(trials):
            data_error, ancilla_error = self.create_error()
            correction = self.decoder.single_run(ancilla_error)
            resulting_operator = correction ^ data_error
            ancilla_error_cords = self.layout.data_qubits_to_ancilla_qubits(correction)
            syndrome_correction = {
                "blue": ancilla_error_cords.intersection(
                    self.layout.ancilla_qubits_blue
                ),
                "red": ancilla_error_cords.intersection(self.layout.ancilla_qubits_red),
                "green": ancilla_error_cords.intersection(
                    self.layout.ancilla_qubits_green
                ),
            }

            if ancilla_error != syndrome_correction:
                print(data_error, "data error")
                print(ancilla_error, "ancilaa")
                print(syndrome_correction, "syndrome")
                stop
                #                print(error)

            logical_errors = self.check_for_logical_error(
                resulting_operator, logical_errors
            )
            if logical_errors == max_logicals:
                break
            # if len(resulting_operator) % 2 == 1:
            #    logical_errors += 1

        return (logical_errors, n_runs)

    def check_for_logical_error(self, resulting_operator, logical_errors):
        """
        error + correction is one of the following three assuming we are in the codespace:
        1. Logical error
        2. Stabilizer (identity)
        3. Null

        stabilizer operators always have even parity, therefore we
        can check for a logical error by just checking the parity
        """
        difference = self.layout.logical_operator - resulting_operator

        if len(difference) % 2 == 0:
            logical_errors += 1
            # print('logical?')
            # print(resulting_operator, 'r op')
            # print(difference, 'difference')
        return logical_errors

    def create_error(self):
        """
        calls the error model and splits the error into colours!
        """
        data_error_cords = self.error_model.create_random_error()
        ancilla_error_cords = self.layout.data_qubits_to_ancilla_qubits(
            data_error_cords
        )
        ancilla_error = {
            "blue": ancilla_error_cords.intersection(self.layout.ancilla_qubits_blue),
            "red": ancilla_error_cords.intersection(self.layout.ancilla_qubits_red),
            "green": ancilla_error_cords.intersection(self.layout.ancilla_qubits_green),
        }

        return (data_error_cords, ancilla_error)
