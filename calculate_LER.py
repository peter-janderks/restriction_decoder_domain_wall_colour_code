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
        logical_errors_x = 0
        logical_errors_z = 0
        total_logical_errors = 0
        for n_runs in range(1, trials):
            (
                data_error_X,
                data_error_Z,
                ancilla_error_X,
                ancilla_error_Z,
            ) = self.create_error()

            correction_X = self.decoder.single_run(
                self.decoder.matching_graph_X,
                self.decoder.coords_matching_graph_X,
                ancilla_error_X,
            )
            resulting_operator_X = correction_X ^ data_error_X
            ancilla_error_cords_X = self.layout.data_qubits_to_ancilla_qubits(
                correction_X
            )
            syndrome_correction_X = {
                "blue": ancilla_error_cords_X.intersection(
                    self.layout.ancilla_qubits_blue
                ),
                "red": ancilla_error_cords_X.intersection(
                    self.layout.ancilla_qubits_red
                ),
                "green": ancilla_error_cords_X.intersection(
                    self.layout.ancilla_qubits_green
                ),
            }

            if ancilla_error_X != syndrome_correction_X:
                print(data_error_X, "data error")
                print(ancilla_error_X, "ancilaa")
                print(syndrome_correction_X, "syndrome")

            correction_Z = self.decoder.single_run(
                self.decoder.matching_graph_Z,
                self.decoder.coords_matching_graph_Z,
                ancilla_error_Z,
            )
            resulting_operator_Z = correction_Z ^ data_error_Z
            ancilla_error_cords_Z = self.layout.data_qubits_to_ancilla_qubits(
                correction_Z
            )
            syndrome_correction_Z = {
                "blue": ancilla_error_cords_Z.intersection(
                    self.layout.ancilla_qubits_blue
                ),
                "red": ancilla_error_cords_Z.intersection(
                    self.layout.ancilla_qubits_red
                ),
                "green": ancilla_error_cords_Z.intersection(
                    self.layout.ancilla_qubits_green
                ),
            }

            if ancilla_error_Z != syndrome_correction_Z:
                print(data_error, "data error")
                print(ancilla_error, "ancilaa")
                print(syndrome_correction, "syndrome")

            # to update check for logical error

            lx_bool = self.check_for_logical_error(resulting_operator_X)
            lz_bool = self.check_for_logical_error(resulting_operator_Z)
            if lx_bool == True or lz_bool == True:

                total_logical_errors += 1
                if lx_bool == True:

                    logical_errors_x += 1

                if lz_bool == True:
                    logical_errors_z += 1

            if total_logical_errors == max_logicals:
                break
            # if len(resulting_operator) % 2 == 1:
            #    logical_errors += 1
        return (total_logical_errors, logical_errors_x, logical_errors_z, n_runs)

    def check_for_logical_error(self, resulting_operator):
        """
        error + correction is one of the following three assuming we are in the codespace:
        1. Logical error
        2. Stabilizer (identity)
        3. Null

        stabilizer operators always have even parity, therefore we
        can check for a logical error by just checking the parity

        TODO: implement qecsim check
        """
        #        difference = self.layout.logical_operator - resulting_operator
        if len(resulting_operator) % 2 == 1:
            return True
        else:
            return False
            # print('logical?')
            # print(resulting_operator, 'r op')
            # print(difference, 'difference')

    def create_error(self):
        """
        calls the error model and splits the error into colours!
        """
        data_error_cords_X, data_error_cords_Z = self.error_model.create_random_error()
        ancilla_error_cords_Z = self.layout.data_qubits_to_ancilla_qubits(
            data_error_cords_Z
        )

        ancilla_error_cords_X = self.layout.data_qubits_to_ancilla_qubits(
            data_error_cords_X
        )
        ancilla_error_X = {
            "blue": ancilla_error_cords_X.intersection(self.layout.ancilla_qubits_blue),
            "red": ancilla_error_cords_X.intersection(self.layout.ancilla_qubits_red),
            "green": ancilla_error_cords_X.intersection(
                self.layout.ancilla_qubits_green
            ),
        }

        ancilla_error_Z = {
            "blue": ancilla_error_cords_Z.intersection(self.layout.ancilla_qubits_blue),
            "red": ancilla_error_cords_Z.intersection(self.layout.ancilla_qubits_red),
            "green": ancilla_error_cords_Z.intersection(
                self.layout.ancilla_qubits_green
            ),
        }
        return (
            data_error_cords_X,
            data_error_cords_Z,
            ancilla_error_X,
            ancilla_error_Z,
        )
