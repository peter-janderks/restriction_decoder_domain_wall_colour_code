class Hexagonal_layout:
    def __init__(self, distance):
        self.ancilla_indexes(distance)

    def ancilla_indexes(self, distance):
        x_max = distance + distance // 2
        if distance > 7:
            ValueError("not imremented yet, need to add corner qubits!")

        self.distance = distance
        self.ancilla_qubits = set()
        self.data_qubits = set()

        self.ancilla_qubits_red = set()
        self.ancilla_qubits_blue = set()
        self.ancilla_qubits_green = set()

        self.ancilla_red_boundary = set()
        self.ancilla_blue_boundary = set()
        self.ancilla_green_boundary = set()

        self.logical_operator = set()

        colour = ["r", "b", "g"]
        y = 0
        while x_max > 0:
            x_row = x_max
            ancilla_colour = colour[y % 3]
            if ancilla_colour == "r":

                self.red_row(x_max, y)

            elif ancilla_colour == "b":

                self.blue_row(x_max, y)
                if y == 1:
                    self.ancilla_green_boundary = self.ancilla_qubits.copy()
            else:
                self.green_row(x_max, y)
            x_max -= 1
            y += 1

        self.ancilla_qubits_red_green = self.ancilla_qubits_red.union(
            self.ancilla_qubits_green
        )
        self.ancilla_qubits_red_blue = self.ancilla_qubits_red.union(
            self.ancilla_qubits_blue
        )
        self.ancilla_qubits_blue_green = self.ancilla_qubits_blue.union(
            self.ancilla_qubits_green
        )

        self.corner_correction_dict()
        self.ancilla_coords_to_index = {
            coords: index for index, coords in enumerate(self.ancilla_qubits)
        }
        self.ancilla_coords_to_index.update({"vG": len(self.ancilla_qubits)})
        self.ancilla_coords_to_index.update({"vR": len(self.ancilla_qubits) + 1})
        self.ancilla_coords_to_index.update({"vB": len(self.ancilla_qubits) + 2})

        self.ancilla_index_to_coords = {
            index: coords for coords, index in self.ancilla_coords_to_index.items()
        }

    def corner_correction_dict(self):
        if self.distance == 7:
            self.corner_corrections = {
                "(8, 8)b": (9, 9),
                "(2, 0)b": (0, 0),
                "(17, 1)b": (18, 0),
            }
        elif self.distance == 5:
            self.corner_corrections = {
                "(11, 1)b": (12, 0),
                "(5, 5)b": (6, 6),
                "(2, 0)b": (0, 0),
            }
        elif self.distance == 3:
            self.corner_corrections = {
                "(2, 0)b": (0, 0),
                "(5, 1)b": (6, 0),
                "(2, 2)b": (3, 3),
            }

    def red_row(self, x_max, y):
        """
        this is also where data qubits are created
        """
        i = 0
        x_row = y

        while i < (x_max):
            data_or_ancilla = i % 3

            if data_or_ancilla == 0 or data_or_ancilla == 2:
                self.data_qubits.add((x_row, y))
                if y == 0:
                    self.logical_operator.add((x_row, 0))
            else:
                self.ancilla_qubits_red.add((x_row, y))
                self.ancilla_qubits.add((x_row, y))
                if i == 1:
                    self.ancilla_blue_boundary.add((x_row, y))
            i += 1
            x_row += 2

    def blue_row(self, x_max, y):
        """
        create blue ancilla qubits
        """
        i = 0
        x_row = y
        while i < x_max:
            data_or_ancilla = i % 3

            if data_or_ancilla == 0 or data_or_ancilla == 1:
                self.data_qubits.add((x_row, y))

            else:
                self.ancilla_qubits_blue.add((x_row, y))
                self.ancilla_qubits.add((x_row, y))
                if i == x_max - 1:
                    self.ancilla_red_boundary.add((x_row, y))
            i += 1
            x_row += 2

    def green_row(self, x_max, y):
        """
        create green ancilla qubits
        """
        i = 0
        x_row = y
        while i < x_max:

            data_or_ancilla = i % 3

            if data_or_ancilla == 1 or data_or_ancilla == 2:
                self.data_qubits.add((x_row, y))
            else:
                self.ancilla_qubits_green.add((x_row, y))
                self.ancilla_qubits.add((x_row, y))
                if i == 0:
                    self.ancilla_blue_boundary.add((x_row, y))
                if i == (x_max - 2):
                    # elif i == (x_max - 2):
                    self.ancilla_red_boundary.add((x_row, y))
            x_row += 2
            i += 1

    def pairs_to_path(self, v0, v1, set_of_ancillas):
        if type(v1) == str:
            path = self.path_node_to_boundary(v0, set_of_ancillas)
        else:
            path = self.path_two_nodes(v0, v1, set_of_ancillas)
        return path

    def add_boundary_edge(self, v0, matching_graph):
        matching_graph.add_edge(v0, str(v0) + "b")
        if v0 in self.ancilla_blue_boundary:
            matching_graph.nodes[str(v0) + "b"]["colour"] = "b"
        elif v0 in self.ancilla_green_boundary:
            matching_graph.nodes[str(v0) + "b"]["colour"] = "g"
        else:
            matching_graph.nodes[str(v0) + "b"]["colour"] = "r"
        return matching_graph

    def data_qubits_to_ancilla_qubits(self, data_cords):
        ancilla_qubits = set()
        for data_qubits in data_cords:
            potential_ancilla_qubits = [
                (data_qubits[0] + 1, data_qubits[1] + 1),
                (data_qubits[0] + 2, data_qubits[1]),
                (data_qubits[0] + 1, data_qubits[1] - 1),
                (data_qubits[0] - 1, data_qubits[1] - 1),
                (data_qubits[0] - 2, data_qubits[1]),
                (data_qubits[0] - 1, data_qubits[1] + 1),
            ]

            new_ancilla_qubits = set(
                qubit
                for qubit in potential_ancilla_qubits
                if qubit in self.ancilla_qubits
            )

            # xor (if ancilla is found even amount of times it does not fire)
            ancilla_qubits ^= new_ancilla_qubits

        return ancilla_qubits

