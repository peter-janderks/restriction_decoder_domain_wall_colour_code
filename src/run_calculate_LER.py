import json
from pathlib import Path
from typing import List
from calculate_LER import Calculate_LER
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
import random


class Calculate_Threshold:
    def __init__(
        self, per_list, distance_array, n_runs, n_logical_errors, bias, cpus, code
    ):
        """Runs Monte Carlo simulations of a QEC code using multiproccessing.

        Data is saved to a json file. See the save function for details.

        Args:
            per_list: List of physical error rates (floats) to simulate.
            distance_array: List of distances (integers) to simulate.
            n_runs: Max number or runs to simulate.
            n_logical_errors: Max number of logical error to simulate.
            bias: The bias to simulate, represented using an integer or float.
            cpus: The number of cpus to use.
            code: The name of the code to simulate, right now either ColourCode or SurfaceCode.
        """
        self.per_list = list(per_list)
        self.distance_array = list(distance_array)
        self.n_runs = n_runs
        self.bias = bias
        self.n_logical_errors = n_logical_errors
        self.code = code
        self.ler_distances = list()
        for d in distance_array:
            ler_a = []
            ler_x_a = []
            ler_z_a = []
            ler_eb_array = []
            for per in per_list:
                mp.freeze_support()
                pool = mp.Pool()
                results = []
                for i in range(cpus):
                    if code == "ColourCode":
                        result = pool.apply_async(
                            self.task,
                            args=(
                                i,
                                d,
                                per,
                                self.bias,
                                self.n_runs // cpus,
                                n_logical_errors // cpus,
                            ),
                        )
                    elif code == "SurfaceCode":
                        result = pool.apply_async(
                            self.surface_code_task,
                            args=(
                                i,
                                d,
                                per,
                                self.bias,
                                self.n_runs // cpus,
                                n_logical_errors // cpus,
                            ),
                        )
                    else:
                        raise ValueError(
                            f"code must be SurfaceCode or ColourCode, here {code} was passed in."
                        )
                    results.append(result)
                pool.close()
                pool.join()
                total_errors = 0
                x_errors = 0
                z_errors = 0
                total_runs = 0
                for result in results:
                    new_result = result.get()
                    total_errors += new_result[0]
                    x_errors += new_result[1]
                    z_errors += new_result[2]
                    total_runs += new_result[3]

                ler = total_errors / total_runs
                ler_x = x_errors / total_runs
                ler_z = z_errors / total_runs
                ler_a.append(ler)
                ler_x_a.append(ler_x)
                ler_z_a.append(ler_z)
                ler_eb_array.append(np.sqrt((1 - ler) * ler / total_runs))

            if type(d) == list:
                ratio = d[1] / d[0]
                d = d[0]

            else:
                ratio = None
            self.ler_distances.append(ler_a)
            self.data = dict()

            self.data[str(d)] = dict()
            self.data[str(d)]["ler"] = list(ler_a)
            self.data[str(d)]["ler_x"] = list(ler_x_a)
            self.data[str(d)]["ler_z"] = list(ler_z_a)
            self.data[str(d)]["per"] = list(per_list)
            self.data[str(d)]["ler_eb"] = list(ler_eb_array)
            self.data[str(d)]["ratio"] = ratio
            self.save_data(d)

    def task(self, pid, d, per, bias, n_runs, n_logicals):
        np.random.seed(random.randint(10000, 20000))
        model = Calculate_LER(d, per, bias)
        l_errors, l_errors_x, l_errors_z, n_runs = model.run(n_runs, n_logicals)
        return l_errors, l_errors_x, l_errors_z, n_runs

    def surface_code_task(self, pid, d, per, bias, n_runs, n_logicals):
        np.random.seed(random.randint(10000, 20000))
        model = SimSurfaceCode(d[0], d[1], per, code=self.code, bias=bias)
        l_errors, l_errors_x, l_errors_z, n_runs = model.run(n_runs, n_logicals)
        return l_errors, l_errors_x, l_errors_z, n_runs

    def save_data(self, d):
        if self.code == "ColourCode":
            file_name = Path(
                "./example_data/6.6.6_DW(23)_-pi4_bias_"
                + str(self.bias)
                + "/_n_runs_"
                + str(self.n_runs)
                + "_n_logical"
                + str(self.n_logical_errors)
                + "_distance"
                + str(d)
                + "_"
                + str(np.random.randint(1, 100000))
                + ".json"
            )
        else:
            file_name = Path(
                f"./data_18_11/{self.code}_bias_"
                + str(self.bias)
                + "/_n_runs_"
                + str(self.n_runs)
                + "_n_logical"
                + str(self.n_logical_errors)
                + "_distance"
                + str(d)
                + "_"
                + str(np.random.randint(1, 100000))
                + ".json"
            )

        file_name.parent.mkdir(exist_ok=True, parents=True)
        f = open(file_name, "w+")

        print(json.dumps(self.__dict__, sort_keys=True, indent=4), file=f)
        f.close()


if __name__ == "__main__":
    # update the following parameters to change the simulation
    n_runs = 100 # maximum number of shots to sample
    n_logical_errors = 50 # maximum number of logical errors to sample
    cpus = 10 # number of cpus to use
    distance_array = [9,11]
    per_lists = []
    bias_list = []
    bias_list.append(30)
    per_lists.append(np.linspace(0.25, 0.4, 21))

    for index, bias in enumerate(bias_list):

        Calculate_Threshold(
            per_lists[index],
            distance_array,
            n_runs,
            n_logical_errors,
            bias,
            cpus,
            "ColourCode",
        )
