import json
from pathlib import Path
from typing import List
from calculate_LER import Calculate_LER
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
import random

from sim_surface_code_dw import SimSurfaceCode


class Calculate_Threshold:
    def __init__(
        self, per_list, distance_array, n_runs, n_logical_errors, bias, cpus, code
    ):

        self.per_list = list(per_list)
        self.distance_array = list(distance_array)
        self.n_runs = n_runs
        self.bias = bias
        self.n_logical_errors = n_logical_errors
        self.code = code
        self.ler_distances = list()
        for d in distance_array:
            print(d, "d")
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
                    else:
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
                "./data_18_11/6.6.6_DW(23)_-pi4_bias_"
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

    n_runs = 100000
    n_logical_errors = 2000
    cpus = 20
    bias_list = [0.5, 1, 3, 10, 30, 100, 300, 1000]
    distance_ratio = [1, 3, 5, 7, 11, 23, 53, 157]

    code = "XZZX"
    per_lists = [
        np.linspace(0.1, 0.2, 11),
        np.linspace(0.1, 0.2, 11),
        np.linspace(0.15, 0.25, 11),
        np.linspace(0.17, 0.32, 16),
        np.linspace(0.25, 0.4, 16),
        np.linspace(0.3, 0.45, 16),
        np.linspace(0.37, 0.47, 11),
        np.linspace(0.38, 0.48, 11),
    ]
    assert len(bias_list) == len(per_lists)
    assert len(bias_list) == len(distance_ratio)

    distance_list = []
    for index,ratio in enumerate(distance_ratio):
        if index < 5:
            distance_list.append(
                [[3, 3 * ratio], [5, 5 * ratio], [7, 7 * ratio], [9, 9 * ratio]]
            )
        else:
            distance_list.append(
                [[3, 3 * ratio], [5, 5 * ratio], [7, 7 * ratio]]
            )

    for index, bias in enumerate(bias_list):
        print(bias, "bias")
        Calculate_Threshold(
            per_lists[index],
            distance_list[index],
            n_runs,
            n_logical_errors,
            bias,
            cpus,
            code,
        )

    n_runs = 1000000
    n_logical_errors = 5000
    cpus = 20
    distance_array = [5, 7,9, 13, 15, 17, 21,25]
    per_lists = []
    bias_list = []

    bias_list.append(0.5)
    per_lists.append(np.linspace(0.08, 0.18, 21))

    bias_list.append(1)
    per_lists.append(np.linspace(0.1, 0.2, 21))

    bias_list.append(3)
    per_lists.append(np.linspace(0.12, 0.22, 21))

    bias_list.append(10)
    per_lists.append(np.linspace(0.16, 0.26, 21))

    bias_list.append(30)
    per_lists.append(np.linspace(0.2, 0.3, 21))

    bias_list.append(100)
    per_lists.append(np.linspace(0.25, 0.35, 21))

    bias_list.append(300)
    per_lists.append(np.linspace(0.28, 0.38, 21))

    bias_list.append(1000)
    per_lists.append(np.linspace(0.3, 0.4, 21))

    bias_list.append(3000)
    per_lists.append(np.linspace(0.34, 0.44, 21))

    bias_list.append(10000)
    per_lists.append(np.linspace(0.38, 0.5, 21))

    bias_list.append(30000)
    per_lists.append(np.linspace(0.38, 0.5, 21))

    bias_list.append(100000)
    per_lists.append(np.linspace(0.38, 0.5, 21))

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
