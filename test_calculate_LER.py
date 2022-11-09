import json
from pathlib import Path
from calculate_LER import Calculate_LER
import matplotlib.pyplot as plt
import numpy as np
import multiprocessing as mp
import random


class Calculate_Threshold:
    def __init__(self, per_list, distance_array, n_runs, n_logical_errors, bias, cpus):

        self.per_list = list(per_list)
        self.distance_array = list(distance_array)
        self.n_runs = n_runs
        self.bias = bias
        self.n_logical_errors = n_logical_errors
        self.ler_distances = list()
        print(self.bias, "bias")
        for d in distance_array:
            print(d, "distance")
            ler_a = []
            ler_eb_array = []
            for per in per_list:
                mp.freeze_support()
                pool = mp.Pool()
                results = []
                for i in range(cpus):
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
                    results.append(result)
                pool.close()
                pool.join()
                total_errors = 0
                total_runs = 0
                for result in results:
                    new_result = result.get()
                    total_errors += new_result[0]
                    total_runs += new_result[1]
                ler = total_errors / total_runs
                ler_a.append(ler)
                ler_eb_array.append(np.sqrt((1 - ler) * ler / total_runs))

            self.ler_distances.append(ler_a)
            self.data = dict()
            self.data[str(d)] = dict()
            self.data[str(d)]["ler"] = list(ler_a)
            self.data[str(d)]["per"] = list(per_list)
            self.data[str(d)]["ler_eb"] = list(ler_eb_array)
            self.save_data(d)

    def task(self, pid, d, per, bias, n_runs, n_logicals):
        np.random.seed(random.randint(10000, 20000))
        model = Calculate_LER(d, per, bias)
        l_errors, n_runs = model.run(n_runs, n_logicals)
        return l_errors, n_runs

    def save_data(self, d):

        file_name = Path(
            "./new_data/6.6.6_DW(23)_-pi4_bias_"
            + str(self.bias)
            + str()
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
    """
    per_list = np.linspace(0.2, 0.4, 20)
    distance_array = [25]
    n_runs = 300000
    n_logical_errors = 1000
    bias = 10
    cpus = 10
    Calculate_Threshold(per_list, distance_array, n_runs, n_logical_errors, bias, cpus)

    per_list = np.linspace(0.25, 0.45, 20)
    distance_array = [7, 9, 11, 13, 15, 17, 19]
    n_runs = 250000
    bias = 30
    cpus = 10
    Calculate_Threshold(per_list, distance_array, n_runs, bias, cpus)
    """
    per_list = np.linspace(0.42, 0.5, 10)
    distance_array = [25]
    n_runs = 300000
    n_logical_errors=1000
    bias = 100
    cpus = 8
    Calculate_Threshold(per_list, distance_array, n_runs,n_logical_errors, bias, cpus)

    per_list = np.linspace(0.45, 0.5, 10)
    distance_array = [13,19,25]
    n_runs = 300000
    n_logical_errors=1000
    bias = 300
    cpus = 8
    Calculate_Threshold(per_list, distance_array, n_runs, n_logical_errors,bias, cpus)

    per_list = np.linspace(0.45, 0.5, 10)
    distance_array = [13,19,25]
    n_runs = 300000
    n_logical_errors=1000
    bias = 1000
    cpus = 8
    Calculate_Threshold(per_list, distance_array, n_runs, n_logical_errors, bias, cpus)

