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
            ler_x_a = []
            ler_z_a = []
            ler_eb_array = []
            for per in per_list:
                print(per, "per")
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

            self.ler_distances.append(ler_a)
            self.data = dict()
            self.data[str(d)] = dict()
            self.data[str(d)]["ler"] = list(ler_a)
            self.data[str(d)]["ler_x"] = list(ler_x_a)
            self.data[str(d)]["ler_z"] = list(ler_z_a)
            self.data[str(d)]["per"] = list(per_list)
            self.data[str(d)]["ler_eb"] = list(ler_eb_array)
            self.save_data(d)

    def task(self, pid, d, per, bias, n_runs, n_logicals):
        np.random.seed(random.randint(10000, 20000))
        model = Calculate_LER(d, per, bias)
        l_errors, l_errors_x, l_errors_z, n_runs = model.run(n_runs, n_logicals)
        return l_errors, l_errors_x, l_errors_z, n_runs

    def save_data(self, d):

        file_name = Path(
            "./data_18_11/6.6.6_DW(23)_-pi4_bias_"
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

    n_runs = 1000000
    n_logical_errors = 1000
    cpus = 20

    # per_list = np.linspace(0.3, 0.4, 11)
    #   bias = 300
    distance_array = [5, 9, 13, 17, 21]
    per_lists = []
    bias_list = []

    bias_list.append(0.5)
    per_lists.append(np.linspace(0.08,0.18,21))

    bias_list.append(1)
    per_lists.append(np.linspace(0.1,0.2,21))

    bias_list.append(3)
    per_lists.append(np.linspace(0.12,0.22,21))

    bias_list.append(10)
    per_lists.append(np.linspace(0.16,0.26,21))

    bias_list.append(30)
    per_lists.append(np.linspace(0.2,0.3,21))


    bias_list.append(100)
    per_lists.append(np.linspace(0.25,0.35,21))

    bias_list.append(300)
    per_lists.append(np.linspace(0.28,0.38,21))
    
    bias_list.append(1000)
    per_lists.append(np.linspace(0.3,0.4,21))

    bias_list.append(3000)
    per_lists.append(np.linspace(0.34,0.44,21))

    bias_list.append(10000)
    per_lists.append(np.linspace(0.38,0.5,21))

    bias_list.append(30000)
    per_lists.append(np.linspace(0.38,0.5,21))

    for index, bias in enumerate(bias_list):    
        Calculate_Threshold(per_lists[index], distance_array, n_runs, n_logical_errors, bias, cpus)
