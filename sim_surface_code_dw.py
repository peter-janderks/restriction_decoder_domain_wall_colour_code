from qec.hgp import hgp
from ldpc.codes import rep_code
import scipy
import numpy as np
from pymatching import Matching


class SimSurfaceCode(object):
    def __init__(
        self,
        n_1,
        n_2,
        p1,
        code="Surface_Code",
        decoding_method="MWPM",
        bias =0,
        bias_type = 'Z'
        
    ):
        self.p1 = p1
        self.n1 = n_1
        self.n2 = n_2
        self.distance = min(n_1, n_2)

        self.initialize_surface_code()
        if code == 'XXXZ':
            self.initialize_error_model_XXXZ(p1, bias, bias_type)
        elif code == '3dw':
            self.initialize_error_model_3dw(p1, bias, bias_type)
        elif code == 'XZZX':
            self.initialize_error_model_XZZX(p1, bias, bias_type)

        else:
            self.initialize_error_model(p1, code, bias)

        if decoding_method == "MWPM":
            self.initialize_decoders_mwpm(manhatten=False)


    def initialize_surface_code(self):
        sc = hgp(rep_code(self.n1), rep_code(self.n2))
        # switched to agree with qecsim
        self.Hx = scipy.sparse.csr_matrix(sc.hz)
        self.Hz = scipy.sparse.csr_matrix(sc.hx)
        self.num_stabilisers, self.num_qubits = self.Hz.shape
        self.lx = sc.lz  # switched to agree with qecsim
        self.lz = sc.lx
        self.N = sc.N

    def surface_code_error_model(self,bias, bias_type):
        if bias_type == "X":
            pz_sc =  1 / (2 * (bias + 1)) 
            py_sc =  1 / (2 * (bias + 1))
            px_sc = bias / (bias + 1) 
        elif bias_type == "Y":
            px_sc =  1 / (2 * (bias + 1)) 
            pz_sc =  1 / (2 * (bias + 1))
            py_sc = bias / (bias + 1) 
        elif bias_type == "Z":
            px_sc =  1 / (2 * (bias + 1)) 
            py_sc =  1 / (2 * (bias + 1))
            pz_sc = bias / (bias + 1) 
        return(px_sc, py_sc, pz_sc)
    
    def initialize_error_model_XZZX(self, per, bias=0, bias_type="Z"):
        self.bias = bias
        self.bias_type=bias_type
        px, py, pz =  np.zeros(self.N), np.zeros(self.N), np.zeros(self.N)
        q_index = 0
        if bias_type in ['X','Y','Z']:
            px_sc, py_sc, pz_sc = self.surface_code_error_model(bias,bias_type)
        elif bias_type == 'Z_inf':
            px_sc, py_sc,pz_sc = 0,0,1

        for i in range(2*self.n1-1):
            if i < self.n1:
                for _ in range(self.n2):
                    px[q_index], py[q_index], pz[q_index] = px_sc, py_sc, pz_sc
                    q_index += 1
            else:
                for _ in range(self.n2-1):
                    pz[q_index], py[q_index], px[q_index] = px_sc, py_sc, pz_sc
                    q_index += 1
        self.px = px * per
        self.pz = pz * per
        self.py = py * per
        self.pxy = self.px+self.py
        self.pzy = self.pz+self.py


    def initialize_error_model_3dw(self, per, bias=0,bias_type="Z"):
        self.bias = bias
        self.bias_type=bias_type
        px, py, pz =  np.zeros(self.N), np.zeros(self.N), np.zeros(self.N)
        q_index = 0
        if bias_type in ['X','Y','Z']:
            px_sc, py_sc, pz_sc = self.surface_code_error_model(bias,bias_type)
        elif bias_type == 'Z_inf':
            px_sc, py_sc,pz_sc = 0,0,1

        for i in range(self.n1):
            if i % 3 == 2:
                for _ in range(self.n2):
                    px[q_index], py[q_index], pz[q_index] = px_sc, py_sc, pz_sc
                    q_index += 1
            else:
                for _ in range(self.n2):
                    pz[q_index], py[q_index], px[q_index] = px_sc, py_sc, pz_sc
                    q_index += 1

        for i in range(self.n1-1):
            if i % 3 == 0:
                for _ in range(self.n2-1):
                    pz[q_index], py[q_index], px[q_index] = px_sc, py_sc, pz_sc
                    q_index += 1
            else:
                for _ in range(self.n2-1):
                    px[q_index], py[q_index], pz[q_index] = px_sc, py_sc, pz_sc
                    q_index += 1


        self.px = px * per
        self.pz = pz * per
        self.py = py * per
        self.pxy = self.px+self.py
        self.pzy = self.pz+self.py

    def initialize_error_model_XXXZ(self, per, bias=0,bias_type="Z"):
        self.bias = bias
        self.bias_type=bias_type
        px, py, pz =  np.zeros(self.N), np.zeros(self.N), np.zeros(self.N)
        q_index = 0
        if bias_type in ['X','Y','Z']:
            px_sc, py_sc, pz_sc = self.surface_code_error_model(bias,bias_type)
        elif bias_type == 'Z_inf':
            px_sc, py_sc,pz_sc = 0,0,1

        for i in range(self.n1):
            if i % 2 == 0:
                for _ in range(self.n2):
                    px[q_index], py[q_index], pz[q_index] = px_sc, py_sc, pz_sc
                    q_index += 1
            else:
                for _ in range(self.n2):
                    pz[q_index], py[q_index], px[q_index] = px_sc, py_sc, pz_sc
                    q_index += 1

        for i in range(self.n1-1):
            if i % 2 == 0:
                for _ in range(self.n2-1):
                    pz[q_index], py[q_index], px[q_index] = px_sc, py_sc, pz_sc
                    q_index += 1
            else:
                for _ in range(self.n2-1):
                    px[q_index], py[q_index], pz[q_index] = px_sc, py_sc, pz_sc
                    q_index += 1


        self.px = px * per
        self.pz = pz * per
        self.py = py * per
        self.pxy = self.px+self.py
        self.pzy = self.pz+self.py


    def initialize_error_model( 
        self, per, code, fixed_seed, non_uniform, std_t=0, std_p=0.5,
        bias=0
    ):
        self.per = per
        self.bias = bias
        qubit_cords_even = np.array(
            [
                [hor_cor, ver_cor]
                for hor_cor in np.arange(0, 2 * self.n1 - 1, 2)
                for ver_cor in np.arange(0, 2 * self.n2 - 1, 2)
            ]
        )
        qubit_cords_uneven = np.array(
            [
                [hor_cor, ver_cor]
                for hor_cor in np.arange(1, 2 * self.n1 - 1, 2)
                for ver_cor in np.arange(1, 2 * self.n2 - 1, 2)
            ]
        )

        

        qubit_cords = np.concatenate((qubit_cords_even, qubit_cords_uneven))
        self.mean = 0.5
        self.std = std_p
        if non_uniform is True:
            self.std_t = std_t

        qecsim_code = LocalCodePlanar(
            self.n1,
            self.n2,
            self.mean,
            self.std,
            seed_h=np.random.randint(1, 10000),
            seed_m=np.random.randint(1, 10000),
            seed_l=np.random.randint(1, 10000),
            seed_n=np.random.randint(1, 10000),
            bias=self.bias,
        )  
        px, py, pz = np.zeros(self.N), np.zeros(self.N), np.zeros(self.N)

        if code == "XZZX":
            q_index = 0
            for q_cords in qubit_cords_even:

                (
                    pz[q_index],
                    py[q_index],
                    px[q_index],
                ) = qecsim_code.qubit_error_XZZX_layout(self.bias)[:, q_cords[0], q_cords[1]]
                q_index += 1
            for q_cords in qubit_cords_uneven:

                (
                    px[q_index],
                    py[q_index],
                    pz[q_index],
                ) = qecsim_code.qubit_error_XZZX_layout(self.bias)[:, q_cords[0], q_cords[1]]

                q_index += 1
        
        elif code == "Surface_Code":
            q_index=0
            for q_index, q_cords in enumerate(qubit_cords):
                    (
                        px[q_index],
                        py[q_index],
                        pz[q_index],
                    ) = qecsim_code.qubit_error_probabilities_biased(self.bias)[
                        :, q_cords[0], q_cords[1]
                    ]

        self.px = px * per
        self.pz = pz * per
        self.py = py * per
        self.pxy = self.px+self.py
        self.pzy = self.pz+self.py

    def initialize_decoders_mwpm(self, manhatten=False):
        weights_x = [np.log((1 - p) / p) if p!=0 else 100 for p in self.pxy ]
        weights_z = [np.log((1 - p) / p) if p!=0 else 100 for p in self.pzy ]

        self.x_decoder = Matching(self.Hz, spacelike_weights=weights_x)

        self.z_decoder = Matching(self.Hx, spacelike_weights=weights_z)

    def run_single_round(self):
        x_error = np.array(
            [1 if np.random.random() < self.px[i] else 0 for i in range(self.N)]
        )
        y_error = np.array(
            [1 if np.random.random() < self.py[i] else 0 for i in range(self.N)]
        )
        z_error = np.array(
            [1 if np.random.random() < self.pz[i] else 0 for i in range(self.N)]
        )

        x_syndrome = self.Hz @ (x_error ^ y_error) % 2
        z_syndrome = self.Hx @ (z_error ^ y_error) % 2

        x_correction = self.x_decoder.decode(x_syndrome, num_neighbours=None)
        z_correction = self.z_decoder.decode(z_syndrome, num_neighbours=None)
        x_residual = (x_error ^ y_error + x_correction) % 2
        z_residual = (z_error ^ y_error + z_correction) % 2
        logical_x_error = (self.lz @ x_residual) % 2
        logical_z_error = (self.lx @ z_residual) % 2

        return logical_x_error[0],logical_z_error[0]

    def run(self, max_runs, max_logical_errors):
        n_logical_errors_x = 0
        n_logical_errors_z = 0
        n_logical_errors_total = 0
        n_runs = 0

        for n_runs in range(max_runs):
            new_lx, new_lz = self.run_single_round()
            n_logical_errors_x += new_lx
            n_logical_errors_z += new_lz            

            if new_lx == 1 or new_lz == 1:
                n_logical_errors_total += 1

            if n_logical_errors_total == max_logical_errors:
                break
        return (n_logical_errors_total, n_logical_errors_x, n_logical_errors_z, n_runs+1)



if __name__ == "__main__":

    decoder = "MWPM"
    #"""
    sim_setup = SimSurfaceCode(
        9,
        3*9,
        0.2,
        decoding_method='MWPM',
        code="3dw",
        bias=100,
        bias_type='Z'
    )
    ler_x, ler_z, n_runs = sim_setup.calculate_logical_error_rate(10)
    print(ler_x, n_runs, 'x')
    print(ler_z, n_runs, 'z')
    
    """
    sim_setup = SimSurfaceCode(
        5,
        11*5,
        0.35,
        decoding_method='MWPM',
        code="XZZX",
        bias=30,
        bias_type='Z'
    )
    ler_x, ler_z, n_runs = sim_setup.calculate_logical_error_rate(10)
    print(ler_x, n_runs, 'x')
    print(ler_z, n_runs, 'z')
    """
