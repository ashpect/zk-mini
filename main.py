import numpy as np
from py_ecc.bn128 import G1, G2, pairing, add, multiply, eq, curve_order, FQ12
from verifier.v0 import verify, verify_qap, verifier_qap_tau_random
from prover.p0 import prover, prover_qap

class R1cs_matrix:
    def __init__(self, matrix):
        self.matrix = matrix
        self.rows_size = matrix.shape[0]  
        self.cols_size = matrix.shape[1]
        
class Witness:
    def __init__(self, witness):
        self.matrix = witness # naming it matrix for consistency
        self.size = witness.shape[0]

    def get_witness_G1(self):
        # Converts from shape (n,) to (n, 2)
        return np.array([multiply(G1, self.matrix[i]) for i in range(self.size)])
    
    def get_witness_G2(self):
        # Converts from shape (n,) to (n, 2)
        return np.array([multiply(G2, self.matrix[i]) for i in range(self.size)])

class R1CS:
    def __init__(self, W: Witness, L: R1cs_matrix, R: R1cs_matrix, OUT: R1cs_matrix):
        # W: witness vector (1D array)
        # L, R, OUT: constraint matrices (2D arrays with same dimensions)

        self.W = W
        self.L = L
        self.R = R
        self.OUT = OUT
        self.sanity_check()

    def sanity_check(self):
        assert(np.all(np.matmul(self.OUT.matrix, self.W.matrix) == np.matmul(self.L.matrix, self.W.matrix) * np.matmul(self.R.matrix, self.W.matrix)))

# Mock Example : 
def mock_example():
    
    # Eqn :  z = x*x*y + 1
    # Constraints:
    # 1) a = x*x
    # 2) -1 + z = a*y

    witness = np.array([1, 19, 3, 2, 9])  # Witness vector [1, z, x, y, a] for x = 3, y = 2, z = 19, a = 9
    W = Witness(witness)

    # OUT = Lw*Rw (element wise - *)
    out = np.array([
        [0, 0, 0, 0, 1],
        [-1, 1, 0, 0, 0],
    ])
    OUT = R1cs_matrix(out)

    l = np.array([
        [0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1],
    ])
    L = R1cs_matrix(l)

    r = np.array([
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0],
    ])
    R = R1cs_matrix(r)

    return R1CS(W, L, R, OUT)

def v0_main():
    # Step 1 : Mock Example
    r1cs = mock_example()

    # Step 2 : Prover
    Ow_g1, Rw, Lw = prover(r1cs)

    # Step 3 : Verify the proof
    verify(Ow_g1, Rw, Lw)

def v1_main():
    r1cs = mock_example()
    tau = verifier_qap_tau_random()
    L, R, O, t, h = prover_qap(r1cs, tau)
    verify_qap(L, R, O, t, h)
    pass

if __name__ == "__main__":
    v1_main()

