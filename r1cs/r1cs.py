import numpy as np
from py_ecc.bn128 import G1, G2, multiply

class R1cs_matrix:
    def __init__(self, matrix):
        self.matrix = matrix
        self.rows_size = matrix.shape[0]  
        self.cols_size = matrix.shape[1]
        
class Witness:
    def __init__(self, witness):
        # W: witness vector (1D array)
        self.matrix = witness # naming it matrix for consistency
        self.size = witness.shape[0]

    def get_witness_G1(self):
        # Converts from shape (n,) to (n, 2)
        return np.array([multiply(G1, self.matrix[i]) for i in range(self.size)])
    
    def get_witness_G2(self):
        # Converts from shape (n,) to (n, 2)
        return np.array([multiply(G2, self.matrix[i]) for i in range(self.size)])

class R1CS:
    def __init__(self, L: R1cs_matrix, R: R1cs_matrix, OUT: R1cs_matrix):
        # L, R, OUT: constraint matrices (2D arrays with same dimensions)
        self.L = L
        self.R = R
        self.OUT = OUT

def sanity_check(L, R, OUT, w):
        assert(np.all(
            np.matmul(OUT, w) ==
            np.multiply(
                np.matmul(L, w),
                np.matmul(R, w)
            )
        ))