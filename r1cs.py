import numpy as np
from py_ecc.bn128 import G1, G2, multiply

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
