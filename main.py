import numpy as np
from py_ecc.bn128 import G1, G2, pairing, add, multiply, eq, curve_order, FQ12
from verifier.v0 import verify
from prover.p0 import prover

class Witness:
    def __init__(self, W):
        self.witness_vector = W
        self.size = len(W)
        
    def get_witness_G1(self):
        return np.array([multiply(G1, self.witness_vector[i]) for i in range(self.size)])
    
    def get_witness_G2(self):
        return np.array([multiply(G2, self.witness_vector[i]) for i in range(self.size)])
    

# Mock Example : 
def mock_example():
    
    # Eqn :  z = x*x*y + 1
    # Constraints:
    # 1) a = x*x
    # 2) -1 + z = a*y

    witness = [1, 19, 3, 2, 9]  # Witness vector [1, z, x, y, a] for x = 3, y = 2, z = 19, a = 9
    W = Witness(witness)

    # OUT = Lw*Rw (element wise - *)
    OUT = np.array([
        [0, 0, 0, 0, 1],
        [-1, 1, 0, 0, 0],
    ])

    L = np.array([
        [0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1],
    ])

    R = np.array([
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0],
    ])

    # Sanity Check :
    assert(np.all(np.matmul(OUT, W.witness_vector) == np.matmul(L, W.witness_vector) * np.matmul(R, W.witness_vector)))

    return W, L, R, OUT

def main():
    # Step 1 : Mock Example
    W, L, R, OUT = mock_example()

    # Step 2 : Prover
    Ow_g1, Rw, Lw = prover(W, L, R, OUT)

    # Step 3 : Verify the proof
    verify(Ow_g1, Rw, Lw)

if __name__ == "__main__":
    main()

