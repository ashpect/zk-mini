import numpy as np
from r1cs import R1CS, R1cs_matrix, Witness

def mock_example_1():
    """
    Get an R1CS example
    """
 
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