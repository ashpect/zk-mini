import numpy as np
import random
from r1cs.r1cs import R1CS, R1cs_matrix, Witness

def get_mock_example(example_type: int):
    """
    Get a mock example based on the provided number
    """
    if example_type == 1:
        return mock_example_1()
    elif example_type == 2:
        return mock_example_2()
    else:
        raise ValueError(f"Unknown example type: {example_type}")

def mock_example_1():
    """
    Get an R1CS example
    """
 
    # Eqn :  z = x*x*y + 1
    # Constraints:
    # 1) a = x*x
    # 2) -1 + z = a*y

    witness = np.array([1, 19, 3, 2, 9])  # Witness vector [1, z, x, y, a] for x = 3, y = 2, z = 19, a = 9

    # OUT = Lw*Rw (element wise - *)
    out = np.array([
        [0, 0, 0, 0, 1],
        [-1, 1, 0, 0, 0],
    ])

    l = np.array([
        [0, 0, 1, 0, 0],
        [0, 0, 0, 0, 1],
    ])

    r = np.array([
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0],
    ])

    return R1CS(R1cs_matrix(l), R1cs_matrix(r), R1cs_matrix(out)), Witness(witness)

def mock_example_2():
    """
    Get an R1CS example
    """
    
    # Eqn : r = x y z * u + 10

    # v1 = x*y
    # v2 = z*u
    # -10 +r = v1*v2

    L = np.array([
                [0,0,1,0,0,0,0,0],
                [0,0,0,0,1,0,0,0],
                [0,0,0,0,0,0,1,0]
            ])

    R = np.array([[0,0,0,1,0,0,0,0],
              [0,0,0,0,0,1,0,0],
              [0,0,0,0,0,0,0,1]])

    O = np.array([[0,0,0,0,0,0,1,0],
              [0,0,0,0,0,0,0,1],
              [-10,1,0,0,0,0,0,0]])
    
    x = random.randint(1,1000)
    y = random.randint(1,1000)
    z = random.randint(1,1000)
    u = random.randint(1,1000)

    # compute the algebraic circuit
    r = x * y * z * u + 10
    v1 = x*y
    v2 = z*u

    # create the witness vector
    witness = np.array([1, r, x, y, z, u, v1, v2])

    return R1CS(R1cs_matrix(L), R1cs_matrix(R), R1cs_matrix(O)), Witness(witness)

# Doesnt work cause we expect the shape of O,L,R to be 2D 
# def mock_example():
#     """
#     Get an R1CS example
#     """
    
#     # Eqn : 

#     O = np.array([[0, 1, 0, 0]])
#     L = np.array([[0, 0, 1, 0]])
#     R = np.array([[0, 0, 0, 1]])

#     # witness vector
#     witness = np.array([1, 4223, 41, 103])    

#     return R1CS(R1cs_matrix(L), R1cs_matrix(R), R1cs_matrix(O)), Witness(witness)
