from py_ecc.bn128 import G1, G2, multiply, add, curve_order
import numpy as np

def prover(r1cs):
    """
    Prover does the following :
    1. Convert witness vector to Elliptic Curve Points E1 and E2
    """ 
    # Step 1 : Converting witness vector to Elliptic Curve Point
    Witness_ec1 = r1cs.W.get_witness_G1()
    Witness_ec2 = r1cs.W.get_witness_G2()


    # Step 2 : Finding LW1, Rw2, OW1 [Lw, Rw, Ow for simplicity]
    W = r1cs.W.matrix
    L = r1cs.L.matrix

    # print(W.shape)
    # print(L.shape)
    # print(Witness_ec1.shape)
    # print(Witness_ec2.shape)

    Lw = np.full(L.shape[0], G1)
    for j in range(L.shape[0]):
        point = multiply(Witness_ec1[0], L[j][0]) # NOTE - curve points comes first, learnt it the hard way lol.
        for i in range(1, W.size): 
            point = add(point, multiply(Witness_ec1[i], L[j][i]))
        Lw[j] = point

    R = r1cs.R.matrix
    Rw = np.full(R.shape[0], G2)
    for j in range(len(R)):
        point = multiply(Witness_ec2[0], R[j][0])
        for i in range(1, W.size):
            point = add(point, multiply(Witness_ec2[i], R[j][i]))
        Rw[j] = point

    OUT = r1cs.OUT.matrix
    Ow_g1 = np.full(OUT.shape[0], G1)
    for j in range(len(OUT)):
        if OUT[j][0] >= 0:
            point = multiply(Witness_ec1[0], OUT[j][0])
        else:
            point = multiply(Witness_ec1[0], (curve_order-1)) # To handle -1 case, generalize it later.

        for i in range(1, W.size):
            if OUT[j][i] >= 0:
                point = add(point, multiply(Witness_ec1[i], OUT[j][i]))
            else:
                point = add(point, multiply(Witness_ec1[i], (curve_order-1))) # To handle -1 case, generalize it later.

        Ow_g1[j] = point

    return Ow_g1, Rw, Lw