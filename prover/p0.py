from py_ecc.bn128 import G1, G2, multiply, add, curve_order
import numpy as np

def prover(W, L, R, OUT):
    # Step 1 : Converting witness vector to Elliptic Curve Point
    Witness_ec1 = W.get_witness_G1()
    Witness_ec2 = W.get_witness_G2()

    # Step 2 : Finding E1(Lw), E2(Rw), E1(Ow)
    Lw = np.array([G1,G1])
    for j in range(2):
        point = multiply(Witness_ec1[0], L[j][0])
        for i in range(1, W.size):
            point = add(point, multiply(Witness_ec1[i], L[j][i]))
        Lw[j] = point

    Rw = np.array([G2,G2])
    for j in range(2):
        point = multiply(Witness_ec2[0], R[j][0])
        for i in range(1, W.size):
            point = add(point, multiply(Witness_ec2[i], R[j][i]))
        Rw[j] = point

    Ow_g1 = np.array([G1,G1])
    for j in range(2):
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