from py_ecc.bn128 import G1, G2, multiply, add, curve_order
import numpy as np
import galois
from utils.matrix import matrix_GF
from utils.lagrange import get_polys
from utils.constants import p, GF

def r1cs_to_qap_polys(r1cs):
    """
    Converts R1CS system to QAP Polynomials
    """
    
    # Convert matrix elements to Galois Field elements
    W_galois = matrix_GF(r1cs.W.matrix)
    L_galois = matrix_GF(r1cs.L.matrix)
    R_galois = matrix_GF(r1cs.R.matrix)
    OUT_galois = matrix_GF(r1cs.OUT.matrix)

    # Sanity Check :
    assert(np.all(np.matmul(OUT_galois, W_galois) == np.matmul(L_galois, W_galois) * np.matmul(R_galois, W_galois)))

    L_polys, R_polys, O_polys = get_polys(L_galois, R_galois, OUT_galois, r1cs.L.rows_size)
    print("R1CS to QAP Polys Done")

    return L_polys, R_polys, O_polys, W_galois

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

    print("Prover Done")
    return Ow_g1, Rw, Lw

def prover_qap(r1cs, tau):
    """
    Prover does the following :
    TODO
    """ 
    L_polys, R_polys, O_polys, W_galois = r1cs_to_qap_polys(r1cs)
    
    # Random poly of deg n-1 in GF(p)
    polynomial = galois.primitive_poly(p, r1cs.L.rows_size-1 , method="random")

    # Final Lw , initiating with random poly and size m(colms) for example
    Lw_qap = galois.Poly([0], field=GF)
    for i in range(r1cs.W.size):
        Lw_qap = L_polys[i]*W_galois[i] + Lw_qap

    Rw_qap = galois.Poly([0], field=GF)
    for i in range(r1cs.W.size):
        Rw_qap = R_polys[i]*W_galois[i] + Rw_qap
    
    Ow_qap = galois.Poly([0], field=GF)
    for i in range(r1cs.W.size):
        Ow_qap = O_polys[i]*W_galois[i] + Ow_qap

    # Calulate h
    t, h = calculate_t_and_h(Lw_qap, Rw_qap, Ow_qap, r1cs.L.rows_size)

    # Get values of all poly at tau
    Lw_qap_tau = Lw_qap(tau)
    Rw_qap_tau = Rw_qap(tau)
    Ow_qap_tau = Ow_qap(tau)
    t_tau = t(tau)
    h_tau = h(tau)

    return Lw_qap_tau, Rw_qap_tau, Ow_qap_tau, t_tau, h_tau

    
def calculate_t_and_h(Lw_qap, Rw_qap, Ow_qap, deg):
    """
    Calculate the h polynomial
    """
    # Sanity Check : Deg of Lw_qap, Rw_qap, Ow_qap should be at most n-1
    assert(Lw_qap.degree <= deg-1 and Rw_qap.degree <= deg-1 and Ow_qap.degree <= deg-1)

    # Verification step will be : 
    # The "underlying" vectors of Lw_qap*Rw_qap and Ow_qap are equal, even though the polynomials that interpolate have diff degree
    # Underlying here kinda the x axis values
    # For 2 vectors, the product polynomial interpolates the Hadamard product

    t = galois.Poly([1], field=GF)
    for i in range(1, deg):
        t = t * galois.Poly([1, p-i], field=GF)

    h = (Lw_qap*Rw_qap - Ow_qap) // t

    # Check cause above doesn't return remainder
    assert(Lw_qap*Rw_qap == Ow_qap + h*t)

    return t,h