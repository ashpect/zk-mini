import galois
from py_ecc.bn128 import multiply, G1, add
from r1cs.qap import calculate_t_and_h
from utils.matrix import matrix_GF
from utils.constants import p, GF

def prover_grothv0(r1cs, alpha_g1, beta_g2, tau_g1, tau_g2, ht_srs, phi_arr, L_polys, R_polys, O_polys, w):
    """
    Input : 
    - r1cs - R1CS system
    - alpha_g1, beta_g2 - elements in G1, G2
    - tau_g1, tau_g2 - tuple of elements in G1, G2
    - ht_srs - tuple of elements in G1
    - phi_arr - tuple of elements in G1
    - L_polys, R_polys, O_polys - tuple of Galois Polynomials
    - w - witness
    Output :
    - A, B, C - tuple of elements in G1, G2, G1
    """
    # Get in the galois field
    W_galois = matrix_GF(w.matrix)

    # QAP polys were calulated by trusted setup, so no need to recompute
    # Calculate Lw_qap, Rw_qap, Ow_qap
    Lw_qap = galois.Poly([0], field=GF)
    for i in range(w.size):
        Lw_qap = L_polys[i]*W_galois[i] + Lw_qap

    Rw_qap = galois.Poly([0], field=GF)
    for i in range(w.size):
        Rw_qap = R_polys[i]*W_galois[i] + Rw_qap
    
    Ow_qap = galois.Poly([0], field=GF)
    for i in range(w.size):
        Ow_qap = O_polys[i]*W_galois[i] + Ow_qap

    # Get h and t
    _, h_qap = calculate_t_and_h(Lw_qap, Rw_qap, Ow_qap, r1cs.L.rows_size)

    A = add(alpha_g1, inner_product(Lw_qap, tau_g1))
    B = add(beta_g2, inner_product(Rw_qap, tau_g2))
    C = inner_product(h_qap, ht_srs) 

    # phi_arr is in G1, W_galois is in GF(p), int'ing for now for faster calculation
    for i in range(len(phi_arr)):
        if multiply(phi_arr[i],  int(W_galois[i])) is not None: # TODO : HANDLE THIS CASE - multiply (G1, 0) returns [None]
            C = add(C, multiply(phi_arr[i],  int(W_galois[i])))

    print("Prover A,B,C generation complete")

    return A, B, C

def inner_product(poly, value_arr):
    """
    Calculate the inner product of a polynomial and an array of values
    """
    # The arr given will be of type [G1, tau*G1 .... tau^n-1*G1]
    coefficients = poly.coeffs # in a degree ascending order
    # tau_g1 has elements in G1 with field order as field_modulus and coeff[i] are in modulus p
    # For now doing int for faster calculation and converting to field_modulus since p is significantly lower
    result = None
    for i in range(len(value_arr)):
        if result is None:
            result = multiply(value_arr[i], int(coefficients[i]))
        else:
            result = add(result, multiply(value_arr[i], int(coefficients[i])))
    # print("type of result", type(result), type(result[0]))
    return result