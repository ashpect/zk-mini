import galois
from py_ecc.bn128 import multiply, G1, add, G2, pairing, neg, FQ12
from r1cs.qap import calculate_t_and_h, r1cs_to_qap
from utils.matrix import matrix_GF
from utils.constants import p, GF
from functools import reduce

def qap_trusted(r1cs, w):

    L_polys, R_polys, O_polys = r1cs_to_qap(r1cs)
    print("L_polys", L_polys)
    print("R_polys", R_polys)
    print("O_polys", O_polys)
    tau = 5 # for testing

    deg = r1cs.L.rows_size
    print(deg)
    tau_g1 = [multiply(G1, tau**i) for i in range(deg-1,-1,-1)]
    tau_g2 = [multiply(G2, tau**i) for i in range(deg-1,-1,-1)] 
    test =  [tau**i for i in range(deg-1,-1,-1)] 
    print("test", test)

    # t = (x-1)(x-2)...(x-n)
    t = galois.Poly([1], field=GF)
    for i in range(1, deg+1):
        t = t * galois.Poly([1, p-i], field=GF)

    t_at_tau = t(tau)
    print("t_at_tau", t_at_tau)
    ht_srs = [multiply(G1, int(t_at_tau * tau**i)) for i in range(deg-2,-1,-1)]
    test2 = [int(t_at_tau * tau**i) for i in range(deg-2,-1,-1)]
    print("test2", test2)

    W_galois = matrix_GF(w.matrix)
    
    Lw_qap = galois.Poly([0], field=GF)
    for i in range(w.size):
        Lw_qap = L_polys[i]*W_galois[i] + Lw_qap

    Rw_qap = galois.Poly([0], field=GF)
    for i in range(w.size):
        Rw_qap = R_polys[i]*W_galois[i] + Rw_qap
    
    Ow_qap = galois.Poly([0], field=GF)
    for i in range(w.size):
        Ow_qap = O_polys[i]*W_galois[i] + Ow_qap

    t, h_qap = calculate_t_and_h(Lw_qap, Rw_qap, Ow_qap, r1cs.L.rows_size)

    print("Lw_qap", Lw_qap)
    print("Rw_qap", Rw_qap)
    print("Ow_qap", Ow_qap)
    print("h_qap", h_qap)
    print("t", t)

    # Sanity Check 
    assert(Lw_qap*Rw_qap == Ow_qap + h_qap*t)

    Lw_qap_coeffs = [int(x) for x in Lw_qap.coeffs]
    Rw_qap_coeffs = [int(x) for x in Rw_qap.coeffs]
    Ow_qap_coeffs = [int(x) for x in Ow_qap.coeffs]
    h_qap_coeffs = [int(x) for x in h_qap.coeffs]

    print("Lw_qap_coeffs", Lw_qap_coeffs)
    print("Rw_qap_coeffs", Rw_qap_coeffs)
    print("Ow_qap_coeffs", Ow_qap_coeffs)
    print("h_qap_coeffs", h_qap_coeffs)

    A = inner_product(tau_g1, Lw_qap_coeffs)
    B = inner_product(tau_g2, Rw_qap_coeffs) 
    C_part1 = inner_product(tau_g1, Ow_qap_coeffs) 
    C_part2 = inner_product(ht_srs, h_qap_coeffs) 

    # Test for inner product
    mockpoly = [multiply(G1, 1) for i in range(deg,-1,-1)]
    assert inner_product(mockpoly, Lw_qap_coeffs) == multiply(G1, sum(Lw_qap_coeffs))
    # TODO : Do i need to handle case when [G1,G1,G1][0,6,14] = 20G1 but the GF is 17, so do i need to reduce ?

    C = add(C_part1, C_part2)

    lhs = pairing(B, neg(A))
    rhs = pairing(G2, C)

    # print("lhs", lhs)
    # print("rhs", rhs)
    assert (lhs * rhs) == FQ12.one(), f" lhs should match rhs"
    print("Done")


def prover_grothv0(r1cs, w):
    """
    Input : 
    - r1cs - R1CS system
    - alpha_g1, beta_g2 - elements in G1, G2
    - tau_g1, tau_g2 - tuple of elements in G1, G2
    - ht_srs - tuple of elements in G1
    - phi_srs - tuple of elements in G1
    - L_polys, R_polys, O_polys - tuple of Galois Polynomials
    - w - witness
    Output :
    - A, B, C - tuple of elements in G1, G2, G1
    """

    L_polys, R_polys, O_polys = r1cs_to_qap(r1cs)
    print("L_polys", L_polys)
    print("R_polys", R_polys)
    print("O_polys", O_polys)
    tau = 5 # for testing

    deg = r1cs.L.rows_size
    print(deg)
    tau_g1 = [multiply(G1, tau**i) for i in range(deg-1,-1,-1)]
    tau_g2 = [multiply(G2, tau**i) for i in range(deg-1,-1,-1)] 
    test =  [tau**i for i in range(deg-1,-1,-1)] 
    print("test", test)

    # t = (x-1)(x-2)...(x-n)
    t = galois.Poly([1], field=GF)
    for i in range(1, deg+1):
        t = t * galois.Poly([1, p-i], field=GF)

    t_at_tau = t(tau)
    print("t_at_tau", t_at_tau)
    ht_srs = [multiply(G1, int(t_at_tau * tau**i)) for i in range(deg-2,-1,-1)]
    test2 = [int(t_at_tau * tau**i) for i in range(deg-2,-1,-1)]
    print("test2", test2)

    W_galois = matrix_GF(w.matrix)
    
    Lw_qap = galois.Poly([0], field=GF)
    for i in range(w.size):
        Lw_qap = L_polys[i]*W_galois[i] + Lw_qap

    Rw_qap = galois.Poly([0], field=GF)
    for i in range(w.size):
        Rw_qap = R_polys[i]*W_galois[i] + Rw_qap
    
    Ow_qap = galois.Poly([0], field=GF)
    for i in range(w.size):
        Ow_qap = O_polys[i]*W_galois[i] + Ow_qap

    t, h_qap = calculate_t_and_h(Lw_qap, Rw_qap, Ow_qap, r1cs.L.rows_size)

    print("Lw_qap", Lw_qap)
    print("Rw_qap", Rw_qap)
    print("Ow_qap", Ow_qap)
    print("h_qap", h_qap)
    print("t", t)

    # Sanity Check 
    assert(Lw_qap*Rw_qap == Ow_qap + h_qap*t)

    Lw_qap_coeffs = [int(x) for x in Lw_qap.coeffs]
    Rw_qap_coeffs = [int(x) for x in Rw_qap.coeffs]
    Ow_qap_coeffs = [int(x) for x in Ow_qap.coeffs]
    h_qap_coeffs = [int(x) for x in h_qap.coeffs]

    print("Lw_qap_coeffs", Lw_qap_coeffs)
    print("Rw_qap_coeffs", Rw_qap_coeffs)
    print("Ow_qap_coeffs", Ow_qap_coeffs)
    print("h_qap_coeffs", h_qap_coeffs)

    # Test for inner product
    mockpoly = [multiply(G1, 1) for i in range(deg,-1,-1)]
    assert inner_product(mockpoly, Lw_qap_coeffs) == multiply(G1, sum(Lw_qap_coeffs))
    # TODO : Do i need to handle case when [G1,G1,G1][0,6,14] = 20G1 but the GF is 17, so do i need to reduce ?

    # -----------

    alpha = 5
    beta = 5
    alpha_g1 = multiply(G1, alpha)
    beta_g2 = multiply(G2, beta)

    # Issue is in this part of code. Hmmmmmm
    phi_srs = []
    for i in range(len(L_polys)):
        phi_element = alpha * L_polys[i](tau) + beta * R_polys[i](tau) + O_polys[i](tau)
        phi_element = multiply(G1, int(phi_element))
        phi_srs.append(phi_element)

    C_mid = [multiply(phi_srs[i],int(W_galois[i])) for i in range(len(phi_srs))]
    C_part1 = None
    for i in range(len(C_mid)):
        C_part1 = add(C_part1,C_mid[i])


    C_part2 = inner_product(ht_srs, h_qap_coeffs) 
    C = add(C_part1,C_part2)

    A = add(alpha_g1, inner_product(tau_g1, Lw_qap_coeffs))
    B = add(beta_g2, inner_product(tau_g2, Rw_qap_coeffs)) 

    print("Prover A,B,C generation complete")

    lhs = pairing(B, neg(A))
    rhs_1 = pairing(beta_g2, alpha_g1)
    rhs_2 = pairing(G2, C)

    # the equivalent of adding in exponent space would be multiply 
    assert (lhs * rhs_1 * rhs_2) == FQ12.one(), f" lhs should match rhs"
    print("Done")

    return A, B, C

# Coeffs and points are in a degree descending order 
# Eg  : poly = 4x^2 + 7x + 8, points = [tau^2, tau, 1]
# Then we have : coeffs = [4, 7, 8]
# and we have : points = [tau^2, tau, 1]
# So the result is : 4*tau^2 + 7*tau + 8*1 


def inner_product(points, coeffs):
    # TODO: Handle when coeffs may be smaller in size than points
    return reduce(add, map(multiply, points, coeffs))

# def inner_product(poly, value_arr):
#     """
#     Calculate the inner product of a polynomial and an array of values
#     """
#     # The arr given will be of type [G1, tau*G1 .... tau^n-1*G1]
#     coefficients = poly.coeffs # in a degree descending order 
#     # tau_g1 has elements in G1 with field order as field_modulus and coeff[i] are in modulus p
#     # For now doing int for faster calculation and converting to field_modulus since p is significantly lower
#     print("coefficients", coefficients)
#     coefficients = coefficients[::-1]  # Reverse the coefficients array
    
#     result = None
#     for i in range(len(value_arr)):
#         if result is None:
#             result = multiply(value_arr[i], int(coefficients[i]))
#         else:
#             result = add(result, multiply(value_arr[i], int(coefficients[i])))
#     # print("type of result", type(result), type(result[0]))
#     return result