import numpy as np
import galois
from utils.matrix import matrix_GF
from utils.lagrange import get_r1cs_interpolated_polys
from utils.constants import p, GF

def r1cs_to_qap(r1cs):
    """
    Converts R1CS system to QAP Polynomials
    """
    
    # Convert matrix elements to Galois Field elements
  
    L_galois = matrix_GF(r1cs.L.matrix)
    R_galois = matrix_GF(r1cs.R.matrix)
    OUT_galois = matrix_GF(r1cs.OUT.matrix)

    L_polys, R_polys, O_polys = get_r1cs_interpolated_polys(L_galois, R_galois, OUT_galois, r1cs.L.rows_size)

    return L_polys, R_polys, O_polys

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
    for i in range(1, deg+1):
        t = t * galois.Poly([1, p-i], field=GF)

    h = (Lw_qap*Rw_qap - Ow_qap) // t

    # Check cause above doesn't return remainder
    assert(Lw_qap*Rw_qap == Ow_qap + h*t)

    return t,h

