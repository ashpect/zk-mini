import numpy as np
import galois
from utils.constants import p, GF

def interpolate_column(col, num_rows):
    xs = GF(np.arange(1, num_rows + 1))
    # print(galois.lagrange_poly(xs, col))
    return galois.lagrange_poly(xs, col)

def get_r1cs_interpolated_polys(L_galois, R_galois, O_galois, num_rows):
    """
    Takes in r1cs matrices and returns the lagrange polynomial
    """
    # axis 0 is the columns. (iterates over columns)
    # apply_along_axis is the same as doing a for loop over the columns 
    # and collecting the results in an array
    U_polys = np.apply_along_axis(lambda col: interpolate_column(col, num_rows), 0, L_galois)
    V_polys = np.apply_along_axis(lambda col: interpolate_column(col, num_rows), 0, R_galois)
    W_polys = np.apply_along_axis(lambda col: interpolate_column(col, num_rows), 0, O_galois)

    return U_polys, V_polys, W_polys
