import galois
from py_ecc.bn128 import curve_order

# Global prime for Galois Field
p = curve_order
GF = galois.GF(p) 