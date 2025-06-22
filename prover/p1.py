# QAP Prover
import galois
from r1cs.qap import r1cs_to_qap, calculate_t_and_h
from utils.matrix import matrix_GF
from utils.constants import p, GF

def prover_qap_no_trusted_setup(r1cs, w, tau):
    """
    Prover does the following :
    1. Get values of all poly at tau
    2. Return the values
    We are trusting the prover to do correct poly operations.
    """ 

    # Calulate W in galois field
    W_galois = matrix_GF(w.matrix)

    # Get the R1CS to QAP Polynomials
    L_polys, R_polys, O_polys = r1cs_to_qap(r1cs)

    # Random poly of deg n-1 in GF(p)
    polynomial = galois.primitive_poly(p, r1cs.L.rows_size-1 , method="random")

    # Final Lw , initiating with random poly and size m(colms) for example
    Lw_qap = galois.Poly([0], field=GF)
    for i in range(w.size):
        Lw_qap = L_polys[i]*W_galois[i] + Lw_qap

    Rw_qap = galois.Poly([0], field=GF)
    for i in range(w.size):
        Rw_qap = R_polys[i]*W_galois[i] + Rw_qap
    
    Ow_qap = galois.Poly([0], field=GF)
    for i in range(w.size):
        Ow_qap = O_polys[i]*W_galois[i] + Ow_qap

    # Get the t and h polynomials
    t, h = calculate_t_and_h(Lw_qap, Rw_qap, Ow_qap, r1cs.L.rows_size)

    # Get values of all poly at tau
    Lw_qap_tau = Lw_qap(tau)
    Rw_qap_tau = Rw_qap(tau)
    Ow_qap_tau = Ow_qap(tau)
    t_tau = t(tau)
    h_tau = h(tau)

    return Lw_qap_tau, Rw_qap_tau, Ow_qap_tau, t_tau, h_tau