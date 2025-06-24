import numpy as np
from utils.constants import p
from py_ecc.bn128 import pairing, G2, eq, add, neg, FQ12

def verify(Ow_g1, Rw, Lw):
    """
    Version 0 : Verify the proof
    """

    # Compute E2(E1(Ow)) 
    Ow_g2 = np.array([FQ12.one(), FQ12.one()])
    for i in range(2):
        Ow_g2[i] = pairing(G2, Ow_g1[i])

    # Compute Bilinear pairing of Lw and Rw
    Lw_Rw = np.array([FQ12.one(), FQ12.one()])
    for i in range(2):
        Lw_Rw[i] = pairing(Rw[i], Lw[i])

    # Assert E2(E1(Ow)) == E2(E1(Lw*Rw))
    # print("Ow_g2:", Ow_g2)
    # print("Lw_Rw:", Lw_Rw)
    assert(eq(Ow_g2[0], Lw_Rw[0]) and eq(Ow_g2[1], Lw_Rw[1]))
    print("Verification successful")

def verify_qap(L, R, O, t, h):
    """
    Verifier verifies the proof
    """
    assert(L*R == O + h*t)
    print("Verification successful")

def verify_grothv0(A, B, C, alpha_g1, beta_g2):
    """
    Verifier verifies the proof
    Input : 
    - A, B, C - Elements in G1, G2, G1
    - alpha_g1, beta_g2 - Elements in G1, G2
    Output :
    - True if proof is valid, False otherwise
    """
    lhs = pairing(B, neg(A))
    rhs = pairing(beta_g2, alpha_g1) * pairing(G2, C)
    
    print("Pairing check:", lhs * rhs == FQ12.one())
