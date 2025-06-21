import numpy as np
from py_ecc.bn128 import G2, pairing, eq, FQ12
import galois
from utils.constants import GF,p

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

def verifier_qap_tau_random():
    """
    Generate a random tau for the QAP verifier
    """
    return np.random.randint(1, p)

def verify_qap(L, R, O, t, h):
    """
    Verify the proof
    """
    assert(L*R == O + h*t)
    print("Verification successful")