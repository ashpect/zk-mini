# QAP Verifier

import numpy as np
from utils.constants import p

def verifier_qap_tau_random():
    """
    Verifier generates a random tau for the prover
    """
    return np.random.randint(1, p)

def verify_qap(L, R, O, t, h):
    """
    Verifier verifies the proof
    """
    assert(L*R == O + h*t)
    print("Verification successful")