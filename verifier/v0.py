# Basic R1CS Verifier based on EC bilinear pairings

import numpy as np
from py_ecc.bn128 import G2, pairing, eq, FQ12

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