from py_ecc.bn128 import pairing, G2, eq

def verify_grothv0(A, B, C, alpha_g1, beta_g2):
    """
    Verifier verifies the proof
    Input : 
    - A, B, C - Elements in G1, G2, G1
    - alpha_g1, beta_g2 - Elements in G1, G2
    Output :
    - True if proof is valid, False otherwise
    """
    lhs = pairing(B, A)
    rhs = pairing(beta_g2, alpha_g1) + pairing(G2, C)

    assert(eq(lhs, rhs))

    print("Verifier A,B,C verification complete")

