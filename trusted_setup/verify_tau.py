from py_ecc.bn128 import pairing, multiply

def verify_tau(alpha_g1, beta_g2, tau_g1, tau_g2, ht_srs, phi_srs, L_polys, R_polys, O_polys, deg):
    # TODO : alpha_g1 is point on curve 
    # For tau_g1 and tau_g2
    if (deg - 1) % 2 == 0:  # Even case
        for i in range((deg - 1) // 2 + 1):
            a = pairing(tau_g2[deg-1-i], tau_g1[i])
            b = pairing(tau_g2[i], tau_g1[deg-1-i])
            assert(a == b)
    else:
        for i in range((deg - 1) // 2 + 1):
            a = pairing(tau_g2[deg-1-i], tau_g1[i])
            b = pairing(tau_g2[i], tau_g1[deg-1-i])
            assert(a == b)

    print("Tau verification complete")