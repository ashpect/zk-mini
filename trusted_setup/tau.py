import random
import galois
from utils.constants import p, GF
from py_ecc.bn128 import G1, G2, multiply
from r1cs.qap import r1cs_to_qap

def grothv0(r1cs):
    """
    Returns alpha_g1, beta_g1, tau_g1, tau_g2, ht_srs, phi_arr
    """
    L_polys, R_polys, O_polys = r1cs_to_qap(r1cs)
    deg = r1cs.L.rows_size

    alpha, beta, tau = generate_random_scalars()
    tau = 25 # for testing

    tau_g1, tau_g2, ht_srs = generate_srs(tau, deg)
    phi_arr = generate_phi(L_polys, R_polys, O_polys, alpha, beta, tau)

    alpha_g1 = multiply(G1, alpha)
    beta_g1 = multiply(G1, beta)

    return alpha_g1, beta_g1, tau_g1, tau_g2, ht_srs, phi_arr

def generate_random_scalars():
    """
    Generate random scalars for the trusted setup
    """
    alpha = random.randint(1, p-1)
    beta = random.randint(1, p-1)
    tau = random.randint(1, p-1)

    return alpha, beta, tau
    
def generate_srs(tau, deg):
    """
    Generate the SRS
    """
    tau_powers = []
    for i in range(deg):
        tau_powers.append(pow(tau, i, p)) # from tau to tau^(deg-1)

    tau_g1 = []
    for i in range(deg):
        tau_g1.append(multiply(G1, tau_powers[i])) # modulus tau^i % p
    
    tau_g2 = []
    for i in range(deg):
        tau_g2.append(multiply(G2, tau_powers[i]))

    ht_srs = [] # srs for h(tau) * t(tau) : deg n-2

    # t = (x-1)(x-2)...(x-(deg-1))
    t = galois.Poly([1], field=GF)
    for i in range(1, deg):
        t = t * galois.Poly([1, p-i], field=GF)

    t_at_tau = t(tau)
    
    for i in range(deg-1):
        ht_srs.append(t_at_tau * tau_powers[i])

    # print(ht_srs)
    # print(t_at_tau)

    return tau_g1, tau_g2, ht_srs

def generate_phi(L_polys, R_polys, O_polys, alpha, beta, tau):
    """ 
    Generate the phi polynomial
    """
    phi_arr = []
    print(alpha, beta)
    print(L_polys[0](tau))
    print(R_polys[0](tau))
    print(O_polys[0](tau))
    print(multiply(G1, 0))
    for i in range(len(L_polys)):
        phi_element = alpha * L_polys[i](tau) + beta * R_polys[i](tau) + O_polys[i](tau)
        phi_element = multiply(G1, phi_element)
        phi_arr.append(phi_element)

    return phi_arr


