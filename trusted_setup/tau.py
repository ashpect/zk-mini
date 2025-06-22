import random
import galois
from utils.constants import p, GF
from py_ecc.bn128 import G1, G2, multiply
from r1cs.qap import r1cs_to_qap

def grothv0(r1cs):
    """
    Input : r1cs
    Output : 
        alpha_g1 - a tuple of elements in G1
        beta_g2 - a tuple of elements in G2
        tau_g1 - a tuple of elements in G1 
        tau_g2 - a tuple of elements in G2 
        ht_srs - a tuple of elements in G1 
        phi_arr - a tuple of elements in G1 
        L_polys - a tuple of Galois Polynomials
        R_polys - a tuple of Galois Polynomials
        O_polys - a tuple of Galois Polynomials
    """
    L_polys, R_polys, O_polys = r1cs_to_qap(r1cs)
    print("R1CS to QAP Polys Done")
    # print(" Type of polys", type(L_polys[0]), type(R_polys), type(O_polys[0])) # Debugging
    deg = r1cs.L.rows_size

    alpha, beta, tau = generate_random_scalars()
    # print(" Type of alpha, beta, tau", type(alpha), type(beta), type(tau)) # Debugging
    tau = GF(25) # for testing

    tau_g1, tau_g2, ht_srs = generate_srs(tau, deg)
    print("SRS Generated")

    phi_arr = generate_phi(L_polys, R_polys, O_polys, alpha, beta, tau)
    print("Phi Generated")

    alpha_g1 = multiply(G1, int(alpha))
    beta_g2 = multiply(G2, int(beta))

    print("Grothv0 trusted setup complete")

    return alpha_g1, beta_g2, tau_g1, tau_g2, ht_srs, phi_arr, L_polys, R_polys, O_polys
    

def generate_random_scalars():
    """
    Generate random scalars for the trusted setup
    """
    alpha = random.randint(1, p-1)
    beta = random.randint(1, p-1)
    tau = random.randint(1, p-1)

    return GF(alpha), GF(beta), GF(tau)
    
def generate_srs(tau, deg):
    """
    Generate the SRS
    Input : tau is in GF(p)
    Output : 
        tau_g1 - a tuple of elements in G1
        tau_g2 - a tuple of elements in G2
        ht_srs - a tuple of elements in G1
    """
    tau_powers = [] # All elements in GF(p)
    for i in range(deg):
        tau_powers.append(pow(tau, i)) # modulus tau^i % p and return value in GF(p) as well
    # if z is in GF(p) and i is in Z, pow(z,i)) is equivalent to pow(z,i,p)

    tau_g1 = []
    for i in range(deg):
        tau_g1.append(multiply(G1, int(tau_powers[i]))) # Need to convert GF to int for multiply op
        # Even if you use p = field_modulus, the mult wont work.

    tau_g2 = []
    for i in range(deg):
        tau_g2.append(multiply(G2, int(tau_powers[i])))

    ht_srs = [] # srs for h(tau) * t(tau) : deg n-2

    # t = (x-1)(x-2)...(x-(deg-1))
    t = galois.Poly([1], field=GF)
    for i in range(1, deg):
        t = t * galois.Poly([1, p-i], field=GF)

    # print(t) - DEBUGGING
    t_at_tau = t(tau) # in GF(p)
    # print("t_at_tau", t_at_tau, type(t_at_tau)) - DEBUGGING

    t_tau_in_G1 = multiply(G1, int(t_at_tau))
    for i in range(deg-1):
        ht_srs.append(multiply(t_tau_in_G1, int(tau_powers[i])))

    # print(ht_srs)
    # print(t_at_tau)

    return tau_g1, tau_g2, ht_srs

def generate_phi(L_polys, R_polys, O_polys, alpha, beta, tau):
    """ 
    Generate the phi polynomial
    Input : 
    - L_polys, R_polys, O_polys are Galois Polynomials
    - alpha, beta, tau are Galois Field Elements
    Output : 
    - phi_arr - a tuple of elements in G1
    """
    phi_arr = []
    # print("Type of L_polys", type(L_polys[0])) # Galois Polynomial
    # print("Type of L_polys[0](tau)", type(L_polys[0](tau))) # Galois Field Element

    for i in range(len(L_polys)):
        phi_element = alpha * L_polys[i](tau) + beta * R_polys[i](tau) + O_polys[i](tau) # All addition, mult is in GF(p)
        phi_element = multiply(G1, int(phi_element)) # Convert GF to int for multiply op
        phi_arr.append(phi_element)

    return phi_arr


