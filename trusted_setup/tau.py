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
        phi_srs - a tuple of elements in G1 
        L_polys - a tuple of Galois Polynomials
        R_polys - a tuple of Galois Polynomials
        O_polys - a tuple of Galois Polynomials
    """
    L_polys, R_polys, O_polys = r1cs_to_qap(r1cs)
    print("R1CS to QAP Polys Done")
    # print("L_polys", L_polys)
    # print("R_polys", R_polys)
    # print("O_polys", O_polys)
    # print(" Type of polys", type(L_polys[0]), type(R_polys), type(O_polys[0])) # Debugging
    deg = r1cs.L.rows_size

    alpha, beta, tau = generate_random_scalars()
    # print(" Type of alpha, beta, tau", type(alpha), type(beta), type(tau)) # Debugging
    tau = GF(5) # for testing
    beta = GF(10) # for testing
    alpha = GF(15) # for testing

    tau_g1, tau_g2, ht_srs = generate_srs(tau, deg)
    print("SRS Generated")

    phi_srs = generate_phi(L_polys, R_polys, O_polys, alpha, beta, tau)
    print("Phi Generated")

    alpha_g1 = multiply(G1, int(alpha))
    beta_g2 = multiply(G2, int(beta))

    print("Grothv0 trusted setup complete")

    return alpha_g1, beta_g2, tau_g1, tau_g2, ht_srs, phi_srs, L_polys, R_polys, O_polys
    
def generate_srs(tau, deg):
    """
    Generate the SRS
    Input : tau is in GF(p)
    Output : 
        tau_g1 - a tuple of elements in G1
        tau_g2 - a tuple of elements in G2
        ht_srs - a tuple of elements in G1
    # We are taking deg as input but generally the srs are calultaed over sufficient number of powers of tau
    """
    # # Not needed 
    # tau_powers = [] # All elements in GF(p)
    # for i in range(deg):
    #     tau_powers.append(pow(tau, i)) # modulus tau^i % p and return value in GF(p) as well
    # # if z is in GF(p) and i is in Z, pow(z,i)) is equivalent to pow(z,i,p)

    tau = int(tau) # Cause multiply func needs it to be int, TODO: maybe int not needed if same curve order, see later

    # (tau^n-1)G1, (tau^n-2)G1, ... , 1G1
    tau_g1 = [multiply(G1, tau**i) for i in range(deg-1,-1,-1)] # deg-1 to 0 in descending order
    # for deg = 3, we have [tau^2, tau, 1]

    # (tau^n-1)G2, (tau^n-2)G2, ... , 1G2
    tau_g2 = [multiply(G2, tau**i) for i in range(deg-1,-1,-1)] 

    # t = (x-1)(x-2)...(x-n)
    t = galois.Poly([1], field=GF)
    for i in range(1, deg+1):
        t = t * galois.Poly([1, p-i], field=GF)

    # print("t",t) # DEBUGGING
    t_at_tau = t(tau) # in GF(p)
    # print("t_at_tau", t_at_tau, type(t_at_tau)) # DEBUGGING
   
    # (tau^n-2)(t(tau))G1, (tau^n-3)(t(tau))G1, ... , 1G1
    ht_srs = [multiply(G1, int(t_at_tau * tau**i)) for i in range(deg-2,-1,-1)] # deg-2 to 0 in descending order

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

# TODO : normal generate and then GF return and is GF necessary ?
def generate_random_scalars():
    """
    Generate random scalars for the trusted setup
    """
    alpha = random.randint(1, p-1)
    beta = random.randint(1, p-1)
    tau = random.randint(1, p-1)

    return GF(alpha), GF(beta), GF(tau)