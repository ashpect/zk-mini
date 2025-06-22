from verifier.v0 import verify
from verifier.v1 import verify_qap, verifier_qap_tau_random
from verifier.v2 import verify_grothv0
from prover.p0 import prover
from prover.p1 import prover_qap_no_trusted_setup
from prover.p2 import prover_grothv0
from trusted_setup.tau import grothv0
from mock_examples import get_mock_example

def v0_main(example_number):
    # Step 1 : Mock Example
    r1cs, witness = get_mock_example(example_number)

    # Step 2 : Prover
    Ow_g1, Rw, Lw = prover(r1cs, witness)

    # Step 3 : Verify the proof
    verify(Ow_g1, Rw, Lw)

def v1_main(example_number):
    # Step 0 : Mock Example
    r1cs, witness = get_mock_example(example_number)

    # Step 1 : Generate a random tau
    tau = verifier_qap_tau_random()

    # Step 2 : Prover
    L, R, O, t, h = prover_qap_no_trusted_setup(r1cs, witness, tau)

    # Step 3 : Verify the proof
    verify_qap(L, R, O, t, h)

def test_all():
    for example_type in range(1, 3):
        v0_main(example_type)
        v1_main(example_type)

def v2_main(example_number):

    # Step 0 : Mock Example
    r1cs, witness = get_mock_example(example_number)

    # Step 1 : Trusted setup generates public parameters
    alpha_g1, beta_g2, tau_g1, tau_g2, ht_srs, phi_arr, L_polys, R_polys, O_polys = grothv0(r1cs)

    # Prover computation to generate proof
    A, B, C = prover_grothv0(r1cs, alpha_g1, beta_g2, tau_g1, tau_g2, ht_srs, phi_arr, L_polys, R_polys, O_polys, witness)

    # Step 2 : Verify the proof
    # verify_grothv0(A, B, C, alpha_g1, beta_g2)

if __name__ == "__main__":
    #test_all()
    v2_main(1)
