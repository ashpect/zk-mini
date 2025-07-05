from verifier.verifier import verify, verify_qap, verify_grothv0
from prover.p0 import prover
from prover.p1 import prover_qap_no_trusted_setup
from prover.p2 import prover_grothv0, qap_trusted
from trusted_setup.tau import grothv0, generate_random_scalars
from trusted_setup.verify_tau import verify_tau
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
    _,_,tau = generate_random_scalars()

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
    alpha_g1, beta_g2, tau_g1, tau_g2, ht_srs, phi_srs, L_polys, R_polys, O_polys = grothv0(r1cs)

    # Step 2 : Verify the trusted outputs (ideally done inside prover)
    verify_tau(alpha_g1, beta_g2, tau_g1, tau_g2, ht_srs, phi_srs, L_polys, R_polys, O_polys, r1cs.L.rows_size)

    # Prover computation to generate proof
    A, B, C = prover_grothv0(r1cs, witness)

    # Step 2 : Verify the proof
    # verify_grothv0(A, B, C, alpha_g1, beta_g2)

def poc_main(example_number):
    r1cs, witness = get_mock_example(example_number)
    print("POC STARTING")
    qap_trusted(r1cs, witness)
    print("POC DONE")

if __name__ == "__main__":
    # test_all()    
    # poc_main(1)
    v2_main(1)
