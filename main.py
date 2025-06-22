from verifier.v0 import verify
from verifier.v1 import verify_qap, verifier_qap_tau_random
from prover.p0 import prover
from prover.p1 import prover_qap
from mock_examples import get_mock_example

def v0_main(example_number):
    # Step 1 : Mock Example
    r1cs = get_mock_example(example_number)

    # Step 2 : Prover
    Ow_g1, Rw, Lw = prover(r1cs)

    # Step 3 : Verify the proof
    verify(Ow_g1, Rw, Lw)

def v1_main(example_number):
    # Step 0 : Mock Example
    r1cs = get_mock_example(example_number)

    # Step 1 : Generate a random tau
    tau = verifier_qap_tau_random()

    # Step 2 : Prover
    L, R, O, t, h = prover_qap(r1cs, tau)

    # Step 3 : Verify the proof
    verify_qap(L, R, O, t, h)

def test_all():
    for example_type in range(1, 3):
        v0_main(example_type)
        v1_main(example_type)

if __name__ == "__main__":
    test_all()

