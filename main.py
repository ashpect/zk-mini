import numpy as np
from py_ecc.bn128 import G1, G2, pairing, add, multiply, eq, curve_order, FQ12
from verifier.v0 import verify
from verifier.v1 import verify_qap, verifier_qap_tau_random
from prover.p0 import prover
from prover.p1 import prover_qap
from r1cs import R1CS, R1cs_matrix, Witness
from mock_examples import mock_example_1

def v0_main():
    # Step 1 : Mock Example
    r1cs = mock_example_1()

    # Step 2 : Prover
    Ow_g1, Rw, Lw = prover(r1cs)

    # Step 3 : Verify the proof
    verify(Ow_g1, Rw, Lw)

def v1_main():
    # Step 0 : Mock Example
    r1cs = mock_example_1()

    # Step 1 : Generate a random tau
    tau = verifier_qap_tau_random()

    # Step 2 : Prover
    L, R, O, t, h = prover_qap(r1cs, tau)

    # Step 3 : Verify the proof
    verify_qap(L, R, O, t, h)

if __name__ == "__main__":
    v0_main()
    v1_main()

