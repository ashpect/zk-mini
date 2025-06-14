import numpy as np
from py_ecc.bn128 import G1, G2, pairing, add, multiply, eq, curve_order
from util import matrix_multiply

# Eqn :  z = x*x*y + 1
# Constraints:
# 1) a = x*x
# 2) -1 + z = a*y

# OUT = Lw*Rw (element wise - *)
OUT = np.array([
    [0, 0, 0, 0, 1],
    [-1, 1, 0, 0, 0],
])

L = np.array([
    [0, 0, 1, 0, 0],
    [0, 0, 0, 0, 1],
])

R = np.array([
    [0, 0, 1, 0, 0],
    [0, 0, 0, 1, 0],
])

# Witness vector is [1, z , x, y, a] ( for x = 3, y = 2, z = 18, a = 9)
W = [1, 19, 3, 2 , 9]

# verify once
assert(np.all(np.matmul(OUT, W) == np.matmul(L, W) * np.matmul(R, W)))

# Step 1 : Converting witness vector to Elliptic Curve Point
Witness_ec1 = np.array([
    multiply(G1, W[0]),
    multiply(G1, W[1]),
    multiply(G1, W[2]),
    multiply(G1, W[3]),
    multiply(G1, W[4]),
])

Witness_ec2 = np.array([
    multiply(G2, W[0]),
    multiply(G2, W[1]),
    multiply(G2, W[2]),
    multiply(G2, W[3]),
    multiply(G2, W[4]),
])

# Step 2 : Finding E1(Lw), E2(Rw), E1(Ow)
Lw = np.array([G1,G1])
for j in range(2):
    point = G1
    for i in range(5):
        point = add(point, multiply(Witness_ec1[i], L[j][i]))
    Lw[j] = point

Rw = np.array([G2,G2])
for j in range(2):
    point = G2
    for i in range(5):
        point = add(point, multiply(Witness_ec2[i], R[j][i]))
    Rw[j] = point

Ow_g1 = np.array([G1,G1])
for j in range(2):
    point = G1
    for i in range(5):
        if OUT[j][i] >= 0:
            point = add(point, multiply(Witness_ec1[i], OUT[j][i]))
        else:
            point = add(point, multiply(Witness_ec1[i], (curve_order-1))) # To handle -1 case, generalize it later.

    Ow_g1[j] = point


# Verifier Steps : 
# Step 1 : Compute E2(E1(Ow)) 
Ow_g2 = np.array([G2,G2])
for i in range(2):
    Ow_g2[i] = pairing(G2, Ow_g1[i])

# Step 2 : Computer Bilinear pairing of Lw and Rw
Lw_Rw = np.array([pairing(Rw[0],Lw[0]), pairing(Rw[1],Lw[1])])

# Step 3 : Assert E2(E1(Ow)) == E2(E1(Lw*Rw))
print("Lw_Rw:", Lw_Rw)
print("Ow_g2:", Ow_g2)
# assert(np.all(Ow_g2 == Lw_Rw))


