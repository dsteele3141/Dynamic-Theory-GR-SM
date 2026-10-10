#!/usr/bin/env sage
# 5D Dirac/Clifford algebra (Eqs 108-109)
#
# Source: P. E. Williams, LA-8370-MS (LASL, 1980), eqs (108)-(109). The Eq-109
#   literal A-matrix fails anticommutation; the corrected A is [[0,iI2],[-iI2,0]].
# Full derivation, context and provenance:
#   Supplementary_Material.pdf — §S2, Appendix A [C6].

from sage.all import *

def ok(label, cond):
    print(f"  {label}: {'PASS' if cond else 'FAIL'}")
    assert cond

def check_clifford(gammas, g_diag):
    I8 = identity_matrix(SR, 8)
    for j in range(5):
        for k in range(j, 5):
            anticomm = (gammas[j]*gammas[k] + gammas[k]*gammas[j]).simplify_full()
            expected = 2 * g_diag[j] * I8 if j == k else zero_matrix(SR, 8)
            if anticomm != expected:
                return False, f"{{gamma{j}, gamma{k}}} failed"
    return True, ""

# 1. BASIC BLOCKS (2x2)
I2 = identity_matrix(SR, 2)
Z2 = zero_matrix(SR, 2)
sigma1 = matrix(SR, [[0, 1], [1, 0]])
sigma2 = matrix(SR, [[0, -I], [I, 0]])
sigma3 = matrix(SR, [[1, 0], [0, -1]])

# 2. STANDARD 4x4 DIRAC MATRICES
beta = block_matrix([[I2, Z2], [Z2, -I2]], subdivide=False)
alpha1 = block_matrix([[Z2, sigma1], [sigma1, Z2]], subdivide=False)
alpha2 = block_matrix([[Z2, sigma2], [sigma2, Z2]], subdivide=False)
alpha3 = block_matrix([[Z2, sigma3], [sigma3, Z2]], subdivide=False)

# 3. THE 5th GENERATOR (A-matrix)
A_source  = block_matrix([[Z2, I2], [I2, Z2]], subdivide=False)
A_correct = block_matrix([[Z2, I*I2], [-I*I2, Z2]], subdivide=False)

# 4. FULL 8x8 MATRICES
Z4 = zero_matrix(SR, 4)
betaTilde = block_matrix([[beta, Z4], [Z4, beta]], subdivide=False)
alphaTilde1 = block_matrix([[Z4, alpha1], [alpha1, Z4]], subdivide=False)
alphaTilde2 = block_matrix([[Z4, alpha2], [alpha2, Z4]], subdivide=False)
alphaTilde3 = block_matrix([[Z4, alpha3], [alpha3, Z4]], subdivide=False)

alphaTilde4_source  = block_matrix([[A_source,  Z4], [Z4,  A_source]], subdivide=False)
alphaTilde4_correct = block_matrix([[A_correct, Z4], [Z4,  A_correct]], subdivide=False)

def get_gammas(aT4):
    aT = [alphaTilde1, alphaTilde2, alphaTilde3, aT4]
    return [betaTilde] + [-betaTilde * aT[i] for i in range(4)]

gammas_source  = get_gammas(alphaTilde4_source)
gammas_correct = get_gammas(alphaTilde4_correct)

g_signature = [1, -1, -1, -1, -1]

print("="*70)
print("02_dirac_algebra_adjudicated.sage")
print("="*70)

print("\n(A) Diagnostic: Literal Manuscript A-matrix [[0,I],[I,0]]")
res_source, msg_source = check_clifford(gammas_source, g_signature)
print(f"  Clifford Algebra PASS: {res_source}")
if not res_source:
    print(f"  Note: Literal source failed as expected - {msg_source}")

print("\n(B) Verification: Corrected Complex A-matrix [[0,iI],[-iI,0]]")
res_correct, msg_correct = check_clifford(gammas_correct, g_signature)
ok("Corrected A-matrix satisfies 5D Clifford algebra", res_correct)

print("\n(C) A-matrix Properties (Corrected)")
ok("A^2 = I", A_correct**2 == identity_matrix(SR, 4))
ok("A is Hermitian", A_correct == A_correct.conjugate().transpose())
ok("{A, beta} = 0", (A_correct*beta + beta*A_correct) == zero_matrix(SR, 4))
ok("{A, alpha1} = 0", (A_correct*alpha1 + alpha1*A_correct) == zero_matrix(SR, 4))

print("\n(D) Bianchi Setup Check (Eq 68)")
var('E1 E2 E3 B1 B2 B3 V0 V1 V2 V3')
F = matrix(SR, [
    [ 0,    I*E1,  I*E2,  I*E3,  I*V0],
    [-I*E1, 0,     B3,   -B2,    V1  ],
    [-I*E2,-B3,    0,     B1,    V2  ],
    [-I*E3, B2,   -B1,   0,     V3  ],
    [-I*V0,-V1,   -V2,  -V3,    0   ]
])
ok("Field Tensor F is antisymmetric", (F + F.transpose()) == zero_matrix(SR, 5))

print("\nSummary: Corrected 5D Dirac algebra verified with complex A-matrix.")
print("="*70)
