#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE (primary): LA-8370-MS (1980) Eqs 109-112: Williams' 8x8 Dirac
#   representation (two inequivalent 4x4 blocks, kappa = +1 and -1) and the
#   5D equation (i gamma^j d_j + 1) psi = 0 on the mass axis x^4 = gamma/a0.
# Context: noflux_boundary_derivation (chiral edge, ledger Q4),
#   weak_neutral_current (N_nu), strong_sector_su3 (commutants, Q14).
# ---------------------------------------------------------------------------
#
# su2_from_5d.sage                                     [A] (assertion-hardened)
#
# PURPOSE.  Can SU(2) be DERIVED from DT's 5D structure, and does it survive
# the switch u4 = 0 (i.e. is it in the recovered 4D sector, not a relaxation)?
#
#   [A] geometry: one extra dimension cannot carry SU(2) (Kaluza-Klein route
#       closed)
#   [B] one 4x4 irrep: internal symmetry U(1) only (psi treated as REAL, so
#       particle-antiparticle mixings are included); DT's mass kills the
#       massless-only extras
#   [C] Williams' 8x8: internal symmetry su(2) + u(1) -- NATIVE, non-abelian
#   [D] it survives the switch: it acts on the zero modes; the edge condition
#       invariant under it and 4D Lorentz is unique up to the hand; the
#       doublet is (block-A zero mode, conjugate of block-B zero mode), both
#       LEFT-handed: a chiral doublet
#   [E] charge: Q = T3 + Y with Y the native u(1); one neutral member and
#       integer charge (Ze) force |Y| = 1/2, i.e. (nu_L, e_L).  Williams' own
#       common phase is pure T3 (no neutral member)
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("su2_from_5d.sage -- is SU(2) in DT's 5D structure?")
print("=" * 74)

K = CyclotomicField(4); i_ = K.gen()
I2, Z2, Z4 = identity_matrix(K, 2), zero_matrix(K, 2), zero_matrix(K, 4)
s = [matrix(K, [[0, 1], [1, 0]]), matrix(K, [[0, -i_], [i_, 0]]),
     matrix(K, [[1, 0], [0, -1]])]
B = lambda a, b, c, d: block_matrix([[a, b], [c, d]], subdivide=False)
g = [B(I2, Z2, Z2, -I2)] + [B(Z2, si, -si, Z2) for si in s]
g5 = i_ * g[0] * g[1] * g[2] * g[3]
G4 = g + [i_ * g5]                                     # one 4x4 irrep (Q3)
beta = B(I2, Z2, Z2, -I2); alph = [B(Z2, si, si, Z2) for si in s]
Amat = B(Z2, i_ * I2, -i_ * I2, Z2); bT = B(beta, Z4, Z4, beta)
aT = [B(Z4, a, a, Z4) for a in alph] + [B(Amat, Z4, Z4, Amat)]
G8 = [bT] + [-bT * a for a in aT]                      # Williams Eq 111
eta = [1, -1, -1, -1, -1]
for GG, nm in ((G4, "4x4"), (G8, "8x8")):
    n = GG[0].nrows()
    check(all(GG[j] * GG[k] + GG[k] * GG[j] == (2 * eta[j] if j == k else 0) * identity_matrix(K, n)
              for j in range(5) for k in range(5)), "%s: 5D Clifford algebra" % nm)

def realify(A):
    n = A.nrows()
    R = matrix(QQ, n, n, lambda a, b: A[a, b].list()[0])
    S = matrix(QQ, n, n, lambda a, b: A[a, b].list()[1])
    return block_matrix([[R, -S], [S, R]], subdivide=False)
def rcomm(gens, sign=1):
    n = gens[0].nrows(); In = identity_matrix(QQ, n)
    M = block_matrix([[x.transpose().tensor_product(In) - sign * In.tensor_product(x)] for x in gens],
                     subdivide=False)
    return [matrix(QQ, n, n, list(v)).transpose() for v in M.right_kernel().basis()]
def norm_preserving(basis):          # antisymmetric part: preserves psi^dag psi
    if not basis: return []
    n = basis[0].nrows()
    M = matrix(QQ, [[(b + b.transpose()).list()[t] for b in basis] for t in range(n * n)])
    return [sum(v[k] * basis[k] for k in range(len(basis))) for v in M.right_kernel().basis()]
def span_rank(ms): return matrix(QQ, [m.list() for m in ms]).rank() if ms else 0
def derived(L):
    D = [a * b - b * a for a in L for b in L]
    if span_rank(D) == 0: return []
    n = L[0].nrows()
    return [matrix(QQ, n, n, list(v)) for v in matrix(QQ, [m.list() for m in D]).row_space().basis()]
def symmetry(GG): return norm_preserving(rcomm([realify(i_ * x) for x in GG]))

# ============================================================ [A]
print("\n[A] Geometry (Kaluza-Klein route)")
iso = lambda d: d * (d + 1) // 2          # max isometry-group dimension of a d-manifold
check(iso(1) == 1 and iso(1) < 3,
      "one extra dimension: isometry group of dim <= 1 < dim SU(2) = 3")
print("    (DT's mass axis is a half-line with an edge: not even the U(1) of a")
print("    circle.  DT's U(1) comes from the Weyl phase, not an isometry.)")

# ============================================================ [B]
print("\n[B] One 4x4 irrep (ledger Q3)")
L4 = symmetry(G4)
check(len(L4) == 1, "massive 5D Dirac, one irrep: internal symmetry u(1) only (the charge phase)")
extra = norm_preserving(rcomm([realify(i_ * x) for x in G4], sign=-1))
check(len(extra) == 2, "massless-only extras exist (Pauli-Guersey type), but they flip "
      "the kinetic term relative to the mass: DT's '+1' (Eq 110) removes them")

# ============================================================ [C]
print("\n[C] Williams' 8x8 (Eqs 109-111)")
L8 = symmetry(G8); D8 = derived(L8)
check(len(L8) == 4 and len(D8) == 3 and len(derived(D8)) == 3,
      "internal symmetry dim 4 = u(1) + su(2): NON-ABELIAN, from DT's own representation")
check(all(T + T.transpose() == 0 for T in L8),
      "generators antisymmetric (norm-preserving): compact, so su(2), not sl(2,R)")
Ii8 = realify(i_ * identity_matrix(K, 8))
check(any(T * Ii8 != Ii8 * T for T in D8),
      "su(2) does not commute with i: it mixes one block with the CONJUGATE of the other")
check(len(norm_preserving(rcomm([realify(i_ * x) for x in G8] + [Ii8]))) == 2,
      "complex-linear symmetry alone is U(1) x U(1) (strong_sector_su3 [A]); su(2) needs the conjugation")

# ============================================================ [D]
print("\n[D] Does it survive the switch?  (acts on the zero modes)")
I8 = identity_matrix(K, 8)
Gam = G8[0] * G8[1] * G8[2] * G8[3] * G8[4]
check(Gam * Gam == I8 and all(Gam * x == x * Gam for x in G8),
      "Gamma = g0g1g2g3g4 is central, square 1: its eigenspaces are the blocks A, B")
PA, PB = (I8 + Gam) / 2, (I8 - Gam) / 2
iG4 = i_ * G8[4]
Pp, Pm = (I8 + iG4) / 2, (I8 - iG4) / 2               # V+ = zero-mode space (M > 0)
check(all(T * realify(iG4) == realify(iG4) * T for T in L8),
      "the symmetry commutes with i*gamma^4: it preserves the zero-mode space V+")
Sig = [realify((G8[a] * G8[b] - G8[b] * G8[a]) / 4) for a in range(4) for b in range(a + 1, 4)]
check(len(rcomm(L8 + Sig)) == 4,
      "commutant of (symmetry + 4D Lorentz) = 4 = 2 + 2: exactly two invariant edge "
      "subspaces, V+ and V-")
a4 = G8[0] * G8[4]
check((Pp.conjugate_transpose() * a4 * Pp).is_zero() and (Pm.conjugate_transpose() * a4 * Pm).is_zero(),
      "both are self-adjoint (J^4 = 0) edge conditions: the choice is the hand (Q3 sign)")
g58 = i_ * G8[0] * G8[1] * G8[2] * G8[3]
hand = {}
for nm, P in (("A", PA), ("B", PB)):
    V = (Pp * P).column_space()
    hand[nm] = [c for c in (1, -1) if all(g58 * v == c * v for v in V.basis())]
    print("    block %s zero modes: dim %d, gamma^5 = %s" % (nm, V.dimension(), hand[nm]))
check(hand["A"] == [-1] and hand["B"] == [1],
      "block A zero mode LEFT, block B zero mode RIGHT (weak_chirality_u4 [C])")
RA, RB = realify(Pp * PA), realify(Pp * PB)
check(any(not (RB * T * RA).is_zero() for T in D8),
      "su(2) maps the A zero mode into (the conjugate of) the B zero mode")
print("    conjugate of a right-handed Weyl field is left-handed, so the doublet")
print("    (psi_A,L , (psi_B,R)^c) is TWO LEFT-HANDED fields: a chiral doublet.")

# ============================================================ [E]
print("\n[E] Where does electric charge sit?")
check(span_rank(D8 + [Ii8]) == span_rank(D8),
      "Williams' common phase (same e^{i alpha} on both blocks) lies INSIDE su(2): "
      "Q = T3, Y = 0, doublet charges (+q, -q) -- no neutral member")
QB = realify(i_ * PB)                                   # charge on block B only
rows, rhs = [], []
for Dj in D8:
    tgt = (QB * Dj - Dj * QB).list()
    cols = [(Dk * Dj - Dj * Dk).list() for Dk in D8]
    for t in range(len(tgt)):
        rows.append([c[t] for c in cols]); rhs.append(tgt[t])
sol = matrix(QQ, rows).solve_right(vector(QQ, rhs))
T3 = sum(sol[k] * D8[k] for k in range(3)); Y = QB - T3
check(all(Y * T == T * Y for T in L8) and span_rank(L8 + [Y]) == span_rank(L8),
      "charge 0 on A, 1 on B: Q = T3 + Y with T3 in su(2) and Y the native central u(1)")
check(Y == realify(i_ * (PB - PA) / 2),
      "Y = (1/2) x relative phase: Y(psi_A) = -1/2, Y(psi_B) = +1/2")
print("    doublet (psi_A, psi_B^c): Q = (0, -1), T3 = (+1/2, -1/2), Y = -1/2:")
print("    exactly (nu_L, e_L).  Given T3 = +-1/2, integer charge (q = Ze) with")
print("    one neutral member FORCES |Y| = 1/2 -- the SM lepton assignment.")

print("\n" + "=" * 74)
print("""VERDICT
  1. Geometry cannot supply SU(2): one extra dimension carries at most U(1)
     ([A]); a single 4x4 irrep carries only the U(1) of charge ([B]).
  2. Williams' OWN 8x8 representation (Eqs 109-111) carries a native
     su(2) + u(1) internal symmetry ([C]).  It is NOT imported: it exists
     because the two inequivalent blocks are related by charge conjugation.
  3. It SURVIVES THE SWITCH ([D]): it acts on the zero modes, the edge
     condition invariant under it (and 4D Lorentz) is unique up to the hand,
     and the doublet it acts on is two LEFT-HANDED fields.  So this is
     structure DT already has at u4 = 0, not a relaxation.
  4. Charge ([E]): with Q = T3 + Y and Y the native u(1), one neutral member
     plus integer charge gives (nu_L, e_L) with Y = -1/2.  Williams' common
     phase puts Q inside SU(2) (Y = 0, no neutral member): a WILLIAMS
     DISCREPANCY, the data-selected reading being the SM's.
  STANDING.  The GLOBAL SU(2)_L x U(1)_Y of one lepton doublet, its
  chirality (left-handed doublet) and |Y| = 1/2 are DERIVED from Williams'
  8x8 with DT's chiral edge.  The 8x8 is then one lepton doublet, not a
  doubled neutrino: the N_nu ~ 6 exclusion (weak_neutral_current [C])
  applied only to the reading where both blocks are neutral.
  Bracket: (i) GAUGING SU(2) x U(1) (W, Z) and sin^2 thW remain imported,
  just as the SM assumes them (DT derives only its abelian Weyl gauge);
  (ii) e_R needs a separate singlet field and its mass a Yukawa (Q13), as
  in the SM; (iii) the charge assignment (0 on A, 1 on B) is SELECTED, not
  derived -- the forcing is of |Y| given a neutral member; (iv) the symmetry
  is of the free field; interactions must respect it (as in the SM).
""")
print("=" * 74)
print("ALL ASSERTIONS PASSED")
print("=" * 74)
