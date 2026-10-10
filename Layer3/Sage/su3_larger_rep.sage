#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE (primary): LA-8370-MS (1980) Eqs 109-112: the 5D Dirac equation
#   (i gamma^j d_j + 1) psi = 0.  Any matrices obeying the 5D Clifford algebra
#   (Eq 112) are a DT representation; Williams' own choice is the 8x8 (4 + 4').
# Context: su2_from_5d (8x8 -> native su(2)+u(1)), strong_sector_su3 (N_c odd
#   derived, N_c = 3 selected; ledger Q14, Q15), noflux_boundary_derivation.
# ---------------------------------------------------------------------------
#
# su3_larger_rep.sage                                  [A] (assertion-hardened)
#
# PURPOSE.  su2_from_5d found SU(2) native in Williams' 8x8.  Does SU(3) sit
# in a larger DT representation the same way?
#
#   [A] every DT field is a copies of 4 (kappa = +1) plus b copies of 4'
#       (kappa = -1); its internal symmetry (psi treated as real) is EXACTLY
#       u(N), N = a + b -- computed for N = 1..4
#   [B] every such symmetry survives u4 = 0: it acts on N LEFT-handed zero
#       modes (4' zero modes after conjugation)
#   [C] N = 3: su(3) present; with all blocks alike the charge phase commutes
#       with it (neutral gluons); the triplet carries equal charges
#   [D] what fixes N?  Clifford theory allows every N; Williams' factorisation
#       stops at N = 2.  So SU(3) is AVAILABLE, not singled out
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("su3_larger_rep.sage -- SU(3) in a larger DT representation?")
print("=" * 74)

K = CyclotomicField(4); i_ = K.gen()
I2, Z2 = identity_matrix(K, 2), zero_matrix(K, 2)
s = [matrix(K, [[0, 1], [1, 0]]), matrix(K, [[0, -i_], [i_, 0]]),
     matrix(K, [[1, 0], [0, -1]])]
B = lambda a, b, c, d: block_matrix([[a, b], [c, d]], subdivide=False)
g = [B(I2, Z2, Z2, -I2)] + [B(Z2, si, -si, Z2) for si in s]
g5 = i_ * g[0] * g[1] * g[2] * g[3]
irrep = {+1: g + [i_ * g5], -1: g + [-i_ * g5]}       # 4 and 4'
eta = [1, -1, -1, -1, -1]

def field(kappas):
    G = [irrep[kappas[0]][j] for j in range(5)]
    for k in kappas[1:]:
        G = [G[j].block_sum(irrep[k][j]) for j in range(5)]
    return G
def realify(A):
    n = A.nrows()
    R = matrix(QQ, n, n, lambda a, b: A[a, b].list()[0])
    S = matrix(QQ, n, n, lambda a, b: A[a, b].list()[1])
    return block_matrix([[R, -S], [S, R]], subdivide=False)
def rcomm(gens):
    n = gens[0].nrows(); In = identity_matrix(QQ, n)
    M = block_matrix([[x.transpose().tensor_product(In) - In.tensor_product(x)] for x in gens],
                     subdivide=False)
    return [matrix(QQ, n, n, list(v)).transpose() for v in M.right_kernel().basis()]
def norm_preserving(basis):
    n = basis[0].nrows()
    M = matrix(QQ, [[(b + b.transpose()).list()[t] for b in basis] for t in range(n * n)])
    return [sum(v[k] * basis[k] for k in range(len(basis))) for v in M.right_kernel().basis()]
def span_rank(ms): return matrix(QQ, [m.list() for m in ms]).rank() if ms else 0
def derived(L):
    D = [a * b - b * a for a in L for b in L]
    if span_rank(D) == 0: return []
    n = L[0].nrows()
    return [matrix(QQ, n, n, list(v)) for v in matrix(QQ, [m.list() for m in D]).row_space().basis()]

# ============================================================ [A]
print("\n[A] Internal symmetry of every DT field  a*4 + b*4'")
res = {}
for kap in [(1,), (-1,), (1, 1), (1, -1), (1, 1, 1), (1, 1, -1), (1, -1, -1), (1, 1, -1, -1)]:
    G = field(list(kap)); n = G[0].nrows()
    assert all(G[j] * G[k] + G[k] * G[j] == (2 * eta[j] if j == k else 0) * identity_matrix(K, n)
               for j in range(5) for k in range(5))
    L = norm_preserving(rcomm([realify(i_ * x) for x in G])); D = derived(L)
    N = len(kap); res[kap] = (G, L, D)
    print("    a=%d b=%d (N=%d, %2dx%-2d): symmetry dim %2d, derived dim %2d  -> u(%d)"
          % (kap.count(1), kap.count(-1), N, n, n, len(L), len(D), N))
check(all(len(L) == len(k)**2 and len(D) == len(k)**2 - 1 for k, (G, L, D) in res.items()),
      "symmetry = u(N), N = a + b, for every mix: 4' counts as a conjugate 4 "
      "(the 8x8 su(2) of su2_from_5d is the N = 2 case)")

# ============================================================ [B]
print("\n[B] Survives the switch: N left-handed zero modes")
ok = True
for kap, (G, L, D) in res.items():
    n = G[0].nrows()
    iG4 = realify(i_ * G[4])
    ok &= all(T * iG4 == iG4 * T for T in L)
    V = ((identity_matrix(K, n) + i_ * G[4]) / 2).column_space()
    ok &= len(V.basis()) == 2 * len(kap)
check(ok, "every u(N) commutes with i*gamma^4 (preserves the zero-mode space); "
      "one 2-component zero mode per block (handedness per block: su2_from_5d [D])")

# ============================================================ [C]
print("\n[C] N = 3: colour")
G3, L3, D3 = res[(1, 1, 1)]
check(len(D3) == 8 and len(derived(D3)) == 8, "su(3) (dim 8, simple part of u(3)) acts on the triplet")
Ii = realify(i_ * identity_matrix(K, 12))
check(all(T * Ii == Ii * T for T in L3),
      "a = 3, b = 0: the charge phase commutes with su(3) -- gluons neutral, triplet "
      "members carry EQUAL charge")
G21, L21, D21 = res[(1, 1, -1)]
check(any(T * realify(i_ * identity_matrix(K, 12)) != realify(i_ * identity_matrix(K, 12)) * T for T in D21),
      "a = 2, b = 1: same su(3), but it mixes particle and antiparticle -- equal "
      "charges again require the conjugate assignment on the 4' block (as for Y in su2_from_5d)")
print("    The triplet charge itself (fractional, Q_u = (N+1)/2N, Q_d = -(N-1)/2N)")
print("    is fixed by the anomaly conditions, strong_sector_su3 [B].")

# ============================================================ [D]
print("\n[D] What fixes N?")
print("    Eq 112 is satisfied by a*4 + b*4' for EVERY a, b: Clifford theory does")
print("    not fix N.  Williams' factorisation (Eqs 106-111) yields the 8x8, N = 2,")
print("    which is why SU(2) came out singled out: U(2) is all there is.")
check(all(len(res[k][1]) == len(k)**2 for k in res),
      "the symmetry depends ONLY on N: nothing in the 5D structure prefers N = 3")
print("    For N = 3 DT supplies: odd N (Witten, given doublets) and fractional")
print("    charge with confinement as a selection rule -- DERIVED (strong_sector_su3).")
print("    N = 3 itself: SELECTED by data (nucleon charges, R = 2.21, R_l = 20.82).")

print("\n" + "=" * 74)
print("""VERDICT
  1. The internal symmetry of ANY DT Dirac field is exactly u(N), N = number
     of 4x4 blocks ([A]), and it always survives u4 = 0, acting on N
     left-handed zero modes ([B]).  So the FORM of colour is native: any DT
     field with three blocks carries su(3), with neutral gluons and equal
     charges across the triplet ([C]).
  2. But N is not fixed by the 5D structure ([D]).  SU(2) was singled out
     because Williams' own factorisation gives N = 2; nothing in DT gives a
     12x12.  SU(3) is therefore AVAILABLE in DT, not SINGLED OUT: N = 3 is
     selected by data, with N odd and fractional charge derived.
  STANDING (ledger Q14).  Colour: the group FORM is DERIVED (u(N) on N zero
  modes); the multiplicity N = 3 is SELECTED (odd N derived); the gauging
  and the SU(3) dynamics remain SHARED-SM by import (Q15), just as the SM
  assumes them.  One level weaker than SU(2), whose N came from Williams.
  Bracket: (i) u(N) is the generic symmetry of N identical free Weyl
  fields -- DT's contribution is that its 5D fields deliver exactly this
  and nothing larger, and that it survives the switch; (ii) a quark
  doublet would be N = 6, u(6) containing su(3) x su(2) x u(1) as the
  (3,2) subgroup -- choosing that subgroup is an import; (iii) free-field
  symmetry only; interactions must respect the gauged subgroup.
""")
print("=" * 74)
print("ALL ASSERTIONS PASSED")
print("=" * 74)
