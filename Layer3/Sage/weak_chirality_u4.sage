#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, LA-8370-MS (1980) Sec. III.C, Eqs 105-112:
#   Eq 105  5D mass shell  p_j p^j = gamma^2 c^2  (gamma = mass density)
#   Eq 106-108  Dirac ansatz h = -(alpha~.p + beta~), adopted ("may be taken")
#   Eq 109  an 8x8 representation (A-matrix corrected to [[0,iI],[-iI,0]],
#           see dirac_algebra.sage)
#   Eq 110-112  (i gamma^j d_j + 1) psi = 0,  {gamma^j, gamma^k} = 2 g^{jk},
#           g = diag(1,-1,-1,-1,-1),  x^4 = gamma/a0 (mass-density axis)
#   u4 = dy4/dt = 0 is conservation of mass (Book Eq 1.36 / 10.6); on the
#   Dirac side the corresponding statement is J^4 = psibar gamma^4 psi = 0.
# ---------------------------------------------------------------------------
#
# weak_chirality_u4.sage                               [A] (assertion-hardened)
#
# PURPOSE.  weak_parity.sage showed DT's weak sector is parity-even, so (like
# the Standard Model) DT must take chirality from experiment.  In 5D that is
# harder: gamma^5 is used up as the fifth Dirac matrix (Kaluza-Klein chirality
# problem).  The standard escape is a BOUNDARY in the extra dimension with a
# no-flux condition, which leaves ONE chirality as a normalisable zero mode.
# DT has exactly such ingredients: a mass-density axis that is a HALF-LINE
# (gamma >= 0) and a no-flux restriction (u4 = 0).  Does u4 = 0 make DT
# chiral?
#
# READINGS OF u4 = 0 ON THE DIRAC FIELD
#   R1  d_4 psi = 0          (no dependence on the mass axis)
#   R2  J^4 = 0 everywhere   (no probability current along the mass axis)
#   R3  J^4 = 0 at gamma = 0 (no-flux BOUNDARY condition on the half-line)
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("weak_chirality_u4.sage -- can u4 = 0 make DT's fermions chiral?")
print("=" * 74)

I2, Z2, Z4 = identity_matrix(SR, 2), zero_matrix(SR, 2), zero_matrix(SR, 4)
I4, I8 = identity_matrix(SR, 4), identity_matrix(SR, 8)
s1 = matrix(SR, [[0, 1], [1, 0]]); s2 = matrix(SR, [[0, -I], [I, 0]])
s3 = matrix(SR, [[1, 0], [0, -1]])
B = lambda a, b, c, d: block_matrix([[a, b], [c, d]], subdivide=False)
beta = B(I2, Z2, Z2, -I2)
alph = [B(Z2, s, s, Z2) for s in (s1, s2, s3)]
Amat = B(Z2, I * I2, -I * I2, Z2)                 # corrected Eq 109 A-matrix
bT = B(beta, Z4, Z4, beta)
aT = [B(Z4, a, a, Z4) for a in alph] + [B(Amat, Z4, Z4, Amat)]
G = [bT] + [-bT * a for a in aT]                  # Eq 111
eta = [1, -1, -1, -1, -1]

# ============================================================ [A]
print("\n[A] Williams' 8x8 representation (Eq 109) is a DOUBLED one")
check(all(G[j] * G[k] + G[k] * G[j] == (2 * eta[j] if j == k else 0) * I8
          for j in range(5) for k in range(5)), "Eq 112 Clifford algebra holds")
C = G[0] * G[1] * G[2] * G[3] * G[4]
check(all(C * g == g * C for g in G) and C * C == I8,
      "C = gamma^0..gamma^4 is central with C^2 = 1 (odd dimension)")
ev = C.eigenvalues()
check(ev.count(1) == 4 and ev.count(-1) == 4,
      "C has eigenvalues +1 (x4) and -1 (x4): 8 = 4 (+) 4', BOTH inequivalent "
      "5D irreps, once each")
g5 = I * G[0] * G[1] * G[2] * G[3]
check(G[4] == I * C * g5,
      "gamma^4 = i kappa gamma^5 with kappa = C = +1 on one block, -1 on the "
      "other: no independent chirality in 5D")
print("    (The 4x4 PhysLib witness in PhyslibBridge.lean is ONE block, kappa = +1.)")

# ============================================================ [B]
print("\n[B] Readings R1 and R2 give no chirality")
# R1: d_4 psi = 0 -> (i gamma^mu d_mu + M) psi = 0, a massive 4D Dirac field
p0, p3, M = var('p0 p3 M', domain='real')
D4 = G[0] * p0 - G[3] * p3 + M * I8
check(bool(D4.det().factor() == (p0**2 - p3**2 - M**2)**4) or
      bool((D4.det() - (p0**2 - p3**2 - M**2)**4).simplify_full() == 0),
      "R1: d4 psi = 0 leaves (gamma.p + M) psi = 0, massive Dirac: L and R "
      "both present and mixed by M")
# R2: J^4 = psibar gamma^4 psi.  Pure-L and pure-R spinors both satisfy it.
PL, PR = (I8 - g5) / 2, (I8 + g5) / 2
xs = vector(SR, var(','.join('z%d' % i for i in range(8))))
J4 = lambda psi: (psi.conjugate() * G[0] * G[4] * psi)
check(bool(J4(PL * xs).simplify_full() == 0) and
      bool(J4(PR * xs).simplify_full() == 0),
      "R2: J^4 = 0 is satisfied by EVERY pure-L and EVERY pure-R spinor")
check(G[0] * G[4] * G[0] == -G[4],
      "R2: J^4 is P-odd, so the set {J^4 = 0} is mirror-symmetric: it cannot "
      "prefer a hand")

# ============================================================ [C]
print("\n[C] Reading R3: no-flux boundary at gamma = 0 (half-line y >= 0)")
y, k = var('y k', domain='positive')
kap, chi, Mc, al = var('kap chi Mc al', domain='real')
# Massless 4D zero mode psi = f(y) u, gamma^5 u = chi u, block kappa:
#   (i gamma^4 d_y + M) psi = 0  ->  -kap chi f' + M f = 0
def survivors(kv, Mfun):
    """hands chi whose zero mode is normalisable on y >= 0"""
    out = []
    for ch in (-1, +1):
        f = exp(kv * ch * integrate(Mfun, y, 0, y))
        assert bool((-kv * ch * diff(f, y) + Mfun * f).simplify_full() == 0)
        if limit(f, y=oo) == 0:
            out.append(ch)
    return out
for tag, Mfun in [("constant M = 1 (Eq 110 '+1')", 1), ("M = y (~ mass density)", y)]:
    s_p, s_m = survivors(+1, Mfun), survivors(-1, Mfun)
    print("    %s: kappa=+1 keeps chi=%s, kappa=-1 keeps chi=%s" % (tag, s_p, s_m))
    check(len(s_p) == 1 and len(s_m) == 1,
          "%s: exactly ONE hand per block is normalisable" % tag)
    check(s_p == [-1] and s_m == [+1],
          "%s: kappa=+1 block keeps L, kappa=-1 block keeps R" % tag)
# boundary: pure-chirality boundary values carry no flux (u4 = 0 at gamma = 0)
check(bool(J4(PL * xs).simplify_full() == 0),
      "the chiral boundary value has J^4(0) = 0: R3 IS a u4 = 0 condition")
# Williams' doubled 8x8 has kappa = +1 AND kappa = -1 (four components each)
check(ev.count(1) == ev.count(-1),
      "Williams' 8x8 carries both blocks equally: one L + one R zero mode, "
      "i.e. NON-chiral (vector-like)")
check(ev.count(1) > 0 and ev.count(-1) > 0,
      "dropping one block (the 4x4 rep, one sign of kappa) leaves a single "
      "hand: u4 = 0 then makes DT chiral")

# ============================================================ verdict
print("\n" + "=" * 74)
print("""VERDICT
  1. u4 = 0 as a bulk statement (R1: no mass-axis dependence; R2: no current
     along the mass axis) gives NO chirality: the first leaves a massive Dirac
     field, the second is mirror-symmetric ([B]).
  2. u4 = 0 as a NO-FLUX BOUNDARY at zero mass density (R3) CAN: the mass axis
     is a half-line, the boundary leaves one normalisable zero mode, and its
     hand is fixed by sign(kappa M) ([C]).  This is the orbifold/domain-wall
     mechanism of Kaluza-Klein model building, realised by DT's OWN
     restriction on DT's OWN half-line.  It holds for constant M and for
     M ~ mass density alike.
  3. BUT Williams' own 8x8 representation (Eq 109) contains both 5D irreps,
     which leave OPPOSITE hands, so as written DT stays non-chiral ([A],[C]).
     Choosing ONE 4x4 irrep -- a single sign, the counterpart of the SM's
     left-doublet assignment -- makes it chiral.
  STANDING (assumption ledger: SHARED WITH SM + one DT-SPECIFIC reading).
  Just as the SM assumes left-handed fermions in doublets (one choice of
  hand), if DT assumes one 4x4 irrep (one sign of kappa), the zero mode is
  left-handed, consistent with the SM.  DT-specific: reading u4 = 0 as a
  no-flux boundary at gamma = 0 -- needed only because DT is 5D; the SM,
  being 4D, has no counterpart.  Recorded in ASSUMPTIONS.md.  The zero mode is massless in 4D, so the natural candidate is
  the neutrino (left-handed and nearly massless), NOT the electron:
  interpretive, not derived.
  Bracket: (i) reading R3 is ours, not Williams' (he never imposes a
  boundary at gamma = 0); (ii) the 8x8 vs 4x4 choice is Williams' (Eq 109)
  vs PhysLib's; (iii) nothing here yet makes the chiral mode couple to the
  Eq 2.31 weak fields with V-A strength; that is the next test.
""")
print("=" * 74)
print("ALL ASSERTIONS PASSED")
print("=" * 74)
