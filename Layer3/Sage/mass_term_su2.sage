#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE (primary): LA-8370-MS (1980) Eq 105 (5D mass shell p_j p^j =
#   gamma^2 c^2), Eqs 110-112 ((i gamma^j d_j + 1) psi = 0, x^4 = gamma/a0):
#   in DT a fermion's mass comes from the 5D operator -- the bulk term "+1"
#   and the mass-axis derivative gamma^4 d_4.
# Context: weak_chirality_u4 (zero modes, kappa), weak_neutral_current
#   (ledger Q6 SU(2)_L x U(1)_Y couplings, Q12 doublet Higgs).
# ---------------------------------------------------------------------------
#
# mass_term_su2.sage                                   [A] (assertion-hardened)
#
# PURPOSE.  Once DT imports SU(2)_L (Q6), is DT's own 5D mass term compatible
# with it, or must the charged-fermion masses come from a Higgs Yukawa?
#
#   [A] both 5D mass sources (bulk M psibar psi, mass-axis psibar gamma^4 d4 psi)
#       FLIP chirality: they link psi_L to psi_R
#   [B] a 5D field carries ONE SU(2) rep for both chiralities, so DT's mass is
#       SU(2)-invariant only for a VECTOR-LIKE electron (L and R both doublets)
#   [C] that option against data: Z-pole electron couplings
#   [D] the chiral option (L doublet, R singlet = two 5D fields): DT's mass
#       gives the zero modes NOTHING; only a doublet scalar (Yukawa) works
#
# Measured (imported): g_A^e = -0.503(17) (CHARM II); A_e = 0.1515(19)
# (LEP/SLD); beta-decay electron helicity -v/c (Frauenfelder 1957).
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("mass_term_su2.sage -- DT's 5D mass term vs imported SU(2)_L")
print("=" * 74)

I2, Z2 = identity_matrix(SR, 2), zero_matrix(SR, 2)
s = [matrix(SR, [[0, 1], [1, 0]]), matrix(SR, [[0, -I], [I, 0]]),
     matrix(SR, [[1, 0], [0, -1]])]
Bk = lambda a, b, c, d: block_matrix([[a, b], [c, d]], subdivide=False)
g = [Bk(I2, Z2, Z2, -I2)] + [Bk(Z2, si, -si, Z2) for si in s]
g5 = I * g[0] * g[1] * g[2] * g[3]
I4 = identity_matrix(SR, 4)
PL, PR = (I4 - g5) / 2, (I4 + g5) / 2
g4 = I * g5                                          # kappa = +1 block (weak_chirality_u4 [A])

# ============================================================ [A]
print("\n[A] Both 5D mass sources flip chirality")
# psibar Gamma psi = psi^dag g0 Gamma psi ; its L-L and R-R parts are
# psi^dag P_L^dag g0 Gamma P_L psi etc.  Show they vanish for Gamma = 1, g4.
for name, Gm in [("bulk mass  psibar psi", I4), ("mass axis  psibar g^4 d4 psi", g4)]:
    LL = PL.conjugate_transpose() * g[0] * Gm * PL
    RR = PR.conjugate_transpose() * g[0] * Gm * PR
    LR = PL.conjugate_transpose() * g[0] * Gm * PR
    check(LL.is_zero() and RR.is_zero() and not LR.is_zero(),
          "%s = psibar_L (.) psi_R + h.c. only (no L-L, no R-R part)" % name)

# ============================================================ [B]
print("\n[B] SU(2) invariance of an L-R bilinear")
T = [si / 2 for si in s]
def invariant_couplings(dimL, repL, dimR, repR):
    """dimension of the space of C with psibar_L C psi_R invariant:
       repL[a]^dag C - C repR[a] = 0 for all generators a  (Hermitian gens)"""
    cs = var(','.join('c%d' % i for i in range(dimL * dimR)))
    C = matrix(SR, dimL, dimR, cs)
    eqs = [e for a in range(3) for e in (repL[a] * C - C * repR[a]).list()]
    M = matrix(SR, [[e.coefficient(c) for c in cs] for e in eqs]) if eqs else matrix(SR, 0, len(cs))
    return len(cs) - M.rank()
singlet = [zero_matrix(SR, 1)] * 3
n_vec = invariant_couplings(2, T, 2, T)
n_chi = invariant_couplings(2, T, 1, singlet)
print("    L doublet, R doublet (vector-like): %d invariant mass coupling(s)" % n_vec)
print("    L doublet, R singlet (chiral, SM) : %d invariant mass coupling(s)" % n_chi)
check(n_vec == 1, "vector-like: DT's 5D mass IS SU(2)-invariant (one coupling, the identity)")
check(n_chi == 0, "chiral: NO SU(2)-invariant mass term exists for L doublet x R singlet")
print("    A 5D field has ONE SU(2) rep for both chiralities: DT's native mass is")
print("    therefore the vector-like case.")

# ============================================================ [C]
print("\n[C] The vector-like electron against data")
s2 = 0.2315
gL_e, gR_chi, gR_vec = -0.5 + s2, s2, -0.5 + s2
gA = lambda gl, gr: gl - gr
Ae = lambda gl, gr: (gl**2 - gr**2) / (gl**2 + gr**2)
gAm, gAs, Aem, Aes = -0.503, 0.017, 0.1515, 0.0019
for tag, gr in [("chiral (SM)", gR_chi), ("vector-like (DT mass)", gR_vec)]:
    print("    %-22s g_A^e = %+.3f (pull %+5.1f)   A_e = %.4f (pull %+6.1f)"
          % (tag, gA(gL_e, gr), (gA(gL_e, gr) - gAm) / gAs, Ae(gL_e, gr), (Ae(gL_e, gr) - Aem) / Aes))
check(abs(gA(gL_e, gR_chi) - gAm) < 2 * gAs and abs(Ae(gL_e, gR_chi) - Aem) < 3 * Aes,
      "chiral assignment fits g_A^e (<2 sigma) and A_e (<3 sigma; the 2.3 sigma is the "
      "SM's own SLD-vs-LEP sin^2 tension, not DT's)")
check(abs(gA(gL_e, gR_vec) - gAm) / gAs > 25 and abs(Ae(gL_e, gR_vec) - Aem) / Aes > 75,
      "vector-like electron EXCLUDED: g_A^e = 0 (~30 sigma), A_e = 0 (~80 sigma)")
print("    (and a W coupling to e_R would make beta electrons helicity 0, not -v/c)")

# ============================================================ [D]
print("\n[D] Chiral option: two 5D fields, zero modes of opposite hand")
y = var('y', domain='positive'); M = var('M', domain='positive')
fD, fS = exp(-M * y), exp(-M * y)          # L zero mode of doublet, R zero mode of singlet
uL = PL * vector(SR, [1, 0, 0, 0]) * 2; uR = PR * vector(SR, [1, 0, 0, 0]) * 2
for name, Gm in [("bulk mass", I4), ("mass axis", g4)]:
    within = [(u.conjugate() * g[0] * Gm * u).simplify_full() for u in (uL, uR)]
    check(all(w == 0 for w in within), "%s gives each zero mode ZERO mass (within one field)" % name)
cross = (uL.conjugate() * g[0] * uR)
check(cross != 0 and n_chi == 0,
      "the only L-R link is ACROSS the two fields, and it is SU(2)-forbidden ([B])")
# Yukawa psibar_L Phi psi_R: Phi psi_R transforms as a doublet, so it is the
# vector-like count again with Phi supplying the doublet index
n_yuk = invariant_couplings(2, T, 2, T)
check(n_yuk == 1, "psibar_L Phi psi_R with a DOUBLET scalar Phi: one invariant (the Yukawa)")
overlap = integrate(fD * fS, y, 0, oo)
print("    Yukawa strength in 4D ~ overlap of the two zero-mode profiles = %s" % overlap)
check(bool(overlap > 0), "with the Higgs doublet the zero modes get m = y_e v x overlap")

print("\n" + "=" * 74)
print("""VERDICT
  1. DT's 5D mass (bulk "+1" and the mass-axis derivative) links psi_L to
     psi_R ([A]).  Inside one 5D field both carry the same SU(2) charge, so
     the mass is SU(2)-invariant ONLY for a vector-like electron ([B]).
  2. The vector-like electron is EXCLUDED: g_A^e = 0 vs -0.503(17) (~30 sigma),
     A_e = 0 vs 0.1515(19) (~80 sigma) ([C]).
  3. The chiral electron (L doublet, R singlet) must be two 5D fields with
     zero modes of opposite hand.  DT's 5D mass then gives each zero mode
     ZERO mass, and the only L-R link is SU(2)-forbidden: the electron's mass
     must come from a Yukawa coupling to the Higgs doublet ([D]).
  So the tension flagged in Q12 resolves AGAINST DT's native mass for
  elementary charged fermions.  DT's 5D mass survives only as the scale of
  the heavy (vector-like) mass-axis modes and the zero-mode profiles.
  STANDING (assumption ledger: SHARED WITH SM, by import).  Just as the SM
  assumes Yukawa couplings of fermions to the Higgs doublet, if DT assumes
  the same, the charged-lepton masses and couplings are consistent with the
  SM.  Neutral, but plain: DT's own mass mechanism (mass as the fifth
  coordinate, Report Eq 105) does NOT give the electron its mass once
  SU(2)_L is imported.  Ledger Q13.
  Bracket: (i) this follows from importing SU(2)_L (Q6), which the neutral
  current forced; DT without SU(2) keeps its native mass but loses the
  weak sector; (ii) the electron mass value (y_e) is an input, as in the
  SM; (iii) the hydrogen tests used a massive non-chiral electron (reading
  R1) -- unaffected in QED, where e_L and e_R carry the same charge.
""")
print("=" * 74)
print("ALL ASSERTIONS PASSED")
print("=" * 74)
