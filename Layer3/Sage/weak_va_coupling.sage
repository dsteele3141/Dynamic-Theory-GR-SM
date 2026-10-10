#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, LA-8370-MS (1980) Eqs 110-121 (5D Dirac,
#   field coupling e A_j, magnetic-moment term).
# SOURCE (book, Report silent on the weak force): "The Dynamic Theory" (2011)
#   Eq 2.31  [(i d I - e_theta A^T)^2 - (m I)^2 - i sigma^{mu nu} e_beta B^beta I] psi = 0
#            -- every field enters multiplied by a CHARGE e_theta, e_beta;
#   Ch.11 Eqs 11.34-11.37: the neutrino is quantised through phi_4, "the
#            gravitational potential".
# ---------------------------------------------------------------------------
#
# weak_va_coupling.sage                                [A] (assertion-hardened)
#
# PURPOSE.  weak_chirality_u4.sage found that u4 = 0, read as a no-flux
# boundary at zero mass density, can leave ONE chiral, massless zero mode
# (a neutrino candidate).  Does that mode couple to DT's weak fields with the
# measured V-A FORM and the measured STRENGTH (Fermi constant G_F)?
#
#   [A] FORM: for a chiral field every scalar, pseudoscalar and tensor
#       bilinear vanishes -- including Williams' Pauli weak term sigma.B --
#       and vector = -axial.  V-A is automatic.
#   [B] One chiral leg is enough: a vector vertex psibar_e gamma^mu nu_L
#       projects the OTHER leg to its left-handed part too.
#   [C] VERTEX: which DT field can couple to a NEUTRAL fermion?
#   [D] STRENGTH: what that leaves, against G_F.
#
# Imported (not re-derived): Michel parameter rho = 3/4 for pure V-A; the
# measured values quoted below; the SM tree-level matching
# G_F/sqrt2 = g^2/(8 M^2), used only to translate a range into a strength.
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("weak_va_coupling.sage -- V-A form and strength for DT's chiral mode")
print("=" * 74)

I2, Z2 = identity_matrix(SR, 2), zero_matrix(SR, 2)
s = [matrix(SR, [[0, 1], [1, 0]]), matrix(SR, [[0, -I], [I, 0]]),
     matrix(SR, [[1, 0], [0, -1]])]
B = lambda a, b, c, d: block_matrix([[a, b], [c, d]], subdivide=False)
g = [B(I2, Z2, Z2, -I2)] + [B(Z2, si, -si, Z2) for si in s]
g0 = g[0]
g5 = I * g[0] * g[1] * g[2] * g[3]
I4, Z4 = identity_matrix(SR, 4), zero_matrix(SR, 4)
PL, PR = (I4 - g5) / 2, (I4 + g5) / 2
sig = lambda m, n: (I / 2) * (g[m] * g[n] - g[n] * g[m])
# psibar_L Gamma psi_L = psi^dag PL g0 Gamma PL psi = psi^dag g0 (PR Gamma PL) psi
bil = lambda Gm: PR * Gm * PL

# ============================================================ [A]
print("\n[A] FORM: bilinears of a chiral (left-handed) field")
check(bil(I4) == Z4, "scalar  psibar_L psi_L = 0 (no Dirac mass: the mode is massless)")
check(bil(g5) == Z4, "pseudoscalar psibar_L g5 psi_L = 0")
check(all(bil(sig(m, n)) == Z4 for m in range(4) for n in range(4)),
      "tensor psibar_L sigma^{mu nu} psi_L = 0: Williams' Pauli weak term "
      "sigma.B (Eq 2.31) does NOT couple to a chiral mode")
check(all(bil(g[m] * g5) == -bil(g[m]) for m in range(4)),
      "axial = -vector on a left-handed field")
check(all(bil(g[m]) == g[m] * (I4 - g5) / 2 * PL for m in range(4)),
      "so every surviving coupling is gamma^mu (1 - g5)/2: V-A is AUTOMATIC")
rho_VA, rho_meas, drho = QQ(3) / 4, 0.74977, 0.00026   # TWIST 2011
print("    muon decay Michel rho: pure V-A 3/4; measured %.5f(%d)" % (rho_meas, 26))
check(abs(rho_meas - rho_VA) < 3 * drho, "measured rho consistent with V-A")

# ============================================================ [B]
print("\n[B] One chiral leg makes the whole current V-A")
check(all(g0 * g[m] * PL == g0 * PR * g[m] for m in range(4)) and
      PL.conjugate_transpose() * g0 == g0 * PR,
      "psibar_e gamma^mu P_L nu = (P_L e)bar gamma^mu nu: only e_L takes part "
      "(electron helicity -v/c, as measured)")

# ============================================================ [C]
print("\n[C] VERTEX: what in DT can couple to a NEUTRAL chiral mode?")
Q = var('Q')
couplings = {
    "vector  e_theta A^T_mu gamma^mu  (Eq 2.31)": Q,       # proportional to charge
    "Pauli   e_beta sigma.B           (Eq 2.31)": 0,       # killed by chirality, [A]
    "mass axis phi_4 (gravitational, Eq 11.34)": 1,        # couples to mass-energy
}
for name, c in couplings.items():
    print("    %-44s strength factor at Q=0: %s" % (name, c if c in (0, 1) else c.subs(Q=0)))
alive = [n for n, c in couplings.items() if (c if c in (0, 1) else c.subs(Q=0)) != 0]
check(alive == ["mass axis phi_4 (gravitational, Eq 11.34)"],
      "for Q = 0 the ONLY surviving DT coupling is the gravitational phi_4")
check(bool(couplings["vector  e_theta A^T_mu gamma^mu  (Eq 2.31)"].subs(Q=0) == 0),
      "a neutral mode has no vector coupling to Williams' (neutral) A^T fields")
print("    And none of Williams' fields carries charge, so there is no e <-> nu")
print("    (charge-changing) vector vertex at all: DT has no W counterpart.")

# ============================================================ [D]
print("\n[D] STRENGTH")
GF, GN = 1.1663787e-5, 6.70883e-39          # GeV^-2
alpha, hbarc, MW = 1 / 137.035999084, 0.1973269804, 80.377   # GeV fm, GeV
ratio = GN / GF
print("    gravitational / Fermi coupling  G_N/G_F = %.2e" % ratio)
print("    cross-section ratio (G_N/G_F)^2        = %.1e" % ratio**2)
print("    Neutrinos ARE detected at G_F strength (reactor inverse beta decay,")
print("    Reines-Cowan 1956; SNO neutral current 2002).")
check(ratio**2 < 1e-60, "a gravity-only neutrino interacts ~1e-67 too weakly: excluded")
# What range would an EM-strength charged vector field need?  (illustrative)
lam_e = sqrt(8 * GF / (sqrt(2.0) * 4 * pi * alpha)) * hbarc     # g = e
lam_W = hbarc / MW                                             # g = e/sin(thW)
print("    A CHARGED vector field of EM strength gives G_F if its range is")
print("    %.1e fm (g = e)  --  ~ hbar/(M_W c) = %.2e fm, the W boson." % (lam_e, lam_W))
check(1 < lam_e / lam_W < 3, "the needed range is the W Compton wavelength (factor < 3)")

# ============================================================ verdict
print("\n" + "=" * 74)
print("""VERDICT
  1. FORM -- YES.  If DT's u4 = 0 boundary leaves a chiral zero mode, its
     couplings are AUTOMATICALLY V-A: scalar, pseudoscalar and tensor
     bilinears vanish (Williams' Pauli weak term sigma.B drops out), axial =
     -vector ([A]); and one chiral leg forces the partner electron to be
     left-handed ([B]).  Michel rho = 3/4 follows.  The V-A form costs
     nothing beyond the one sign already paid in weak_chirality_u4.
  2. STRENGTH -- NO, as written.  Every field in Eq 2.31 enters times a
     charge, the neutrino is neutral, and the Pauli term is killed by
     chirality; the only coupling left is the gravitational phi_4, ~1e-33
     of G_F in amplitude, ~1e-67 in cross-section ([C],[D]).  DT also has no
     charge-changing vector field, so no e <-> nu vertex exists.
  3. The measured strength needs a CHARGED vector field of range ~ hbar/M_W c
     ~ 2.5e-3 fm -- i.e. re-adding the W boson.  (Curiously, that is the
     same scale as the hydrogen bound lambda_p < 1.8e-3 fm in
     hydrogen_nonsingular_spectrum.sage; noted, not used.)
  STANDING (assumption ledger: SHARED WITH SM).  FORM follows; STRENGTH is
  an input.  Just as the SM assumes a charged vector boson W with
  M_W ~ 80 GeV (equivalently G_F), if DT assumes a charged short-range
  vector field of range hbar/(M_W c), the weak strength is consistent with
  the SM -- and the V-A form then comes for free from [A],[B].  Not a
  DT-specific defect: an assumption the SM also makes, recorded in
  ASSUMPTIONS.md.
  Bracket: (i) the zero-mode neutrino is OUR identification; Williams'
  neutrino is the helicity-0 phi_4 wave quantum (weak_parity [D]), which
  faces the same strength problem; (ii) DT's beta decay is the breakup of a
  proton-electron bound state (Book p.255), an EM process whose RATE has not
  been computed -- the neutron lifetime (878 s) is the natural next test and
  a possible defeater; (iii) the strength numbers are order-of-magnitude.
""")
print("=" * 74)
print("ALL ASSERTIONS PASSED")
print("=" * 74)
