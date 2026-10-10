#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE: Report (LA-8370-MS) Eqs 105-112 (5D Dirac); book (2011) Eq 2.31
#   (weak/EM fields, each entering times a charge) and Ch.11 Eq 11.34 (phi_4).
#   Williams never discusses a neutral weak current (no match for "neutral
#   current", "Z boson", "Weinberg" or "electroweak" in the book).
# ---------------------------------------------------------------------------
#
# weak_neutral_current.sage                            [A] (assertion-hardened)
#
# PURPOSE.  Ledger Q6.  The neutral current (Gargamelle 1973) is measured in
# four independent ways.  What does DT give natively, what must be added, and
# does DT's one native weak ingredient -- the chiral zero-mode neutrino
# (weak_chirality_u4, neutrino_helicity) -- survive the Z-pole data?
#
#   [A] DT natively: no neutral weak coupling (nu: gravity only; e, nucleons:
#       parity-even) -> NC/CC ~ 1e-67, atomic weak charge Q_W = 0
#   [B] the SM's neutral sector: ONE extra parameter (sin^2 thW) fixes all
#       four observables, and rho = 1
#   [C] DT's native test: Z invisible width counts neutrino CHIRALITIES --
#       one 4x4 irrep (left only) vs Williams' doubled 8x8 (left + right)
#   [D] parameter count if DT adds a neutral field WITHOUT the SM structure
#
# Measured (imported, rounded; tree-level comparisons, ~2% radiative scale):
#   R_nu = sigma_NC/sigma_CC (nu, isoscalar) ~ 0.31        CDHS/CHARM 1980s
#   g_V^e = -0.035(17), g_A^e = -0.503(17)                 CHARM II 1994
#   Q_W(Cs-133) = -72.8(4)                                 Wood 1997 + atomic theory
#   Gamma_inv/Gamma_ll = 499.0/83.98 MeV, rho_0 = 1.00038(20)   LEP / PDG
#   sin^2 thW (Z-pole leptonic effective) = 0.2315        LEP/SLD
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("weak_neutral_current.sage -- DT and the neutral weak current")
print("=" * 74)

Rnu_m, gVe_m, gAe_m, sg_e = 0.31, -0.035, -0.503, 0.017
QW_m, QW_s = -72.8, 0.4
Ginv, Gll, rho_m = 499.0, 83.984, 1.00038
s2 = 0.2315
Zc, Nc = 55, 78                                     # Cs-133

# ============================================================ [A]
print("\n[A] DT's own fields")
GF, GN = 1.1663787e-5, 6.70883e-39
RNC_DT = (GN / GF)**2                               # nu couples via phi_4 only
print("    nu neutral coupling: gravitational phi_4 only (weak_va_coupling [C])")
print("    NC/CC from DT's own fields ~ (G_N/G_F)^2 = %.1e  vs measured %.2f" % (RNC_DT, Rnu_m))
check(RNC_DT / Rnu_m < 1e-60, "DT natively has NO neutrino neutral current: excluded (Gargamelle)")
QW_DT = 0                                           # parity-even sector (weak_parity)
print("    atomic parity violation: DT's weak sector is parity-even -> Q_W = %d" % QW_DT)
check(abs(QW_DT - QW_m) / QW_s > 100, "Q_W(Cs) = 0 is excluded at > 100 sigma")
print("    So a neutral, parity-violating vector field must be ADDED (as for the W, Q5).")

# ============================================================ [B]
print("\n[B] The SM's neutral sector: one parameter, sin^2 thW = %.4f" % s2)
gL = lambda T3, Q: T3 - Q * s2
gR = lambda Q: -Q * s2
gV = lambda T3, Q: T3 - 2 * Q * s2
Rnu = 0.5 - s2 + RR(20) / 27 * s2**2                # Llewellyn Smith, isoscalar
QW = Zc * (1 - 4 * s2) - Nc
rows = [("R_nu (DIS NC/CC)", Rnu, Rnu_m, 0.03 * Rnu_m),
        ("g_V^e (nu_mu e)", gV(-0.5, -1), gVe_m, sg_e),
        ("g_A^e (nu_mu e)", -0.5, gAe_m, sg_e),
        ("Q_W(Cs) tree", QW, QW_m, 0.03 * abs(QW_m))]
for n, p, m, s in rows:
    print("    %-18s SM tree %+8.4f   measured %+8.4f   pull %+.1f" % (n, p, m, (p - m) / s))
check(all(abs(p - m) <= 2 * s for n, p, m, s in rows),
      "one parameter reproduces all four NC observables (within 2 sigma / 3% tree)")
check(abs(rho_m - 1) < 1e-3, "rho_0 = M_W^2/(M_Z^2 cos^2 thW) = 1 to 4e-4: the doublet-Higgs prediction")

# ============================================================ [C]
print("\n[C] DT's native ingredient: neutrino chirality, counted at the Z")
Gff = lambda gl, gr: gl**2 + gr**2                  # partial width / common factor
r_lep = Gff(gL(0.5, 0), 0) / Gff(gL(-0.5, -1), gR(-1))
Nnu = (Ginv / Gll) / r_lep
print("    Gamma(nu nu)/Gamma(l l) per family, left-only nu: %.3f" % r_lep)
print("    => N_nu = %.3f  (3 families: an input, as in the SM)" % Nnu)
check(abs(Nnu - 3) < 0.03, "ONE 4x4 irrep (left-handed zero mode): 3 families fit the Z width")
r_vec = Gff(gL(0.5, 0), 0.5) / Gff(gL(-0.5, -1), gR(-1))
Nvec = 3 * r_vec / r_lep
sigN = 0.0082
print("    Williams' 8x8 (both blocks: an L AND an R zero mode, same coupling):")
print("    3 families would look like N_nu = %.1f  vs %.4f(82): %.0f sigma"
      % (Nvec, Nnu, (Nvec - Nnu) / sigN))
check((Nvec - Nnu) / sigN > 100, "the doubled 8x8 (vector-like neutrino) is EXCLUDED by the Z width")

# ============================================================ [D]
print("\n[D] If DT adds a neutral field WITHOUT the SM's SU(2)xU(1) structure")
free = {"nu": 1, "e": 2, "u": 2, "d": 2, "M_Z (rho)": 1}    # nu: shape fixed by [C], norm free
nDT, nSM = sum(free.values()), 1
print("    free couplings: %s  -> %d inputs (SM: %d, plus rho = 1 predicted)"
      % (free, nDT, nSM))
check(nDT > len(rows), "DT would need more inputs than the observables used here: no prediction")

print("\n" + "=" * 74)
print("""VERDICT
  1. DT has NO neutral weak current of its own: the neutrino couples only
     gravitationally (NC/CC ~ 1e-67 vs 0.31) and the charged sector is
     parity-even (Q_W(Cs) = 0 vs -72.8) ([A]).  As for the W (Q5), a
     neutral, parity-violating vector field must be added.
  2. The SM adds it with ONE parameter: sin^2 thW = 0.2315 reproduces R_nu,
     g_V^e, g_A^e and Q_W(Cs) at tree level, and the doublet Higgs gives
     rho = 1 ([B]).  Without that structure DT would need ~8 inputs for 4
     observables ([D]) -- consistent but empty.  So DT must adopt the SM's
     SU(2)_L x U(1)_Y couplings T3 - Q sin^2 thW, not merely "a Z".
  3. DT's one native weak ingredient PASSES a precision test: the chiral
     zero mode (one 4x4 irrep, left-handed) gives N_nu = 2.99 with three
     families, while Williams' doubled 8x8 (an L and an R zero mode) would
     show as N_nu ~ 6, excluded at > 300 sigma ([C]).  The data select the
     4x4 irrep (Q3) over Williams' Eq 109 representation.
  STANDING (assumption ledger: SHARED WITH SM, by import).  Just as the SM
  assumes SU(2)_L x U(1)_Y with one mixing angle and a doublet Higgs, if
  DT assumes the same neutral couplings (T3 - Q sin^2 thW, rho = 1), its
  neutral current is consistent with the SM.  Neutral, but say plainly:
  this IMPORTS the SM's electroweak structure; DT supplies none of it.
  DT's own contribution here is the neutrino's chirality (Q3 + Q4), which
  the Z width confirms (selected, not derived).  Ledger Q6 (couplings)
  and Q12 (doublet Higgs, rho = 1 -- which brings the 125 GeV boson and
  forbids a bare charged-fermion mass: open tension with DT's 5D mass).
  Bracket: (i) tree level, ~2% radiative corrections (full SM Q_W(Cs) =
  -73.2); (ii) three families are an input in both theories; (iii) DT has
  no doublet pairing nu with e -- the zero mode and the electron are
  separate DT fields -- so the SU(2) assignment is wholly imported.
""")
print("=" * 74)
print("ALL ASSERTIONS PASSED")
print("=" * 74)
