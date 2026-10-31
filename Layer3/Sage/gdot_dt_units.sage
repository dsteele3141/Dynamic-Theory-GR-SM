#!/usr/bin/env sage
# --------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, LA-8370-MS (LASL, 1980) and Book Eqs
#   3.2/3.3, 3.54-3.55 (the secular (1 - H0 t) drift and the variable unit of
#   action).  The SN candle model and its rate xi: dark_energy_standard_candle,
#   pantheon_dark_energy_fit (DT distance d_L = (c/H0)(1+z) ln(1+z), candle
#   Delta m = (15/4) log10(G_past/G0) = (15/4) xi ln(1+z)/ln10).  The joint
#   drift-Lambda valley: dt_lambda_joint_fit.  The Chandrasekhar-screening
#   no-go: density_screening_nogo, local_gdot_escape_candidates.  alpha and
#   m_e/m_p fixed (selected by clocks): hbar_potential_clocks.
# Data: D. B. Guenther, L. M. Krauss & P. Demarque, ApJ 498 (1998) 871
#   (helioseismology, |Gdot/G| < 1.6e-12 /yr).  Solar age 4.57 Gyr
#   (meteoritic, radiometric = atomic/nuclear time).  [A. Bonanno &
#   H.-E. Froehlich 2020 report Gdot/G = +(1.25 +- 0.30)e-13 /yr, a claimed
#   non-zero value, contested; NOT used.]
# Supplementary Materials: step [P5], Sec. S5 (code index: Sec. S6).
# --------------------------------------------------------------------------

# ============================================================================
# gdot_dt_units.sage                               [A]  (assertion-hardened)
#
# Gdot/G IN DT'S OWN UNITS: does the drift the SN fit needs survive locally?
#
# The open question (future-work 1g(i)): the LLR bound on Gdot/G was imported
# from a GR-framework reduction, not re-derived through DT's drifting units.
# This script removes the unit question instead of answering it piecemeal.
#
#   [A] Units cannot matter.  A local gravitational observable read with
#       atomic clocks and rulers depends on G only through the dimensionless
#       alpha_G = G m_p^2/(hbar c) (with alpha and m_e/m_p, fixed by clocks).
#       Checked for the orbital frequency in atomic units, with invariance
#       under an arbitrary rescaling of the units of M, L and T.
#   [B] The SN candle is the SAME number: M_Ch/m_p = alpha_G^(-3/2), so the
#       fit's xi is a drift of alpha_G, however DT splits it among G, hbar
#       and m.  A drift that cancels in units also removes the candle effect.
#   [C] Same epochs: DT's distance law gives ln(1+z) = H0 t_lb, so the drift
#       is ln alpha_G(now)/alpha_G(past) = -xi H0 t_lb.  The Sun formed at
#       H0 t = 0.33, inside the SN range.  Other t(z) maps (Milne, EdS,
#       LCDM) change the inferred rate by < 30%.
#   [D] Helioseismology (stellar structure in atomic units, i.e. GR + SM,
#       which DT recovers in the solar interior at u4 = 0 to ~1e-31): the
#       fitted xi = 0.23 is excluded by ~6-8x; Williams' xi = 1 by ~30x.
#   [E] Escapes: unit drift (removes the candle too), density screening
#       (NO-GO, density_screening_nogo), a different clock (< 30%).
#   [F] Consequence (dt_lambda_joint_fit): with xi ~ 0 the SN data need a
#       NEGATIVE Lambda (~ -0.5 x Planck) in DT's distance law, with chi2 as
#       good as the drift fit.
#
# RUN: sage gdot_dt_units.sage
# ============================================================================

from sage.all import *


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


def z(e):
    return bool(SR(e).canonicalize_radical().simplify_full() == 0)


print("=" * 72)
print("gdot_dt_units.sage -- Gdot/G in DT's own units")
print("=" * 72)

# ---------------------------------------------------------------- [A]
print("\n[A] Local observables see only alpha_G")
G, hb, c, mp, me, al, Nb = var('G hb c mp me al Nb', domain='positive')
aB = hb / (me * c * al)                      # Bohr radius
tA = hb / (me * c ** 2)                      # atomic time unit
A_, nu = var('A_ nu', domain='positive')     # orbit radius, frequency in atomic units
# Kepler: n^2 a^3 = G M, with a = A_ aB, n = nu / tA, M = Nb mp
nu2 = (G * Nb * mp / (A_ * aB) ** 3) * tA ** 2
alG = G * mp ** 2 / (hb * c)
target = Nb * al ** 3 * alG * (me / mp) / A_ ** 3
ok("A1: orbital frequency in atomic units: nu^2 = N alpha^3 alpha_G (m_e/m_p) / A^3",
   z(nu2 - target))
mu, l, t = var('mu l t', domain='positive')
dims = {G: (-1, 3, -2), hb: (1, 2, -1), c: (0, 1, -1), mp: (1, 0, 0), me: (1, 0, 0)}
sub = {s: s * mu ** d[0] * l ** d[1] * t ** d[2] for s, d in dims.items()}
ok("A2: nu^2 is invariant under any rescaling of the units of M, L, T",
   z(nu2.subs(sub) / nu2 - 1))
ok("A3: alpha_G is invariant too", z(alG.subs(sub) / alG - 1))
print("      a drift that is only a unit change is invisible to every local")
print("      dimensionless observable -- and, by [B], to the SN candle as well")

# ---------------------------------------------------------------- [B]
print("\n[B] The Chandrasekhar candle is the same number")
MCh = (hb * c / G) ** (QQ(3) / 2) / mp ** 2
ok("B1: M_Ch / m_p = alpha_G^(-3/2)", z(MCh / mp - alG ** (-QQ(3) / 2)))
ok("B2: invariant under unit rescaling", z((MCh / mp).subs(sub) / (MCh / mp) - 1))
print("      so the fit's xi is an alpha_G drift, however DT shares it among")
print("      G, hbar and m (Williams' hbar(Phi) is a potential effect,")
print("      hbar_potential_clocks, ~1e-8, and cannot cancel it)")

# ---------------------------------------------------------------- [C]
print("\n[C] The Sun samples the SN epochs")
H0 = 70 * 1e3 / 3.0857e22 * 3.15576e7          # /yr
tsun = 4.57e9
x = H0 * tsun
zsun = RR(exp(x) - 1)
print("      H0 = %.3e /yr ; H0 t_sun = %.3f ; DT redshift of the Sun's birth z = %.3f"
      % (H0, x, zsun))
ok("C1: the Sun's birth lies inside the Pantheon range (0.01 < z < 2.3)", 0.01 < zsun < 2.3)
milne = zsun / (1 + zsun)
eds = (2.0 / 3) * (1 - (1 + zsun) ** (-1.5))
lcdm = numerical_integral(lambda q: 1 / ((1 + q) * sqrt(0.3 * (1 + q) ** 3 + 0.7)), 0, zsun)[0]
rat = [RR(r / log(1 + zsun)) for r in (milne, eds, lcdm)]
print("      H0 t_lb / ln(1+z) at that z: Milne %.2f, EdS %.2f, LCDM %.2f (DT 1)"
      % tuple(rat))
fmin = min(rat)
ok("C2: any standard t(z) changes the inferred rate by < 30%", fmin > 0.7)

# ---------------------------------------------------------------- [D]
print("\n[D] Helioseismology")
bound = 1.6e-12                                  # /yr, Guenther et al. 1998
ximax = bound / H0 / fmin                        # most generous mapping
print("      |Gdot/G| < %.1e /yr -> |xi| < %.3f (DT clock), < %.3f (most generous t(z))"
      % (bound, bound / H0, ximax))
for name, xi in [("fitted, Lambda = 0", 0.229), ("Williams (1 - H0 t)", 1.0)]:
    print("      xi = %.3f (%s): excluded by %.1fx; ln alpha_G change over the Sun's life %.3f"
          % (xi, name, xi / ximax, xi * x))
ok("D1: the fitted drift xi = 0.229 is excluded by > 5x", 0.229 / ximax > 5)
ok("D2: Williams' xi = 1 is excluded by > 25x", 1.0 / ximax > 25)
dL = exp(7 * 0.229 * x)
print("      (homology L ~ alpha_G^7: the young Sun %.1fx brighter than today," % dL)
print("       against the standard ~30% fainter)")
ok("D3: the fitted drift would make the young Sun brighter than today", dL > 1)
print("      the inference uses DT's own recovered sector (GR + SM at u4 = 0 to")
print("      ~1e-31 in the solar interior); no GR-framework bound is imported")

# ---------------------------------------------------------------- [E]
print("\n[E] Escapes")
print("      unit drift only   -> alpha_G fixed -> no candle effect either ([A],[B])")
print("      density screening -> screens the white dwarf before the Sun")
print("                           (NO-GO, density_screening_nogo)")
ok("E1: a white dwarf (~1e6 g/cc) is denser than the solar core (150 g/cc) by > 1e3",
   1e6 / 150 > 1e3)
print("      a different clock -> covered by [C] (< 30%)")

# ---------------------------------------------------------------- [F]
print("\n[F] What the SN data then need (dt_lambda_joint_fit, full covariance)")
Lslope, Lplanck = 1.49, 0.320
Lmax = (ximax - 0.229) / Lslope
Lbest, chi2_0, chi2_fit = -0.158, 47.3, 47.1
print("      valley xi = 0.229 + %.2f L: xi < %.3f -> L < %.2f; at xi = 0, L = %.3f"
      % (Lslope, ximax, Lmax, Lbest))
print("      Lambda = %.2f x Planck; chi2 %.1f vs %.1f for the drift fit"
      % (Lbest / Lplanck, chi2_0, chi2_fit))
ok("F1: with the drift bounded the SN data need a NEGATIVE Lambda", Lmax < 0 and Lbest < 0)
ok("F2: and fit as well as the drift did (Delta chi2 < 1)", abs(chi2_0 - chi2_fit) < 1)

print("\n" + "=" * 72)
print("SUMMARY")
print("  Local Gdot/G tests read with atomic clocks measure alpha_G, and so")
print("  does the SN candle: DT's units cannot hide one from the other.")
print("  Helioseismology, in DT's recovered sector, gives |xi| < %.3f: the" % ximax)
print("  fitted drift (0.23) is excluded by ~%.0fx, Williams' (1 - H0 t) by ~%.0fx."
      % (0.229 / ximax, 1 / ximax))
print("  The drift is not DT's dark energy.  Without it the SN data need DT's")
print("  own Lambda, NEGATIVE (~ -0.5 Planck) in DT's distance law.")
print("=" * 72)
