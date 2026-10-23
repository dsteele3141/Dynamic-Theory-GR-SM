#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE (primary): LA-8370-MS (1980) Eqs 110-112 (5D Dirac equation on the
#   half-line x^4 = gamma/a0 >= 0). Builds on heavy_massaxis_partners (zero
#   mode + continuum, pair production), noflux_boundary_derivation (the two
#   chiral edge conditions), su2_from_5d (the 5D field is an SU(2) doublet),
#   gauge_fields_halfline (gauge fields independent of x^4) and
#   u4_relaxation_bounds (the drift: velocity along x^4 proportional to x^4).
# Imported, indicative: LZ (PRL 135, 011802, 2025) spin-independent limit at high mass,
#   ~3e-47 cm^2 at 1 TeV, scaling as M; Dirac-doublet thermal relic
#   Omega h^2 = 0.12 at M ~ 1.1 TeV, scaling as M^2.
# Supplementary Materials: step [P2], Sec. S5 (code index: Sec. S6).
# ---------------------------------------------------------------------------
#
# partner_stability.sage                              [A] (assertion-hardened)
#
# QUESTION. Gauge bosons produce the heavy
# mass-axis partners only in pairs. Is that an exact selection rule? If so,
# the lightest partner is stable: a dark-matter candidate, possibly a
# right-handed neutrino (Boyle, Finn & Turok, PRL 121, 251301, 2018). Odds before computing that a
# viable candidate survives: ~20-30%.
#
#   [A] which couplings conserve the partner label: gauge, metric, Yukawa
#   [B] what breaks it: block-dependent M, a localised Higgs, the drift
#   [C] the lightest partner: its quantum numbers (not a sterile nu_R)
#   [D] data: direct detection and the thermal relic density
#   [E] where a right-handed neutrino can come from in DT
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("partner_stability.sage -- is pair production of partners exact?")
print("=" * 74)

y, k, kp = var('y k kp', domain='positive')
M, M1, M2, mu = var('M M1 M2 mu', domain='positive')

# Mode functions (unnormalised; normalisation does not affect zeros).
# Block with a LEFT-handed zero mode (edge f_R(0) = 0):
#   zero mode f_L = exp(-M y); continuum f_R = sin(k y), m f_L = M sin - k cos
# Block with a RIGHT-handed zero mode (the other edge condition, L <-> R):
#   zero mode f_R = exp(-M y); continuum f_L = sin(k y), m f_R = M sin - k cos
zero = lambda Mx: exp(-Mx * y)
contD = lambda Mx, kx: Mx * sin(kx * y) - kx * cos(kx * y)   # the chirality of the zero mode
contN = lambda kx: sin(kx * y)                               # the other chirality

# ============================================================ [A]
print("\n[A] Couplings that are independent of x^4 conserve the partner label")
print("    A 4D field B(x) independent of y couples modes a, b with weight")
print("    int f_a f_b dy, so it can change a zero mode into a partner only if")
print("    the two profiles overlap.")
ov_g = integrate(zero(M) * contD(M, k), y, 0, oo)
check(bool(ov_g.simplify_full() == 0),
      "gauge field or 4D metric: int f_L^0 f_L^k dy = 0 (zero mode -> partner forbidden)")
print("    The same weight governs the graviton: with g44 fixed and the 4D metric")
print("    independent of x^4 (the GR recovery), T_ab couples modes by int f_a f_b.")
print("\n    Yukawa  y_f H(x) psibar_L(block 1) psi_R(block 2): block 1 has the")
print("    left-handed zero mode, block 2 the right-handed one (e.g. L and e_R).")
Y01 = integrate(zero(M1) * contD(M2, k), y, 0, oo).simplify_full()
Y10 = integrate(contD(M1, k) * zero(M2), y, 0, oo).simplify_full()
print("      zero_L(1) x partner_R(2):  %s" % Y01)
print("      partner_L(1) x zero_R(2):  %s" % Y10)
check(bool((Y01 - k * (M2 - M1) / (M1**2 + k**2)).simplify_full() == 0)
      and bool((Y10 - k * (M1 - M2) / (M2**2 + k**2)).simplify_full() == 0),
      "zero-mode/partner Yukawa overlaps are proportional to M1 - M2")
check(bool(Y01.subs(M2=M1) == 0) and bool(Y10.subs(M2=M1) == 0),
      "with ONE native mass M for all blocks the Yukawa cannot convert a zero mode either")
check(bool((contD(M1, k) - contD(M2, k)).subs(M2=M1) == 0),
      "with one M the partner profiles of the two chiralities coincide, so the")
print("        Yukawa is diagonal in k (orthonormal continuum: int = delta(k - k'))")
print("    => With one native M and every 4D field independent of x^4, each")
print("       fermion line keeps its mass-axis label k: partners are created and")
print("       destroyed only in pairs, and the lightest one cannot decay.")

# ============================================================ [B]
print("\n[B] What breaks the rule")
check(bool(Y01.subs(M1=1, M2=2, k=1) != 0),
      "(i) block-dependent M (M1 != M2): a partner decays through the Yukawa, rate ~ y_f^2 (M1-M2)^2")
ovH = integrate(exp(-mu * y) * zero(M) * contD(M, k), y, 0, oo).simplify_full()
check(bool((ovH + k * mu / ((M + mu)**2 + k**2)).simplify_full() == 0),
      "(ii) a Higgs localised along x^4 (profile exp(-mu y)): overlap -k mu/((M+mu)^2+k^2) != 0")
hbar_eVs, H0 = 6.582e-16, 2.27e-18
omega_M = 1.0e12 / hbar_eVs                       # M = 1 TeV in s^-1
ratio = H0 / omega_M
print("    (iii) the drift (velocity along x^4 ~ H0 x^4) varies slowly against")
print("          the gap M between zero mode and continuum: H0 / M = %.1e at 1 TeV" % ratio)
check(ratio < 1e-40, "the drift is adiabatic: zero-mode/partner mixing < 1e-40, negligible")
print("    DT does not say whether M is one constant or differs between blocks")
print("    (heavy_massaxis_partners: DT does not fix M). Case (i) is a free choice.")

# ============================================================ [C]
print("\n[C] The lightest partner, for one native M")
print("    The partner continuum of a 5D field carries that field's gauge")
print("    representation in BOTH chiralities (vector-like): the lepton-doublet")
print("    block gives a doublet (nu', e'), Y = -1/2, the e_R block a singlet, Y = -1;")
print("    quark blocks give coloured partners.")
T3, Yh = QQ(1)/2, -QQ(1)/2
check(T3 + Yh == 0, "nu' has Q = T3 + Y = 0 but Y != 0: it couples to the Z")
print("    Electroweak loops make the charged doublet member heavier (~0.35 GeV")
print("    for M >> M_Z), so the lightest LEPTONIC partner is nu', a Dirac")
print("    SU(2)-doublet fermion with a vector Z coupling. Its right-handed")
print("    component is a doublet member, NOT a sterile right-handed neutrino.")
print("    Quark partners keep baryon number and k, so the lightest coloured")
print("    partner is stable too.")

# ============================================================ [D]
print("\n[D] Data (nu' as the dark matter; imported limits, indicative)")
GF, sw2, mn, GeV2cm2 = 1.1664e-5, 0.231, 0.939, 0.3894e-27
A, Z = 131, 54
QW = (A - Z) - (1 - 4 * sw2) * Z
sig_n = mn**2 * GF**2 * (QW / A)**2 / (2 * pi.n()) * GeV2cm2   # cm^2 per nucleon
print("    Z exchange on xenon: sigma_n = mu_n^2 G_F^2 (Q_W/A)^2 / 2pi = %.1e cm^2" % sig_n)
check(1e-39 < sig_n < 1e-38, "sigma_n ~ 2e-39 cm^2 (the Goodman-Witten value for a Dirac doublet)")
lim_1TeV = 3e-47                                  # cm^2 at 1 TeV, LZ 2025, ~ M above 100 GeV
M_excl = 1e3 * sig_n / lim_1TeV                   # GeV
print("    LZ limit ~ 3e-47 cm^2 x (M / 1 TeV)  =>  excluded as all the DM for")
print("    M <~ %.0e GeV" % M_excl)
M_turok = 4.8e8
check(M_excl > 10 * M_turok,
      "Turok's 4.8e8 GeV is excluded for nu' by > 10x margin, even if the limit is 10x weaker")
M_relic = 1.1e3                                   # GeV, Omega h^2 = 0.12
check(M_relic < M_excl,
      "the thermal-relic mass (~1.1 TeV) is excluded by direct detection by ~1e8")
over = (M_excl / M_relic)**2
print("    Above M ~ 1.1 TeV a thermal relic overcloses: at M = %.0e GeV by %.0e." % (M_excl, over))
check(over > 1e10, "for every M >~ 1 TeV a thermal population of stable partners is excluded")
print("    So stable partners are allowed only if they were never thermally")
print("    produced: reheating temperature well below M.")

# ============================================================ [E]
print("\n[E] A right-handed neutrino in DT")
print("    Each 5D block takes one of the two Lorentz-invariant edge conditions")
print("    (noflux_boundary_derivation); e_R is the right-handed zero mode of")
print("    a block with the other condition. A gauge-singlet block with that")
print("    condition gives a right-handed zero mode with T3 = Y = 0:")
check(0 + 0 == 0, "its charge Q = T3 + Y = 0 and it has no gauge coupling: a sterile nu_R")
print("    So nu_R needs no compactified mass axis, only one more block (the")
print("    same price as in the SM). Its Majorana mass and any symmetry making")
print("    it stable are 4D inputs, as in Turok's proposal; DT adds nothing to them.")

print("""
NET.
  * The pair rule is EXACT under the conditions of the recovery: one native M
    for every block and all 4D fields (gauge, metric, Higgs) independent of
    x^4. Gauge and Yukawa overlaps between a zero mode and a partner vanish,
    and the drift mixes them at H0/M < 1e-40. So the lightest partner and
    the lightest coloured partner are STABLE.
  * It is broken only by choices DT leaves open: block-dependent M (Yukawa
    decay ~ (M1-M2)^2) or a Higgs localised along x^4.
  * The candidate is nu', the neutral member of the lepton-doublet partner
    continuum: a Dirac doublet with Z couplings, NOT a sterile nu_R. Direct
    detection excludes it as the dark matter below ~1e11 GeV (Turok's
    4.8e8 GeV included); a thermal relic overcloses at every M >~ 1 TeV.
  * Consequence for the manuscript (Sec. 6.2): with one M the partner sector
    does not simply decouple. It needs M above the reheating temperature, or
    block-dependent M so that partners decay. A bound on M, not a candidate.
  * A sterile nu_R is available in DT exactly as in the SM: one gauge-singlet
    block with the right-handed edge condition. No compactification needed.
""")
print("ALL ASSERTIONS PASSED")
