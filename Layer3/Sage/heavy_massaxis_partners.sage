#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE (primary): LA-8370-MS (1980) Eq 105 (mass shell p_j p^j = gamma^2 c^2),
#   Eqs 110-112 ((i gamma^j d_j + 1) psi = 0, x^4 = gamma/a0 >= 0). Builds on
#   weak_chirality_u4 / noflux_boundary_derivation (zero mode exp(-M y)),
#   mass_term_su2 (the native mass survives only in the heavy modes) and
#   gauge_fields_halfline (gauge fields independent of y). a0 = 3.65e7 kg/m^4
#   (Layer-3 script 18).
# Collider bounds (imported, indicative): vector-like leptons excluded to
#   ~1 TeV (CMS 2022), vector-like quarks to ~1.5 TeV (ATLAS/CMS 2022-23).
# ---------------------------------------------------------------------------
#
# heavy_massaxis_partners.sage                        [A] (assertion-hardened)
#
# QUESTION (future-work 4). Above each chiral zero mode, DT's 5D field has
# heavy vector-like partners whose scale is the native mass M. Does DT fix M?
# With a scale this is a collider prediction; without one it is untestable.
# Odds before computing: a usable scale ~10%.
#
#   [A] spectrum on the half-line: ONE zero mode + a CONTINUUM from m = M
#       (no discrete tower, no other bound state)
#   [B] y-independent gauge fields cannot turn a zero mode into a partner
#       (orthogonality): partners are produced in pairs, above 2M
#   [C] candidate scales: Williams' "+1" in Compton units (M = m_f) and DT's
#       own density 1/K_gamma are excluded; M = a0 / gamma* for other
#       densities spans every M, so it selects nothing
#   [D] data: M >~ 1 TeV, so the zero mode lives within gamma* = a0/M
#       <~ 7e-12 kg/m^3 of the edge
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("heavy_massaxis_partners.sage -- does DT fix the heavy-partner scale M?")
print("=" * 74)

y, M, m, k = var('y M m k', domain='positive')
hbarc = 1.97327e-7                    # eV m
a0 = 3.65e7                           # kg / m^4
inv_m = lambda E_eV: E_eV / hbarc     # mass (eV) -> 1/m

# ============================================================ [A]
print("\n[A] Spectrum of (i g^mu d_mu + i g^4 d_4 - M) psi = 0 on y >= 0")
print("    Chiral split psi = f_L(y) psi_L(x) + f_R(y) psi_R(x), 4D mass m:")
print("      (d_y + M) f_L = m f_R,  (-d_y + M) f_R = m f_L,  edge f_R(0) = 0")
fL, fR = function('fL')(y), function('fR')(y)
# zero mode: m = 0
z = exp(-M * y)
check(bool(diff(z, y) + M * z == 0) and bool(integrate(z**2, y, 0, oo) == 1 / (2 * M)),
      "m = 0: f_L = exp(-M y), normalisable; f_R = 0 meets the edge condition")
# massive modes: f_R obeys -f'' + M^2 f = m^2 f with f_R(0) = 0
fRk = sin(k * y)
check(bool((-diff(fRk, y, 2) + M**2 * fRk - (M**2 + k**2) * fRk).simplify_full() == 0)
      and bool(fRk.subs(y=0) == 0),
      "f_R = sin(k y) with m^2 = M^2 + k^2: bounded for every real k > 0")
print("    For m < M the solutions grow or decay as exp(+-sqrt(M^2 - m^2) y);")
q = var('q', domain='positive')                  # q = sqrt(M^2 - m^2) > 0
check(bool((-diff(sinh(q * y), y, 2) + M**2 * sinh(q * y)
            - (M**2 - q**2) * sinh(q * y)).simplify_full() == 0),
      "sinh(q y) solves the f_R equation with m^2 = M^2 - q^2 < M^2")
check(bool(sinh(q * y).subs(y=0) == 0) and bool(limit(sinh(q * y), y=oo) == oo),
      "the only m < M solution meeting f_R(0) = 0 is sinh, which diverges: no bound state")
print("    => one chiral zero mode, then a CONTINUUM of vector-like (L+R) states")
print("       from m = M. Not a Kaluza-Klein tower: the half-line is not compact.")

# ============================================================ [B]
print("\n[B] Gauge couplings between zero mode and continuum")
fLk = M * sin(k * y) - k * cos(k * y)            # m f_L = (-d_y + M) f_R
check(bool((-diff(sin(k * y), y) + M * sin(k * y) - fLk).simplify_full() == 0)
      and bool((diff(fLk, y) + M * fLk - (M**2 + k**2) * sin(k * y)).simplify_full() == 0),
      "continuum: m f_L = M sin(k y) - k cos(k y), and (d_y + M) f_L = m f_R closes")
ov = integrate(exp(-M * y) * fLk, y, 0, oo)
check(bool(ov.simplify_full() == 0),
      "int f_L^0 f_L^k dy = 0: a y-independent gauge field cannot turn a zero mode into a partner")
print("    So W, Z, photon and gluons produce partners only in PAIRS (threshold 2M);")
print("    single production needs Yukawa mixing, suppressed by (m_f / M)^2.")

# ============================================================ [C]
print("\n[C] Candidate scales for M")
mf = {"e": 0.511e6, "mu": 105.66e6, "tau": 1776.9e6}
print("    (i) Williams' '+1' in Compton units: M = m_f for each fermion")
for n, v in mf.items():
    print("        %-4s partners from %.3g GeV" % (n, v / 1e9))
check(mf["e"] < 45.6e9,
      "M = m_e gives charged partners below M_Z/2: excluded by the Z width and LEP pair searches")
print("    (ii) M = a0 / gamma* for a physical mass density gamma* (x^4 = gamma/a0)")
H0, c, G = 2.27e-18, 2.998e8, 6.674e-11
rho_c = 3 * H0**2 / (8 * pi.n() * G)
invK = a0 * c / H0                                 # 1/K_gamma, K_gamma = H0/(a0 c)
dens = [("1/K_gamma (DT's own)", invK), ("Planck", 5.16e96), ("nucleus", 2.3e17),
        ("water", 1.0e3), ("air", 1.2), ("interstellar", 1.0e-21), ("critical", rho_c)]
Mev = {}
for n, rho in dens:
    Mev[n] = a0 / rho * hbarc                      # eV
    print("        %-22s rho = %.2e kg/m^3  ->  M = %.2e eV" % (n, rho, Mev[n]))
check(Mev["1/K_gamma (DT's own)"] < 1e-30,
      "DT's own density scale 1/K_gamma gives M = hbar H0 / c^2 ~ 1e-33 eV: excluded")
check(min(Mev.values()) < 1e12 < max(Mev.values()),
      "physical densities give M from 1e-96 eV to 1e27 eV: 1 TeV is inside, nothing selects it")
print("        (M at the critical density, ~8e17 GeV, is near the Planck scale;")
print("         a coincidence of H0, G and a0, not a derivation)")
print("    (iii) a0 with G, hbar, c: a0 enters only as a density per length, and")
print("          any further scale needs a second input; none in the source.")

# ============================================================ [D]
print("\n[D] Data")
Mmin = 1.0e12                                      # eV, indicative (VLL ~1 TeV)
gstar = a0 / inv_m(Mmin)
print("    M >~ 1 TeV  =>  zero-mode width 1/M <~ %.1e m along x^4" % (1 / inv_m(Mmin)))
print("               =>  gamma* = a0 / M <~ %.1e kg/m^3" % gstar)
check(gstar < 1e-11, "the zero mode sits within gamma* < 1e-11 kg/m^3 of the edge")
print("    So x^4 of a field mode cannot be the ambient density of the matter")
print("    the particle sits in (water is at x^4 = 2.7e-5 m, where exp(-M y) ~")
print("    exp(-1e14)). The edge reading (Q4) needs only gamma >= 0, not this.")
check(inv_m(Mmin) * (1.0e3 / a0) > 1e13, "M y(water) > 1e13 for M = 1 TeV")
print("""
NET.
  * DT's 5D field on the half-line has one chiral zero mode and, above it, a
    CONTINUUM of vector-like partners from m = M. Gauge bosons produce them
    only in pairs. These are DT's heavy mass-axis states.
  * DT does not fix M. Williams' own normalisation (M = m_f) puts charged
    partners at the fermion mass, excluded by LEP and the Z width (a Williams
    discrepancy under the recovered, Yukawa-mass reading). DT's own density
    scale 1/K_gamma gives M ~ hbar H0/c^2: excluded. Any other density
    gamma* gives M = a0/gamma*, and physical densities span every M, so
    picking one is choosing M.
  * Data bound M >~ 1 TeV (indicative: the LHC bounds are for discrete
    vector-like states; a continuum without a mass peak needs a recast).
  * Consequence for the reading: x^4 of a field mode is not the ambient mass
    density; the zero mode lives within ~1e-11 kg/m^3 of the edge.
  * Verdict: NO PREDICTION. M is a free parameter with a lower bound. The
    qualitative signature, a continuum rather than a tower, would distinguish
    DT from compact extra dimensions if partners were ever seen.
""")
print("ALL ASSERTIONS PASSED")
