#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, LA-8370-MS (LASL, 1980), Section III,
#   Eqs (71)-(77): the 5D continuity d_i J^i = 0, unpacked as Eq (72)
#   d rho/dt + div J + a0 dJ4/d gamma = 0 (prediction P6, T2_2_charge_
#   conservation.sage), with x^4 = gamma/a0. Book (2011) Eq (1.36)/(10.6):
#   u^4 = 0 is conservation of mass; the (1 - H0 t) secular drift of the
#   mass potential (Book Ch. 3) is DT's only native u^4 != 0.
# Builds on: weak_chirality_u4 / noflux_boundary_derivation (u^4 = 0 is the
#   no-flux edge J^4(0) = 0, ledger Q4), heavy_massaxis_partners (zero mode
#   e^{-M y}, M >~ 1 TeV indicative, ledger Q17), kk_radion_g44 (the free
#   radion gives gamma_PPN = 1/2), gravitational_radiation (the breathing
#   mode; the gauge-sector coupling K_g gamma ~ 1e-31 in the laboratory).
# Data (imported): Borexino, PRL 115 (2015) 231802: tau(e -> nu gamma)
#   > 6.6e28 yr (90% CL). Cassini, Bertotti et al., Nature 425 (2003) 374:
#   gamma_PPN - 1 = (2.1 +/- 2.3)e-5. LIGO/Virgo GW170814, PRL 119 (2017)
#   141101: first polarisation test; pure tensor favoured over pure scalar
#   and pure vector. Mixed-polarisation amplitudes are bounded only at
#   order 0.1-1 by current analyses (used here as an ORDER, not a number).
# ---------------------------------------------------------------------------
#
# u4_relaxation_bounds.sage                           [A] (assertion-hardened)
#
# QUESTION (manuscript item 5). Relaxing u^4 = 0 is where DT should predict.
# Two consequences can be confronted with data:
#   (i)  P6: charge may move along the mass axis. Does the electron's
#        stability (tau > 6.6e28 yr) bound it, and what does DT itself give?
#   (ii) the breathing (radion) polarisation of gravitational waves that the
#        restriction removes. What do LIGO's polarisation tests bound, and
#        what does DT itself give?
# Odds before computing: (i) a null result plus a bound ~90%, something
# sharp ~5%; (ii) a bound with a native value far below reach ~90%, with
# Cassini stronger than LIGO on the same parameter ~80%.
#
#   [A] 4D charge from the 5D continuity: Q changes ONLY by flux through the
#       edge gamma = 0; DT's native drift (a dilation of gamma) has zero flux
#       there, so 4D charge is exactly conserved
#   [B] a leak at the edge, if one were added: electron stability bounds it
#   [C] breathing mode: one parameter eps (eps = 1 is the free radion),
#       gamma_PPN = 1 - 2 alpha^2/(1+alpha^2), alpha^2 = eps/3; Cassini vs LIGO
#   [D] DT's native eps: K_g gamma, laboratory and neutron star
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("u4_relaxation_bounds.sage -- relaxing u4 = 0: charge leak and breathing mode")
print("=" * 74)

hbar, c, eV = 1.054571817e-34, 2.99792458e8, 1.602176634e-19
yr = 3.15576e7
H0 = 2.2e-18                    # s^-1 (~68 km/s/Mpc)
a0 = 3.65e7                     # kg/m^4 (Layer-3 script 18)

# ============================================================ [A]
print("\n[A] 4D charge from Eq (72), integrated over space: d rho/dt + d J4/d y = 0")
print("    (y = x^4 = gamma/a0 >= 0, so a0 d/d gamma = d/d y). Q(t) = int_0^oo rho dy.")
t, y, H, v0 = var('t y H v0', domain='real')
assume(H > 0); assume(y >= 0)
f = function('f')
# DT's native u4 != 0: the secular drift rescales mass density, gamma -> gamma (1 - H t),
# i.e. the 5D velocity along y is  dy/dt = -H y  (proportional to y: zero at the edge)
vel = -H * y
rho_dil = exp(H * t) * f(y * exp(H * t))          # transported profile
J4_dil = rho_dil * vel
cont = diff(rho_dil, t) + diff(J4_dil, y)
check(bool(cont.simplify_full() == 0),
      "the dilated profile rho = e^{Ht} f(y e^{Ht}) solves Eq (72) with u = -H y")
check(bool(J4_dil.subs(y=0) == 0),
      "its flux through the edge y = 0 vanishes identically (the velocity is prop. to y)")
# Q(t) = int_0^oo e^{Ht} f(y e^{Ht}) dy = int_0^oo f(s) ds  (s = y e^{Ht}): t-independent
s = var('s')
check(bool((exp(H * t) * diff(s * exp(-H * t), s)).simplify_full() == 1),
      "Q(t) = int f(s) ds after s = y e^{Ht}: 4D charge EXACTLY conserved under the drift")
# general statement: dQ/dt = -[J4]_0^oo = J4(0) for a profile that decays
print("    In general dQ/dt = -int_0^oo dJ4/dy dy = J4(0) - J4(oo) = J4(0).")
print("    So 4D charge can change ONLY through the edge. The no-flux edge (u4 = 0,")
print("    ledger Q4) closes it; DT's own u4 != 0 (the drift) has no edge flux either,")
print("    because zero mass density stays zero under a rescaling of mass density.")
print("    => P6 is a statement about the 5D distribution of charge along the mass")
print("       axis; it predicts NO 4D charge non-conservation.")

# ============================================================ [B]
print("\n[B] A leak at the edge (velocity u4_edge * c at y = 0), if one were added")
M_GeV = 1000.0                                   # 1 TeV, indicative (Q17)
M_inv_m = M_GeV * 1e9 * eV / (hbar * c)          # 1/m
# zero mode |psi|^2 ~ e^{-2 M y}: int rho dy = rho(0)/(2M); Gamma = J4(0)/Q = 2 M c u4
tau_e = 6.6e28 * yr
u4_edge_max = 1 / (2 * M_inv_m * c * tau_e)
print("    zero mode rho ~ e^{-2My}: Gamma = J4(0)/Q = 2 M c u4_edge")
print("    M = 1 TeV -> M = %.3e /m; tau_e > 6.6e28 yr -> u4_edge < %.1e" % (M_inv_m, u4_edge_max))
check(u4_edge_max < 1e-60, "an edge leak is bounded below 1e-60 (tighter for larger M)")
u4_drift_lab = 2e-31
print("    For comparison DT's drift has u4 ~ 2e-31 in the laboratory, but zero AT the")
print("    edge ([A]); a constant leak at the edge is NOT something DT supplies.")
check(u4_drift_lab / u4_edge_max > 1e30,
      "a leak of the drift's size at the edge would be excluded by >1e30")

# ============================================================ [C]
print("\n[C] Breathing mode: relax h_A4 = 0 by a fraction eps (eps = 1: the free radion)")
al2 = var('al2', domain='positive'); eps = var('eps', domain='positive')
gam = 1 - 2 * al2 / (1 + al2)
omega = (1 / al2 - 3) / 2
check(bool((gam - (1 + omega) / (2 + omega)).simplify_full() == 0),
      "gamma_PPN = 1 - 2 alpha^2/(1+alpha^2) is Brans-Dicke (1+w)/(2+w)")
check(gam.subs(al2=QQ(1) / 3) == QQ(1) / 2,
      "alpha^2 = 1/3 (w = 0, the KK radion) gives gamma = 1/2 (kk_radion_g44)")
gam_eps = gam.subs(al2=eps / 3)
cassini = 2.1e-5 + 2 * 2.3e-5                   # |gamma - 1| at ~2 sigma
eps_cassini = find_root(lambda e: abs(gam_eps.subs(eps=e) - 1) - cassini, 1e-12, 1)
print("    gamma - 1 = -2 eps/(3 + eps);  Cassini |gamma - 1| < %.1e -> eps < %.1e"
      % (cassini, eps_cassini))
eps_ligo = 0.1                                  # order of current mixed-mode limits
print("    LIGO-Virgo polarisation tests: pure breathing excluded; mixed scalar")
print("    amplitude bounded only at O(0.1-1). Taking eps_LIGO ~ %.1f (optimistic):" % eps_ligo)
check(eps_ligo / eps_cassini > 500,
      "on the SAME parameter Cassini is ~1e3 tighter than the polarisation tests")

# ============================================================ [D]
print("\n[D] DT's native eps: the mass-coordinate coupling K_g gamma, K_g = H0/(a0 c)")
Kg = H0 / (a0 * c)
for name, g_ in (("laboratory (1e3 kg/m^3)", 1e3), ("neutron star (5e17 kg/m^3)", 5e17)):
    print("    %-28s K_g gamma = %.1e" % (name, Kg * g_))
eps_ns = Kg * 5e17
check(Kg * 1e3 < 1e-29 and eps_ns < 1e-15,
      "native eps: ~1e-31 in the laboratory, ~1e-16 in a neutron star")
check(eps_ns / eps_cassini < 1e-10,
      "even for a neutron-star merger the native breathing fraction is >1e10 below Cassini")

print("""
NET.
  * P6 (Williams' modified continuity, Eq 72) moves charge along the mass
    axis only. The 4D charge changes solely by flux through the edge
    gamma = 0. The no-flux edge (u4 = 0) closes it, and DT's one native
    u4 != 0, the secular drift, rescales gamma and so has zero velocity at
    the edge. 4D charge is EXACTLY conserved: a null result, consistent with
    the electron's stability. Any edge leak one might add is bounded by the
    electron lifetime below ~1e-64 (M = 1 TeV): a BOUND, not a prediction.
  * The breathing mode is the relaxed restriction h_A4 = 0. On one
    parameter eps (eps = 1 the free radion, gamma_PPN = 1/2), Cassini gives
    eps < ~1e-4 and is ~1e3 tighter than LIGO-Virgo polarisation tests. DT's
    own eps is the mass-coordinate coupling, ~1e-31 in the lab and ~1e-16 in
    a neutron star: far below both. No prediction within reach.
  * Both relaxations are consistent with data and neither yields a
    measurable signal: what relaxing u4 = 0 leaves in 4D is too small to see.
""")
print("ALL ASSERTIONS PASSED")
