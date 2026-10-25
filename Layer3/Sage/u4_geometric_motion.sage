#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE (primary): LA-8370-MS (1980): x^4 = gamma/a0, u^i = dx^i/ds, the
#   equation of motion delta int gamma ds = 0 (Eq 62); signature (+,-,-,-,-).
# Context: u4_definition_reconciliation.sage (F5), test_body_motion.sage,
#   gr_projection_audit.sage (G2), u4_rest_mass_nonconservation.sage (F5).
# Supplementary Materials: step [F5], Sec. S2 (code index: Sec. S6).
# ---------------------------------------------------------------------------
#
# u4_geometric_motion.sage                             [A] (assertion-hardened)
#
# REVIEW POINT.  Take u4 = dx4/ds literally (the GEOMETRIC reading): a fluid
# element that is compressed moves along x4 although its mass is conserved.
# Two questions:
#
#   [A] Does motion along x4 disturb the 4D motion?  On the product metric
#       (g_a4 = 0, g_44 = -1, g_ab independent of x4), for an ARBITRARY
#       history x4(s), the space-time components of the 5D equation leave a
#       residual proportional to u_a only: a reparametrisation.  The 4D path
#       is an exact geodesic of g_ab (checked on a Schwarzschild circular
#       orbit with x4(s) left free).  Any reaction that drives x4 is
#       orthogonal to the 5D velocity, so its space-time part is also along
#       u_a.
#   [B] Mutation: if g_tt depends on x4, the radial residual is not zero:
#       the 4D path is no longer a geodesic.  So it is d_4 g_ab = 0 (K = 0
#       with g_a4 = 0), not u4 = 0, that protects the 4D motion.  (A purely
#       conformal warp f(x4) g_ab only reparametrises too, so K = 0 is
#       sufficient here, not the only possibility.)
#   [C] The cost.  A 5D timelike path has |dx4/dt| < c, so on the geometric
#       reading the COMPRESSION RATE of matter is capped:
#         |D gamma/Dt| = gamma |div v| < a0 c ~ 1.1e16 kg m^-3 s^-1.
#       Order-of-magnitude rates: lab planar shocks ~1e-4 of the cap,
#       laser shocks ~0.1, cosmological expansion ~1e-50, stellar core
#       bounce ~1e5 OVER the cap.
#
#   [D] The weight of an element at rest on the mass axis: with gamma = a0 x4
#       u4 = 0 needs the reaction a0 (manuscript Sec. 3.4); with the
#       element's fixed mass m as the weight, u4 = 0 is a free solution.
#       The 4D path is a geodesic in both cases.
#
# RESULT: the GR test-body recovery does not need u4 = 0 under either reading.
# But the geometric reading turns the c-cap into a bound on compression rates
# that core collapse exceeds by ~1e5 (order of magnitude; with a0 as
# calibrated).  So x4 is taken to be the comoving mass coordinate
# (u4_definition_reconciliation.sage), for which no such bound arises; by [A]
# the choice cannot be seen in 4D.
#
# RUN: sage u4_geometric_motion.sage
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("u4_geometric_motion.sage -- motion along x4 and the compression cap")
print("=" * 74)

s = var('s')
M, r0, a0, eps = var('M r0 a0 eps')
assume(M > 0, r0 > 3 * M, a0 > 0)
t_, r_, th_, ph_, w_ = var('t_ r_ th_ ph_ w_')
X = [t_, r_, th_, ph_, w_]

def EL_residual(gmat, curve, weight=lambda q: a0 * q):
    """F_k = d/ds(gamma u_k) - (gamma/2) d_k g_ij u^i u^j - d_k gamma,
    gamma = weight(x4) (default a0 x4), on the given curve x^i(s)."""
    sub = dict(zip(X, curve))
    u = [diff(c, s) for c in curve]
    gam = weight(curve[4])
    F = []
    for k in range(5):
        uk = sum(gmat[k, l].subs(sub) * u[l] for l in range(5))
        dg = sum(diff(gmat[i, j], X[k]).subs(sub) * u[i] * u[j]
                 for i in range(5) for j in range(5))
        dgam = diff(weight(w_), w_).subs({w_: curve[4]}) if k == 4 else 0
        F.append(diff(gam * uk, s) - gam * dg / 2 - dgam)
    return F, u

def schw(warp):
    Fr = 1 - 2 * M / r_
    return diagonal_matrix(SR, [Fr * warp, -1 / Fr, -r_**2,
                                -r_**2 * sin(th_)**2, -1])

# circular geodesic: dphi/dt = sqrt(M/r0^3), 4D-normalised in sigma
A = 1 / sqrt(1 - 3 * M / r0)
Om = A * sqrt(M / r0**3)
w = function('w')(s)                       # x4(s): ARBITRARY history
sig = function('sig')(s)
curve = [A * sig, r0, pi / 2, Om * sig, w]

print("\n[A] Product metric, x4(s) free")
g = schw(1)
sub = dict(zip(X, curve))
u = [diff(c, s) for c in curve]
norm5 = sum(g[i, i].subs(sub) * u[i]**2 for i in range(5))
sigp = sqrt(1 + diff(w, s)**2)             # 5D arc length fixes sigma'
fix = {diff(sig, s): sigp}
check((norm5.subs(fix) - 1).simplify_full() == 0,
      "5D arc length: sigma' = sqrt(1 + (dx4/ds)^2)")
F, u = EL_residual(g, curve)
F = [f.subs({diff(sig, s, 2): diff(sigp, s)}).subs(fix).simplify_full() for f in F]
uf = [ui.subs(fix) for ui in u]
u_low = [(g[k, k].subs(sub) * uf[k]).simplify_full() for k in range(5)]
check(F[1] == 0 and F[2] == 0,
      "radial and polar residuals vanish for any x4(s)")
lam = (F[0] / u_low[0]).simplify_full()
check((F[3] - lam * u_low[3]).simplify_full() == 0,
      "space-time residual F_a = lambda(s) u_a: a reparametrisation only")
print("    lambda =", lam)
check(lam != 0 and lam.subs({w: s, diff(w, s): 1, diff(w, s, 2): 0}) != lam, "lambda is driven by the x4 history (non-trivial)")
orth = sum(F[k] * uf[k] for k in range(5)).simplify_full()
check(orth == 0, "the total reaction is orthogonal to the 5D velocity")
print("    so the 4D path is an exact geodesic of g_ab whatever x4 does")

print("\n[B] Mutation: g_tt depends on x4 (g_tt -> g_tt e^{2 eps x4})")
gw = schw(exp(2 * eps * w_))
Fw, _ = EL_residual(gw, curve)
Fr_w = Fw[1].subs(fix).simplify_full()
check(Fr_w.subs(eps=0).simplify_full() == 0 and Fr_w != 0,
      "radial residual != 0 when d_4 g_ab != 0: the 4D path is no geodesic")

print("\n[C] The c-cap on compression under the geometric reading")
a0n, cn, Gn = 3.65e7, 2.998e8, 6.674e-11
cap = a0n * cn
print("    cap a0 c = %.2e kg m^-3 s^-1" % cap)
rho_nuc = 3e17
t_bounce = 1 / sqrt(Gn * rho_nuc)
cases = [
    ("cosmological expansion, 3 H rho_crit", 3 * 2.3e-18 * 9e-27),
    ("planar shock, Al 2.7e3 -> 5e3 kg/m^3 in 1 ns", 2.3e3 / 1e-9),
    ("laser shock, 2.7e3 -> 6e3 kg/m^3 in 3 ps", 3.3e3 / 3e-12),
    ("core bounce, rho_nuc / (G rho_nuc)^-1/2", rho_nuc / t_bounce),
]
ratios = {}
for name, rate in cases:
    ratios[name] = rate / cap
    print("    %-48s %.1e x cap" % (name, rate / cap))
check(ratios[cases[0][0]] < 1e-40 and ratios[cases[1][0]] < 1e-3,
      "cosmology and ordinary shocks sit far below the cap")
check(ratios[cases[3][0]] > 1e4,
      "core bounce exceeds the cap by > 1e4 (order of magnitude)")
print("    With x4 the comoving mass coordinate (u4_definition_reconciliation)")
print("    every mass-conserving process has u4 = 0 and no cap applies to")
print("    compression.")

print("\n[D] The weight of an element, at rest on the mass axis (u4 = 0)")
w0, m = var('w0 m')
assume(w0 > 0, m > 0)
curve0 = [A * s, r0, pi / 2, Om * s, w0]
for name, wt in [("gamma = a0 x4 (weight a function of x4)", lambda q: a0 * q),
                 ("fixed element mass m", lambda q: m + 0 * q)]:
    F0, _ = EL_residual(g, curve0, wt)
    F0 = [f.simplify_full() for f in F0]
    check(all(F0[k] == 0 for k in range(4)),
          "%s: 4D path a geodesic" % name)
    print("    %s: x4 residual = %s" % (name, F0[4]))
    if "fixed" in name:
        check(F0[4] == 0, "fixed mass: u4 = 0 is a free solution, no reaction")
    else:
        check((F0[4] + a0).simplify_full() == 0,
              "gamma = a0 x4: u4 = 0 needs the reaction gamma f_4 = a0")

print("""
VERDICT
  Motion along x4 leaves the 4D path an exact geodesic on the product metric:
  the u4 terms are a reparametrisation, and a warp (d_4 g_ab != 0) is what
  would spoil it.  So the GR test-body recovery does not need u4 = 0.  Read
  geometrically (x4 the density field), u4 would cap compression at
  a0 c ~ 1e16 kg m^-3 s^-1, which a stellar core bounce exceeds by ~1e5
  (order of magnitude, a0 as calibrated).  We therefore take x4 to be the
  comoving mass coordinate (u4_definition_reconciliation), for which u4 is
  the continuity residual and the cap bounds only mass non-conservation;
  that choice is invisible in 4D because of [A].  At rest on the mass axis
  the 4D motion is a geodesic for either weight; u4 = 0 needs the reaction
  a0 if the weight is a0 x4, and none if it is the element's fixed mass.
""")
print("=" * 74)
print("ALL ASSERTIONS PASSED")
print("=" * 74)
