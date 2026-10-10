#!/usr/bin/env sage
# --------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, "The Dynamic Theory" (2011), Ch.13 "Black
#   Holes" (p.452) and Ch.10 (perihelion, Eq 10.36): the non-singular
#   gravitational potential V_g = -(GM/r) e^{-lambda/r}, lambda = GM/c^2, as
#   used by nonsingular_black_holes.sage.  The vacuum surface equations
#   G_mu nu = 0 are DT's recovered GR sector (gr_on_4d_surface,
#   kk_radion_g44, gravitational_radiation).
# Data: Bertotti, Iess & Tortora, Nature 425 (2003) 374 (gamma - 1 =
#   (2.1 +/- 2.3)e-5); planetary ephemerides |beta - 1| < ~1e-4 (Mercury
#   perihelion); LIGO/Virgo, PRL 116 (2016) 061102 and PRL 116 (2016) 221101
#   (GW150914: remnant 62 Msun source-frame, z = 0.09, spin 0.67; ringdown
#   consistent with the Kerr 220 mode); Berti, Cardoso & Will, PRD 73 (2006)
#   064030 (Kerr QNM fit); Event Horizon Telescope, ApJL 875 (2019) L1 and
#   ApJL 930 (2022) L12 (M87*, Sgr A* rings consistent with Kerr shadows).
# Full derivation, context and provenance:
#   Supplementary_Material.pdf -- §S3, Appendix A.
# --------------------------------------------------------------------------

# ============================================================================
# black_hole_ringdown.sage                         [A]  (assertion-hardened)
#
# THE RINGDOWN TEST OF DT'S NON-SINGULAR BLACK HOLE.
#
# After two black holes merge, the remnant "rings" at frequencies set by its
# LIGHT RING, the radius where light can orbit (r = 3GM/c^2 for
# Schwarzschild).  The real part of each mode is the orbital frequency of
# that ring, and the damping rate is how fast light leaks off it (its
# Lyapunov exponent).  So the question for DT is simple: does Williams'
# non-singular potential produce a light ring, and if so where?
#
# Williams gives only the Newtonian potential, so the metric must be closed.
# The script promotes V_g to g_00 = 1 + 2V_g/c^2 and lets the solar system
# choose the closure: PPN gamma = beta = 1 (Cassini, Mercury).
#
# WHAT THE SCRIPT ESTABLISHES
#   (A) Controls: Schwarzschild has a light ring at 3M, and the eikonal
#       ringdown Omega_c = lambda_L = 1/(3 sqrt3 M); the Kerr fit gives
#       f_220 ~ 250 Hz for GW150914's remnant, as observed.
#   (B) Closure: read with an AREAL radius, Williams' g_00 has beta = 2,
#       so the perihelion advance is 2/3 of GR's -- excluded by Mercury.  Read
#       with an ISOTROPIC radius it has beta = 1 (the e^{-lambda/r} factor
#       supplies exactly GR's second-order term), and gamma = 1 fixes the
#       spatial part to first order.  Three such closures are kept.
#   (C) No horizon: g_00 >= 1 - 2/e = 0.264 everywhere (closure-independent).
#   (D) No light ring: in every solar-system-allowed closure the orbital
#       frequency of light, Omega(r), is monotonic -- there is no photon
#       orbit at any radius.  The potential is too shallow (|V_g| <= c^2/e)
#       to trap light.  (Even the excluded areal closure has none.)
#   (E) Not a vacuum solution: the closed metric has G^t_t ~ -(lambda/r)^4
#       outside the mass, i.e. an exterior source of NEGATIVE energy density
#       (weak energy condition violated; sign fixed by a control; computed
#       for the closure B = 1/g_00 only).  The vacuum surface equations
#       that DT recovers force Schwarzschild instead (Birkhoff).
#
# VERDICT (honest bracket).  DT contains a FORK at strong field.
#   * GR reading (the vacuum surface equations recovered in gaps #1-#3):
#     Schwarzschild/Kerr, horizon, light ring at 3M, the Kerr ringdown --
#     consistent with GW150914 and the EHT rings.
#   * Williams' non-singular reading (e^{-lambda/r} promoted to the metric):
#     no horizon AND no light ring, in every closure the solar system allows.
#     With no light ring there are no light-ring modes, so the remnant
#     cannot ring like Kerr, and there is no photon ring for EHT to see.  The
#     ringdown it would show is set by an interior DT does not specify, so
#     it cannot even be predicted -- but it is not the observed Kerr one.
#   Ringdown and EHT therefore select the GR reading, as Cassini selected
#   the restricted reading in kk_radion_g44.  This is a strong NEGATIVE for
#   the non-singular-black-hole claim of nonsingular_black_holes.sage, but
#   not a falsification of DT, whose own vacuum sector predicts the observed
#   ringdown.  Scope: static, non-rotating, eikonal (geometric-optics)
#   analysis; the metric closure is ours (Williams gives only V_g); the
#   2PN-and-beyond spatial part is not fixed by the solar system.
#
# RUN: /home/dsteele/miniconda3/envs/sage/bin/sage black_hole_ringdown.sage
# ============================================================================

from sage.all import *


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


print("=" * 74)
print("black_hole_ringdown.sage -- does DT's non-singular black hole ring like Kerr?")
print("=" * 74)

r = var('r', domain='positive')       # radius in units of lambda = GM/c^2
u = var('u')                          # u = lambda / r


def light_rings(Om2, lo=0.02, hi=60, N=6000):
    """radii where Omega^2(r) = f/R^2 has an extremum (sign change of d/dr)"""
    d = fast_callable(diff(Om2, r), vars=[r], domain=RDF)
    xs = [lo + (hi - lo) * k / N for k in range(N + 1)]
    ds = [d(x) for x in xs]
    return [xs[k] for k in range(N) if ds[k] * ds[k + 1] < 0]


def eikonal(f, g, R, r0):
    """Omega_c and Lyapunov exponent at the light ring r0, metric
       ds^2 = f dt^2 - dr^2/g - R^2 dOmega^2 (tortoise d/dr* = sqrt(f g) d/dr)"""
    V = f / R**2
    D = lambda E: sqrt(f * g) * diff(E, r)
    V2 = D(D(V))
    Om = sqrt(V.subs(r=r0))
    lam = sqrt(-V2.subs(r=r0) / (2 * V.subs(r=r0)))
    return Om.n(), lam.n()


# ---------------------------------------------------------------------------
print("\n[A] Controls: Schwarzschild and the GW150914 remnant")
fS = 1 - 2 / r
lrS = light_rings(fS / r**2)
print("    Schwarzschild light ring(s) at r/M = %s" % [round(x, 3) for x in lrS])
ok("A1: Schwarzschild has one light ring, at r = 3M", len(lrS) == 1 and abs(lrS[0] - 3) < 0.02)
OmS, lamS = eikonal(fS, fS, r, 3)
print("    eikonal: Omega_c M = %.5f, lambda_L M = %.5f  (1/(3 sqrt3) = %.5f)"
      % (OmS, lamS, 1 / (3 * sqrt(3.0))))
ok("A2: eikonal Omega_c = lambda_L = 1/(3 sqrt3 M) (Cardoso et al. 2009)",
   abs(OmS - 1 / (3 * sqrt(3.0))) < 1e-9 and abs(lamS - 1 / (3 * sqrt(3.0))) < 1e-9)
Tsun = 4.925490947e-6
Mdet, a = 62.0 * 1.09, 0.67
f220 = (1.5251 - 1.1568 * (1 - a)**0.1292) / (2 * pi.n() * Mdet * Tsun)
print("    Kerr f_220 for GW150914 (M_det = %.1f Msun, a = %.2f): %.0f Hz" % (Mdet, a, f220))
ok("A3: GR reading predicts the observed ~250 Hz ringdown", 235 < f220 < 265)

# ---------------------------------------------------------------------------
print("\n[B] Closing Williams' potential into a metric: what the solar system allows")
gW = 1 - 2 * u * exp(-u)                          # g_00 = 1 + 2 V_g / c^2
s00 = gW.series(u, 4).truncate()
print("    g_00 = %s + O(u^4)" % s00)
c2 = s00.coefficient(u, 2)
# areal radius, g_rr = 1/g_00 (Schwarzschild-like): g_00 = 1 - 2u + 2(beta-gamma)u^2
gam_areal = 1
beta_areal = c2 / 2 + gam_areal
perih = (2 + 2 * gam_areal - beta_areal) / 3
print("    AREAL reading:     gamma = %s, beta = %s, perihelion = %s x GR" % (gam_areal, beta_areal, perih))
ok("B1: areal reading has beta = 2 -> perihelion 2/3 of GR (excluded by Mercury)",
   beta_areal == 2 and perih == QQ(2) / 3)
beta_iso = c2 / 2                                  # isotropic: g_00 = 1 - 2u + 2 beta u^2
print("    ISOTROPIC reading: beta = %s  (GR's u^2 term supplied by e^{-u})" % beta_iso)
ok("B2: isotropic reading has beta = 1", beta_iso == 1)
gS_iso = ((1 - u / 2) / (1 + u / 2))**2
dev3 = (s00.coefficient(u, 3) - gS_iso.series(u, 4).truncate().coefficient(u, 3))
print("    first departure from Schwarzschild: u^3 coefficient differs by %s (2PN)" % dev3)
ok("B3: Williams' g_00 departs from GR only at 2PN (u^3)", dev3 != 0)
closures = {
    "B = 1/g00":          1 / gW,
    "B = 1 + 2u e^-u":    1 + 2 * u * exp(-u),
    "B = exp(2u e^-u)":   exp(2 * u * exp(-u)),
}
for name, B in closures.items():
    gam = B.series(u, 2).truncate().coefficient(u, 1) / 2
    ok("B4: closure %-18s has PPN gamma = 1 (Cassini)" % name, gam == 1)

# ---------------------------------------------------------------------------
print("\n[C] Horizon")
umin = solve(diff(gW, u) == 0, u)[0].rhs()
gmin = gW.subs(u=umin)
print("    min g_00 = 1 - 2/e = %.4f at r = lambda" % gmin.n())
ok("C1: g_00 > 0 everywhere -> no surface of infinite redshift (no horizon)",
   umin == 1 and gmin.n() > 0.26)

# ---------------------------------------------------------------------------
print("\n[D] Light ring (the thing that rings)")
A = gW.subs(u=1 / r)
for name, B in closures.items():
    Bt = B.subs(u=1 / r)
    lr = light_rings(A / (Bt * r**2))
    print("    isotropic, %-18s light rings: %s" % (name, lr if lr else "none"))
    ok("D1: no light ring for closure %s" % name, len(lr) == 0)
lrA = light_rings(A / r**2)
print("    areal (excluded anyway)          light rings: %s" % (lrA if lrA else "none"))
ok("D2: even the areal reading has no light ring", len(lrA) == 0)
# analytic: isotropic B = 1/A light-ring condition 2u(2-u)e^{-u} = 1 has no root
hmax = find_local_maximum(2 * u * (2 - u) * exp(-u), 0, 2)[0]
print("    isotropic B=1/A: max of 2u(2-u)e^-u = %.4f < 1, so no root" % hmax)
ok("D3: analytic check -- the light-ring equation has no solution", hmax < 1)
print("    (Schwarzschild in the same isotropic gauge has its ring at r = %.3f M)"
      % ((2 + sqrt(3.0)) / 2))
lrSi = light_rings(gS_iso.subs(u=1 / r) / ((1 + 1 / (2 * r))**4 * r**2), lo=0.55)  # outside the horizon rho = M/2
ok("D4: control -- the same finder locates Schwarzschild's isotropic ring",
   len(lrSi) == 1 and abs(lrSi[0] - (2 + sqrt(3.0)) / 2) < 0.02)

# ---------------------------------------------------------------------------
print("\n[E] Is the non-singular metric a vacuum solution of the surface equations?")
t, th, ph = var('t th ph')
X = [t, r, th, ph]
n = 4


def Gtt_of(g):
    gi = g.inverse()
    Gam = [[[sum(gi[a, d] * (diff(g[d, b], X[c]) + diff(g[d, c], X[b]) - diff(g[b, c], X[d]))
                 for d in range(n)) / 2 for c in range(n)] for b in range(n)] for a in range(n)]
    Ric = matrix(SR, n, n, lambda b, d: sum(diff(Gam[a][d][b], X[a]) - diff(Gam[a][a][b], X[d])
                                            + sum(Gam[a][a][f] * Gam[f][d][b] - Gam[a][d][f] * Gam[f][a][b]
                                                  for f in range(n)) for a in range(n)))
    Rs = sum(gi[a, b] * Ric[a, b] for a in range(n) for b in range(n))
    return (gi[0, 0] * Ric[0, 0] - Rs / 2).simplify_full()      # G^t_t = 8 pi G rho


# sign control: areal metric with mass function m(r) = k r^3 has rho = 3k/(8 pi G) > 0 for k > 0
k = var('k')
Gc = Gtt_of(diagonal_matrix(SR, [1, -1 / (1 - 2 * k * r**2), -r**2, -r**2 * sin(th)**2]))
ok("E0: sign control -- G^t_t = +6k for a positive-density ball (k > 0)",
   bool((Gc - 6 * k).simplify_full() == 0))
Bt = (1 / gW).subs(u=1 / r)
Gtt = Gtt_of(diagonal_matrix(SR, [A, -Bt, -Bt * r**2, -Bt * r**2 * sin(th)**2]))
lead = (Gtt.subs(r=1 / u)).series(u, 6).truncate()
print("    G^t_t (closure B = 1/g00) = %s + ...   (u = lambda/r)" % lead)
ok("E1: G^t_t != 0 outside the mass -> needs an exterior source (not vacuum)",
   not bool(lead == 0))
ok("E2: the leading exterior density is NEGATIVE (-u^4): violates the weak energy condition",
   lead.coefficient(u, 4) < 0)
print("    (vacuum surface equations + spherical symmetry -> Schwarzschild, Birkhoff)")

# ---------------------------------------------------------------------------
print("\n[VERDICT]")
print("  * GR reading of DT (vacuum surface equations): Schwarzschild/Kerr, light")
print("    ring at 3M, f_220 ~ %.0f Hz for GW150914 -- matches LIGO and EHT." % f220)
print("  * Williams' non-singular reading: no horizon (g_00 >= 0.264) and NO light")
print("    ring in any closure the solar system allows -> no Kerr ringdown, no")
print("    photon ring.  It also needs an exterior source of NEGATIVE energy")
print("    density, so it is not a solution of DT's own vacuum surface equations.")
print("  * Ringdown + EHT select the GR reading: a strong negative for the")
print("    non-singular black hole, not for DT.")
print("\n" + "=" * 74)
print("black_hole_ringdown COMPLETE.")
print("=" * 74)
