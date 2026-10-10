#!/usr/bin/env sage
# --------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, "The Dynamic Theory" (2011) Ch. 3,
#   Eqs (3.24)-(3.36), pp. 68-71 (time-dependent gravitational field,
#   rotation curves, the MOND scale).  Solar-system inputs (IMPORTED, not
#   re-derived): Mercury perihelion from MESSENGER ranging (Park et al.,
#   AJ 153 (2017) 121: beta-1 = (-2.7 +- 3.9)e-5); Gdot/G from lunar laser
#   ranging (Hofmann & Mueller, CQG 35 (2018) 035015: (7.1 +- 7.6)e-14 /yr --
#   a GR-framework reduction, NOT re-derived through DT's drifting units).
# Builds on dark_matter_rotation_curves.sage, G_variation_constraints.sage,
#   black_hole_ringdown.sage (strong-field fork), matter_source_tab.sage.
# Full derivation, context and provenance:
#   Supplementary_Material.pdf -- §S3, Appendix A.
# --------------------------------------------------------------------------

# ============================================================================
# gr_limit_domain.sage                             [A]  (assertion-hardened)
#
# GAP #5: THE DOMAIN OF DT'S GR (NEWTONIAN) LIMIT.
#
# Williams' weak-field mass sector (Eqs 3.24-3.27) is Newton's law times two
# factors:  (1 - lambda/r)  and  (1 - H0 tau).
#   * lambda = GM/c^2 is the strong-field scale.  Its domain question was
#     settled by black_hole_ringdown.sage: the data select DT's metric (GR)
#     reading, which has no lambda-departure at all.  Not repeated here.
#   * (1 - H0 tau) is the TIME-DEPENDENT field.  tau obeys
#     tau'' = -H0 lambda / r (Eq 3.25), so tau = tau0dot t + O(H0 lambda t^2/r),
#     and ONE number, tau0dot, sets every departure from Newton/GR.
#
# WHAT THE SCRIPT ESTABLISHES
#   (A) tau0dot = 0 gives Newton exactly (to O(H0^2 lambda t^2 / r) ~ 1e-24 per century).
#   (B) ONE DIAL.  The rotation-curve anomaly (Eq 3.30: kappa = H0 tau0dot/c)
#       and the fractional rate of change of the source strength,
#       d ln(G_eff M)/dt = -H0 tau0dot, are the SAME parameter:
#       kappa = |d ln G_eff/dt| / c.  The spatial anomaly IS the secular drift,
#       read across the light-travel time r/c (Williams' own step 3.29->3.30).
#   (C) Perihelion: a = a_N(1 + kappa r) precesses by -pi kappa r per orbit
#       (exact apsidal formula, cross-checked by direct numerical integration).
#   (D) Flat curves need kappa = 1/r_c, so tau0dot = c/(r_c H0).  For the
#       Milky Way this is ~5e5, and then the galaxy's gravitational strength
#       must fall by its own size every r_c/c ~ 3e4 yr -- and the linear
#       factor (1 - H0 tau) crosses ZERO that soon after.  Disks are ~1e10 yr
#       old.  The flat part (r >~ r_c) is also exactly where |H0 tau| >~ 1,
#       outside the regime the expansion (1 - H0 tau) was written for.
#   (E) Eq 3.33 drops v^2/c^2 as small, but on a circular orbit
#       v^2/c^2 = lambda/r0 EXACTLY, the very term kept.  Corrected:
#       tau0dot^2 = 2E0/(m c^2).  tau0dot is a free initial condition, and
#       tau0dot ~ 5e5 needs an initial tau-kinetic energy ~1e11 rest masses.
#   (F) Solar-system domain.  The universal-a0 reading (a0 the same for the
#       Sun) gives Mercury ~1.5e4 arcsec/century of extra precession and a
#       solar G e-folding in ~41 days: EXCLUDED.  Mercury alone allows
#       tau0dot_sun < ~1.4e4; LLR (imported bound) < ~1.4e-3.
#   (G) Domain map: Newton/GR holds where tau0dot H0 r/c and tau0dot H0 t are
#       below the test precision.  Solar system: tau0dot <~ 1e-3 (LLR,
#       imported) / 1.4e4 (Mercury, DT-native kinematics).  Galaxies (to be
#       flat): ~5e5.
#   (H) No compromise.  v_flat^2 x T_fade = GM/c for every tau0dot, so speed
#       is bought only with faster fading.  The natural reading (Eq 3.19,
#       tau0dot = 1) gives the Milky Way a ~0.25 km/s floor.  The shortfall
#       (= tau0dot needed) falls as M^-1/2: ~1e7 dwarfs, ~5e5 Milky Way, ~4e3
#       clusters; it closes only at ~the mass inside the Hubble radius (the
#       a0 ~ c H0 coincidence restated).  Any fixed tau0dot gives v^2 ~ M, not
#       Tully-Fisher: BTF itself needs the per-system tuning.
#
# VERDICT (honest bracket).  DT's GR limit is EXACT at tau0dot = 0, and its
# domain is controlled by that single initial condition (plus lambda/r, which
# the data already settle in favour of the metric reading).  The solar
# system bounds it.  But Williams' flat-rotation-curve mechanism, taken at
# face value, needs tau0dot ~ 5e5, i.e. each galaxy's gravity changing by
# order unity on its light-crossing time (~3e4 yr) and the linearised drift
# factor (1 - H0 tau) leaving its own validity range on the flat part.  This
# is an INTERNAL-CONSISTENCY FLAG on the dark-matter POSTDICTION (a reading
# of Ch. 3), NOT a falsification of DT, and NOT a negative for DT's GR
# sector, which survives untouched at tau0dot -> 0.  Defeater (why it is not
# a kill under the strict bar): Williams flags "two times" in 3.28-3.30 and
# sets tau0 = t0 = 0 "when the light left the star"; a reading in which the
# retardation term is not a secular decline of the source would escape it,
# but the source supplies no such reading, and on the face-value reading
# the surviving GR sector does NOT carry the dark-matter result.  Unlike
# Cassini or ringdown, this is a cost, not only a selection.
#
# RUN: /home/dsteele/miniconda3/envs/sage/bin/sage gr_limit_domain.sage
# ============================================================================

from sage.all import *


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


print("=" * 74)
print("gr_limit_domain.sage -- gap #5: where does DT reduce to Newton/GR?")
print("=" * 74)

G, M, r, H0, c, t, td, lam, kap = var('G M r H0 c t td lam kap', domain='positive')

# ---------------------------------------------------------------------------
print("\n[A] Eqs 3.24-3.27: Newton x (1 - H0 tau) x (1 - lambda/r)")
tau = td * t - H0 * lam / (2 * r) * t**2              # Eq 3.26, tau0 = t0 = 0
ok("A1: tau solves Eq 3.25, tau'' = -H0 lambda / r",
   bool((diff(tau, t, 2) + H0 * lam / r).simplify_full() == 0))
v2 = (1 - H0 * tau) * G * M / r                         # Eq 3.27, lambda << r
dev = (v2 / (G * M / r) - 1).subs(td=0).simplify_full()
print("    with tau0dot = 0 the fractional departure is %s" % dev)
ok("A2: tau0dot = 0 leaves only H0^2 lambda t^2 / (2r)", bool((dev - H0**2 * lam * t**2 / (2 * r)).simplify_full() == 0))
H0n, cn, GMs = 70e3 / 3.0857e22, 2.99792458e8, 1.32712440018e20
yr = 3.15576e7
rM = 5.7909e10
resid = float(H0n**2 * (GMs / cn**2) * (100 * yr)**2 / (2 * rM))
print("    size at Mercury over a century: %.2e" % resid)
ok("A3: that residual is negligible (< 1e-20) -- Newton/GR is exact at tau0dot = 0", resid < 1e-20)

# ---------------------------------------------------------------------------
print("\n[B] One dial: the rotation-curve anomaly IS the secular drift")
v2_330 = v2.subs(t=-r / c).series(lam, 1).truncate()   # Eq 3.29 -> 3.30
kappa_expr = ((v2_330 / (G * M / r) - 1) / r).simplify_full()
print("    Eq 3.30 anomaly per unit radius kappa = %s" % kappa_expr)
ok("B1: kappa = H0 tau0dot / c (Eq 3.30, lambda -> 0)", bool((kappa_expr - H0 * td / c).simplify_full() == 0))
dlnG = (diff(1 - H0 * tau, t) / (1 - H0 * tau)).subs(t=0).simplify_full()
print("    d ln(G_eff M)/dt at emission = %s" % dlnG)
ok("B2: d ln G_eff/dt = -H0 tau0dot, so kappa = |d ln G_eff/dt| / c", bool((dlnG + c * kappa_expr).simplify_full() == 0))

# ---------------------------------------------------------------------------
print("\n[C] Perihelion advance of a = a_N (1 + kappa r)")
x = var('x')
F = G * M / r**2 * (1 + kap * r)
n_eff = (r * diff(F, r) / F).simplify_full()            # d ln F / d ln r
Psi = pi / sqrt(3 + n_eff)                               # apsidal angle (near-circular)
dw = (2 * Psi - 2 * pi).subs(kap=x / r).series(x, 2).truncate().subs(x=kap * r)
print("    precession per orbit = %s" % dw.simplify_full())
ok("C1: delta omega = -pi kappa r per orbit (retrograde) to first order",
   bool((dw + pi * kap * r).simplify_full() == 0))
# cross-check: direct integration of an eccentric orbit, GM = r0 = 1
def apsidal_numeric(k, rp, ra):
    U = lambda s: -1.0 / s + k * log(s)
    L2 = 2 * (U(ra) - U(rp)) / (1 / rp**2 - 1 / ra**2)
    E = U(rp) + L2 / (2 * rp**2)
    m_, d_ = (ra + rp) / 2, (ra - rp) / 2
    def f(phi):
        s = m_ - d_ * cos(phi)
        rad = 2 * (E - U(s)) - L2 / s**2
        return sqrt(L2) / s**2 * d_ * sin(phi) / sqrt(max(rad, 1e-300))
    val, _ = numerical_integral(f, 0, float(pi))
    return 2 * val - 2 * float(pi)
k_test = 1e-3
num0 = apsidal_numeric(0.0, 0.98, 1.02)
num = apsidal_numeric(k_test, 0.98, 1.02)
print("    numeric: Kepler %.2e rad, kappa=1e-3 gives %.4e rad (formula %.4e)"
      % (num0, num, -float(pi) * k_test))
ok("C2: numerical orbit reproduces -pi kappa r to 2% (Kepler control closes to 1e-6)",
   abs(num0) < 1e-6 and abs(num / (-float(pi) * k_test) - 1) < 0.02)

# ---------------------------------------------------------------------------
print("\n[D] What flat rotation curves cost: the Milky Way")
a0n = 1.2e-10
Msun = 1.98892e30
GMmw = 6.674e-11 * 6e10 * Msun
rc_mw = sqrt(GMmw / a0n).n()
td_mw = (cn / (rc_mw * H0n)).n()
tef_mw = (rc_mw / cn / yr).n()
print("    r_c = %.3g m (%.1f kpc), tau0dot = c/(r_c H0) = %.3g" % (rc_mw, rc_mw / 3.0857e19, td_mw))
print("    G_eff M changes by its own size every r_c/c = %.3g yr" % tef_mw)
print("    (1 - H0 tau) reaches zero %.3g yr after the tau = 0 epoch" % tef_mw)
ok("D1: flat curves need tau0dot >> 1 for a galaxy (~5e5)", 1e5 < td_mw < 1e6)
ok("D2: the galaxy's gravity must change by O(1) within ~3e4 yr, vs disk ages ~1e10 yr",
   tef_mw < 1e5 and 1e10 / tef_mw > 1e5)
ok("D3: on the flat part (r = 3 r_c) |H0 tau| = r/r_c = 3, outside |H0 tau| << 1",
   abs(float(H0n * td_mw * 3 * rc_mw / cn) - 3) < 1e-9)

# ---------------------------------------------------------------------------
print("\n[E] Eq 3.32-3.33: what sets tau0dot")
E0, m, r0 = var('E0 m r0', domain='positive')
v2circ = G * M / r0                                      # circular, lambda << r0
eq332 = 2 * E0 / (m * c**2) - (x**2 + v2circ / c**2 - (G * M / c**2) / r0)
sol = solve(eq332 == 0, x)
tdsq = [s.rhs()**2 for s in sol][0].simplify_full()
print("    Eq 3.32 on a circular orbit: tau0dot^2 = %s" % tdsq)
ok("E1: v^2/c^2 = lambda/r0 exactly, so the dropped term cancels the kept one: tau0dot^2 = 2E0/(m c^2)",
   bool((tdsq - 2 * E0 / (m * c**2)).simplify_full() == 0))
print("    -> tau0dot is a FREE initial condition (the initial tau-kinetic energy)")
print("    -> the Milky Way value needs 2E0/(m c^2) = %.2g rest energies" % float(td_mw**2))
ok("E2: the galactic value needs an initial tau-energy > 1e10 x the rest energy", td_mw**2 > 1e10)

# ---------------------------------------------------------------------------
print("\n[F] Solar-system domain")
rc_sun = sqrt(GMs / a0n).n()
kap_univ = 1 / rc_sun
orb_cy = 36525 / 87.969
as_rad = 180 * 3600 / float(pi)
dw_univ = float(pi * kap_univ * rM * orb_cy * as_rad)
print("    universal-a0 reading: r_c(Sun) = %.3g m, Mercury extra precession %.3g arcsec/cy" % (rc_sun, dw_univ))
print("    and solar G_eff e-folds every r_c/c = %.1f days" % float(rc_sun / cn / 86400))
bound_park = 3 * 42.98 * 3.9e-5 / 3                      # 3 sigma of (2+2g-b)/3 * 42.98"
ok("F1: universal-a0 reading is excluded by Mercury (> 1e6 x the MESSENGER bound)",
   dw_univ / bound_park > 1e6)
kap_max = bound_park / orb_cy / as_rad / (float(pi) * rM)
td_merc = kap_max * cn / H0n
gdot_llr = 1e-13
td_llr = gdot_llr / (H0n * yr)
print("    Mercury (|d omega| < %.4f\"/cy, imported): tau0dot_sun < %.2g" % (bound_park, td_merc))
print("    LLR |Gdot/G| < 1e-13/yr -> tau0dot < %.2g   [IMPORTED bound: a GR-framework"
      % td_llr)
print("      reduction, not re-derived through DT's drifting units]")
ok("F2: both solar bounds sit below the galactic requirement (Mercury x%.0f, LLR x%.1e)"
   % (td_mw / td_merc, td_mw / td_llr), td_merc < td_mw and td_llr < td_mw)

# ---------------------------------------------------------------------------
print("\n[G] Domain map of the Newton/GR limit (one dial, tau0dot)")
print("    %-34s %s" % ("regime", "tau0dot"))
print("    %-34s %s" % ("exact Newton/GR", "0"))
print("    %-34s < %.1g  (DT-native kinematics)" % ("solar system, Mercury perihelion", td_merc))
print("    %-34s < %.1g  (imported bound)" % ("solar system, LLR Gdot/G", td_llr))
print("    %-34s ~ %.1g  (Williams Eqs 3.35-3.36)" % ("galaxy, flat rotation curve", td_mw))
ok("G1: the GR limit and the MOND regime sit at opposite ends of the same dial",
   td_llr < 1 < td_merc < td_mw)

# ---------------------------------------------------------------------------
print("\n[H] No compromise: the trade-off rule and the mass dependence")
T_fade = 1 / (H0 * td)                                   # time for G_eff M to change by itself
v2flat = G * M * (H0 * td / c)                           # flat floor, Eq 3.30
ok("H1: v_flat^2 x T_fade = GM/c for EVERY tau0dot (and independent of H0)",
   bool((v2flat * T_fade - G * M / c).simplify_full() == 0))
v_nat = sqrt(GMmw * H0n / cn)                            # Eq 3.19: tau0dot = 1
print("    natural reading (Eq 3.19, tau0dot = 1): Milky Way v_flat = %.3f km/s (observed ~200)" % (v_nat / 1e3))
ok("H2: the natural reading gives a Milky Way floor of ~0.25 km/s", 0.1e3 < v_nat < 0.5e3)
print("    %-24s %10s %14s %16s" % ("system", "M [Msun]", "tau0dot needed", "r_c/c [yr]"))
need = []
for name, Mx in [("dwarf galaxy", 1e8), ("Milky Way", 6e10), ("galaxy cluster", 1e15)]:
    GMx = 6.674e-11 * Mx * Msun
    rcx = sqrt(GMx / a0n)
    need.append(float(cn / (H0n * rcx)))
    print("    %-24s %10.0e %14.2g %16.2g" % (name, Mx, need[-1], float(rcx / cn / yr)))
ok("H3: every real system falls short, the smallest worst (shortfall falls as M^-1/2, > 1e3 even for clusters)",
   need[0] > need[1] > need[2] > 1e3)
M_match = a0n * (cn / H0n)**2 / 6.674e-11 / Msun         # r_c = c/H0  <=>  tau0dot = 1
M_hub = (4 * float(pi) / 3) * (3 * H0n**2 / (8 * float(pi) * 6.674e-11)) * (cn / H0n)**3 / Msun
print("    the shortfall closes (tau0dot = 1) only at M = %.2g Msun" % M_match)
print("    mass inside the Hubble radius (critical density) = %.2g Msun" % M_hub)
ok("H4: that closing mass is of order the mass inside the Hubble radius (a0 ~ c H0 restated)",
   0.1 < M_match / M_hub < 10)
p = var('p')
ok("H5: with any FIXED tau0dot, v_flat^2 ~ M^1, not the Tully-Fisher M^(1/2): BTF needs tau0dot ~ M^(-1/2)",
   bool(diff(log(v2flat), M) * M == 1)
   and bool((diff(log(sqrt(G * M * a0n)), M) * M - QQ(1) / 2).simplify_full() == 0))

# ---------------------------------------------------------------------------
print("\n[VERDICT]")
print("  * DT's Newton/GR limit is exact at tau0dot = 0; one initial condition")
print("    controls its domain (lambda/r is already settled by the metric reading).")
print("  * The galactic anomaly and the secular drift are the SAME parameter.")
print("  * Flat curves, taken at face value, need each galaxy's gravity to change by")
print("    O(1) every ~3e4 yr: an internal-consistency FLAG on the Ch. 3 dark-matter")
print("    postdiction, not a falsification of DT and not a negative for its GR")
print("    sector.  Defeater: Williams' 'two times' bookkeeping (Eqs 3.28-3.30).")
print("  * No compromise: v_flat^2 x T_fade = GM/c.  The natural reading gives the")
print("    right SHAPE at 1e-3..1e-7 of the size (worst for dwarfs); it is the")
print("    right size only for the universe as a whole.")
print("\n" + "=" * 74)
print("gr_limit_domain COMPLETE.")
print("=" * 74)
