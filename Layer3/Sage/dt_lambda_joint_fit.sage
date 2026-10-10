#!/usr/bin/env sage
# --------------------------------------------------------------------------
# SOURCE: model = P. E. Williams, "The Dynamic Theory" Ch. 3 (DT redshift
#   relation Eq 3.6 + the evolving Chandrasekhar candle, Eqs 3.52-3.55), as fitted
#   in pantheon_dark_energy_fit.sage / pantheon_systematic_covariance.sage, PLUS a
#   cosmological-constant term (Report Eq 67 on Weyl's constant reading of lambda,
#   |Lambda| = lambda/4 != 0: weyl_lambda_reading.sage).
#   DATA = Pantheon binned SN-Ia (Scolnic et al. 2018, ApJ 859:101):
#   lcparam_DS17f.txt (40 bins) + sys_DS17f.txt (40x40 systematic covariance),
#   github.com/dscolnic/Pantheon; local copy pantheon_sys_DS17f.txt.
# --------------------------------------------------------------------------

# ============================================================================
# dt_lambda_joint_fit.sage                         [A]  (assertion-hardened)
#
# QUESTION. DT's secular drift (log-redshift distance + evolving candle, rate xi)
# already absorbs ~98% of the SN dark-energy signal. Weyl's procedure requires
# lambda != 0, so Lambda is nonzero -- but how LARGE may it be once the drift
# does the dark-energy work? Fit drift (xi) and Lambda JOINTLY and report the
# allowed range of Lambda.
#
# MODEL (DT is non-expanding, so Lambda enters through its STATIC metric, not a
# Friedmann expansion). Weak-field static patch of G + Lambda g = kappa T:
#     g_tt = 1 - Lambda r^2 / 3   ->   1 + z_Lambda = (1 - Lambda r^2/3)^(-1/2)
# (a static source at areal distance r, observer at r = 0). With r = (c/H0) x,
# where x = ln(1 + z_DT) is DT's own distance variable (Eq 3.6), and
#     L := Lambda c^2 / (6 H0^2)          (dimensionless; Lambda = 6 L H0^2/c^2)
# the observed redshift composes multiplicatively:
#     ln(1 + z) = x - (1/2) ln(1 - 2 L x^2)
# and the magnitude shape is DT's, evaluated at the solved x:
#     shape(z; xi, L) = 5 log10[(1+z) x] + (15/4) xi x / ln 10
# At L = 0 this is exactly the shape of pantheon_dark_energy_fit (checked).
#
# RESULT. Drift and Lambda are DEGENERATE (valley xi ~ 0.23 + 1.5 L). Full
# covariance leans mildly to Lambda > 0, ~Planck size (Lambda = 0: Delta chi2 6.9),
# but the lean rests on the two calibration-limited high-z bins (dropped: 3.1) and
# the best fit sits at the ansatz edge. "Lambda small" follows only once the drift
# rate is fixed independently (xi ~ 0.23); Gdot/G would fix it.
#
# NOTE ON "Lambda = 0". Weyl's constant reading REQUIRES lambda != 0 (else no
# Einstein term), so Lambda exactly 0 is not a DT option. Below, "Lambda = 0" is
# the reference point for the LIMIT Lambda -> 0: any |Lambda| << 1e-53 m^-2 gives
# a fit indistinguishable from it, so it stands for "negligible Lambda".
#
# HONEST BRACKET. The Lambda coupling is a LEADING-ORDER static-patch ansatz
# (gravitational redshift of the Lambda term composed with DT's drift redshift);
# Williams gives no Lambda term in his redshift law. Lambda's lensing / areal-
# distance corrections are the same order and are NOT included, so the bound is
# order-of-magnitude, not precision cosmology. The SIGN of Lambda is scanned both
# ways. The offset (absolute magnitude + H0) is marginalised, so the fit is
# H0-independent in L; converting to Lambda in m^-2 uses H0 = 70 km/s/Mpc.
#
# RUN: /home/dsteele/miniconda3/envs/sage/bin/sage dt_lambda_joint_fit.sage
# ============================================================================

from sage.all import *
import os


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


LN10 = float(log(10))
DATA = [
 (0.0140,14.57001926,0.0311),(0.0194,15.2279603507,0.02645),(0.0264,15.934102317,0.0231),
 (0.0329,16.4217829558,0.0235),(0.0396,16.836418956,0.02845),(0.0475,17.2346439906,0.0334),
 (0.0560,17.5758459622,0.03215),(0.0640,17.9104543501,0.0439),(0.0721,18.1780074589,0.03865),
 (0.0811,18.5543855764,0.04225),(0.0889,18.7006554037,0.03635),(0.1001,19.0388242428,0.03515),
 (0.1071,19.1720185353,0.03355),(0.1195,19.3711295152,0.0260),(0.1278,19.5555026238,0.02635),
 (0.1396,19.8185930766,0.02335),(0.1519,19.9445025196,0.0237),(0.1635,20.0642537602,0.02585),
 (0.1778,20.3389008923,0.01965),(0.1906,20.5047456137,0.02205),(0.2067,20.6854007995,0.0214),
 (0.2216,20.8610322507,0.0233),(0.2405,21.0518440051,0.0222),(0.2558,21.2037930719,0.02025),
 (0.2762,21.3636137766,0.0221),(0.2972,21.5790014334,0.0216),(0.3215,21.813277279,0.0214),
 (0.3453,21.9665751024,0.0234),(0.3708,22.1428940056,0.02145),(0.4049,22.3799192408,0.03225),
 (0.4355,22.5579272526,0.0254),(0.4738,22.797376642,0.02935),(0.5174,23.0011594448,0.02685),
 (0.5742,23.3004345026,0.0245),(0.6299,23.5036280324,0.031),(0.7240,23.8666293039,0.027),
 (0.8210,24.2445519795,0.0248),(0.9511,24.6411323707,0.0276),(1.2336,25.3039940191,0.05635),
 (1.6123,25.9259729107,0.0735),
]
NALL = len(DATA)

# --- systematic covariance (same file as pantheon_systematic_covariance) ----
HERE = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "."
vals = None
for p in [os.path.join(HERE, "pantheon_sys_DS17f.txt"), "Sage/pantheon_sys_DS17f.txt", "Layer3/Sage/pantheon_sys_DS17f.txt",
          "pantheon_sys_DS17f.txt"]:
    try:
        toks = open(p).read().split()
        dim = int(float(toks[0]))
        vals = [float(t) for t in toks[1:1 + dim * dim]]
        assert dim == NALL and len(vals) == NALL * NALL
        break
    except Exception:
        vals = None
assert vals is not None, "could not load pantheon_sys_DS17f.txt"
CSYS = matrix(RDF, NALL, NALL, vals)


class Fit:
    """chi2 for a subset of bins; cov = 'full' (sys + stat) or 'stat'."""

    def __init__(self, idx, cov):
        self.z = [DATA[i][0] for i in idx]
        self.mb = vector(RDF, [DATA[i][1] for i in idx])
        st = [DATA[i][2]**2 for i in idx]
        C = diagonal_matrix(RDF, st)
        if cov == "full":
            C = C + CSYS.matrix_from_rows_and_columns(idx, idx)
        self.Ci = C.inverse()
        self.one = vector(RDF, [1.0] * len(idx))
        self.A = self.one * (self.Ci * self.one)

    def chi2(self, shape):
        r0 = self.mb - vector(RDF, shape)
        Cr = self.Ci * r0
        return float(r0 * Cr - (self.one * Cr)**2 / self.A)


def solve_x(z, L):
    """Solve ln(1+z) = x - 0.5 ln(1 - 2 L x^2) for DT's distance variable x."""
    target = float(log(1 + z))
    x = target
    for _ in range(60):
        g = 1 - 2 * L * x * x
        if g <= 0:
            return None                       # beyond the static-patch horizon
        f = x - 0.5 * log(g) - target
        fp = 1 + 2 * L * x / g
        dx = f / fp
        x -= dx
        if abs(dx) < 1e-14:
            break
    return float(x)


def shape(zs, xi, L):
    out = []
    for z in zs:
        x = solve_x(z, L)
        if x is None or x <= 0:
            return None
        out.append(5.0 * float(log((1 + z) * x)) / LN10 + 3.75 * xi * x / LN10)
    return out


LS = [-1.0 + 0.002 * i for i in range(0, 751)]    # -1.00 .. +0.50 (step 0.002)


import numpy as np


def profile(fit):
    """min over xi of chi2, for each L; returns list (L, chi2min, xi_at_min).
    chi2 is quadratic in xi (shape is linear in xi), so the xi minimum is exact:
    shape = s0 + xi*u  ->  chi2(xi) = a - 2 b xi + d xi^2 (offset marginalised)."""
    Ci = np.array(fit.Ci); one = np.ones(len(fit.z)); A = one @ Ci @ one
    P = Ci - np.outer(Ci @ one, one @ Ci) / A          # offset-marginalised metric
    mb = np.array(fit.mb)
    prof = []
    for L in LS:
        xs = [solve_x(z, L) for z in fit.z]
        if any(x is None or x <= 0 for x in xs):
            continue
        xs = np.array(xs); zz = np.array(fit.z)
        s0 = 5.0 * np.log((1 + zz) * xs) / LN10
        u = 3.75 * xs / LN10
        r = mb - s0
        b = u @ P @ r; d = u @ P @ u
        xi = max(b / d, 0.0)          # DT: G weakens with time, so xi >= 0
        c2 = float(r @ P @ r - 2 * xi * b + xi * xi * d)
        prof.append((L, c2, float(xi)))
    return prof


def interval(prof, dchi2):
    cmin = min(p[1] for p in prof)
    inside = [p[0] for p in prof if p[1] - cmin <= dchi2]
    return min(inside), max(inside)


H0 = 70.0e3 / 3.0857e22          # s^-1
c = 2.99792458e8
UNIT = 6 * H0**2 / c**2          # Lambda = UNIT * L   [m^-2]
LAM_OBS = 1.1e-52                # Planck 2018 LambdaCDM value [m^-2]

print("=" * 72)
print("JOINT FIT: DT drift (xi) + cosmological constant Lambda -- Pantheon SN-Ia")
print("=" * 72)

# [0] machinery checks ------------------------------------------------------
print("\n[0] Machinery")
ok("x solver: L = 0 returns x = ln(1+z) exactly", abs(solve_x(1.0, 0.0) - float(log(2))) < 1e-14)
xL = solve_x(1.0, 0.3)
ok("x solver: residual of ln(1+z) = x - 0.5 ln(1-2Lx^2) < 1e-12 at L=0.3",
   abs(xL - 0.5 * log(1 - 0.6 * xL**2) - log(2.0)) < 1e-12)
ok("Lambda > 0 adds redshift (x < ln(1+z)); Lambda < 0 removes it",
   solve_x(1.0, 0.3) < log(2.0) < solve_x(1.0, -0.3))
full = Fit(list(range(NALL)), "full")
stat = Fit(list(range(NALL)), "stat")
ok("L = 0 reproduces pantheon_dark_energy_fit (stat-only chi2 79.8 at xi=0.235)",
   abs(stat.chi2(shape(stat.z, 0.235, 0.0)) - 79.8) < 0.1)
ok("L = 0 reproduces pantheon_systematic_covariance (full-cov chi2 47.1 at xi=0.230)",
   abs(full.chi2(shape(full.z, 0.230, 0.0)) - 47.1) < 0.1)
print("    unit: Lambda = %.3e * L  m^-2   (H0 = 70 km/s/Mpc)" % UNIT)
print("    Planck-2018 Lambda = %.2e m^-2  <->  L = %.3f" % (LAM_OBS, LAM_OBS / UNIT))
print("    NB 'Lambda = 0' below = the limit Lambda -> 0 (negligible Lambda);")
print("       exactly 0 is excluded by Weyl's constant reading (lambda != 0).")

# [1]-[3] profiles (xi >= 0 profiled out exactly) ----------------------------
cases = [
    ("FULL covariance, all 40 bins (primary)", full),
    ("STAT-only errors, all 40 bins", stat),
    ("FULL covariance, top-2 high-z bins dropped (calibration-limited tail)",
     Fit(list(range(NALL - 2)), "full")),
]
RES = {}
for k, (name, fit) in enumerate(cases, 1):
    prof = profile(fit)
    Lb, cb, xb = min(prof, key=lambda p: p[1])
    c0 = [p for p in prof if abs(p[0]) < 1e-9][0]
    lo2, hi2 = interval(prof, 3.84)
    RES[k] = dict(prof=prof, Lb=Lb, cb=cb, xb=xb, c0=c0, lo2=lo2, hi2=hi2)
    print("\n[%d] %s" % (k, name))
    print("    L      Lambda[m^-2]   best xi   chi2")
    for p in prof:
        if abs(p[0] * 10 - round(p[0] * 10)) < 1e-6 and p[0] > -0.55:
            print("   %+.1f   %+.2e     %.3f   %6.1f" % (p[0], p[0] * UNIT, p[2], p[1]))
    print("    best: L = %+.2f (xi = %.3f, chi2 = %.1f) -- at the scan EDGE L = 0.5,"
          % (Lb, xb, cb))
    print("          where 2 L x^2 -> 1 and the weak-field ansatz fails")
    print("    Lambda = 0: chi2 = %.1f (xi = %.3f);  Delta chi2 = %.1f"
          % (c0[1], c0[2], c0[1] - cb))
    print("    95%% (dchi2 = 3.84): %+.2f <= L <= %+.2f  ->  %+.1e .. %+.1e m^-2"
          % (lo2, hi2, lo2 * UNIT, hi2 * UNIT))

# [4] the degeneracy ----------------------------------------------------------
print("\n[4] The drift-Lambda DEGENERACY (primary case)")
P1 = RES[1]["prof"]
pa = [p for p in P1 if abs(p[0]) < 1e-9][0]
pb = [p for p in P1 if abs(p[0] - 0.3) < 1e-9][0]
slope = (pb[2] - pa[2]) / 0.3
print("    along the chi2 valley the best drift rate rises with Lambda:")
print("    xi ~ %.2f + %.2f L   (a faster G-drift dims more; Lambda > 0 brightens" % (pa[2], slope))
print("    by adding redshift at fixed distance -- the two trade off)")
ok("valley: xi rises monotonically with L for L in [0, 0.5]",
   all(P1[i][2] <= P1[i + 1][2] + 1e-12 for i in range(len(P1) - 1) if P1[i][0] >= 0))
ok("valley is FLAT: chi2 varies by < 8 from Lambda = 0 to the Planck value (L=0.32)",
   abs(pa[1] - [p for p in P1 if abs(p[0] - 0.32) < 1e-9][0][1]) < 8)
ok("full covariance, all bins: Lambda = 0 is DISFAVOURED at 95% (Delta chi2 > 3.84)",
   RES[1]["c0"][1] - RES[1]["cb"] > 3.84)
ok("... but that rests on the calibration-limited tail: drop the top-2 bins and "
   "Lambda = 0 is inside the 95% interval", RES[3]["lo2"] <= 0 <= RES[3]["hi2"])
ok("the Planck value (L = 0.32) is inside the 95% interval in both full-cov cases",
   all(RES[k]["lo2"] <= 0.32 <= RES[k]["hi2"] for k in (1, 3)))
ok("no upper bound inside the ansatz: every case's 95% interval reaches the L=0.5 "
   "edge, where the weak-field coupling fails", all(RES[k]["hi2"] >= 0.498 for k in RES))
ok("but large NEGATIVE Lambda is excluded (L < -0.5 => Delta chi2 > 3.84)",
   all(p[1] - RES[1]["cb"] > 3.84 for p in P1 if p[0] < -0.5))

# [5] what fixes Lambda: an independent drift rate ---------------------------
print("\n[5] Lambda is fixed only once the drift rate xi is fixed independently")
def best_L_at_xi(fit, xi):
    out = []
    for L in LS:
        s = shape(fit.z, xi, L)
        if s is not None:
            out.append((fit.chi2(s), L))
    return min(out)
for xi_fix, label in [(0.229, "fitted drift alone (Lambda=0 best fit)"),
                      (0.0, "no drift (e.g. if the LLR Gdot/G bound holds in DT units)"),
                      (1.0, "linearised Gdot/G = H0 (Williams' (1-H0 t), first order)")]:
    c2, Lx = best_L_at_xi(full, xi_fix)
    print("    xi = %.3f [%s]:" % (xi_fix, label))
    print("        best L = %+.3f -> Lambda = %+.2e m^-2 (%.2f x Planck), chi2 = %.1f"
          % (Lx, Lx * UNIT, Lx * UNIT / LAM_OBS, c2))
c2a, La = best_L_at_xi(full, 0.229)
c2b, Lb0 = best_L_at_xi(full, 0.0)
c2c, Lc = best_L_at_xi(full, 1.0)
ok("fitted drift (xi=0.229) => Lambda consistent with ZERO (|L| < 0.02)", abs(La) < 0.02)
ok("no drift (xi=0) => Lambda NEGATIVE (L < 0)", Lb0 < 0)
ok("Gdot/G = H0 (xi=1) => Lambda pushed to the weak-field edge (L >= 0.45, ~Planck or more)",
   Lc >= 0.45)

print("\n  NET. The SN Hubble diagram does NOT separate DT's drift from Lambda: along")
print("  the valley xi ~ 0.23 + 1.5 L the fit changes by only a few in chi2. Within it")
print("  the data lean mildly TOWARD a positive Lambda of Planck size with a faster")
print("  drift (full covariance: Delta chi2 = 6.9 for Lambda = 0, ~2.6 sigma), but")
print("  that lean rests on the two calibration-limited high-z bins (drop them:")
print("  Delta chi2 = 3.1, Lambda = 0 allowed), and the best fit sits at the edge")
print("  where the weak-field coupling fails. So 'Lambda small' is NOT a fit result;")
print("  it follows only if the drift rate is fixed at its drift-alone value")
print("  (xi ~ 0.23). Lambda's size is tied to the drift rate, so an independent")
print("  measurement of xi -- Gdot/G, the open LLR question -- would fix Lambda.")
print("  Order-of-magnitude only (leading-order static-patch coupling; areal and")
print("  lensing terms omitted).")

print("\n" + "=" * 72)
print("dt_lambda_joint_fit COMPLETE.")
print("=" * 72)
