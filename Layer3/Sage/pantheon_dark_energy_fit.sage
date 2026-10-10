#!/usr/bin/env sage
# --------------------------------------------------------------------------
# SOURCE: model = P. E. Williams, "The Dynamic Theory" Ch. 3 (DT redshift
#   relation Eq 3.6 + the evolving Chandrasekhar candle, Eqs 3.52-3.55).
#   DATA = Pantheon binned Type-Ia supernovae (Scolnic et al. 2018,
#   ApJ 859:101; Binned_data/lcparam_DS17f.txt), 40 bins (zcmb, mb, dmb),
#   fetched from the public repository github.com/dscolnic/Pantheon.
# --------------------------------------------------------------------------

# ============================================================================
# 27_pantheon_dark_energy_fit.sage              [A]  (assertion-hardened)
#
# THE LIVE, FALSIFIABLE TEST: does DT (log-redshift distance + evolving candle)
# fit the REAL SN-Ia Hubble diagram as well as LambdaCDM, and far better than a
# no-dark-energy decelerating (Einstein-de Sitter) universe?
#
# Method: SN data fix the SHAPE m_b(z), not the absolute scale, so for each model
# the constant offset (= absolute magnitude M + 5 log10(c/H0) + const) is
# marginalised analytically (inverse-variance weighted). chi^2 = sum w (m_b -
# shape - offset)^2 then measures shape agreement only (H0-independent).
#
#   shape_LCDM(z)   = 5 log10[(1+z) INT_0^z dz'/E],  E=sqrt(0.3(1+z)^3+0.7)
#   shape_EdS(z)    = 5 log10[(1+z) INT_0^z dz'/(1+z')^{3/2}]      (Omega_M=1)
#   shape_DT(z,xi)  = 5 log10[(1+z) ln(1+z)]  +  (15/4) xi ln(1+z)/ln10
#                     \___ DT redshift distance (Eq 3.6) ___/  \__ evolving candle __/
#
# RUN: /home/dsteele/miniconda3/envs/sage/bin/sage 27_pantheon_dark_energy_fit.sage
# ============================================================================

from sage.all import *


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


LN10 = float(log(10))

# --- Pantheon binned data: (zcmb, mb, dmb), Scolnic+ 2018 -------------------
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
zs = [d[0] for d in DATA]
mb = [d[1] for d in DATA]
dmb = [d[2] for d in DATA]
w = [1.0 / s**2 for s in dmb]
N = len(DATA)


def chi2(shape):
    # marginalise the constant offset analytically (weighted mean of residuals)
    resid = [mb[i] - shape[i] for i in range(N)]
    off = sum(w[i] * resid[i] for i in range(N)) / sum(w)
    return float(sum(w[i] * (resid[i] - off)**2 for i in range(N)))


def shape_FRW(Om):
    out = []
    for z in zs:
        I = numerical_integral(lambda zp: 1.0 / sqrt(Om * (1 + zp)**3 + (1 - Om)), 0, z)[0]
        out.append(5.0 * log((1 + z) * I) / LN10)
    return out


def shape_DT(xi):
    return [5.0 * log((1 + z) * log(1 + z)) / LN10 + (15.0 / 4.0) * xi * log(1 + z) / LN10
            for z in zs]


print("=" * 70)
print("PANTHEON SN-Ia HUBBLE DIAGRAM -- DT vs LambdaCDM vs Einstein-de Sitter")
print("=" * 70)
print("    %d binned SNe, 0.014 < z < 1.61 (Scolnic+ 2018)" % N)

# --- LambdaCDM (Om=0.3) and EdS (Om=1) --------------------------------------
chi2_LCDM = chi2(shape_FRW(0.3))
chi2_EdS = chi2(shape_FRW(1.0))
print("\n[1] LambdaCDM (Om=0.3) vs Einstein-de Sitter (Om=1, NO dark energy)")
print("    chi2_LCDM = %.1f  (dof=%d, chi2/dof=%.2f)" % (chi2_LCDM, N - 1, chi2_LCDM / (N - 1)))
print("    chi2_EdS  = %.1f  (dof=%d, chi2/dof=%.2f)" % (chi2_EdS, N - 1, chi2_EdS / (N - 1)))
print("    Delta chi2 (EdS - LCDM) = %.1f  <- the dark-energy detection significance"
      % (chi2_EdS - chi2_LCDM))
ok("LambdaCDM is a good fit (chi2/dof < 1.5).", chi2_LCDM / (N - 1) < 1.5)
ok("Einstein-de Sitter (no dark energy) is strongly disfavoured (Delta chi2 > 25).",
   (chi2_EdS - chi2_LCDM) > 25)

# --- DT: scan xi (evolving-candle rate) -------------------------------------
print("\n[2] DT (log-z distance + evolving candle): scan the candle rate xi")
best = None
xi = 0.0
while xi <= 1.0001:
    c2 = chi2(shape_DT(xi))
    if best is None or c2 < best[1]:
        best = (xi, c2)
    xi += 0.005
xi_best, chi2_DT = best
chi2_DT_nocandle = chi2(shape_DT(0.0))
print("    DT, no candle (xi=0, pure log-z): chi2 = %.1f (chi2/dof=%.2f)"
      % (chi2_DT_nocandle, chi2_DT_nocandle / (N - 1)))
print("    DT best fit: xi = %.3f -> chi2 = %.1f (dof=%d, chi2/dof=%.2f)"
      % (xi_best, chi2_DT, N - 2, chi2_DT / (N - 2)))
ok("DT best-fit candle rate xi is order 0.1-0.5 (consistent with scripts 24/25).",
   0.1 < xi_best < 0.5)
ok("the evolving candle is NEEDED: it improves DT's fit over pure log-z by >25 in chi2.",
   (chi2_DT_nocandle - chi2_DT) > 25)

# --- Verdict: is DT competitive with LambdaCDM, and far better than EdS? ----
print("\n[3] Verdict on the Hubble diagram")
print("    Delta chi2 (DT - LCDM)  = %+.1f  (mild residual preference for LambdaCDM)"
      % (chi2_DT - chi2_LCDM))
print("    Delta chi2 (EdS - DT)   = %+.1f" % (chi2_EdS - chi2_DT))
absorbed = (chi2_EdS - chi2_DT) / (chi2_EdS - chi2_LCDM)
print("    DT absorbs %.1f%% of the dark-energy chi2 excess (EdS->LCDM)." % (100 * absorbed))
ok("DT fits FAR better than no-dark-energy EdS (by >1000 in chi2).",
   (chi2_EdS - chi2_DT) > 1000)
ok("DT absorbs >90% of the dark-energy chi2 excess (most of the acceleration signal).",
   absorbed > 0.90)
ok("but a RESIDUAL mild preference for LambdaCDM remains (DT not a perfect fit).",
   (chi2_DT - chi2_LCDM) > 10)

print("\n  NET: on the REAL Pantheon Hubble diagram, DT (its log-redshift distance")
print("  plus an evolving Chandrasekhar candle, xi~%.2f) ABSORBS ~98%% of the" % xi_best)
print("  dark-energy signal: it crushes the no-dark-energy Einstein-de Sitter fit")
print("  (chi2 1469 -> 80) and lands close to LambdaCDM (chi2 51). But a residual")
print("  Delta chi2 ~ %.0f favours LambdaCDM, so dark energy does NOT vanish" % (chi2_DT - chi2_LCDM))
print("  completely -- DT is a near-miss, not an exact substitute, on current SN data.")
print("  And even this near-fit needs the |Gdot/G| question settled (script 26):")
print("  the candle rate xi~0.24 is a local-G drift (excluded) unless it is a gauge")
print("  effect (allowed, but then degenerate with the very data it fits here).")

print("\n" + "=" * 70)
print("27 COMPLETE.")
print("=" * 70)
