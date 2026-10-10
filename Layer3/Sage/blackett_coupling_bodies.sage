#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, The Dynamic Theory (2011), Ch. 10:
#   Eq 10.17-10.19 (beta fixed by turning the vector-field force into Newton's
#   law; printed sqrt(4 pi eps0 G) but valued 2.4296e-11 = sqrt(eps0 G)),
#   Eq 10.53 (beta = sqrt(4 pi eps G) = 8.6127e-11), the Earth moment
#   mu = (q_eff / 2M) I omega = 8.6e22 A m^2 (Book p.295), and Ch. 6 / p.437:
#   the moment of a neutral spinning body comes from "the fifth dimensional
#   field, being the gravitational field". Builds on Layer-3 script 20
#   (20_earth_magnetic_moment.sage, frozen corpus).
# Blackett, Nature 159, 658 (1947); Phil. Trans. A 245, 309 (1952).
# Data (imported; limits taken CONSERVATIVELY, i.e. looser than published):
#   moments: Earth 8.0e22, Mercury 2.8e19 (MESSENGER), Jupiter 1.56e27 (Juno),
#   Saturn 4.6e25 (Cassini) A m^2; Sun dipole B_pole ~ 1 G; millisecond pulsar
#   B ~ 3e8 G (mu ~ 3e23 A m^2). Upper limits: Moon < 1e16, Mars < 2e18,
#   Venus < 1e19 A m^2. Earth paleointensity falls below ~20% of today's
#   during reversals.
# ---------------------------------------------------------------------------
#
# blackett_coupling_bodies.sage                       [A] (assertion-hardened)
#
# QUESTION (manuscript item 4). Williams gives a spinning neutral mass the
# magnetic moment mu = (beta/2) L with beta = sqrt(eps0 G): Blackett's 1947
# conjecture, Earth to ~6%. Is it a real magnetic moment, and does it survive
# the bodies that magnetometers have measured?
# Odds before computing: as a magnetometer field, excluded ~85%; the
# gravitational (gravitomagnetic) reading survives as recovery ~80%.
#
#   [A] the two values of beta: only sqrt(4 pi eps0 G) maps Coulomb onto
#       Newton; the Earth match uses the other one
#   [B] Williams' law is Blackett's with beta_B = 1/sqrt(4 pi); the 1952
#       laboratory test does not reach a Maxwellian version
#   [C] read as a magnetometer field: Moon, Mars, Venus, millisecond pulsars,
#       the Sun and Earth's reversals exclude it
#   [D] read as Williams' own source says (the gravitational field): converted
#       back with beta it is the gravitomagnetic dipole field, a GR quantity
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("blackett_coupling_bodies.sage -- Williams' mu = (beta/2) L against the bodies")
print("=" * 74)

eps0, mu0, G, c = 8.8541878128e-12, 1.25663706212e-6, 6.674e-11, 2.99792458e8
b1 = sqrt(eps0 * G).n()                    # Eq 10.19 value
b2 = sqrt(4 * pi * eps0 * G).n()           # Eq 10.53

# ============================================================ [A]
print("\n[A] beta: Eq 10.19 value %.4e, Eq 10.53 %.4e C/kg" % (b1, b2))
check(abs(b1 - 2.4296e-11) / 2.4296e-11 < 2e-3 and abs(b2 - 8.6127e-11) / 8.6127e-11 < 2e-3,
      "the book carries both values, a factor sqrt(4 pi) = 3.545 apart")
# effective charge q = beta M; Coulomb field q/(4 pi eps0 r^2) must give, on a
# test charge beta m, Newton's force G M m / r^2
bb, MM, mm, rr = var('bb MM mm rr', domain='positive')
F_C = (bb * MM) * (bb * mm) / (4 * pi * eps0 * rr**2)
sol = solve(F_C == G * MM * mm / rr**2, bb)
check(any(abs(s.rhs().n() - b2) / b2 < 1e-9 for s in sol),
      "Coulomb(beta M, beta m) = Newton  <=>  beta = sqrt(4 pi eps0 G) (Eq 10.53)")
print("    So the self-consistent ratio is Eq 10.53; the Earth estimate used Eq 10.19.")

# ============================================================ [B]
print("\n[B] Blackett's law, P = beta_B sqrt(G) U / (2c) (Gaussian), in SI:")
print("    mu = beta_B sqrt(4 pi eps0 G) L / 2. Williams: mu = sqrt(eps0 G) L / 2.")
betaB = b1 / b2
check(abs(betaB - 1 / sqrt(4 * pi).n()) < 1e-12,
      "Williams' law = Blackett's with beta_B = 1/sqrt(4 pi) = 0.282 (Blackett fitted beta_B ~ 0.3)")
print("    Blackett's 1952 null result tested a LOCAL law B ~ beta sqrt(G) rho v / c")
print("    for a gold cylinder carried by the Earth (prediction ~1e-13 T).")
rho_au, v_E, R_cyl, om_E = 19300.0, 340.0, 0.05, 7.292e-5
B_maxwell = mu0 * b1 * rho_au * om_E * R_cyl**2 / 3        # rotating uniformly charged sphere, centre
print("    A Maxwellian version (rho_q = beta rho, Ampere's law) gives for the")
print("    cylinder's own rotation with the Earth B ~ %.1e T." % B_maxwell)
check(B_maxwell < 1e-18,
      "Maxwellian prediction < 1e-18 T: far below the 1952 sensitivity; that test does not reach it")

# ============================================================ [C]
print("\n[C] Read as a magnetometer field: mu_DT = (beta/2) C M R^2 omega")
day, hr = 86400.0, 3600.0
#  name         M (kg)     R (m)      C      period (s)        measured  kind
bodies = [
    ("Moon",     7.342e22, 1.7374e6, 0.393, 27.32 * day,      1.0e16,  "lim"),
    ("Mars",     6.417e23, 3.3895e6, 0.3645, 1.026 * day,     2.0e18,  "lim"),
    ("Venus",    4.867e24, 6.0518e6, 0.330, 243.0 * day,      1.0e19,  "lim"),
    ("Mercury",  3.301e23, 2.4397e6, 0.346, 58.65 * day,      2.8e19,  "obs"),
    ("Earth",    5.972e24, 6.371e6,  0.3307, 0.99727 * day,   8.0e22,  "obs"),
    ("Saturn",   5.683e26, 5.8232e7, 0.210, 10.56 * hr,       4.6e25,  "obs"),
    ("Jupiter",  1.898e27, 6.9911e7, 0.254, 9.925 * hr,       1.56e27, "obs"),
]
Rsun = 6.957e8
mu_sun = 2 * pi.n() * Rsun**3 * 1e-4 / mu0                  # B_pole ~ 1 G
L_sun = 1.92e41
ratio = {}
print("    %-9s %10s %10s %10s %8s" % ("body", "L", "mu(10.19)", "measured", "ratio"))
for n, M_, R_, C_, P_, meas, kind in bodies:
    L_ = C_ * M_ * R_**2 * 2 * pi.n() / P_
    mu_ = b1 / 2 * L_
    ratio[n] = mu_ / meas
    print("    %-9s %10.2e %10.2e %9s%.1e %8.2g" % (n, L_, mu_, "<" if kind == "lim" else " ", meas, ratio[n]))
mu_s = b1 / 2 * L_sun; ratio["Sun"] = mu_s / mu_sun
print("    %-9s %10.2e %10.2e %10.1e %8.2g" % ("Sun", L_sun, mu_s, mu_sun, ratio["Sun"]))
L_msp = 1e38 * 2 * pi.n() / 3e-3; mu_msp = b1 / 2 * L_msp; ratio["MSP"] = mu_msp / 3e23
print("    %-9s %10.2e %10.2e %10.1e %8.2g" % ("ms pulsar", L_msp, mu_msp, 3e23, ratio["MSP"]))
check(0.7 < ratio["Earth"] < 1.2, "Earth (modern C = 0.3307): within ~10% with beta = sqrt(eps0 G)")
check(ratio["Earth"] * (b2 / b1) > 3, "with the self-consistent beta (Eq 10.53) Earth is off by >3x")
check(ratio["Moon"] > 100 and ratio["Mars"] > 500,
      "Moon (>100x) and Mars (>500x) exceed their limits even with conservative limits")
check(ratio["Venus"] > 10, "Venus exceeds its limit by >10x")
check(ratio["MSP"] > 1e6, "millisecond pulsars: predicted moment >1e6 x the observed")
print("    Moon and Mars have no active dynamo, so nothing can cancel an always-on")
print("    contribution. The Sun's dipole reverses every ~11 yr, and Earth's")
print("    paleointensity drops below ~20% in reversals; a fixed contribution of")
print("    ~90% of today's moment, parallel to the spin, forbids both.")
check(ratio["Earth"] > 0.2, "Earth's fixed part (~%d%% of today) exceeds the <20%% reversal minima" % round(100 * ratio["Earth"]))
print("    => As a field that magnetometers measure, the law is EXCLUDED, with")
print("       margins of 1e2-1e7 and no dynamo assumption (Moon, Mars, pulsars).")

# ============================================================ [D]
print("\n[D] Read as the source says: the gravitational field of the 5D gauge field")
J, r = var('J r', domain='positive')
# Williams' field of a dipole mu = beta J / 2, acting on a test 'charge' beta m
# moving at v: F = beta m v B, with B = mu0 mu / (4 pi r^3) (equatorial scale)
for bname, bval in (("10.19", b1), ("10.53", b2)):
    acc = bval * mu0 * (bval * J / 2) / (4 * pi * r**3)       # per unit v
    ratio_gm = (acc / (G * J / (c**2 * r**3))).n()
    print("    beta (Eq %s): acceleration per unit v = %.4g x G J /(c^2 r^3)" % (bname, ratio_gm))
acc2 = (b2 * mu0 * (b2 * J / 2) / (4 * pi * r**3))
check(abs((acc2 / (G * J / (c**2 * r**3))).n() - 0.5) < 1e-6,
      "with the self-consistent beta the field is G J / (2 c^2 r^3) per unit v: a gravitomagnetic dipole")
print("    No electromagnetic constant survives: it is gravitomagnetism (frame")
print("    dragging), the field GR gives for a spinning mass (Lense-Thirring),")
print("    whose normalisation the recovered metric (section 3) fixes and Gravity")
print("    Probe B and LAGEOS/LARES measure. It is not seen by a magnetometer.")
print("""
NET.
  * Williams' moment mu = (beta/2) L is Blackett's 1947 law with beta_B =
    1/sqrt(4 pi). Its ~6% Earth match uses beta = sqrt(eps0 G) (Eq 10.19);
    the ratio that maps Coulomb onto Newton is sqrt(4 pi eps0 G) (Eq 10.53),
    which puts Earth off by 3.5x. The book holds both: an internal slip.
  * READ AS A MAGNETOMETER FIELD (Williams' comparison): EXCLUDED. The Moon
    (>280x) and Mars (>1000x) have no dynamo to cancel it, millisecond
    pulsars miss by ~1e7, and solar and geomagnetic reversals forbid a fixed
    spin-aligned part. A WILLIAMS DISCREPANCY at the criminal standard.
  * READ AS THE SOURCE SAYS (the fifth-dimensional, gravitational field):
    with the self-consistent beta the field is G J/(2 c^2 r^3), the
    gravitomagnetic dipole of a spinning mass. That is GR's frame dragging,
    recovered, not a new coupling: NO PREDICTION.
  * The 1952 laboratory null tested Blackett's local law; a Maxwellian
    version predicts < 1e-18 T there, so that experiment decides nothing.
""")
print("ALL ASSERTIONS PASSED")
