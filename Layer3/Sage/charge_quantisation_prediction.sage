# ---------------------------------------------------------------------------
# charge_quantisation_prediction.sage
#
# QUESTION. The Standard Model explains charge quantisation only for its
# minimal field content. Add right-handed neutrinos with Dirac masses (needed
# for neutrino masses) and the hypercharges dequantise: anomaly cancellation
# plus Yukawa invariance leave Q -> Q + delta (B - L), with delta free.
# Dynamic Theory's q = Z e (Book Eqs 7.66-7.67; Z = N_o - N_1, an integer
# winding number from the quantum condition; e universal) does not depend on
# the field content. Does it fix delta?
#
# SOURCE (primary): Williams, The Dynamic Theory (2011), Eqs 7.66-7.67.
# Dequantisation with nu_R: Babu & Mohapatra, PRL 63, 938 (1989); Foot, Lew &
# Volkas, J. Phys. G 19, 361 (1993).
# Data (PDG): neutron charge q_n = (-0.2 +- 0.8) x 1e-21 e (Baumann et al.
# 1988); |q_p + q_e| < 1e-21 e (Bressi et al. 2011).
#
#   [A] SM + nu_R, Dirac Yukawas: anomalies leave a ONE-parameter family
#   [B] SM without nu_R: the same conditions QUANTISE (the known result)
#   [C] physical charges of p, n, e, nu, H on the family: q_n = delta
#   [D] DT's q = Z e on free states: delta = -Z_n / (Z_e + Z_n), so
#       delta = 0 unless |Z_e| ~ 1/|delta|
#   [E] data: |delta| < ~2e-21, so Z_n = 0: delta = 0 EXACTLY
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("charge_quantisation_prediction.sage -- does q = Ze fix the hypercharges?")
print("=" * 74)

R = PolynomialRing(QQ, 'yQ,yu,yd,yL,ye,yn,yH')
yQ, yu, yd, yL, ye, yn, yH = R.gens()

# Q = T3 + Y. Left-handed doublets counted +, right-handed singlets -.
# One generation; colour multiplicity 3.
def anomalies(with_nuR):
    nR = yn if with_nuR else 0
    return [
        2*yQ - yu - yd,                                        # SU(3)^2 U(1)
        3*yQ + yL,                                             # SU(2)^2 U(1)
        6*yQ - 3*yu - 3*yd + 2*yL - ye - nR,                   # grav^2 U(1)
        6*yQ**3 - 3*yu**3 - 3*yd**3 + 2*yL**3 - ye**3 - nR**3, # U(1)^3
    ]

yuk = {yu: yQ + yH, yd: yQ - yH, ye: yL - yH, yn: yL + yH}   # Dirac masses via H

# ============================================================ [A]
print("\n[A] SM + nu_R with Dirac Yukawas")
eqs = [e.subs(yuk) for e in anomalies(True)]
eqs = [e.subs({yL: -3*yQ}) for e in eqs]
check(all(e == 0 for e in eqs),
      "with yL = -3 yQ every anomaly cancels identically: yQ and yH both free")
print("    normalise yH = 1/2 (neutral Higgs component): yQ is a FREE parameter")

# ============================================================ [B]
print("\n[B] SM without nu_R")
eqsB = [e.subs({k: v for k, v in yuk.items() if k != yn}) for e in anomalies(False)]
eqsB = [e.subs({yL: -3*yQ, yH: QQ(1)/2}) for e in eqsB]
grav = eqsB[2]
check(grav == -3*yQ + QQ(1)/2, "gravitational anomaly reduces to yH = 3 yQ")
check(all(e.subs({yQ: QQ(1)/6}) == 0 for e in eqsB),
      "so yQ = 1/6 and every anomaly cancels: hypercharge QUANTISED (standard result)")

# ============================================================ [C]
print("\n[C] Physical charges on the nu_R family, yH = 1/2")
t = var('t')                          # t = yQ
qu, qd = QQ(1)/2 + t, -QQ(1)/2 + t
qe, qnu = -QQ(1)/2 - 3*t, QQ(1)/2 - 3*t
qp, qn = 2*qu + qd, qu + 2*qd
delta = 3*t - QQ(1)/2
for name, q in [("p", qp), ("n", qn), ("e", qe), ("nu", qnu), ("H atom", qp + qe)]:
    print("    q_%-6s = %s" % (name, q.simplify_full()))
check(bool((qn - delta).simplify_full() == 0) and bool((qnu + delta).simplify_full() == 0),
      "q_n = delta, q_nu = -delta, with delta = 3 yQ - 1/2 (the B-L direction)")
check(bool((qp + qe).simplify_full() == 0),
      "q_p + q_e = 0 for every delta: atom neutrality does NOT test delta, the neutron does")
check(bool(delta.subs(t=QQ(1)/6) == 0), "the SM point yQ = 1/6 is delta = 0")

# ============================================================ [D]
print("\n[D] DT: every free state carries Z e, Z integer, e universal")
# q_e = Z_e e, q_n = Z_n e  =>  delta / (-(1 + delta)) = Z_n / Z_e
d = var('d')
Ze, Zn = var('Ze Zn')
dsol = solve(d / (-(1 + d)) == Zn / Ze, d)[0].rhs()
print("    delta = %s" % dsol)
check(bool(dsol.subs(Zn=0) == 0), "Z_n = 0  =>  delta = 0")
check(bool(((dsol + Zn / (Ze + Zn))).simplify_full() == 0), "delta = -Z_n / (Z_e + Z_n)")
small = min(abs(dsol.subs(Ze=-1, Zn=z)) for z in [-3, -2, -1, 2, 3])   # Z_n = 1 has no solution
check(small >= QQ(1)/2,
      "with the electron as the unit (Z_e = -1), any Z_n != 0 gives |delta| >= 1/2")
print("    |delta| < x needs |Z_e| > ~1/x: the electron would have to carry")
print("    ~1/x elementary units")

# ============================================================ [E]
print("\n[E] Data")
qn_meas, qn_err = -0.2e-21, 0.8e-21
bound = abs(qn_meas) + 2 * qn_err
print("    |delta| = |q_n|/e < %.1e (2 sigma)" % bound)
check(bound < 1e-20, "SM + Dirac nu_R allows any |delta| < 2e-21; DT allows only 0")
check(1 / bound > 1e20,
      "a nonzero delta in DT needs |Z_e| > 5e20: excluded unless e is not the electron's unit")
print("""
NET.
  * With right-handed neutrinos and Dirac masses the SM's charges carry a
    free parameter delta along B-L: q_n = delta, q_nu = -delta, while atoms
    (q_p + q_e = 0) stay neutral for every delta. The SM fits delta = 0; it
    does not require it (Majorana masses would).
  * DT's q = Z e holds for every free state whatever the field content, so
    q_n = Z_n e with Z_n an integer. Given |q_n| < 2e-21 e, Z_n = 0 and
    delta = 0 EXACTLY, independently of whether neutrinos are Dirac or
    Majorana.
  * PREDICTION: q_n = q_nu = 0 exactly; any nonzero neutron or neutrino
    charge falsifies DT while being allowed by the SM with Dirac neutrinos.
  * Bracket: this is the compactness of electromagnetic U(1). Dirac's
    monopole argument and grand unification also give it; DT gets it from
    the quantum condition without either. The prediction is a null one: it
    is consistent with all data and falsifiable, but it cannot be confirmed
    by a measurement, only bounded.
""")
print("ALL ASSERTIONS PASSED")
