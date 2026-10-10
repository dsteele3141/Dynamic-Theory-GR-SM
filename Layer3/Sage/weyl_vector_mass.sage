#!/usr/bin/env sage
# --------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, LA-8370-MS (LASL, 1980), Eq (39) (4D) and
#   Eq (67) (5D):  delta INT [R + (1/4) alpha F^ F^ - lambda(1/2 - k phi^.phi^)]
#   sqrt(-g) = 0, k = 6 / 12, "Weyl's procedure in his unified field theory".
#   The EM reading of the gauge vector is the corpus' London coupling
#   phi = i (e/hbar) A  (gateway_london_quantisation.sage: alpha1 = i e^2/(hbar c);
#   T1_8_weyl_integrability.sage part B).  London (1927), Weyl (1929).
# Full derivation, context and provenance:
#   Supplementary_Material.pdf -- §S1 (T1.8, T1.9), Appendix A.
# Companions: weyl_lambda_reading.sage (lambda is a constant, Lambda != 0),
#   einstein_from_dt_action.sage (metric variation of Eq 67).
# --------------------------------------------------------------------------

# ============================================================================
# weyl_vector_mass.sage                         [A]  (assertion-hardened)
#
# THE GAUGE-VECTOR MASS IN EQS 39 / 67, AND WHAT IT DOES TO THE GR SECTOR.
#
# weyl_lambda_reading left two items open: the mass the constant lambda gives
# the gauge vector (it needs the F^2 coefficient alpha), and the effect of the
# imaginary London gauge.  Both are fixed here, convention-free:
#
#   [1] The phi.phi term is GEOMETRIC: for constant phi on flat space
#       R~ = -(n-1)(n-2) phi.phi (k = 6 in 4D, 12 in 5D), so its coefficient
#       relative to R is fixed -- it does not depend on lambda or its sign.
#   [2] Dispersion of L = R - k phi.phi + (1/4) a F^2 (a = alpha in the
#       physical, un-hatted field): Omega^2 = K^2 + 2k/a.  m^2 = 2k/a.
#   [3] Which sign of a is physical.  The Ricci convention here gives the
#       weak-field G_00 = 2 lap(Phi), so attractive gravity is G = +kappa T;
#       metric variation gives G^ab = a S^ab (S = F-stress), and positive
#       EM energy then needs a > 0 for a REAL gauge vector, a < 0 for an
#       IMAGINARY one (phi = i q A flips S^00).
#   [4] The London coupling q = e/hbar fixes |a| = 2 kappa/(mu0 q^2), hence
#             m^2 = -k mu0 e^2/(kappa hbar^2) = -k alpha_fs/(2 l_P^2):
#       a TACHYONIC mass of order the Planck mass (0.15 / 0.21 M_P for
#       k = 6 / 12), independent of lambda.  So the imaginary EM phase cannot
#       live inside the Eq 39/67 connection: there R~ is complex and the
#       "photon" is a Planck-scale tachyon.  It must be the Weyl-1929 form --
#       a phase on matter, Maxwell's F^2 separate, no phi_I in R~.
#   [5] The REAL, Weyl-integrable (gravitational) sector -- phi = grad of
#       ln f, F = 0 (T1.8 part A; Report Eq 65) -- is forced to phi = 0 on
#       shell by its own field equation (2k phi = 0), so Weyl geometry
#       collapses to Riemannian and the metric equations are EXACTLY
#       G_ab + Lambda g_ab = kappa T_ab.  Mutation: a real NON-integrable
#       vector (F != 0) is a healthy Proca field with m^2 = +k mu0 q^2/kappa,
#       q unfixed -- Weyl's 1918 situation.
#
# VERDICT (honest bracket).  On Weyl's procedure the gauge vector in Eqs
# 39/67 carries a Planck-scale mass whatever lambda is.  For DT's gravitational
# sector this is good news: the integrable real Weyl vector vanishes on shell
# and the u4 = 0 surface carries exactly Einstein + Lambda -- the Weyl
# "second-clock" structure drops out dynamically, not by assumption.  For the
# EM sector it is a constraint on how DT must be read: the London phase
# phi = i(e/hbar)A cannot be the same phi as in Eq 39/67's R~ (that gives a
# complex action and a Planck-mass tachyonic photon); it is consistent only
# as the Weyl-1929 phase on matter fields, with Maxwell's F^2 added
# separately.  The corpus' T1.8 split (real integrable / imaginary EM) is
# compatible with this IF the imaginary part is kept out of Eq 39/67.  Not
# settled: matter that carries Weyl weight would source the real phi; the
# 5D (k = 12) number uses the 4D kappa of the u4 = 0 surface.
#
# RUN: /home/dsteele/miniconda3/envs/sage/bin/sage weyl_vector_mass.sage
# ============================================================================

from sage.all import *


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


def z(e):
    return bool(SR(e).simplify_full() == 0)


t, x, y, w, v = var('t x y w v')
XALL = [t, x, y, w, v]


def inverse(g):
    return [list(r) for r in matrix(SR, g).inverse()]


def levi_civita(g, gi, X):
    n = len(X)
    dg = [[[diff(g[a][b], X[c]) for c in range(n)] for b in range(n)]
          for a in range(n)]
    return [[[sum(gi[k][s] * (dg[s][j][i] + dg[s][i][j] - dg[i][j][s])
                  for s in range(n)) / 2
              for j in range(n)] for i in range(n)] for k in range(n)]


def ricci(G, X):
    n = len(X)
    return [[sum(diff(G[r][vv][s], X[r]) - diff(G[r][r][s], X[vv])
                 + sum(G[r][r][l] * G[l][vv][s] - G[r][vv][l] * G[l][r][s]
                       for l in range(n))
                 for r in range(n))
             for vv in range(n)] for s in range(n)]


def weyl_scalar(g, ph, X):
    """Scalar curvature of Gamma + (delta phi + delta phi - g phi)."""
    n = len(X)
    gi = inverse(g)
    phu = [sum(gi[k][m] * ph[m] for m in range(n)) for k in range(n)]
    Gam = levi_civita(g, gi, X)
    d = lambda a, b: 1 if a == b else 0
    Gw = [[[Gam[k][i][j] + d(k, i) * ph[j] + d(k, j) * ph[i] - g[i][j] * phu[k]
            for j in range(n)] for i in range(n)] for k in range(n)]
    Ric = ricci(Gw, X)
    return sum(gi[a][b] * Ric[a][b] for a in range(n) for b in range(n))


def eta(n):
    return [[(1 if i == 0 else -1) if i == j else 0 for j in range(n)]
            for i in range(n)]


print("=" * 72)
print("weyl_vector_mass: the gauge-vector mass in Eqs 39 / 67")
print("=" * 72)

# ==========================================================================
# [1] The phi.phi coefficient is geometric
# ==========================================================================
print("\n[1] R~ for constant phi on flat space")
c0, c1, c2, c3, c4 = var('c0 c1 c2 c3 c4', domain='real')
C = [c0, c1, c2, c3, c4]
for n, k in [(4, 6), (5, 12)]:
    X = XALL[:n]
    g = eta(n)
    ph = C[:n]
    pp = sum(g[i][i] * ph[i]**2 for i in range(n))      # eta^ii = eta_ii
    Rw = weyl_scalar(g, ph, X)
    ok("%dD: R~ = -%d phi.phi  (k = (n-1)(n-2), sign fixed by geometry)"
       % (n, k), z(Rw + k * pp))
    phI = [I * ci for ci in C[:n]]
    RwI = weyl_scalar(g, phI, X)
    ok("%dD: imaginary phi = i c flips it: R~ = +%d c.c" % (n, k),
       z(RwI - k * pp))
# non-constant imaginary phi makes R~ complex (Im R~ = -2(n-1) div phi_I)
a1f = function('A1')(t, x)
phc = [0, I * a1f, 0, 0]
Rc = weyl_scalar(eta(4), phc, XALL[:4])
ok("4D: imaginary non-constant phi = (0, i A1, 0, 0) gives the complex "
   "R~ = -6 A1^2 + 6 i d_x A1  (Im R~ = -6 div phi_I)",
   z(Rc - (-6 * a1f**2 + 6 * I * diff(a1f, x))))

# ==========================================================================
# [2] Dispersion relation of L = R - k phi.phi + (1/4) a F^2
# ==========================================================================
print("\n[2] Dispersion of the gauge vector (flat space, transverse mode)")
a, kk, P, Pt, Px, Om, K = var('a kk P Pt Px Om K')
et = eta(4)
phi = [0, 0, P, 0]
dphi = {(0, 2): Pt, (1, 2): Px}                       # d_t phi_2, d_x phi_2
Fl = [[dphi.get((i, j), 0) - dphi.get((j, i), 0) for j in range(4)]
      for i in range(4)]
F2 = sum(et[i][i] * et[j][j] * Fl[i][j]**2 for i in range(4) for j in range(4))
Lphi = -kk * sum(et[i][i] * phi[i]**2 for i in range(4)) + a * F2 / 4
# EL for a constant-coefficient quadratic L on P = exp(i(K x - Om t)):
#   dL/dP - D_t dL/dPt - D_x dL/dPx  ->  [dL/dP / P + Om^2 L_tt + K^2 L_xx] P
EL = (diff(Lphi, P) / P + Om**2 * diff(Lphi, Pt, 2) + K**2 * diff(Lphi, Px, 2))
Om2 = solve(EL == 0, Om**2)[0].rhs()
ok("Omega^2 = K^2 + 2k/a, i.e. m^2 = 2k/a", z(Om2 - (K**2 + 2 * kk / a)))

# ==========================================================================
# [3] Which sign of a is physical
# ==========================================================================
print("\n[3] Physical sign of a (convention fixed by the Newtonian limit)")
eps = var('eps')
Phi = function('Phi')(x, y, w)
X4 = [t, x, y, w]
gN = [[(1 + 2 * eps * Phi) if i == j == 0 else
       (-(1 - 2 * eps * Phi) if i == j else 0) for j in range(4)]
      for i in range(4)]
giN = inverse(gN)
RicN = ricci(levi_civita(gN, giN, X4), X4)
RN = sum(giN[i][j] * RicN[i][j] for i in range(4) for j in range(4))
G00 = diff((RicN[0][0] - gN[0][0] * RN / 2), eps).subs({eps: 0})
lapPhi = diff(Phi, x, 2) + diff(Phi, y, 2) + diff(Phi, w, 2)
ok("weak field: G_00 = 2 lap(Phi) in this Ricci convention, so Newton "
   "(lap Phi = 4 pi G rho) is G = +kappa T", z(G00 - 2 * lapPhi))
# F-stress from metric variation of sqrt|g| F^2/4 (einstein_from_dt_action
# [3]): S^ab = -1/2 F^a_j F^bj + 1/8 g^ab F^2.  Electric field E along x.
E = var('E', domain='real')


def S00(F01):
    Fd = [[0] * 4 for _ in range(4)]
    Fd[0][1], Fd[1][0] = F01, -F01
    Fu = [[et[i][i] * et[j][j] * Fd[i][j] for j in range(4)] for i in range(4)]
    Fmix = [[et[i][i] * Fd[i][j] for j in range(4)] for i in range(4)]  # F^i_j
    FF = sum(Fd[i][j] * Fu[i][j] for i in range(4) for j in range(4))
    return -sum(Fmix[0][j] * Fu[0][j] for j in range(4)) / 2 + FF / 8


q = var('q', domain='positive')
ok("real field: S^00 = +E^2/4 > 0, so G^00 = a S^00 = kappa T^00 needs a > 0",
   z(S00(E) - E**2 / 4))
ok("imaginary field F = i q F_A: S^00 = -q^2 E^2/4 < 0, so a < 0",
   z(S00(I * q * E) + q**2 * E**2 / 4))
ok("hence m^2 = 2k/a is TACHYONIC for the imaginary (London) gauge vector "
   "and healthy for a real one (k = 6 > 0 fixed by [1]): Omega^2(K=0) < 0 "
   "for a < 0, > 0 for a > 0",
   bool(Om2.subs({K: 0, kk: 6, a: -1}) < 0) and bool(Om2.subs({K: 0, kk: 6, a: 1}) > 0))

# ==========================================================================
# [4] London coupling q = e/hbar: the magnitude
# ==========================================================================
print("\n[4] London coupling q = e/hbar fixes a, hence the mass")
mu0, kap, hb, e_, G_, c_ = var('mu0 kap hb e_ G_ c_', domain='positive')
# Maxwell T^ab = (1/mu0)(-F F + 1/4 g F^2) = (2/mu0) S_A^ab;
# G = a S_phi = a (-q^2) S_A  =  kappa (2/mu0) S_A
a_L = -2 * kap / (mu0 * q**2)
ok("a = -2 kappa/(mu0 q^2) reproduces kappa T_Maxwell",
   z(a_L * (-q**2) - kap * 2 / mu0))
m2 = (2 * kk / a_L).subs({q: e_ / hb})
afs = mu0 * e_**2 * c_ / (4 * pi * hb)
lP2 = hb * G_ / c_**3
ok("m^2 = -k mu0 e^2/(kappa hbar^2) = -k alpha_fs/(2 l_P^2)  "
   "(kappa = 8 pi G/c^4), independent of lambda",
   z(m2.subs({kap: 8 * pi * G_ / c_**4}) + kk * afs / (2 * lP2))
   and z(m2 + kk * mu0 * e_**2 / (kap * hb**2)))
alpha_fs = RR(7.2973525643e-3)
lP = RR(1.616255e-35)
EP_GeV = RR(1.220890e19)
for kv, lab in [(6, 'Eq 39 (4D)'), (12, 'Eq 67 on the u4 = 0 surface')]:
    ratio = sqrt(kv * alpha_fs / 2)
    print("    %-28s k = %2d: |m| = %.3f / l_P,  M c^2 = %.3f E_P = %.2e GeV,"
          " range %.2e m" % (lab, kv, ratio, ratio, ratio * EP_GeV, lP / ratio))
ok("both are Planck-scale (0.1 < |m| l_P < 0.3)",
   all(0.1 < sqrt(kv * alpha_fs / 2) < 0.3 for kv in (6, 12)))
print("    So the London phase cannot sit inside the Eq 39/67 connection: the")
print("    'photon' would be a Planck-mass tachyon and R~ complex ([1]).  It is")
print("    consistent only as the Weyl-1929 phase on matter fields.")

# ==========================================================================
# [5] The real, integrable (gravitational) sector: phi = 0 on shell
# ==========================================================================
print("\n[5] Real Weyl-integrable sector: the vector vanishes on shell")
# Field equation from varying phi_b in R~-based L (the 2(n-1) div phi term is
# a total derivative):  -2k phi^b - a d_a F^ab = 0.
f = function('f')(t, x, y, w)
phiR = [diff(log(f) / 2, xi) for xi in X4]           # Report Eq 65 form
FR = [[diff(phiR[j], X4[i]) - diff(phiR[i], X4[j]) for j in range(4)]
      for i in range(4)]
ok("phi = d ln f^(1/2) (Report Eq 65) has F = 0", all(z(Fij) for r_ in FR for Fij in r_))
ok("with F = 0 the phi equation reads 2k phi = 0, so phi = 0 for k != 0: "
   "Weyl geometry collapses to Riemannian", z(diff(Lphi, P).subs({Pt: 0, Px: 0}) - 2 * kk * P))
print("    With phi = 0 the metric equations of Eq 39/67 are exactly")
print("    G_ab + Lambda g_ab = kappa T_ab (weyl_lambda_reading: Lambda != 0,")
print("    free; einstein_from_dt_action: the operator; gr_on_4d_surface: the")
print("    u4 = 0 restriction).")
mR = 2 * kk / (2 * kap / (mu0 * q**2))
ok("mutation: a real NON-integrable vector is a healthy Proca field, "
   "m^2 = +k mu0 q^2/kappa > 0 (q unfixed: Weyl 1918)",
   z(mR - kk * mu0 * q**2 / kap))

print("""
  NET: the gauge vector of Eqs 39/67 is Planck-heavy whatever lambda is.
  Gravity: the real integrable vector vanishes on shell, so the u4 = 0
  surface carries exactly G + Lambda g = kappa T -- Weyl's structure drops
  out dynamically.  EM: the London phase i(e/hbar)A cannot be the same phi
  (complex R~, Planck-mass tachyonic photon); it must be Weyl-1929's phase
  on matter, with Maxwell's F^2 separate.""")
print("=" * 72)
print("weyl_vector_mass COMPLETE.")
print("=" * 72)
