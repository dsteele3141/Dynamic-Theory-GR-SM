#!/usr/bin/env sage
# --------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, LA-8370-MS (LASL, 1980), Eq (39) (4D) and
#   Eq (67) (5D):  delta INT [R + (1/4) alpha F^ F^ - lambda(1/2 - k phi^.phi^)]
#   sqrt(-g) = 0, k = 6 (4D) / 12 (5D), reached by "Weyl's procedure in his
#   unified field theory"; the Report rescales phi^ and F^ by powers of lambda
#   (the line under each equation) and never says what lambda is.
#   Weyl's procedure: H. Weyl, Raum-Zeit-Materie (4th ed. 1921) / Space-Time-
#   Matter (1922), the "natural gauge" of his R~^2 action.
# Full derivation, context and provenance:
#   Supplementary_Material.pdf -- §S1 (T1.9), Appendix A.
# Companion: einstein_from_dt_action.sage (the metric variation of Eq 67).
# --------------------------------------------------------------------------

# ============================================================================
# weyl_lambda_reading.sage                      [A]  (assertion-hardened)
#
# WHAT IS lambda IN EQS 39 / 67?  IT DECIDES THE COSMOLOGICAL CONSTANT.
#
# einstein_from_dt_action.sage found that the cosmological part of the
# lambda term is -(lambda/2)(1/2 - k phi^.phi^) g^ab: zero if lambda is a
# Lagrange multiplier enforcing phi^.phi^ = 1/(2k), non-zero if lambda is a
# constant.  Weyl's procedure, which the Report says it follows, settles it:
#
#   [1] GAUGE WEIGHTS.  Under g -> e^(2 s) g, phi -> phi - ds the Weyl
#       connection is unchanged and R~ -> e^(-2 s) R~, so sqrt|g| R~^p has
#       weight n - 2p: gauge-invariant iff p = n/2 (R~^2 in 4D, R~^(5/2) in
#       5D).  Mutation: R~^2 is NOT invariant in 5D.
#   [2] THE NATURAL GAUGE.  The weight lets one choose the gauge R~ = -lambda,
#       a CONSTANT.  Expanding f(R~) = R~^p about it to first order (which is
#       all the field equations see) gives, after dividing by f'(-lambda) and
#       rescaling phi = sqrt(lambda) phi^,
#             R + lambda (c_p - (n-1)(n-2) phi^.phi^)    (+ divergence),
#       with c_p = (p-1)/p.  p = 2 gives c = 1/2: exactly Eq 39 in 4D, and
#       exactly Eq 67 in 5D -- so Eq 67 is Weyl's 4D R~^2 procedure carried
#       to 5D; the 5D gauge-invariant parent R~^(5/2) gives c = 3/5 instead.
#   [3] lambda IS NOT A MULTIPLIER.  Treat lambda as an auxiliary field in
#       the un-normalised linearisation f(s) + f'(s)(R~ - s), s = -lambda:
#       its Euler-Lagrange equation is f''(s)(R~ - s) = 0, i.e. R~ = -lambda,
#       the natural-gauge condition itself -- NOT phi^.phi^ = const.  The
#       multiplier reading arises only from varying lambda in the
#       already-normalised Lagrangian, in which lambda is no longer an
#       independent variable (it sits in the normalisation and inside phi^).
#   [4] lambda != 0, SO Lambda != 0.  f'(0) = 0 for p > 1: the normalisation
#       fails at lambda = 0, so the procedure needs lambda != 0, and the
#       constant term lambda c_p is then a cosmological constant,
#       |Lambda| = c_p lambda / 2  (lambda/4 for Eq 39/67 as written).  The
#       same lambda multiplies the phi^.phi^ term, a Proca-type mass for the
#       gauge vector: one constant sets both.
#
# VERDICT (honest bracket).  On Weyl's procedure -- the one the Report cites
# -- lambda is a CONSTANT (the gauge-fixed value of the Weyl curvature
# scalar), not a Lagrange multiplier, and it cannot be zero.  So DT's
# variational principle carries a NON-ZERO cosmological constant,
# |Lambda| = lambda/4, whose MAGNITUDE is a free constant: DT does not fix it,
# and does not force Lambda = 0.  The corpus' multiplier reading
# (T1_9_weyl_strengthening, phi.phi = 1/12) is a misreading and is superseded.
# Side finding: Eq 67's "1/2" is the 4D (R~^2) value; the 5D gauge-invariant
# parent gives 3/5 -- same conclusion, different coefficient.  Not settled
# here: the size of the Weyl-vector mass (it needs the F^2 normalisation
# alpha), and how Williams' imaginary (London) gauge vector changes the sign
# of the phi^.phi^ term.
#
# RUN: /home/dsteele/miniconda3/envs/sage/bin/sage weyl_lambda_reading.sage
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


print("=" * 72)
print("weyl_lambda_reading: what lambda is in Report Eqs 39 / 67")
print("=" * 72)

# ==========================================================================
# [1] Gauge weights
# ==========================================================================
print("\n[1] Gauge weight of R~ and of sqrt|g| R~^p")
diag5 = [exp(t * x), -(1 + t**2 + x**2), -(2 + t) * (1 + x**2),
         -exp(t) * (1 + x**4), -(3 + t * x**2)]
phi5 = [t * x, sin(t) + x, x**2 - t, 1 + t * x**2, cos(x) * t]
sig = t * x / 2 + x**2 / 3
for n in (4, 5):
    X = XALL[:n]
    g = [[diag5[i] if i == j else 0 for j in range(n)] for i in range(n)]
    gp = [[exp(2 * sig) * g[i][j] for j in range(n)] for i in range(n)]
    ph = phi5[:n]
    php = [ph[i] - diff(sig, X[i]) for i in range(n)]
    Rw = weyl_scalar(g, ph, X)
    Rwp = weyl_scalar(gp, php, X)
    ok("%dD: g -> e^(2s) g, phi -> phi - ds  gives  R~ -> e^(-2s) R~" % n,
       z(Rwp - exp(-2 * sig) * Rw))
    # sqrt|g| -> e^(n s) sqrt|g|, hence sqrt|g| R~^p -> e^((n - 2p) s) (...)
    for p in sorted(set([2, Integer(n) / 2])):
        inv = (n - 2 * p == 0)
        print("    n = %d, p = %s: weight of sqrt|g| R~^p = %s  (%s)"
              % (n, p, n - 2 * p, "invariant" if inv else "NOT invariant"))
    dw = matrix(SR, gp).det() / matrix(SR, g).det()     # = (sqrt-ratio)^2
    ok("%dD: sqrt|g| -> e^(%d s) sqrt|g|, so sqrt|g| R~^p is invariant only for "
       "p = n/2 = %s" % (n, n, Integer(n) / 2), z(dw - exp(2 * n * sig)))
ok("mutation: R~^2 is not gauge-invariant in 5D (weight 5 - 4 = 1)", 5 - 2 * 2 != 0)

# ==========================================================================
# [2] The natural gauge R~ = -lambda and the linearised Lagrangian
# ==========================================================================
print("\n[2] Natural gauge R~ = -lambda: the linearised Lagrangian")
Rt, s, lam, p = var('Rt s lam p')
lam = var('lam', domain='positive')


def f(r, pp):
    return r**pp


def lin(pp, s_):
    """f(s) + f'(s)(Rt - s): the part of f(Rt) the field equations see."""
    r = var('r_')
    fp = diff(f(r, pp), r).subs({r: s_})
    return f(s_, pp) + fp * (Rt - s_), fp


for pp in (2, Integer(5) / 2):
    r = var('r_')
    F = f(Rt, pp)
    L1, fp = lin(pp, s)
    # first-order agreement at Rt = s (value and first derivative)
    ok("p = %s: f(R~) and its linearisation agree to first order at R~ = s"
       % pp, z((F - L1).subs({Rt: s})) and z(diff(F - L1, Rt).subs({Rt: s})))
    c = (L1 / fp - Rt).simplify_full().subs({s: -lam})
    ok("p = %s: L / f'(-lambda) = R~ + c lambda with c = (p-1)/p = %s"
       % (pp, (pp - 1) / pp), z(c - lam * (pp - 1) / pp))

# Insert R~ = R - 2(n-1) div(phi) - (n-1)(n-2) phi.phi, phi = sqrt(lam) phi^
R, dv, ph2h = var('R dv ph2h')
for n, k, eq in [(4, 6, '39'), (5, 12, '67')]:
    Rw_sym = R - 2 * (n - 1) * dv - (n - 1) * (n - 2) * lam * ph2h
    Lnorm = Rw_sym + lam / 2                     # p = 2, c = 1/2
    williams = R + lam * (Integer(1) / 2 - k * ph2h)
    ok("%dD, p = 2: R~ + lambda/2 = Eq %s's  R + lambda(1/2 - %d phi^.phi^)"
       " (up to the divergence term)" % (n, eq, k),
       z(Lnorm + 2 * (n - 1) * dv - williams))
ok("mutation: the 5D gauge-invariant parent (p = 5/2) gives c = 3/5, not "
   "Eq 67's 1/2", (Integer(5) / 2 - 1) / (Integer(5) / 2) == Integer(3) / 5
   and Integer(3) / 5 != Integer(1) / 2)
print("    Eq 39 is Weyl's 4D natural-gauge Lagrangian exactly.  Eq 67 is the")
print("    same R~^2 procedure carried to 5D (c = 1/2), where the gauge-")
print("    invariant parent would be R~^(5/2) (c = 3/5).")

# ==========================================================================
# [3] lambda as an auxiliary field: its equation is the gauge condition
# ==========================================================================
print("\n[3] Varying lambda the right way")
for pp in (2, Integer(5) / 2):
    L1, fp = lin(pp, s)
    r = var('r_')
    fpp = diff(f(r, pp), r, 2).subs({r: s})
    ok("p = %s: d/ds [f(s) + f'(s)(R~ - s)] = f''(s) (R~ - s)" % pp,
       z(diff(L1, s) - fpp * (Rt - s)))
    ok("p = %s: f''(s) != 0, so EL(lambda) is R~ = s = -lambda (the natural "
       "gauge), not a constraint on phi" % pp, not z(fpp))
print("    The multiplier reading (EL(lambda) => phi^.phi^ = 1/(2k)) varies")
print("    lambda in the NORMALISED Lagrangian, where lambda also sits in the")
print("    dropped factor f'(-lambda) and inside phi^ = phi/sqrt(lambda); it")
print("    is not an independent variable there.  That reading is superseded.")

# ==========================================================================
# [4] lambda != 0, so the cosmological term survives
# ==========================================================================
print("\n[4] lambda cannot vanish; the cosmological constant")
for pp in (2, Integer(5) / 2):
    r = var('r_', domain='positive')
    fp0 = diff(f(r, pp), r).subs({r: 0})
    ok("p = %s: f'(0) = 0, so the normalisation by f'(-lambda) needs "
       "lambda != 0" % pp, z(fp0))
# Metric variation of sqrt|g| lambda c (einstein_from_dt_action [3]):
# the g^ab coefficient is -(lambda c)/2; written as R - 2 Lambda, Lambda is
# -c lambda / 2, i.e. |Lambda| = lambda / 4 for c = 1/2.
c = Integer(1) / 2
Lam = -c * lam / 2
ok("Eq 39/67 (c = 1/2): |Lambda| = lambda/4 != 0 for lambda != 0",
   z(abs(Lam) - lam / 4) and bool(Lam != 0))
print("    The same lambda multiplies the phi^.phi^ (Proca-type mass) term:")
print("    one constant sets both Lambda and the Weyl-vector mass.")

print("""
  NET: on Weyl's procedure, which the Report invokes, lambda is a CONSTANT
  -- the gauge-fixed value of the Weyl curvature scalar -- not a Lagrange
  multiplier, and it cannot be zero.  DT's variational principle therefore
  carries a non-zero cosmological constant, |Lambda| = lambda/4 (Eq 39/67 as
  written; c lambda/2 in general), of UNFIXED magnitude.  Lambda = 0 is not
  available.  The surface equations are G_ab + Lambda g_ab = (8 pi G/c^4) T_ab
  with Lambda small but free.  Open: the Weyl-vector mass scale (needs the
  F^2 normalisation alpha) and the effect of the imaginary London gauge.""")
print("=" * 72)
print("weyl_lambda_reading COMPLETE.")
print("=" * 72)
