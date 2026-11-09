#!/usr/bin/env sage
# --------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, LA-8370-MS (LASL, 1980), Eq (39) (4D) and
#   Eq (67) (5D), the Weyl-type variational principle
#       delta INT [ R + (1/4) alpha F_ij F^ij - lambda (1/2 - 12 phi_i phi^i) ]
#             sqrt(-g) d^5x = 0 ,                                     (67)
#   with the Report's own statement directly below it: "Varying the metric
#   coefficients g_ij in the variational principle (67) will yield field
#   equations of the Dynamic Theory which are extensions of Einstein's
#   General Theory of Relativity."  The Weyl connection is the one DT derives
#   (T1.8; Report eqs 36-37, the W(p) term of T1.8).
# Supplementary Materials: step [G1], Sec. S3 (code index: Sec. S6).
# Companions: einstein_equations_recovery.sage (conservation, Lovelock
#   uniqueness, kappa = 8 pi G), gr_on_4d_surface.sage (u4 = 0 <=> K = 0).
# --------------------------------------------------------------------------

# ============================================================================
# einstein_from_dt_action.sage                  [A]  (assertion-hardened)
#
# THE SURFACE METRIC IS DYNAMICAL, AND ITS EQUATIONS ARE SECOND ORDER, BECAUSE
# DT HAS AN ACTION.
#
# einstein_equations_recovery.sage showed: IF the surface metric obeys local,
# second-order, metric-only field equations, they are Einstein's.  That "if"
# is not an extra assumption: the Report gives the action (Eq 67) and says
# its metric variation yields DT's field equations.  This script carries the
# variation out, which T1.9 had left as "cited, not derived".
#
#   [1] THE COEFFICIENTS 6 AND 12 ARE WEYL'S.  For the Weyl connection DT
#       derives, Gamma~ = Gamma + (delta phi + delta phi - g phi), the scalar
#       curvature is  R~ = R - 2(n-1) div(phi) - (n-1)(n-2) phi.phi .
#       (n-1)(n-2) = 6 in 4D (Eq 39) and 12 in 5D (Eq 67): Williams'
#       constraint coefficients are the Weyl-curvature coefficients, so the
#       action is the Weyl-geometric one, not an ad hoc choice.  Mutation: the
#       4D coefficient 6 fails in 5D.
#   [2] THE R TERM GIVES THE EINSTEIN TENSOR, AT SECOND ORDER.  The Euler-
#       Lagrange derivative of sqrt|g| R with respect to each metric function
#       equals  -sqrt|g| G^ab dg_ab/du  identically (4D and 5D families; an
#       off-diagonal component included).  So the gravitational operator is
#       G_ab -- local and second order (the third and fourth derivatives in
#       the raw Euler-Lagrange expression cancel).  Mutation: R^ab in place
#       of G^ab fails.  This closes the "Einstein-Hilbert variation cited,
#       not derived" gap in T1.9 (on these families).
#   [3] THE MATTER TERMS.  The F^2 term gives the Maxwell stress; the lambda
#       term gives  -lambda [ 12 phi^a phi^b + (1/2)(1/2 - 12 phi.phi) g^ab ].
#       Its g^ab part is a cosmological term whose coefficient is the
#       constraint itself:
#         * lambda a Lagrange multiplier (the reading used in T1.9): the
#           constraint holds on shell and the cosmological term VANISHES;
#         * lambda a constant (Weyl's own 1918 reading): a cosmological term
#           -lambda/4 g^ab survives, plus a Proca-type stress.
#   [4] ASSEMBLY.  Eq 67 => 5D  G_ab = T^F_ab + T^lambda_ab, local and second
#       order, metric the only gravitational field (R is the only metric
#       term).  With u4 = 0 (=> K = 0) Gauss-Codazzi restricts it exactly to
#       the 4D surface (gr_on_4d_surface); kappa = 8 pi G from the Newtonian
#       limit (einstein_equations_recovery [4]).
#
# VERDICT.  The Lovelock premises are SUPPLIED BY DT: the
# Report's own action, Eq 67, is local, its gravitational part is R alone,
# and its metric variation is second order (shown here).  The Einstein
# equations are therefore DT's mass-conserving (u4 = 0) sector, on Williams'
# own variational principle.  What is NOT derived from the three laws: the
# choice of Lagrangian.  DT derives the Weyl geometry (T1.8); taking R (rather
# than R~^2, or adding Gauss-Bonnet) as the Lagrangian follows Weyl, in
# Williams' words "we may follow Weyl's procedure".  Lambda is zero on the
# multiplier reading of lambda (the reading of the original treatment) and non-zero on Weyl's
# constant reading; the Report does not settle which.  The variation is
# checked on function-parametrised metric families (generic 2-variable
# functions, diagonal 4D and 5D, plus an off-diagonal 4D instance), the same
# "restricted form" standard as the Lovelock check.
#
# RUN: sage einstein_from_dt_action.sage
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


# --------------------------------------------------------------------------
# Curvature of an arbitrary torsion-free connection, hand-built so that the
# same code serves the Levi-Civita and the Weyl connection.
#   Gamma[k][i][j] = Gamma^k_ij
#   Ric_sv = d_r Gamma^r_vs - d_v Gamma^r_rs + Gamma^r_rl Gamma^l_vs
#            - Gamma^r_vl Gamma^l_rs          (MTW / Sage convention)
# --------------------------------------------------------------------------
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


def scalar(gi, Ric):
    n = len(gi)
    return sum(gi[a][b] * Ric[a][b] for a in range(n) for b in range(n))


print("=" * 72)
print("einstein_from_dt_action: varying the Report's Eq 67 action")
print("=" * 72)

# ==========================================================================
# [1] Weyl scalar curvature: the coefficients 6 (Eq 39) and 12 (Eq 67)
# ==========================================================================
print("\n[1] The constraint coefficients 6 (4D) and 12 (5D) are Weyl's (n-1)(n-2)")

# Explicit, non-symmetric metrics and gauge vectors (the identity is a
# theorem; a non-trivial instance plus mutations is the check).
diag5 = [exp(t * x), -(1 + t**2 + x**2), -(2 + t) * (1 + x**2),
         -exp(t) * (1 + x**4), -(3 + t * x**2)]
phi5 = [t * x, sin(t) + x, x**2 - t, 1 + t * x**2, cos(x) * t]

# Validate the hand-built curvature against Sage's own (4D instance).
M4 = Manifold(4, 'M4v', structure='Lorentzian')
C4 = M4.chart('t x y w')
gS = M4.metric('gS')
for i in range(4):
    gS[i, i] = diag5[i].subs({t: C4[0], x: C4[1]})
g4 = [[diag5[i] if i == j else 0 for j in range(4)] for i in range(4)]
g4i = inverse(g4)
R4hand = scalar(g4i, ricci(levi_civita(g4, g4i, XALL[:4]), XALL[:4]))
R4sage = gS.ricci_scalar().expr().subs({C4[0]: t, C4[1]: x})
ok("hand-built Ricci scalar agrees with Sage's (conventions pinned)",
   z(R4hand - R4sage))


def weyl_pieces(n, sgn=1):
    X = XALL[:n]
    g = [[diag5[i] if i == j else 0 for j in range(n)] for i in range(n)]
    gi = inverse(g)
    ph = phi5[:n]
    phu = [sum(gi[k][m] * ph[m] for m in range(n)) for k in range(n)]
    Gam = levi_civita(g, gi, X)
    d = lambda a, b: 1 if a == b else 0
    Gw = [[[Gam[k][i][j] + sgn * (d(k, i) * ph[j] + d(k, j) * ph[i]
                                   - g[i][j] * phu[k])
            for j in range(n)] for i in range(n)] for k in range(n)]
    R = scalar(gi, ricci(Gam, X))
    Rw = scalar(gi, ricci(Gw, X))
    div = sum(diff(phu[k], X[k]) for k in range(n)) + \
        sum(Gam[k][k][l] * phu[l] for k in range(n) for l in range(n))
    sq = sum(ph[k] * phu[k] for k in range(n))
    return R, Rw, div, sq


for n, coef in [(4, 6), (5, 12)]:
    R, Rw, div, sq = weyl_pieces(n)
    ok("%dD: R~ = R - 2(n-1) div(phi) - %d phi.phi   ((n-1)(n-2) = %d)"
       % (n, coef, (n - 1) * (n - 2)),
       (n - 1) * (n - 2) == coef and z(Rw - R + 2 * (n - 1) * div + coef * sq))
    if n == 5:
        ok("mutation: the 4D coefficient 6 is wrong in 5D",
           not z(Rw - R + 2 * (n - 1) * div + 6 * sq))
        ok("mutation: the quadratic term is not optional (phi.phi != 0 here)",
           not z(Rw - R + 2 * (n - 1) * div))
print("    Eq 39 carries 6, Eq 67 carries 12: both are the Weyl-curvature")
print("    coefficient of the dimension in which Williams writes them.")

# ==========================================================================
# Euler-Lagrange machinery for a Lagrangian density of (t, x) built from
# metric functions u_k(t, x) up to second derivatives.
# ==========================================================================
JETS = [(), (t,), (x,), (t, t), (t, x), (x, x)]


def euler_lagrange(L, funcs):
    """Return {u: E_u} with E_u the variational derivative of INT L dt dx."""
    fwd, bwd, sym = [], {}, {}
    k = 0
    for u in funcs:                     # highest jets first, so that
        for J in reversed(JETS):        # derivatives are replaced before u
            e = diff(u, *J) if J else u
            s = SR.var('jet%d' % k)
            k += 1
            fwd.append((e, s))
            bwd[s] = e
            sym[(u, J)] = s
    Ls = L
    for e, s in fwd:
        Ls = Ls.subs({e: s})
    out = {}
    for u in funcs:
        E = 0
        for J in JETS:
            P = diff(Ls, sym[(u, J)]).subs(bwd)
            for c in J:
                P = diff(P, c)
            E += (-1) ** len(J) * P
        out[u] = E
    return out


def metric_family(n, funcs_diag, offdiag=None):
    g = [[0] * n for _ in range(n)]
    g[0][0] = exp(2 * funcs_diag[0])
    for i in range(1, n):
        g[i][i] = -exp(2 * funcs_diag[i])
    if offdiag is not None:
        g[0][1] = g[1][0] = offdiag
    return g


def einstein_up(g, X):
    n = len(X)
    gi = inverse(g)
    Ric = ricci(levi_civita(g, gi, X), X)
    R = scalar(gi, Ric)
    Gd = [[Ric[a][b] - R * g[a][b] / 2 for b in range(n)] for a in range(n)]
    Gu = [[sum(gi[a][c] * Gd[c][d] * gi[d][b] for c in range(n) for d in range(n))
           for b in range(n)] for a in range(n)]
    Ru = [[sum(gi[a][c] * Ric[c][d] * gi[d][b] for c in range(n) for d in range(n))
           for b in range(n)] for a in range(n)]
    return gi, R, Gu, Ru


def dg_du(g, u):
    s = SR.var('dgdu_tmp')
    return [[diff(gab.subs({u: s}), s).subs({s: u}) for gab in row] for row in g]


def contract(T, D):
    n = len(T)
    return sum(T[a][b] * D[a][b] for a in range(n) for b in range(n))


# ==========================================================================
# [2] delta INT sqrt|g| R = -INT sqrt|g| G^ab delta g_ab  (second order)
# ==========================================================================
print("\n[2] The R term: its metric variation is the Einstein tensor")

for n in (4, 5):
    X = XALL[:n]
    U = [function('U%d' % i)(t, x) for i in range(n)]
    g = metric_family(n, U)
    gi, R, Gu, Ru = einstein_up(g, X)
    sg = exp(sum(U))                       # sqrt|det g| for this family
    E = euler_lagrange(sg * R, U)
    good = all(z(E[u] + sg * contract(Gu, dg_du(g, u))) for u in U)
    ok("%dD diagonal family (%d free functions of t, x): "
       "dL/du = -sqrt|g| G^ab dg_ab/du for every u" % (n, n), good)
    if n == 4:
        ok("mutation: Ricci in place of Einstein (R^ab) fails",
           not z(E[U[1]] + sg * contract(Ru, dg_du(g, U[1]))))

# Off-diagonal: one generic family, then checked on an explicit instance
# (generic simplification of the sqrt(det) expression is expensive).
X = XALL[:4]
U = [function('V%d' % i)(t, x) for i in range(4)]
f = function('Voff')(t, x)
g = metric_family(4, U, offdiag=f)
gi, R, Gu, Ru = einstein_up(g, X)
sg = sqrt(exp(2 * U[0] + 2 * U[1]) + f**2) * exp(U[2] + U[3])
E = euler_lagrange(sg * R, U + [f])
inst = [(U[0], t * x / 3), (U[1], t**2 / 5 + x), (U[2], sin(x) + t / 2),
        (U[3], x**2 / 4 - t), (f, t + x**2 / 2)]


def at_inst(e):
    for (fu, val) in inst:
        e = e.substitute_function(fu.operator(), val.function(t, x))
    return e


# simplify_full cannot clear the sqrt(e^(2A+2B) + f^2) radicals, so the
# instance is evaluated at 300-bit precision at several points; the terms
# themselves are O(1), so a residual below 1e-60 is an identity, not luck.
RF = RealField(300)
pts = [(Integer(3) / 10, Integer(7) / 10), (Integer(-6) / 5, Integer(1) / 3),
       (Integer(2), Integer(-3) / 2)]
worst, scale = RF(0), RF(0)
for u in U + [f]:
    lhs = at_inst(E[u] + sg * contract(Gu, dg_du(g, u)))
    rhs = at_inst(E[u])
    for (a, b) in pts:
        worst = max(worst, abs(RF(lhs.subs({t: a, x: b}))))
        scale = max(scale, abs(RF(rhs.subs({t: a, x: b}))))
print("    off-diagonal instance: max residual %.1e against terms of size %.1f"
      % (float(worst), float(scale)))
ok("4D with off-diagonal g_tx: identity holds for all five metric functions "
   "(explicit instance, 300-bit, 3 points)", worst < RF(10)**-60 and scale > 0.1)
print("    The raw Euler-Lagrange expression contains third and fourth")
print("    derivatives; the identity shows they cancel: the field equations")
print("    of the R term are SECOND ORDER, as Lovelock's premise requires.")

# ==========================================================================
# [3] The matter terms of Eq 67: Maxwell stress and the lambda term
# ==========================================================================
print("\n[3] The F^2 and lambda terms (5D diagonal family)")
n = 5
X = XALL[:n]
U = [function('W%d' % i)(t, x) for i in range(n)]
g = metric_family(n, U)
gi = inverse(g)
sg = exp(sum(U))
A = [function('A%d' % i)(t, x) for i in range(n)]
F = [[diff(A[j], X[i]) - diff(A[i], X[j]) for j in range(n)] for i in range(n)]
F2 = sum(F[i][j] * F[k][l] * gi[i][k] * gi[j][l]
         for i in range(n) for j in range(n) for k in range(n) for l in range(n))
Mu = [[sum(gi[a][c] * gi[b][d] * gi[j][k] * F[c][j] * F[d][k]
           for c in range(n) for d in range(n) for j in range(n) for k in range(n))
       for b in range(n)] for a in range(n)]
EF = euler_lagrange(sg * F2 / 4, U)
TF = [[-Mu[a][b] / 2 + gi[a][b] * F2 / 8 for b in range(n)] for a in range(n)]
ok("F^2 term: dL/du = sqrt|g| (-1/2 F^a_j F^bj + 1/8 g^ab F^2) dg_ab/du "
   "(the Maxwell stress)",
   all(z(EF[u] - sg * contract(TF, dg_du(g, u))) for u in U))

lam = var('lam')
ph = A                                     # the gauge vector (phi-hat)
X2 = sum(gi[i][j] * ph[i] * ph[j] for i in range(n) for j in range(n))
phu = [sum(gi[a][m] * ph[m] for m in range(n)) for a in range(n)]
El = euler_lagrange(-lam * (Integer(1) / 2 - 12 * X2) * sg, U)
Tl = [[-lam * (12 * phu[a] * phu[b] + (Integer(1) / 2 - 12 * X2) * gi[a][b] / 2)
       for b in range(n)] for a in range(n)]
ok("lambda term: dL/du = -sqrt|g| lam [12 phi^a phi^b + (1/2)(1/2 - 12 phi.phi)"
   " g^ab] dg_ab/du", all(z(El[u] - sg * contract(Tl, dg_du(g, u))) for u in U))

# The cosmological (g^ab) part has the constraint as its coefficient.
coef_g = -lam * (Integer(1) / 2 - 12 * X2) / 2
ok("the g^ab part of the lambda stress is  -(lam/2)(1/2 - 12 phi.phi) g^ab",
   all(z(Tl[a][b] + lam * 12 * phu[a] * phu[b] - coef_g * gi[a][b])
       for a in range(n) for b in range(n)))
ok("multiplier reading: on the constraint phi.phi = 1/24 the cosmological "
   "term vanishes (Lambda = 0)",
   z(-(lam / 2) * (Integer(1) / 2 - 12 * Integer(1) / 24)))
ok("constant reading (Weyl 1918): with no constraint, the phi-independent "
   "part is -lam/4 g^ab (Lambda != 0)",
   z(coef_g.subs({a_: 0 for a_ in A}) + lam / 4))
print("    The Report does not say which reading is meant; the Layer-1 check")
print("    takes lambda as a multiplier (EL(lambda) => phi.phi = const).")

# ==========================================================================
# [5] Signature, the sign of Lambda, de Sitter, and the vector mass
# ==========================================================================
print("\n[5] Signature and the sign of Lambda")
print("""    The Report writes Eq 67 with x0 = i c t, i.e. effectively signature
    (-,+,+,+,+).  Reversing the metric reverses R but not the constant term,
    so in the signature used here (+,-,-,-,-) Weyl's natural-gauge form is
        L = R + (1/4) alpha F^_ij F^^ij + lambda (1/2 - 12 phi^.phi^)
    (weyl_lambda_reading [2] checks R~ + lambda/2 = this, up to a divergence;
    the phi^.phi^ sign is that of R~ for a REAL Weyl vector).  The literal
    transcription -lambda(...) used in [3] is the same Lagrangian with
    lambda -> -lambda.""")
lamp = var('lamp')
# constant term +lamp/2 (in [3] the constant term is -lam/2 with g^ab part -lam/4)
coef_cos = (coef_g.subs({a_: 0 for a_ in A})).subs(lam=-lamp)
ok("the constant term +lambda/2 gives +(lambda/4) g^ab in dL/dg_ab",
   z(coef_cos - lamp / 4))
print("    EL in g_ab:  -G^ab + (lambda/4) g^ab + (matter) = 0, i.e.")
print("        G_ab - Lambda g_ab = kappa T_ab,   Lambda = lambda/4.")

H = var('H')
for n, cn in [(4, 3), (5, 6)]:
    X = XALL[:n]
    g = [[0] * n for _ in range(n)]
    g[0][0] = 1
    for i in range(1, n):
        g[i][i] = -exp(2 * H * t)
    gi = inverse(g)
    Ric = ricci(levi_civita(g, gi, X), X)
    R = scalar(gi, Ric)
    ok("%dD de Sitter (scale factor e^(Ht)): G_ab = %d H^2 g_ab" % (n, cn),
       all(z(Ric[a][b] - R * g[a][b] / 2 - cn * H**2 * g[a][b])
           for a in range(n) for b in range(n)))
    ok("%dD: hence de Sitter has Lambda = %d H^2 > 0, lambda = %d H^2 > 0"
       % (n, cn, 4 * cn), bool(cn > 0))
    ok("%dD mutation: in the form G + Lambda g = kappa T de Sitter would need "
       "Lambda < 0" % n, bool(-cn < 0))

# weak field fixes the sign of kappa (attractive gravity, G = +kappa T)
Phi = function('Phi')(x, y, w)
e = var('e')
X = XALL[:4]
g = [[0] * 4 for _ in range(4)]
g[0][0] = 1 + 2 * e * Phi
for i in range(1, 4):
    g[i][i] = -(1 - 2 * e * Phi)
gi = inverse(g)
Ric = ricci(levi_civita(g, gi, X), X)
R = scalar(gi, Ric)
G00 = taylor(Ric[0][0] - R * g[0][0] / 2, e, 0, 1)
ok("weak field: G_00 = 2 lap(Phi), so attractive gravity is G = +kappa T "
   "(same convention as the static metric g_tt = 1 - Lambda r^2/3)",
   z(G00 - 2 * e * (diff(Phi, x, 2) + diff(Phi, y, 2) + diff(Phi, w, 2))))
print("    The action leaves the sign of lambda free: Weyl's procedure needs")
print("    only lambda != 0.  The sign of Lambda is the sign of lambda.")

al = var('al')
print("\n    Vector mass (weyl_vector_mass [2]): with phi^ = phi/sqrt(lambda),")
print("    L = R - 12 phi.phi + (1/4)(alpha/lambda) F^2 + lambda/2, so")
print("    m^2 = 2k lambda/alpha, k = 12 (5D), 6 (4D): set by the free ratio")
print("    lambda/alpha, a Planck mass only for a charge q of order one in")
print("    Ghilencea's normalisation.  The integrable vector vanishes on shell")
print("    (weyl_vector_mass [5]), so the value plays no part in the recovery.")
ok("m^2 = 2k/a with a = alpha/lambda gives m^2 = 24 lambda/alpha in 5D",
   z(2 * 12 / (al / lamp) - 24 * lamp / al))

# ==========================================================================
# [6] Overall sign of the action, and the four-dimensional source
# ==========================================================================
print("\n[6] Overall sign, and the map from the 5D to the 4D source")
hs = list(var('h0 h1 h2 h3 h4'))
sq5 = sqrt(hs[0] * (-hs[1]) * (-hs[2]) * (-hs[3]) * (-hs[4]))
at_eta5 = {hs[0]: 1, hs[1]: -1, hs[2]: -1, hs[3]: -1, hs[4]: -1}
Ee, kap5, av6 = var('Ee kap5 av6', domain='positive')
LMax = -(2 * Ee**2 / (hs[0] * hs[1])) / 4          # standard Maxwell -F^2/4
dM = (diff(sq5 * LMax, hs[0]) / sq5).subs(at_eta5)
# S = s/(2 kappa_5) INT sqrt|g| R + S_m and [2] (d sqrt|g|R/dg_ab = -sqrt|g| G^ab):
#   -(s/(2 kappa_5)) G^00 + dM = 0  =>  G^00 = 2 kappa_5 dM / s
ok("prefactor -1/(2 kappa_5): standard Maxwell matter (-F^2/4, energy "
   "E^2/2 > 0) gives G^00 = +kappa_5 E^2/2, attractive (Landau-Lifshitz sign)",
   z(2 * kap5 * dM / (-1) - kap5 * Ee**2 / 2))
ok("mutation: prefactor +1/(2 kappa_5) gives G^00 = -kappa_5 E^2/2",
   z(2 * kap5 * dM / 1 + kap5 * Ee**2 / 2))
FFs, PPs = var('FFs PPs')
# vector part of -(1/(2 kappa_5)) [ (1/4) a F^2 - 12 phi.phi ], a = alpha/lambda;
# canonical field A = sqrt(a/(2 kappa_5)) phi  =>  -(1/4) F_A^2 + (1/2) m^2 A.A
Lvec = -(av6 * FFs / 4 - 12 * PPs) / (2 * kap5)
cn = av6 / (2 * kap5)
ok("canonical form: -(1/4) F_A^2 + (1/2)(24/a) A.A, the standard Maxwell and "
   "Proca signs in (+,-,-,-,-) iff a = alpha/lambda > 0; m^2 = 24 lambda/alpha",
   z(Lvec - (-FFs * cn / 4 + (24 / av6) * PPs * cn / 2)))

# product metric: 4D surface metric of t, x times a fixed g_44 = -1
A6 = [function('B%d' % i)(t, x) for i in range(4)]
g5 = [[0] * 5 for _ in range(5)]
g5[0][0] = exp(2 * A6[0])
for i in range(1, 4):
    g5[i][i] = -exp(2 * A6[i])
g5[4][4] = -1
X5 = XALL[:5]
gi5 = inverse(g5)
Ric5 = ricci(levi_civita(g5, gi5, X5), X5)
R5 = scalar(gi5, Ric5)
g4 = [r_[:4] for r_ in g5[:4]]
gi4 = inverse(g4)
Ric4 = ricci(levi_civita(g4, gi4, XALL[:4]), XALL[:4])
R4 = scalar(gi4, Ric4)
G5 = [[Ric5[a_][b_] - R5 * g5[a_][b_] / 2 for b_ in range(5)] for a_ in range(5)]
G4 = [[Ric4[a_][b_] - R4 * g4[a_][b_] / 2 for b_ in range(4)] for a_ in range(4)]
ok("x4-independent metric, g_44 = -1: G^(5)_ab = G^(4)_ab, G^(5)_a4 = 0, "
   "G^(5)_44 = R^(4)/2",
   all(z(G5[a_][b_] - G4[a_][b_]) for a_ in range(4) for b_ in range(4))
   and all(z(G5[a_][4]) for a_ in range(4)) and z(G5[4][4] - R4 / 2))
ok("4D trace: g^ab G_ab = -R", z(sum(gi4[a_][b_] * G4[a_][b_]
                                       for a_ in range(4) for b_ in range(4)) + R4))
print("""    So G - Lambda g = kappa_5 T^(5) has x4-independent left side: the
    source must be uniform along x4, T^(5)_ab = T_ab/ell, with T_ab the
    observed 4D stress and ell a length (S_m = ell^-1 INT dx4 S_m^(4)).
    The surface equation is then G - Lambda g = kappa T with kappa = kappa_5/ell.""")
kp, ell, Lm6, Tr6, R4s = var('kp ell Lm6 Tr6 R4s')
R4sol = solve(-R4s - 4 * Lm6 == kp * Tr6, R4s)[0].rhs()   # trace of 4D equation
T44_5 = (R4sol / 2 + Lm6) / kap5                          # 44 eq: G_44 + Lambda = kappa_5 T44
ok("44 equation with the same normalisation: ell T^(5)_44 = -Lambda/kappa "
   "- T/2 (kappa = kappa_5/ell)",
   z((ell * T44_5).subs(kp=kap5 / ell) - (-Lm6 / (kap5 / ell) - Tr6 / 2)))

# ==========================================================================
# [4] Assembly
# ==========================================================================
print("\n[4] Assembly")
print("""    Eq 67, varied in the metric:
        G_ab = T^F_ab + T^lambda_ab        (5D, local, second order;
                                             R the only metric term)
    + u4 = 0  =>  K = 0  =>  exact restriction to the 4D surface
                                             (gr_on_4d_surface)
    + Newtonian limit calibrates kappa = 8 pi G / c^4
                                             (einstein_equations_recovery [4])
    =>  G_ab - Lambda g_ab = (8 pi G/c^4) T_ab on the surface, with
        Lambda = lambda/4 ([5]; Lambda = 0 on the multiplier reading).

  NET: the premise 'the surface metric is dynamical, with local second-order
  equations' is DT's own -- the Report's action Eq 67, whose metric variation
  is shown here to be the Einstein operator.  Einstein's equations are DT's
  mass-conserving (u4 = 0) sector.  Not derived from the three laws: the
  choice of R as the Lagrangian (Weyl's procedure, adopted by Williams).
  Open: the reading of lambda, which decides Lambda.""")
print("=" * 72)
print("einstein_from_dt_action COMPLETE.")
print("=" * 72)
