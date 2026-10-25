#!/usr/bin/env sage
# --------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, "The Dynamic Theory" (2011), Eqs (1.34)-
#   (1.39) / (10.1)-(10.9); u4 = dy4/dt = 0 is Eq 1.36 / 10.6 (conservation
#   of mass).  Report LA-8370-MS eqs 163-164 (induced metric).
# Context: gr_on_4d_surface.sage (Gauss-Codazzi in Gaussian-normal form,
#   g_a4 = 0 and g_44 = -1 built in), einstein_equations_recovery.sage (G2),
#   kk_radion_g44.sage (Cassini selects g_44 fixed).
# Supplementary Materials: step [G2], Sec. S3 (code index: Sec. S6).
# --------------------------------------------------------------------------
#
# ============================================================================
# gr_projection_audit.sage                        [A]  (assertion-hardened)
#
# AUDIT OF THE PROJECTION STEP IN G2.  Does the tangential projection of the
# 5D Einstein tensor equal the 4D Einstein tensor of the surface x4 = const
# for EVERY metric DT admits, or only for the Gaussian-normal ansatz used in
# gr_on_4d_surface?  In a general 5D geometry K_ab = 0 is not enough: the
# projection also carries the normal ("electric") curvature
#       E_ab = R5(n, e_a, n, e_b) ,
# and that contains d_n K_ab and the Hessian of the lapse.  Here the slices
# are the PHYSICAL surfaces x4 = y = const, and the metric is general:
#       ds^2 = q_ab(x, y) dx^a dx^b + 2 g_a4 dx^a dy + g_44 dy^2 .
#
#   [A] Gauss (contracted) for a general lapse N (g_44 = eps N^2), checked on
#       an abstract non-diagonal family of functions of (t, x, y):
#         G5_ab - G4_ab = eps[K K_ab - K_ac K^c_b]
#                         + eps E_ab ... (coefficients fixed by the code)
#       i.e. the residual is EXACTLY the K-terms plus the E-terms.
#   [B] the normal curvature:  E_ab = (1/N) d_y K_ab + K K-terms
#       + (1/N) D_a D_b N  (coefficients fixed by the code).
#   [C] SUFFICIENCY: K = 0 on EVERY slice, zero shift, constant g_44 =>
#       residual = 0 identically, for arbitrary q_ab(x).  The class DT
#       admits under its stated conditions is exactly this one.
#   [D] NECESSITY, one counter-example per dropped condition:
#       (D1) K = 0 on one slice only: d_y K survives;
#       (D2) g_44 = -N(x)^2 non-constant (a radion): the Hessian term
#            (D_a D_b N - q_ab box N)/N survives; it vanishes iff
#            D_a D_b N = 0;
#       (D3) g_a4 != 0: a generic shift already gives K != 0 (free fall
#            leaves the slice, so u4 = 0 excludes it); a Killing shift with
#            dA != 0 keeps K = 0 but leaves a residual.  So g_a4 = 0 is an
#            INDEPENDENT condition, supplied by the restriction.
#   [E] the remaining components under the conditions: G5_a4 = 0 and
#       G5_44 = -(eps/2) R4 (the T_a4, T_44 constraints).
# ============================================================================

from sage.all import *


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


def is_zero(e):
    return bool(SR(e).simplify_full() == 0)


def mat_zero(Mx):
    return all(is_zero(e) for e in Mx.list())


eps = -1
t, x, u, v, y = var('t x u v y')
crd4 = [t, x, u, v]
_n = [0]


def geometry(qm, g4col=None, N=1):
    """5D metric with blocks qm (4x4), g_a4 = g4col, g_44 = eps N^2 (lapse
    N, default 1).  Slices y = const; e_a = d_a are tangent to them, so the tangential
    projection of a 5D tensor is its (a, b) block.  Returns the residual
    Delta = G5_ab - G4_ab[q] and, for zero shift, the Gauss pieces."""
    _n[0] += 1
    k = str(_n[0])
    M5 = Manifold(5, 'M5_' + k); M5.chart('t x u v y')
    M4 = Manifold(4, 'M4_' + k); M4.chart('t x u v')
    g4col = g4col or [0, 0, 0, 0]
    g5 = M5.metric('g5_' + k, signature=-3)
    q = M4.metric('q_' + k, signature=-2)
    for i in range(4):
        for j in range(i, 4):
            g5[i, j] = qm[i, j]
            q[i, j] = qm[i, j]
        g5[i, 4] = g4col[i]
    g5[4, 4] = -SR(N)**2                    # eps N^2 with eps = -1 (spacelike)
    R5 = g5.ricci(); R5s = g5.ricci_scalar().expr()
    R4 = q.ricci(); R4s = q.ricci_scalar().expr()
    G5 = matrix(SR, 4, 4, lambda i, j: R5[i, j].expr() - qm[i, j] * R5s / 2)
    G4 = matrix(SR, 4, 4, lambda i, j: R4[i, j].expr() - qm[i, j] * R4s / 2)
    out = dict(g5=g5, q=q, R5=R5, R5s=R5s, R4s=R4s, Delta=G5 - G4)
    if any(c != 0 for c in g4col):
        return out
    # zero shift: n = (1/N) d_y, n.n = eps, N^2 = eps g_44
    N = SR(N)
    qi = qm.inverse()
    K = qm.apply_map(lambda e: diff(e, y)) / (2 * N)
    Kup = qi * K
    trK = Kup.trace()
    KK = K * qi * K
    KdK = (Kup * Kup).trace()
    Riem = g5.riemann()
    # E_ab = R(n, e_a, n, e_b) = g_44 R^4_{a4b} / N^2 = eps R^4_{a4b}
    E = matrix(SR, 4, 4, lambda a, b: eps * Riem[4, a, 4, b].expr())
    Gam = q.christoffel_symbols()
    DDN = matrix(SR, 4, 4, lambda a, b: diff(N, crd4[a], crd4[b])
                 - sum(Gam[c, a, b].expr() * diff(N, crd4[c]) for c in range(4)))
    out.update(N=N, qi=qi, K=K, trK=trK, KK=KK, KdK=KdK, E=E, Etr=(qi * E).trace(),
               dyK=K.apply_map(lambda e: diff(e, y)) / N, DDN=DDN,
               boxN=(qi * DDN).trace())
    return out


def gauss_basis(G):
    q = G['q']
    qm = matrix(SR, 4, 4, lambda i, j: q[i, j].expr())
    return [G['trK'] * G['K'] - G['KK'], G['E'],
            qm * (G['trK']**2 - G['KdK']), qm * G['Etr']]


def normal_basis(G):
    return [G['dyK'], G['trK'] * G['K'], G['KK'], G['DDN'] / G['N']]


def fit(target, basis, pts):
    """Exact coefficients c with target = sum c_k basis_k, fitted at several
    rational points on a polynomial instance (then re-checked symbolically).
    Several points are needed: q_ab X and q_ab Y are parallel at any one."""
    rows, rhs = [], []
    for pt in pts:
        for a in range(4):
            for b in range(a, 4):
                rows.append([QQbar(SR(B_[a, b]).subs(pt)) for B_ in basis])
                rhs.append(QQbar(SR(target[a, b]).subs(pt)))
    A = matrix(QQbar, rows); r = vector(QQbar, rhs)
    c = A.solve_right(r)
    ok("fit is determined (basis independent over the sample points)",
       A.rank() == len(basis))
    ok("fit is consistent (residual zero on every component at every point)", (A * c - r).is_zero())
    return [QQ(ci) for ci in c]


print("=" * 72)
print("gr_projection_audit -- is the 5D -> 4D projection complete?")
print("=" * 72)

# ----------------------------------------------------------------------------
# Polynomial instance (exact rationals) for fixing coefficients
qpoly = matrix(SR, [[1 + t*y/3 + x*y**2/5, t*x*y/11 + y/13, 0, 0],
                    [t*x*y/11 + y/13, -(1 + x*y/4 + t**2*y/7), 0, 0],
                    [0, 0, -(1 + t*y/6 + x*y**2/9), 0],
                    [0, 0, 0, -(1 + t*y/6 + x*y**2/9)]])
Npoly = 1 + t*x*y/5 + x/7 + y**2/3
pt = {t: QQ(1)/3, x: QQ(1)/2, u: 0, v: 0, y: QQ(1)/5}
pts = [pt, {t: QQ(2)/7, x: QQ(3)/4, u: 0, v: 0, y: QQ(2)/9}, {t: QQ(5)/6, x: QQ(1)/8, u: 0, v: 0, y: QQ(3)/7}]

print("\nComputing the polynomial instance (general lapse, y-dependent q) ...")
P = geometry(qpoly, N=Npoly)
ok("residual G5_ab - G4_ab is non-zero in general (the question is real)",
   not mat_zero(P['Delta'].subs(pt)))

print("\n[A] Contracted Gauss equation with general lapse")
cG = fit(P['Delta'], gauss_basis(P), pts)
print("    G5_ab - G4_ab = (%s)[K K_ab - K_ac K^c_b] + (%s) E_ab"
      "\n                    + q_ab[(%s)(K^2 - K.K) + (%s) E]" % tuple(cG))
print("    (eps = %d)" % eps)

print("\n[B] Normal curvature E_ab with general lapse")
cN = fit(P['E'], normal_basis(P), pts)
print("    E_ab = (%s)(1/N) d_y K_ab + (%s) K K_ab + (%s) K_ac K^c_b"
      " + (%s)(1/N) D_a D_b N" % tuple(cN))
ok("the lapse Hessian enters E_ab with non-zero weight", cN[3] != 0)
ok("d_y K enters E_ab with non-zero weight", cN[0] != 0)

# Abstract family: both identities hold identically, not just at a point
print("\nAbstract family (functions of t, x, y; non-diagonal; general lapse) ...")
A_, B_, C_, f_, N_ = [function(n)(t, x, y) for n in ('A', 'B', 'C', 'f', 'N')]
qabs = matrix(SR, [[exp(2*A_), f_, 0, 0],
                   [f_, -exp(2*B_), 0, 0],
                   [0, 0, -exp(2*C_), 0],
                   [0, 0, 0, -exp(2*C_)]])
Gab = geometry(qabs, N=exp(N_))
ok("[A] holds identically on the abstract family: the residual is EXACTLY "
   "K-terms + normal-curvature terms",
   mat_zero(Gab['Delta'] - sum(c * Bm for c, Bm in zip(cG, gauss_basis(Gab)))))
ok("[B] holds identically on the abstract family",
   mat_zero(Gab['E'] - sum(c * Bm for c, Bm in zip(cN, normal_basis(Gab)))))
mut = list(cG); mut[1] = -mut[1]
ok("mutation: flipping the sign of the E_ab term fails",
   not mat_zero(Gab['Delta'] - sum(c * Bm for c, Bm in zip(mut, gauss_basis(Gab)))))

# Combined: residual in terms of K, d_y K and the lapse
c1, c2, c3, c4 = cG
print("\n    => with K = 0 at a point, the residual there is")
print("       (%s)(1/N)[d_y K_ab + ...] + Hessian terms -- E_ab, not K, decides."
      % (c2 * cN[0]))

# ----------------------------------------------------------------------------
print("\n[C] Sufficiency: K = 0 on every slice, zero shift, g_44 constant")
A0, B0, C0, f0 = [function(n)(t, x) for n in ('A0', 'B0', 'C0', 'f0')]
qcyl = qabs.subs({A_: A0, B_: B0, C_: C0, f_: f0})
Gc = geometry(qcyl)                              # g_44 = eps
ok("K = 0 on every slice (q independent of y)", mat_zero(Gc['K']))
ok("E_ab = 0 (no normal curvature: d_y K = 0, lapse constant)", mat_zero(Gc['E']))
ok("G5_ab = G4_ab identically for ARBITRARY q_ab(x): the projection is complete",
   mat_zero(Gc['Delta']))

# ----------------------------------------------------------------------------
print("\n[D] Necessity: drop one condition at a time")

# (D1) K = 0 on one slice y = 0 only
y0 = 0
qD1 = matrix(SR, [[1 + t*x/5 + y**2*(1 + x/3), y**2*t/7, 0, 0],
                  [y**2*t/7, -(1 + x**2/9) - y**2*t/2, 0, 0],
                  [0, 0, -(1 + t/4), 0],
                  [0, 0, 0, -(1 + t/4) - y**2*x/3]])
D1 = geometry(qD1)
p0 = {t: QQ(1)/3, x: QQ(1)/2, u: 0, v: 0, y: y0}
ok("(D1) K_ab = 0 on the slice y = 0", mat_zero(D1['K'].subs({y: y0})))
ok("(D1) but d_y K != 0 there, so E_ab != 0", not mat_zero(D1['E'].subs(p0)))
ok("(D1) and the residual G5_ab - G4_ab != 0 on that slice: K = 0 on ONE "
   "slice is not enough", not mat_zero(D1['Delta'].subs(p0)))
ok("(D1) residual there = (%s) E_ab + (%s) q_ab E, as [A] predicts" % (c2, c4),
   mat_zero((D1['Delta'] - c2 * D1['E']
             - c4 * qD1 * D1['Etr']).subs({y: y0})))

# (D2) radion: g_44 = eps N(x)^2, q independent of y
Nr = function('Nr')(t, x)
D2 = geometry(qcyl, N=Nr)
ok("(D2) K = 0 on every slice, yet", mat_zero(D2['K']))
hess = cN[3] * (c2 * D2['DDN'] + c4 * qcyl * D2['boxN']) / Nr
ok("(D2) residual = (%s)(D_a D_b N)/N + (%s) q_ab (box N)/N exactly"
   % (c2 * cN[3], c4 * cN[3]), mat_zero(D2['Delta'] - hess))
ok("(D2) ... which is non-zero for a generic N(x): a free g_44 spoils the "
   "reduction", not mat_zero(D2['Delta']))
print("    (it vanishes iff D_a D_b N = 0: tracing (DDN - q box N) gives box N = 0.)")
print("    DT holds g_44 fixed (kk_radion_g44: Cassini selects it), so N = 1.")

# (D3) shift g_a4 = A_a(x), q and g_44 independent of y
eta4 = diagonal_matrix(SR, [1, -1, -1, -1])
Ag = [function('Ag%d' % i)(t, x) for i in range(4)]
D3g = geometry(eta4, g4col=Ag)
Gam5 = D3g['g5'].christoffel_symbols()
ok("(D3) a generic shift g_a4 = A_a(x): Gamma^4_ab != 0, so free fall leaves "
   "the slice (u4 = 0 violated, K != 0)",
   not all(is_zero(Gam5[4, a, b].expr()) for a in range(4) for b in range(a, 4)))
# Killing shift on flat q: rotation in the x-u plane, dA != 0
Akill = [0, -u, x, 0]          # A_a = eta_ab xi^b for xi = (0, u, -x, 0)
D3k = geometry(eta4, g4col=Akill)
Gk = D3k['g5'].christoffel_symbols()
ok("(D3) a Killing shift (rotation, dA != 0): Gamma^4_ab = 0, u4 = 0 holds",
   all(is_zero(Gk[4, a, b].expr()) for a in range(4) for b in range(a, 4)))
ok("(D3) ... but G5_ab - G4_ab != 0: g_a4 = 0 is an INDEPENDENT condition",
   not mat_zero(D3k['Delta']))
D3t = geometry(eta4, g4col=[QQ(1)/3, 0, 0, 0])
ok("(D3) a constant shift (dA = 0) leaves no residual: the field strength, "
   "not the gauge, is what survives", mat_zero(D3t['Delta']))

# ----------------------------------------------------------------------------
print("\n[E] The other components under the conditions (class [C])")
R5c = Gc['R5']
ok("G5_a4 = R5_a4 = 0: the Codazzi constraint T_a4 = 0",
   all(is_zero(R5c[a, 4].expr()) for a in range(4)))
ok("G5_44 = -(eps/2) R4: the T_44 constraint",
   is_zero(R5c[4, 4].expr() - eps * Gc['R5s'] / 2 + eps * Gc['R4s'] / 2))
ok("R5 = R4 on the class", is_zero(Gc['R5s'] - Gc['R4s']))

print("""
VERDICT
  1. For a general 5D metric the projected Einstein tensor differs from the
     surface's by K-terms AND the normal curvature E_ab ([A], identically on
     an abstract non-diagonal family with general lapse).  E_ab carries
     d_y K and the lapse Hessian ([B]).  So the review's point is right:
     K_ab = 0 by itself does not give the reduction.
  2. The reduction is complete for EVERY metric with
       (i)   K_ab = 0 on EVERY slice x4 = const (u4 = 0 for all free motion
             at every mass density, not on one surface),
       (ii)  g_a4 = 0,
       (iii) g_44 constant (or, more weakly, D_a D_b N = 0),
     with q_ab(x) arbitrary ([C]).  These are the conditions the recovery
     states: (i) from u4 = 0, (ii)-(iii) from the restriction that holds
     g_a4 and g_44 fixed.  G2 therefore holds for the full class DT admits
     under them, not only for an ansatz.
  3. Each condition is necessary ([D]): K = 0 on one slice leaves d_y K;
     a free g_44 leaves a radion Hessian; a Killing g_a4 with dA != 0
     satisfies u4 = 0 yet leaves a residual, so g_a4 = 0 is not implied by
     K = 0 and must come from the restriction.
  4. Wording: the manuscript's 'for all free motion' should say 'on every
     surface x4 = const' (condition (i) is about all slices), and the
     reduction should name (ii)-(iii) where Gauss is invoked.
  Bracket: the 5D field equations are taken Einstein-type (as in
  gr_on_4d_surface); the identities are checked on abstract families of
  functions of (t, x, y) (two coordinates kept symbolic for cost), and the
  general Gauss-Codazzi theorem is standard, not re-proven.
""")
print("=" * 72)
print("ALL ASSERTIONS PASSED")
print("=" * 72)
