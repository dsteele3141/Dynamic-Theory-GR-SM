#!/usr/bin/env sage
# --------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, LA-8370-MS (LASL, 1980), Eq (39) (4D) and
#   Eq (67) (5D):  delta INT [R + (1/4) alpha F^ F^ - lambda(1/2 - k phi^.phi^)]
#   sqrt(-g) = 0, k = 6 / 12, "Weyl's procedure in his unified field theory".
# COMPARISON (literature): D. M. Ghilencea, "Spontaneous breaking of Weyl
#   quadratic gravity to Einstein action and Higgs potential", JHEP 03 (2019)
#   049, arXiv:1812.08613 -- Eqs (5), (8)-(10), (25); and "Stueckelberg
#   breaking of Weyl conformal geometry with applications to gravity",
#   Phys. Rev. D 101 (2020) 045010, arXiv:1904.06596.  Precedent: L. Smolin,
#   Nucl. Phys. B 160 (1979) 253 (Planck-mass Weyl field coupled to GR).
# Full derivation, context and provenance:
#   Supplementary_Material.pdf -- §S1 (T1.9), Appendix A.
# Companions: weyl_lambda_reading.sage (natural gauge, |Lambda| = lambda/4),
#   weyl_vector_mass.sage (m^2 = 2k/a), einstein_from_dt_action.sage.
# --------------------------------------------------------------------------

# ============================================================================
# weyl_ghilencea_comparison.sage                [A]  (assertion-hardened)
#
# IS DT's ROUTE TO EINSTEIN + LAMBDA THE PUBLISHED WEYL-QUADRATIC-GRAVITY ONE?
#
# Ghilencea (2019) shows that Weyl's R~^2 gravity, after gauge fixing, IS the
# Einstein-Hilbert action with a positive cosmological constant plus a Proca
# action for the Weyl vector, whose mass is Planck-scale; below that mass the
# vector decouples and Weyl geometry becomes Riemannian.  The corpus reached
# the same three statements from Williams' Eqs 39/67.  This script checks,
# term by term and in one convention, that the two are the SAME calculation:
#
#   [1] Weyl scalar curvature, computed here from the Weyl connection in
#       signature (+,-,...):  R~ = R - 2(n-1) div(phi) - (n-1)(n-2) phi.phi,
#       i.e. k = 6 (4D), 12 (5D) -- the corpus' k.
#   [2] Dictionary phi = (q/2) omega maps it onto Ghilencea Eq (5),
#       R~ = R - 3q D.omega - (3/2) q^2 omega.omega.
#   [3] Linearisation.  Ghilencea Eq (10): (xi0/4!)[-2 phi0^2 R~ - phi0^4];
#       its phi0 equation gives phi0^2 = -R~ and returns (xi0/4!) R~^2.
#       Normalised by its R~ coefficient it is R~ + phi0^2/2: Williams'
#       "R - lambda/2" form with  lambda = phi0^2  (natural gauge R~ = -lambda
#       is Ghilencea's phi0^2 = -R~).  Only p = 2 gives Williams' "1/2".
#   [4] Cosmological constant.  M^2 = xi0 phi0^2/6 reproduces Ghilencea's
#       V = 3 M^4/(2 xi0), and Lambda = V/M^2 = phi0^2/4 = lambda/4 -- the
#       corpus' |Lambda| = lambda/4, reached independently.
#   [5] Weyl-vector mass.  The corpus' m^2 = 2k/a, with Ghilencea's field
#       normalisation (-1/4 F_omega^2), gives m^2 = (3/2) q^2 M^2 -- his
#       m_omega exactly; the 5D (k = 12) analogue is 3 q^2 M^2.
#   [6] Stueckelberg mode = DT's integrable sector.  For phi = d ln f
#       (Report Eq 65, F = 0) the Weyl gauge s = ln f removes phi and turns
#       (eta, d ln f) into the Riemannian (f^2 eta, 0):  R~ = f^2 R[f^2 eta],
#       checked for an arbitrary f.  So the vector the corpus finds vanishing
#       on shell is Ghilencea's "eaten" dilaton; both end in Riemannian GR.
#   [7] Mutations: phi = q omega breaks [2]; p = 3 gives 2/3, not 1/2.
#
# VERDICT (honest bracket).  The GR step DT takes (Eq 39/67 -> Einstein +
# Lambda + a heavy, decoupled Weyl vector) is the gauge-fixed Weyl-quadratic-
# gravity result published by Ghilencea, with the same Lambda = lambda/4 and
# the same mass formula.  So DT's GR recovery inherits a peer-reviewed
# result -- and its novelty here is correspondingly narrow: Weyl's procedure
# is Weyl's, the Planck-mass Weyl field is Smolin's (1979), the Stueckelberg
# reading is Ghilencea's.  What is DT's own: Eq 67 carries the procedure to
# 5D (k = 12, the 4D "1/2" kept), and u4 = 0 then gives the 4D surface.
# Differences, NOT checked here: Ghilencea's F^2 comes with the Weyl-tensor
# term C~^2 = C^2 + (3/2) q^2 F^2 (his Eq 9), so his action also carries the
# Riemannian C^2; Williams has F^2 alone.  Ghilencea's matter (a Higgs) has
# Weyl weight; the corpus assumes matter without it.  In both, the SIZE of
# Lambda is set by a free constant (xi0 there, lambda here); Ghilencea's
# xi0 > 0 makes Lambda > 0, a sign not re-derived in Williams' conventions.
#
# RUN: /home/dsteele/miniconda3/envs/sage/bin/sage weyl_ghilencea_comparison.sage
# ============================================================================

from sage.all import *


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


def z(e):
    return bool(SR(e).simplify_full() == 0)


def weyl_scalar(n, g, phi, X):
    """Scalar curvature of the Weyl connection
    G~^a_bc = G^a_bc + d^a_b phi_c + d^a_c phi_b - g_bc phi^a
    (g -> e^{2s} g, phi -> phi - ds leaves it invariant);
    R^a_bcd = d_c G^a_db - d_d G^a_cb + G^a_ce G^e_db - G^a_de G^e_cb,
    R_bd = R^a_bad.  With phi = 0 this is the Riemannian R."""
    gi = g.inverse()
    G = [[[sum(gi[a, e] * (diff(g[e, b], X[c]) + diff(g[e, c], X[b])
                            - diff(g[b, c], X[e])) for e in range(n)) / 2
           for c in range(n)] for b in range(n)] for a in range(n)]
    up = [sum(gi[a, e] * phi[e] for e in range(n)) for a in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                G[a][b][c] += ((a == b) * phi[c] + (a == c) * phi[b]
                               - g[b, c] * up[a])

    def Riem(a, b, c, d):
        return (diff(G[a][d][b], X[c]) - diff(G[a][c][b], X[d])
                + sum(G[a][c][e] * G[e][d][b] - G[a][d][e] * G[e][c][b]
                      for e in range(n)))
    Ric = [[sum(Riem(a, b, a, d) for a in range(n)) for d in range(n)]
           for b in range(n)]
    return sum(gi[b, d] * Ric[b][d] for b in range(n) for d in range(n))


print("=" * 72)
print("weyl_ghilencea_comparison.sage -- Williams Eq 39/67 vs Ghilencea 2019")
print("=" * 72)

# ---------------------------------------------------------------- [1] -----
print("\n[1] Weyl scalar curvature: R~ = R - 2(n-1) div phi - (n-1)(n-2) phi.phi")
for n in (4, 5):
    X = var(' '.join('x%d' % i for i in range(n)))
    eta = diagonal_matrix(SR, [1] + [-1] * (n - 1))
    cs = var(' '.join('c%d' % i for i in range(n)))
    A = matrix(SR, n, n, lambda i, j: var('a%d%d' % (i, j)))
    phi = [cs[m] + sum(A[m, v] * X[v] for v in range(n)) for m in range(n)]
    Rt = weyl_scalar(n, eta, phi, X)
    div = sum(eta[m, m] * A[m, m] for m in range(n))          # eta^{mm} d_m phi_m
    pp = sum(eta[m, m] * phi[m] ** 2 for m in range(n))
    ok("n=%d: R~ = -%d div(phi) - %d phi.phi on flat space (k = %d)"
       % (n, 2 * (n - 1), (n - 1) * (n - 2), (n - 1) * (n - 2)),
       z(Rt - (-2 * (n - 1) * div - (n - 1) * (n - 2) * pp)))

# ---------------------------------------------------------------- [2] -----
print("\n[2] Dictionary phi = (q/2) omega  ->  Ghilencea Eq (5)")
q, Rr, Dw, ww = var('q Rr Dw ww')        # R, D.omega, omega.omega (4D)
ours = Rr - 6 * (q / 2) * Dw - 6 * (q / 2) ** 2 * ww
ghil = Rr - 3 * q * Dw - QQ(3) / 2 * q ** 2 * ww
ok("R - 6 D.phi - 6 phi.phi  ==  R - 3q D.omega - (3/2) q^2 omega.omega",
   z(ours - ghil))

# ---------------------------------------------------------------- [3] -----
print("\n[3] Linearisation of R~^2: Ghilencea Eq (10) = Williams' 'R - lambda/2'")
xi0, Rt, P = var('xi0 Rt P')             # P = phi0^2
L10 = xi0 / 24 * (-2 * P * Rt - P ** 2)
sol = solve(diff(L10, P) == 0, P)
ok("phi0 equation gives phi0^2 = -R~", len(sol) == 1 and z(sol[0].rhs() + Rt))
ok("on shell it returns (xi0/4!) R~^2", z(L10.subs(P == -Rt) - xi0 / 24 * Rt ** 2))
coefR = L10.coefficient(Rt)
ok("normalised by its R~ coefficient: R~ + phi0^2/2",
   z(L10 / coefR - (Rt + P / 2)))
lam = var('lam')
ok("Williams/corpus natural gauge R~ + lambda/2 (R~ = -lambda) <=> lambda = phi0^2",
   z((Rt + lam / 2).subs(lam == P) - L10 / coefR))
for p, cp in ((2, QQ(1) / 2), (3, QQ(2) / 3)):
    s = var('s')
    f = s ** p
    lin = f + diff(f, s) * (Rt - s)             # linearised about R~ = s
    norm = (lin / diff(f, s)).subs(s == -lam).expand()
    ok("R~^%d linearised and normalised: R~ + %s lambda%s"
       % (p, cp, "" if p == 2 else "  ([7] mutation: not Williams' 1/2)"),
       z(norm - (Rt + cp * lam)) and ((cp == QQ(1) / 2) == (p == 2)))

# ---------------------------------------------------------------- [4] -----
print("\n[4] Cosmological constant: Lambda = phi0^2/4 = lambda/4")
M = var('M')
assume(xi0 > 0, P > 0)
M2 = -2 * coefR                             # coefficient of R~ is -M^2/2
ok("coefficient of R~ is -M^2/2 with M^2 = xi0 phi0^2/6", z(M2 - xi0 * P / 6))
V = xi0 / 24 * P ** 2                        # minus the constant term of L10
Psol = 6 * M ** 2 / xi0
ok("V = (xi0/4!) phi0^4 = 3 M^4/(2 xi0)  (Ghilencea Eq 25)",
   z(V.subs(P == Psol) - 3 * M ** 4 / (2 * xi0)))
Lam = V / M2                                 # -M^2/2 (R + 2 Lambda) form
ok("Lambda = V/M^2 = phi0^2/4", z(Lam - P / 4))
ok("with lambda = phi0^2: Lambda = lambda/4  (corpus |Lambda| = lambda/4)",
   z(Lam.subs(P == lam) - lam / 4))
ok("in terms of M: Lambda = 3 M^2/(2 xi0), > 0 for xi0 > 0 (Ghilencea's sign)",
   z(Lam.subs(P == Psol) - 3 * M ** 2 / (2 * xi0)))

# ---------------------------------------------------------------- [5] -----
print("\n[5] Weyl-vector mass: corpus m^2 = 2k/a vs Ghilencea m_omega^2")
# Ghilencea: -M^2/2 R~ - 1/4 F_omega^2.  Divide by -M^2/2:
#   R~ + (1/(2 M^2)) F_omega^2,  F_omega = (2/q) F_phi
#   => (1/4) a F_phi^2 with a = 4 * (1/(2 M^2)) * (2/q)^2
a = 4 * (1 / (2 * M ** 2)) * (2 / q) ** 2
ok("F^2 normalisation in corpus form: a = 8/(q^2 M^2)", z(a - 8 / (q ** 2 * M ** 2)))
m2_ours = {k: 2 * k / a for k in (6, 12)}
# direct: -M^2/2 * (-(3/2) q^2 omega^2) = +(1/2) m^2 omega^2
m2_ghil = 2 * (-M ** 2 / 2) * (-QQ(3) / 2 * q ** 2)
ok("Ghilencea's mass term gives m_omega^2 = (3/2) q^2 M^2",
   z(m2_ghil - QQ(3) / 2 * q ** 2 * M ** 2))
ok("corpus m^2 = 2k/a (k = 6) == Ghilencea m_omega^2", z(m2_ours[6] - m2_ghil))
ok("5D analogue (k = 12): m^2 = 3 q^2 M^2", z(m2_ours[12] - 3 * q ** 2 * M ** 2))
ok("[7] mutation: dictionary phi = q omega would give m^2 = 6 q^2 M^2 (no match)",
   not z(2 * 6 / (4 * (1 / (2 * M ** 2)) * (1 / q) ** 2) - m2_ghil))

# ---------------------------------------------------------------- [6] -----
print("\n[6] Integrable sector (phi = d ln f) = Ghilencea's Stueckelberg mode")
n = 4
X = var('t x y w')
sf = function('sf')(*X)                     # s = ln f, arbitrary
eta = diagonal_matrix(SR, [1, -1, -1, -1])
phi = [diff(sf, X[m]) for m in range(n)]
Fmn = [diff(phi[j], X[i]) - diff(phi[i], X[j]) for i in range(n) for j in range(n)]
ok("phi = d ln f is integrable: F = 0 (Report Eq 65)", all(z(e) for e in Fmn))
Rt_W = weyl_scalar(n, eta, phi, X)
gR = exp(2 * sf) * eta
R_riem = weyl_scalar(n, gR, [0] * n, X)
ok("(eta, d ln f) ~ (f^2 eta, 0):  R~ = f^2 R[f^2 eta]  for arbitrary f",
   z(Rt_W - exp(2 * sf) * R_riem))
ok("so after the gauge s = ln f the vector is zero and the geometry Riemannian",
   all(z(phi[m] - diff(sf, X[m])) for m in range(n)))

print("""
NET.  Williams' Eqs 39/67, read by Weyl's procedure, and Ghilencea's
spontaneously broken Weyl quadratic gravity are the same calculation:
same R~ (k = 6), same linearisation (lambda = phi0^2), same cosmological
constant (Lambda = lambda/4), same Weyl-vector mass (m^2 = (3/2) q^2 M_P^2),
and the same fate for the integrable vector (gauged away -> Riemannian).
DT's own additions: Eq 67 carries the procedure to 5D (k = 12) and u4 = 0
gives the 4D Einstein surface.  Not checked: Ghilencea's C~^2 term (absent
in Williams), Weyl-weighted matter, and the sign of Lambda in Williams'
conventions.""")
