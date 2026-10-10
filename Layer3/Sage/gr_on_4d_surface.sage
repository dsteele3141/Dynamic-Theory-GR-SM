#!/usr/bin/env sage
# --------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, "The Dynamic Theory" (2011), Eqs (1.34)-
#   (1.39) and their Ch.10 repeat (10.1)-(10.9): "conservation of the fifth
#   dimension embeds a four-dimensional hyper-surface within the five
#   dimensional manifold and ... the equations giving the curvature of this
#   hyper-surface are Einstein's field equations."  u4 = dy4/dt = 0 is Eq 1.36
#   / 10.6.  Report LA-8370-MS eqs 163-164 give the induced metric (T2.11).
# Full derivation, context and provenance:
#   Supplementary_Material.pdf -- §S3, §6.1 (layer boundary), Appendix A.
# --------------------------------------------------------------------------

# ============================================================================
# gr_on_4d_surface.sage                           [A]  (assertion-hardened)
#
# "GR ON A 4D SURFACE" FROM u4 = 0 -- WHAT THE STEP ACTUALLY NEEDS.
#
# Williams' own text does not derive the 4D field equations: after building
# the surface energy-momentum tensor (Eq 1.37 / 10.7) he "suggests writing"
# G = C T and adds "when Equations (1.38) are taken to be Einstein's field
# equations".  The geometric content he needs is the Gauss-Codazzi reduction
# of the 5D curvature onto the hypersurface y4 = const.  This script supplies
# it and asks exactly what u4 = 0 contributes.
#
# Setting: Gaussian-normal 5D metric (the general local form near any
# hypersurface whose normal is the mass direction)
#       ds^2 = g_ab(x, y) dx^a dx^b + eps dy^2 ,   eps = -1  (y = y4 spacelike,
#       signature (+,-,-,-,-) as everywhere in the corpus),
# extrinsic curvature  K_ab = (1/2) d_y g_ab.
#
# WHAT THE SCRIPT ESTABLISHES
#   (A) Gauss-Codazzi, verified on an explicit non-diagonal family with
#       functions of (t, x, y):
#         R5_ab = R4_ab + (d_y K_ab - 2 K_am K^m_b + K K_ab)          (Gauss)
#         R5_44 = -(d_y K + K_ab K^ab)
#         R5_a4 = nabla_b K^b_a - d_a K                              (Codazzi)
#   (B) u4 = 0 under free fall.  The 5D geodesic equation gives
#         d^2 y / ds^2 = -Gamma^4_ab u^a u^b = eps K_ab u^a u^b ,
#       so a free particle STAYS on y = const (mass conserved, u4 = 0 kept)
#       iff K_ab u^a u^b = 0.  Requiring it for EVERY timelike u forces
#       K_ab = 0 (a quadratic form vanishing on an open cone vanishes).
#   (C) With K = 0 the reduction is exact: R5_ab = R4_ab, R5_44 = 0,
#       R5_a4 = 0, R5 = R4, hence G5_ab = G4_ab.  So 5D field equations
#       G5_AB = kappa T_AB restrict to the 4D Einstein equations
#       G4_ab = kappa T_ab (Williams' "surface energy-momentum tensor found
#       within the space tensor", Eq 1.37 / 10.7), plus two constraints:
#       T_a4 = 0, and T_44 = (eps/2) T^a_a -- for dust T_44 = -rho c^2/2,
#       the same mass-direction source weak_field_poisson.sage finds.
#   (D) With K != 0 the 4D equations acquire an induced source built from
#       K.  Explicit case: the 5D Ricci-flat metric
#         ds^2 = (y/L)^2 [dt^2 - e^{2t/L} dx_i dx_i] - dy^2
#       has K != 0 and induces on y = y0 a 4D G_ab proportional to the
#       metric with |Lambda| = 3/y0^2: curvature with no 4D matter.
#
# VERDICT (honest bracket).  Williams' sentence holds exactly under two
# conditions, of which u4 = 0 supplies one:
#   * u4 = 0 for all free motion  <=>  K_ab = 0 on the surface;
#   * the 5D field equations are Einstein-type (G5 = kappa T5).  DT does not
#     postulate these (weak_field_poisson.sage: the mass sector is a GAUGE
#     sector); the Ricci-flat vacuum case is the one DT uses (T2.14).
# Given both, GR on the 4D surface is DERIVED, not identified.  If mass is
# conserved only by a non-geodesic (fifth-force) constraint, K may be non-zero
# and the surface equations are Einstein's plus K-terms (D).  Type-I
# structural result about the reduction step; it does not test DT's physics.
#
# RUN: /home/dsteele/miniconda3/envs/sage/bin/sage gr_on_4d_surface.sage
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

print("=" * 72)
print("GR ON A 4D SURFACE from u4 = 0 -- Gauss-Codazzi reduction (5D -> 4D)")
print("=" * 72)

# ----------------------------------------------------------------------------
# Generic Gaussian-normal family: non-diagonal, functions of (t, x, y)
# ----------------------------------------------------------------------------
M5 = Manifold(5, 'M5')
X5 = M5.chart('t x u v y')
t, x, u, v, y = X5[:]
M4 = Manifold(4, 'M4')
X4 = M4.chart('t x u v')          # same symbolic names: y is a parameter here

A = function('A')(t, x, y)
B = function('B')(t, x, y)
C = function('C')(t, x, y)
f = function('f')(t, x, y)

g4m = matrix(SR, [[exp(2*A), f,          0,           0],
                  [f,        -exp(2*B),  0,           0],
                  [0,        0,          -exp(2*C),   0],
                  [0,        0,          0,   -exp(2*C)]])

g5 = M5.metric('g5', signature=-3)
g4 = M4.metric('g4', signature=-2)
for i in range(4):
    for j in range(i, 4):
        g5[i, j] = g4m[i, j]
        g4[i, j] = g4m[i, j]
g5[4, 4] = -1   # = eps, y4 spacelike

print("\nComputing 5D and 4D Ricci tensors (generic family) ...")
R5 = g5.ricci()
R4 = g4.ricci()
R5m = matrix(SR, 5, 5, lambda i, j: R5[i, j].expr())
R4m = matrix(SR, 4, 4, lambda i, j: R4[i, j].expr())

dy = lambda Mx: Mx.apply_map(lambda e: diff(e, y))
g4i = g4m.inverse()
K = dy(g4m) / 2                       # K_ab
Kmix = g4i * K                        # K^a_b
trK = Kmix.trace()
KK = K * g4i * K                      # K_am K^m_b
KabKab = (Kmix * Kmix).trace()

print("\n(A) Gauss-Codazzi")
gauss = dy(K) - 2 * KK + trK * K
ok("Gauss:  R5_ab = R4_ab + (d_y K_ab - 2 K_am K^m_b + K K_ab)",
   mat_zero(R5m[:4, :4] - R4m - gauss))
ok("mutation: the opposite sign of the K-terms fails",
   not mat_zero(R5m[:4, :4] - R4m + gauss))
ok("mutation: dropping -2 K_am K^m_b fails",
   not mat_zero(R5m[:4, :4] - R4m - (dy(K) + trK * K)))
ok("R5_44 = -(d_y K + K_ab K^ab)",
   is_zero(R5m[4, 4] + diff(trK, y) + KabKab))

# Codazzi: nabla_b K^b_a - d_a K with the 4D Levi-Civita connection
Kt = M4.tensor_field(1, 1, name='K')
for i in range(4):
    for j in range(4):
        Kt[i, j] = Kmix[i, j]
DK = g4.connection()(Kt)              # DK[b, a, c] = nabla_c K^b_a
coords4 = [t, x, u, v]
cod = [sum(DK[b, a, b].expr() for b in range(4)) - diff(trK, coords4[a])
       for a in range(4)]
sgn = None
for s in (1, -1):
    if all(is_zero(R5m[a, 4] - s * cod[a]) for a in range(4)):
        sgn = s
ok("Codazzi term is non-trivial on this family (sign test is meaningful)",
   not all(is_zero(c) for c in cod))
ok("Codazzi: R5_a4 = %s(nabla_b K^b_a - d_a K)" % ("+" if sgn == 1 else "-"),
   sgn is not None)

# ----------------------------------------------------------------------------
print("\n(B) u4 = 0 under free fall")
Gam5 = g5.christoffel_symbols()
ok("Gamma^4_ab = -eps K_ab  (so d^2y/ds^2 = eps K_ab u^a u^b)",
   all(is_zero(Gam5[4, a, b].expr() + eps * K[a, b])
       for a in range(4) for b in range(4)))
ok("Gamma^4_a4 = Gamma^4_44 = 0 (no other source of y-acceleration)",
   all(is_zero(Gam5[4, a, 4].expr()) for a in range(5)))

# K_ab u^a u^b = 0 for all u in an open (timelike) cone => K = 0.
# Local orthonormal frame: u = (1, w1, w2, w3), |w| < 1, an open set.
ks = [[var('k%d%d' % (min(i, j), max(i, j))) for j in range(4)] for i in range(4)]
Kg = matrix(SR, ks)
w = [var('w%d' % i) for i in range(1, 4)]
uu = vector(SR, [1] + w)
quad = (uu * Kg * uu).expand()
coeffs = [quad.subs({wi: 0 for wi in w})]
for i in range(3):
    coeffs.append(diff(quad, w[i]).subs({wi: 0 for wi in w}))
    for j in range(i, 3):
        coeffs.append(diff(quad, w[i], w[j]))
sol = solve([c == 0 for c in coeffs], [ks[i][j] for i in range(4) for j in range(i, 4)])
ok("K_ab u^a u^b = 0 on an open cone of timelike u forces all 10 K_ab = 0",
   all(eq.rhs() == 0 for eq in sol[0]) and len(sol[0]) == 10)

# ----------------------------------------------------------------------------
print("\n(C) K = 0 (mass conserved under free fall): the reduction is exact")
# the same family with the y-dependence removed (cylinder condition)
A0f, B0f, C0f, f0f = [function(n)(t, x) for n in ('A0', 'B0', 'C0', 'f0')]
g0m = g4m.subs({A: A0f, B: B0f, C: C0f, f: f0f})
M5c = Manifold(5, 'M5c'); X5c = M5c.chart('t x u v y')
M4c = Manifold(4, 'M4c'); X4c = M4c.chart('t x u v')
g5c = M5c.metric('g5c', signature=-3); g4c = M4c.metric('g4c', signature=-2)
for i in range(4):
    for j in range(i, 4):
        g5c[i, j] = g0m[i, j]
        g4c[i, j] = g0m[i, j]
g5c[4, 4] = -1   # = eps, y4 spacelike
R5c = g5c.ricci(); R4c = g4c.ricci()
ok("K = 0: R5_ab = R4_ab",
   all(is_zero(R5c[i, j].expr() - R4c[i, j].expr())
       for i in range(4) for j in range(4)))
ok("K = 0: R5_a4 = 0 and R5_44 = 0",
   all(is_zero(R5c[a, 4].expr()) for a in range(5)))
R5s = g5c.ricci_scalar().expr(); R4s = g4c.ricci_scalar().expr()
ok("K = 0: R5 = R4, hence G5_ab = G4_ab (5D Einstein eqs restrict to 4D)",
   is_zero(R5s - R4s))
ok("K = 0: G5_44 = -(eps/2) R4, so T_44 = (eps/2) T^a_a (constraint, not a 4D eq)",
   is_zero((R5c[4, 4].expr() - eps * R5s / 2) + eps * R4s / 2))

# ----------------------------------------------------------------------------
print("\n(D) K != 0: an induced source.  Ricci-flat 5D with de Sitter slices")
L, y0 = var('L y0', domain='positive')


def canonical(H):
    Mc = Manifold(5, 'Mc'); Xc = Mc.chart('t x u v y')
    gc = Mc.metric('gc', signature=-3)
    gc[0, 0] = (y / L)**2
    for i in (1, 2, 3):
        gc[i, i] = -(y / L)**2 * exp(2 * H * t)
    gc[4, 4] = -1   # = eps, y4 spacelike
    return gc


Rc = canonical(1 / L).ricci()
ok("ds^2 = (y/L)^2[dt^2 - e^{2t/L} dx^2] - dy^2 is 5D Ricci-flat",
   all(is_zero(Rc[i, j].expr()) for i in range(5) for j in range(i, 5)))
Rbad = canonical(2 / L).ricci()
ok("mutation: Hubble rate 2/L instead of 1/L is not Ricci-flat",
   not all(is_zero(Rbad[i, j].expr()) for i in range(5) for j in range(i, 5)))
ok("its extrinsic curvature is non-zero: K_tt = y/L^2",
   is_zero(diff((y / L)**2, y) / 2 - y / L**2))

Ms = Manifold(4, 'Ms'); Xs = Ms.chart('t x u v')
q = Ms.metric('q', signature=-2)
q[0, 0] = (y0 / L)**2
for i in (1, 2, 3):
    q[i, i] = -(y0 / L)**2 * exp(2 * t / L)
Gq = q.ricci() - q.ricci_scalar() * q / 2
lam_sgn = None
for s in (1, -1):
    if all(is_zero(Gq[i, j].expr() - s * 3 / y0**2 * q[i, j].expr())
           for i in range(4) for j in range(4)):
        lam_sgn = s
ok("induced metric on y = y0 has G_ab = %s(3/y0^2) q_ab: a cosmological term "
   "with no 4D matter" % ("+" if lam_sgn == 1 else "-"), lam_sgn is not None)

print("""
  NET: 'GR on a 4D surface' is DERIVED (not merely identified, as in Eq
  1.38/10.8) under two stated conditions:
    (i)  u4 = 0 for all free motion, which is equivalent to K_ab = 0;
    (ii) Einstein-type 5D field equations (Ricci-flat in vacuum, as T2.14).
  Then G5_ab = G4_ab exactly; T_a4 = 0 and T_44 = (eps/2) T^a_a are
  constraints (dust: T_44 = -rho c^2/2, as in weak_field_poisson.sage).
  If mass is conserved only under an imposed (non-geodesic) constraint,
  K may be non-zero and the surface equations are Einstein's plus K-terms.""")
print("=" * 72)
print("gr_on_4d_surface COMPLETE.")
print("=" * 72)
