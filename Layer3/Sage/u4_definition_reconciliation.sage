#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE (primary): LA-8370-MS (1980): x^4 = gamma/a0, u^i = dx^i/ds, and the
#   five-force density.  Williams (2001) p.106 and Book
#   Eq 1.36 / 10.6: u4 = dy4/dt = d_t y4 + div(y4 u) = 0, "the conservation
#   of the fifth dimension".
# Context: u4_geometric_motion.sage (motion along x4 leaves the 4D path a
#   geodesic), u4_rest_mass_nonconservation.sage (F5: the SIZE of u4),
#   gr_projection_audit.sage (G2), matter_source_tab.sage (G5).
# Precedents for the definition: the Lagrangian mass coordinate of stellar
#   structure (Kippenhahn, Weigert & Weiss) and the matter-space (label)
#   formulation of relativistic fluids (Andersson & Comer, Living Rev.
#   Relativ.), in which the conserved current is built from comoving labels.
# Supplementary Materials: step [F5], Sec. S2 (code index: Sec. S6).
# ---------------------------------------------------------------------------
#
# u4_definition_reconciliation.sage                    [A] (assertion-hardened)
#
# REVIEW POINT.  With x4 = gamma(x)/a0 read as the density FIELD, dx4/ds along
# a path is (1/a0)(dt/ds) D gamma/Dt, while the continuity residual is
# D gamma/Dt + gamma div v.  These are different quantities, and no
# substitution makes them equal.  This script therefore DEFINES the mass
# coordinate of a continuum element as its mass per unit REFERENCE volume,
#
#     x4 = gamma0 / a0,    gamma0 = gamma J     (J = det dx/dX, X = labels),
#
# a comoving (Lagrangian) mass coordinate, and checks:
#
#   [A] residual - D gamma/Dt = gamma div v identically: the field reading and
#       the continuity residual agree only for incompressible flow.
#   [B] the reviewer's test: a compressible 3D flow with shear (div v != 0)
#       that conserves mass.  Field reading: dx4/dt = D gamma/Dt / a0 != 0.
#       Comoving coordinate: dx4/dt = 0.  Continuity residual: 0.  In label
#       form, gamma = gamma0(X(t,x)) det(dX/dx) satisfies continuity
#       identically and x4 = gamma0(X)/a0 is advected with the flow.
#   [C] for ANY density field (mass conserved or not)
#         d/dt (gamma J) |_X = J * [d_t gamma + div(gamma v)],
#       so u4 = (J/a0)(dt/ds) * (continuity residual): an identity, not an
#       assertion.  Williams' conservative form (Book Eq 1.36) is J^-1 times
#       the rate of the comoving coordinate, at every instant.
#   [D] the two definitions coincide where they can be compared: J = 1 at
#       the reference instant, and for any element whose volume does not
#       change (point bodies, rigid bodies, static stars).
#   [E] compatibility with the 5D geometry.  On the recovered product metric
#       (g_a4 = 0, g_44 = -1, g_ab independent of x4) d_4 is a symmetry
#       direction, and motion along x4 leaves the 4D path a geodesic
#       (u4_geometric_motion.sage).  So where an element sits on the mass
#       axis is invisible to 4D physics, and the definition is a choice of
#       how matter is placed in the 5D space.  Control: the same relabelling
#       done as a COORDINATE change x4 -> x4 J(x) would bring in g_a4 != 0 and
#       g_44 != -1, which the GR projection excludes.  It is a relabelling of
#       the matter, not of the coordinates.
#
# RESULT: with x4 the comoving mass coordinate, u4 is the continuity residual
# by an identity, every mass-conserving process (compression, shocks,
# collapse) has u4 = 0, and the c-cap bounds only mass non-conservation.  The
# field reading differs by gamma div v.  Bracket: this is OUR definition,
# consistent with the source's conservative form (Eq 1.36) and equal to its
# x4 = gamma/a0 for single bodies; it is not the source's literal wording for
# a continuum.
#
# RUN: sage u4_definition_reconciliation.sage
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("u4_definition_reconciliation.sage -- the comoving mass coordinate")
print("=" * 74)

t, x, y, z = var('t x y z')
X, Y, Z = var('X Y Z')
a, b, c, a0 = var('a b c a0')
assume(a > 0, b > 0, a0 > 0)

def div(F):
    return diff(F[0], x) + diff(F[1], y) + diff(F[2], z)

def residual(g, v):
    return diff(g, t) + div([g * vi for vi in v])

def material(g, v):
    return diff(g, t) + v[0] * diff(g, x) + v[1] * diff(g, y) + v[2] * diff(g, z)

# ---------------------------------------------------------------------------
print("\n[A] General identity (3D, arbitrary density and velocity fields)")
G = function('G')(t, x, y, z)
V = [function('V%d' % k)(t, x, y, z) for k in range(3)]
check((residual(G, V) - material(G, V) - G * div(V)).expand() == 0,
      "residual - D gamma/Dt = gamma div v identically")
print("    so the two agree only when div v = 0 (incompressible flow)")

# explicit, invertible, compressive 3D Lagrangian map with shear
chi = [X * exp(a * t), Y * exp(-b * t) + c * t * X, Z * (1 + t**2)]
Jmat = matrix(SR, [[diff(ci, Q) for Q in (X, Y, Z)] for ci in chi])
J = Jmat.det().simplify_full()
Xs = x * exp(-a * t)
Ys = (y - c * t * Xs) * exp(b * t)
Zs = z / (1 + t**2)
inv = {X: Xs, Y: Ys, Z: Zs}
check(all((ci.subs(inv) - w).simplify_full() == 0
          for ci, w in zip(chi, (x, y, z))), "inverse map correct")
vC = [diff(ci, t).subs(inv).simplify_full() for ci in chi]
on_path = {x: chi[0], y: chi[1], z: chi[2]}

# ---------------------------------------------------------------------------
print("\n[B] Reviewer's test: compressible flow that conserves mass")
check(div(vC).simplify_full() != 0, "div v != 0 (compressible, with shear)")
rho0 = 1 + X**2 + Y * Z                       # mass per reference volume
dXdx = matrix(SR, [[diff(q, w) for w in (x, y, z)] for q in (Xs, Ys, Zs)])
gB = (rho0.subs(inv) * dXdx.det()).simplify_full()   # label form of gamma
check(residual(gB, vC).simplify_full() == 0,
      "label form gamma = gamma0(X) det(dX/dx): continuity holds identically")
field_rate = (material(gB, vC) / a0).simplify_full()
check(field_rate != 0, "field reading: dx4/dt = (D gamma/Dt)/a0 != 0")
como = (diff((gB.subs(on_path) * J).simplify_full(), t) / a0).simplify_full()
check(como == 0, "comoving coordinate gamma J/a0: dx4/dt = 0")
x4f = rho0.subs(inv) / a0
check(material(x4f, vC).simplify_full() == 0,
      "x4 = gamma0(X)/a0 is advected with the flow (d_t x4 + v.grad x4 = 0)")
pt = {t: 1, x: 1, y: 2, z: 3, a: 1, b: QQ(1)/2, c: 1, a0: 1}
print("    at a sample point: field rate %.3f, comoving rate 0, residual 0"
      % field_rate.subs(pt).n())

# ---------------------------------------------------------------------------
print("\n[C] Any density field, mass conserved or not")
gC = (1 + x**2 + y * z) * exp(-t) + t * x * y
lhs = diff(gC.subs(on_path) * J, t)
rhs = J * residual(gC, vC).subs(on_path)
check((lhs - rhs).simplify_full() == 0,
      "d/dt (gamma J)|_X = J * [d_t gamma + div(gamma v)] identically")
check(residual(gC, vC).simplify_full() != 0,
      "non-trivial: this field does not conserve mass, and the identity still holds")
check(((lhs / J) - residual(gC, vC).subs(on_path)).simplify_full() == 0,
      "Book Eq 1.36 = J^-1 x (rate of the comoving coordinate), at every t")

# ---------------------------------------------------------------------------
print("\n[D] Where the two definitions coincide")
check(J.subs(t=0) == 1, "reference instant: J = 1, so x4 = gamma/a0 there")
rigid = [X * cos(t) - Y * sin(t), X * sin(t) + Y * cos(t), Z]   # rotation
Jr = matrix(SR, [[diff(ci, Q) for Q in (X, Y, Z)] for ci in rigid]).det()
check(Jr.simplify_full() == 1,
      "rigid motion: J = 1 for all t, so gamma J = gamma (point/rigid/static bodies)")

# ---------------------------------------------------------------------------
print("\n[E] Compatibility with the 5D geometry")
r_, w_, yy = var('r_ w_ yy')
f, h = function('f')(t, r_), function('h')(t, r_)
g5 = diagonal_matrix(SR, [f, -h, -1])           # (t, r, x4) block, product form
check(all(diff(g5[i, j], w_) == 0 for i in range(3) for j in range(3)),
      "product metric: no component depends on x4, so d_4 is a Killing direction")
Jf = function('Jf')(t, r_)                      # a relabelling factor J(x)
xmap = [t, r_, yy / Jf]                         # x4 = y/J(x): y = x4 J
Jac = matrix(SR, [[diff(q, s) for s in (t, r_, yy)] for q in xmap])
gnew = (Jac.transpose() * g5 * Jac).simplify_full()
check(gnew[0, 2] != 0 and gnew[1, 2] != 0,
      "control: x4 -> x4 J(x) as a coordinate change gives g_a4 != 0")
check((gnew[2, 2] + 1).simplify_full() != 0,
      "control: and g_44 != -1, both excluded by the GR projection (G2, G3)")
print("    so the comoving coordinate relabels where matter sits on the mass")
print("    axis, which 4D physics cannot see (u4_geometric_motion.sage); it is")
print("    not a change of the 5D coordinates.")

# ---------------------------------------------------------------------------
print("""
VERDICT
  The field reading x4 = gamma(x)/a0 and the continuity residual differ by
  gamma div v.  With x4 DEFINED as the comoving mass coordinate gamma J/a0
  (mass per reference volume), u4 = (J/a0)(dt/ds) x (continuity residual)
  identically, so every mass-conserving process has u4 = 0 and the c-cap
  bounds only mass non-conservation.  The definition equals x4 = gamma/a0
  for single bodies and at the reference instant, and it is compatible with
  the 5D geometry because the mass axis is a symmetry direction of the
  recovered metric.  Bracket: our definition, consistent with Eq 1.36, not
  the source's literal wording for a continuum.
""")
print("=" * 74)
print("ALL ASSERTIONS PASSED")
print("=" * 74)
