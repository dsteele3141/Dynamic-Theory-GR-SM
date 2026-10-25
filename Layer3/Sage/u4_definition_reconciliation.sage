#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE (primary): LA-8370-MS (1980): x^4 = gamma/a0, u^i = dx^i/ds, and the
#   five-force density.  Williams (2001) p.106 and Book
#   Eq 1.36 / 10.6: u4 = dy4/dt = d_t y4 + div(y4 u) = 0, "the conservation
#   of the fifth dimension".
# Context: u4_rest_mass_nonconservation.sage (F5: the SIZE of u4),
#   matter_source_tab.sage (G5), weak_chirality_u4.sage (S1).
# Supplementary Materials: step [F5], Sec. S2 (code index: Sec. S6).
# ---------------------------------------------------------------------------
#
# u4_definition_reconciliation.sage                    [A] (assertion-hardened)
#
# REVIEW POINT.  The text defines u4 = dx4/ds with x4 = gamma/a0, and
# also calls u4 the residual of the continuity equation.  Read naively these
# differ: the rate of change of the density FIELD along a trajectory is the
# material derivative D gamma/Dt, which is non-zero under compression, while
# the continuity residual d_t gamma + div(gamma v) is zero whenever mass is
# conserved.  Checked here:
#
#   [A] residual - D gamma/Dt = gamma div v, identically (3D): the two agree
#       only for incompressible flow.
#   [B] a compression that conserves mass: residual = 0, D gamma/Dt != 0.
#       On the field reading, compression would violate u4 = 0.
#   [C] the reconciliation.  Follow a comoving element of a general 3D flow
#       (explicit Lagrangian map, Jacobian J) and take as its mass coordinate
#       the mass it carries per unit reference volume, gamma*J/a0.  Then
#         d/dt (gamma J) |_X  =  J * [d_t gamma + div(gamma v)]
#       identically.  Choosing the reference volume to be the current one at
#       the instant considered (J = 1), x4 = gamma/a0 there and
#         dx4/dt = (1/a0) (d_t gamma + div(gamma v)),
#       which is Williams' own conservative form (Book Eq 1.36).
#   [D] mutation: the field reading (D gamma/Dt) does not reproduce [C].
#
# RESULT: the two statements are consistent once dx4/ds is read as the rate
# of change along a COMOVING ELEMENT (Reynolds transport), not as the rate of
# change of the field gamma(x) along a path.  For a single body the two
# readings coincide (its mass per its own volume).  Bracket: this is a stated
# DEFINITION, supported by the source's conservative form, not a derivation.
#
# RUN: sage u4_definition_reconciliation.sage
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("u4_definition_reconciliation.sage -- dx4/ds vs the continuity residual")
print("=" * 74)

t, x, y, z = var('t x y z')
X, Y, Z = var('X Y Z')
a, b, c = var('a b c')
assume(a > 0, b > 0)

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

# ---------------------------------------------------------------------------
print("\n[B] Mass-conserving compression")
g0 = var('g0'); assume(g0 > 0)
L = exp(-a * t)                         # uniform 1D compression, x -> x L(t)
gB = g0 / L
vB = [x * diff(L, t) / L, 0, 0]
check(residual(gB, vB).simplify_full() == 0,
      "continuity residual = 0 (mass conserved)")
check(material(gB, vB).simplify_full() != 0,
      "but D gamma/Dt = a gamma != 0: the density of each element rises")
print("    On the field reading, compression would count as u4 != 0;")
print("    the source says it must not.")

# ---------------------------------------------------------------------------
print("\n[C] Comoving element: mass per unit reference volume")
# explicit, invertible, compressive 3D Lagrangian map with shear
chi = [X * exp(a * t), Y * exp(-b * t) + c * t * X, Z * (1 + t**2)]
Jmat = matrix(SR, [[diff(ci, Q) for Q in (X, Y, Z)] for ci in chi])
J = Jmat.det().simplify_full()
check(J != 1, "the map compresses/expands (J != 1), with shear")
# invert: X, Y, Z as functions of (t, x, y, z)
Xs = x * exp(-a * t)
Ys = (y - c * t * Xs) * exp(b * t)
Zs = z / (1 + t**2)
inv = {X: Xs, Y: Ys, Z: Zs}
check(all((ci.subs(inv) - w).simplify_full() == 0
          for ci, w in zip(chi, (x, y, z))), "inverse map correct")
# Eulerian velocity v(t, x) = d_t chi evaluated at X(t, x)
vC = [diff(ci, t).subs(inv).simplify_full() for ci in chi]
# a non-trivial density field (any function would do)
gC = (1 + x**2 + y * z) * exp(-t) + t * x * y
on_path = {x: chi[0], y: chi[1], z: chi[2]}
lhs = diff(gC.subs(on_path) * J, t)
rhs = J * residual(gC, vC).subs(on_path)
check((lhs - rhs).simplify_full() == 0,
      "d/dt (gamma J)|_X = J * [d_t gamma + div(gamma v)] identically")
at0 = {t: 0}
check(J.subs(at0) == 1, "reference = current volume at t = 0 (J = 1)")
check((lhs.subs(at0) - residual(gC, vC).subs(on_path).subs(at0)).simplify_full() == 0,
      "there x4 = gamma/a0 and a0 dx4/dt = continuity residual (Book Eq 1.36)")

# ---------------------------------------------------------------------------
print("\n[D] Mutation: the field reading")
field_rate = material(gC, vC).subs(on_path)
check((lhs.subs(at0) - field_rate.subs(at0)).simplify_full() != 0,
      "D gamma/Dt does NOT reproduce the comoving rate (differs by gamma div v)")

# ---------------------------------------------------------------------------
print("""
VERDICT
  u4 = dx4/ds is the continuity residual, and so vanishes for every
  mass-conserving process however its density changes, when dx4/ds is the
  rate of change along a comoving element: the mass it carries per unit
  (reference) volume, with the reference taken as the current volume.  That
  is the conservative form the source writes (Book Eq 1.36).  Read as the
  rate of change of the density field along a path (D gamma/Dt) it is not,
  and compression would then violate u4 = 0.  For a single body the readings
  coincide.  Bracket: a stated definition, not a derivation.
""")
print("=" * 74)
print("ALL ASSERTIONS PASSED")
print("=" * 74)
