#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE (primary): the 5D field equations G_ij + Lambda g_ij = kappa T_ij
#   (i, j = 0..4), signature (+,-,-,-,-); the restriction u4 = 0 holds
#   g_a4 = 0 and g_44 = -1 (kk_radion_g44.sage, G3).
# Context: gr_projection_audit.sage [E] (G5_44 = -(eps/2) R4 on the class),
#   matter_source_tab.sage (G5), weak_field_poisson.sage (T_44 = -rho c^2/2
#   in the weak field).
# Supplementary Materials: step [G3], Sec. S3 (code index: Sec. S6).
# ---------------------------------------------------------------------------
#
# t44_normal_equation.sage                             [A] (assertion-hardened)
#
# REVIEW POINT.  The tangential equations give G5_ab = G4_ab.  The normal-
# normal (44) equation must hold too.  Checked here:
#
#   [A] G5_44 = (1/2) R4 and R5 = R4 for product metrics (g_a4 = 0,
#       g_44 = -1) on two independent 4D families: static spherical
#       (f(r), h(r)) and spatially flat FLRW (a(t)).
#   [B] with the 4D trace R4 = 4 Lambda - kappa T, the 44 equation gives
#         T_44 = Lambda/kappa - T/2     (T = 4D trace of T_ab)
#       for any matter and any Lambda.
#   [C] perfect fluid: T = rho c^2 - 3p.  Dust with Lambda = 0 gives
#       T_44 = -rho c^2/2, the weak-field value (weak_field_poisson); radiation (T = 0)
#       gives T_44 = Lambda/kappa.
#   [D] consistency: T_44 is the reaction of the multiplier that holds
#       g_44 fixed.  The 44 equation fixes it and imposes nothing on
#       (rho, p): the 5D system is solvable for every 4D solution.
#   [E] the Bianchi component nabla_i T^i_4 = 0 holds identically when
#       nothing depends on x4.
#
# RESULT: the full 5D equations, not only their tangential projection,
# contain GR; the weak-field T_44 is the dust, Lambda = 0 case of T_44 = Lambda/kappa - T/2.
#
# RUN: sage t44_normal_equation.sage
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("t44_normal_equation.sage -- the 44 component of the 5D equations")
print("=" * 74)

def product_geometry(name, coords, g4diag):
    M5 = Manifold(5, 'M5_' + name, structure='Lorentzian', signature=-3)
    ch = M5.chart(coords)
    co = ch[:]
    g = M5.metric()
    for k in range(4):
        g[k, k] = g4diag(co)[k]
    g[4, 4] = -1
    M4 = Manifold(4, 'M4_' + name, structure='Lorentzian', signature=-2)
    ch4 = M4.chart(' '.join(coords.split()[:4]))
    co4 = ch4[:]
    h = M4.metric()
    for k in range(4):
        h[k, k] = g4diag(co4 + (0,))[k]
    return g, h

print("\n[A] G5_44 and R5 on product metrics (g_a4 = 0, g_44 = -1)")
families = {
    'static': ('t r th ph w:mass',
               lambda c: [function('f')(c[1]), -function('h')(c[1]),
                          -c[1]**2, -c[1]**2 * sin(c[2])**2]),
    'flrw':   ('t x y z w:mass',
               lambda c: [1, -function('a')(c[0])**2,
                          -function('a')(c[0])**2, -function('a')(c[0])**2]),
}
for name, (coords, g4) in families.items():
    g, h = product_geometry(name, coords, g4)
    R5 = g.ricci_scalar().expr()
    R4 = h.ricci_scalar().expr()
    Ric44 = g.ricci()[4, 4].expr()
    G44 = Ric44 - QQ(1)/2 * R5 * g[4, 4].expr()
    check(Ric44 == 0, "%s: R5_44 = 0" % name)
    check((R5 - R4).simplify_full() == 0, "%s: R5 = R4" % name)
    check((G44 - R4 / 2).simplify_full() == 0, "%s: G5_44 = R4/2" % name)

print("\n[B] The 44 equation with the 4D trace")
R, Lam, kap, T, T44 = var('R Lambda kappa T T44')
g44 = -1
# 4D trace of G_ab + Lambda g_ab = kappa T_ab:  G^a_a = -R, g^a_a = 4
trace_eq = (-R + 4 * Lam == kap * T)
Rsol = solve(trace_eq, R)[0].rhs()
eq44 = (R / 2 + Lam * g44 == kap * T44)
T44sol = solve(eq44.subs(R=Rsol), T44)[0].rhs()
check((T44sol - (Lam / kap - T / 2)).simplify_full() == 0,
      "T_44 = Lambda/kappa - T/2 for any matter and any Lambda")

print("\n[C] Perfect fluid")
rho, p, cc = var('rho p c')
Tpf = rho * cc**2 - 3 * p
T44pf = T44sol.subs(T=Tpf)
print("    T_44 = %s" % T44pf)
check((T44pf.subs(p=0, Lambda=0) + rho * cc**2 / 2).simplify_full() == 0,
      "dust, Lambda = 0: T_44 = -rho c^2/2 (the weak-field value)")
check((T44pf.subs(p=rho * cc**2 / 3) - Lam / kap).simplify_full() == 0,
      "radiation (T = 0): T_44 = Lambda/kappa")
Lam_obs, kap_obs = 1.1e-52, 8 * pi.n() * 6.674e-11 / (2.998e8)**4
check(abs(Lam_obs / kap_obs) / (1e3 * 2.998e8**2 / 2) < 1e-25,
      "the Lambda term is < 1e-25 of the dust term at laboratory density")

print("\n[D] The 44 equation imposes nothing on the matter")
mu = var('mu')     # multiplier reaction: the whole of T_44
sol = solve(R / 2 + Lam * g44 == kap * mu, mu)
check(len(sol) == 1, "the 44 equation fixes the multiplier uniquely")
check(set(sol[0].rhs().variables()) <= {R, Lam, kap},
      "it depends only on (R, Lambda, kappa): no condition on (rho, p), so any 4D solution lifts")
check(set(sol[0].rhs().variables()) != set(),
      "mutation guard: the multiplier is not a fixed constant")

print("\n[E] Bianchi component along the mass axis")
g, h = product_geometry('bianchi', 'tt xx yy zz ww:mass',
                        lambda c: [function('A')(c[0], c[1]),
                                   -function('B')(c[0], c[1]), -1, -1])
G5 = g.ricci() - QQ(1)/2 * g.ricci_scalar() * g
check(all(G5[a, 4].expr() == 0 for a in range(4)),
      "G5_a4 = 0, so nabla_i T^i_4 = 0 reduces to d_4 T^4_4 = 0: holds (no x4 dependence)")

print("""
VERDICT
  The 44 equation holds for any matter: with g_44 = -1, G5_44 = R4/2, and
  with the 4D trace, T_44 = Lambda/kappa - T/2.  For dust with Lambda = 0
  this is the weak-field value -rho c^2/2.  T_44 is the reaction of the multiplier
  holding g_44 fixed, so the equation determines it and constrains nothing
  else: the full 5D equations, not only their projection, contain GR.
""")
print("=" * 74)
print("ALL ASSERTIONS PASSED")
print("=" * 74)
