#!/usr/bin/env sage
# --------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, LA-8370-MS (LASL, 1980), Eqs (58)-(63):
#   the entropy principle delta Int gamma sqrt(g_ij dq^i dq^j) = 0 and the
#   5-force density gamma f_k = d_k gamma - a0 u4 (g.u)_k (Eq 62; manuscript
#   Eq 6, eq:5force).  Book (Williams 2011) Eqs (6.2), (6.5)-(6.6): "if the
#   mass density is considered to be conserved such that u4 = 0 ... the volume
#   integral of Equation (6.5) becomes the force-mass-acceleration relation";
#   and Ch.10 / p.~14935: "within the Dynamic Theory gravitational and inertial
#   mass are equal".
# Comparison: the space-time-matter "fifth force" f^mu = -(1/m)(dm/ds) u^mu of
#   induced-matter Kaluza-Klein theory (Mashhoon, Wesson & Liu, Gen. Rel.
#   Grav. 30 (1998) 555; Wesson, "Space-Time-Matter", World Scientific 1999).
# Full derivation, context and provenance:
#   Supplementary_Material.pdf -- §S2 [C1], §6.1 (layer boundary), Appendix A.
# --------------------------------------------------------------------------

# ============================================================================
# test_body_motion.sage                            [A]  (assertion-hardened)
#
# HOW TEST BODIES MOVE: DOES DT's 5-FORCE LAW GIVE 4D GEODESICS?
#
# GR's motion postulate is that free bodies follow 4D geodesics of the
# surface metric, independent of their composition (weak equivalence
# principle).  DT's motion comes instead from its entropy principle, whose
# Lagrangian gamma ds is WEIGHTED by the mass density.  This script asks
# what that weighting does on the 4D surface of gr_on_4d_surface (K = 0,
# g_a4 = 0, g_44 = -1 held as in kk_radion_g44).
#
# WHAT THE SCRIPT ESTABLISHES
#   (A) Eq 62 DERIVED, not just quoted: Euler-Lagrange of L = gamma sqrt(Q)
#       on a fully general 5D metric (15 free functions of all 5 coordinates)
#       gives, at arc-length (Q = 1, dQ/ds = 0),
#         gamma (g_kl du^l/ds + Gamma_k,ij u^i u^j) = d_k gamma - u_k dgamma/ds,
#       which with dgamma/ds = a0 u4 is Eq 62.  (T2_1 only described it.)
#   (B) The 4D part of the extra term is EXACTLY the space-time-matter fifth
#       force, -(1/gamma)(dgamma/ds) u_mu, with m -> gamma.  It vanishes
#       identically under u4 = 0 (conservation of mass).
#   (C) On the K = 0 surface the 5D Christoffels reduce to the 4D ones, and
#       every u4-carrying term vanishes at u4 = 0.  With the d_mu gamma term
#       removed (D) the 4D equation is the GEODESIC equation with no body
#       parameter left in it: universal free fall (weak equivalence).
#   (D) The remaining 4D term d_mu gamma is the gradient of the body's own
#       density.  Integrated over a compact body it vanishes (spatial parts
#       are surface terms; the time part is dM/dq0 = 0 under conservation of
#       mass).  So the centre of mass moves on a 4D geodesic, to the usual
#       test-body order (tidal and spin terms not treated).
#   (E) The fifth (mass) equation is NOT satisfied by u4 = 0 on its own.
#       Along a path gamma = a0 q4, so d_4 gamma = a0 and the q4 equation
#       demands du4/ds = -1/q4 = -a0/gamma, ~1e4 per metre for ordinary
#       matter.  u4 = 0 is therefore held by a reaction along the mass axis,
#       exactly as in kk_radion_g44 (E): conservation of mass is a
#       RESTRICTION in the particle equations too, not a consequence.
#
# VERDICT (honest bracket).  Given conservation of mass imposed as a
# restriction, DT's own motion law reduces to GR's: compact bodies move on 4D
# geodesics of the surface metric, composition-independently, and the
# space-time-matter fifth force that 5D theories generically carry is
# switched off.  Without the restriction, the mass coordinate is driven at a
# huge rate, so the restriction is not optional.  The same restriction
# carried the field equations in kk_radion_g44; this is one assumption, used
# twice, not two.  Not treated: the strong equivalence principle (Nordtvedt
# effect, self-gravitating bodies), spin/tidal corrections, and the gauge-
# field view of gravity (Williams' universal beta = sqrt(4 pi eps0 G) is set
# by matching Newton, so EP there is by construction).
#
# RUN: /home/dsteele/miniconda3/envs/sage/bin/sage test_body_motion.sage
# ============================================================================

from sage.all import *


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


def is_zero(E):
    return bool(SR(E).expand().simplify_full() == 0)


print("=" * 74)
print("test_body_motion.sage -- DT's 5-force law vs 4D geodesic motion")
print("=" * 74)

n = 5
X = list(var('x0 x1 x2 x3 x4'))
U = list(var('u0 u1 u2 u3 u4'))
ACC = list(var('w0 w1 w2 w3 w4'))         # du^i/ds
a_0 = var('a_0', domain='positive')


def dds(F):
    """total s-derivative along the path, by the chain rule"""
    return (sum(diff(F, X[j]) * U[j] for j in range(n))
            + sum(diff(F, U[j]) * ACC[j] for j in range(n)))


# ---------------------------------------------------------------------------
print("\n[A] Euler-Lagrange of L = gamma sqrt(g_ij u^i u^j), general 5D metric")
gf = {}
for i in range(n):
    for j in range(i, n):
        gf[(i, j)] = function('g%d%d' % (i, j))(*X)
        gf[(j, i)] = gf[(i, j)]
g = matrix(SR, n, n, lambda i, j: gf[(i, j)])
gam = function('gam')(*X)
Q = sum(g[i, j] * U[i] * U[j] for i in range(n) for j in range(n))
L = gam * sqrt(Q)


def chris_low(k, i, j, gm):
    return (diff(gm[k, j], X[i]) + diff(gm[k, i], X[j]) - diff(gm[i, j], X[k])) / 2


dgam = sum(diff(gam, X[j]) * U[j] for j in range(n))
dQ = dds(Q)
allok = True
for k in range(n):
    EL = dds(diff(L, U[k])) - diff(L, X[k])
    uk = sum(g[k, l] * U[l] for l in range(n))
    geo = (sum(g[k, l] * ACC[l] for l in range(n))
           + sum(chris_low(k, i, j, g) * U[i] * U[j] for i in range(n) for j in range(n)))
    target = (Q**(-QQ(1) / 2) * (gam * geo + uk * dgam - Q * diff(gam, X[k]))
              - gam * uk * Q**(-QQ(3) / 2) * dQ / 2)
    allok = allok and is_zero((EL - target) * Q**(QQ(3) / 2))
ok("A1: EL_k = Q^-1/2 [gamma(g a + Gamma uu)_k + u_k dgamma/ds - Q d_k gamma]"
   " - (gamma u_k/2) Q^-3/2 dQ/ds, for all k, general metric", allok)
print("    At arc-length (Q = 1, dQ/ds = 0) EL = 0 reads")
print("      gamma f_k = d_k gamma - u_k dgamma/ds,  f_k = (g a + Gamma uu)_k")
print("    and dgamma/ds = a0 u4 (q4 = gamma/a0) gives Eq 62 / manuscript Eq 6.")

# ---------------------------------------------------------------------------
print("\n[B] The u4 term is the space-time-matter fifth force")
gs, ul = var('gs ul')
dgam_ds = a_0 * U[4]
fifth_DT = -a_0 * U[4] * ul / gs               # 4D part of -a0 u4 u_mu / gamma
fifth_STM = -(dgam_ds / gs) * ul               # -(1/m)(dm/ds) u_mu, m -> gamma
ok("B1: -a0 u4 u_mu / gamma = -(1/gamma)(dgamma/ds) u_mu (STM form, m -> gamma)",
   is_zero(fifth_DT - fifth_STM))
ok("B2: it vanishes identically at u4 = 0 (conservation of mass)",
   is_zero(fifth_DT.subs(u4=0)))

# ---------------------------------------------------------------------------
print("\n[C] The K = 0 surface: 5D Christoffels reduce to 4D ones")
X4 = X[:4]
hf = [function('h%d' % a)(*X4) for a in range(4)]
g5 = diagonal_matrix(SR, hf + [-1])
g4 = diagonal_matrix(SR, hf)


def chris_up(gm, m):
    gi = gm.inverse()
    return [[[sum(gi[a, d] * chris_low(d, b, c, gm) for d in range(m))
              for c in range(m)] for b in range(m)] for a in range(m)]


G5 = chris_up(g5, 5)
G4 = chris_up(g4, 4)
ok("C1: Gamma5^mu_nu,lam = Gamma4^mu_nu,lam (mu,nu,lam = 0..3)",
   all(is_zero(G5[a][b][c] - G4[a][b][c]) for a in range(4) for b in range(4) for c in range(4)))
ok("C2: Gamma5^mu_nu4 = Gamma5^mu_44 = Gamma5^4_ij = 0 (no mass-axis coupling)",
   all(is_zero(G5[a][b][4]) for a in range(5) for b in range(5))
   and all(is_zero(G5[4][b][c]) for b in range(5) for c in range(5)))
# 4D equation at u4 = 0 and without the d_mu gamma term:
f_mu = [ACC[a] + sum(G5[a][b][c] * U[b] * U[c] for b in range(5) for c in range(5))
        for a in range(4)]
geo4 = [ACC[a] + sum(G4[a][b][c] * U[b] * U[c] for b in range(4) for c in range(4))
        for a in range(4)]
ok("C3: at u4 = 0 the 5D acceleration f^mu IS the 4D geodesic operator",
   all(is_zero(f_mu[a].subs(u4=0) - geo4[a]) for a in range(4)))
ok("C4: the 4D equation f^mu = 0 contains no gamma, a0 or composition parameter"
   " (universal free fall)",
   all(not geo4[a].has(gam) and not geo4[a].has(a_0) for a in range(4)))

# ---------------------------------------------------------------------------
print("\n[D] The d_mu gamma term integrates to zero over a compact body")
xx, yy, zz, tt, Rb, rho0 = var('xx yy zz tt Rb rho0', domain='positive')
xx, yy, zz = var('xx yy zz', domain='real')
blob = rho0 * exp(-(xx**2 + yy**2 + zz**2) / Rb**2)
Fx = integral(integral(integral(diff(blob, xx), xx, -oo, oo), yy, -oo, oo), zz, -oo, oo)
ok("D1: Int d_x gamma d^3x = 0 for a compact (Gaussian) body", is_zero(Fx))
M0 = var('M0', domain='positive')
Rr, Rdot = var('Rr Rdot', domain='positive')  # body breathes: R(t), dR/dt = Rdot
blob_t = M0 / (pi**(QQ(3) / 2) * Rr**3) * exp(-(xx**2 + yy**2 + zz**2) / Rr**2)
dblob_dt = diff(blob_t, Rr) * Rdot            # chain rule d/dt = Rdot d/dR
Ft = integral(integral(integral(dblob_dt, xx, -oo, oo), yy, -oo, oo), zz, -oo, oo)
ok("D2: Int d_0 gamma d^3x = dM/dq0 = 0 for a body whose density changes"
   " (breathing) while its mass is conserved", is_zero(Ft))
print("    => Int gamma f_mu d^3x = 0: the centre of mass follows the 4D geodesic.")

# ---------------------------------------------------------------------------
print("\n[E] The mass equation: u4 = 0 needs a reaction")
qq4 = var('qq4', domain='positive')
gam_path = a_0 * qq4                          # gamma = a0 q4 on the path
d4gam = diff(gam_path, qq4)
u4lower = -U[4]                               # g_44 = -1
rhs4 = d4gam - u4lower * a_0 * U[4]           # d_4 gamma - u_4 dgamma/ds
f4_required = (rhs4 / gam_path).subs(u4=0)   # f_4 = -du4/ds on the K=0 surface
du4ds = (-f4_required).simplify_full()
print("    q4 equation at u4 = 0: du4/ds = %s" % du4ds)
ok("E1: du4/ds = -1/q4 = -a0/gamma != 0: u4 = 0 is NOT preserved freely",
   is_zero(du4ds + 1 / qq4))
a0num, gAl = 3.65e7, 2.70e3                   # kg m^-4 (Layer-3 script 18); Al density
rate = a0num / gAl
print("    aluminium: |du4/ds| = a0/gamma = %.2e per metre of proper path" % rate)
ok("E2: the free drive is enormous (> 1e3 per metre): the restriction is not optional",
   rate > 1e3)
print("    => conservation of mass is held by a reaction gamma f_4 = a0 along the")
print("       mass axis -- the particle analogue of the T_44 reaction in kk_radion_g44.")

# ---------------------------------------------------------------------------
print("\n[VERDICT]")
print("  * Eq 62 derived from delta Int gamma ds on a general 5D metric.")
print("  * Its u4 term is the STM fifth force; u4 = 0 switches it off.")
print("  * On the K = 0 surface, compact bodies follow 4D geodesics with no")
print("    composition parameter: GR's motion law and weak equivalence.")
print("  * u4 = 0 is not dynamically preserved; it is a restriction with a")
print("    reaction along the mass axis -- the same one the field equations need.")
print("  * Not treated: strong EP / Nordtvedt, spin and tidal terms.")
print("\n" + "=" * 74)
print("test_body_motion COMPLETE.")
print("=" * 74)
