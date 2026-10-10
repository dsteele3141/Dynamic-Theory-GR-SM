#!/usr/bin/env sage
# --------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, "The Dynamic Theory" (2011), Eqs (1.33)-
#   (1.39) and their Ch.10 repeat: the induced surface metric g_ab = a_ab +
#   h_ab (Eq 1.33), the space energy-momentum tensor T^ij_sp = gamma u^i u^j +
#   field terms (Eq 1.34), u4 = dy4/dt = 0 as "the statement required by the
#   conservation of the fifth dimension" (Eq 1.36), and ONLY the surface field
#   equations G^ab = C T^ab (Eq 1.38).  Conservation of mass is called a
#   "restriction".  Report LA-8370-MS eqs 163-164 (T2.11).
# Observational bound: Bertotti, Iess & Tortora, Nature 425 (2003) 374
#   (Cassini): gamma - 1 = (2.1 +/- 2.3) x 10^-5.
# Full derivation, context and provenance:
#   Supplementary_Material.pdf -- §S3, §6.1 (layer boundary), Appendix A.
# --------------------------------------------------------------------------

# ============================================================================
# kk_radion_g44.sage                               [A]  (assertion-hardened)
#
# THE EXTRA 5D FIELDS g_44 (RADION) AND g_a4 (KK VECTOR) -- DO THEY SPOIL GR?
#
# Every 5D theory that recovers 4D gravity must dispose of the two extra
# metric fields.  A free g_44 is a Brans-Dicke scalar coupled to matter; it
# gives the PPN parameter gamma = 1/2, excluded by Cassini by ~10^4.  The
# corpus so far held g_44 = -1 fixed (weak_field_poisson, gr_on_4d_surface)
# without justification.  This script asks what in DT fixes it.
#
# Setting: linearised 5D metric about eta = diag(+,-,-,-,-), static source,
#   g_00 = 1 + 2 Phi,  g_ij = -(1 - 2 Psi) d_ij,  g_44 = -(1 + 2 Chi),
#   g_04 = a  (KK vector, static),
# PPN gamma = Psi / Phi.  Matter: Williams' Eq 1.34 dust with u4 = 0.
#
# WHAT THE SCRIPT ESTABLISHES
#   (A) Linearised 5D Einstein tensor with a free radion Chi:
#         G_00 = 2 lap Psi - lap Chi,   G_44 = lap Phi - 2 lap Psi,
#         G_ij = -(d_i d_j - d_ij lap)(Phi + Chi - Psi),   G_04 ~ lap a.
#   (B) u4 = 0 dust (T^44 = T^a4 = T^ij = 0) with g_44 DYNAMICAL:
#         Psi = Phi/2, Chi = -Phi/2  =>  gamma = 1/2, G_eff = (4/3) G.
#       Excluded by Cassini (|gamma - 1| = 0.5 vs 2.3e-5).
#   (C) u4 = 0 does NOT fix the radion: for ds^2 = g_ab(x) dx dx - N(x)^2 dy^2
#       the 5D geodesic keeps dy/ds = 0 for ANY N(x).  Conservation of mass
#       along worldlines is compatible with a free g_44.
#   (D) The lambda term cannot rescue it: the only scale it supplies is
#       Lambda = lambda/4, so any radion mass is O(sqrt Lambda), range ~1e26 m,
#       ~1e15 AU -- no screening in the solar system (order of magnitude).
#   (E) Conservation of mass as a CONSTRAINT (a "restriction", imposed by a
#       Lagrange multiplier on g_44 = -1 rather than varied): the surface
#       equations give Psi = Phi, gamma = 1 and Newton's coefficient exactly;
#       the 44 equation no longer constrains the metric and instead fixes the
#       multiplier, whose value is T_44 = -rho c^2/2 -- precisely the "mass-
#       direction source" that gr_on_4d_surface (C) and weak_field_poisson
#       (C3) found was needed.  That source IS the reaction of the restriction.
#       With K = 0 (gr_on_4d_surface) the surface equations are exactly 4D
#       GR, so beta = 1 as well.
#   (F) KK vector: T_04 = rho u_0 u_4 = 0 under u4 = 0, so g_04 obeys an
#       unsourced Laplace equation and vanishes for decaying boundary data.
#   (G) Correction to weak_field_poisson branch B: with g_44 frozen, imposing
#       G_44 = 0 (Psi = Phi/2) leaves G_ij != 0, i.e. it needs anisotropic
#       stress; it is not a consistent dust solution.  The consistent
#       "unsourced fifth direction" result is (B), gamma = 1/2.
#
# VERDICT (honest bracket).  The radion does not go away by itself, and u4 = 0
# along worldlines does not remove it.  DT recovers GR's PPN gamma = 1 iff
# conservation of mass is a restriction on the variation (Williams' word, and
# his practice: he writes only the surface equations 1.38), not a dynamical
# 5D field equation for g_44.  On that reading the extra fields carry no
# propagating scalar, the earlier T_44 = -rho c^2/2 is explained as the
# constraint reaction, and the 4D surface obeys GR exactly.  On the dynamical
# reading DT is excluded by Cassini.  So Cassini SELECTS the reading; it does
# not derive it.  Linear order for gamma; beta via the exact K = 0 reduction.
#
# RUN: /home/dsteele/miniconda3/envs/sage/bin/sage kk_radion_g44.sage
# ============================================================================

from sage.all import *


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


def is_zero(E):
    return bool(SR(E).simplify_full() == 0)


print("=" * 74)
print("kk_radion_g44.sage -- the extra 5D fields g_44 and g_a4 vs PPN gamma")
print("=" * 74)

t, x, y, z, w = var('t x y z w')
X = [t, x, y, z, w]
n = 5
e = var('e')
Phi = function('Phi')(x, y, z)
Psi = function('Psi')(x, y, z)
Chi = function('Chi')(x, y, z)
A = function('A')(x, y, z)


def lap(F):
    return diff(F, x, 2) + diff(F, y, 2) + diff(F, z, 2)


def christoffel(g):
    gi = g.inverse()
    return [[[sum(gi[a, d] * (diff(g[d, b], X[c]) + diff(g[d, c], X[b])
                              - diff(g[b, c], X[d])) for d in range(n)) / 2
              for c in range(n)] for b in range(n)] for a in range(n)]


def ricci(g):
    G = christoffel(g)

    def Ric(b, d):
        return sum(diff(G[a][d][b], X[a]) - diff(G[a][a][b], X[d])
                   + sum(G[a][a][f] * G[f][d][b] - G[a][d][f] * G[f][a][b]
                         for f in range(n)) for a in range(n))
    return matrix(SR, n, n, lambda i, j: Ric(i, j))


# ---------------------------------------------------------------------------
print("\n[A] Linearised 5D Einstein tensor with radion Chi and KK vector a")
eta = diagonal_matrix(SR, [1, -1, -1, -1, -1])
h = diagonal_matrix(SR, [2 * Phi, 2 * Psi, 2 * Psi, 2 * Psi, -2 * Chi])
h[0, 4] = A
h[4, 0] = A
g = eta + e * h
lin = lambda E: diff(E, e).subs(e=0).simplify_full()
Rl = ricci(g).apply_map(lin)
Rs = sum(eta[i, i] * Rl[i, i] for i in range(n))
Gl = (Rl - eta * Rs / 2).apply_map(lambda E: E.simplify_full())

S = Phi + Chi - Psi
ok("A1: G_00 = 2 lap Psi - lap Chi", is_zero(Gl[0, 0] - (2 * lap(Psi) - lap(Chi))))
ok("A2: G_44 = lap Phi - 2 lap Psi", is_zero(Gl[4, 4] - (lap(Phi) - 2 * lap(Psi))))
ok("A3: G_ij = -(d_i d_j - d_ij lap)(Phi + Chi - Psi), all i,j",
   all(is_zero(Gl[i, j] + diff(S, X[i], X[j]) - (lap(S) if i == j else 0))
       for i in range(1, 4) for j in range(1, 4)))
ok("A4: G_04 = (1/2) lap a  (KK vector decoupled at linear order)",
   is_zero(Gl[0, 4] - lap(A) / 2))
ok("A5: no mixing of a into the scalar sector (G_00, G_44, G_ij free of a)",
   all(not Gl[i, j].has(A) for (i, j) in [(0, 0), (4, 4), (1, 1), (1, 2)]))

# Work with Laplacians (fields decay at infinity, so G_ij = 0 => S = 0)
kap, rho, c, Gn = var('kappa rho c G_N', domain='positive')
LPhi, LPsi, LChi = var('LPhi LPsi LChi')
T00 = rho * c**2

# ---------------------------------------------------------------------------
print("\n[B] Williams' Eq 1.34 dust with u4 = 0 and g_44 DYNAMICAL")
print("    u4 = 0  =>  T^44 = gamma u4 u4 = 0,  T^a4 = 0;  dust  =>  T^ij = 0")
solB = solve([LPsi == LPhi + LChi,                 # G_ij = 0  => S = 0
              LPhi - 2 * LPsi == 0,                # G_44 = kappa T_44 = 0
              2 * LPsi - LChi == kap * T00],       # G_00 = kappa T_00
             [LPhi, LPsi, LChi], solution_dict=True)[0]
gammaB = (solB[LPsi] / solB[LPhi]).simplify_full()
chiB = (solB[LChi] / solB[LPhi]).simplify_full()
print("    lap Phi = %s, Psi/Phi = %s, Chi/Phi = %s" % (solB[LPhi], gammaB, chiB))
ok("B1: Psi = Phi/2  =>  PPN gamma = 1/2", gammaB == QQ(1) / 2)
ok("B2: Chi = -Phi/2  (the radion is SOURCED by ordinary matter)", chiB == -QQ(1) / 2)
kap4 = 8 * pi * Gn / c**4
Geff = (solB[LPhi].subs(kappa=kap4) * c**2 / (4 * pi * rho)).simplify_full()
print("    G_eff = %s" % Geff)
ok("B3: G_eff = (4/3) G (absorbed into the constant; gamma is the physics)",
   is_zero(Geff - 4 * Gn / 3))
cassini = 2.3e-5
ok("B4: |gamma - 1| = 0.5 exceeds the Cassini bound 2.3e-5 by > 2e4",
   abs(gammaB - 1) / cassini > 2e4)

# ---------------------------------------------------------------------------
print("\n[C] u4 = 0 along worldlines does not fix g_44")
N = function('N')(x, y, z)
gC = diagonal_matrix(SR, [1 + 2 * Phi, -(1 - 2 * Psi), -(1 - 2 * Psi),
                          -(1 - 2 * Psi), -N**2])
GC = christoffel(gC)
ok("C1: Gamma^4_ab = 0 for a,b in 0..3 (no 4D velocity drives dy/ds)",
   all(is_zero(GC[4][a][b]) for a in range(4) for b in range(4)))
ok("C2: Gamma^4_a4 = d_a N / N (y-equation is homogeneous in dy/ds)",
   all(is_zero(GC[4][a][4] - diff(N, X[a]) / N) for a in range(4)))
print("    => d2y/ds2 = -2 (d_a N/N) u^a dy/ds : dy/ds = 0 stays 0 for ANY N(x)")
print("    Conservation of mass along worldlines is compatible with a free radion.")

# ---------------------------------------------------------------------------
print("\n[D] The lambda term: any radion mass is O(sqrt Lambda)")
Lam = 1.1e-52                         # m^-2, Lambda = lambda/4 (Planck-size)
AU = 1.496e11
rng = 1 / sqrt(Lam)
print("    range ~ 1/sqrt(Lambda) = %.2e m = %.1e AU" % (rng, rng / AU))
ok("D1: range exceeds the solar system by > 1e14: no Yukawa screening",
   rng / AU > 1e14)
print("    (dimensional argument: Lambda is the only scale in the gravitational")
print("     action besides M_P, which multiplies every term.)")

# ---------------------------------------------------------------------------
print("\n[E] Conservation of mass as a RESTRICTION: g_44 = -1 held by a multiplier")
print("    Chi = 0 imposed; the 44 equation fixes the multiplier mu, not the metric.")
solE = solve([LPsi == LPhi,                        # G_ij = 0 with Chi = 0
              2 * LPsi == kap * T00],              # G_00 = kappa T_00
             [LPhi, LPsi], solution_dict=True)[0]
ok("E1: Psi = Phi  =>  PPN gamma = 1 (GR)", (solE[LPsi] / solE[LPhi]) == 1)
lapphiN = (solE[LPhi].subs(kappa=kap4) * c**2).simplify_full()
ok("E2: lap(phi_N) = 4 pi G rho, Newton's coefficient exactly",
   is_zero(lapphiN - 4 * pi * Gn * rho))
mu = ((solE[LPhi] - 2 * solE[LPsi]) / kap).simplify_full()   # G_44 / kappa
print("    multiplier (reaction) T_44 = G_44 / kappa = %s" % mu)
ok("E3: the reaction is T_44 = -rho c^2/2 = (eps/2) T^a_a, the source found by "
   "gr_on_4d_surface (C) and weak_field_poisson (C3)",
   is_zero(mu + rho * c**2 / 2))
print("    Nonlinear: with K = 0 (gr_on_4d_surface) R5_ab = R4_ab, so the surface")
print("    equations are exactly 4D GR and beta = 1 as well.")

# ---------------------------------------------------------------------------
print("\n[F] KK vector g_04")
print("    T_04 = rho u_0 u_4 = 0 under u4 = 0, so G_04 = lap(a)/2 = 0")
ok("F1: G_04 contains no Phi, Psi, Chi: the scalar sector does not drive a",
   not any(Gl[0, 4].has(F) for F in (Phi, Psi, Chi)))
print("    => lap(a) = 0 with decaying data gives a = 0 (no gravito-vector)")

# ---------------------------------------------------------------------------
print("\n[G] Check of weak_field_poisson branch B (g_44 frozen, G_44 = 0 imposed)")
SB = (Phi - Phi / 2)                       # Phi + Chi - Psi with Chi=0, Psi=Phi/2
Gxy = -diff(SB, x, y)
ok("G1: branch B leaves G_xy = -(1/2) d_x d_y Phi != 0 (needs anisotropic stress)",
   not is_zero(Gxy))
print("    => branch B is not a dust solution; the consistent unsourced-44 result")
print("       is [B] (gamma = 1/2 with a dynamical radion).")

# ---------------------------------------------------------------------------
print("\n[VERDICT]")
print("  * Dynamical g_44 + u4 = 0 dust: gamma = 1/2. Excluded by Cassini.")
print("  * u4 = 0 along worldlines does not fix g_44; Lambda cannot screen it.")
print("  * Conservation of mass as a restriction (Williams' word; he writes only")
print("    the surface equations 1.38): gamma = 1, Newton exact, beta = 1, and")
print("    the T_44 = -rho c^2/2 found earlier is the restriction's reaction.")
print("  * Cassini SELECTS the constraint reading; it does not derive it.")
print("\n" + "=" * 74)
print("kk_radion_g44 COMPLETE.")
print("=" * 74)
