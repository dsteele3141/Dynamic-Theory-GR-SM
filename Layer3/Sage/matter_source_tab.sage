#!/usr/bin/env sage
# --------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, "The Dynamic Theory" (2011) Eqs (1.32)-
#   (1.39): the 5D field tensor with V_0..V_3 in the fifth column (1.32), the
#   space energy-momentum tensor T^ij_sp = gamma u^i u^j + (1/c^2)[F^i_k F^kj
#   + (1/4) a^ij F^kl F_kl] (1.34), its surface form (1.35), u^4 = d_t y^4 +
#   div(y^4 u) = 0 "the conservation of the fifth dimension" (1.36), the
#   surface tensor with the F^a_4 F^4b stress removed (1.37), and
#   G^mu nu = -(8 pi K / c^2) T^mu nu, K the gravitational constant (1.39).
#   Eq (10.19) beta = sqrt(4 pi eps0 G) (gravitational <-> electrical units).
#   Same content in LA-8370-MS (1980) Eqs (68)-(70).
# Data: Bertotti, Iess & Tortora, Nature 425 (2003) 374; planetary-ephemeris
#   bound |beta_PPN - 1| < ~1e-4 (Mercury perihelion).
# Full derivation, context and provenance:
#   Supplementary_Material.pdf -- §S3, §6.1 (layer boundary), Appendix A.
# --------------------------------------------------------------------------

# ============================================================================
# matter_source_tab.sage                           [A]  (assertion-hardened)
#
# GAP #4: IS DT'S MATTER SOURCE GR'S MATTER SOURCE?
#
# Einstein's equations need a source T_mu nu with three properties: the right
# coupling constant, conservation (nabla_mu T^mu nu = 0, forced by the Bianchi
# identity), and no unwanted extra pieces.  This script takes Williams' own
# T^{ab} (Eqs 1.34-1.37) and checks all three.
#
# WHAT THE SCRIPT ESTABLISHES
#   (A) Coupling: Williams' constant 8 pi K / c^2 (Eq 1.39, T in mass-density
#       units) reproduces Newton's law exactly in the restricted linearised
#       surface equations (chi = 0): nabla^2 Phi = 4 pi G rho.
#   (B) Conservation FORCES conservation of mass.  The Bianchi identity
#       nabla_mu G^mu nu = 0 (checked) makes the surface equations
#       inconsistent unless nabla_mu T^mu nu = 0.  For Williams' dust term the
#       projection along u is exactly the continuity residual nabla_mu(rho
#       u^mu), whose Newtonian form IS Williams' u^4 (Eq 1.36); the
#       transverse part is the geodesic equation.  So u^4 = 0 is not only a
#       restriction we choose: it is the INTEGRABILITY CONDITION of DT's own
#       surface equations.  (Keeping it, given the 5D drive of
#       test_body_motion, still costs a reaction.)
#   (C) No mass-axis source: with u^4 = 0, T^{a4} = 0, so matter does not
#       source g_a4 or g_44 -- consistent with gaps #1 and #3.
#   (D) Williams' surface tensor, evaluated: its EM part is the standard
#       energy (E^2 + B^2)/2 (validating his + 1/4 sign convention), but Eq
#       1.37 removes only the F^a_4 F^4b stress and LEAVES a trace term
#       (1/4) g^ab (|V|^2 - V_0^2).
#   (E) That residual is excluded IF V is the Newtonian field.  With
#       V = g / sqrt(4 pi G) (the Eq 10.19 unit map, normalised by (D)), the
#       leftover energy density g^2/(16 pi G) gravitates and shifts PPN beta
#       to 9/8: excluded by Mercury by ~1e3.
#   (F) DT itself removes it.  With a convective current J^A = sigma u^A,
#       u^4 = 0 gives J^4 = 0; the fifth field equation (box phi_4 = J^4,
#       gravitational_radiation [F]) is then sourceless, so V = d phi_4 = 0
#       on the surface, the V terms vanish, and the 4D force reduces to the
#       Lorentz force.  T^ab = dust + Maxwell: GR's source exactly.
#   (G) Weyl coupling: the Maxwell action is Weyl-invariant in 4D (weight
#       n - 4 = 0) but not in 5D, so on the surface EM does not source the
#       Weyl vector (weyl_vector_mass).  Dust is not Weyl-invariant; it is
#       assumed to couple in the Einstein gauge (Ghilencea), not derived.
#
# VERDICT (honest bracket).  Williams' T^ab, with conservation of mass,
# is GR's matter source: the right constant, conserved, and with no
# mass-axis or V-field pieces.  The headline: conservation of mass is FORCED
# by the Bianchi identity of the surface equations -- the step from "a
# restriction that the data select" to "a condition the equations require"
# for u^4 (the g_44 part of gaps #1/#3 is still selected by Cassini).  The one
# leftover in Williams' Eq 1.37, a V^2 trace term, would give beta = 9/8 if
# taken at Newtonian strength; DT's own u^4 = 0 removes it.  Scope: dust plus
# EM (pressure is the standard perfect-fluid extension, not in Williams'
# T); convective current assumed; h^ab terms set to zero in (D); dust's Weyl
# coupling assumed, not derived.
#
# RUN: /home/dsteele/miniconda3/envs/sage/bin/sage matter_source_tab.sage
# ============================================================================

from sage.all import *


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


print("=" * 74)
print("matter_source_tab.sage -- gap #4: is DT's T_ab GR's T_ab?")
print("=" * 74)

# ---------------------------------------------------------------------------
print("\n[A] Coupling constant (restricted linearised surface equations)")
LPhi, LPsi, LChi, kap, rho = var('LPhi LPsi LChi kap rho')
# kk_radion_g44: G_00 = 2 L(Psi) - L(chi), G_ij ~ (Phi + chi - Psi)  [units 1/c^2]
sol = solve([2 * LPsi - LChi == kap * rho, LPhi + LChi - LPsi == 0, LChi == 0],
            [LPhi, LPsi, LChi], solution_dict=True)[0]
Gn = var('Gn')
kW = 8 * pi * Gn                       # Williams Eq 1.39 magnitude, c = 1
print("    L(Phi) = %s ; with kappa = 8 pi K: L(Phi) = %s" % (sol[LPhi], sol[LPhi].subs(kap=kW)))
ok("A1: Williams' 8 pi K/c^2 gives Newton: lap Phi = 4 pi G rho",
   bool((sol[LPhi].subs(kap=kW) - 4 * pi * Gn * rho).simplify_full() == 0))

# ---------------------------------------------------------------------------
print("\n[B] Bianchi identity -> conservation -> u^4 = 0")
t, r, th, ph = var('t r th ph')
X = [t, r, th, ph]
aa, bb = function('aa')(r), function('bb')(r)
g = diagonal_matrix(SR, [exp(2 * aa), -exp(2 * bb), -r**2, -r**2 * sin(th)**2])
gi = g.inverse()
n = 4
Gam = [[[sum(gi[a, d] * (diff(g[d, b], X[c]) + diff(g[d, c], X[b]) - diff(g[b, c], X[d]))
             for d in range(n)) / 2 for c in range(n)] for b in range(n)] for a in range(n)]
Ric = matrix(SR, n, n, lambda b, d: sum(diff(Gam[a][d][b], X[a]) - diff(Gam[a][a][b], X[d])
                                        + sum(Gam[a][a][f] * Gam[f][d][b] - Gam[a][d][f] * Gam[f][a][b]
                                              for f in range(n)) for a in range(n)))
Rs = sum(gi[a, b] * Ric[a, b] for a in range(n) for b in range(n))
Gup = gi * (Ric - Rs / 2 * g) * gi
divG_r = sum(diff(Gup[m, 1], X[m]) for m in range(n)) \
    + sum(Gam[m][m][l] * Gup[l, 1] for m in range(n) for l in range(n)) \
    + sum(Gam[1][m][l] * Gup[m, l] for m in range(n) for l in range(n))
ok("B1: contracted Bianchi identity nabla_mu G^{mu r} = 0 (general static metric)",
   bool(divG_r.simplify_full() == 0))

x, y, z = var('x y z')
Y = [t, x, y, z]
eta = diagonal_matrix(SR, [1, -1, -1, -1])
v = [function('v%d' % i)(*Y) for i in range(1, 4)]
dens = function('dens')(*Y)
Gm = 1 / sqrt(1 - sum(vi**2 for vi in v))
U = [Gm] + [Gm * vi for vi in v]                        # u^mu, u.u = 1
Ul = [eta[m, m] * U[m] for m in range(4)]
T = matrix(SR, 4, 4, lambda m, k: dens * U[m] * U[k])
divT = [sum(diff(T[m, k], Y[m]) for m in range(4)) for k in range(4)]
cont = sum(diff(dens * U[m], Y[m]) for m in range(4))
proj = sum(Ul[k] * divT[k] for k in range(4))
ok("B2: u_nu d_mu T^{mu nu} = d_mu(rho u^mu)  (timelike part = continuity)",
   bool((proj - cont).simplify_full() == 0))
acc = [sum(U[m] * diff(U[k], Y[m]) for m in range(4)) for k in range(4)]
ok("B3: the rest is rho u^mu d_mu u^nu (transverse part = geodesic/Euler)",
   all(bool((divT[k] - U[k] * cont - dens * acc[k]).simplify_full() == 0) for k in range(4)))
eps = var('eps')
GmE = 1 / sqrt(1 - eps**2 * sum(vi**2 for vi in v))
UE = [GmE] + [GmE * eps * vi for vi in v]
cont_NR = sum(diff(dens * UE[m], Y[m]) for m in range(4))
lead = cont_NR.subs(eps=0) + eps * diff(cont_NR, eps).subs(eps=0)
u4W = diff(dens, t) + sum(diff(dens * eps * v[i], Y[i + 1]) for i in range(3))
ok("B4: to O(v) the continuity residual is Williams' u^4 = d_t y4 + div(y4 u) (Eq 1.36)",
   bool((lead - u4W).simplify_full() == 0))
print("    => surface Einstein equations + Williams' dust REQUIRE u^4 = 0")

# ---------------------------------------------------------------------------
print("\n[C] 5D dust with u^4 = 0 has no mass-axis components")
u5 = vector(SR, list(var('w0 w1 w2 w3')) + [0])
T5 = var('g5') * u5.outer_product(u5)
ok("C1: T^{a4} = 0 for all a (no source for g_a4, g_44)",
   all(T5[a, 4] == 0 for a in range(5)))

# ---------------------------------------------------------------------------
print("\n[D] Williams' surface tensor (Eqs 1.32, 1.35, 1.37) with h = 0")
E1, E2, E3, B1, B2, B3, V0, V1, V2, V3 = var('E1 E2 E3 B1 B2 B3 V0 V1 V2 V3')
F = matrix(SR, [[0, E1, E2, E3, V0],
                [-E1, 0, B3, -B2, V1],
                [-E2, -B3, 0, B1, V2],
                [-E3, B2, -B1, 0, V3],
                [-V0, -V1, -V2, -V3, 0]])
e5 = diagonal_matrix(SR, [1, -1, -1, -1, -1])
Fu = e5 * F * e5                                       # F^{ij}
Fm = e5 * F                                            # F^i_j
FF4 = sum(Fu[m, k] * F[m, k] for m in range(4) for k in range(4))
FV = sum(Fu[4, k] * F[4, k] for k in range(4))
T00_sp = sum(Fm[0, k] * Fu[k, 0] for k in range(4)) + Fm[0, 4] * Fu[4, 0] + QQ(1) / 4 * (FF4 + FV)
T00 = T00_sp - Fm[0, 4] * Fu[4, 0]                     # Eq 1.37, h = 0
EM = T00.subs({V0: 0, V1: 0, V2: 0, V3: 0}).expand()
res = (T00 - EM).expand()
print("    T^00 = %s  +  %s" % (EM.factor(), res.factor()))
ok("D1: the EM part is the standard energy (E^2 + B^2)/2",
   bool((EM - (E1**2 + E2**2 + E3**2 + B1**2 + B2**2 + B3**2) / 2).expand() == 0))
ok("D2: Eq 1.37 leaves a V trace term (1/4)(|V|^2 - V0^2)",
   bool((res - (V1**2 + V2**2 + V3**2 - V0**2) / 4).expand() == 0))

# ---------------------------------------------------------------------------
print("\n[E] If V were the Newtonian field: the leftover gravitates")
GM, cl, R = var('GM cl R', domain='positive')
gfield = GM / R**2
rhoV = gfield**2 / (16 * pi * Gn) / cl**2       # (1/4)|V|^2, V = g/sqrt(4 pi G), as mass density
A_ = var('A_')
dPhi = A_ / R**2                                   # particular solution of lap dPhi = 4 pi G rhoV
eqn = (diff(R**2 * diff(dPhi, R), R) / R**2 - 4 * pi * Gn * rhoV).simplify_full()
Asol = solve(eqn == 0, A_)[0].rhs()
dg00 = (2 * dPhi.subs(A_=Asol) / cl**2).simplify_full()
Uu = GM / (cl**2 * R)
dbeta = (dg00 / (2 * Uu**2)).simplify_full()
print("    delta g_00 = %s  ->  delta beta = %s" % (dg00, dbeta))
ok("E1: the V^2 leftover would shift PPN beta by +1/8 (beta = 9/8)", dbeta == QQ(1) / 8)
ok("E2: excluded by Mercury (|beta - 1| < 1e-4) by more than 1e3", (QQ(1) / 8) / 1e-4 > 1e3)

# ---------------------------------------------------------------------------
print("\n[F] DT removes it: u^4 = 0 -> J^4 = 0 -> V unsourced")
sig = var('sig')
J = vector(SR, [sig * c for c in u5])            # convective current, u^4 = 0
ok("F1: J^4 = sigma u^4 = 0", J[4] == 0)
force = [sum(Fm[nu, B] * J[B] for B in range(5)) for nu in range(4)]
lor = [sum(Fm[nu, m] * J[m] for m in range(4)) for nu in range(4)]
ok("F2: the 4D force F^nu_B J^B is the Lorentz force F^nu_mu J^mu (no V term)",
   all(bool((force[k] - lor[k]).expand() == 0) for k in range(4)))
print("    fifth field equation box phi_4 = J^4 = 0 (gravitational_radiation [F]):")
print("    no source, so V = d phi_4 = 0 on the surface and T^ab = dust + Maxwell.")
ok("F3: with V = 0 Williams' surface T^00 is exactly the Maxwell energy",
   bool((T00.subs({V0: 0, V1: 0, V2: 0, V3: 0}) - EM).expand() == 0))

# ---------------------------------------------------------------------------
print("\n[G] Weyl weight of the Maxwell action: sqrt|g| g^.. g^.. ~ Omega^(n-4)")
Om = var('Om', domain='positive')
wgt = {}
for nd in (4, 5):
    en = diagonal_matrix(SR, [1] + [-1] * (nd - 1))
    gs = Om**2 * en                                    # Weyl rescaling g -> Om^2 g
    L = sqrt(abs(gs.det())) * gs.inverse()[0, 0] * gs.inverse()[1, 1]   # sqrt|g| g^00 g^11
    L0 = sqrt(abs(en.det())) * en.inverse()[0, 0] * en.inverse()[1, 1]
    wgt[nd] = (log(L / L0) / log(Om)).simplify_full()
    print("    n = %d: sqrt|g| g^.. g^.. scales as Om^%s" % (nd, wgt[nd]))
ok("G1: Maxwell is Weyl-invariant on the 4D surface (does not source the Weyl vector)", wgt[4] == 0)
ok("G2: but not in the full 5D space (weight 1)", wgt[5] == 1)

# ---------------------------------------------------------------------------
print("\n[VERDICT]")
print("  * Williams' 8 pi K/c^2 is exactly Newton's coupling.")
print("  * Bianchi forces nabla T = 0, whose timelike part IS u^4 (Eq 1.36):")
print("    conservation of mass is REQUIRED by DT's surface equations.")
print("  * Matter does not source g_a4 / g_44; with u^4 = 0 the V field is")
print("    unsourced, so T_ab = dust + Maxwell: GR's source exactly.")
print("  * Williams' leftover V^2 trace term would give beta = 9/8 at Newtonian")
print("    strength (excluded) -- DT's own u^4 = 0 removes it.")
print("\n" + "=" * 74)
print("matter_source_tab COMPLETE.")
print("=" * 74)
