#!/usr/bin/env sage
# --------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, "The Dynamic Theory: A New View of
#   Space, Time, and Matter," LA-8370-MS (LASL, 1980). The 5D space-time-mass
#   manifold with signature (+,-,-,-,-) and the mass coordinate q4 = gamma/a0.
#   The weak-field question addressed here is posed in the Layer-3 mass-sector
#   addendum (05_GR-QM/mass_sector_note.tex, sec.1, check tag [C-wf]), which
#   claims the linearisation "yields" Poisson's equation with Newton's sign and
#   coefficient, in contrast to the Schwarzschild_4 x R embedding (T2.14) which
#   only "admits" the general-relativistic solution.
# --------------------------------------------------------------------------

# ============================================================================
# 39_weak_field_poisson.sage                     [A]  (assertion-hardened)
#
# THE "YIELD" HALF OF THE GR CORRESPONDENCE, CHECKED.  [C-wf]
#
# Linearise the 5D metric about flat eta = diag(+,-,-,-,-) on the static
# Newtonian ansatz
#       h_00 = 2 Phi(x,y,z),   h_ij = 2 Psi(x,y,z) delta_ij,   h_4mu = 0,
# and ask what the field equations give.
#
# SCOPE (stated up front, load-bearing).  DT does NOT postulate 5D Einstein
# equations: in the mass sector the potential phi_m is a GAUGE object of the
# eight extended Maxwell equations (Book Eq. 2.9/2.10), and "no Einstein
# equation enters" (mass_sector_note sec.2).  What is checked here is therefore
# the GEOMETRIC correspondence only: IF the 5D system is governed by sourced
# Einstein equations G_ij = kappa T_ij, THEN what does the weak field give?
# The Ricci-flat vacuum case IS DT-operative (it is the equation the T2.14
# Schwarzschild_4 x R embedding solves); the sourced case is an added
# assumption, flagged as such in the verdict.
#
# WHAT THE SCRIPT ESTABLISHES
#   (A) The three linearised identities claimed by [C-wf]:
#         R_00 = lap(Phi),  R_44 = 0,  G_00 = 2 lap(Psi).
#   (B) The full linearised set the note does not give: R_ij, R, G_44.
#   (C) Branch A (impose Psi = Phi, the 4D GR relation): Newton's sign AND
#       coefficient are recovered, lap(phi_N) = 4 pi G rho.  BUT the SAME
#       ansatz then forces a mass-direction source T_44 = -rho c^2/2, so
#       Psi = Phi is an EXTRA assumption, not a consequence of the 5D
#       equations with T_00-only matter.
#   (D) Branch B (impose T_44 = 0, a genuinely unsourced fifth direction):
#       the equations force Psi = Phi/2 and lap(phi_N) = 8 pi G rho, TWICE
#       Newton -- the standard 5D Kaluza-Klein deficit.
#   (E) In VACUUM (rho = 0) the two branches coincide and R_44 = 0 is
#       consistent, which is exactly the "admit = vacuum + spectator limit of
#       yield" statement of mass_sector_note sec.1 (T4').
#
# NET: the Newtonian limit is recovered from the 5D geometry, but only with
# Psi = Phi imposed by hand; it is NOT a free-standing consequence of a
# T_00-only source in five dimensions.  The correspondence is therefore
# conditional on that relation, and that condition should be stated wherever
# the weak-field result is cited.
#
# RUN: /home/dsteele/miniconda3/envs/sage/bin/sage 39_weak_field_poisson.sage
# ============================================================================

from sage.all import *


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


print("=" * 74)
print("39_weak_field_poisson.sage  --  [C-wf] the yield half of the GR bridge")
print("=" * 74)

# ---------------------------------------------------------------------------
# (0) the 5D manifold, the flat background, and the static Newtonian ansatz
# ---------------------------------------------------------------------------
print("\n[0] 5D manifold, signature (+,-,-,-,-), static Newtonian ansatz")

M = Manifold(5, 'M', structure='Lorentzian')
X = M.chart('t x y z q')
t, x, y, z, q = X[:]

eps = var('eps')                      # linearisation parameter
Phi = function('Phi')(x, y, z)        # h_00 = 2 Phi   (time potential)
Psi = function('Psi')(x, y, z)        # h_ij = 2 Psi d_ij (space potential)

g = M.metric()
g[0, 0] = 1 + 2 * eps * Phi
g[1, 1] = -1 + 2 * eps * Psi
g[2, 2] = -1 + 2 * eps * Psi
g[3, 3] = -1 + 2 * eps * Psi
g[4, 4] = -1                          # fifth row/column strictly flat, h_4mu = 0

print("    g_00 = 1 + 2 eps Phi ; g_ij = -(1 - 2 eps Psi) d_ij ; g_44 = -1")
print("    Phi, Psi = Phi(x,y,z), Psi(x,y,z)  (static: no t, no q dependence)")


def lin(expr):
    """the O(eps) coefficient: the linearised part of expr."""
    return SR(expr).taylor(eps, 0, 1).coefficient(eps, 1).expand()


def lap(f):
    """flat 3D Laplacian."""
    return (diff(f, x, 2) + diff(f, y, 2) + diff(f, z, 2)).expand()


Ric = g.ricci()
Rs = g.ricci_scalar()

# Einstein tensor G_ij = R_ij - (1/2) R g_ij, linearised
def G_lin(i, j):
    return lin(Ric[i, j].expr() - Rs.expr() / 2 * g[i, j].expr())


R00 = lin(Ric[0, 0].expr())
R44 = lin(Ric[4, 4].expr())
R11 = lin(Ric[1, 1].expr())
R12 = lin(Ric[1, 2].expr())
Rsc = lin(Rs.expr())
G00 = G_lin(0, 0)
G44 = G_lin(4, 4)

# ---------------------------------------------------------------------------
# (A) the three identities claimed by [C-wf]
# ---------------------------------------------------------------------------
print("\n[A] The three linearised identities claimed by [C-wf]")
print("    R_00 = %s" % R00)
ok("A1: R_00 = lap(Phi).", (R00 - lap(Phi)).is_trivial_zero()
   or (R00 - lap(Phi)).simplify_full() == 0)
print("    R_44 = %s" % R44)
ok("A2: R_44 = 0 (the fifth direction is an unsourced spectator at O(eps)).",
   R44 == 0)
print("    G_00 = %s" % G00)
ok("A3: G_00 = 2 lap(Psi).", (G00 - 2 * lap(Psi)).simplify_full() == 0)

# ---------------------------------------------------------------------------
# (B) the rest of the linearised set (not given in the note)
# ---------------------------------------------------------------------------
print("\n[B] The remaining linearised components")
print("    R_11 = %s" % R11)
ok("B1: R_ij = d_i d_j (Psi - Phi) + delta_ij lap(Psi)  [i=j=1 component].",
   (R11 - (diff(Psi - Phi, x, 2) + lap(Psi))).simplify_full() == 0)
print("    R_12 = %s" % R12)
ok("B2: off-diagonal R_12 = d_1 d_2 (Psi - Phi).",
   (R12 - diff(Psi - Phi, x, y)).simplify_full() == 0)
print("    R    = %s" % Rsc)
ok("B3: Ricci scalar R = 2 lap(Phi) - 4 lap(Psi).",
   (Rsc - (2 * lap(Phi) - 4 * lap(Psi))).simplify_full() == 0)
print("    G_44 = %s" % G44)
ok("B4: G_44 = lap(Phi) - 2 lap(Psi)  (NOT zero unless lap(Phi) = 2 lap(Psi)).",
   (G44 - (lap(Phi) - 2 * lap(Psi))).simplify_full() == 0)

# ---------------------------------------------------------------------------
# (C) branch A: impose Psi = Phi (the 4D GR relation) -> Newton, plus a cost
# ---------------------------------------------------------------------------
print("\n[C] Branch A: impose Psi = Phi (the relation 4D GR itself satisfies)")

G_newton, c_light, rho = var('G_newton c_light rho', domain='positive')
kappa = 8 * pi * G_newton / c_light**4
T00 = rho * c_light**2                     # energy density
print("    kappa = 8 pi G / c^4 ; T_00 = rho c^2")

# G_00 = 2 lap(Psi) -> with Psi = Phi -> 2 lap(Phi) = kappa T_00
lapPhi_A = solve(2 * var('L') == kappa * T00, var('L'))[0].rhs()
print("    lap(Phi) = %s" % lapPhi_A)
ok("C1: Psi=Phi gives lap(Phi) = 4 pi G rho / c^2.",
   (lapPhi_A - 4 * pi * G_newton * rho / c_light**2).simplify_full() == 0)
# metric potential Phi is the Newtonian potential over c^2 (g_00 = 1 + 2 phi_N/c^2)
lap_phiN_A = (lapPhi_A * c_light**2).simplify_full()
print("    => lap(phi_N) = %s   (phi_N = c^2 Phi)" % lap_phiN_A)
ok("C2: NEWTON's sign AND coefficient recovered: lap(phi_N) = 4 pi G rho.",
   (lap_phiN_A - 4 * pi * G_newton * rho).simplify_full() == 0)

# the cost: G_44 no longer vanishes where rho != 0
G44_A = (lapPhi_A - 2 * lapPhi_A).simplify_full()        # lap(Psi)=lap(Phi)
T44_A = (G44_A / kappa).simplify_full()
print("    but G_44 = lap(Phi) - 2 lap(Phi) = %s" % G44_A)
print("    => required mass-direction source T_44 = G_44/kappa = %s" % T44_A)
ok("C3: Psi=Phi FORCES a non-zero mass-direction source T_44 = -rho c^2/2.",
   (T44_A + rho * c_light**2 / 2).simplify_full() == 0)
ok("C4: that source is non-zero wherever rho > 0 (so Psi=Phi is an EXTRA "
   "assumption, not a consequence of T_00-only matter in 5D).",
   not bool(T44_A.subs(rho=1, c_light=1) == 0))

# ---------------------------------------------------------------------------
# (D) branch B: impose T_44 = 0 (genuinely unsourced fifth direction)
# ---------------------------------------------------------------------------
print("\n[D] Branch B: impose T_44 = 0 (a genuinely unsourced fifth direction)")
# T_44 = 0  <=>  G_44 = 0  <=>  lap(Phi) = 2 lap(Psi)
# with G_00 = 2 lap(Psi) = kappa T_00
lapPsi_B = (kappa * T00 / 2).simplify_full()
lapPhi_B = (2 * lapPsi_B).simplify_full()
print("    G_44 = 0  =>  lap(Phi) = 2 lap(Psi)  =>  Psi = Phi/2 (not Psi = Phi)")
print("    lap(Psi) = %s" % lapPsi_B)
print("    lap(Phi) = %s" % lapPhi_B)
lap_phiN_B = (lapPhi_B * c_light**2).simplify_full()
print("    => lap(phi_N) = %s" % lap_phiN_B)
ok("D1: T_44=0 gives lap(phi_N) = 8 pi G rho, TWICE the Newtonian coefficient.",
   (lap_phiN_B - 8 * pi * G_newton * rho).simplify_full() == 0)
ok("D2: the two branches disagree by exactly a factor 2 in the sourced region.",
   (lap_phiN_B / lap_phiN_A).simplify_full() == 2)

# ---------------------------------------------------------------------------
# (E) vacuum: the branches coincide -> the admit slice is the vacuum limit
# ---------------------------------------------------------------------------
print("\n[E] Vacuum (rho = 0): the two branches coincide")
ok("E1: rho=0 sends both branches to lap(phi_N) = 0 (Laplace, exterior Newton).",
   bool(lap_phiN_A.subs(rho=0) == 0) and bool(lap_phiN_B.subs(rho=0) == 0))
ok("E2: in vacuum G_44 = 0 is consistent with R_44 = 0 and R = 0.",
   bool((lapPhi_A.subs(rho=0) - 2 * lapPhi_A.subs(rho=0)) == 0))
print("    => the Schwarzschild_4 x R embedding (T2.14) is the vacuum +")
print("       spectator corner of this problem, as mass_sector_note sec.1 (T4')")
print("       states; the sourced interior is where the derivational content is,")
print("       and that is where the extra relation Psi = Phi has to be supplied.")

# ---------------------------------------------------------------------------
print("\n[VERDICT]")
print("  [C-wf] CONFIRMED as stated: R_00 = lap(Phi), R_44 = 0, G_00 = 2 lap(Psi),")
print("  and with Psi = Phi the weak field reproduces Poisson's equation with")
print("  Newton's sign and coefficient, lap(phi_N) = 4 pi G rho.  The note's chain")
print("  is correct.")
print()
print("  SCOPE the note does not state, established here:")
print("   (i) Psi = Phi is an ASSUMPTION, not a consequence.  In 5D it forces a")
print("       mass-direction source T_44 = -rho c^2/2, so a T_00-only source is")
print("       inconsistent with it.")
print("  (ii) Imposing T_44 = 0 instead forces Psi = Phi/2 and lap(phi_N) =")
print("       8 pi G rho, twice Newton (the 5D Kaluza-Klein deficit).")
print(" (iii) The two agree only in VACUUM, which is precisely the regime the")
print("       Schwarzschild x R embedding already covers.  So 'yield' improves on")
print("       'admit' only under the added relation Psi = Phi.")
print("  (iv) Throughout, sourced 5D Einstein equations are an ADDED assumption:")
print("       DT's mass sector derives phi_m as a gauge object of the extended")
print("       Maxwell set, with no Einstein equation.  This is a check of the")
print("       geometric correspondence, not a DT-native derivation.")

print("\n" + "=" * 74)
print("39_weak_field_poisson COMPLETE.")
print("=" * 74)
