#!/usr/bin/env sage
# --------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, "The Dynamic Theory" (2011), Eqs (1.34)-
#   (1.39) and their Ch.10 repeat (10.1)-(10.9).  Williams builds the surface
#   energy-momentum tensor (Eq 1.37/10.7) and then "suggests writing"
#   G = C T (Eq 1.38/10.8), "when Equations (1.38) are taken to be Einstein's
#   field equations".  u4 = dy4/dt = 0 (Eq 1.36/10.6) is conservation of mass.
#   The force law used is DT's (Report eqs 66-67; force density C5).
# Full derivation, context and provenance:
#   Supplementary_Material.pdf -- §S3, §6.1 (layer boundary), Appendix A.
# Companion: gr_on_4d_surface.sage (u4 = 0 <=> K = 0, Gauss-Codazzi).
# --------------------------------------------------------------------------

# ============================================================================
# einstein_equations_recovery.sage               [A]  (assertion-hardened)
#
# RECOVERING THE EINSTEIN EQUATIONS ON THE MASS HYPERSURFACE.
#
# Williams' step G = C T is an identification motivated by the form of the
# surface energy-momentum tensor.  The implicit argument is: the source is
# conserved, and G is the conserved curvature tensor.  Lovelock's theorem
# makes that argument exact: in FOUR dimensions the only local, symmetric,
# divergence-free tensor built from the metric and its derivatives up to
# second order is  a G_ab + Lambda g_ab.  The chain checked here:
#
#   [1] CONSERVATION.  u4 = 0 makes the surface source conserved.
#       In a local frame with no mass-direction dependence (cylinder, the
#       field version of K = 0), the 5D gauge stress splits into the 4D
#       Maxwell stress plus the stress of a massless scalar phi = A_4.  With
#       the current carried by matter (J^k ~ u^k), u4 = 0 gives J^4 = 0,
#       which is exactly the wave equation box(phi) = 0, so the scalar stress
#       is separately conserved; the Maxwell stress balances the 4D Lorentz
#       force, which DT's force law hands to matter.  Total surface T is
#       conserved.  With u4 != 0 the surface tensor leaks (mutation).
#   [2] UNIQUENESS (restricted Lovelock).  Among a R_ab + b R g_ab + l g_ab,
#       divergence-free for a metric with non-constant R iff b = -a/2, with
#       l (= Lambda) free.  Contracted Bianchi checked explicitly.
#   [3] WHY u4 = 0 MATTERS FOR UNIQUENESS.  In 5D Lovelock also allows the
#       Gauss-Bonnet tensor H_ab, so 5D field equations are NOT unique.
#       H_ab vanishes identically in 4D, is non-zero in 5D, and on a K = 0
#       product (4D surface x mass line) its 4D block vanishes: whatever
#       second-order 5D equations DT has, the surface sees Einstein + Lambda.
#       The Gauss-Bonnet term survives only in the 44 constraint.
#   [4] COUPLING.  Linearised 4D equations: the off-diagonal spatial
#       equations force Psi = Phi (in 4D this is DERIVED; in 5D it was an
#       extra assumption, C12), and G_00 = kappa rho gives Newton's
#       lap(Phi) = 4 pi G rho iff kappa = 8 pi G/c^4.
#
# VERDICT (honest bracket).  Given u4 = 0 (=> K = 0 and J^4 = 0), DT's force
# law, and the Lovelock premises -- field equations LOCAL, SECOND ORDER, with
# the metric as the only gravitational field -- the surface field equations
# are G_ab + Lambda g_ab = (8 pi G/c^4) T_ab.  Imported, not re-proven:
# Lovelock's theorem (only its restricted form is checked here).  Not fixed:
# Lambda (DT's cosmic acceleration comes from its (1 - H0 t) time-dependence,
# so DT would set Lambda = 0).  Part [1] is a local flat-frame check; the
# curved version follows by minimal coupling.  Structural result about the
# recovery step; it does not test DT's physics.
#
# RUN: /home/dsteele/miniconda3/envs/sage/bin/sage einstein_equations_recovery.sage
# ============================================================================

from sage.all import *


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


def z(e):
    return bool(SR(e).simplify_full() == 0)


def zx(e):                      # polynomial in derivatives: expand suffices
    return bool(SR(e).expand() == 0)


def find_sign(test):
    """return s in (1, -1) with test(s) True, or None"""
    for s in (1, -1):
        if test(s):
            return s
    return None


print("=" * 72)
print("RECOVERING THE EINSTEIN EQUATIONS on the mass hypersurface (u4 = 0)")
print("=" * 72)

# ============================================================================
# [1] CONSERVATION OF THE SURFACE SOURCE
# ============================================================================
print("\n[1] Conservation of the surface source (flat local frame, d_4 = 0)")
X = list(var('x0 x1 x2 x3'))
eta5 = diagonal_matrix(SR, [1, -1, -1, -1, -1])
eta4 = diagonal_matrix(SR, [1, -1, -1, -1])
A = [function('A%d' % j)(*X) for j in range(5)]


def d(e, i):
    return diff(e, X[i]) if i < 4 else SR(0)      # cylinder: nothing depends on y4


def maxwell(n, eta, Acomp):
    F = matrix(SR, n, n, lambda i, j: d(Acomp[j], i) - d(Acomp[i], j))
    Fup = eta * F * eta                            # F^{ij}
    Fmix = Fup * eta                               # F^i_k
    FF = sum(F[k, l] * Fup[k, l] for k in range(n) for l in range(n))
    T = Fmix * Fup + (1 / 4) * eta * FF           # Rpt/Book Eq 1.34 field part
    J = [sum(d(Fup[i, k], i) for i in range(n)) for k in range(n)]
    return F, Fup, Fmix, T, J


F5, Fup5, Fmix5, T5, J5 = maxwell(5, eta5, A)
div5 = [sum(d(T5[i, j], i) for i in range(5)) for j in range(5)]
force5 = [sum(Fmix5[j, k] * J5[k] for k in range(5)) for j in range(5)]
s_f = find_sign(lambda s: all(zx(div5[j] - s * force5[j]) for j in range(5)))
ok("5D: d_i T^ij = %sF^j_k J^k (force-density identity, C5)"
   % ("" if s_f == 1 else "-"), s_f is not None)

F4, Fup4, Fmix4, T4, J4 = maxwell(4, eta4, A[:4])
phi = A[4]
dphi_up = [eta4[a, a] * diff(phi, X[a]) for a in range(4)]
dphi2 = sum(eta4[a, a] * diff(phi, X[a])**2 for a in range(4))
TS = matrix(SR, 4, 4, lambda a, b: dphi_up[a] * dphi_up[b] - eta4[a, b] * dphi2 / 2)
sig = find_sign(lambda s: all(zx(T5[a, b] - T4[a, b] - s * TS[a, b])
                              for a in range(4) for b in range(4)))
ok("cylinder: 4D block of the 5D stress = Maxwell stress(F_ab) %s scalar stress(phi = A_4)"
   % ("+" if sig == 1 else "-"), sig is not None)
ok("scalar energy density T_S^00 = (1/2)(phi_t^2 + |grad phi|^2) >= 0 (sign +)",
   sig == 1)

box = sum(eta4[a, a] * diff(phi, X[a], 2) for a in range(4))
tau = find_sign(lambda s: zx(J5[4] - s * box))
ok("J^4 = %sbox(phi): the mass-direction current is the scalar's wave operator"
   % ("" if tau == 1 else "-"), tau is not None)
ok("4D current J^a (5D, d_4 = 0) equals the 4D Maxwell current",
   all(zx(J5[a] - J4[a]) for a in range(4)))

divS = [sum(diff(TS[a, b], X[a]) for a in range(4)) for b in range(4)]
ok("d_a T_S^ab = box(phi) d^b phi  ->  conserved iff J^4 = 0, i.e. iff u4 = 0",
   all(zx(divS[b] - box * dphi_up[b]) for b in range(4)))
ok("mutation: with u4 != 0 (J^4 != 0) the scalar stress is NOT conserved",
   not all(zx(divS[b]) for b in range(4)))
div4 = [sum(diff(T4[a, b], X[a]) for a in range(4)) for b in range(4)]
ok("d_a T_EM^ab = %sF^b_a J^a: the 4D Lorentz force, balanced by matter (DT force law)"
   % ("" if s_f == 1 else "-"),
   all(zx(div4[b] - s_f * sum(Fmix4[b, a] * J4[a] for a in range(4)))
       for b in range(4)))
print("    => with J^k ~ u^k, u4 = 0 gives J^4 = 0: the A_4 sector is source-free on")
print("       the surface, so the Newtonian tail must come from the METRIC (curvature),")
print("       consistent with C10's assignment of -GM/r to curvature.")

# ============================================================================
# [2] RESTRICTED LOVELOCK: which combination is divergence-free?
# ============================================================================
print("\n[2] Uniqueness among a R_ab + b R g_ab + l g_ab (4D, explicit diagonal metric)")
t, x = var('t x')


def diag_metric(name, n, entries):
    Mn = Manifold(n, 'M_' + name)
    ch = Mn.chart(' '.join(['t', 'x', 'u', 'v', 'w'][:n]))
    g = Mn.metric('g_' + name, signature=2 - n)
    for i in range(n):
        g[i, i] = entries[i]
    return Mn, ch, g


# Explicit, non-symmetric diagonal metric of (t, x).  Fully generic functions
# exhaust memory (Bianchi divergence here, curvature-squared contractions in
# [3]); the identities checked are known theorems, so a non-trivial instance
# plus a mutation test is the check.
Pe, Qe, Re, Se = exp(t * x), 1 + t**2 + x**2, (2 + t) * (1 + x**2), exp(t) * (1 + x**4)
M4, ch4, g4 = diag_metric('a', 4, [Pe, -Qe, -Re, -Se])
Ric = g4.ricci()
Rsc = g4.ricci_scalar()
ginv = [1 / g4[i, i].expr() for i in range(4)]
DRic = g4.connection()(Ric)                      # DRic[a,b,c] = nabla_c R_ab
co = ch4[:]
dR = [diff(Rsc.expr(), co[b]) for b in range(4)]
divRic = [sum(ginv[a] * DRic[a, b, a].expr() for a in range(4)) for b in range(4)]
ok("contracted Bianchi: nabla^a R_ab = (1/2) d_b R", all(z(divRic[b] - dR[b] / 2) for b in range(4)))
ok("R is not constant on this family (the test is non-trivial)", not z(dR[0]))
al, be, lm = var('al be lm')
divE = [al * divRic[b] + be * dR[b] for b in range(4)]   # nabla g = 0: lm drops
ok("divergence of a R_ab + b R g_ab + l g_ab = (a/2 + b) d_b R for any l",
   all(z(divE[b] - (al / 2 + be) * dR[b]) for b in range(4)))
ok("divergence-free  <=>  b = -a/2 (the Einstein tensor), Lambda = l free",
   solve(al / 2 + be == 0, be)[0].rhs() == -al / 2)
ok("mutation: b = -a/4 is not divergence-free",
   not z(divE[0].subs(be=-al / 4)))

# ============================================================================
# [3] GAUSS-BONNET: allowed in 5D, absent in 4D, absent on a K = 0 surface
# ============================================================================
print("\n[3] The Lovelock (Gauss-Bonnet) tensor H_ab in 4D, 5D, and on the K = 0 surface")


def lovelock_H(g, n, c_mut=-2):
    gd = [g[i, i].expr() for i in range(n)]
    gi = [1 / e for e in gd]
    Rm = g.riemann()
    Ru = [[[[Rm[m, b, c, dd].expr() for dd in range(n)] for c in range(n)]
           for b in range(n)] for m in range(n)]
    Rd = [[[[gd[a] * Ru[a][b][c][dd] for dd in range(n)] for c in range(n)]
           for b in range(n)] for a in range(n)]
    Rc = [[sum(Ru[m][a][m][b] for m in range(n)) for b in range(n)] for a in range(n)]
    Rs = sum(gi[a] * Rc[a][a] for a in range(n))
    RcRc = sum(gi[a] * gi[b] * Rc[a][b]**2 for a in range(n) for b in range(n))
    RR = sum(gi[a] * gi[b] * gi[c] * gi[dd] * Rd[a][b][c][dd]**2
             for a in range(n) for b in range(n) for c in range(n) for dd in range(n))
    GB = Rs**2 - 4 * RcRc + RR
    H = matrix(SR, n, n)
    for a in range(n):
        for b in range(a, n):
            e = (Rs * Rc[a][b]
                 - 2 * sum(Rc[a][c] * gi[c] * Rc[c][b] for c in range(n))
                 + c_mut * sum(Rd[a][c][b][dd] * gi[c] * gi[dd] * Rc[c][dd]
                               for c in range(n) for dd in range(n))
                 + sum(gi[c] * gi[dd] * gi[ee] * Rd[a][c][dd][ee] * Rd[b][c][dd][ee]
                       for c in range(n) for dd in range(n) for ee in range(n)))
            H[a, b] = 2 * e - (gd[a] if a == b else 0) * GB / 2
            H[b, a] = H[a, b]
    return H, GB


g4e = g4
H4, GB4 = lovelock_H(g4e, 4)
print("    4D curvature-squared contractions done")
ok("the test metric has non-zero Gauss-Bonnet density (instance is non-trivial)",
   not z(GB4))
ok("4D: H_ab = 0 identically (Gauss-Bonnet adds nothing to 4D field equations)",
   all(z(H4[a, b]) for a in range(4) for b in range(a, 4)))
H4bad, _ = lovelock_H(g4e, 4, c_mut=-1)
ok("mutation: a wrong coefficient in H_ab is not identically zero in 4D",
   not all(z(H4bad[a, b]) for a in range(4) for b in range(a, 4)))

aF = function('aF')(t)
M5f, ch5f, g5f = diag_metric('f', 5, [1] + [-aF**2] * 4)
H5f, _ = lovelock_H(g5f, 5)
ok("5D: H_ab != 0 (FRW-type 5D metric) -> 5D second-order equations are NOT unique",
   not z(H5f[0, 0]))
Ht = M5f.tensor_field(0, 2, name='H', sym=(0, 1))
for a in range(5):
    for b in range(a, 5):
        Ht[a, b] = H5f[a, b]
DH = g5f.connection()(Ht)
gi5 = [1 / g5f[i, i].expr() for i in range(5)]
ok("5D: nabla^a H_ab = 0 (so it is a legitimate Lovelock term there)",
   all(z(sum(gi5[a] * DH[a, b, a].expr() for a in range(5))) for b in range(5)))

M5p, ch5p, g5p = diag_metric('p', 5, [Pe, -Qe, -Re, -Se, -1])   # 4D surface x mass line, K = 0
H5p, _ = lovelock_H(g5p, 5)
ok("K = 0 product: the 4D block of the 5D Gauss-Bonnet tensor vanishes",
   all(z(H5p[a, b]) for a in range(4) for b in range(a, 4)))
ok("K = 0 product: H_a4 = 0 and H_44 = -(1/2) g_44 GB(4D) -- only the 44 constraint feels it",
   all(z(H5p[a, 4]) for a in range(4)) and z(H5p[4, 4] - GB4 / 2))

# ============================================================================
# [4] THE COUPLING: Newtonian limit on the surface
# ============================================================================
print("\n[4] Coupling constant from the Newtonian limit (4D, linearised)")
y_, zc, ep = var('y_ zc ep')
Phi = function('Phi')(x, y_, zc)
Psi = function('Psi')(x, y_, zc)
Mw = Manifold(4, 'Mw')
chw = Mw.chart('t x y_ zc')
gw = Mw.metric('gw', signature=-2)
gw[0, 0] = 1 + 2 * ep * Phi
for i in (1, 2, 3):
    gw[i, i] = -(1 - 2 * ep * Psi)
Gw = gw.ricci() - gw.ricci_scalar() * gw / 2


def lin(e):
    return diff(e, ep).subs(ep=0)


s_xy = find_sign(lambda s: z(lin(Gw[1, 2].expr()) - s * diff(Phi - Psi, x, y_)))
ok("G_xy = %sd_x d_y (Phi - Psi): vacuum spatial equations force Psi = Phi (DERIVED in 4D)"
   % ("" if s_xy == 1 else "-"), s_xy is not None)
lapPhi = sum(diff(Phi, c_, 2) for c_ in (x, y_, zc))
G00 = lin(Gw[0, 0].expr()).substitute_function(Psi.operator(), Phi.operator())
s00 = find_sign(lambda s: z(G00 - s * 2 * lapPhi))
ok("with Psi = Phi: G_00 = %s2 lap(Phi)" % ("" if s00 == 1 else "-"), s00 is not None)
kap, Gn, rho = var('kap Gn rho')
ksol = solve(s00 * 2 * (4 * pi * Gn * rho) == kap * rho, kap)[0].rhs()
ok("G_00 = kappa rho reproduces lap(Phi) = 4 pi G rho iff kappa = 8 pi G (c = 1)",
   bool(ksol == 8 * pi * Gn))
ok("mutation guard: G_00 = %s2 lap(Psi) depends on Psi ONLY -- perturbing g_00 alone "
   "(Psi = 0) sources nothing" % ("" if s00 == 1 else "-"),
   z(lin(Gw[0, 0].expr()) - s00 * 2 * sum(diff(Psi, c_, 2) for c_ in (x, y_, zc))))
print("    (G_00 depends on Psi only: the g_00 perturbation alone sources nothing;")
print("     the 4D spatial equations supply Psi = Phi, and with it Newton.)")
print("    Williams writes G = -(8 pi K/c^2) T; the sign tracks the curvature")
print("    convention, the magnitude 8 pi G (per c^4 in SI) is what is checked.")

print("""
  NET: given u4 = 0 (=> K = 0, J^4 = 0), DT's force law, and Lovelock's
  premises (local, second order, metric-only), the surface field equations are
        G_ab + Lambda g_ab = (8 pi G / c^4) T_ab ,
  with T_ab conserved, Psi = Phi derived, and any 5D Gauss-Bonnet term
  confined to the 44 constraint.  Imported: Lovelock's theorem (restricted
  form checked).  Open: Lambda (DT would set it to 0).""")
print("=" * 72)
print("einstein_equations_recovery COMPLETE.")
print("=" * 72)
