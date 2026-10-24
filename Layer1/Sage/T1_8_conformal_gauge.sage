#!/usr/bin/env sage
# T1.8 - Weyl conformal gauge: no real rescaling of physical lengths
#
# Source: P. E. Williams, LA-8370-MS (LASL, 1980), eq 37.
# Supplementary Materials: step [L7], Sec. S2 (code index: Sec. S6).
#
# Weyl objection residual check: the conformal gauge introduces NO real
# (observable) rescaling of physical lengths.
#
# Weyl geometry has a gauge freedom (eq 37): ĝ -> f·ĝ together with
# phi_k -> phi_k + (1/2) d_k ln f, under which the intrinsic geometry is
# unchanged and only the local length UNIT changes (ĥ̂_ell^2 = f·ell^2). The
# worry for the second-clock objection: does this rescaling produce a real,
# observable, history-dependent effect? It must not. We verify the four facts
# that establish it does not:
#
#   (1) the field strength F_ij = d_i phi_j - d_j phi_i is GAUGE-INVARIANT
#       (phi -> phi + dLambda leaves F unchanged) — the observable field is fixed;
#   (2) the closed-loop HOLONOMY ∮ phi is GAUGE-INVARIANT (∮ dLambda = 0) — the
#       second-clock-relevant quantity is fixed;
#   (3) hence the gauge-dependent part of the OPEN-path ∫ phi is PATH-INDEPENDENT
#       (= Lambda(P) - Lambda(P0)): a change of local length UNIT, not a
#       history-dependent rescaling;
#   (4) physical length RATIOS at a point are invariant (ell1^2/ell2^2 unchanged
#       under ell^2 -> f ell^2); and the Weyl CONNECTION is invariant under the
#       combined transformation (eq 37), via linearity of the Weyl term + the
#       conformal shift = Weyl-term(1/2 d ln f) (Lean T2_10/weyl_connection).
# ============================================================================

# Output is captured via stdout redirection by the runner.

_lines = []
def out(s=''):
    print(s); _lines.append(str(s))

out('=' * 70)
out('Weyl objection residual — conformal gauge: no real length rescaling')
out('=' * 70)
out('')

var('t x y z')
coords = [t, x, y, z]

# ----------------------------------------------------------------------------
# (1) Field strength F_ij is gauge-invariant under phi -> phi + dLambda
# ----------------------------------------------------------------------------
out('-' * 70)
out('(1) F_ij = d_i phi_j - d_j phi_i is GAUGE-INVARIANT (phi -> phi + dLambda)')
out('-' * 70)
phi = [function('phi%d' % i)(t, x, y, z) for i in range(4)]
Lam = function('Lambda')(t, x, y, z)
phi_g = [phi[i] + diff(Lam, coords[i]) for i in range(4)]      # gauge-shifted
allinv = True
for i in range(4):
    for j in range(4):
        F  = diff(phi[j],   coords[i]) - diff(phi[i],   coords[j])
        Fg = diff(phi_g[j], coords[i]) - diff(phi_g[i], coords[j])
        if (Fg - F).simplify_full() != 0:
            allinv = False
            out('  F_%d%d changed!' % (i, j))
out('  F[phi+dLambda] - F[phi] = 0 for all i,j: %s' % allinv)
assert allinv
out('  => the observable field strength does NOT depend on the conformal gauge.')
out('')

# ----------------------------------------------------------------------------
# (2) Closed-loop holonomy ∮ phi is gauge-invariant (∮ dLambda = 0)
# ----------------------------------------------------------------------------
out('-' * 70)
out('(2) Holonomy ∮ phi is GAUGE-INVARIANT: ∮ dLambda = 0 on a closed loop')
out('-' * 70)
var('s')
# unit circle in (x,y); the gauge change of ∮phi is ∮ dLambda.
# Representative smooth scalar (the result holds for ANY Lambda: along the loop
# the integrand is d/ds[Lambda(loop(s))], a total derivative, so the closed-loop
# integral is Lambda(end) - Lambda(start) = 0).
Lc = x^2*y - y^3 + x
Lx = diff(Lc, x); Ly = diff(Lc, y)
sub_loop = {x: cos(s), y: sin(s)}
dLam_loop = (Lx.subs(sub_loop) * (-sin(s)) + Ly.subs(sub_loop) * cos(s))
loopL = integrate(dLam_loop.trig_reduce(), s, 0, 2*pi).simplify_full()
out('  ∮ dLambda (unit circle, representative Lambda) = %s' % loopL)
assert loopL == 0
out('  => the gauge change of the holonomy vanishes: ∮phi is gauge-invariant')
out('     (and equals ∬F by Stokes — fixed by the observable field).')
out('')

# ----------------------------------------------------------------------------
# (3) The gauge-dependent part of the OPEN-path ∫phi is PATH-INDEPENDENT
#     (= Lambda(P) - Lambda(P0)): a change of local UNIT, not history-dependent
# ----------------------------------------------------------------------------
out('-' * 70)
out('(3) Open-path gauge piece ∫ dLambda is PATH-INDEPENDENT (= unit change)')
out('-' * 70)
# same endpoints A=(0,0) -> B=(1,1), two different paths; ∫dLambda equal both.
Lxy = x^2*y - y^3 + x                # same representative scalar
gx = diff(Lxy, x); gy = diff(Lxy, y)
var('u')
# path 1: straight line x=u, y=u, u:0..1
p1 = {x: u, y: u}
i1 = integrate((gx.subs(p1)*1 + gy.subs(p1)*1), u, 0, 1)
# path 2: parabola x=u, y=u^2, u:0..1
p2 = {x: u, y: u^2}
i2 = integrate((gx.subs(p2)*1 + gy.subs(p2)*(2*u)), u, 0, 1)
diff_paths = (i1 - i2).simplify_full()
out('  ∫dLambda(straight) - ∫dLambda(parabola), same endpoints = %s' % diff_paths)
assert diff_paths == 0
endpts = (Lxy.subs({x: 1, y: 1}) - Lxy.subs({x: 0, y: 0}))
out('  ∫dLambda = Lambda(B) - Lambda(A) = %s (depends only on endpoints)' % endpts)
assert (i1 - endpts).simplify_full() == 0
out('  => the gauge rescaling factor exp(∫dLambda) depends only on endpoints:')
out('     it is a position-dependent change of the local length UNIT, NOT a')
out('     path/history-dependent (second-clock) effect.')
out('')

# ----------------------------------------------------------------------------
# (4a) Physical length RATIOS at a point are gauge-invariant (ell^2 -> f ell^2)
# ----------------------------------------------------------------------------
out('-' * 70)
out('(4a) Length RATIOS at a point invariant under ell^2 -> f·ell^2')
out('-' * 70)
var('f l1 l2', domain='positive')
ratio_before = l1^2 / l2^2
ratio_after = (f*l1^2) / (f*l2^2)
out('  (f l1^2)/(f l2^2) - l1^2/l2^2 = %s' % (ratio_after - ratio_before).simplify_full())
assert (ratio_after - ratio_before).simplify_full() == 0
out('  => only the overall unit changes; relative (measured) lengths are fixed.')
out('')

# ----------------------------------------------------------------------------
# (4b) Weyl CONNECTION invariant under (g->f g, phi->phi+1/2 d ln f), eq 37,
#      via linearity of the Weyl term (the conformal shift = Weyl-term(1/2 d ln f)
#      is the Lean result weyl_connection; linearity closes the cancellation).
# ----------------------------------------------------------------------------
out('-' * 70)
out('(4b) Weyl term is LINEAR in phi => connection invariant under eq 37')
out('-' * 70)
nn = 3
gm = matrix(SR, nn, nn, lambda a, b: var('g_%d%d' % (min(a,b), max(a,b))))  # symmetric g
gi = matrix(SR, nn, nn, lambda a, b: var('gi_%d%d' % (min(a,b), max(a,b))))  # symmetric ginv
phis = [var('ph_%d' % a) for a in range(nn)]
psis = [var('ps_%d' % a) for a in range(nn)]
def Wterm(p, k, i, j):
    return sum(gi[k, s]*(gm[s, i]*p[j] + gm[s, j]*p[i] - gm[i, j]*p[s]) for s in range(nn))
sumpsi = [phis[a] + psis[a] for a in range(nn)]
lin_ok = True
for k in range(nn):
    for i in range(nn):
        for j in range(nn):
            d = Wterm(sumpsi, k, i, j) - Wterm(phis, k, i, j) - Wterm(psis, k, i, j)
            if d.simplify_full() != 0:
                lin_ok = False
out('  W(phi+psi) - W(phi) - W(psi) = 0 for all k,i,j: %s' % lin_ok)
assert lin_ok
out('  => combined transform: Gamma_Weyl[f g, phi + (1/2)d ln f]')
out('       = (Gamma2[g] + W((1/2)d ln f))  -  (W(phi) + W((1/2)d ln f))')
out('       = Gamma2[g] - W(phi) = Gamma_Weyl[g, phi].   (eq 37)')
out('     [conformal shift Gamma2[f g] = Gamma2[g] + W((1/2)d ln f): Lean')
out('      T2_… weyl_connection; here we supply the linearity that closes it.]')
out('')

out('=' * 70)
out('SUMMARY — conformal gauge introduces NO real rescaling (eq 37):')
out('  (1) field strength F gauge-invariant.                                 [OK]')
out('  (2) closed-loop holonomy ∮phi gauge-invariant.                        [OK]')
out('  (3) open-path gauge piece path-independent (= local unit change only). [OK]')
out('  (4) length ratios invariant; Weyl connection invariant (eq 37).       [OK]')
out('  => the conformal/Weyl gauge changes only the local length UNIT; all')
out('     observables (F, holonomy, ratios, connection/parallel transport) are')
out('     gauge-invariant. No real, history-dependent second-clock effect from')
out('     the gauge freedom. This closes the last residual of the Weyl escape')
out('     (the residual check); together with T1_8_weyl_integrability.sage the')
out('     two-sector escape is fully verified.')
out('=' * 70)

# (result captured via stdout redirection by the runner)
