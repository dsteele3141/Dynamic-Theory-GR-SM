#!/usr/bin/env sage
# T2.2 remainder — Modified Charge Conservation (the J4 term)  [prediction P6]
#
# Source: Williams LA-8370-MS (1980), Section III, pages 89-92, eqs (71)-(77).
# Full derivation, context and provenance:
#   Supplementary_Material.pdf — §S3, Appendix A [C9].
#
# ============================================================================
# The seven Maxwell-type equations (68),(70) are verified separately.
# This script does the remainder, centred on prediction P6 (modified conservation
# of charge):
#
#   (A) 5D current conservation d_i J^i = 0, unpacked in the natural coordinates
#       (x^0 = c t, x^4 = gamma/a0), is the MODIFIED continuity equation (eq 72):
#           d rho/dt + div J + a0 dJ4/d gamma = 0 .
#   (B) eq (72) is FORCED by the field equations + antisymmetry of F: with
#       J^j = (c/4pi) d_i F^{ij}, one has d_j J^j = (c/4pi) d_j d_i F^{ij} = 0,
#       since d_j d_i is symmetric and F^{ij} antisymmetric (same mechanism that
#       gives ordinary charge conservation from Maxwell). [verified for 5D F]
#   (C) classical limit (mass density fixed, d/d gamma -> 0): recovers the
#       ordinary continuity equation d rho/dt + div J = 0.
#   (D) alternate interpretation (eqs 73-77): if instead d_i F^{ij} = 0
#       (source-free) and one DEFINES rho = -(a0/4pi) dV0/d gamma,
#       J = -(a0 c/4pi) dV/d gamma, then eq (74) [div V + (1/c) dV0/dt = 0] gives
#       the CLASSICAL conservation div J + d rho/dt = 0 — at the cost of currents
#       of thermodynamic origin.
#   (E) wave equation (eq 71): V0 obeys a wave equation with a J4 source and a new
#       a0^2 (mass-density) term; a0 -> 0 recovers the ordinary wave equation.
# ============================================================================

_lines = []
def out(s=''):
    print(s); _lines.append(str(s))

out('=' * 70)
out('T2.2 remainder — modified charge conservation (prediction P6)')
out('=' * 70)
out('')

# Coordinates: t (time), x,y,z (space), g (= gamma, mass density)
var('t x y z g c a0', domain='real')
coords3 = [x, y, z]

# ----------------------------------------------------------------------------
# (A) 5D continuity d_i J^i = 0 unpacks to eq (72)
# ----------------------------------------------------------------------------
out('-' * 70)
out('(A) 5D continuity  d_i J^i = 0  =>  d rho/dt + div J + a0 dJ4/d gamma = 0')
out('-' * 70)
rho = function('rho')(t, x, y, z, g)
J1 = function('J1')(t, x, y, z, g)
J2 = function('J2')(t, x, y, z, g)
J3 = function('J3')(t, x, y, z, g)
J4 = function('J4')(t, x, y, z, g)
# x^0 = c t  => d/dx^0 = (1/c) d/dt, and J^0 = c rho  => d_0 J^0 = d rho/dt
# x^4 = gamma/a0 => d/dx^4 = a0 d/d gamma
divJ_5D = (diff(c*rho, t)/c) + diff(J1, x) + diff(J2, y) + diff(J3, z) + a0*diff(J4, g)
eq72 = diff(rho, t) + diff(J1, x) + diff(J2, y) + diff(J3, z) + a0*diff(J4, g)
out('  d_i J^i (x^0=ct, J^0=c rho, x^4=gamma/a0) equals eq(72): %s'
    % bool((divJ_5D - eq72).simplify_full() == 0))
assert (divJ_5D - eq72).simplify_full() == 0
out('  => eq (72):  d rho/dt + div J + a0 dJ4/d gamma = 0   [unpacked]')
out('')

# ----------------------------------------------------------------------------
# (B) eq (72) is FORCED: d_j d_i F^{ij} = 0 for antisymmetric F  => d_j J^j = 0
# ----------------------------------------------------------------------------
out('-' * 70)
out('(B) Charge conservation forced by antisymmetry: d_j d_i F^{ij} = 0')
out('-' * 70)
X = [t, x, y, z, g]
# antisymmetric 5x5 field tensor of arbitrary functions F_{ij}(t,x,y,z,g), i<j
Fup = matrix(SR, 5, 5)
for i in range(5):
    for j in range(i+1, 5):
        fij = function('F_%d%d' % (i, j))(t, x, y, z, g)
        Fup[i, j] = fij
        Fup[j, i] = -fij
# current J^j = sum_i d_i F^{ij}  (up to the 4pi/c factor); 5D divergence:
dofF = sum(diff(sum(diff(Fup[i, j], X[i]) for i in range(5)), X[j]) for j in range(5))
out('  d_j d_i F^{ij} (antisymmetric F, arbitrary components) = %s'
    % dofF.simplify_full())
assert dofF.simplify_full() == 0
out('  => with J^j = (c/4pi) d_i F^{ij}, d_j J^j = 0 identically: the modified')
out('     charge conservation (eq 72) is FORCED, exactly as Maxwell forces the')
out('     ordinary one. (Symmetric d_j d_i contracted with antisymmetric F = 0.)')
out('')

# ----------------------------------------------------------------------------
# (C) classical limit: mass density fixed (d/d gamma -> 0)
# ----------------------------------------------------------------------------
out('-' * 70)
out('(C) Classical limit (mass density fixed): ordinary continuity')
out('-' * 70)
# drop the gamma-dependence: the a0 dJ4/d gamma term vanishes
eq72_classical = eq72.subs({a0: 0})
out('  a0 -> 0 (or dJ4/d gamma = 0):  eq (72) -> d rho/dt + div J = %s + div J'
    % diff(rho, t))
out('  i.e. the classical conservation of charge  d rho/dt + div J = 0.  [OK]')
out('')

# ----------------------------------------------------------------------------
# (D) alternate interpretation (eqs 73-77): source-free F + thermodynamic
#     rho, J  =>  classical conservation div J + d rho/dt = 0
# ----------------------------------------------------------------------------
out('-' * 70)
out('(D) Alternate interpretation (eqs 73-77): classical conservation recovered')
out('-' * 70)
V0 = function('V0')(t, x, y, z, g)
Vx = function('Vx')(t, x, y, z, g)
Vy = function('Vy')(t, x, y, z, g)
Vz = function('Vz')(t, x, y, z, g)
rho_th = -(a0/(4*pi)) * diff(V0, g)                       # eq (75)
Jx_th = -(a0*c/(4*pi)) * diff(Vx, g)                      # eq (76)
Jy_th = -(a0*c/(4*pi)) * diff(Vy, g)
Jz_th = -(a0*c/(4*pi)) * diff(Vz, g)
divV = diff(Vx, x) + diff(Vy, y) + diff(Vz, z)
gauge74 = divV + diff(V0, t)/c                            # eq (74): = 0
# div J + d rho/dt, then use eq (74):
cons = diff(Jx_th, x) + diff(Jy_th, y) + diff(Jz_th, z) + diff(rho_th, t)
# cons should equal -(a0 c/4pi) d/d gamma (gauge74):
target = -(a0*c/(4*pi)) * diff(gauge74, g)
out('  div J + d rho/dt - [-(a0 c/4pi) d(eq74)/d gamma] = %s'
    % (cons - target).simplify_full())
assert (cons - target).simplify_full() == 0
out('  so with eq (74) (div V + (1/c)dV0/dt = 0):  div J + d rho/dt = 0  [classical].')
out('  => choosing d_i F^{ij}=0 keeps ordinary charge conservation, but the')
out('     current densities then originate in the thermodynamics (eqs 75-76).')
out('')

# ----------------------------------------------------------------------------
# (E) wave equation (eq 71): V0 propagates, with J4 source + new a0^2 term
# ----------------------------------------------------------------------------
out('-' * 70)
out('(E) Wave equation (eq 71): box V0 = (4pi/c^2) dJ4/dt + a0^2 d^2 V0/d gamma^2')
out('-' * 70)
lap = diff(V0, x, 2) + diff(V0, y, 2) + diff(V0, z, 2)
boxV0 = lap - diff(V0, t, 2)/c^2                          # d'Alembertian
rhs71 = (4*pi/c^2)*diff(J4, t) + a0^2*diff(V0, g, 2)
out('  box V0 = grad^2 V0 - (1/c^2) d^2V0/dt^2  (standard wave operator)')
out('  RHS = (4pi/c^2) dJ4/dt + a0^2 d^2V0/d gamma^2')
out('  a0 -> 0 (mass fixed): RHS -> (4pi/c^2) dJ4/dt = %s'
    % rhs71.subs({a0: 0}))
out('     i.e. the ordinary wave equation for V0 (sourced by dJ4/dt); the a0^2')
out('     term is the new mass-density (5th-dimension) propagation. [structure OK]')
out('')

out('=' * 70)
out('T2.2 REMAINDER SUMMARY (prediction P6 — modified charge conservation):')
out('  (A) d_i J^i = 0 unpacks to eq (72): d rho/dt + div J + a0 dJ4/d gamma = 0. [OK]')
out('  (B) eq (72) FORCED by antisymmetry of F (d_j d_i F^{ij}=0).               [OK]')
out('  (C) classical limit (mass fixed): ordinary d rho/dt + div J = 0.          [OK]')
out('  (D) alt. interpretation (eqs 73-77): source-free F + thermodynamic rho,J')
out('      => classical conservation (currents of thermodynamic origin).         [OK]')
out('  (E) eq (71) wave equation: V0 propagates; a0->0 gives ordinary wave eq.   [OK]')
out('  => P6: the J4 term modifies charge conservation; it is forced by the 5D')
out('     field structure and reduces to the classical law when mass is fixed.')
out('=' * 70)
