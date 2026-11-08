#!/usr/bin/env sage
# --------------------------------------------------------------------------
# tolman_ehrenfest_recovery - thermodynamic equilibrium in a gravitational field
# (case P != 0, F != 0): T sqrt(g00) = const
#
# Source (primary): P. E. Williams, LA-8370-MS (LASL, 1980), Eq (1) (dE = dU + P dV
#   - F dq: thermodynamic and mechanical work terms in ONE relation), p.9 ("A system
#   with a combination of thermodynamic and mechanical work terms will be considered
#   later"), p.76 ("p dV itself is just another F dq type term"). The Report does not
#   state the temperature of a system in a gravitational field; this script supplies
#   it. Inputs: the integrating factor of a purely thermodynamic system is the
#   absolute temperature (Report pp.23-24: "In the case of a purely thermodynamic
#   system the absolute quantity is the absolute zero temperature"), an
#   identification; positive heat capacities from the positive-definite stability
#   form (T1_5_stability_conditions); E = m c^2, so all energy, heat included, has
#   weight (the GR matter source, matter_source_tab); the static metric
#   g00 = e^(2 nu) and hydrostatic equilibrium from Div T = 0 on the recovered 4D
#   surface (gr_on_4d_surface, matter_source_tab).
#
# Supplementary Materials: step [G8], Sec. S3 (code index: Sec. S6).
# --------------------------------------------------------------------------
# Claim checked: in the case with both forces nonzero, equilibrium in a static
# gravitational field has T sqrt(g00) = const (Tolman 1930; Tolman & Ehrenfest 1930),
# and mu sqrt(g00) = const (Klein 1949). A recovery, not a prediction.
#   A. Weak field from the mixed first law: heat dQ lifted through height h does
#      work against F = -(dQ/c^2) g, so it arrives as dQ (1 - g h/c^2); no entropy
#      change on reversible transfer gives T_top = T_bot (1 - g h/c^2), i.e.
#      T (1 + Phi/c^2) = const to first order.
#   B. Exact, static metric g00 = e^(2 nu): the conserved (Killing) energy of a cell
#      is e^nu E_local; maximising S1 + S2 at fixed e^nu1 E1 + e^nu2 E2 gives
#      T1 e^nu1 = T2 e^nu2; the stationary point is a MAXIMUM when the heat
#      capacities are positive (the stability conditions).
#   C. Consistency with hydrostatic equilibrium dP/dr = -(rho + P) nu' (the
#      Tolman-Oppenheimer-Volkoff form of Div T = 0): with the Gibbs-Duhem relation
#      dP = s dT + n dmu and the Euler relation rho + P = T s + mu n, T e^nu = const
#      and mu e^nu = const give it identically.
#   D. Weak-field reduction e^nu = 1 + Phi/c^2 + O(c^-4); sizes: 1 km on Earth
#      Delta T/T = g h/c^2 ~ 1.1e-13 (unmeasurable); a neutron star, surface
#      redshift ~0.2-0.3, a temperature ratio of tens of percent (used in
#      neutron-star cooling models).
#   E. CONTROLS: weightless heat (uniform T) is NOT an entropy maximum (moving heat
#      down raises S) and contradicts hydrostatic equilibrium unless s = 0; the
#      opposite sign T e^(-nu) = const fails both tests.
#
# Inputs: CONTEXT P != 0, F != 0, static field, equilibrium; IDENTIFICATIONS the
#   integrating factor is the absolute temperature (Report pp.23-24), energy has
#   inertia and weight (E = m c^2). No new premise.

from sage.all import *


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


def zero(e):
    return bool(SR(e).simplify_full() == 0)


c, g, h, dQ, Tb = var('c g h dQ Tb', domain='positive')
eps = var('eps')                         # bookkeeping for 1/c^2

print("=" * 70)
print("tolman_ehrenfest_recovery.sage")
print("=" * 70)

# ---------------------------------------------------------------------------
print("\n(A) Weak field from the mixed first law dE = dU + P dV - F dq")
# heat dQ has mass dQ/c^2 (E = m c^2), so in field g the force on it is F = -(dQ/c^2) g;
# lifting it through h does work (dQ/c^2) g h, taken from the heat itself.
dQ_top = dQ - (dQ / c ** 2) * g * h
ok("heat dQ lifted through h arrives as dQ (1 - g h/c^2)", zero(dQ_top - dQ * (1 - g * h / c ** 2)))
Tt = var('Tt', domain='positive')
dS = -dQ / Tb + dQ_top / Tt
Tsol = solve(dS == 0, Tt)[0].rhs()
ok("reversible transfer dS = 0 gives T_top = T_bot (1 - g h/c^2)", zero(Tsol - Tb * (1 - g * h / c ** 2)))
Phi_b, Phi_t = 0, g * h
expr = Tsol * (1 + Phi_t / c ** 2) - Tb * (1 + Phi_b / c ** 2)
ok("so T (1 + Phi/c^2) is the same at both heights up to O(c^-4)",
   zero(taylor(expr.subs(c=1 / sqrt(eps)), eps, 0, 1)))

# ---------------------------------------------------------------------------
print("\n(B) Exact: entropy maximum at fixed Killing energy in g00 = e^(2 nu)")
nu1, nu2, E1, E2, lam, dE = var('nu1 nu2 E1 E2 lam dE', domain='real')
S1 = function('S1')
S2 = function('S2')
T1 = 1 / diff(S1(E1), E1)
T2 = 1 / diff(S2(E2), E2)
# move dE of Killing energy from cell 2 to cell 1: E1 += e^-nu1 dE, E2 -= e^-nu2 dE
Stot = S1(E1 + exp(-nu1) * dE) + S2(E2 - exp(-nu2) * dE)
first = diff(Stot, dE).subs(dE=0)
ok("dS/dE_Killing = e^-nu1/T1 - e^-nu2/T2: stationary iff T1 e^nu1 = T2 e^nu2",
   zero(first - (exp(-nu1) / T1 - exp(-nu2) / T2)))
second = diff(Stot, dE, 2).subs(dE=0)
C1, C2, t1, t2 = var('C1 C2 t1 t2', domain='positive')
# S'' = -1/(C T^2) with C = dE/dT > 0 (stability conditions)
sec = second.subs({diff(S1(E1), E1, 2): -1 / (C1 * t1 ** 2), diff(S2(E2), E2, 2): -1 / (C2 * t2 ** 2)})
ok("second variation = -e^(-2nu1)/(C1 T1^2) - e^(-2nu2)/(C2 T2^2) < 0: a MAXIMUM when C > 0",
   zero(sec + exp(-2 * nu1) / (C1 * t1 ** 2) + exp(-2 * nu2) / (C2 * t2 ** 2)))

# ---------------------------------------------------------------------------
print("\n(C) Consistency with hydrostatic equilibrium (TOV form of Div T = 0)")
r = var('r', domain='positive')
nu = function('nu')(r)
s_ = function('s_')(r)
n_ = function('n_')(r)
T0, mu0 = var('T0 mu0', domain='positive')
Tr = T0 * exp(-nu)                         # Tolman
mur = mu0 * exp(-nu)                       # Klein
dPdr = s_ * diff(Tr, r) + n_ * diff(mur, r)   # Gibbs-Duhem
rho_plus_P = Tr * s_ + mur * n_             # Euler relation
ok("T e^nu, mu e^nu const + Gibbs-Duhem give dP/dr = -(rho + P) nu' identically",
   zero(dPdr + rho_plus_P * diff(nu, r)))

# ---------------------------------------------------------------------------
print("\n(D) Weak-field reduction and sizes")
Phi = var('Phi')
ok("e^nu = sqrt(1 + 2 Phi/c^2) = 1 + Phi/c^2 + O(c^-4): (B) reduces to (A)",
   zero(taylor(sqrt(1 + 2 * Phi * eps), eps, 0, 1) - (1 + Phi * eps)))
cv, gv = 299792458, 9.80665
ratio = gv * 1000 / cv ** 2
ok("1 km on Earth: Delta T/T = g h/c^2 = %.2e (about 1e-13, unmeasurable)" % ratio,
   bool(1.0e-13 < ratio < 1.2e-13))
zs = 0.25
ok("neutron star, surface redshift z ~ 0.25: local T at the surface vs at infinity, 1/(1+z) = %.2f,"
   " a 20%% effect (as in cooling models)" % (1 / (1 + zs)), bool(abs(1 / (1 + zs) - 0.8) < 1e-12))

# ---------------------------------------------------------------------------
print("\n(E) CONTROLS")
Tu = var('Tu', domain='positive')
fu = first.subs({diff(S1(E1), E1): 1 / Tu, diff(S2(E2), E2): 1 / Tu})
ok("weightless heat (T1 = T2 = T) with nu1 < nu2: dS/dE = (e^-nu1 - e^-nu2)/T > 0, not a maximum"
   " (moving heat down raises S)",
   bool(fu.subs(nu1=0, nu2=QQ(1) / 10, Tu=1) > 0))
dPdr_u = s_ * 0 + n_ * diff(mur, r)
ok("uniform T contradicts hydrostatic equilibrium: residual = -T s nu' != 0 unless s = 0",
   zero(dPdr_u + (T0 * s_ + mur * n_) * diff(nu, r) - T0 * s_ * diff(nu, r))
   and not zero(dPdr_u + (T0 * s_ + mur * n_) * diff(nu, r)))
Tw = T0 * exp(nu)
muw = mu0 * exp(nu)
res_w = s_ * diff(Tw, r) + n_ * diff(muw, r) + (Tw * s_ + muw * n_) * diff(nu, r)
ok("wrong sign T e^(-nu) = const: residual 2 (rho + P) nu' != 0",
   zero(res_w - 2 * (Tw * s_ + muw * n_) * diff(nu, r)) and not zero(res_w))

print("\nAll checks passed.")
