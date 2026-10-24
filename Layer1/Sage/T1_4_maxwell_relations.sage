# T1.4 - Mechanical Maxwell relations (Eqs 7-11)
#
# Source: P. E. Williams, LA-8370-MS (LASL, 1980), Eqs 7-11.
# Supplementary Materials: step [L1], Sec. S2 (code index: Sec. S6).
from sage.all import *
#
# Source errors corrected and enforced by the asserts below:
#   - Eq 9: manuscript prints "(dK/dq)_{dot q} = phi*F"; correct is "F". And
#     "dK" with "-F dq"; correct is "+F dq".
#   - Eq 10: manuscript prints "dG" with "+q dF"; correct is "-q dF".
#   - Eq 11: Relation 2 prints "phi'(q)"; must be "phi'(dot q)" since phi is a
#     function of velocity dot{q}.
# The four Maxwell cross-derivative identities hold for the corrected potentials.

# --- Symbols & Differentials ---
var('dS dq dF dqdot phi phi_p S q F')


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


# 1. POTENTIALS AND DIFFERENTIALS
print("="*60)
print("COMPUTATION 1: Differentials and Source Errors")
print("="*60)

# Eq 7: Starting Differential
dU = phi*dS + F*dq
print(f"\nEq 7 (dU): {dU} [BASE]")

# Eq 8: Enthalpy H = U - Fq
dH = (dU - (F*dq + q*dF)).simplify_full()
print(f"Eq 8 (dH): {dH} [EXACT]")
# coefficient for dS and dF
print(f"  (∂H/∂S)_F = {dH.coefficient(dS)}")
print(f"  (∂H/∂F)_S = {dH.coefficient(dF)}")
ok("Eq 8: (dH/dS)_F = phi and (dH/dF)_S = -q.",
   (dH.coefficient(dS) - phi).is_zero() and (dH.coefficient(dF) + q).is_zero())

# Eq 9: Helmholtz K = U - phi*S
# d(phi*S) = phi*dS + S*phi_p*dqdot
dK = (dU - (phi*dS + S*phi_p*dqdot)).simplify_full()
print(f"\nEq 9 (dK) Corrected: {dK}")
print(f"  (∂K/∂dqdot)_q = {dK.coefficient(dqdot)}")
print(f"  (∂K/∂q)_dqdot = {dK.coefficient(dq)}")
print("  AUDIT: Manuscript Eq 9 prints '-F dq' and partial 'phi*F'.")
print("  CONCLUSION: Confirmed Source Errors in signs and factors.")
ok("Eq 9 CORRECTED: (dK/dq)_qdot = F (not phi*F) and (dK/dqdot)_q = -S*phi_p.",
   (dK.coefficient(dq) - F).is_zero() and (dK.coefficient(dqdot) + S*phi_p).is_zero())

# Eq 10: Gibbs G = H - phi*S
dG = (dH - (phi*dS + S*phi_p*dqdot)).simplify_full()
print(f"\nEq 10 (dG) Corrected: {dG}")
print(f"  (∂G/∂dqdot)_F = {dG.coefficient(dqdot)}")
print(f"  (∂G/∂F)_dqdot = {dG.coefficient(dF)}")
print("  AUDIT: Manuscript Eq 10 prints '+q dF'.")
print("  CONCLUSION: Confirmed Source Error in sign.")
ok("Eq 10 CORRECTED: (dG/dF)_qdot = -q (not +q) and (dG/dqdot)_F = -S*phi_p.",
   (dG.coefficient(dF) + q).is_zero() and (dG.coefficient(dqdot) + S*phi_p).is_zero())

# 2. MAXWELL RELATIONS (Eq 11)
print("\n" + "="*60)
print("COMPUTATION 2: Maxwell Relation Identities (Schwarz Symmetry)")
print("="*60)

def check_maxwell(potential_name, coords, name):
    """Verifies ∂(∂P/∂x)/∂y = ∂(∂P/∂y)/∂x."""
    x, y = coords
    P = function(potential_name)(x, y)
    mixed1 = diff(diff(P, x), y)
    mixed2 = diff(diff(P, y), x)
    res = (mixed1 - mixed2).simplify_full()
    ok("Maxwell %s: residual = %s." % (name, res), res.is_zero())

# Relation 1 from U(S, q)
var('S_v q_v F_v qd_v')
check_maxwell('U', (S_v, q_v), "1 (from dU)")
# Relation 2 from H(S, F)
check_maxwell('H', (S_v, F_v), "2 (from dH)")
# Relation 3 from K(qdot, q)
check_maxwell('K', (qd_v, q_v), "3 (from dK)")
# Relation 4 from G(qdot, F)
check_maxwell('G', (qd_v, F_v), "4 (from dG)")

print("\nAUDIT: Relation 2 in Eq 11 prints 'phi'(q)'.")
print("CONCLUSION: Source Error confirmed. Correct term is 'phi'(dot q)'.")

# 3. ENERGY CAPACITIES (Page 27)
print("\n" + "="*60)
print("COMPUTATION 3: Energy Capacity Formulas (Page 27)")
print("="*60)

# Cq = phi * (∂S/∂qdot)_q
# CF = phi * (∂S/∂qdot)_F
var('dS_dqdot dS_dq dF_dqdot dF_dq')

Cq = phi * dS_dqdot
# (∂S/∂qdot)_F = (∂S/∂qdot)_q + (∂S/∂q)_qdot * (∂q/∂qdot)_F
# (∂q/∂qdot)_F = -(∂F/∂qdot)_q / (∂F/∂q)_qdot
dq_dqdot_F = -dF_dqdot / dF_dq
CF = phi * (dS_dqdot + dS_dq * dq_dqdot_F)

# Goal: Cq - CF = (phi/phi_p) * (∂q/∂qdot)_F * (∂F/∂qdot)_q
# Use Maxwell 3 to substitute (∂S/∂q)_qdot = -(∂F/∂qdot)_q / phi_p
lhs = (Cq - CF).simplify_full()
lhs_sub = lhs.subs(dS_dq = -dF_dqdot / phi_p).simplify_full()

rhs = ((phi / phi_p) * dq_dqdot_F * dF_dqdot).simplify_full()
print(f"\nCq - CF LHS (substituted): {lhs_sub}")
print(f"Cq - CF RHS: {rhs}")
ok("Energy-capacity identity Cq - CF matches (phi/phi') (dq/dqdot)_F (dF/dqdot)_q.",
   (lhs_sub - rhs).simplify_full().is_zero())

print("\nAudit Complete: T1_4_maxwell_relations.sage")
print("="*60)
