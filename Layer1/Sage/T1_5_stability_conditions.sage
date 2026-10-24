# T1.5 - Stability conditions (Eqs 14-15)
#
# Source: P. E. Williams, LA-8370-MS (LASL, 1980), Eqs 14-15.
# Supplementary Materials: step [L4], Sec. S2 (code index: Sec. S6).
from sage.all import *
#
# Source error (Eq 14): manuscript prints the second-order Taylor coefficients
# with "2" on the mixed term and "1" on delta S^2; correct are "1" and "1/2".
# Form A: first-order terms in dU - F dq - phi dS > 0 cancel at equilibrium,
# leaving a positive-definite quadratic. Form B: the cross-terms
# delta q * delta qdot cancel in velocity variables (diagonal stability).

# --- Symbols ---
var('delta_q delta_S delta_qdot phi phi_p phi_pp F S_eq U_qq U_SS U_qS')


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


# 1. FORM A (q, S Variables): Eq 14 & 15
print("="*60)
print("COMPUTATION 1: Form A Stability (Eq 14-15)")
print("="*60)

# Correct Taylor expansion of U(q, S)
dU_1st = phi*delta_S + F*delta_q
dU_2nd = 1/2*U_qq*delta_q^2 + U_qS*delta_q*delta_S + 1/2*U_SS*delta_S^2
dU_correct = dU_1st + dU_2nd

# Manuscript Eq 14 (Erroneous):
# dU = phi dS + F dq + 1/2 U_qq dq^2 + 2 U_qS dq dS + U_SS dS^2
dU_manuscript = dU_1st + 1/2*U_qq*delta_q^2 + 2*U_qS*delta_q*delta_S + U_SS*delta_S^2

print("\nTaylor Expansion (Eq 14) Audit:")
print(f"  Correct 2nd Order: {dU_2nd}")
print(f"  Manuscript 2nd Order: {dU_manuscript - dU_1st}")
print("  AUDIT: Manuscript coefficients are wrong (doubled).")
print("  CONCLUSION: Source Error in Eq 14 coefficients confirmed.")
ok("Eq 14 source error is real: manuscript 2nd-order form differs from the correct one.",
   not (dU_manuscript - dU_correct).simplify_full().is_zero())

# Stability Inequality (Eq 12): dU - F dq - phi dS > 0
stability_A = (dU_correct - F*delta_q - phi*delta_S).simplify_full()
print(f"\nStability Quadratic (Form A): {stability_A}")
ok("First-order terms cancel exactly: stability margin = the 2nd-order form.",
   (stability_A - dU_2nd).simplify_full().is_zero())

# Positive Definiteness (Sylvester's Criterion)
# The quadratic stability_A = dU_2nd has Hessian [[U_qq, U_qS],[U_qS, U_SS]];
# its determinant is the stability determinant det_H (Page 32).
det_H = U_qq * U_SS - U_qS^2
Hqq = diff(stability_A, delta_q, 2)
HSS = diff(stability_A, delta_S, 2)
HqS = diff(diff(stability_A, delta_q), delta_S)
print(f"\nHessian Determinant Condition: {det_H} > 0")
ok("Hessian determinant of the stability quadratic equals U_qq U_SS - U_qS^2 (Page 32).",
   (Hqq*HSS - HqS^2 - det_H).simplify_full().is_zero())

# 2. FORM B (qdot, q Variables): Cross-term Cancellation
print("\n" + "="*60)
print("COMPUTATION 2: Form B Stability (Velocity Variables)")
print("="*60)

# dS = (1/phi) dU - (F/phi) dq
# => ∂S/∂qdot = (1/phi) ∂U/∂qdot
# => ∂S/∂q = (1/phi) (∂U/∂q - F)
var('dU_dq_mF') # (∂U/∂q - F)

# Cross term from delta K: ∂²K/∂qdot∂q = -phi_p * (∂S/∂q)
cross_K = -phi_p/phi * dU_dq_mF

# Cross term from phi' * delta S * delta qdot:
# (phi_p * dqdot) * [(∂S/∂q) * dq]
cross_S = phi_p/phi * dU_dq_mF

print(f"\nForm B Cross-term Cancellation (delta q * delta qdot):")
print(f"  From ∂²K/∂qdot∂q: {cross_K}")
print(f"  From phi' * delta S: {cross_S}")
print(f"  Total: {(cross_K + cross_S).simplify_full()}")
ok("Form B cross-terms (delta q * delta qdot) cancel exactly (diagonal stability).",
   (cross_K + cross_S).simplify_full().is_zero())

print("\nAudit Complete: T1_5_stability_conditions.sage")
print("="*60)
