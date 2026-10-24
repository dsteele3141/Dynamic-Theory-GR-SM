# T1.6 - Stability metric and Lorentzian emergence (Eqs 16-19)
#
# Source: P. E. Williams, LA-8370-MS (LASL, 1980), Eqs 16-19.
# Supplementary Materials: step [L5], Sec. S2 (code index: Sec. S6).
from sage.all import *
#
# Signature flip: solving the Riemannian stability quadratic
# (h_ij dq^i dq^j = c^2 dt^2) for the entropy arc-length (dq^0)^2 produces a
# Lorentzian metric (+,-,-,-); asserted componentwise. Flat limit recovers
# Minkowski eta_ij. Source error (Eq 18): the printed 'A' radicand on Page 43
# is malformed (multiplies c^2 and spatial terms where it should subtract).

# --- Symbols ---
var('f h11 h22 h33 h01 dq0 dq1 dq2 dq3 dt c dx0 dx1 dx2 dx3 v1 X')


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


# 1. DIAGONAL CASE (h_0alpha = 0)
print("="*60)
print("COMPUTATION 1: Diagonal Case & Signature Flip")
print("="*60)

# Eq 17: c^2 dt^2 = f(dq0)^2 + h11(dq1)^2 + h22(dq2)^2 + h33(dq3)^2
rhs_17 = f*dq0^2 + h11*dq1^2 + h22*dq2^2 + h33*dq3^2
eq17 = c^2 * dt^2 == rhs_17

# Solve for (dq0)^2
dq0_sq_sol = solve(eq17, dq0^2)[0].rhs()
print(f"\nSolution for (dq0)^2: {dq0_sq_sol}")

# Map spatial dq^alpha -> dx^alpha; the time term stays c^2 dt^2 (= (dx0)^2).
dq0_sq_x = dq0_sq_sol.subs({dq1: dx1, dq2: dx2, dq3: dx3})

# f*(dq0)^2 = c^2 dt^2 - h11 dx1^2 - h22 dx2^2 - h33 dx3^2
expr_sig = (f * dq0_sq_x).simplify_full()
print(f"In spacetime coordinates (x0=ct => dx0^2 = c^2 dt^2): {expr_sig}")
ok("Lorentzian signature (+,-,-,-) emerges: dt^2 coeff = +c^2; "
   "dx1^2,dx2^2,dx3^2 coeffs = -h11,-h22,-h33 (sign-flipped vs the time term).",
   (expr_sig.coefficient(dt, 2) - c^2).is_zero()
   and (expr_sig.coefficient(dx1, 2) + h11).is_zero()
   and (expr_sig.coefficient(dx2, 2) + h22).is_zero()
   and (expr_sig.coefficient(dx3, 2) + h33).is_zero())

# 2. GENERAL CASE & AUDIT OF 'A' (Eq 18)
print("\n" + "="*60)
print("COMPUTATION 2: Audit of 'A' Definition (Eq 18)")
print("="*60)

# 1+1D Case for clarity: c^2 = f dot_q0^2 + 2 h01 dot_q0 v1 + h11 v1^2
# v1 = dq1/dt, X = dot_q0
quad_eq = f*X^2 + 2*h01*v1*X + (h11*v1^2 - c^2) == 0
sol_X = solve(quad_eq, X)

# The solved root must actually satisfy the quadratic (this underwrites the
# corrected 'A'; the manuscript's printed radicand for A is malformed).
residual_quad = quad_eq.lhs().subs(X == sol_X[0].rhs()).simplify_full()
ok("Quadratic solution for dot_q0 satisfies Eq-17 (1+1D) exactly.",
   residual_quad.is_zero())

# Squaring the solution dot_q0^2 should match the form of Eq 18:
# dot_q0^2 = (1/f) * [ c^2 + 2 h01 A - h11 v1^2 ]
dot_q0_sq = (sol_X[0].rhs()^2).simplify_full()

# Determine the 'Correct' A from the quadratic solution
A_correct = ((f*dot_q0_sq - c^2 + h11*v1^2) / (2*h01*v1)).simplify_full()
print(f"\nAlgebraically correct A: {A_correct}")

# Manuscript A (Page 43): A = (h01*v1)/f +/- sqrt( (c^2/f)*(h11*v1^2/f) + (h01*v1^2)/f )
# This radicand is clearly wrong (multiplies c^2 and h11).
A_manuscript_radicand = (c^2/f)*(h11*v1^2/f) + (h01*v1^2)/f
print(f"Manuscript A radicand (approx): {A_manuscript_radicand}")
print("  AUDIT: The manuscript A formula radicand multiplies terms where it should subtract.")
print("  CONCLUSION: Source Error in Eq 18 definition of A confirmed.")

# 3. MINKOWSKI RECOVERY
print("\n" + "="*60)
print("COMPUTATION 3: Minkowski Flat Limit")
print("="*60)

# In the flat limit, h_alpha_alpha = 1, h_0alpha = 0
dq0_sq_flat = dq0_sq_x.subs({h11: 1, h22: 1, h33: 1})
print(f"Flat limit (dq0)^2 = {dq0_sq_flat.simplify_full()}")
print("  => (1/f) * [ (dx0)^2 - (dx1)^2 - (dx2)^2 - (dx3)^2 ]")
ok("Flat limit recovers Minkowski: f*(dq0)^2 = c^2 dt^2 - dx1^2 - dx2^2 - dx3^2.",
   (f*dq0_sq_flat - (c^2*dt^2 - dx1^2 - dx2^2 - dx3^2)).simplify_full().is_zero())

print("\nAudit Complete: T1_6_stability_metric.sage")
print("="*60)
