#!/usr/bin/env sage
# T2.10 — London Quantisation (the L2->L3 gateway / quantum condition)
#   London condition ∮ phi_k dx^k = 2pi i n; Bohr radii r = n^2 hbar^2/(m e^2).
#
# Source: Williams LA-8370-MS (1980), Section II.F, pages 64-72, eqs (55)-(57).
# Full derivation, context and provenance:
#   Supplementary_Material.pdf — §S3, Appendix A [C1].

from sage.all import *

def ok(label, cond):
    print(f"  {label}: {'PASS' if cond else 'FAIL'}")
    assert cond

m, e, r, v = var('m e r v', domain='real')
assume(m > 0, e > 0, r > 0, v > 0)
hbar, h_planck, n = var('hbar h n', domain='real')
assume(hbar > 0, h_planck > 0, n > 0)
alpha, alpha1 = var('alpha alpha1')
c = var('c', domain='real')
assume(c > 0)

print("="*70)
print("15_london_quantisation_adjudicated.sage")
print("="*70)

print("\n(A) London Condition (Single-Valuedness)")
# Verify that exp(2*pi*i*n) is 1 for any integer n.
# Sage symbolic engine doesn't always automatically simplify exp(2*pi*I*n) == 1
# unless n is explicitly an integer and we use the right simplification.
n_val = 5
ok(f"exp(2*pi*i*n) = 1 for integer n (tested n={n_val})", bool(exp(2*pi*I*n_val) == 1))
# SCOPE (finding L3-01): this checks the FORWARD holonomy only -- the quantised
# form 2pi*i*n gives trivial holonomy. It does NOT prove the converse
# (exp(integral)=1 => integral=2pi*i*n) nor the DT bridge from length-transport
# to the loop integral. So single-valuedness is shown CONSISTENT with the London
# condition, not derived from first principles.
print("  Forward holonomy: the quantised form 2pi*i*n gives single-valuedness.")

print("\n(B) Bohr Radii Derivation")
v_sol = sqrt(e^2 / (m*r))
T = 2*pi*r / v_sol
phi_0 = alpha1 / r
integral = (phi_0 * c * T).canonicalize_radical()

forget(r > 0)
eq_london = (integral^2 == (2*pi*I*n)^2)
r_sol = solve(eq_london, r)[0].rhs()

alpha1_def = I * e^2 / (hbar * c)
r_final = r_sol.subs(alpha1=alpha1_def).simplify_full()

expected = n^2 * hbar^2 / (m * e^2)   # = n^2 a_B, a_B = Bohr radius (NOT the DT a0)
ok("Radius matches Bohr n^2*hbar^2/(m*e^2) = n^2 a_B (a_B != DT a0)",
   (r_final - expected).simplify_full() == 0)

print("\nSummary: forward holonomy + Bohr-radius algebra (a_B), conditional on the")
print("London orbit relation; the converse and the length-transport bridge are open.")
print("="*70)
