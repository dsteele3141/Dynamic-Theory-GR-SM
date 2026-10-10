# 5D Minkowski Metric & Force Density (T2.1, Eqs 58-63)
#
# Source: P. E. Williams, LA-8370-MS (LASL, 1980), eqs (58)-(63).
# Full derivation, context and provenance:
#   Supplementary_Material.pdf — §S2, Appendix A [C1].

from sage.all import *
from sage.manifolds.all import *

# --- Setup SageManifolds ---
M = Manifold(5, 'M', structure='Lorentzian')
c_coords = M.chart('x0:q0 x1:x x2:y x3:z x4:q4')
q0, x, y, z, q4 = c_coords[:]


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


# Symbols for calculations
c_sym = var('c', domain='real')
assume(c_sym > 0)
a0 = var('a0', domain='real')
assume(a0 > 0)
gamma_f = function('gamma_f')(q0, x, y, z, q4)
beta = var('beta', domain='real')

# 1. 5D METRIC & NORMALIZATION
print("="*60)
print("COMPUTATION 1: 5D Metric and Velocity Normalization")
print("="*60)

# 5D Minkowski g_ij = diag(1, -1, -1, -1, -1)
# Note: g44 = -1 as requested.
g = M.metric()
g[0,0] = 1
g[1,1] = -1
g[2,2] = -1
g[3,3] = -1
g[4,4] = -1

print(f"\n5D Metric g_ij:\n{g[:]}")

u0, u1, u2, u3, u4 = var('u0 u1 u2 u3 u4', domain='real')
u_vec = [u0, u1, u2, u3, u4]
norm = sum(sum(g[i,j] * u_vec[i] * u_vec[j] for j in range(5)) for i in range(5))
print(f"\nVelocity Normalization g_ij u^i u^j = {norm.expr().simplify_full()}")
ok("(+,-,-,-,-) normalization: g_ij u^i u^j = u0^2 - u1^2 - u2^2 - u3^2 - u4^2.",
   (norm.expr() - (u0^2 - u1^2 - u2^2 - u3^2 - u4^2)).simplify_full().is_zero())
ok("Mass dimension is spacelike: g44 = -1.", (g[4,4].expr() + 1).is_zero())

# 2. EULER-LAGRANGE DERIVATION (Eq 62)
print("\n" + "="*60)
print("COMPUTATION 2: 5D Force Density (Eq 62)")
print("="*60)

print("\nEuler-Lagrange structure (NOT symbolically derived here):")
print("  LHS: gamma * f_i")
print("  RHS: d_gamma/dx^i - a0 * u4 * g_ij u^j")
print("  NOTE: structural description only -- no [PASS]. The Eq-62 Euler-Lagrange")
print("        derivation is verified rigorously in the Lean layer")
print("        (T2_4_ForceDensity) and the assertion-driven Layer-2 scripts.")

# 3. PAGE 86 AUDIT (Time Conversion)
print("\n" + "="*60)
print("COMPUTATION 3: Audit of Page 86 Time Conversion Factors")
print("="*60)

dq0_dt = c_sym * sqrt(1 - beta^2)
gamma_sym = var('gamma_sym', domain='real')

# Corrected factor for d/dt[v_alpha/sqrt(1-beta^2)]: gamma / (c^2 sqrt(1-beta^2)).
correct_factor = gamma_sym / (dq0_dt * c_sym)
# Manuscript (Page 86) printed factor: c^2 gamma / sqrt(1-beta^2).
manuscript_factor = c_sym^2 * gamma_sym / sqrt(1 - beta^2)

print(f"\nCorrect conversion factor: {correct_factor}")
print("  (Factor should be gamma / (c^2 * sqrt(1-beta^2)))")
print(f"Manuscript (Page 86) factor: {manuscript_factor}")
print("  CONCLUSION: Source Error (dimensional / factor-of-c mismatch).")
ok("Page-86 source error is real: corrected factor != manuscript factor.",
   not (correct_factor - manuscript_factor).simplify_full().is_zero())

# 4. SR REDUCTION
print("\n" + "="*60)
print("COMPUTATION 4: Special Relativity Limit (u4 = 0)")
print("="*60)

t = var('t')
v_alpha_func = function('v_alpha_func')(t)
gamma_rel = gamma_sym / sqrt(1 - beta^2)
sr_force = diff(gamma_rel * v_alpha_func, t)

print(f"\nStandard SR force density d/dt(gamma_rel v) = {sr_force}")
print("  NOTE: structural description only -- no [PASS]. The 5D force density is")
print("        not constructed here, so the u4->0 reduction is not symbolically")
print("        performed; the SR correspondence is documented in SI-6.2.3.")

print("\nAudit Complete: 11_5D_minkowski.sage")
print("="*60)
