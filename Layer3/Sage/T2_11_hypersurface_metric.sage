# T2.11 — Hypersurface Metric Decomposition (Eqs 163-164)
#
# Source: Williams LA-8370-MS (LASL, 1980), eqs 163-164.
# Full derivation, context and provenance:
#   Supplementary_Material.pdf — §S3, Appendix A [C2].

from sage.all import *


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


# 1. SETUP MANIFOLDS
# Base 4D manifold M
M = Manifold(4, 'M')
U.<u0, u1, u2, u3> = M.chart()
u_coords = U[:]

# Ambient 5D manifold N
N = Manifold(5, 'N')
X.<x0, x1, x2, x3, x4> = N.chart()

# 2. DEFINE EMBEDDING Phi: M -> N
# x^a = u^a, x^4 = f(u)
f = function('f')(*u_coords)
phi = M.diff_map(N, (u0, u1, u2, u3, f), name='phi')

# 3. DEFINE 5D METRIC gh ON N
gh = N.tensor_field(0, 2, 'gh', sym=[(0,1)])
# We define components as functions of u because the pullback is on M
for i in range(5):
    for j in range(i, 5):
        # Using a trick to make gh components appear as functions of u
        gh[i,j] = function(f'gh{i}{j}')(*u_coords)

# 4. COMPUTE INDUCED METRIC A = phi^* gh
print("="*60)
print("COMPUTATION 1: Induced Metric Pullback (Eq 131)")
print("="*60)

A = phi.pullback(gh)

# Define the CORRECT formula for verification
f_derivs = [diff(f, u_coords[a]) for a in range(4)]
def correct_pullback_formula(a, b):
    # gh components are functions on N; .expr() pulls the symbolic expression
    # (in the shared u-coordinates) so it can combine with the M-chart pullback.
    return (gh[a,b].expr() + gh[a,4].expr()*f_derivs[b]
            + gh[b,4].expr()*f_derivs[a] + gh[4,4].expr()*f_derivs[a]*f_derivs[b])

print("\nVerifying A_ab matches the correct pullback formula (Eq 131 expansion):")
all_match = all((A[a,b].expr() - correct_pullback_formula(a,b)).simplify_full().is_zero()
                for a in range(4) for b in range(4))
ok("All A_ab components match the correct Eq-131 pullback.", all_match)

# 5. AUDIT MANUSCRIPT EQ (163)
print("\n" + "="*60)
print("COMPUTATION 2: Audit of Manuscript Eq (163) Literal Entries")
print("="*60)

# The manuscript's erroneous A_01 entry: gh01 + 2*gh04*f_1 + gh44*f_0*f_1
A01_manuscript = gh[0,1].expr() + 2*gh[0,4].expr()*f_derivs[1] + gh[4,4].expr()*f_derivs[0]*f_derivs[1]
diff_01 = (A[0,1].expr() - A01_manuscript).simplify_full()

print(f"\nManuscript Eq (163) entry A_01: {A01_manuscript}")
print(f"Pullback result A_01: {A[0,1]}")
print(f"  Residual (correct - manuscript) = {diff_01}")
print("  (Printed entry duplicates gh_04*f_1 instead of using gh_04*f_1 + gh_14*f_0)")
ok("Eq-163 source error is real: manuscript A_01 differs from the correct pullback.",
   not diff_01.is_zero())

# 6. DECOMPOSITION A_ab = gh_ab + h_ab (Eq 164)
print("\n" + "="*60)
print("COMPUTATION 3: Decomposition Analysis (Eq 164)")
print("="*60)

# Define h_ab from Eq 164 by computing components directly
# (A is the pullback on M; gh[a,b] are the user-supplied functions of u)
h_comp = [[(A[a,b].expr() - gh[a,b].expr()).simplify_full() for b in range(4)] for a in range(4)]

sym_ok = all((h_comp[a][b] - h_comp[b][a]).simplify_full().is_zero()
             for a in range(4) for b in range(4))
ok("h_ab is symmetric.", sym_ok)

# Verify h_ab structure
def expected_h(a, b):
    return (gh[a,4].expr()*f_derivs[b] + gh[b,4].expr()*f_derivs[a]
            + gh[4,4].expr()*f_derivs[a]*f_derivs[b])

all_h_pass = all((h_comp[a][b] - expected_h(a,b)).simplify_full().is_zero()
                 for a in range(4) for b in range(4))
ok("h_ab = A_ab - gh_ab matches the expected f-derivative structure (Eq 164).", all_h_pass)

print("\nAudit Complete: 05_hypersurface_metric.sage")
print("="*60)
