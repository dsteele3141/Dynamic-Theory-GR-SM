# T1.10 - Equilibrium conditions (Eqs 12-13)
#
# Source: P. E. Williams, LA-8370-MS (LASL, 1980), Eqs 12-13.
# Supplementary Materials: step [L3], Sec. S2 (code index: Sec. S6).
from sage.all import *
#
# Verifies the equilibrium inequality dU - F dq - phi dS >= 0 (Eq 12), the
# free-energy maximum A_max and least-work principle (Eq 13), and the Helmholtz
# connection: equilibrium minimizes (K - Fq) at constant velocity (dot q),
# consistent with the corrected Eq 9 from T1.4.

# --- Symbols ---
var('phi S U q F delta_S delta_U delta_q Delta_E_rev')


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


# 1. CLAUSIUS INEQUALITY & EQUILIBRIUM (Eq 12)
print("="*60)
print("COMPUTATION 1: Equilibrium Inequality (Eq 12)")
print("="*60)

# Clausius: phi * dS >= Delta_E
# First Law: Delta_E = dU - F*dq
# For restoration (reverse of virtual displacement):
# phi * (-delta_S) >= (-delta_U) - F*(-delta_q)
restoration_lhs = phi * (-delta_S)
restoration_rhs = -delta_U + F*delta_q

print("\nRestoration Process (Spontaneous):")
print(f"  LHS: {restoration_lhs}")
print(f"  RHS: {restoration_rhs}")
print("  Condition: LHS >= RHS")

# Multiply by -1 to get the virtual displacement AWAY from equilibrium:
# -phi * (-delta_S) <= -(-delta_U) + (-F*delta_q) => phi * delta_S <= delta_U - F*delta_q
# delta_U - F*delta_q - phi * delta_S >= 0
eq12_lhs = delta_U - F*delta_q - phi*delta_S
print(f"\nRearranged Virtual Displacement form (Eq 12 LHS): {eq12_lhs}")
ok("Eq 12 LHS is the (sign-flipped) restoration margin (restoration_lhs - restoration_rhs).",
   (restoration_lhs - restoration_rhs - eq12_lhs).simplify_full().is_zero())
print("  NOTE: Williams uses strict '>' for true equilibrium (irreversibility).")

# 2. FREE ENERGY (Eq 13)
print("\n" + "="*60)
print("COMPUTATION 2: Free Energy Maximum (Eq 13)")
print("="*60)

# A = Delta_E_rev - Delta_U + F*Delta_q
# At maximum (reversible), Delta_E_rev = phi * Delta_S
A_max = phi*delta_S - delta_U + F*delta_q
print(f"\nFree Energy A_max (Eq 13): {A_max}")
ok("Eq 13: A_max = -(Eq 12 LHS).",
   (A_max + eq12_lhs).simplify_full().is_zero())

# Least work to DISPLACE = work released on the (reversible) restoration.
# Restoration quantities are (-delta_S, -delta_U, -delta_q).
# The correct relation is delta_A_min = A_restoration = Eq 12 LHS (a
# non-negative least-work cost).
A_restoration = phi*(-delta_S) - (-delta_U) + F*(-delta_q)
delta_A_min = A_restoration.simplify_full()
print(f"\nLeast work delta_A_min: {delta_A_min}")
ok("Least-work form delta_A_min equals Eq 12 LHS.",
   (delta_A_min - eq12_lhs).simplify_full().is_zero())

# 3. HELMHOLTZ CONNECTION (Eq 9 Context)
print("\n" + "="*60)
print("COMPUTATION 3: Connection to Helmholtz function K")
print("="*60)

# From T1.4 (Corrected Eq 9): K = U - phi*S
# Differential at constant velocity (dqdot=0):
# delta_K = delta_U - phi*delta_S
delta_K = delta_U - phi*delta_S

# Condition delta_A_min >= 0 becomes:
# (delta_U - phi*delta_S) - F*delta_q >= 0
# delta_K - F*delta_q >= 0
print(f"\nEquilibrium condition: {delta_K} - F*delta_q >= 0")
print("  => delta(K - Fq) >= 0 at constant dot{q}")
ok("Mechanical equilibrium minimizes (K - Fq): delta_K - F*delta_q = Eq 12 LHS.",
   (delta_K - F*delta_q - eq12_lhs).simplify_full().is_zero())

print("\nAudit Complete: T1_10_equilibrium.sage")
print("="*60)
