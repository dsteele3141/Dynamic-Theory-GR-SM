# velocity_factor_form - the FORM of the velocity factor phi(v) is derived
#
# Source: P. E. Williams, The Dynamic Theory (2011), Eqs 5.110-5.120 (pp. 136-141).
#   The 1980 Report (LA-8370-MS, p.64) takes phi = sqrt(1 - v^2/c^2) from
#   mass-spectrometer data; Williams 1997 (MPLA 12, 2725) derives only the
#   constancy of c and defers the form "to the equations of motion". The book
#   supplies that step; this script checks it.
# Supplementary Materials: step [L2], Sec. S2 (code index: Sec. S6).
#
# Chain checked:
#   1. Flat Lorentzian metric with arc length measured by c (the stability-form
#      identification, T1_6): d(sigma)/dt = sqrt(c^2 - v^2)        [Eq 5.111]
#   2. Equation of motion d/dsigma[m0 u] = F (geodesic of the entropy principle)
#      gives, in coordinate time, d/dt[m0 v/sqrt(1-b^2)] = Phi(v) Fhat(x) with
#      Phi = sqrt(1 - b^2)                                          [Eqs 5.112-5.114]
#   3. Exactness of dS = (dU - F dx)/phi with U = U(v) gives (Phi/phi)' = 0, so
#      phi = Phi = sqrt(1 - v^2/c^2) after phi(0) = 1               [Eqs 5.118-5.120]
#
# Premises (named, not derived here): the stability form as metric with arc
# length c dt (T1_6); U = U(v), i.e. all position dependence is carried by the
# work term F dx of the first law. CONTROL: if U may also depend on x, exactness
# no longer fixes phi (any phi is admissible), so that premise is load-bearing.
from sage.all import *

var('v x t c m0', domain='real')
assume(c > 0, m0 > 0, v > 0, v < c)
b = v / c
phi = function('phi')
Fhat = function('Fhat')


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


print("=" * 60)
print("STEP 1: d(sigma)/dt from the flat Lorentzian metric (Eq 5.111)")
print("=" * 60)
# d sigma^2 = c^2 dt^2 - dx^2  =>  (dsigma/dt)^2 = c^2 - v^2
dsig_dt = sqrt(c**2 - v**2)
ok("(dsigma/dt)^2 = c^2 - v^2 on the metric c^2 dt^2 - dx^2.",
   (dsig_dt**2 - (c**2 - v**2)).simplify_full() == 0)

print("\n" + "=" * 60)
print("STEP 2: coordinate-time force law (Eqs 5.112-5.114)")
print("=" * 60)
# u = dx/dsigma = v/(dsigma/dt); d/dsigma = (dt/dsigma) d/dt.
# Treat v as a function of t to differentiate.
vt = function('vt')(t)
u_t = vt / sqrt(c**2 - vt**2)
lhs = diff(m0 * c * u_t, t) / sqrt(c**2 - vt**2) * c   # c^2 d/dsigma[m0 u], sigma in length
# Williams: F^l = (1/sqrt(1-b^2)) d/dt[m0 v / sqrt(1-b^2)]
rhs = (1 / sqrt(1 - vt**2 / c**2)) * diff(m0 * vt / sqrt(1 - vt**2 / c**2), t)
ok("c^2 d/dsigma[m0 dx/dsigma] = (1/sqrt(1-b^2)) d/dt[m0 v/sqrt(1-b^2)] (Eq 5.112).",
   (lhs - rhs).simplify_full() == 0)
Phi = sqrt(1 - b**2)
print("  => d/dt[m0 v/sqrt(1-b^2)] = Phi(v) Fhat(x), Phi(v) = %s  (Eq 5.114)" % Phi)
ok("Newtonian limit: Phi(0) = 1.", Phi.subs(v=0) == 1)

print("\n" + "=" * 60)
print("STEP 3: exactness of dS fixes phi (Eqs 5.118-5.120)")
print("=" * 60)
# dS = U'(v)/phi dv - Phi Fhat/phi dx ; exactness: d/dv[-Phi Fhat/phi] = d/dx[U'/phi] = 0
U = function('U')
M = -Phi * Fhat(x) / phi(v)
N = diff(U(v), v) / phi(v)
cond = (diff(M, v) - diff(N, x))
ok("With U = U(v), exactness reduces to Fhat(x) d/dv[Phi/phi] = 0.",
   (cond + Fhat(x) * diff(Phi / phi(v), v)).simplify_full() == 0)
# Solve (Phi/phi)' = 0 => phi = Phi / k ; normalisation phi(0) = 1 => k = 1.
k = var('k')
phi_sol = Phi / k
ok("phi = Phi/k solves d/dv[Phi/phi] = 0.",
   diff(Phi / phi_sol, v).simplify_full() == 0)
k_val = solve(phi_sol.subs(v=0) == 1, k)[0].rhs()
phi_final = (Phi / k_val).simplify_full()
print("  => phi(v) = %s" % phi_final)
ok("phi(0) = 1 fixes k = 1, so phi(v) = sqrt(1 - v^2/c^2).",
   (phi_final - sqrt(1 - v**2 / c**2)).simplify_full() == 0)
ok("Consistency with T1_2: phi vanishes at the absolute velocity c.",
   limit(phi_final, v=c, dir='-') == 0)
# Energy integral (Eq 5.123): dU/dt = F v gives U = m0 c^2/sqrt(1-b^2) + const
Ugamma = m0 * c**2 / sqrt(1 - b**2)
p = m0 * v / sqrt(1 - b**2)
ok("Energy integral: dU/dv = v dp/dv for U = m0 c^2/sqrt(1-b^2) (Eq 5.123).",
   (diff(Ugamma, v) - v * diff(p, v)).simplify_full() == 0)

print("\n" + "=" * 60)
print("CONTROL: premise U = U(v) is load-bearing")
print("=" * 60)
# Allow U = h(v) + G(v) Vhat(x), Vhat' = Fhat. Exactness condition from T1_1-3:
#   C = phi F_v + (U_x - F) phi' = 0 with F = Phi Fhat.
# Choose G = Phi - Phi' phi/phi' : then C = 0 for ANY phi.
psi = sqrt(1 - b**2)**3      # an arbitrary wrong factor
G = Phi - diff(Phi, v) * psi / diff(psi, v)
C = psi * diff(Phi, v) * Fhat(x) + (G * Fhat(x) - Phi * Fhat(x)) * diff(psi, v)
ok("With x-dependent U, a different phi (here (1-b^2)^(3/2)) also passes exactness.",
   C.simplify_full() == 0)
ok("... so phi = sqrt(1-b^2) needs U = U(v) (the work term carries x).",
   (psi - Phi).simplify_full() != 0)

print("\nAudit Complete: velocity_factor_form.sage")
print("=" * 60)
