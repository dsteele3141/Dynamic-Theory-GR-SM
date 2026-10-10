#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, LA-8370-MS (1980) Sec. II.F, pp.70-73:
#   hydrogen from London's quantum condition \oint phi_k dx^k = 2 pi i n with a
#   PURE Coulomb field phi_0 = e/r  ->  the Bohr radii.  The Report has no
#   non-singular potential, so it predicts exactly Bohr/Schroedinger hydrogen.
# SOURCE (book, needed because the Report stops at Bohr): P. E. Williams,
#   "The Dynamic Theory" (2011) Ch.9 "Nuclear Physics", pp.234-242:
#     Eq 9.2   phi = (k/r) e^{-lambda/r}, lambda different for each particle
#     Eq 9.6/7 F_1 uses lambda_2, F_2 uses lambda_1  (Newton III broken)
#     Eq 9.16  V   = -(mu k/r)[e^{-l1/r}/m2 + e^{-l2/r}/m1]
#              V_c = -(k/r)[e^{-l1/r} - e^{-l2/r}]          (energy OF the COM)
#     Eq 9.19  first-order Bohr shift; "a lambda_2 of 1 Fermi ... approximately
#              2x10^-5 different from Bohr's model"
#     p.241    radial Schroedinger eq with V + V_c, lambda_i << r:
#              s(s+1) = l(l+1) + lambda,
#              lambda = (2 mu Z e^2/4 pi eps0 hbar^2)[(2m1+m2) l1 - m1 l2]/(m1+m2)
#              delta  = -(l+1/2)[1 - sqrt(1 + lambda/(l+1/2)^2)],
#              E - E_c = E_n/(1 + 2 delta/n)
#   VALUES (Williams' own: text of Figs 2.3 and 9.3, and near Eq 11.170):
#     lambda_p ~ 1 fm  (from proton-proton scattering)
#     lambda_e ~ 1e-3 fm (assumed, "consistent with" e-e scattering)
# ---------------------------------------------------------------------------
#
# hydrogen_nonsingular_spectrum.sage                [A]  (assertion-hardened)
#
# PURPOSE. Does DT's non-singular potential give a hydrogen spectrum that
# matches experiment?  The object is to match MEASUREMENT with DT's own
# (non-standard) physics, not to reproduce QED.
#
# KEY POINT. For r >> lambda, e^{-lambda/r}/r = 1/r - lambda/r^2 + ...  The
# 1/r^2 tail is LONG-RANGE: it is not a contact term like the proton's finite
# size, it reaches every l.  It shifts the effective l, so it lifts the
# Coulomb degeneracy -- a Lamb-type shift by a non-QED route, with
#     Delta E(n,l) = Lt / (n^3 (2l+1)),   Lt = 4 Lambda Ry (mu/m_e)^2 / a0,
# where Lambda (a length) is the 1/r^2 coefficient divided by k.
#
# TWO READINGS (Williams never says which, same species as the phat photon):
#   REPLACES   -- DT has no QED radiative corrections; the Lambda term must
#                 supply the WHOLE Lamb shift.
#   SUPPLEMENTS -- QED is kept and the Lambda term is extra; the data then
#                 BOUND Lambda.
#
# What is imported (not re-derived here):  measured frequencies and the
# standard QED Lamb-shift values L(nlj) (to ~0.1 MHz, far finer than any
# residual found below), and the published statement that Dirac + QED + recoil
# + hyperfine reproduce the measured hydrogen lines to the kHz level.
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("hydrogen_nonsingular_spectrum.sage -- DT non-singular potential vs H")
print("=" * 74)

# ---------------------------------------------------------------- constants
Ry_Hz  = 3.2898419602508e15      # R_inf c        (CODATA 2018)
a0     = 5.29177210903e-11       # Bohr radius, m
fm     = 1e-15
me_mp  = 1 / 1836.15267343       # m_e / m_p
mmu_me = 206.7682830             # m_mu / m_e
h_eVs  = 4.135667696e-15         # Planck, eV s

LAM_E, LAM_P = 1e-3 * fm, 1.0 * fm     # Williams' values

# ============================================================ [A] Eq 9.16
print("\n[A] The 1/r^2 tail of Williams' V + V_c (Eqs 9.16)")
r, k, l1, l2, m1, m2 = var('r k l1 l2 m1 m2', domain='positive')
M  = m1 + m2
mu = m1 * m2 / M
x  = var('x')                                    # x = 1/r
V  = -(mu * k * x) * (exp(-l1 * x) / m2 + exp(-l2 * x) / m1)
Vc = -(k * x) * (exp(-l1 * x) - exp(-l2 * x))
def coeffs(expr):
    s = expr.taylor(x, 0, 2)
    return s.coefficient(x, 1).simplify_full(), s.coefficient(x, 2).simplify_full()
c1_V, c2_V   = coeffs(V)
c1_T, c2_T   = coeffs(V + Vc)
check(bool(c1_T == -k), "V + V_c -> -k/r at large r (Coulomb recovered)")
W_C2 = k * ((2 * m1 + m2) * l1 - m1 * l2) / M          # Williams p.241
check(bool((c2_T - W_C2).simplify_full() == 0),
      "1/r^2 coefficient of V+V_c = k[(2m1+m2)l1 - m1 l2]/M  (Williams' lambda)")
check(bool((c2_V - k * (m1 * l1 + m2 * l2) / M).simplify_full() == 0),
      "V alone: 1/r^2 coefficient = k(m1 l1 + m2 l2)/M")
# relabelling 1<->2 flips V_c, so the tail is NOT label-symmetric
swapped = (k * ((2 * m2 + m1) * l2 - m2 * l1) / M)
check(bool((swapped - W_C2).simplify_full() != 0),
      "tail is NOT symmetric under relabelling: V_c is odd in (1<->2)")

# ============================================================ [B] p.241
print("\n[B] Frobenius solution (p.241) re-derived")
ell, lam, n = var('ell lam n')
s = -QQ(1)/2 + sqrt((ell + QQ(1)/2)**2 + lam)        # root of s(s+1)=l(l+1)+lam
check(bool((s * (s + 1) - ell * (ell + 1) - lam).simplify_full() == 0),
      "s = -1/2 + sqrt((l+1/2)^2 + lam) solves s(s+1) = l(l+1) + lam")
delta_W = -(ell + QQ(1)/2) * (1 - sqrt(1 + lam / (ell + QQ(1)/2)**2))
check(bool((s - ell - delta_W).simplify_full() == 0) or
      all(abs(N((s - ell - delta_W).subs(ell=L_, lam=0.01))) < 1e-14 for L_ in [0,1,2]),
      "Williams' delta = s - l  (exact)")
check(bool(delta_W.taylor(lam, 0, 1) == lam / (2 * ell + 1)),
      "delta = lam/(2l+1) + O(lam^2)")
E_exact = n**2 / (n + var('d'))**2
E_W     = 1 / (1 + 2 * var('d') / n)
check(bool((E_exact - E_W).taylor(var('d'), 0, 1) == 0),
      "E_n n^2/(n+delta)^2 and Williams' E_n/(1+2 delta/n) agree to 1st order")
print("    => Delta E(n,l) = -2 delta E_n / n = Lt/(n^3 (2l+1)),"
      " Lt = 4 Lambda Ry_red / a_red")

def Lt_Hz(Lam, mred_over_me):
    # 4 Lambda Ry (mu/me) / (a0 me/mu)
    return 4 * Lam * Ry_Hz * mred_over_me**2 / a0
def dE(Lt, n_, l_):
    return Lt / (n_**3 * (2 * l_ + 1))

# ============================================================ [C] Eq 9.19
print("\n[C] Williams' Eq 9.19 number ('1 fm -> ~2e-5 from Bohr')")
mu_H = 1 / (1 + me_mp)                                # mu / m_e for H
corr_m2_e = mu_H * mu_H * (1.0 * fm) / a0             # body 2 = electron
corr_m2_p = mu_H * mu_H * me_mp * (1.0 * fm) / a0     # body 2 = proton
print("    body 2 = electron, lambda_2 = 1 fm : %.2e" % corr_m2_e)
print("    body 2 = proton,   lambda_2 = 1 fm : %.2e" % corr_m2_p)
check(1.5e-5 < corr_m2_e < 2.5e-5,
      "Williams' 2e-5 is reproduced only with lambda = 1 fm on the ELECTRON")
check(corr_m2_p < 1e-7, "with lambda = 1 fm on the proton it is ~1e-8, not 2e-5")
print("    SOURCE INCONSISTENCY: elsewhere lambda_e = 1e-3 fm and lambda_p = 1 fm."
      "\n    The 2e-5 illustration puts the proton's 1 fm on the electron.")

# ============================================================ [D] readings
print("\n[D] Lambda for hydrogen under each reading (Williams' lambda values)")
mp, me_ = 1.0, me_mp
Mh = mp + me_
readings = {
    "R1 e = body 1 (V+V_c)": (LAM_E * (2 * me_ + mp) - LAM_P * me_) / Mh,
    "R2 p = body 1 (V+V_c)": (LAM_P * (2 * mp + me_) - LAM_E * mp) / Mh,
    "R3 V only (no COM)   ": (LAM_E * me_ + LAM_P * mp) / Mh,
}
LAMB_MEAS = 1057.83e6                     # H 2S1/2-2P1/2, Hz (Bezginov 2019)
Lstar = LAMB_MEAS * 12 / Lt_Hz(1.0, mu_H)  # Lambda giving the WHOLE Lamb shift
print("    Lambda needed for the whole Lamb shift: %.4f fm" % (Lstar / fm))
for name, Lam in readings.items():
    lamb = dE(Lt_Hz(Lam, mu_H), 2, 0) - dE(Lt_Hz(Lam, mu_H), 2, 1)
    print("    %s Lambda = %+.3e fm   2S-2P = %+10.2f MHz  (%.3g x Lamb)"
          % (name, Lam / fm, lamb / 1e6, lamb / LAMB_MEAS))
check(readings["R1 e = body 1 (V+V_c)"] > 0 and
      readings["R1 e = body 1 (V+V_c)"] < 0.05 * Lstar,
      "R1: right sign, but < 5% of the Lamb shift")
check(readings["R2 p = body 1 (V+V_c)"] > 20 * Lstar,
      "R2: > 20x the Lamb shift")
check(readings["R3 V only (no COM)   "] > 10 * Lstar,
      "R3: > 10x the Lamb shift")

# ============================================================ [E] REPLACES
print("\n[E] REPLACES reading: can a free Lambda (no QED) fit hydrogen?")
# Truth (to ~kHz): E = Dirac(R) + L_QED.  DT model: Dirac(R') + Lt/(n^3(2l+1)).
# Dirac(R') - Dirac(R) = -eps/n^2 at the precision needed.  Residual per level:
#   res(n,l) = -eps/n^2 + Lt/(n^3(2l+1)) - L_QED(n,l)
L_QED = {(1,0): 8172.9e6, (2,0): 1045.0e6, (2,1): -12.8e6, (4,1): -1.6e6}
def res(eps, Lt, n_, l_):
    return -eps / n_**2 + dE(Lt, n_, l_) - L_QED[(n_, l_)]
eps_, Lt_ = var('eps_ Lt_')
sol = solve([res(eps_, Lt_, 2, 0) - res(eps_, Lt_, 2, 1) == 0,      # Lamb
             res(eps_, Lt_, 2, 0) - res(eps_, Lt_, 1, 0) == 0],     # 1S-2S
            eps_, Lt_, solution_dict=True)[0]
eps_f, Lt_f = float(sol[eps_]), float(sol[Lt_])
miss = abs(res(eps_f, Lt_f, 4, 1) - res(eps_f, Lt_f, 2, 0))          # 2S-4P
print("    fit to Lamb + 1S-2S: Lt = %.0f MHz, Rydberg shift eps = %.0f MHz"
      % (Lt_f / 1e6, eps_f / 1e6))
print("    predicted 2S-4P misses by %.0f MHz; measured to 2.3 kHz (Beyer 2017)"
      % (miss / 1e6))
check(miss > 1e8, "REPLACES: 2S-4P misses by > 100 MHz (~1e5 sigma)")
print("    Why: QED's Lamb shift lives almost entirely in S states and scales")
print("    ~1/n^3; the 1/r^2 tail shifts P (and D) states by 1/(2l+1) of the S")
print("    shift, so the l-pattern is wrong whatever Lambda is.")

# ============================================================ [F] SUPPLEMENTS
print("\n[F] SUPPLEMENTS reading: QED kept, Lambda bounded")
# Conservative: a free proton radius can absorb an S-only shift, so bound only
# the P-state part Lt/24, against ~10 kHz QED-vs-experiment agreement on 2S-2P.
TOL_H = 10e3
Lt_max_H = 24 * TOL_H
Lam_max_H = Lt_max_H / Lt_Hz(1.0, mu_H)
print("    |Lambda| < %.2e fm  (from |Lt|/24 < 10 kHz)" % (Lam_max_H / fm))
for name, Lam in readings.items():
    print("    %s exceeds the bound by %.0fx" % (name, abs(Lam) / Lam_max_H))
check(all(abs(L_) > 100 * Lam_max_H for L_ in readings.values()),
      "every reading with Williams' lambda values is excluded by > 100x")
lp_max_R1 = Lam_max_H * Mh / me_          # R1 with lambda_e -> 0
print("    R1 with lambda_e -> 0: lambda_p < %.1e fm (Williams: 1 fm)" % (lp_max_R1 / fm))

# mass-scaled reading lambda_i proportional to m_i (lambda_e = lambda_p me/mp)
Lam_ms_H = (LAM_P * me_ * (2 * me_ + mp) - LAM_P * me_) / Mh
print("\n    Mass-scaled lambda (lambda_i ~ m_i), R1: Lambda_H = %.2e fm, "
      "P-part %.1f kHz" % (Lam_ms_H / fm, Lt_Hz(Lam_ms_H, mu_H) / 24 / 1e3))
check(abs(Lam_ms_H) < Lam_max_H, "mass-scaled lambda SURVIVES ordinary hydrogen")

# ============================================================ [G] muonic H
print("\n[G] Muonic hydrogen (mu-p), Lamb shift 202.37 meV (Antognini 2013)")
mmu = mmu_me * me_mp
Mmu = mp + mmu
mu_mu = mmu_me * mp / Mmu                            # mu_red / m_e
TOL_MU = 0.12e-3                                     # eV: r_p free within PRad
for tag, lam_mu in [("lambda_mu = lambda_e (lepton-universal)", LAM_E),
                    ("lambda_mu = lambda_p m_mu/m_p (mass-scaled)", LAM_P * mmu)]:
    Lam = (lam_mu * (2 * mmu + mp) - LAM_P * mmu) / Mmu
    Lt_eV = Lt_Hz(Lam, mu_mu) * h_eVs
    print("    %s: Lambda = %+.3e fm, 2S-2P extra = %+.1f meV, P-part %.1f meV"
          % (tag, Lam / fm, Lt_eV / 12 * 1e3, Lt_eV / 24 * 1e3))
    check(abs(Lt_eV / 24) > 100 * TOL_MU,
          "muonic H excludes %s by > 100x" % tag.split(' (')[0])
Lam_max_mu = TOL_MU * 24 / (Lt_Hz(1.0, mu_mu) * h_eVs)
print("    muonic bound |Lambda_mu-p| < %.1e fm" % (Lam_max_mu / fm))

# ============================================================ verdict
print("\n" + "=" * 74)
print("""VERDICT
  1. Williams' algebra checks: the 1/r^2 tail of V + V_c, the Frobenius root,
     delta and the first-order energy are all re-derived exactly ([A],[B]).
     His 2e-5 illustration (Eq 9.19) uses 1 fm on the ELECTRON, inconsistent
     with his own lambda_e = 1e-3 fm ([C]).  The tail also depends on which
     particle is called '1' (V_c is odd under relabelling): an unresolved
     ambiguity in the source ([A],[D]).
  2. The non-singular potential gives a Lamb-TYPE shift by a non-QED route:
     long-range 1/r^2, right sign on reading R1.  That is a genuine,
     non-standard mechanism.
  3. REPLACES reading: it cannot be the Lamb shift.  Its l-pattern (P states
     shifted 1/3 as much as S) is wrong for any Lambda: fitting the Lamb shift
     and 1S-2S, 2S-4P misses by ~%.0f MHz against a 2.3 kHz measurement ([E]).
  4. SUPPLEMENTS reading: with Williams' lambda values every labelling is
     excluded by > 100x; hydrogen needs |Lambda| < %.1e fm.  Mass-scaled
     lambda (lambda_i ~ m_i) survives ordinary hydrogen but muonic hydrogen
     excludes it, as it does the lepton-universal choice ([F],[G]).
  STANDING: ordinary plus muonic hydrogen exclude every reading tried with
  Williams' lambda_p ~ 1 fm (from p-p scattering).  The 1/r^2 tail survives
  only if |Lambda| < ~1e-6 fm (e-p) and < ~1e-4 fm (mu-p).
  Bracket (criminal bar): a NEGATIVE CANDIDATE, not a falsification of DT.
  Defeaters: (i) Williams fixes lambda only roughly ('assumed'); (ii) the
  COM energy V_c and its labelling are not settled in the source; (iii) the
  e^{-lambda/r} form may hold only at nuclear range (a lambda that runs with
  separation would evade [E]-[G]).  What survives: the Report's own content
  -- London's condition -> Bohr/Schroedinger hydrogen -- is untouched.
""" % (miss / 1e6, Lam_max_H / fm))
print("=" * 74)
print("ALL ASSERTIONS PASSED")
print("=" * 74)
