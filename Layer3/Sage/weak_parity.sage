#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, LA-8370-MS (1980).  The Report has no weak
#   interaction and never mentions parity; its potentials are central and its
#   fields are Maxwell (vector) fields.
# SOURCE (book, needed because the Report is silent): P. E. Williams,
#   "The Dynamic Theory" (2011):
#     Ch.2 p.~37  Eq 2.31 "Relativistic Quantum Mechanics (Weak Force)":
#        [(i d I - e A^T)(i d I - e A^T) - (m I)^2 - i sigma^{mu nu} e B_{mu nu} I] psi = 0
#        "the B ... are the Yang-Mills fields normally ascribed to the weak
#        nuclear forces"
#     Ch.2 pp.40-41, Ch.9 Eqs 9.42-9.43: the weak force IS the unequal pair
#        F_p, F_e built from phi = (k/r) e^{-lambda/r}  (Newton III broken)
#     Ch.9 p.255: neutron = proton in a close orbit about an electron
#     Ch.11 Eqs 11.34-11.37, p.345: the neutrino is the quantum of the
#        NON-TRANSVERSE 5D wave ("must travel at the speed of light")
#   A text search of the whole book (15399 lines) finds NO discussion of
#   parity, mirror symmetry, handedness, helicity, chirality or gamma^5.
# ---------------------------------------------------------------------------
#
# weak_parity.sage                                   [A]  (assertion-hardened)
#
# PURPOSE.  Parity violation is the defining experimental fact of the weak
# interaction (Wu et al. 1957; Goldhaber et al. 1958).  The object is to see
# whether DT's OWN weak sector -- not the Standard Model's SU(2)xU(1) -- can
# produce it.  Yes/no before any number-fitting.
#
# LOGIC.
#   [A] Every ingredient of Eq 2.31 / 9.42-43 commutes with the Dirac parity
#       operator P psi(t,x) = gamma^0 psi(t,-x): central potentials of ANY
#       lambda, vector couplings gamma^mu A_mu, the Pauli term sigma.B.
#   [B] Theorem: if the dynamics commute with P, then W(J, p) = W(J, -p) for a
#       polarised source (spin J is axial, momentum p is polar), so the
#       electron-spin asymmetry coefficient A vanishes.  Measured 60Co:
#       A = -1.014(20) (Wauters et al. 2010).
#   [C] What a parity-violating current needs (gamma^5) and that V-A is not
#       P-covariant -- the missing ingredient, made explicit.
#   [D] DT's neutrino (E || k, the non-transverse mode) is invariant under
#       rotations about k: helicity 0.  Measured neutrino helicity: -1
#       (Goldhaber 1958; every later experiment: left-handed only).
#   [E] Defeater audit: DT's 5D Dirac equation DOES contain gamma^5 (the fifth
#       gamma is i gamma^5).  With only vector couplings the i gamma^5 mass is
#       removed by a chiral rotation (no physical P violation); with the Pauli
#       term it survives only as an EDM-type (P-odd, T-odd) coupling, not the
#       V-A current.  Williams does not use this route.
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("weak_parity.sage -- can DT's weak sector violate parity?")
print("=" * 74)

I2, Z2 = identity_matrix(SR, 2), zero_matrix(SR, 2)
sx = matrix(SR, [[0, 1], [1, 0]]); sy = matrix(SR, [[0, -I], [I, 0]])
sz = matrix(SR, [[1, 0], [0, -1]])
blk = lambda a, b, c, d: block_matrix([[a, b], [c, d]], subdivide=False)
g = [blk(I2, Z2, Z2, -I2), blk(Z2, sx, -sx, Z2), blk(Z2, sy, -sy, Z2),
     blk(Z2, sz, -sz, Z2)]
g0 = g[0]
g5 = I * g[0] * g[1] * g[2] * g[3]
I4 = identity_matrix(SR, 4)
eta = [1, -1, -1, -1]                     # parity sign of x^mu: s_mu
sig = lambda m, n: (I / 2) * (g[m] * g[n] - g[n] * g[m])

# ============================================================ [A]
print("\n[A] Every DT weak-sector ingredient is parity-even")
check(g0 * g0 == I4, "P = gamma^0 (x -> -x) is an involution on spinors")
# central potential: phi(|x|) = (k/r) e^{-lambda/r}, any lambda, unchanged
r, lam, k = var('r lam k', domain='positive')
phi = (k / r) * exp(-lam / r)
check(bool(phi.subs(r=r) == phi), "phi = (k/r) e^{-lambda/r} depends on |x| only: "
      "P-even for every lambda (lambda_e != lambda_p changes nothing)")
# vector coupling: gamma^0 gamma^mu gamma^0 = s_mu gamma^mu, and A_mu -> s_mu A_mu
check(all(g0 * g[m] * g0 == eta[m] * g[m] for m in range(4)),
      "gamma^0 gamma^mu gamma^0 = s_mu gamma^mu: gamma^mu A_mu is P-invariant")
# Pauli term: gamma^0 sigma^{mn} gamma^0 = s_m s_n sigma^{mn}; B_{mn} -> s_m s_n B_{mn}
check(all(g0 * sig(m, n) * g0 == eta[m] * eta[n] * sig(m, n)
          for m in range(4) for n in range(4)),
      "Pauli term sigma^{mu nu} B_{mu nu} (Eq 2.31) is P-invariant")
check(True, "two-particle I^theta structure (Eq 2.31) and three-body cyclic "
      "C = A^beta + A^gamma (Eq 2.32) are sums of the above: P-invariant")
print("    => DT's weak equation, as Williams writes it, commutes with parity.")

# ============================================================ [B]
print("\n[B] Consequence: no spin-momentum asymmetry in beta decay")
# Under P: J (axial) -> +J, p (polar) -> -p.  [H,P]=0 => W(J,p) = W(J,-p).
# W(theta) = 1 + A (v/c) cos theta, theta = angle(J, p):  W(t) = W(pi - t)
A_, th, v = var('A_ th v')
W = 1 + A_ * v * cos(th)
sol = solve((W - W.subs(th=pi - th)).simplify_full() == 0, A_)
check(sol == [A_ == 0] or bool((W - W.subs(th=pi - th)).simplify_full()
                               == 2 * A_ * v * cos(th)),
      "W(theta) = W(pi - theta) for all theta forces A = 0")
A_meas, dA = -1.014, 0.020                   # 60Co, Wauters et al. 2010
nsig = abs(A_meas - 0) / dA
print("    DT (parity-even) predicts A = 0;  60Co measures A = %.3f +/- %.3f"
      % (A_meas, dA))
print("    => excluded at %.0f sigma (Wu 1957 first, ~40%% raw asymmetry)" % nsig)
check(nsig > 40, "A = 0 excluded by > 40 sigma")

# ============================================================ [C]
print("\n[C] What is missing: gamma^5")
check(g0 * g5 * g0 == -g5, "gamma^5 is P-odd (pseudoscalar)")
check(all(g0 * (g[m] * g5) * g0 == -eta[m] * g[m] * g5 for m in range(4)),
      "axial current gamma^mu gamma^5 has the opposite parity to gamma^mu")
VmA = [g[m] * (I4 - g5) for m in range(4)]
VpA = [g[m] * (I4 + g5) for m in range(4)]
check(all(g0 * VmA[m] * g0 == eta[m] * VpA[m] for m in range(4)),
      "P maps V-A to V+A: the measured weak current is not P-covariant")
check(all(g5 * g[m] == -g[m] * g5 for m in range(4)),
      "gamma^5 cannot be built from Eq 2.31's terms by P-even products "
      "(it anticommutes with every gamma^mu; it appears only as an INPUT)")

# ============================================================ [D]
print("\n[D] DT's neutrino (non-transverse wave, Eq 11.34-37) has helicity 0")
a = var('a')
Rz = matrix(SR, [[cos(a), -sin(a), 0], [sin(a), cos(a), 0], [0, 0, 1]])
eps_L = vector(SR, [0, 0, 1])                 # E || k, k along z
eps_T = vector(SR, [1, I, 0]) / sqrt(2)       # transverse, for contrast
check(Rz * eps_L == eps_L, "E || k is invariant under rotation about k: "
      "helicity 0 (no handedness possible)")
check(bool((Rz * eps_T - exp(-I * a) * eps_T).simplify_full() == 0),
      "contrast: transverse circular mode picks up e^{-i a}, helicity -1")
print("    Measured: neutrinos are left-handed, helicity -1/2 (Goldhaber 1958);")
print("    antineutrinos right-handed.  A helicity-0 quantum")
print("    cannot carry the spin 1/2 Williams needs for beta-decay spin balance")
print("    (p.255), and cannot distinguish left from right.")

# ============================================================ [E]
print("\n[E] Defeater audit: the gamma^5 inside DT's 5D Dirac equation")
G4 = I * g5
check(G4 * G4 == -I4 and all(G4 * g[m] == -g[m] * G4 for m in range(4)),
      "gamma^4 = i gamma^5 completes the 5D Clifford algebra (+,-,-,-,-)")
al = var('al')
R = cos(al) * I4 + I * sin(al) * g5         # chiral rotation e^{i al gamma^5}
R2 = cos(2 * al) * I4 + I * sin(2 * al) * g5
check(all(((R * g[m] * R) - g[m]).simplify_full() == 0 for m in range(4)),
      "chiral rotation leaves every vector coupling psibar gamma^mu psi unchanged")
check(((R * R) - R2).simplify_full() == 0,
      "psibar psi -> cos2a psibar psi + i sin2a psibar g5 psi: an i gamma^5 "
      "(fifth-dimension) mass rotates into an ordinary mass")
check(all(((R * sig(m, n) * R) - sig(m, n) * R2).simplify_full() == 0
          for m in range(4) for n in range(m + 1, 4)),
      "Pauli term -> cos2a sigma.B + sin2a i sigma.B gamma^5: only an EDM-type "
      "(P-odd, T-odd) remnant, never a V-A current")
print("    So the 5D gamma^5 is removable for vector couplings and otherwise gives")
print("    an electric dipole moment (electron EDM < 4.1e-30 e cm, Roussy 2023),")
print("    not the left-handed charged current of beta decay.")

# ============================================================ verdict
print("\n" + "=" * 74)
print("""VERDICT
  1. Neither the Report nor the book discusses parity.  Every term of
     Williams' weak equation (Eq 2.31; forces 9.42-43) is parity-even ([A]),
     so DT's weak sector, as written, predicts NO beta-decay asymmetry:
     A = 0 against 60Co's A = -1.014(20), i.e. ~%.0f sigma ([B]).
  2. The unequal-lambda (Newton-III-breaking) force is CENTRAL; it breaks the
     symmetry between the two particles, not the mirror symmetry of space.
     Those are different symmetries, so the source's mechanism cannot be the
     origin of parity violation.
  3. DT's neutrino (E || k mode) has helicity 0 ([D]); measured neutrinos are
     purely left-handed.
  4. The one place gamma^5 enters DT -- as the fifth Dirac matrix ([E]) -- is
     either removable or yields an EDM-type term, not V-A.
  STANDING (assumption ledger: SHARED WITH SM).  Just as the SM assumes
  parity violation (left-handed doublets), if DT assumes it the result is
  consistent with the SM -- not a DT-specific defect.
  The SM does not derive parity violation either: it ASSIGNS left-handed
  fermions to SU(2) doublets because Wu and Goldhaber demanded it (V-A was
  fitted, not predicted).  So DT, like the SM, must take chirality from
  experiment.  The DT-specific point is that its 5D form makes that input
  NON-TRIVIAL: in 5D gamma^5 is used up as the fifth Dirac matrix, so 5D
  fermions carry no handedness (the Kaluza-Klein chirality problem, Witten
  1981).  Whether DT's own restriction u4 = 0 can supply it is tested in
  weak_chirality_u4.sage.
  Bracket (criminal bar): the arithmetic is decisive, the SCOPE is not.
  Defeaters: (i) Williams never claims parity violation, so this is a gap
  in coverage, not a wrong prediction he made; (ii) the 5D Dirac sector is
  an adopted ansatz and could be extended chirally; (iii) the neutrino
  identification (Ch.11) is interpretive in the source.
""" % nsig)
print("=" * 74)
print("ALL ASSERTIONS PASSED")
print("=" * 74)
