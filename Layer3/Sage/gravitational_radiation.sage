#!/usr/bin/env sage
# --------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, LA-8370-MS (LASL, 1980), Eqs (39)/(67)
#   (the Weyl action whose gauge-fixed form is Einstein-Hilbert + Lambda,
#   weyl_ghilencea_comparison.sage) and Eqs (68)-(70) (the 5D gauge field
#   F_ij = d_i phi_j - d_j phi_i and its field equations); Book (2011) Eqs
#   (1.32)-(1.39) and Ch.11 "Five Dimensional Waves".  Williams' own remark
#   on gravity waves, "Against the Odds" (2009) p.170, as quoted in
#   nontransverse_pressure_rederivation.sage: "there is no gravity wave
#   independent of electromagnetic waves ... the gravity wave does not exist
#   by itself ... a good reason the gravity wave experiments have come up
#   negative."
# Data: Weisberg & Huang, ApJ 829 (2016) 55 (PSR B1913+16: m1 = 1.438,
#   m2 = 1.390 Msun, Pb = 27906.98 s, e = 0.6171334, intrinsic Pbdot =
#   -(2.398 +/- 0.004)e-12); LIGO/Virgo/Fermi/INTEGRAL, ApJL 848 (2017) L13
#   (GW170817: -3e-15 < (c_gw - c)/c < +7e-16); GW150914, PRL 116 (2016)
#   061102 (binary black hole, no electromagnetic counterpart).
# Full derivation, context and provenance:
#   Supplementary_Material.pdf -- §S3, §6.1 (layer boundary), Appendix A.
# --------------------------------------------------------------------------

# ============================================================================
# gravitational_radiation.sage                     [A]  (assertion-hardened)
#
# GRAVITATIONAL RADIATION IN DT: WHICH WAVES, HOW MANY, HOW FAST, HOW STRONG.
#
# A 5D metric theory radiates more than GR does.  This script counts DT's
# metric wave modes, applies the restriction that carried gaps #1 and #2
# (conservation of mass: g_44 and g_a4 held, kk_radion_g44 /
# test_body_motion), and checks what survives against the binary pulsar,
# GW170817 and the polarisation content.  It then asks what Williams' OWN
# gauge-field picture of gravity can radiate.
#
# WHAT THE SCRIPT ESTABLISHES
#   (A) Mode count (linearised Riemann on the Ricci-flat plane-wave space,
#       modulo gauge): 4D gives 2 (GR); free 5D gives 5 = 2 tensor + 2 vector
#       (g_a4) + 1 scalar (g_44).
#   (B) The 5D scalar mode is seen by a 4D detector as a BREATHING
#       polarisation (trace of the transverse tidal tensor != 0): the
#       signature of the dynamical radion, already excluded by Cassini.
#   (C) Under the restriction (h_A4 = 0) exactly 2 modes survive and they
#       are GR's plus and cross: traceless, transverse, no breathing.
#   (D) They travel on the null cone of the surface metric: the Ricci
#       condition forces omega^2 = c^2 k^2, so c_gw = c exactly, inside the
#       GW170817 window.
#   (E) Strength: the surface equations are 4D GR (gr_on_4d_surface, K = 0),
#       so the quadrupole formula holds.  Peters-Mathews for PSR B1913+16
#       gives Pbdot = -2.402e-12, against -(2.398 +/- 0.004)e-12 observed.
#   (F) Williams' GAUGE picture of gravity cannot carry these waves.  With
#       d_4 = 0 (no dependence on the mass coordinate) the 5D field equations
#       split exactly: the electromagnetic F_mu nu sector contains no phi_4,
#       and the gravitational V_mu = F_mu4 = d_mu phi_4 is the gradient of a
#       SCALAR.  A scalar has no tensor polarisation (and, like Nordstrom's
#       theory, gives no light bending), so the tensor waves come from the
#       metric sector only.  Coupling between the electric and the V sector
#       needs d_4 != 0; through DT's static gauge factor e^{-K_g gamma} it is
#       of order K_g gamma ~ 1e-31 for laboratory densities.
#   (G) Williams' remark that "the gravity wave does not exist by itself"
#       describes the gauge sector at most.  The recovered metric sector
#       carries free tensor waves with no electromagnetic partner, as
#       GW150914 (a binary black-hole merger with no counterpart) requires.
#       So the GR recovery SUPERSEDES that remark; it is not a prediction of
#       DT once the Weyl action is taken seriously.
#
# VERDICT (honest bracket).  Under the same single restriction that fixed
# light bending and free fall, DT's gravitational radiation is GR's: two
# tensor polarisations, speed c, quadrupole power, matching the binary pulsar
# to 0.2%.  The extra 5D modes (two vector, one breathing scalar) are exactly
# what the restriction removes, and their absence in LIGO-Virgo-KAGRA
# polarisation tests is a third, independent check of that restriction.
# Williams' "no gravity wave by itself" is superseded by the metric sector.
# What stays DT-specific is the gauge-sector non-transverse ("electrogravitic")
# mode, which needs dependence on the mass coordinate (d_4 != 0) and is a
# separate, electromagnetic-band question (nontransverse_pressure_
# rederivation.sage), not gravitational radiation in GR's sense.  Linear
# theory; the radiation reaction beyond quadrupole order is inherited from
# GR, not recomputed.
#
# RUN: /home/dsteele/miniconda3/envs/sage/bin/sage gravitational_radiation.sage
# ============================================================================

from sage.all import *


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


print("=" * 74)
print("gravitational_radiation.sage -- DT's wave modes vs GR, pulsars, LIGO")
print("=" * 74)


def setup(n, kup):
    eta = diagonal_matrix(QQ, [1] + [-1] * (n - 1))
    k = eta * vector(kup)                     # k_A (lower)
    pairs = [(i, j) for i in range(n) for j in range(i, n)]
    return eta, k, pairs


def eps_of(vec, n, pairs):
    E = matrix(QQ, n, n)
    for c, (i, j) in zip(vec, pairs):
        E[i, j] = c
        E[j, i] = c
    return E


def riemann(E, k, n):
    """linearised plane-wave Riemann (overall factor dropped)"""
    return {(A, B, C, D): (k[B] * k[C] * E[A, D] + k[A] * k[D] * E[B, C]
                           - k[A] * k[C] * E[B, D] - k[B] * k[D] * E[A, C])
            for A in range(n) for B in range(n) for C in range(n) for D in range(n)}


def ricci(Rm, eta, n):
    ei = eta.inverse()
    return matrix(QQ, n, n, lambda B, D: sum(ei[A, C] * Rm[(A, B, C, D)]
                                             for A in range(n) for C in range(n)))


def mode_count(n, kup, allowed=None):
    eta, k, pairs = setup(n, kup)
    basis = [vector(QQ, [1 if m == p else 0 for m in range(len(pairs))])
             for p in range(len(pairs))]
    if allowed is not None:
        basis = [b for b, pr in zip(basis, pairs) if allowed(pr)]
    # linear maps: eps -> Ricci, eps -> Riemann
    ric_cols, riem_cols = [], []
    for b in basis:
        E = eps_of(b, n, pairs)
        Rm = riemann(E, k, n)
        ric_cols.append(vector(QQ, ricci(Rm, eta, n).list()))
        riem_cols.append(vector(QQ, [Rm[key] for key in sorted(Rm)]))
    Mric = matrix(QQ, ric_cols).transpose()
    Mriem = matrix(QQ, riem_cols).transpose()
    ker = Mric.right_kernel().basis_matrix().transpose()   # Ricci-flat amplitudes
    if ker.ncols() == 0:
        return 0, eta, k, pairs, basis, ker
    return (Mriem * ker).rank(), eta, k, pairs, basis, ker


# ---------------------------------------------------------------------------
print("\n[A] Physical polarisations of a null plane wave along x1")
N4, *_ = mode_count(4, [1, 1, 0, 0])
N5, eta5, k5, pairs5, basis5, ker5 = mode_count(5, [1, 1, 0, 0, 0])
print("    4D: %d   free 5D: %d" % (N4, N5))
ok("A1: 4D metric gravity has 2 polarisations (GR)", N4 == 2)
ok("A2: free 5D metric gravity has 5 = n(n-3)/2 polarisations", N5 == 5 * (5 - 3) / 2)
Nnovec, *_ = mode_count(5, [1, 1, 0, 0, 0], allowed=lambda pr: pr[1] != 4 or pr == (4, 4))
ok("A3: dropping g_a4 (a < 4) removes exactly 2 modes (the KK vector)", N5 - Nnovec == 2)

# ---------------------------------------------------------------------------
print("\n[B] The g_44 mode is a breathing polarisation in 4D")
idx44 = pairs5.index((4, 4))
found = None
for col in ker5.columns():
    if col[idx44] != 0:
        found = col
        break
E = eps_of(found, 5, pairs5)
Rm = riemann(E, k5, 5)
tidal = matrix(QQ, 2, 2, lambda i, j: Rm[(0, i + 2, 0, j + 2)])   # transverse plane
print("    a Ricci-flat amplitude with eps_44 != 0; transverse tidal tensor:")
print("    %s" % str(tidal).replace("\n", "\n    "))
# remove any plus/cross admixture: breathing = trace part
ok("B1: the g_44 mode has a non-zero TRACE in the detector plane (breathing)",
   tidal.trace() != 0)

# ---------------------------------------------------------------------------
print("\n[C] Under the restriction h_A4 = 0: only GR's plus and cross")
Nr, eta_r, k_r, pairs_r, basis_r, ker_r = mode_count(5, [1, 1, 0, 0, 0],
                                                     allowed=lambda pr: pr[1] != 4)
print("    restricted 5D: %d" % Nr)
ok("C1: exactly 2 polarisations survive the restriction", Nr == 2)
traces = []
nonzero_tidal = []
for col in ker_r.columns():
    Ec = matrix(QQ, 5, 5)
    for c, b in zip(col, basis_r):
        Ec += c * eps_of(b, 5, pairs_r)
    Rc = riemann(Ec, k_r, 5)
    T = matrix(QQ, 2, 2, lambda i, j: Rc[(0, i + 2, 0, j + 2)])
    traces.append(T.trace())
    if T != 0:
        nonzero_tidal.append(T)
ok("C2: every surviving mode is traceless in the detector plane (no breathing)",
   all(tr == 0 for tr in traces))
ok("C3: the surviving tidal tensors span plus and cross (rank-2 traceless space)",
   matrix(QQ, [T.list() for T in nonzero_tidal]).rank() == 2)

# ---------------------------------------------------------------------------
print("\n[D] Propagation speed of the surviving modes")
w, kx = var('w kx')
eta4 = diagonal_matrix(SR, [1, -1, -1, -1, -1])
ksym = eta4 * vector(SR, [w, kx, 0, 0, 0])
Eplus = matrix(SR, 5, 5)
Eplus[2, 2] = 1
Eplus[3, 3] = -1
Rp = {(A, B, C, D): (ksym[B] * ksym[C] * Eplus[A, D] + ksym[A] * ksym[D] * Eplus[B, C]
                     - ksym[A] * ksym[C] * Eplus[B, D] - ksym[B] * ksym[D] * Eplus[A, C])
      for A in range(5) for B in range(5) for C in range(5) for D in range(5)}
Ric22 = sum(eta4.inverse()[A, C] * Rp[(A, 2, C, 2)] for A in range(5) for C in range(5))
print("    Ricci_22 of the plus mode = %s" % Ric22.factor())
ok("D1: Ricci-flatness requires w^2 = kx^2 (c_gw = c exactly)",
   bool((Ric22 / (w**2 - kx**2)).simplify_full().is_constant())
   and not bool(Ric22.is_zero()))
print("    so (c_gw - c)/c = 0, inside the GW170817 window [-3e-15, +7e-16]")

# ---------------------------------------------------------------------------
print("\n[E] Quadrupole power: PSR B1913+16 (Peters-Mathews)")
Tsun = 4.925490947e-6                     # G Msun / c^3  [s]
m1, m2, Pb, ecc = 1.438, 1.390, 27906.98161, 0.6171334
fe = (1 + 73 / 24.0 * ecc**2 + 37 / 96.0 * ecc**4) / (1 - ecc**2)**3.5
Pdot = -(192 * pi.n() / 5) * Tsun**(5 / 3.0) * (Pb / (2 * pi.n()))**(-5 / 3.0) \
    * fe * m1 * m2 * (m1 + m2)**(-1 / 3.0)
obs, sig = -2.398e-12, 0.004e-12
print("    GR/restricted-DT Pbdot = %.4e ; observed = (%.3f +/- %.3f)e-12"
      % (Pdot, obs * 1e12, sig * 1e12))
print("    observed / predicted = %.4f" % (obs / Pdot))
ok("E1: quadrupole prediction matches the binary pulsar within 3 sigma",
   abs(obs - Pdot) < 3 * sig)

# ---------------------------------------------------------------------------
print("\n[F] Williams' gauge picture of gravity with d_4 = 0")
X = list(var('x0 x1 x2 x3 x4'))
ph = [function('p%d' % i)(*X[:4]) for i in range(5)]     # no x4 dependence
F = matrix(SR, 5, 5, lambda i, j: diff(ph[j], X[i]) - diff(ph[i], X[j]))
ei = diagonal_matrix(SR, [1, -1, -1, -1, -1])
Fup = ei * F * ei
divF = [sum(diff(Fup[i, j], X[j]) for j in range(5)) for i in range(5)]
ok("F1: the four EM equations (i = 0..3) contain no phi_4",
   all(not divF[i].has(ph[4]) for i in range(4)))
ok("F2: the fifth equation contains only phi_4",
   all(not divF[4].has(ph[m]) for m in range(4)))
box = diff(ph[4], X[0], 2) - sum(diff(ph[4], X[i], 2) for i in range(1, 4))
ok("F3: the fifth equation is the scalar wave operator on phi_4",
   bool((divF[4] - box).simplify_full() == 0) or bool((divF[4] + box).simplify_full() == 0))
ok("F4: V_mu = F_mu4 = d_mu phi_4 -- gravity in the gauge picture is a gradient field",
   all(bool((F[m, 4] - diff(ph[4], X[m])).simplify_full() == 0) for m in range(4)))
H0, a0n, cl = 2.2e-18, 3.65e7, 2.998e8
Kg = H0 / (a0n * cl)                       # Book Eq 2.10
mix = Kg * 1.0e3
print("    EM-V coupling via e^{-K_g gamma}: K_g gamma ~ %.1e (gamma = 1e3 kg/m^3)" % mix)
ok("F5: the static EM-gravity mixing is < 1e-29 at laboratory density", mix < 1e-29)

# ---------------------------------------------------------------------------
print("\n[G] Williams: 'the gravity wave does not exist by itself'")
print("    The restricted metric sector [C] has free tensor waves sourced by")
print("    T_mu nu alone, with no electromagnetic partner required.")
print("    GW150914 (binary black holes, no counterpart) is such a wave, so the")
print("    remark is superseded by DT's own GR sector; read it as a statement")
print("    about the gauge sector at most.")

# ---------------------------------------------------------------------------
print("\n[VERDICT]")
print("  * Free 5D: 5 polarisations (2 tensor + 2 vector + 1 breathing scalar).")
print("  * Conservation of mass as a restriction leaves GR's plus and cross only,")
print("    at speed c, with quadrupole power matching PSR B1913+16 to 0.2%.")
print("  * The absent vector/breathing modes are a third test of the restriction.")
print("  * Williams' gauge gravity is a scalar gradient with d_4 = 0: no tensor")
print("    waves.  His 'no gravity wave by itself' is superseded by the metric sector.")
print("\n" + "=" * 74)
print("gravitational_radiation COMPLETE.")
print("=" * 74)
