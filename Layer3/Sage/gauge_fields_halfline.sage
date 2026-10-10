#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, LA-8370-MS (1980) Sec. III.C, Eqs 105-112
#   (5D Dirac field, g = diag(1,-1,-1,-1,-1), x^4 = gamma/a0 >= 0: a HALF-LINE);
#   u4 = 0 as conservation of mass (Book Eq 1.36 / 10.6). Builds on
#   weak_chirality_u4.sage and noflux_boundary_derivation.sage (J^4(0) = 0,
#   one normalisable left-handed zero mode psi ~ exp(-M y)).
# Gauge fields on an interval / half-line: standard boundary-condition pairing
#   (e.g. Csaki, Grojean, Murayama, Pilo & Terning, PRD 69, 055006 (2004)).
# ---------------------------------------------------------------------------
#
# gauge_fields_halfline.sage                          [A] (assertion-hardened)
#
# QUESTION (the "A4 edge check", future-work 2a). DT adopts the SM's gauge
# fields. Do they sit consistently on DT's 5D half-line y = x^4 >= 0? What
# does their fifth component A_4 do at y = 0, and are the 4D zero modes the
# massless W, Z, photon and gluons? This is the one place where adopting the
# SM's dynamics could clash with DT. Odds before computing: consistent with
# no new prediction ~70%; a clash ~20%; a light scalar ~10%.
#
#   [A] boundary-condition pairing: gauge invariance at y = 0 allows only
#       (A_mu Neumann, A_4 Dirichlet) or (A_mu Dirichlet, A_4 Neumann)
#   [B] zero modes of each pairing; data select (N, D): no A_4 zero mode
#   [C] a BULK gauge field on the flat semi-infinite half-line is excluded:
#       non-normalisable zero mode, a continuum from m = 0, 5D Coulomb law
#   [D] the consistent reading: gauge fields independent of the mass axis
#       (d_4 A_mu = 0, alpha = alpha(x)), coupling to the y-integrated current;
#       charges universal, zero mode and edge condition untouched
#   [E] the abelian A_4 (Williams' V_4): J^4 = 0 on a chiral zero mode, so it
#       decouples from light matter; the non-abelian A_4 must vanish
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("gauge_fields_halfline.sage -- do the SM gauge fields fit DT's half-line?")
print("=" * 74)

y, L, m, M, K = var('y L m M K', domain='positive')

# ============================================================ [A]
print("\n[A] Boundary variation of -1/4 F_MN F^MN on y >= 0")
print("    delta S|_edge = - int d^4x  F^{mu 4} delta A_mu |_{y=0}")
print("    F_{mu 4} = d_mu A_4 - d_4 A_mu (+ commutator). Gauge parameter alpha:")
print("    delta A_mu = D_mu alpha, delta A_4 = D_4 alpha.")
print("    (N,D): d_4 A_mu(0) = 0, A_4(0) = 0  -- residual alpha(x, 0) free")
print("    (D,N): A_mu(0) = 0, d_4 A_4(0) = 0  -- alpha(x, 0) = 0 at the edge")
dA4, dmuA4, d4Amu, dAmu = var('dA4 dmuA4 d4Amu dAmu')      # values at y = 0
edge = lambda F, dA: F * dA                                # F_{mu4} delta A_mu
F_ND = dmuA4 - d4Amu                                       # abelian part suffices
check(bool(edge(F_ND, dAmu).subs({dmuA4: 0, d4Amu: 0}) == 0) and
      bool(edge(F_ND, dAmu).subs({dAmu: 0}) == 0),
      "the edge term vanishes for (N,D) [d4 A_mu = A_4 = 0 so F_{mu4} = 0] and (D,N) [delta A_mu = 0]")
check(not bool(edge(F_ND, dAmu).subs({d4Amu: 0}) == 0),
      "Neumann A_mu alone, with A_4 free at the edge, leaves the edge term: A_4 must be Dirichlet")

# ============================================================ [B]
print("\n[B] Kaluza-Klein modes on [0, L] (L -> oo later)")
# A_mu = sum f_n(y) A_mu^n(x), -f'' = m^2 f; A_4 modes g_n = f_n'/m_n (eaten)
f = function('f')(y)
fN = cos(m * y)         # Neumann at 0
fD = sin(m * y)         # Dirichlet at 0
check(bool(diff(fN, y).subs(y=0) == 0) and bool(fD.subs(y=0) == 0),
      "f_N = cos(m y) is Neumann, f_D = sin(m y) Dirichlet at y = 0")
check(bool(diff(fN, y, 2) + m**2 * fN == 0) and bool(diff(fD, y, 2) + m**2 * fD == 0),
      "both solve -f'' = m^2 f")
print("    (N,D): m = 0 gives f = 1 (a massless 4D vector) and g = f'/m = 0:")
print("           NO A_4 zero mode; the massive A_4 modes are eaten")
print("    (D,N): f(0) = 0 kills the m = 0 vector; A_4 has g = 1, a massless")
print("           4D scalar in the ADJOINT (gauge-Higgs)")
check(bool(fD.subs(m=0) == 0) and bool(fN.subs(m=0) == 1),
      "massless vector only for (N,D); massless adjoint scalar only for (D,N)")
print("    The W, Z, photon and gluons exist, and a massless colour-octet or")
print("    weak-triplet scalar does not  =>  (N,D) is SELECTED (no light scalar).")

# ============================================================ [C]
print("\n[C] Can the (N,D) zero mode be a BULK field on DT's half-line?")
norm_flat = integrate(1, y, 0, L)
check(bool(limit(norm_flat, L=oo) == oo),
      "flat half-line: int_0^L |f_0|^2 dy = L -> oo; g_4^2 = g_5^2 / L -> 0")
print("    The fermion zero mode is localised by its native mass:")
norm_psi = integrate(exp(-2 * M * y), y, 0, oo)
check(bool(norm_psi == 1 / (2 * M)), "int_0^oo exp(-2 M y) dy = 1/(2M): fermion normalisable")
print("    The gauge field has no such mass. Spacing of the bulk modes on [0, L]:")
spacing = pi / L
check(bool(limit(spacing, L=oo) == 0),
      "m_n = n pi / L: the spectrum becomes a continuum from m = 0")
print("    A continuum of light vectors gives a 5D force law, V(r) ~ 1/r^2, at")
print("    every distance; Coulomb's 1/r is measured to ~1e-16 in the exponent")
print("    (Williams, Faller & Hill 1971). A bulk SM gauge field on DT's FLAT")
print("    semi-infinite half-line is EXCLUDED.")
# A weighted measure would localise it: e.g. a weight exp(-K y)
norm_w = integrate(exp(-K * y), y, 0, oo)
check(bool(norm_w == 1 / K),
      "a measure weight exp(-K y) would give g_4^2 = g_5^2 K (finite)")
print("    DT's own gauge factor exp(-K_gamma gamma) has this form, but its")
print("    entry into the gauge kinetic term is NOT established: not adopted.")

# ============================================================ [D]
print("\n[D] The consistent reading: mass-axis-independent gauge fields")
print("    Take d_4 A_mu = 0, A_4 = 0, alpha = alpha(x) (u4 = 0 read as R1 for")
print("    the bosons). The fermion sees D_mu = d_mu - i q A_mu(x), D_4 = d_4.")
I4 = identity_matrix(SR, 4)
s1 = matrix(SR, [[0, 1], [1, 0]]); s2 = matrix(SR, [[0, -I], [I, 0]])
s3 = matrix(SR, [[1, 0], [0, -1]]); Z2 = zero_matrix(SR, 2); I2 = identity_matrix(SR, 2)
B = lambda a, b, c, d: block_matrix([[a, b], [c, d]], subdivide=False)
g0 = B(Z2, I2, I2, Z2)                                  # chiral basis
gs = [B(Z2, s, -s, Z2) for s in (s1, s2, s3)]
g5 = B(-I2, Z2, Z2, I2)
g4 = I * g5                                             # (g4)^2 = -1
G = [g0] + gs + [g4]
eta = diagonal_matrix(SR, [1, -1, -1, -1, -1])
check(all((G[a] * G[b] + G[b] * G[a]) == 2 * eta[a, b] * I4
          for a in range(5) for b in range(5)), "Clifford algebra in 5D, g = diag(1,-1,-1,-1,-1)")
print("    The edge condition J^4(0) = psibar g4 psi |_0 = 0 is gauge invariant")
print("    (psi -> e^{i q alpha(x)} psi leaves psibar g4 psi unchanged), and D_4 = d_4")
print("    leaves the zero mode exp(-M y) and its chirality untouched.")
fz = sqrt(2 * M) * exp(-M * y)
check(bool(integrate(fz**2, y, 0, oo) == 1),
      "normalised zero-mode profile f = sqrt(2M) exp(-M y)")
print("    Coupling: int dy  q psibar g^mu psi A_mu(x) = q (int |f|^2 dy) psi4bar g^mu psi4 A_mu")
check(bool(integrate(fz**2, y, 0, oo).subs(M=3) == integrate(fz**2, y, 0, oo).subs(M=7)),
      "the 4D charge equals q for EVERY M: charges universal across fermions")
print("    (A gauge field confined to the edge alone would give couplings")
print("    proportional to |f(0)|^2 = 2M, non-universal, and would not be gauge")
print("    covariant in the bulk: excluded by charge universality.)")

# ============================================================ [E]
print("\n[E] The fifth component A_4")
PL = (I4 - g5) / 2; PR = (I4 + g5) / 2
a = vector(SR, var('a1 a2 a3 a4'))
psiL = PL * a
J4L = (psiL.conjugate() * g0 * g4 * psiL).simplify_full()
check(bool(J4L == 0), "J^4 = psibar g4 psi = 0 identically on a chiral (left) spinor")
print("    So an abelian A_4 couples to J^4 only, and J^4 vanishes on every")
print("    light (chiral zero-mode) fermion: it DECOUPLES from light matter.")
print("    Abelian A_4 is gauge invariant (D_4 alpha = 0) and neutral: a sterile")
print("    scalar. This is where Williams' V_4 sits; it is not seen through EM.")
print("    Non-abelian A_4 transforms in the adjoint and couples to W and gluons")
print("    through D_mu A_4 = d_mu A_4 - i g [A_mu, A_4]: a light colour octet or")
print("    weak triplet would be produced at colliders and shift the rho")
print("    parameter, so A_4 = 0 (the (N,D) choice) is required for SU(2), SU(3).")
print("""
NET.
  * The adopted SM gauge fields fit DT's half-line under one reading: they are
    independent of the mass axis (d_4 A_mu = 0, gauge parameters alpha(x)),
    with A_4 = 0 for SU(2) and SU(3). This is the (Neumann, Dirichlet) pairing
    at y = 0, which data select. The W, Z, photon and gluons are the massless
    4D modes; there is NO light adjoint scalar.
  * The edge condition J^4(0) = 0, the zero mode exp(-M y) and its chirality
    are untouched, and every fermion's 4D charge is the same q whatever its M.
  * A BULK gauge field on the flat semi-infinite half-line is excluded (no
    normalisable zero mode, a continuum, a 1/r^2 force law). So, unlike the
    fermions, DT's gauge fields do not live in the mass dimension. This is a
    statement of the reading, adopted, not a derivation; a measure weight such
    as exp(-K_gamma gamma) could change it and is not established.
  * The abelian A_4 (Williams' V_4) decouples from light matter because J^4 = 0
    on chiral zero modes. It is sterile, not excluded, and gives no EM test.
  * Verdict: CONSISTENT. The owed check closes with no clash and no new light
    scalar. No prediction (the ~70% outcome).
""")
print("ALL ASSERTIONS PASSED")
