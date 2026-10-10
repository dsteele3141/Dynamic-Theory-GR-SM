#!/usr/bin/env sage
# --------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, "The Dynamic Theory" (2011) Eqs (1.33)-
#   (1.39) (the 4D surface of the 5D space-time-mass manifold and its
#   Einstein equations under conservation of mass).  Exact plane-fronted
#   waves: Brinkmann, Math. Ann. 94 (1925) 119; Ehlers & Kundt, in
#   "Gravitation" (ed. Witten, 1962); Stephani et al., "Exact Solutions of
#   Einstein's Field Equations" (CUP, 2nd ed. 2003) Ch.24.
# Builds on gravitational_radiation.sage (linear mode count) and
#   kk_radion_g44.sage (the g_44 restriction).
# Full derivation, context and provenance:
#   Supplementary_Material.pdf -- §S3, Appendix A.
# --------------------------------------------------------------------------

# ============================================================================
# exact_ppwave_lift.sage                           [A]  (assertion-hardened)
#
# GRAVITATIONAL WAVES IN DT, EXACTLY -- NOT JUST TO FIRST ORDER.
#
# gravitational_radiation.sage counted DT's wave modes in linear theory.
# Here the same statement is made with no approximation, using the exact
# plane-fronted (pp) waves of GR:
#     ds^2 = H(u,x,y) du^2 + 2 du dv - dx^2 - dy^2       (4D)
# with arbitrary amplitude and waveform.  Its Ricci tensor is EXACTLY
# R_uu = +(1/2) lap_perp H (signature +,-,-,-), so it solves the full nonlinear vacuum
# equations whenever H is harmonic in the wave front.
#
# WHAT THE SCRIPT ESTABLISHES (all exact, no linearisation)
#   (A) The 4D pp-wave: the only Ricci component is R_uu = (1/2) lap_perp H.
#       The plus/cross wave H = A(u)(x^2 - y^2) + 2 B(u) x y, with ARBITRARY
#       waveforms A(u), B(u) of any amplitude, is an exact vacuum solution
#       with non-zero curvature (a real wave, not a coordinate artefact).
#   (B) Lift to DT: adding the mass coordinate with g_44 = -1 (the
#       restriction of gaps #1-#4) gives an exact 5D vacuum solution.  The
#       5D Ricci tensor equals the 4D one, with no mass-axis components.
#   (C) Test bodies: on the lifted wave Gamma^4 = 0 identically, so u^4 = 0
#       is preserved exactly -- conservation of mass survives the full
#       nonlinear wave.
#   (D) Exact speed c: the wave depends on u alone and d_v is a covariantly
#       constant NULL vector, so the wave travels exactly on the light cone,
#       at any amplitude.
#   (E) Exact polarisation content: the tidal tensor felt by test masses is
#       (1/2) d_i d_j H (up to sign convention), which the vacuum equation makes TRACELESS: plus and
#       cross only, never breathing, at any amplitude.
#   (F) What the restriction removes, exactly: if the wave may depend on the
#       mass coordinate, H = A(u)(x^2 + y^2 - 2 w^2) is ALSO an exact 5D
#       vacuum wave -- and on the surface w = 0 it is a pure BREATHING wave
#       (trace != 0).  This is the exact form of the linear scalar mode.
#       Without w-dependence it is not a solution.
#   (G) Exact nonlinear superposition: two different waveforms travelling
#       the same way add exactly (Ricci is linear in H for this class).
#
# VERDICT (honest bracket).  The wave sector of DT's GR recovery holds
# EXACTLY, not only to first order: every exact GR plane gravitational wave
# lifts to an exact 5D solution of DT's space-time-mass manifold under the
# conservation-of-mass restriction, travels at exactly c, carries only plus
# and cross, and keeps u^4 = 0 along every test-body path.  The exact
# breathing wave exists only when the field depends on the mass coordinate
# -- the reading the restriction (and Cassini, and LIGO polarisation tests)
# excludes.  Scope: pp-waves are the exact PLANE-wave class, not a binary
# merger; merger waveforms are inherited because the restricted surface
# equations are the full Einstein equations (gr_on_4d_surface), not
# re-solved here.
#
# RUN: /home/dsteele/miniconda3/envs/sage/bin/sage exact_ppwave_lift.sage
# ============================================================================

from sage.all import *


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


print("=" * 74)
print("exact_ppwave_lift.sage -- exact gravitational waves in DT's 5D manifold")
print("=" * 74)

u, v, x, y, w = var('u v x y w')


def geometry(g, X):
    n = len(X)
    gi = g.inverse()
    Gam = [[[sum(gi[a, d] * (diff(g[d, b], X[c]) + diff(g[d, c], X[b]) - diff(g[b, c], X[d]))
                 for d in range(n)) / 2 for c in range(n)] for b in range(n)] for a in range(n)]
    Ric = matrix(SR, n, n, lambda b, d: sum(diff(Gam[a][d][b], X[a]) - diff(Gam[a][a][b], X[d])
                                            + sum(Gam[a][a][f] * Gam[f][d][b] - Gam[a][d][f] * Gam[f][a][b]
                                                  for f in range(n)) for a in range(n)))
    return Gam, Ric.apply_map(lambda E: E.simplify_full())


def riemann_uiuj(g, X, Gam, i, j):
    """R^i_{u j u} -> tidal tensor component (index positions: a=i, b=u, c=j, d=u)"""
    n = len(X)
    a, b, c, d = i, 0, j, 0
    return (diff(Gam[a][d][b], X[c]) - diff(Gam[a][c][b], X[d])
            + sum(Gam[a][c][f] * Gam[f][d][b] - Gam[a][d][f] * Gam[f][c][b] for f in range(n))
            ).simplify_full()


def pp(H, X, extra=None):
    n = len(X)
    g = matrix(SR, n, n)
    g[0, 0] = H
    g[0, 1] = g[1, 0] = 1
    for k in range(2, n):
        g[k, k] = -1
    if extra:
        for (a, b), val in extra.items():
            g[a, b] = g[b, a] = val
    return g


# ---------------------------------------------------------------------------
print("\n[A] Exact 4D pp-wave")
X4 = [u, v, x, y]
Hg = function('Hg')(u, x, y)
G4, R4 = geometry(pp(Hg, X4), X4)
lap = diff(Hg, x, 2) + diff(Hg, y, 2)
ok("A1: the ONLY Ricci component is R_uu = (1/2) lap_perp H (exact, no linearisation)",
   bool((R4[0, 0] - lap / 2).simplify_full() == 0)
   and all(bool(R4[a, b] == 0) for a in range(4) for b in range(4) if (a, b) != (0, 0)))
A, B = function('A')(u), function('B')(u)
Hpc = A * (x**2 - y**2) + 2 * B * x * y
_, Rpc = geometry(pp(Hpc, X4), X4)
ok("A2: plus/cross with ARBITRARY waveforms A(u), B(u) is an exact vacuum solution",
   Rpc == 0)
Gpc, _ = geometry(pp(Hpc, X4), X4)
Rxx = riemann_uiuj(pp(Hpc, X4), X4, Gpc, 2, 2)
print("    curvature R^x_{u x u} = %s" % Rxx)
ok("A3: the curvature is non-zero (a physical wave, not a gauge artefact)",
   not bool(Rxx == 0))

# ---------------------------------------------------------------------------
print("\n[B] Lift to DT's 5D manifold with g_44 = -1 (the restriction)")
X5 = [u, v, x, y, w]
g5 = pp(Hpc, X5)
G5, R5 = geometry(g5, X5)
ok("B1: the lifted plus/cross wave is an exact 5D vacuum solution", R5 == 0)
G5g, R5g = geometry(pp(Hg, X5), X5)
ok("B2: for any H, the 5D Ricci equals the 4D Ricci (no mass-axis components)",
   all(bool((R5g[a, b] - (R4[a, b] if a < 4 and b < 4 else 0)).simplify_full() == 0)
       for a in range(5) for b in range(5)))

# ---------------------------------------------------------------------------
print("\n[C] Test bodies keep u^4 = 0 exactly")
ok("C1: Gamma^4_{ab} = 0 identically, so d u^4/ds = 0 on every geodesic",
   all(bool(G5g[4][b][c] == 0) for b in range(5) for c in range(5)))
print("    (u^4 = 0 initially stays 0: conservation of mass survives the wave)")

# ---------------------------------------------------------------------------
print("\n[D] Exact propagation at c")
kvec = vector(SR, [0, 1, 0, 0, 0])                      # d_v
norm = (kvec * g5 * kvec).simplify_full()
covd = [[(sum(G5[a][b][c] * kvec[c] for c in range(5))).simplify_full()
         for b in range(5)] for a in range(5)]          # nabla_b k^a (k constant)
ok("D1: d_v is null (g_vv = 0)", norm == 0)
ok("D2: d_v is covariantly constant -> the wave front moves exactly on the light cone",
   all(bool(covd[a][b] == 0) for a in range(5) for b in range(5)))

# ---------------------------------------------------------------------------
print("\n[E] Exact polarisation content (tidal tensor on test masses)")
tid = matrix(SR, 2, 2, lambda i, j: riemann_uiuj(g5, X5, G5, i + 2, j + 2))
print("    tidal tensor = %s" % str(tid).replace("\n", "\n                   "))
ok("E1: traceless at any amplitude -> plus and cross only, no breathing",
   bool(tid.trace().simplify_full() == 0))
tidg = matrix(SR, 2, 2, lambda i, j: riemann_uiuj(pp(Hg, X5), X5, G5g, i + 2, j + 2))
ok("E2: in general the tidal trace is +-(1/2) lap_perp H, zero exactly by the vacuum equation",
   bool((tidg.trace() + lap / 2).simplify_full() == 0)
   or bool((tidg.trace() - lap / 2).simplify_full() == 0))

# ---------------------------------------------------------------------------
print("\n[F] The exact breathing wave needs dependence on the mass coordinate")
Cb = function('Cb')(u)
Hbr = Cb * (x**2 + y**2 - 2 * w**2)
gbr = pp(Hbr, X5)
Gbr, Rbr = geometry(gbr, X5)
ok("F1: H = C(u)(x^2 + y^2 - 2w^2) is an exact 5D vacuum wave", Rbr == 0)
tbr = matrix(SR, 2, 2, lambda i, j: riemann_uiuj(gbr, X5, Gbr, i + 2, j + 2))
print("    tidal tensor on the surface = %s" % str(tbr).replace("\n", "\n                                  "))
ok("F2: on the surface it is a pure BREATHING wave (trace != 0)",
   not bool(tbr.trace().simplify_full() == 0) and bool(tbr[0, 1] == 0)
   and bool((tbr[0, 0] - tbr[1, 1]).simplify_full() == 0))
_, Rnow = geometry(pp(Cb * (x**2 + y**2), X5), X5)
ok("F3: without w-dependence the same breathing profile is NOT a vacuum solution",
   not bool(Rnow[0, 0] == 0))
ok("F4: and it carries a tidal force ALONG THE MASS AXIS (Gamma^4 != 0 off w = 0), changing u^4 of neighbours",
   any(not bool(Gbr[4][b][c] == 0) for b in range(5) for c in range(5)))

# ---------------------------------------------------------------------------
print("\n[G] Exact nonlinear superposition")
A2, B2 = function('A2')(u), function('B2')(u)
Hsum = Hpc + A2 * (x**2 - y**2) + 2 * B2 * x * y
_, Rsum = geometry(pp(Hsum, X5), X5)
ok("G1: two arbitrary co-moving waves superpose exactly", Rsum == 0)

# ---------------------------------------------------------------------------
print("\n[VERDICT]")
print("  * Every exact GR plane wave (any amplitude, any waveform) lifts to an")
print("    exact 5D solution of DT under the conservation-of-mass restriction.")
print("  * It moves at exactly c, carries plus/cross only, and keeps u^4 = 0.")
print("  * The exact breathing wave exists only with mass-coordinate dependence,")
print("    and a tidal force along the mass axis (u^4 != 0) -- the excluded reading.")
print("  * Gap #3 holds exactly in the wave sector, not just to first order.")
print("\n" + "=" * 74)
print("exact_ppwave_lift COMPLETE.")
print("=" * 74)
