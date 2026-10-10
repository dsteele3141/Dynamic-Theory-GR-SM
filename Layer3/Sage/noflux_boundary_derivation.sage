#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE (primary): LA-8370-MS (1980) Eqs 105-112: (i gamma^j d_j + 1) psi = 0,
#   x^4 = gamma/a0 with gamma the mass density (so x^4 >= 0: a HALF-LINE).
#   u4 = dy4/dt = 0 is conservation of mass (Book Eq 1.36 / 10.6).
# Context: weak_chirality_u4.sage (readings R1-R3; R3 = J^4 = 0 at gamma = 0
#   gives one chiral zero mode), ledger Q3 (one 4x4 irrep) and Q4 (R3, so far
#   an ASSUMED reading).
# ---------------------------------------------------------------------------
#
# noflux_boundary_derivation.sage                      [A] (assertion-hardened)
#
# PURPOSE.  Ledger Q4 reads u4 = 0 as a no-flux boundary at gamma = 0.  Can
# that be DERIVED instead of assumed?
#
#   [A] the Dirac velocity along the mass axis is alpha^4 = gamma^0 gamma^4
#       (Heisenberg), and psi^dag alpha^4 psi = J^4: J^4 IS the Dirac form of u4
#   [B] 5D current conservation; on the half-line the 4D norm is conserved
#       iff J^4(0) = 0  (u4 = 0 at the edge = conservation of mass)
#   [C] unitarity: H is self-adjoint on y >= 0 only with a boundary condition
#       making B(phi,psi) = phi^dag alpha^4 psi vanish -- a U(2) family
#   [D] 4D Lorentz invariance of the boundary condition leaves exactly TWO
#       members: psi(0) left-handed or psi(0) right-handed (chiral)
#   [E] the chiral condition keeps the normalisable zero mode for one sign:
#       the remaining choice is the hand (Q3, the SM's own choice)
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("noflux_boundary_derivation.sage -- is Q4 (no-flux edge) derivable?")
print("=" * 74)

K = CyclotomicField(4); i_ = K.gen()
I2, Z2, I4 = identity_matrix(K, 2), zero_matrix(K, 2), identity_matrix(K, 4)
s = [matrix(K, [[0, 1], [1, 0]]), matrix(K, [[0, -i_], [i_, 0]]),
     matrix(K, [[1, 0], [0, -1]])]
B = lambda a, b, c, d: block_matrix([[a, b], [c, d]], subdivide=False)
g = [B(I2, Z2, Z2, -I2)] + [B(Z2, si, -si, Z2) for si in s]
g5 = i_ * g[0] * g[1] * g[2] * g[3]
kappa = 1
G = g + [i_ * kappa * g5]                          # one 4x4 irrep (Q3)
eta = [1, -1, -1, -1, -1]
check(all(G[j] * G[k] + G[k] * G[j] == (2 * eta[j] if j == k else 0) * I4
          for j in range(5) for k in range(5)), "5D Clifford algebra (Eq 112)")
H_ = lambda X: X.conjugate_transpose()

# ============================================================ [A]
print("\n[A] What is u4 for a Dirac particle?")
# (i g^j d_j + M) psi = 0  ->  i d_t psi = H psi,
# H = -i alpha^k d_k + beta' ,  alpha^k = g^0 g^k (k = 1..4), beta' = -M g^0
al = [g[0] * G[k] for k in range(1, 5)]
check(all(H_(a) == a for a in al) and all(a * a == I4 for a in al),
      "alpha^1..alpha^4 Hermitian with square 1: H is formally Hermitian")
y = var('y'); f = function('f')(y)
# i[H, y] acting on f: only the -i alpha^4 d_y term fails to commute with y
comm = -I * (diff(y * f, y) - y * diff(f, y))             # [-i d_y, y] f
check(bool((comm + I * f).simplify_full() == 0),
      "[ -i d_y , y ] = -i  =>  dy/dt = i[H, y] = alpha^4  (Heisenberg)")
check(g[0] * G[4] == al[3],
      "so psi^dag (dy/dt) psi = psi^dag g^0 g^4 psi = psibar g^4 psi = J^4")
print("    u4 = dy4/dt for a Dirac particle is the operator alpha^4; its")
print("    density is J^4.  J^4 = 0 IS u4 = 0 on the Dirac side (derived).")

# ============================================================ [B]
print("\n[B] Current conservation on the half-line y >= 0")
check(all(g[0] * H_(G[j]) * g[0] == G[j] for j in range(5)),
      "g^0 (g^j)^dag g^0 = g^j for all five: the adjoint equation follows and "
      "d_j (psibar g^j psi) = 0 for any solution")
print("    integrate d_0 J^0 + d_k J^k + d_y J^4 = 0 over 3-space and y in [0, oo):")
print("    dN/dt = J^4(y=0) - J^4(oo) = J^4(0)   (normalisable: J^4(oo) = 0)")
print("    => the 4D particle number is conserved iff J^4(0) = 0, i.e. iff")
print("       u4 = 0 AT THE EDGE.  Conservation of mass, read on the half-line.")

# ============================================================ [C]
print("\n[C] Unitarity: self-adjoint boundary conditions")
# <phi, H psi> - <H phi, psi> = -i int d_y (phi^dag alpha^4 psi) dy
#                             = i phi(0)^dag alpha^4 psi(0)
a4 = al[3]
ev = a4.charpoly().roots(multiplicities=True)
check(sorted(ev) == [(-1, 2), (1, 2)],
      "boundary form B = alpha^4 has signature (2,2): self-adjoint BCs exist")
print("    self-adjoint <=> psi(0) confined to a 2-dim subspace W with B(W,W) = 0;")
print("    these form a U(2) family (graphs of unitaries between the +-1 eigenspaces).")
print("    Every member has J^4(0) = psi(0)^dag alpha^4 psi(0) = 0.")
print("    => J^4(0) = 0 is FORCED by unitarity on DT's half-line, independent")
print("       of [B]: conservation of mass and unitarity give the same condition.")

# ============================================================ [D]
print("\n[D] Which members respect 4D Lorentz invariance?")
Sig = [(g[m] * g[n] - g[n] * g[m]) / 4 for m in range(4) for n in range(m + 1, 4)]
In = I4
Mk = block_matrix([[X.transpose().tensor_product(In) - In.tensor_product(X)] for X in Sig],
                  subdivide=False)
comm_basis = [matrix(K, 4, 4, list(v)).transpose() for v in Mk.right_kernel().basis()]
print("    commutant of the 4D Lorentz generators: dimension %d" % len(comm_basis))
check(len(comm_basis) == 2, "commutant = span{1, gamma^5}: Dirac = (1/2,0) + (0,1/2), inequivalent")
PL, PR = (I4 - g5) / 2, (I4 + g5) / 2
check(all(X * Y == Y * X for X in comm_basis for Y in (PL, PR)) and
      all((PL * X * PR).is_zero() for X in comm_basis),
      "so the only invariant 2-dim subspaces are L = range(P_L) and R = range(P_R)")
check((H_(PL) * a4 * PL).is_zero() and (H_(PR) * a4 * PR).is_zero(),
      "both are isotropic for B (J^4 = 0 on each): both are self-adjoint BCs")
# a non-chiral member exists (so Lorentz invariance is doing the work)
eP = (a4 - I4).right_kernel().basis(); eM = (a4 + I4).right_kernel().basis()
W = matrix(K, [eP[0] + eM[0], eP[1] + eM[1]]).transpose()
check((H_(W) * a4 * W).is_zero() and not (PR * W).is_zero() and not (PL * W).is_zero(),
      "e.g. W = graph(identity) is self-adjoint but NOT chiral -- and not "
      "Lorentz invariant, so it is excluded")
print("    => unitarity + 4D Lorentz invariance leave EXACTLY two boundary")
print("       conditions: P_R psi(0) = 0 or P_L psi(0) = 0.  The chiral BC is derived.")

# ============================================================ [E]
print("\n[E] The zero mode")
Mv, yy = var('Mv yy', domain='positive')
# massless 4D zero mode psi = f(y) u, g5 u = chi u: (i g^4 d_y + M) psi = 0
#   -> -kappa chi f' + M f = 0 -> f = exp(kappa chi M y)   (weak_chirality_u4 [C])
for chi, name, Pbc in [(-1, "left", PR), (+1, "right", PL)]:
    fz = exp(kappa * chi * Mv * yy)
    normal = bool(limit(fz, yy=oo) == 0)
    u = [c for c in [(PL if chi == -1 else PR).column(j) for j in range(4)] if not c.is_zero()][0]
    allowed = (Pbc * u).is_zero()
    print("    %-5s zero mode: normalisable %s; allowed by BC 'P_%s psi(0)=0': %s"
          % (name, normal, "R" if chi == -1 else "L", allowed))
check(bool(limit(exp(-Mv * yy), yy=oo) == 0) and not bool(limit(exp(Mv * yy), yy=oo) == 0),
      "kappa M > 0: only the LEFT zero mode is normalisable")
uL = [c for c in [PL.column(j) for j in range(4)] if not c.is_zero()][0]
check((PR * uL).is_zero() and not (PL * uL).is_zero(),
      "BC P_R psi(0) = 0 keeps it (one massless left-handed fermion); "
      "BC P_L psi(0) = 0 removes it (no massless mode)")
print("    The remaining freedom is one sign (kappa, equivalently which BC):")
print("    the choice of hand, already ledger Q3 -- the SM makes the same one.")

print("\n" + "=" * 74)
print("""VERDICT
  1. J^4 IS the Dirac form of u4: by Heisenberg, dy/dt = alpha^4 = g^0 g^4,
     whose density is psibar g^4 psi ([A]).  The identification is derived.
  2. On DT's half-line (x^4 = gamma/a0, gamma >= 0) the 4D particle number is
     conserved iff J^4(0) = 0 ([B]); unitarity INDEPENDENTLY forces the same
     ([C]).  Conservation of mass and unitarity coincide at the edge.
  3. Of the U(2) family of unitary boundary conditions, 4D Lorentz invariance
     leaves exactly two, both CHIRAL: psi(0) left- or right-handed ([D]).
  4. Together with the bulk mass (Eq 110 '+1') this leaves one massless
     chiral fermion; its hand is the single remaining sign ([E]) -- ledger
     Q3, the same one choice the SM makes.
  STANDING.  Q4 is no longer an assumption: the no-flux chiral boundary at
  gamma = 0 follows from DT's half-line + conservation of mass, and
  equivalently from unitarity + 4D Lorentz invariance, both SHARED with
  the SM.  Ledger Q4: DT-SPECIFIC -> DERIVED.
  Bracket: (i) it rests on the fermion field living on the mass-density
  axis with gamma >= 0 -- Williams' x^4 = gamma/a0 read literally, a
  half-line with a physical edge; (ii) normalisability in y is assumed
  (as for any localised mode); (iii) a LOCAL, linear boundary condition is
  assumed (standard; non-local conditions are not considered); (iv) the
  distinctive consequences -- heavy mass-axis partner modes, profile
  overlaps -- remain untested, no scale for M being fixed yet.
""")
print("=" * 74)
print("ALL ASSERTIONS PASSED")
print("=" * 74)
