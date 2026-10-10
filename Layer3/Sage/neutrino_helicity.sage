#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE: Report (LA-8370-MS) Eqs 105-112 (5D Dirac); book (2011):
#   Ch.11 Eqs 11.34-11.37  the neutrino is the quantum of the non-transverse
#                  phi_4 wave ("This particle has all of the attributes we
#                  currently assign to a neutrino."); its spin is never stated
#   Eq 2.57 / 9.44 (1/2)(+-h_e +- h_p) + h_e + h_p = h : the neutron's
#                  orbital + intrinsic angular momentum set equal to that of
#                  free p + e, which Williams takes as h/2 + h/2 = h (p.49)
#   p.257          "[Newton III] was found to be violated during beta decay
#                  until the hypothesis of the neutrino reinstated the
#                  summation of particle spins"
#   Eq 2.60        h'_e = 8.0517e-4 h, h'_p = 0.66586 h
# ---------------------------------------------------------------------------
#
# neutrino_helicity.sage                               [A] (assertion-hardened)
#
# PURPOSE.  beta_spectrum_two_body.sage showed neutron decay needs a neutral
# third body.  Which DT object can it be, judged by SPIN and HELICITY?
#
#   [A] angular momentum: n -> p + e + X forces X to be a FERMION
#   [B] Williams' phi_4 quantum: integer spin, helicity 0
#   [C] DT's other candidate, the u4 = 0 chiral zero mode of the 5D Dirac
#       field (weak_chirality_u4.sage): spin 1/2, helicity computed here
#   [D] Williams' own neutron: Eq 2.57 gives it angular momentum h, not h/2
#
# Measured (imported): Goldhaber-Grodzins-Sunyar 1958, nu_e helicity
# negative ("100% within errors"); pi -> mu nu, xi P_mu delta/rho > 0.99682
# (TWIST 2011) => nu_mu helicity -1 to < 0.3%; neutron decay neutrino
# asymmetry B = 0.98 (antineutrino right-handed); neutron spin 1/2 (neutron
# Stern-Gerlach, two beams: Sherwood, Stephenson, Bernstein 1954).
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("neutrino_helicity.sage -- spin and helicity of DT's neutrino candidates")
print("=" * 74)

# ============================================================ [A]
print("\n[A] Angular momentum in n -> p + e + X (all in units of hbar)")
half = QQ(1) / 2
def allowed_J(spins, L_max=4):
    """total J reachable from the given spins plus orbital L (as 2J)"""
    tot = {0}
    for s in spins:
        tot = {abs(a - int(2 * s)) + 2 * k for a in tot
               for k in range(0, min(a, int(2 * s)) + 1)}
    return {j + 2 * l for j in tot for l in range(L_max + 1)} | \
           {abs(j - 2 * l) for j in tot for l in range(L_max + 1)}
def half_int_ok(sX):
    return 1 in allowed_J([half, half, sX])           # 2J = 1: spin 1/2 neutron
check(not half_int_ok(0) and not half_int_ok(1) and not half_int_ok(2),
      "X with spin 0, 1, 2 (any boson): p + e + X cannot have J = 1/2")
check(half_int_ok(half), "X with spin 1/2: J = 1/2 reachable -- X must be a fermion")
check(not (1 in allowed_J([half, half])),
      "no X at all (Williams' p + e): J = 1/2 unreachable")

# ============================================================ [B]
print("\n[B] Williams' neutrino: the phi_4 wave quantum (Eqs 11.34-11.37)")
# phi_4 is the FIFTH component of the 5D gauge vector: under 4D rotations
# (indices 0..3) it does not transform, so its quantum is a 4D SCALAR.
Rz = matrix(SR, 5, 5, lambda i, j: 0); th = var('th')
Rz[0, 0] = Rz[3, 3] = Rz[4, 4] = 1
Rz[1, 1] = Rz[2, 2] = cos(th); Rz[1, 2] = -sin(th); Rz[2, 1] = sin(th)
phi = vector(SR, [0, 0, 0, 0, var('B4')])            # only phi_4 nonzero (Eq 11.36)
check(Rz * phi == phi, "a phi_4-only gauge vector is invariant under 4D rotations: spin 0")
s_phi4, h_phi4 = 0, 0
print("    spin: 0 (4D scalar); helicity along k: %d" % h_phi4)
check(s_phi4 == 0 and not half_int_ok(s_phi4),
      "phi_4 quantum is a boson: EXCLUDED by angular momentum in neutron decay")
check(h_phi4 != -1, "and its helicity 0 contradicts the measured -1 (Goldhaber; pi decay)")

# ============================================================ [C]
print("\n[C] DT's other candidate: the 5D-Dirac chiral zero mode (u4 = 0 boundary)")
I2, Z2 = identity_matrix(SR, 2), zero_matrix(SR, 2)
s = [matrix(SR, [[0, 1], [1, 0]]), matrix(SR, [[0, -I], [I, 0]]),
     matrix(SR, [[1, 0], [0, -1]])]
Bk = lambda a, b, c, d: block_matrix([[a, b], [c, d]], subdivide=False)
g = [Bk(I2, Z2, Z2, -I2)] + [Bk(Z2, si, -si, Z2) for si in s]
g5 = I * g[0] * g[1] * g[2] * g[3]
I4 = identity_matrix(SR, 4)
PL = (I4 - g5) / 2
Sz = Bk(s[2], Z2, Z2, s[2])                         # spin along z (x2)
E = var('E', domain='positive')
slash = lambda p: g[0] * p[0] - g[3] * p[3]
# particle: massless u(p), p along +z, left-chiral
Nu = (slash([E, 0, 0, E])).right_kernel().basis_matrix().T
uL = [PL * Nu.column(i) for i in range(Nu.ncols())]
uL = [v for v in uL if not v.is_zero()]
check(len(uL) == 1, "massless left-chiral solution moving along +z is unique")
u = uL[0]
hel_nu = [ (Sz * u)[i] / u[i] for i in range(4) if u[i] != 0 ][0]
print("    neutrino (left-chiral particle): helicity = %s" % hel_nu)
check(hel_nu == -1, "zero-mode NEUTRINO has helicity -1, as measured (Goldhaber; pi decay)")
# antiparticle.  A left-chiral FIELD creates antiparticles through a v-spinor
# with P_L v = v; massless u and v solve the same equation p-slash = 0, and
# the antiparticle's PHYSICAL spin is opposite to the spinor's (hole theory).
# Cross-check by charge conjugation: v_L = i gamma^2 (u_R)^*.
uR = [w for w in [(I4 - PL) * Nu.column(i) for i in range(Nu.ncols())] if not w.is_zero()][0]
v = vector(SR, I * g[2] * uR.conjugate())
check((PL * v - v).is_zero() and (slash([E, 0, 0, E]) * v).is_zero(),
      "v_L = i gamma^2 u_R^* is left-chiral and solves the +z Dirac equation")
hel_nubar = [ -(Sz * v)[i] / v[i] for i in range(4) if v[i] != 0 ][0]   # physical spin = -Sz
print("    antineutrino (in neutron decay): helicity = %s" % hel_nubar)
check(hel_nubar == +1, "zero-mode ANTINEUTRINO has helicity +1 (neutron-decay B ~ 1: right-handed)")
check(half_int_ok(half), "and spin 1/2 satisfies [A]")
print("    NB: the zero mode is massless; oscillations need m_nu > 0 -- the")
print("    minimal SM is also massless, so this is shared, not DT-specific.")

# ============================================================ [D]
print("\n[D] Williams' neutron: what angular momentum does Eq 2.57 give it?")
he, hp = QQ(80517) / 10**8, QQ(66586) / 10**5
combos = {(a, b): half * (a * he + b * hp) + he + hp for a in (1, -1) for b in (1, -1)}
best = min(combos, key=lambda k: abs(combos[k] - 1))
print("    signs %s give J_n = %.5f hbar (Eq 2.57 target: hbar = h/2 + h/2)" % (best, combos[best]))
check(abs(combos[best] - 1) < 1e-4, "Eq 2.57 sets the neutron's angular momentum to hbar")
Jn = round(combos[best])
check(2 * Jn + 1 == 3 and 2 * half + 1 == 2,
      "J = hbar means 3 magnetic substates; the measured neutron (J = hbar/2) "
      "shows 2 Stern-Gerlach beams: Williams' neutron has the WRONG SPIN")
check(not (1 in allowed_J([half, half])) and (2 in allowed_J([half, half], 0)),
      "and two free spin-1/2 particles give J = 0 or hbar, never hbar/2: "
      "his p + e final state cannot come from a spin-1/2 neutron")
print("    (Internally consistent with his p + e decay -- J = 1 -> p + e needs no")
print("    neutrino -- which is exactly why the measured spin 1/2 forces one.)")

print("\n" + "=" * 74)
print("""VERDICT
  1. Williams' neutrino (the phi_4 wave quantum) is EXCLUDED twice: it is a
     boson, and angular momentum in n -> p + e + X demands a fermion ([A],
     [B]); its helicity 0 contradicts the measured -1.  So the neutral
     fermion beta_spectrum_two_body requires (ledger Q9) is NOT his.
  2. DT's 5D Dirac equation supplies a candidate that PASSES: the chiral zero
     mode left by the u4 = 0 no-flux boundary (weak_chirality_u4) is spin
     1/2 with neutrino helicity -1 and antineutrino helicity +1, exactly the
     measured values ([C]).  It needs the same choices already logged: one
     4x4 irrep (Q3, shared with the SM) and the boundary reading (Q4, ours).
  3. Williams' own neutron (Eq 2.57) has J = hbar.  That is why his p + e
     decay looked self-consistent; the measured neutron spin 1/2 breaks it
     ([D]) -- a second, independent reason the no-neutrino reading fails.
  STANDING (assumption ledger: SHARED WITH SM, SELECTED by data).  Just as
  the SM assumes a massless left-handed Weyl neutrino, if DT's neutrino is
  its 5D Dirac zero mode, its spin and helicity are consistent with the SM.
  The data select the zero mode over Williams' phi_4 quantum (selected,
  not derived).  Ledger Q10 neutrino identity; Q11 Williams'
  neutron spin (DT-SPECIFIC model, negative candidate).
  Bracket: (i) the zero-mode neutrino is OUR identification, built on our
  boundary reading; (ii) "integer spin" for the phi_4 quantum follows from
  it being a classical-field quantum -- Williams never states a spin; (iii)
  [D]'s defeater is Williams' effective unit of action: he could argue J is
  counted in h' inside the neutron, but a FREE neutron shows 2 states, so
  the defeater must explain why it is spin 1/2 outside.
""")
print("=" * 74)
print("ALL ASSERTIONS PASSED")
print("=" * 74)
