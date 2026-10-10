#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE: the Report (LA-8370-MS) has no neutron model; book only.
#   "The Dynamic Theory" (2011):
#   Eqs 2.33-2.34  unequal central forces F_p = f_p(r), F_e = f_e(r) (no Newton III)
#   Eq 2.41-2.42   E = K + V + K_c + V_c  (internal + centre-of-mass energies)
#   p.42 / p.257   "the energy of the electron emerging after decay is
#                  inconsistent with the equal and opposite coulombic forces ...
#                  Here ... the neo-coulombic forces are unequal ... The above
#                  picture of the neutron produces no need to postulate the
#                  existence of neutrinos."
#   p.253          "the decay of a neutron is accompanied by motion of the
#                  center of mass. It must necessarily move in the decay."
# ---------------------------------------------------------------------------
#
# beta_spectrum_two_body.sage                          [A] (assertion-hardened)
#
# PURPOSE.  In DT the neutron decays to a proton and an electron only.  A
# two-body decay gives a single electron energy; the measured spectrum is
# continuous from 0 to Q = 782 keV.  Williams' rescue is the Newton-III
# violation.  Can it work?
#
#   [A] identity: for ANY central forces f_p(r), f_e(r) the INTERNAL energy is
#       exactly conserved; Newton-III violation only moves the centre of mass
#   [B] so the relative speed at infinity is fixed; the electron's lab energy
#       can vary only through the COM boost, which shows up as proton recoil
#   [C] against data: measured proton recoil (endpoint 751 eV) and the
#       measured continuous electron spectrum
#   [D] the escape route: something else carries the energy (a neutrino)
#
# Imported (not re-derived): allowed beta spectrum shape with the Fermi
# function (Fierz b = 0.017(21), Perkeo III 2019, confirms the allowed
# shape); proton recoil endpoint 751 eV (kinematic, confirmed by aSPECT);
# RDK II radiative branching ratio.
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("beta_spectrum_two_body.sage -- can DT's 2-body decay give a continuous spectrum?")
print("=" * 74)

# ============================================================ [A]
print("\n[A] Identity: Newton-III violation moves the COM, nothing else")
t = var('t'); mp_, me = var('m_p m_e', domain='positive')
xp, xe = function('x_p')(t), function('x_e')(t)
fp, fe = function('f_p'), function('f_e')
M = mp_ + me; mu = mp_ * me / M
r = xp - xe
# U'(r) = -(m_e f_p - m_p f_e)/M : the internal force is a gradient in r alone
U = function('U')
eom = {diff(xp, t, 2): fp(r) / mp_, diff(xe, t, 2): fe(r) / me}
Up = -(me * fp(r) - mp_ * fe(r)) / M              # U'(r)
rdd = diff(r, t, 2).subs(eom)                       # from Eqs 2.33-2.34
dE_explicit = mu * diff(r, t) * rdd + Up * diff(r, t)   # d/dt[mu rdot^2/2 + U]
check(bool(dE_explicit.simplify_full() == 0),
      "d/dt [ mu rdot^2/2 + U(r) ] = 0 for ANY f_p, f_e: internal energy conserved")
V = (mp_ * diff(xp, t) + me * diff(xe, t)) / M
work = (fp(r) * diff(xp, t) + fe(r) * diff(xe, t))
Ptot2 = M / 2 * V**2
check(bool((work - diff(Ptot2, t).subs(eom) + Up * diff(r, t)).simplify_full() == 0),
      "total power = d/dt(M V^2/2) - dU/dt: the ONLY non-conserved energy is "
      "the centre-of-mass kinetic energy P^2/2M")
print("    (1D shown; for central forces the transverse parts are identical.)")

# ============================================================ [B]
print("\n[B] Fixed relative speed => electron energy set by the COM boost alone")
mec2, mpc2, Q = 0.51099895, 938.27208816, 0.78233          # MeV
Mc2 = mpc2 + mec2
# electron in the COM frame carries ~ all of the internal energy (proton heavy)
gu = 1 + Q * mpc2 / Mc2 / mec2
u = sqrt(1 - 1 / gu**2)
Tp_max_meas = 751.4e-6                                     # MeV (recoil endpoint)
Vmax = sqrt(2 * Tp_max_meas / mpc2)                        # COM speed it allows
Te = lambda Vb, c: mec2 * (gu / sqrt(1 - Vb**2) * (1 + u * Vb * c) - 1)
lo, hi = Te(Vmax, -1), Te(Vmax, +1)
print("    COM-frame electron: gamma = %.4f, u = %.4f c" % (gu, u))
print("    proton recoil <= 751 eV  =>  COM speed <= %.2e c" % Vmax)
print("    => DT electron energy confined to [%.1f, %.1f] keV (width %.1f keV)"
      % (1e3 * lo, 1e3 * hi, 1e3 * (hi - lo)))
check(1e3 * (hi - lo) < 5, "DT's electron is a LINE a few keV wide, for any Newton-III violation")
# what proton recoil would a typical measured electron energy need?
def Tp_needed(Te_target):
    Vb = find_root(lambda Vb: Te(Vb, -1) - Te_target, 0, 0.999)
    return mpc2 * (1 / sqrt(1 - Vb**2) - 1)

# ============================================================ [C]
print("\n[C] Against the measured spectra")
def dN(T):                                    # allowed shape x Fermi function (Z=1)
    E = T + mec2; p = sqrt(E**2 - mec2**2)
    eta = E / p / 137.035999
    return p * E * (Q - T)**2 * 2 * pi * eta / (1 - exp(-2 * pi * eta))
norm = numerical_integral(dN, 0, Q)[0]
Tmean = numerical_integral(lambda T: T * dN(T), 0, Q)[0] / norm
frac = numerical_integral(dN, lo, Q)[0] / norm
print("    measured: continuous 0..782 keV, mean T_e = %.0f keV" % (1e3 * Tmean))
print("    fraction of measured decays inside DT's line: %.1e" % frac)
check(frac < 1e-6, "fewer than 1 decay in a million lands where DT's 2-body decay puts them")
Tp_mean = Tp_needed(Tmean)
print("    to give the MEAN electron (%.0f keV) DT needs proton recoil %.1f MeV"
      % (1e3 * Tmean, Tp_mean))
print("    = %.0e x the measured 751 eV endpoint" % (Tp_mean / Tp_max_meas))
check(Tp_mean / Tp_max_meas > 1e4, "spreading the line via the COM needs recoil > 1e4 x measured")
Emiss = Q - Tmean
print("    mean energy unaccounted for by p + e: %.0f keV per decay" % (1e3 * Emiss))
BR_rad = 3.35e-3                              # RDK II, 0.4-14 keV photons
check(BR_rad < 1e-2, "photons carry it in < 1% of decays (RDK II): not radiation")

# ============================================================ [D]
print("\n" + "=" * 74)
print("""VERDICT
  1. Williams' rescue cannot work, and the reason is general.  For ANY pair
     of unequal central forces the internal energy is exactly conserved;
     Newton-III violation only changes the centre-of-mass kinetic energy
     ([A]).  So the electron's energy is fixed up to the COM boost, which
     must appear as proton recoil.  The measured recoil (<= 751 eV) confines
     DT's electron to a line ~3 keV wide at ~780 keV ([B]); 6e-8 of the
     measured decays fall there, and the mean electron (~ 300 keV) would
     need a proton recoil ~ 7e6 x the measured endpoint ([C]).
  2. So a p + e final state is EXCLUDED by DT's own mechanics plus two
     measured spectra.  ~ 0.48 MeV per decay must leave in something
     neutral, not a photon, and nearly invisible.  "No need to postulate
     neutrinos" is refuted -- by DT's own internal-energy conservation.
  3. DT HAS candidates.  Williams' own (the phi_4 wave quantum, Ch.11 Eqs
     11.34-37) is spin 0 and is EXCLUDED by angular momentum; the chiral
     zero mode of DT's 5D Dirac field passes (neutrino_helicity.sage).
     With it, neutron decay is p + e + (DT neutrino), and the spectrum
     argument passes over to the SM's.
  STANDING (assumption ledger: SHARED WITH SM, by data-selected reading).
  Just as Pauli (1930) assumed a light neutral particle carries the missing
  energy, if DT's neutron decay emits its 5D-Dirac zero mode, the spectrum is
  consistent with the SM.  Williams' no-neutrino reading is EXCLUDED; DT
  is not -- the data select the reading with a neutrino (selected, not
  derived).  Ledger Q9.  Its spin and helicity: neutrino_helicity.sage
  (Q10); its strength: Q5.
  Bracket: (i) [A] is non-relativistic in the relative coordinate; the
  proton is non-relativistic (v < 2e-3 c) and the electron enters only
  through the exact radial term; (ii) the allowed shape is imported but
  measured (Fierz b = 0.017(21)); (iii) if the bound state is not an energy
  eigenstate the line broadens, but the neutron's width hbar/tau ~ 1e-18 eV
  makes it sharp.
""")
print("=" * 74)
print("ALL ASSERTIONS PASSED")
print("=" * 74)
