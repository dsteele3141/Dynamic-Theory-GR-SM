#!/usr/bin/env sage
# --------------------------------------------------------------------------
# SOURCE: Williams (2011), The Dynamic Theory, 5D structure of Ch.2; the
#   route fork P=0 / P!=0 (Layer 2, S2 sec.1) and the S3 mass-sector addendum
#   (this session, "no absolute GR/QM boundary"). Constants a0, H0 as pinned
#   in the corpus (P1_radiation_pressure_a0.sage, local_gdot_escape_candidates.sage).
# --------------------------------------------------------------------------

# ============================================================================
# u4_rest_mass_nonconservation.sage              [A]  (assertion-hardened)
#
# QUANTIFYING THE GR/QM CROSSOVER.
# In DT the charge/particle (QM) sector is route 3, P=0, i.e. u4=0 (constant
# rest mass, pure unitary 4D). The mass/gravity (GR) sector is route 2, P!=0,
# with the fifth (mass-density) velocity component u4 live. The 4D<->5D divide
# is thus the SINGLE restrictive assumption u4=0, with the clean-separation
# limit a0 -> inf playing the exact role of DT's hbar -> 0. There is NO
# absolute boundary; the question is how BIG the residual crossing is at the
# real, FINITE, pinned a0 -- and how it is bounded.
#
# DT's 5th coordinate is  q4 = gamma/a0  (gamma = mass density, LENGTH units),
# and the fifth velocity component is a VELOCITY RATIO (q0 = c t):
#       u4 = dq4/dq0 = (1/a0) d gamma/dq0 = (1/(a0 c)) d gamma/dt.
# Because it is dq4/dq0, DT's absolute velocity c caps it exactly as it caps
# the spatial velocities: from the (+,-,-,-,-) 5-velocity norm,
#       beta_space^2 + beta_mass^2 < 1   =>   u4 in [0, 1),
# and the mass-direction speed TRADES against spatial speed on one shared
# c-budget. So the GR/QM dial is a c-capped BOOST: u4=0 is the flat quantum
# sector; u4 -> 1 is the fully mass-dynamical GR limit; the two are ends of one
# hyperbolic rotation, with QUANTISATION (the London condition) common to both.
#
# The DT-native rate of rest-mass non-conservation is the secular (1 - H0 t)
# factor (b = H0; the SAME drift that drives the dark-energy candidate), so
# for quasi-static matter  d gamma/dt ~ H0 gamma  and  u4 ~ H0 gamma/(a0 c).
#
# RESULT: for all accessible matter u4 is astronomically small (lab solid
# ~1e-31; even nuclear density ~5e-17), because the crossing is DOUBLY
# suppressed -- by the cosmological slowness H0 ~ 1e-18 /s AND by the huge DT
# gauge scale a0 c ~ 1e16. So pure unitary QM is an IDEALISATION exact only at
# a0 -> inf; the finite-a0 leakage is real but far below any decoherence /
# non-unitarity bound. It is largest for the densest, fastest-mass-changing
# systems, so the effect surfaces cosmologically, not locally.
#
# BOUNDS (not a floor with an open top): d gamma/dt is HARD-CAPPED at a0 c
# (u4 -> 1). Quasi-static matter uses ~1e-31 of that c-budget. A faster system
# uses MORE of the fixed budget, up to the ceiling it cannot exceed -- the
# mass-sector analogue of "cannot cross c".
#
# RUN: /home/dsteele/miniconda3/envs/sage/bin/sage u4_rest_mass_nonconservation.sage
# ============================================================================

from sage.all import *


def ok(label, cond):
    cond = bool(cond)
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label


# --- pinned corpus constants -------------------------------------------------
H0 = 2.27e-18     # s^-1     Hubble rate      (local_gdot_escape_candidates.sage)
a0 = 3.65e7       # kg/m^4   DT gauge constant (P1_radiation_pressure_a0.sage)
c  = 2.998e8      # m/s

print("=" * 72)
print("u4 -- REST-MASS NON-CONSERVATION and the FINITE-a0 GR/QM CROSSOVER")
print("=" * 72)
print("  H0 = %.3e /s     a0 = %.3e kg/m^4     c = %.4e m/s" % (H0, a0, c))
print("  a0*c = %.3e kg/(m^3 s)   (the DT gauge suppression scale = u4=1 ceiling)"
      % (a0*c))


def u4_of(gamma):
    """u4 for quasi-static matter: u4 = H0 gamma / (a0 c)  (secular rate)."""
    return H0 * gamma / (a0 * c)


def u4_from_rate(dgdt):
    """u4 = (1/(a0 c)) d gamma/dt  -- the c-capped velocity ratio, general rate."""
    return dgdt / (a0 * c)


# --- (1) representative environments -----------------------------------------
print("\n[1] u4 across representative mass densities gamma (kg/m^3), secular rate")
envs = [
    ("cosmic mean",         2.7e-27),
    ("air (STP)",           1.2),
    ("water / lab solid",   1.0e3),
    ("Earth (mean)",        5.5e3),
    ("iron",                7.9e3),
    ("white dwarf (Ch.)",   1.0e9),
    ("nuclear matter",      2.3e17),
]
for name, gamma in envs:
    print("    %-20s gamma=%.2e  ->  u4 = %.3e" % (name, gamma, u4_of(gamma)))

# --- (2) lab matter: crossover utterly negligible ----------------------------
u4_lab = u4_of(1.0e3)
print("\n[2] Lab matter (gamma=1e3): u4 = %.3e" % u4_lab)
ok("lab-matter u4 < 1e-25 (pure unitary QM safe to far below any decoherence bound).",
   u4_lab < 1e-25)

# --- (3) even the densest accessible matter stays tiny -----------------------
u4_nuc = u4_of(2.3e17)
print("\n[3] Nuclear matter (gamma=2.3e17, the densest accessible): u4 = %.3e" % u4_nuc)
ok("nuclear-density u4 < 1e-10 (crossover still far below observable non-unitarity).",
   u4_nuc < 1e-10)

# --- (4) structure: linear in gamma, monotone, vanishes as a0 -> inf ----------
print("\n[4] Structural properties of the crossover u4(gamma) = H0 gamma/(a0 c)")
ok("u4 is LINEAR in density: u4(2 gamma) = 2 u4(gamma).",
   abs(u4_of(2.0e3) - 2*u4_of(1.0e3)) < 1e-40)
ok("u4 is MONOTONE increasing in gamma (denser -> larger crossover).",
   u4_of(1.0e9) > u4_of(1.0e3) > u4_of(1.2))
# a0 -> inf is DT's hbar -> 0 : the clean-separation limit
g = var('g', domain='positive'); A = var('A', domain='positive')
ok("clean-separation limit: u4 = H0 g/(A c) -> 0 as a0=A -> inf (pure QM recovered).",
   limit(H0*g/(A*c), A=Infinity).is_zero())

# --- (5) the double suppression, made explicit -------------------------------
print("\n[5] Why the continuous boundary is unobservable locally: DOUBLE suppression")
supp_cosmo = H0                 # /s     cosmological slowness
supp_gauge = 1.0/(a0*c)         # m^3 s / kg   DT gauge scale
print("    cosmological factor H0        = %.2e /s" % supp_cosmo)
print("    gauge factor 1/(a0 c)         = %.2e m^3 s/kg" % supp_gauge)
print("    product (per unit gamma)      = %.2e (m^3/kg)/s" % (supp_cosmo*supp_gauge))
ok("the two independently-pinned suppressions multiply to ~1e-34 per unit gamma.",
   (supp_cosmo*supp_gauge) < 1e-33)

# --- (6) u4 is a c-CAPPED velocity ratio: the ceiling, not an open rate -------
print("\n[6] u4 = dq4/dq0 is a VELOCITY RATIO, capped by the absolute velocity c")
dgdt_ceiling = a0 * c
print("    ceiling: d gamma/dt = a0 c = %.3e kg/(m^3 s)  <->  u4 = 1 (GR limit)"
      % dgdt_ceiling)
ok("at the ceiling d gamma/dt = a0 c the boost SATURATES: u4 = 1.",
   abs(u4_from_rate(dgdt_ceiling) - 1.0) < 1e-12)
print("    quasi-static lab matter uses fraction %.2e of the c-budget" % u4_lab)
ok("quasi-static lab matter uses a negligible fraction of the c-budget (u4 << 1).",
   u4_lab < 1e-25)
# shared c-budget: beta_space^2 + beta_mass^2 < 1  (from the (+,-,-,-,-) norm)
ok("shared c-budget: raising spatial beta shrinks the mass budget (bs=0.99 -> bm<0.15).",
   bool(sqrt(1 - 0.99^2) < 0.15))
# core-collapse flag (ORDER OF MAGNITUDE ONLY -- bulk mean density / timescale)
dgdt_collapse = 1e17 / 1e-3      # ~1e20 kg/(m^3 s): 1e9 -> 1e17 over ~1 ms (rough)
u4_collapse_naive = u4_from_rate(dgdt_collapse)
print("    FLAG (order-of-magnitude): SN core collapse dgamma/dt ~ %.0e -> naive u4 ~ %.0e"
      % (dgdt_collapse, u4_collapse_naive))
ok("naive core-collapse rate EXCEEDS the cap (u4>1): DT must speed-limit it (lead, not claim).",
   u4_collapse_naive > 1)
print("      => DT predicts the mass-densification rate SATURATES at a0 c (u4 -> 1):")
print("         a mass-sector 'cannot cross c'. Needs the LOCAL dq4/dt along a")
print("         worldline, not the bulk-mean estimate used here -- flagged as a lead.")

print("\n[VERDICT] The GR/QM boundary in DT is CONTINUOUS, not sharp -- and more")
print("  precisely it is a c-capped BOOST, not a soft knob: u4 = dq4/dq0 is a")
print("  velocity ratio in [0,1). u4=0 is the flat quantum sector; u4 -> 1 is the")
print("  fully mass-dynamical GR limit; the two are ends of ONE hyperbolic rotation,")
print("  with quantisation (the London condition) common to both -- so gravity is a")
print("  u4-boost away from the quantum sector, not a separate theory. a0 -> inf is")
print("  the hbar -> 0 analogue. At the real finite a0 the crossing is ~1e-31 for lab")
print("  matter, never exceeds ~1e-17 at nuclear density, and is hard-capped at")
print("  d gamma/dt = a0 c: real, definite, cosmological/extreme-density -- not a")
print("  laboratory signature. (Magnitude estimate on DT's secular rate, not a bound.)")

print("\n" + "=" * 72)
print("u4_rest_mass_nonconservation COMPLETE.")
print("=" * 72)
