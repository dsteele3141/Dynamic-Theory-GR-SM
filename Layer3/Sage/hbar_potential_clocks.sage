#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE (primary): P. E. Williams, "The Dynamic Theory" (2011) Ch. 3,
#   Eqs (3.2)-(3.6) (repeated Ch. 10, Eqs 10.54-10.60): the unit of action
#   depends on the gravitational gauge function,
#       hbar_x = hbar * exp[ (W_x/R_x)(1 + b t_x) e^{-lambda_x/R_x} ],
#       W = -GM/c^2, b = -H0,
#   and photon energy is conserved between emitter and receiver (Eq 3.5),
#       hbar_e nu_e = hbar_r nu_r,
#   which Williams uses to obtain the gravitational redshift (Eq 3.6).
#   Not in LA-8370-MS (book-only material).  The source does NOT say how the
#   electric charge or the masses scale with hbar (alpha is never mentioned):
#   that gap is the fork tested here.
# Context: local_unit_of_action_cancellation (Eq 3.3 drift is potential-
#   suppressed); future-work.txt s.1 (relax the GR-QM boundary).
# ---------------------------------------------------------------------------
#
# hbar_potential_clocks.sage                          [A] (assertion-hardened)
#
# QUESTION.  Relaxing the GR-QM boundary: Eq 3.3 couples a QUANTUM constant
# (hbar) to the GRAVITATIONAL potential at full strength (coefficient 1 in
# Phi/c^2).  Does that predict anything for precision atomic clocks?
#
#   [0] KEY RESULT.  A redshift is measured by comparing the received photon
#       with an identical LOCAL clock.  If a clock runs at nu ~ hbar^p and the
#       photon keeps hbar*nu (Eq 3.5), the observed shift is
#           nu_received / nu_local = (hbar_e/hbar_r)^(p+1).
#       Williams' Eq 3.6 is the special case p = 0 (clock frequencies
#       independent of hbar).  If the clocks themselves scale as 1/hbar
#       (p = -1) the hbar shift CANCELS EXACTLY.
#   [A] Reading A: hbar varies, e and masses fixed -> alpha ~ 1/hbar, clocks
#       have different p: non-universal, alpha tracks Phi.  EXCLUDED.
#   [B] LPI (all clocks the same p) forces alpha and m_e/m_p fixed.  Two
#       consistent readings remain:
#         B1  masses fixed     -> p = -1, the hbar shift cancels; the GR
#             redshift comes from DT's metric (u4 = 0 GR sector).
#         B2  masses ~ hbar    -> p =  0, the hbar shift IS the GR redshift
#             (Williams' Eq 3.6); then the metric must not add a second one.
#   [C] B1 and B2 are the same physics in two unit conventions: counting both
#       gives 2x GR (excluded).
#   [D] Relaxation residues (second order, secular H0*Phi/c^2): unmeasurable.
# ---------------------------------------------------------------------------

from sage.all import *

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("hbar_potential_clocks.sage -- does hbar(Phi) predict a clock signal?")
print("=" * 74)

# Imported bounds (approximate; GR-framework reductions, cited not re-derived)
GPA_REDSHIFT  = 1.4e-4   # Gravity Probe A, Vessot et al. 1980
GAL_REDSHIFT  = 2.5e-5   # Galileo 5/6, Delva et al. 2018 / Herrmann et al. 2018
K_ALPHA_BOUND = 3e-8     # c^2 dln(alpha)/dU from annual clock comparisons
                         # (Yb+ E3/E2, Lange et al. 2021: order 1e-8); generous

GM_S, AU, c = 1.32712440018e20, 1.495978707e11, 2.99792458e8
U_sun = GM_S / (AU * c**2)            # Sun at Earth, ~9.9e-9
ECC = 0.0167
dU_annual = 2 * ECC * U_sun           # peak-to-peak annual modulation

Ue, Ur, U, p = var('U_e U_r U p')    # U = GM/(R c^2) > 0, Phi/c^2 = -U
hb = lambda u: exp(-u)               # Eq 3.3 at t = 0, lambda/R -> 0
gr = 1 - Ue + Ur                     # GR first order: nu_r/nu_e = 1 - (U_e - U_r)
lin = lambda f: taylor(f, (Ue, 0), (Ur, 0), 1)

# ============================================================ [0]
print("\n[0] What a clock comparison actually sees")
# emitter clock: nu_e = nu0 * hbar_e^p ; photon keeps hbar*nu (Eq 3.5):
#   nu_arrive = hbar_e nu_e / hbar_r ; local clock: nu_loc = nu0 * hbar_r^p
nu0 = var('nu0')
nu_e = nu0 * hb(Ue)**p
nu_arrive = hb(Ue) * nu_e / hb(Ur)
nu_loc = nu0 * hb(Ur)**p
obs = simplify(nu_arrive / nu_loc)
check(bool(simplify(obs - (hb(Ue) / hb(Ur))**(p + 1)) == 0),
      "observed shift = (hbar_e/hbar_r)^(p+1) for a clock nu ~ hbar^p")
check(bool(lin(obs.subs(p=0)) == gr),
      "p = 0 (clock frequency independent of hbar): GR redshift -- Williams' Eq 3.6")
check(bool(simplify(obs.subs(p=-1)) == 1),
      "p = -1 (clock nu = E/hbar, E fixed): the hbar shift CANCELS exactly")
shift_coef = lambda pp: (pp + 1)      # first-order redshift in units of GR

# ============================================================ [A]
print("\n[A] Reading A: hbar varies, e, m_e, m_p, c fixed  (alpha ~ 1/hbar)")
h = var('h', domain='positive')
def pw(f): return simplify(diff(log(f), h) * h)
e2A, mA = 1, 1
clocksA = {
    "optical gross (Rydberg m e^4/hbar^3)": mA * e2A**2 / h**3,
    "H maser / Cs hyperfine (~alpha^2 Ry)": (e2A / h)**2 * mA * e2A**2 / h**3,
    "Compton-type (m c^2/hbar)":            mA / h,
}
pA = {}
for name, f in clocksA.items():
    pA[name] = pw(f)
    print("    %-38s p = %s -> hbar shift %+d x GR" % (name, pA[name], shift_coef(pA[name])))
check(len(set(pA.values())) > 1, "clocks have different p: redshift is NOT universal")
# GP-A: H maser vs ground.  Total = metric (1) + hbar part (p+1), p = -5 -> 1 - 4
dev_A = abs(shift_coef(-5))
check(dev_A > 1e4 * GPA_REDSHIFT,
      "H-maser redshift off GR by %d x vs GP-A 1.4e-4: excluded by > 1e4" % dev_A)
k_alpha_A = 1                         # alpha ~ 1/hbar ~ exp(U): c^2 dln(alpha)/dU = 1
check(k_alpha_A / K_ALPHA_BOUND > 1e7,
      "alpha tracks the potential, k_alpha = 1 vs bound ~3e-8: excluded by ~%.0e"
      % (k_alpha_A / K_ALPHA_BOUND))
print("    (annual alpha modulation it would give: %.1e peak-to-peak)" % dU_annual)

# ============================================================ [B]
print("\n[B] Universal (LPI) readings: alpha and m_e/m_p fixed")
for label, m_of_h, want in [("B1: masses fixed  ", 1, -1), ("B2: masses ~ hbar ", h, 0)]:
    e2 = h                                         # alpha = e^2/(hbar c) fixed
    alpha = e2 / h
    clocks = [m_of_h * alpha**2 / h,               # Rydberg
              m_of_h * alpha**4 / h,               # hyperfine (m_e/m_p fixed)
              1 / (h**2 / (m_of_h * e2)),          # material cavity c/L, L ~ Bohr radius
              m_of_h / h]                          # Compton
    ps = [pw(f) for f in clocks]
    print("    %s p for Rydberg, hyperfine, cavity, Compton = %s" % (label, ps))
    check(all(bool(q == want) for q in ps),
          "%s every clock has p = %d -> hbar shift %d x GR, universal" % (label.strip(), want, want + 1))
print("    B1: the hbar shift cancels; the GR redshift must come from the METRIC")
print("        of DT's u4 = 0 GR sector (gr recovery scripts) -- exactly GR.")
print("    B2: the hbar shift IS the GR redshift (Williams' Eq 3.6); masses then")
print("        carry the potential (m ~ exp(-U)), and the metric must not add one.")

# ============================================================ [C]
print("\n[C] One effect, not two")
check(abs(2 - 1) > 1e4 * GAL_REDSHIFT,
      "metric + hbar shift together = 2x GR vs Galileo 5/6 2.5e-5: excluded")
print("    B1 and B2 differ only in which quantities are held fixed (units): every")
print("    dimensionless observable agrees.  So hbar(Phi) is a units convention for")
print("    the gravitational redshift, NOT an extra coupling between the sectors,")
print("    and cannot act as an independent lever (cf. local_unit_of_action_")
print("    cancellation, where it cancels nothing).")

# ============================================================ [D]
print("\n[D] What is left when the boundary is relaxed")
U2_dt = taylor(exp(U), U, 0, 2).coefficient(U, 2)
U2_gr = taylor(1 / sqrt(1 - 2 * U), U, 0, 2).coefficient(U, 2)
d2 = abs(U2_dt - U2_gr)
print("    U^2 coefficient: Eq 3.3 exp form %s, Schwarzschild %s (coordinate-" % (U2_dt, U2_gr))
print("    dependent; under B1 the metric, not Eq 3.3, sets this anyway)")
U_E = 6.96e-10
res2_lab = float(d2) * 2 * U_E * 1.1e-16   # clock pair 1 m apart: dU = g h/c^2
H0 = 2.27e-18
drift = H0 * dU_annual * 3.156e7           # secular hbar drift across the orbit
print("    second order, clock pair 1 m apart  : ~%.0e   (clocks ~1e-18)" % res2_lab)
print("    secular H0*Phi/c^2 annual modulation: ~%.0e /yr" % drift)
check(res2_lab < 1e-20 and drift < 1e-19,   # best clock drift bounds ~1e-18/yr
      "both relaxation residues are far below any present clock comparison")

print("\n" + "=" * 74)
print("""VERDICT
  1. A clock comparison sees (hbar_e/hbar_r)^(p+1) for a clock nu ~ hbar^p
     ([0]).  Williams' Eq 3.6 silently assumes p = 0; a clock that scales
     as 1/hbar sees NO hbar shift at all.  The redshift Eq 3.6 claims is
     therefore not fixed by Eq 3.3 alone: it depends on how charge and mass
     scale with hbar, which the source never states.
  2. Reading A (e, masses fixed; alpha ~ 1/hbar): clocks disagree (p = -3,
     -5, -1), the redshift is non-universal and alpha tracks the solar
     potential.  EXCLUDED by GP-A (> 1e4) and annual clock comparisons
     (~3e7).  A WILLIAMS DISCREPANCY if read this way, not a DT negative.
  3. Data select alpha and m_e/m_p fixed.  Then either masses are fixed
     (p = -1, hbar shift cancels, GR redshift from DT's metric) or masses
     scale with hbar (p = 0, Eq 3.6 is the GR redshift, metric adds none).
     Same physics in two unit conventions; counting both = 2x GR, excluded.
  4. Relaxation residues: ~1e-25 (2nd order), ~2e-20/yr (secular).  This relaxation PREDICTS NOTHING
     MEASURABLE: hbar(Phi) is the gravitational redshift re-expressed, not
     a new GR-QM coupling.
  STANDING.  Source gap (how e and m scale with hbar) closed by data.
  Useful by-products: (i) the variable unit of action cannot be an
  independent lever (e.g. on local Gdot/G); (ii) the GR-QM cross-couplings
  worth testing must be dimensionless ones -- this one is not.
  Bracket: imported clock bounds (GP-A, Galileo, Yb+), approximate;
  first-order scaling; nuclear structure (g_I, m_q) not modelled.
""")
print("=" * 74)
print("ALL ASSERTIONS PASSED")
print("=" * 74)
