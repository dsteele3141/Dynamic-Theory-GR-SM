#!/usr/bin/env sage
# ---------------------------------------------------------------------------
# SOURCE: the Report (LA-8370-MS) has no neutron model; book only.
#   "The Dynamic Theory" (2011):
#   Eqs 2.33-2.34  F_p = -k/r^2 (1-lam_e/r) e^{-lam_e/r},  F_e = +k/r^2 (1-lam_p/r) e^{-lam_p/r}
#   Eq 2.42 / 9.94 V_c = -(k/r) [ e^{-lam_p/r} - e^{-lam_e/r} ]   (COM potential)
#   Fig 2.5 / 9.9  "positive energy well", peak ~ 500 MeV, lam_p = 1 F, lam_e = 1e-3 F
#   Eq 2.60        hbar'_e = 8.0517e-4 hbar, hbar'_p = 0.66586 hbar (fitted to
#                  the neutron spin and magnetic moment)
#   p.282          "a prediction of a half-life of something over 600 seconds"
#   p.~520         "This ability to predict the lifetime of the neutron is new."
#   The tunnelling calculation itself is NOT shown in the book.
# ---------------------------------------------------------------------------
#
# neutron_lifetime_tunnelling.sage                     [A] (assertion-hardened)
#
# PURPOSE.  DT's beta decay is the proton-electron bound state tunnelling out
# of the COM barrier V_c.  Reconstruct the calculation Williams describes and
# ask: does it give the measured lifetime tau_n = 878.4(5) s, and how much of
# the answer is fixed by the theory rather than by choices?
#
#   [A] the barrier: shape, peak, turning points (reproduces Fig 2.5)
#   [B] WKB lifetime for every mass / unit-of-action pair the book offers
#   [C] what it takes to hit 878 s; the input count, against the SM
#   [D] the hydrogen bound on lam_p (hydrogen_nonsingular_spectrum.sage)
#
# Imported (not re-derived): relativistic WKB (Klein-Gordon) under-barrier
# momentum, Gamow attempt-frequency estimate, measured constants, and the SM
# lifetime formula with measured g_A and V_ud (for the like-for-like count).
# ---------------------------------------------------------------------------

from sage.all import *
import mpmath as mp
mp.mp.dps = 30

def check(cond, label):
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    assert cond, "ASSERTION FAILED: " + label

print("=" * 74)
print("neutron_lifetime_tunnelling.sage -- DT's neutron lifetime by tunnelling")
print("=" * 74)

hbarc, k = mp.mpf('197.3269804'), mp.mpf('1.43996448')       # MeV fm
me, mp_, Q = mp.mpf('0.51099895'), mp.mpf('938.27208816'), mp.mpf('0.78233')  # MeV
c_fm = mp.mpf('2.99792458e23')                                # fm/s
tau_meas = 878.4
hp_e, hp_p = mp.mpf('8.0517e-4'), mp.mpf('0.66586')

def Vc(r, lp, le):
    return -(k / r) * (mp.exp(-lp / r) - mp.exp(-le / r))

# ============================================================ [A]
print("\n[A] The COM barrier (Eq 2.42), lam_p = 1 F, lam_e = 1e-3 F")
lp0, le0 = mp.mpf(1), mp.mpf('1e-3')
rpk = mp.findroot(lambda r: mp.diff(lambda x: Vc(x, lp0, le0), r), le0)
Vpk = Vc(rpk, lp0, le0)
print("    peak at r = %.3e F, height %.1f MeV (Fig 2.5: ~500 MeV peak)" % (rpk, Vpk))
check(400 < Vpk < 600, "peak height ~ k/(e lam_e): Fig 2.5 reproduced")
check(Vc(mp.mpf('1e-6'), lp0, le0) < 1e-100 and Vc(mp.mpf(1e4), lp0, le0) < Q,
      "V_c -> 0 at r -> 0 and r -> oo: a single positive barrier, the "
      "'well' is the region r < lam_e INSIDE it")

def turning(lp, le, E):
    pk = mp.findroot(lambda r: mp.diff(lambda x: Vc(x, lp, le), r),
                     (le / 20, 20 * lp), solver='bisect')
    rin = mp.findroot(lambda r: Vc(r, lp, le) - E, (mp.mpf('1e-9'), pk), solver='bisect')
    rout = mp.findroot(lambda r: Vc(r, lp, le) - E, (pk, mp.mpf(1e4)), solver='bisect')
    return rin, pk, rout
rin0, _, rout0 = turning(lp0, le0, Q)
print("    trapped energy = Q = %.4f MeV; turning points %.3e F, %.3f F" % (Q, rin0, rout0))
check(rin0 < le0 < rout0, "the decaying state sits at r < lam_e = 1e-3 F")
print("    NB: Ch.9 places the proton orbit at ~10 F, OUTSIDE the barrier; only")
print("    r < %.0e F is trapped.  The two statements are inconsistent." % rin0)

# ============================================================ [B]
print("\n[B] WKB lifetime  tau = exp(G)/nu,  G = (2/hbar' c) Int sqrt((V-E)(V-E+2mc^2)) dr")
def lifetime(m, hfac, lp=lp0, le=le0, E=Q):
    rin, pk, rout = turning(lp, le, E)
    f = lambda r: mp.sqrt(max(Vc(r, lp, le) - E, 0) * (Vc(r, lp, le) - E + 2 * m))
    G = 2 / (hfac * hbarc) * mp.quad(f, [rin, pk, 10 * pk, rout])
    gam = 1 + E / m
    nu = c_fm * mp.sqrt(1 - 1 / gam**2) / (2 * rin)          # Gamow attempt rate
    return G, mp.log10(nu), (G - mp.log(nu)) / mp.log(10)     # log10 tau
masses = [("m_e (~ reduced mass)", me), ("m_p", mp_), ("M = m_p + m_e (COM)", mp_ + me)]
hbars = [("hbar", mp.mpf(1)), ("hbar'_p", hp_p), ("hbar'_e", hp_e)]
table = {}
print("    %-22s %-9s %10s %12s" % ("mass", "action", "G", "log10 tau/s"))
for mn, m in masses:
    for hn, h in hbars:
        G, lnu, lt = lifetime(m, h)
        table[(mn, hn)] = float(lt)
        print("    %-22s %-9s %10.2f %12.1f" % (mn, hn, G, lt))
lts = list(table.values())
print("    measured: log10(878.4) = %.2f" % log(tau_meas, 10).n())
check(max(lts) - min(lts) > 100,
      "the book's own choices span > 100 decades: the lifetime is not fixed "
      "by the theory as written")
near = [kk for kk, v in table.items() if abs(v - log(tau_meas, 10).n()) < 1]
print("    combinations within a factor 10 of 878 s:", near if near else "NONE")
check(not near, "no combination of the book's ingredients gives ~878 s")

# ============================================================ [C]
print("\n[C] What hitting 878 s takes; inputs counted against the SM")
for mn, m in masses:
    hreq = mp.findroot(lambda h: lifetime(m, h)[2] - mp.log10(tau_meas), (mp.mpf('1e-5'), 1), solver='bisect')
    print("    %-22s needs hbar'/hbar = %.4e" % (mn, hreq))
G1 = lifetime(me, hp_e)[0]
print("    sensitivity: d ln tau / d ln hbar' = -G = %.0f for (m_e, hbar'_e)" % (-G1))
print("    i.e. a 1%% change in hbar'_e moves tau by a factor e^%.1f" % (0.01 * G1))
check(G1 > 100, "the lifetime is exponentially sensitive: a fitted hbar' fixes "
      "it, the theory does not")
hreq_e = mp.findroot(lambda h: lifetime(me, h)[2] - mp.log10(tau_meas), (mp.mpf('1e-5'), 1), solver='bisect')
print("    NEAR-MISS: the natural pairing (electron mass, electron hbar'_e) needs")
print("    hbar'/hbar = %.3e; Williams' hbar'_e (fitted to mu_n, NOT to tau) is" % hreq_e)
print("    %.3e -- a factor %.2f apart, 16 decades in tau." % (hp_e, hreq_e / hp_e))
check(1 < hreq_e / hp_e < 2,
      "independently fitted hbar'_e is within a factor 2 of the value 878 s needs")
# SM like-for-like: G_F (muon decay) plus TWO measured inputs g_A, V_ud
gA, Vud, K = 1.2754, 0.97373, 4908.6          # K: 2 pi^3 hbar^7/(m_e^5 c^4 G_F^2 f(1+RC)) in s
tau_SM = K / (Vud**2 * (1 + 3 * gA**2))
print("    SM: tau = %.0f s from G_F + measured g_A, V_ud (+ radiative corr.)" % tau_SM)
check(abs(tau_SM - tau_meas) / tau_meas < 0.01, "SM with measured g_A, V_ud gives 878 s to < 1%")
print("    DT inputs: lam_p, lam_e (assumed), hbar'_p, hbar'_e (fitted to mu_n,")
print("    spin), plus an unspecified mass and attempt rate.  The SM ALSO takes")
print("    g_A and V_ud from experiment -- so the honest comparison is: the SM")
print("    with 2 measured inputs lands at 878 s; DT with its 4 lands anywhere.")

# ============================================================ [D]
print("\n[D] The hydrogen bound on lam_p (hydrogen_nonsingular_spectrum.sage)")
lpH = mp.mpf('1.8e-3')
print("    hydrogen: lam_p < %.1e F (Williams 1 F).  With lam_p at that bound:" % lpH)
rin, pk, rout = turning(lpH, le0, Q)
print("    barrier peak %.1f MeV, width %.3e F (was %.3f F)" % (Vc(pk, lpH, le0), rout - rin, rout0 - rin0))
check(rout - rin < (rout0 - rin0) / 10,
      "the hydrogen bound shrinks the barrier > 10x: the 'something over 600 s' "
      "rests on the lam_p hydrogen excludes")
lt_H = [float(lifetime(m, h, lp=lpH)[2]) for _, m in masses for _, h in hbars]
print("    log10 tau over the same 9 choices: %.1f ... %.1f" % (min(lt_H), max(lt_H)))
check(max(lt_H) < log(tau_meas, 10).n(),
      "at the hydrogen bound EVERY one of the 9 choices decays faster than 878 s "
      "(hbar'_e included)")

print("\n" + "=" * 74)
print("""VERDICT
  1. The mechanism is real and reconstructible: V_c (Eq 2.42) is a single
     positive barrier peaked at r ~ lam_e (~530 MeV, Fig 2.5), and a state
     at r < lam_e with the decay energy Q = 0.782 MeV must tunnel out ([A]).
  2. The number is NOT fixed.  The book never shows the calculation, and its
     own ingredients (which mass, hbar vs hbar'_p vs hbar'_e) give lifetimes
     spread over > 100 decades, none near 878 s ([B]).  The exponent is so
     large that any fitted hbar' can be tuned to 878 s ([C]).  "Something
     over 600 s" is therefore not a prediction we can reproduce.  One
     fair near-miss: the natural pairing (electron mass with the electron's
     hbar'_e, fitted to the magnetic moment, not to tau) needs hbar' only
     1.54x larger than Williams' value -- suggestive, but the exponent turns
     that 1.54x into 10^16 in tau, so it is not agreement.
  3. The barrier rests on lam_p ~ 1 F, which hydrogen excludes (lam_p < 1.8e-3
     F); at the bound the barrier is ~30x narrower and every choice decays
     faster than 878 s, hbar'_e included ([D]).  The neutron model
     also places the trapped state at r < 1e-3 F, not at the ~10 F orbit used
     for the neutron mass (Ch.9) -- an internal inconsistency.
  STANDING (assumption ledger).  The SM is not input-free either: just as the
  SM takes g_A and V_ud from experiment and then gets 878 s, DT may take
  inputs -- but with lam_p, lam_e, hbar'_p, hbar'_e and an unstated mass and
  attempt rate it does not reproduce 878 s; it can be FITTED to it.  So:
  NOT a postdiction (Williams' "this ability ... is new" is unsupported as
  written), and NOT a falsification.  Ledger: Q7 (lam values, DT-SPECIFIC)
  now also carries the neutron; new row Q8 (decay rate: SM 2 inputs, DT
  unfixed).
  Bracket: (i) WKB for a two-body state whose forces violate Newton's third
  law is our reconstruction; Williams' own calculation may differ; (ii) the
  energy is taken as Q (the observable release); (iii) a 2-body p+e final
  state gives a MONOENERGETIC electron unless the field carries momentum --
  the measured continuous beta spectrum is the sharper next test.
""")
print("=" * 74)
print("ALL ASSERTIONS PASSED")
print("=" * 74)
