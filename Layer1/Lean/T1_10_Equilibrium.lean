/-
  T1.10 — Equilibrium Conditions (Mechanical Analogues)

  Source: P. E. Williams, LA-8370-MS (LASL, 1980), §II, eqs (11)–(13).
  Supplementary Materials: step [L3], Sec. S2 (code index: Sec. S6).
-/

import Mathlib.Tactic

namespace DTWilliams.Equilibrium

/-- **Entropy maximum for an E-conservative system (Section A).** With the
    Clausius bound `heat ≤ ΔS` (T1.3) and `đE = 0` (so the heat integral
    `heat = 0`), the entropy cannot decrease: `ΔS ≥ 0`. -/
theorem entropy_max_Econservative (heat ΔS : ℝ)
    (clausius : heat ≤ ΔS) (econs : heat = 0) : 0 ≤ ΔS := by
  rw [econs] at clausius; exact clausius

/-- **Mechanical Clausius inequality (Section B).** Substituting the first law
    `ΔE = ΔU − F·Δq` into Clausius `ΔE ≤ φ·ΔS` gives `ΔU − F·Δq ≤ φ·ΔS`. -/
theorem mechanical_clausius (φ ΔS ΔU F Δq ΔE : ℝ)
    (firstlaw : ΔE = ΔU - F * Δq) (clausius : ΔE ≤ φ * ΔS) :
    ΔU - F * Δq ≤ φ * ΔS := by rw [firstlaw] at clausius; exact clausius

/-- **Free energy (eq 13).** For a reversible process `ΔE_rev = φ·ΔS`, the
    available work is `A_max = φ·ΔS − ΔU + F·Δq`. -/
theorem free_energy_max (φ ΔS ΔU F Δq ΔErev : ℝ) (hrev : ΔErev = φ * ΔS) :
    ΔErev - ΔU + F * Δq = φ * ΔS - ΔU + F * Δq := by rw [hrev]

/-- **Least work / eq (12) identity.** The least work to displace from
    equilibrium equals the maximum free energy recoverable from the (reverse)
    spontaneous restoration: `δA_min = φ·(−δS) − (−δU) + F·(−δq) = δU − F·δq − φ·δS`. -/
theorem least_work_eq12 (φ δS δU F δq : ℝ) :
    φ * (-δS) - (-δU) + F * (-δq) = δU - F * δq - φ * δS := by ring

/-- **Equilibrium criterion / free-energy minimum (eq 12).** The restoration of a
    virtual displacement is spontaneous, so Clausius applies to it:
    `−δU + F·δq ≤ −φ·δS`. Rearranging, `δU − F·δq − φ·δS ≥ 0` — any displacement
    from equilibrium requires work (the mechanical free energy is a minimum). -/
theorem equilibrium_criterion (φ δS δU F δq : ℝ)
    (restore_clausius : -δU + F * δq ≤ -(φ * δS)) :
    0 ≤ δU - F * δq - φ * δS := by linarith

/-- **Connection to the Helmholtz function K (Section E).** At constant q̇,
    `δK = δU − φ·δS` (from T1.4's `dK`). Hence `δA_min = δU − F·δq − φ·δS
    = δK − F·δq`, so the equilibrium criterion `δA_min ≥ 0` is the statement that
    `K − F·q` is minimised at equilibrium (the mechanical Gibbs minimum). -/
theorem equilibrium_via_K (φ δS δU F δq δK : ℝ) (hK : δK = δU - φ * δS) :
    δU - F * δq - φ * δS = δK - F * δq := by rw [hK]; ring

end DTWilliams.Equilibrium
