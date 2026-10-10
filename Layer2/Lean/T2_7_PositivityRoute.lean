/-
  Route (a) — does positivity select the FIRST-ORDER (Dirac) equation over the
  SECOND-ORDER (Klein–Gordon) one?

  Source: P. E. Williams, LA-8370-MS (LASL, 1980), Section III.C, eqs (105), (108).
  Full derivation, context and provenance:
  Supplementary_Material.pdf — §S2, Appendix A [C7].
-/

import Mathlib

namespace DTWilliams.PositivityRoute

open Complex
open scoped BigOperators

/-! ### First-order (Dirac) density: positive-definite -/

/-- The Dirac probability density of an `n`-component spinor: `ρ_D = ψ†ψ = Σ|ψ_a|²`.
    First-order (Dirac) evolution conserves exactly this density. (For the 5D Dirac
    equation `n = 8`, eq 109; the result holds for any `n`.) -/
def diracDensity {n : ℕ} (ψ : Fin n → ℂ) : ℝ := ∑ a, normSq (ψ a)

/-- `ρ_D ≥ 0` always: it is a sum of squared magnitudes. -/
theorem diracDensity_nonneg {n : ℕ} (ψ : Fin n → ℂ) : 0 ≤ diracDensity ψ :=
  Finset.sum_nonneg (fun _ _ => normSq_nonneg _)

/-- `ρ_D = 0` exactly when the state is zero — the density is POSITIVE-DEFINITE. -/
theorem diracDensity_eq_zero {n : ℕ} (ψ : Fin n → ℂ) :
    diracDensity ψ = 0 ↔ ψ = 0 := by
  unfold diracDensity
  rw [Finset.sum_eq_zero_iff_of_nonneg (fun _ _ => normSq_nonneg _)]
  constructor
  · intro h
    funext a
    exact normSq_eq_zero.mp (h a (Finset.mem_univ a))
  · intro h _a _
    simp [h]

/-- Any nonzero Dirac state has strictly positive density. -/
theorem diracDensity_pos {n : ℕ} (ψ : Fin n → ℂ) (h : ψ ≠ 0) :
    0 < diracDensity ψ := by
  rcases (diracDensity_nonneg ψ).lt_or_eq with hlt | heq
  · exact hlt
  · exact absurd ((diracDensity_eq_zero ψ).mp heq.symm) h

/-! ### Second-order (Klein–Gordon) density: indefinite -/

/-- The Klein–Gordon conserved density for a single mode of frequency `ω` and
    complex amplitude `c`. Up to a fixed positive constant the time-component of
    the KG current is `j⁰ = Im(ψ* ∂_t ψ)` with `ψ = c·e^{-iωt}`, which equals
    `ω·|c|²`. Both `ω = ±√(k²+m²)` solve the KG dispersion, so `ω` may have either
    sign — the density is NOT sign-definite. -/
def kgDensity (ω : ℝ) (c : ℂ) : ℝ := ω * normSq c

/-- The KG density is NEGATIVE for a (perfectly admissible) negative-frequency
    mode — Dirac's objection to the second-order equation. -/
theorem kgDensity_neg : ∃ (ω : ℝ) (c : ℂ), c ≠ 0 ∧ kgDensity ω c < 0 := by
  refine ⟨-1, 1, one_ne_zero, ?_⟩
  unfold kgDensity
  rw [normSq_one]; norm_num

/-- The KG density takes BOTH signs: positive for `ω > 0`, negative for `ω < 0`. -/
theorem kgDensity_indefinite :
    (∃ (ω : ℝ) (c : ℂ), 0 < kgDensity ω c) ∧
      (∃ (ω : ℝ) (c : ℂ), kgDensity ω c < 0) := by
  refine ⟨⟨1, 1, ?_⟩, ⟨-1, 1, ?_⟩⟩
  · unfold kgDensity; rw [normSq_one]; norm_num
  · unfold kgDensity; rw [normSq_one]; norm_num

/-! ### Positivity selects the first-order form -/

/-- **Positivity selects first order.** The probability density is non-negative for
    EVERY first-order (Dirac) state, but the second-order (Klein–Gordon) density is
    negative for some admissible state. Hence the requirement of a positive-definite
    probability density — which, in the Dynamic Theory, is the consistency
    condition for the entropy/statistical interpretation of the timelike (entropy)
    coordinate — is satisfied by the Dirac equation and VIOLATED by Klein–Gordon.

    This is the machine-checked kernel of Route (a): positivity does not LOGICALLY
    follow from `(dq⁰)² ≥ 0` (a metric condition), but probability-density
    positivity — the theory's own entropy-positivity — selects the first order. -/
theorem positivity_selects_first_order {n : ℕ} :
    (∀ ψ : Fin n → ℂ, 0 ≤ diracDensity ψ) ∧
      (∃ (ω : ℝ) (c : ℂ), c ≠ 0 ∧ kgDensity ω c < 0) :=
  ⟨fun ψ => diracDensity_nonneg ψ, kgDensity_neg⟩

end DTWilliams.PositivityRoute
