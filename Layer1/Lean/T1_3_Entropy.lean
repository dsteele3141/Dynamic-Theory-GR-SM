/-
  T1.3 — Mechanical Entropy and the Clausius Inequality

  Source: P. E. Williams, LA-8370-MS (LASL, 1980), §II, eq (6) and the Clausius theorem.
  Supplementary Materials: step [L3], Sec. S2 (code index: Sec. S6).
-/

import Mathlib.Tactic
import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus

namespace DTWilliams.Entropy

variable (φ f : ℝ → ℝ)

/-- Coefficient of the reversible energy 1-form: `đE = (φ(q̇)·f(σ)) dσ`
    (from T1.1, λ = φ·f). -/
def dEcoeff (q σ : ℝ) : ℝ := φ q * f σ

/-- **Entropy form is exact (eq 6, doc p.20-21).** Dividing đE by the universal
    integrating-factor velocity part `φ(q̇)` leaves `f(σ)`, a function of σ ALONE.
    This is what makes `dS = đE/φ` an exact differential. -/
theorem entropy_form_eq_f (q σ : ℝ) (hq : φ q ≠ 0) :
    dEcoeff φ f q σ / φ q = f σ := by
  unfold dEcoeff; rw [mul_div_cancel_left₀ _ hq]

/-- **Velocity-independence (exactness, restated).** `đE/φ` is the same function
    of σ at any two velocities — it carries no velocity dependence, so dS is
    well-defined as a state quantity. -/
theorem entropy_form_velocity_indep (q q' σ : ℝ) (hq : φ q ≠ 0) (hq' : φ q' ≠ 0) :
    dEcoeff φ f q σ / φ q = dEcoeff φ f q' σ / φ q' := by
  rw [entropy_form_eq_f φ f q σ hq, entropy_form_eq_f φ f q' σ hq']

/-- **Entropy is unchanged on an E-conservative curve (doc p.20).** A reversible
    E-conservative process keeps σ fixed (dσ = 0), so with `S = F(σ)` the entropy
    change between coincident σ-endpoints is zero. -/
theorem entropy_const_on_Econservative (F : ℝ → ℝ) (σ : ℝ) :
    F σ - F σ = 0 := sub_self _

/-- **Reversible-path entropy (doc p.21).** For a reversible path the entropy
    change equals the σ-integral of `f`:
        `∫_{σ_A}^{σ_B} f dσ = F(σ_B) − F(σ_A) = S(B) − S(A)`,
    where the entropy `S = F` is an antiderivative of `f` (`đE/φ = f dσ`).
    This is the fundamental theorem of calculus. -/
theorem entropy_reversible_path (F : ℝ → ℝ) (σA σB : ℝ)
    (hF : ∀ x ∈ Set.uIcc σA σB, HasDerivAt F (f x) x)
    (hf : IntervalIntegrable f MeasureTheory.volume σA σB) :
    ∫ σ in σA..σB, f σ = F σB - F σA :=
  intervalIntegral.integral_eq_sub_of_hasDerivAt hF hf

/-- **Entropy-increase law for isolated systems (doc p.22) — DERIVED from the
    Clausius bound.** Williams asserts the dynamic second law (Clausius
    inequality) without proof: `heat := ∫_A^B đE/φ ≤ S(B) − S(A)`. For an
    isolated system `đE = 0`, so `heat = 0` and `S(B) ≥ S(A)`: the mechanical
    entropy of an isolated system never decreases. -/
theorem entropy_increase_isolated (heat S_A S_B : ℝ)
    (clausius : heat ≤ S_B - S_A) (isolated : heat = 0) :
    S_A ≤ S_B := by linarith

end DTWilliams.Entropy
