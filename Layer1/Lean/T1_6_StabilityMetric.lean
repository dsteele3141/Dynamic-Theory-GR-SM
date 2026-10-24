/-
  T1.6 — Stability Metric → Space-Entropy Manifold → Minkowski

  Source: P. E. Williams, LA-8370-MS (LASL, 1980), §II, eqs (16)–(19).
  Supplementary Materials: step [L5], Sec. S2 (code index: Sec. S6).
-/

import Mathlib.Tactic
import Mathlib.LinearAlgebra.Matrix.Determinant.Basic

namespace DTWilliams.StabilityMetric

private lemma sq_pos_of_ne (x : ℝ) (hx : x ≠ 0) : 0 < x ^ 2 := by
  rw [← sq_abs]; exact pow_pos (abs_pos.mpr hx) 2

/-! ## The two quadratic forms (diagonal case)

`stabForm` is the stability line element of eq (16) (the metric `h_{ij}`,
positive-definite/Riemannian). `gHatForm` is the spacetime form `ĝ_{ij}` of
eq (19) obtained after solving for (dq⁰)². The coefficients `h1,h2,h3` are the
spatial diagonal entries `h_{αα}`; `f = h_{00}`. -/

/-- Stability Hessian quadratic form `h_{ij} qⁱ qʲ` (diagonal), eq (16). -/
def stabForm (f h1 h2 h3 q0 q1 q2 q3 : ℝ) : ℝ :=
  f * q0 ^ 2 + h1 * q1 ^ 2 + h2 * q2 ^ 2 + h3 * q3 ^ 2

/-- Spacetime quadratic form `ĝ_{ij} xⁱ xʲ = (x⁰)² − h_{αβ} xᵅ xᵝ` (diagonal), eq (19). -/
def gHatForm (h1 h2 h3 x0 x1 x2 x3 : ℝ) : ℝ :=
  x0 ^ 2 - h1 * x1 ^ 2 - h2 * x2 ^ 2 - h3 * x3 ^ 2

/-! ## (1) Solving eq (17) for (dq⁰)²  →  eq (19) -/

/-- **Eq (17) → eq (19) (the solve-for-(dq⁰)² step).**  From the arc-length
    relation `f·(dq⁰)² + Σ h_{αβ}(dqᵅ)² = c²dt²` with `x⁰ = ct` (so `c²dt² = (x⁰)²`)
    and `xᵅ = qᵅ`, one gets `(dq⁰)² = (1/f)·ĝ_{ij} dxⁱ dxʲ`. Pure algebra. -/
theorem solve_for_dq0_sq (f h1 h2 h3 dq0 q1 q2 q3 x0 : ℝ) (hf : f ≠ 0)
    (harc : f * dq0 ^ 2 + (h1 * q1 ^ 2 + h2 * q2 ^ 2 + h3 * q3 ^ 2) = x0 ^ 2) :
    dq0 ^ 2 = (1 / f) * gHatForm h1 h2 h3 x0 q1 q2 q3 := by
  unfold gHatForm
  field_simp
  linarith [harc]

/-! ## (2) The signature flip — the headline of T1.6 -/

/-- The stability form is **positive-definite** (Riemannian, signature +,+,+,+):
    strictly positive off the origin when all diagonal entries are positive. -/
theorem stabForm_posDef (f h1 h2 h3 : ℝ)
    (hf : 0 < f) (hh1 : 0 < h1) (hh2 : 0 < h2) (hh3 : 0 < h3)
    (q0 q1 q2 q3 : ℝ) (hne : q0 ≠ 0 ∨ q1 ≠ 0 ∨ q2 ≠ 0 ∨ q3 ≠ 0) :
    0 < stabForm f h1 h2 h3 q0 q1 q2 q3 := by
  unfold stabForm
  have n0 : 0 ≤ f * q0 ^ 2 := mul_nonneg hf.le (sq_nonneg _)
  have n1 : 0 ≤ h1 * q1 ^ 2 := mul_nonneg hh1.le (sq_nonneg _)
  have n2 : 0 ≤ h2 * q2 ^ 2 := mul_nonneg hh2.le (sq_nonneg _)
  have n3 : 0 ≤ h3 * q3 ^ 2 := mul_nonneg hh3.le (sq_nonneg _)
  rcases hne with h | h | h | h
  · nlinarith [mul_pos hf (sq_pos_of_ne q0 h)]
  · nlinarith [mul_pos hh1 (sq_pos_of_ne q1 h)]
  · nlinarith [mul_pos hh2 (sq_pos_of_ne q2 h)]
  · nlinarith [mul_pos hh3 (sq_pos_of_ne q3 h)]

/-- ĝ is **positive** in the timelike direction `(1,0,0,0)`. -/
theorem gHat_timelike (h1 h2 h3 : ℝ) : 0 < gHatForm h1 h2 h3 1 0 0 0 := by
  simp [gHatForm]

/-- ĝ is **negative** in the spacelike direction `x¹` (and symmetrically x², x³),
    because the positive Hessian entry `h₁ > 0` enters with a minus sign. -/
theorem gHat_spacelike_x (h1 h2 h3 : ℝ) (hh1 : 0 < h1) :
    gHatForm h1 h2 h3 0 1 0 0 < 0 := by simp [gHatForm]; linarith

theorem gHat_spacelike_y (h1 h2 h3 : ℝ) (hh2 : 0 < h2) :
    gHatForm h1 h2 h3 0 0 1 0 < 0 := by simp [gHatForm]; linarith

theorem gHat_spacelike_z (h1 h2 h3 : ℝ) (hh3 : 0 < h3) :
    gHatForm h1 h2 h3 0 0 0 1 < 0 := by simp [gHatForm]; linarith

/-- **Signature flip (T1.6, headline).**  Under the T1.5 stability conditions
    (`f, h₁, h₂, h₃ > 0`), the stability metric is positive-definite (Riemannian)
    yet the spacetime metric `ĝ` is INDEFINITE: positive along the time axis,
    negative along every space axis — signature (+,−,−,−), Lorentzian.
    Thermodynamic stability (positive-definite Hessian) produces a Minkowski-type
    space-time. -/
theorem signature_flip (f h1 h2 h3 : ℝ)
    (hf : 0 < f) (hh1 : 0 < h1) (hh2 : 0 < h2) (hh3 : 0 < h3) :
    -- Riemannian: stability form strictly positive on every axis
    (0 < stabForm f h1 h2 h3 1 0 0 0 ∧ 0 < stabForm f h1 h2 h3 0 1 0 0 ∧
     0 < stabForm f h1 h2 h3 0 0 1 0 ∧ 0 < stabForm f h1 h2 h3 0 0 0 1) ∧
    -- Lorentzian: spacetime form positive on time axis, negative on space axes
    (0 < gHatForm h1 h2 h3 1 0 0 0 ∧ gHatForm h1 h2 h3 0 1 0 0 < 0 ∧
     gHatForm h1 h2 h3 0 0 1 0 < 0 ∧ gHatForm h1 h2 h3 0 0 0 1 < 0) := by
  refine ⟨⟨?_, ?_, ?_, ?_⟩, gHat_timelike h1 h2 h3, gHat_spacelike_x h1 h2 h3 hh1,
    gHat_spacelike_y h1 h2 h3 hh2, gHat_spacelike_z h1 h2 h3 hh3⟩ <;>
  · simp [stabForm]; linarith

/-! ## (2′) Determinant invariants (cross-check with the companion Sage script, SI Appendix A)

The same sign flip seen invariantly: the stability metric `diag(f,h₁,h₂,h₃)` has
positive determinant, the spacetime metric `diag(1,−h₁,−h₂,−h₃)` negative — the
hallmark of a 4D Lorentzian signature (an odd number of negative eigenvalues). -/

/-- `det(h_{ij}) = f·h₁·h₂·h₃ > 0` — Riemannian. -/
theorem stab_det_pos (f h1 h2 h3 : ℝ)
    (hf : 0 < f) (hh1 : 0 < h1) (hh2 : 0 < h2) (hh3 : 0 < h3) :
    0 < (Matrix.diagonal ![f, h1, h2, h3]).det := by
  rw [Matrix.det_diagonal, Fin.prod_univ_four]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons,
    Matrix.cons_val_two, Matrix.cons_val_three, Matrix.tail_cons]
  positivity

/-- `det(ĝ_{ij}) = −h₁·h₂·h₃ < 0` — Lorentzian. -/
theorem gHat_det_neg (h1 h2 h3 : ℝ)
    (hh1 : 0 < h1) (hh2 : 0 < h2) (hh3 : 0 < h3) :
    (Matrix.diagonal ![(1 : ℝ), -h1, -h2, -h3]).det < 0 := by
  rw [Matrix.det_diagonal, Fin.prod_univ_four]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons,
    Matrix.cons_val_two, Matrix.cons_val_three, Matrix.tail_cons]
  nlinarith [mul_pos (mul_pos hh1 hh2) hh3]

/-! ## (3) Minkowski recovery — the flat limit -/

/-- **Minkowski recovery (eq 19, flat limit).**  When the spatial metric is flat,
    `h_{αβ} → δ_{αβ}` (here `h₁=h₂=h₃=1`), the spacetime form `ĝ` becomes exactly
    the special-relativity metric `η = diag(1,−1,−1,−1)`:
        ĝ xⁱ xʲ = (x⁰)² − (x¹)² − (x²)² − (x³)². -/
theorem minkowski_flat (x0 x1 x2 x3 : ℝ) :
    gHatForm 1 1 1 x0 x1 x2 x3 = x0 ^ 2 - x1 ^ 2 - x2 ^ 2 - x3 ^ 2 := by
  unfold gHatForm; ring

/-- The flat-limit metric matrix is literally Minkowski `η = diag(1,−1,−1,−1)`,
    with `det η = −1` (Lorentzian). -/
theorem minkowski_det : (Matrix.diagonal ![(1 : ℝ), -1, -1, -1]).det = -1 := by
  rw [Matrix.det_diagonal, Fin.prod_univ_four]
  simp only [Matrix.cons_val_zero, Matrix.cons_val_one, Matrix.head_cons,
    Matrix.cons_val_two, Matrix.cons_val_three, Matrix.tail_cons]
  norm_num

end DTWilliams.StabilityMetric
