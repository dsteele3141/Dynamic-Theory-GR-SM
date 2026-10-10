/-
  T2.1 — 5D Minkowski Metric and Equations of Motion

  Source: P. E. Williams, LA-8370-MS (LASL, 1980), Section III, eqs (58)-(63).
  Full derivation, context and provenance:
  Supplementary_Material.pdf — §S2, Appendix A [C1].
-/

import Mathlib.Tactic
import Mathlib.Analysis.Calculus.Deriv.Mul
import Mathlib.Analysis.SpecialFunctions.Sqrt

namespace DTWilliams.Minkowski5

/-- The 5D Minkowski quadratic form (eq 60), signature (+,−,−,−,−). -/
def Q5 (x0 x1 x2 x3 x4 : ℝ) : ℝ := x0 ^ 2 - x1 ^ 2 - x2 ^ 2 - x3 ^ 2 - x4 ^ 2

/-! ## (1) 5D metric: eq (58) → eq (60), and signature -/

/-- **Eq (58) → eq (60).** From the 5D Euclidean arc-length
    `(dq⁰)² + Σ(dqᵅ)² + (dq⁴)² = c²dt²` with `x⁰ = ct` (so `c²dt² = (x⁰)²`),
    solving for `(dq⁰)²` gives the 5D Minkowski form. -/
theorem solve_dq0_5D (dq0 dq1 dq2 dq3 dq4 x0 : ℝ)
    (harc : dq0 ^ 2 + (dq1 ^ 2 + dq2 ^ 2 + dq3 ^ 2 + dq4 ^ 2) = x0 ^ 2) :
    dq0 ^ 2 = Q5 x0 dq1 dq2 dq3 dq4 := by
  unfold Q5; linarith

/-- `Q5` is **positive** in the timelike (entropy) direction. -/
theorem Q5_timelike : 0 < Q5 1 0 0 0 0 := by unfold Q5; norm_num

/-- `Q5` is **negative** in each of the three space directions (here x¹). -/
theorem Q5_spacelike_x1 : Q5 0 1 0 0 0 < 0 := by unfold Q5; norm_num
theorem Q5_spacelike_x2 : Q5 0 0 1 0 0 < 0 := by unfold Q5; norm_num
theorem Q5_spacelike_x3 : Q5 0 0 0 1 0 < 0 := by unfold Q5; norm_num

/-- **The mass-density direction `x⁴` is spacelike** (`Q5 < 0`): the 5th
    coordinate enters the metric with a minus sign, like the spatial ones. -/
theorem Q5_mass_spacelike : Q5 0 0 0 0 1 < 0 := by unfold Q5; norm_num

/-- **Determinant of the 5D Minkowski metric = +1** (the product of the diagonal
    `diag(1,−1,−1,−1,−1)`). Contrast 4D Lorentzian (det = −1): the even number
    (four) of spacelike directions makes the 5D determinant positive. -/
theorem minkowski5_det : (1 : ℝ) * (-1) * (-1) * (-1) * (-1) = 1 := by norm_num

/-! ## (2) Velocity normalisation `g_{ij}uⁱuʲ = 1` -/

/-- **Velocity normalisation.** With `uⁱ = dxⁱ/dq⁰` and the metric relation
    `(dq⁰)² = Q5(dx)` (eq 60), the 5-velocity is unit: `Q5(u) = 1`. -/
theorem velocity_normalization (dx0 dx1 dx2 dx3 dx4 dq0 : ℝ) (hdq0 : dq0 ≠ 0)
    (hmetric : dq0 ^ 2 = Q5 dx0 dx1 dx2 dx3 dx4) :
    Q5 (dx0 / dq0) (dx1 / dq0) (dx2 / dq0) (dx3 / dq0) (dx4 / dq0) = 1 := by
  unfold Q5 at hmetric ⊢
  field_simp
  linarith [hmetric]

/-! ## (3) Mass-density coupling `dγ/dq⁰ = a₀u⁴` -/

/-- **Mass-density velocity (Section C).** Since `γ = a₀·q⁴ = a₀·x⁴`, the rate of
    change of mass density along the trajectory is `dγ/dq⁰ = a₀·(dx⁴/dq⁰) = a₀u⁴`.
    (So `u⁴ = 0 ↔ γ` constant — the condition for constant rest mass, the 4D
    limit.) -/
theorem mass_density_velocity (a0 : ℝ) (X : ℝ → ℝ) (s : ℝ)
    (hX : DifferentiableAt ℝ X s) :
    deriv (fun s => a0 * X s) s = a0 * deriv X s :=
  deriv_const_mul a0 hX

/-! ## (4) Force-density equation of motion (eq 62) and (5) SR reduction (eq 63)

The variational principle δ∫γ√(dq⁰)² = 0 has, at the normalisation g_{ij}uⁱuʲ=1,
the Euler–Lagrange momentum `p_k = γ·(g·u)_k` with `d p_k/dq⁰ = ∂γ/∂x^k`. We
write `gu := (g·u)_k` as a function of `q⁰`. -/

/-- **Force-density equation of motion (eq 62).** Expanding the EL momentum
    derivative `d/dq⁰(γ·(g·u)_k) = ∂γ/∂x^k` by the product rule, with the mass
    coupling `dγ/dq⁰ = a₀u⁴`, gives

        γ·f_k = ∂γ/∂x^k − a₀u⁴·(g·u)_k,

    where `f_k = d(g·u)_k/dq⁰`. The second term is the mass-density-coupling
    force density absent from 4D. -/
theorem euler_lagrange_force (γ gu : ℝ → ℝ) (a0 u4 dgamma_dx : ℝ) (s : ℝ)
    (hγ : DifferentiableAt ℝ γ s) (hgu : DifferentiableAt ℝ gu s)
    (hmass : deriv γ s = a0 * u4)
    (hEL : deriv (fun s => γ s * gu s) s = dgamma_dx) :
    γ s * deriv gu s = dgamma_dx - a0 * u4 * gu s := by
  rw [deriv_fun_mul hγ hgu, hmass] at hEL
  linarith [hEL]

/-- **SR reduction at u⁴ = 0 (eq 63).** When the mass density is constant along
    the trajectory (`u⁴ = 0`, so `dγ/dq⁰ = 0`), the coupling term vanishes and the
    equation of motion reduces to the 4D special-relativity force density
    `γ·f_k = ∂γ/∂x^k` (→ `F̄ᵅ = d/dt[γ¹·ẋᵅ]`). -/
theorem SR_reduction (γ gu : ℝ → ℝ) (a0 dgamma_dx : ℝ) (s : ℝ)
    (hγ : DifferentiableAt ℝ γ s) (hgu : DifferentiableAt ℝ gu s)
    (hmass0 : deriv γ s = 0)
    (hEL : deriv (fun s => γ s * gu s) s = dgamma_dx) :
    γ s * deriv gu s = dgamma_dx := by
  have h := euler_lagrange_force γ gu a0 0 dgamma_dx s hγ hgu (by rw [hmass0]; ring) hEL
  simpa using h

/-! ## (6) Relativistic mass density `γ¹ = γ/√(1−β²)` -/

/-- The two `√`-arguments agree up to the `c²` factor: `1 − v²/c² = (c²−v²)/c²`. -/
theorem rel_mass_arg (c v : ℝ) (hc : c ≠ 0) :
    1 - v ^ 2 / c ^ 2 = (c ^ 2 - v ^ 2) / c ^ 2 := by field_simp

/-- **Relativistic mass density (Section F).** `γ¹ ≡ γ/√(1−β²) = γ·c/√(c²−v²)`
    (β = v/c), the mass-density analogue of relativistic mass, arising when
    converting from q⁰-time to observer time t. -/
theorem rel_mass_density (γ c v : ℝ) (hc : 0 < c) (hv : v ^ 2 ≤ c ^ 2) :
    γ / Real.sqrt (1 - v ^ 2 / c ^ 2) = γ * c / Real.sqrt (c ^ 2 - v ^ 2) := by
  rw [rel_mass_arg c v hc.ne',
      Real.sqrt_div (by linarith : (0 : ℝ) ≤ c ^ 2 - v ^ 2), Real.sqrt_sq hc.le,
      div_div_eq_mul_div]

end DTWilliams.Minkowski5
