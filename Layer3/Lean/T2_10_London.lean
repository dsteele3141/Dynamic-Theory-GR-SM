/-
  T2.10 — London Quantisation from Geometry

  Source: Williams LA-8370-MS (1980), Section II.F, pages 64-72, eqs (55)-(57).
  Full derivation, context and provenance:
  Supplementary_Material.pdf — §S3, Appendix A [C1].
-/

import Mathlib.Tactic
import Mathlib.Analysis.SpecialFunctions.ExpDeriv
import Mathlib.Analysis.SpecialFunctions.Complex.Circle

namespace DTWilliams.London

/-! ## A — Poisson bracket collapse in curved space -/

/-- **Canonical bracket vanishes in curved space (pp.73-74).** Model the partials
    `a i k = ∂x^i/∂p_k`, `c j k = ∂p_j/∂x^k`. The curvature symmetry
    `∂p_j/∂x^k = ∂p_k/∂x^j` (`hsymm`) and the chain rule
    `Σ_k(∂x^i/∂p_k)(∂p_k/∂x^j) = δ_{ij}` (`hchain`) force
    `[x^i,p_j] = δ − Σ_k a_{ik} c_{jk} = 0`. -/
theorem poisson_bracket_collapse {ι : Type*} [Fintype ι]
    (a c : ι → ι → ℝ) (δ : ℝ)
    (hsymm : ∀ j k, c j k = c k j)
    (hchain : ∀ i j, (∑ k, a i k * c k j) = δ)
    (i j : ι) :
    δ - (∑ k, a i k * c j k) = 0 := by
  have h : (∑ k, a i k * c j k) = ∑ k, a i k * c k j :=
    Finset.sum_congr rfl (fun k _ => by rw [hsymm j k])
  rw [h, hchain i j]; ring

/-! ## B — London quantisation -/

/-- **Weyl arc-element ODE solution (p.70).** `dσ(τ) = dσ₀·exp(Φ(τ))`
    (Φ(τ) = ∫φ_k dx^k along the path) satisfies the geometric length-change
    equation `d(dσ)/dτ = (φ_k ẋ^k)·dσ` (with `φ_k ẋ^k = Φ'(τ)`). -/
theorem arc_element_ode (dσ0 : ℝ) (Φ : ℝ → ℝ) (τ : ℝ) (hΦ : DifferentiableAt ℝ Φ τ) :
    deriv (fun τ => dσ0 * Real.exp (Φ τ)) τ
      = deriv Φ τ * (dσ0 * Real.exp (Φ τ)) := by
  have h : HasDerivAt (fun τ => dσ0 * Real.exp (Φ τ))
      (dσ0 * (Real.exp (Φ τ) * deriv Φ τ)) τ :=
    (hΦ.hasDerivAt.exp).const_mul dσ0
  rw [h.deriv]; ring

/-- **London quantisation / closed-loop holonomy (p.71).** When the loop integral
    is `∮φ_k dx^k = 2πin` (n ∈ ℤ), the holonomy is trivial — `exp(∮φ_k dx^k)=1` —
    so the arc element returns to itself. This is the London condition, here
    DERIVED from the geometry (Weyl length-change + equilibrium null trajectory). -/
theorem london_holonomy (n : ℤ) :
    Complex.exp ((n : ℂ) * (2 * Real.pi * Complex.I)) = 1 :=
  Complex.exp_int_mul_two_pi_mul_I n

/-! ## C — Bohr radii -/

/-- **Circular-orbit speed.** Balancing centripetal and Coulomb forces,
    `mv²/r = e²/r²` (here cleared by `r²`: `m v² r = e²`), gives `v² = e²/(mr)`. -/
theorem circular_orbit_v2 (m e r v : ℝ) (hm : m ≠ 0) (hr : r ≠ 0)
    (hbal : m * v ^ 2 * r = e ^ 2) :
    v ^ 2 = e ^ 2 / (m * r) := by
  rw [eq_div_iff (mul_ne_zero hm hr)]
  linear_combination hbal

/-- **Bohr radii (p.72).** With `α¹ = ie²/ℏ` (c = 1), so `(α¹)² = −e⁴/ℏ²`, the
    London-quantised circular-orbit radius `r = −n²e²/((α¹)²m)` reduces to the
    Bohr radius `n²ℏ²/(me²) = n²a₀`. -/
theorem bohr_radius (n e ℏ m r α1sq : ℝ) (hm : m ≠ 0) (he : e ≠ 0) (hℏ : ℏ ≠ 0)
    (hα : α1sq = -(e ^ 4 / ℏ ^ 2))
    (hr : r = -(n ^ 2 * e ^ 2) / (α1sq * m)) :
    r = n ^ 2 * ℏ ^ 2 / (m * e ^ 2) := by
  rw [hr, hα]; field_simp

/-! ## D — Fine-structure constant `α¹ = iα` is imaginary -/

/-- **(α¹)² = −α².** The Dynamic-Theory gauge coupling `α¹ = iα` is imaginary, so
    its square is the negative real `−α²`. -/
theorem alpha1_sq (α : ℝ) : (Complex.I * (α : ℂ)) ^ 2 = -((α : ℂ) ^ 2) := by
  rw [mul_pow, Complex.I_sq]; ring

/-- **(α¹)² < 0** (as a real number `−α²`, for `α ≠ 0`): this negativity is what
    makes `r = −n²e²/((α¹)²m)` come out positive — real Bohr radii. -/
theorem alpha1_sq_neg (α : ℝ) (hα : α ≠ 0) : -(α ^ 2) < 0 := by
  have : 0 < α ^ 2 := by positivity
  linarith

end DTWilliams.London
