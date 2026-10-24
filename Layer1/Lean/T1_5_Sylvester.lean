/-
  T1.5 — Stability Conditions → Positive-Definite Quadratic Form

  Source: P. E. Williams, LA-8370-MS (LASL, 1980), §II, eqs (14)–(15).
  Supplementary Materials: step [L5], Sec. S2 (code index: Sec. S6).
-/

import Mathlib.Algebra.Order.Field.Basic
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Tactic

open Real

/-!
## Completing-the-square identity

The algebraic heart of the Sylvester criterion. The identity
  a · Q(x,y) = (ax + by)² + (ac - b²) · y²
shows that a · Q is a sum of two nonneg terms, each zero only when
(x,y) = (0,0) under the criterion conditions.
-/

private lemma completing_square (a b c x y : ℝ) :
    a * (a * x ^ 2 + 2 * b * x * y + c * y ^ 2) =
    (a * x + b * y) ^ 2 + (a * c - b ^ 2) * y ^ 2 := by ring

/-!
## Strict positivity of squares
-/

private lemma sq_pos_of_ne (x : ℝ) (hx : x ≠ 0) : 0 < x ^ 2 := by
  rw [← sq_abs]
  exact pow_pos (abs_pos.mpr hx) 2

/-!
## Main theorem: Sylvester criterion is sufficient for positive definiteness
-/

/-- **Sylvester criterion (sufficiency).**
    For a 2D symmetric quadratic form Q(x,y) = ax² + 2bxy + cy²:
    if a > 0 and ac - b² > 0, then Q is strictly positive definite,
    i.e. Q(x,y) > 0 whenever (x,y) ≠ (0,0).

    This is the key step in T1.5: the stability quadratic form from
    the Hessian of the energy function is positive definite. -/
theorem sylvester_2d {a b c : ℝ} (ha : 0 < a) (hdet : 0 < a * c - b ^ 2) :
    ∀ x y : ℝ, (x ≠ 0 ∨ y ≠ 0) → 0 < a * x ^ 2 + 2 * b * x * y + c * y ^ 2 := by
  intro x y hne
  -- Step 1: The completing-the-square identity
  have key := completing_square a b c x y
  -- Step 2: Both summands of a·Q are nonneg
  have h1 : 0 ≤ (a * x + b * y) ^ 2 := sq_nonneg _
  have h2 : 0 ≤ (a * c - b ^ 2) * y ^ 2 := mul_nonneg hdet.le (sq_nonneg _)
  -- Step 3: Show a · Q > 0 by case analysis on which of x,y is nonzero
  have haQ : 0 < a * (a * x ^ 2 + 2 * b * x * y + c * y ^ 2) := by
    rw [key]
    rcases hne with hx | hy
    · -- Case x ≠ 0
      rcases eq_or_ne y 0 with rfl | hy0
      · -- Subcase y = 0: sum becomes (ax)² = a²x² > 0  (second term vanishes)
        nlinarith [mul_pos ha ha, sq_pos_of_ne x hx]
      · -- Subcase y ≠ 0: second term (ac-b²)y² > 0
        nlinarith [h1, mul_pos hdet (sq_pos_of_ne y hy0)]
    · -- Case y ≠ 0: second term (ac-b²)y² > 0 immediately
      nlinarith [h1, mul_pos hdet (sq_pos_of_ne y hy)]
  -- Step 4: Extract Q > 0 from a · Q > 0, using a > 0
  -- (assume Q ≤ 0; then a·Q ≤ 0 since a > 0, contradicting haQ)
  nlinarith [haQ, ha]

/-- Cleaner version using `nlinarith` for step 3. -/
theorem sylvester_2d' {a b c : ℝ} (ha : 0 < a) (hdet : 0 < a * c - b ^ 2) :
    ∀ x y : ℝ, (x ≠ 0 ∨ y ≠ 0) → 0 < a * x ^ 2 + 2 * b * x * y + c * y ^ 2 := by
  intro x y hne
  have key : a * (a * x ^ 2 + 2 * b * x * y + c * y ^ 2) =
             (a * x + b * y) ^ 2 + (a * c - b ^ 2) * y ^ 2 := by ring
  rcases hne with hx | hy
  · rcases eq_or_ne y 0 with rfl | hy0
    · -- y = 0: Q = ax² > 0 directly (no need for the identity)
      nlinarith [mul_pos ha ha, sq_pos_of_ne x hx]
    · -- y ≠ 0: use the identity and (ac-b²)y² > 0
      nlinarith [sq_nonneg (a * x + b * y), mul_pos hdet (sq_pos_of_ne y hy0)]
  · -- y ≠ 0: (ac-b²)y² > 0 gives the whole RHS > 0
    nlinarith [sq_nonneg (a * x + b * y), mul_pos hdet (sq_pos_of_ne y hy)]

/-!
## Necessity: positive definite → criterion
-/

/-- **Sylvester criterion (necessity).**
    If ax² + 2bxy + cy² > 0 for all (x,y) ≠ (0,0), then a > 0. -/
lemma sylvester_ne_a {a b c : ℝ}
    (hQ : ∀ x y : ℝ, (x ≠ 0 ∨ y ≠ 0) → 0 < a * x ^ 2 + 2 * b * x * y + c * y ^ 2) :
    0 < a := by
  have := hQ 1 0 (Or.inl one_ne_zero)
  simp at this
  exact this

/-- If ax² + 2bxy + cy² > 0 for all (x,y) ≠ (0,0), then ac - b² > 0. -/
lemma sylvester_ne_det {a b c : ℝ}
    (hQ : ∀ x y : ℝ, (x ≠ 0 ∨ y ≠ 0) → 0 < a * x ^ 2 + 2 * b * x * y + c * y ^ 2) :
    0 < a * c - b ^ 2 := by
  -- Evaluate Q at (b, -a): Q(b, -a) = ab² - 2ab² + a²c = a(ac - b²)
  have h1 : 0 < a := sylvester_ne_a hQ
  have h2 := hQ b (-a) (Or.inr (neg_ne_zero.mpr h1.ne'))
  have heq : a * b ^ 2 + 2 * b * b * (-a) + c * (-a) ^ 2 = a * (a * c - b ^ 2) := by ring
  rw [heq] at h2
  -- a > 0 and a·(ac-b²) > 0 ⇒ ac-b² > 0
  nlinarith [h1, h2]

/-!
## Full equivalence
-/

/-- **Theorem T1.5 (Sylvester criterion, full equivalence).**
    The 2D quadratic form Q(x,y) = ax² + 2bxy + cy² is positive definite
    if and only if a > 0 and ac - b² > 0.

    This is used in the Dynamic Theory to show that the stability Hessian
    h_{ij} = ∂²U/∂qⁱ∂qʲ defines a positive-definite metric on the
    space-entropy manifold, which in turn yields the Minkowski metric
    after arc-length parameterisation (T1.6). -/
theorem T1_5_sylvester_iff (a b c : ℝ) :
    (∀ x y : ℝ, (x ≠ 0 ∨ y ≠ 0) → 0 < a * x ^ 2 + 2 * b * x * y + c * y ^ 2) ↔
    (0 < a ∧ 0 < a * c - b ^ 2) :=
  ⟨fun hQ => ⟨sylvester_ne_a hQ, sylvester_ne_det hQ⟩,
   fun ⟨ha, hdet⟩ => sylvester_2d ha hdet⟩

/-!
## Connection to the Dynamic Theory (T1.5)

The stability quadratic form arising from the Hessian of U(q,S) is
  Q = U_qq (δq)² + 2 U_qS (δq)(δS) + U_SS (δS)²
(Williams equation 15). Positive definiteness of this form requires
exactly the Sylvester criterion with a = U_qq, b = U_qS, c = U_SS.
This Lean theorem gives a machine-checked proof of the sufficiency
direction used in T1.5.
-/

/-- The stability quadratic from T1.5 is positive definite under
    Williams' stability conditions:  U_qq > 0  and  U_qq·U_SS - U_qS² > 0. -/
theorem T1_5_stability (U_qq U_qS U_SS : ℝ)
    (h1 : 0 < U_qq)
    (h2 : 0 < U_qq * U_SS - U_qS ^ 2)
    (δq δS : ℝ)
    (hne : δq ≠ 0 ∨ δS ≠ 0) :
    0 < U_qq * δq ^ 2 + 2 * U_qS * δq * δS + U_SS * δS ^ 2 :=
  sylvester_2d h1 h2 δq δS hne
