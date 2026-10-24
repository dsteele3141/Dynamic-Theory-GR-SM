/-
  T1.6 (necessity) — the Lorentzian signature is FORCED, not chosen

  Source: P. E. Williams, LA-8370-MS (LASL, 1980), §II, eqs (16)–(19).
  Supplementary Materials: step [L5], Sec. S2 (code index: Sec. S6).
-/

import Mathlib.Tactic

namespace DTWilliams.SignatureNecessity

variable {n : ℕ}

/-- The spacetime form produced by the T1.6 construction: a single timelike
    (arc-length / absolute-velocity) coordinate `t` entering with `+`, and the
    positive-definite spatial stability Hessian `h` entering with a single overall
    sign flip. `h : (Fin n → ℝ) → ℝ` is the spatial quadratic form. -/
def Q (h : (Fin n → ℝ) → ℝ) (t : ℝ) (x : Fin n → ℝ) : ℝ := t ^ 2 - h x

variable (h : (Fin n → ℝ) → ℝ)

/-- **Timelike direction (one `+`).** The arc-length / absolute-velocity
    coordinate is positive: `Q(1,0) = 1 > 0`. -/
theorem timelike (h0 : h 0 = 0) : 0 < Q h 1 0 := by
  unfold Q; rw [h0]; norm_num

/-- **All spatial directions are spacelike.** For positive-definite `h`, every
    nonzero spatial vector gives `Q(0,x) = −h(x) < 0`: the whole n-dimensional
    spatial subspace is negative-definite. This is the heart of necessity — no
    spatial direction can become timelike, for ANY positive-definite Hessian. -/
theorem spacelike (hpos : ∀ x, x ≠ 0 → 0 < h x) {x : Fin n → ℝ} (hx : x ≠ 0) :
    Q h 0 x < 0 := by
  have hp : 0 < h x := hpos x hx
  have z : (0 : ℝ) ^ 2 = 0 := by norm_num
  unfold Q; rw [z]; linarith

/-- **Euclidean signature is impossible.** Since `Q` is negative on every spatial
    direction (`spacelike`), it is NOT positive-definite — the construction can
    never yield a Euclidean `(+,+,…,+)` signature. -/
theorem not_positive_definite (hpos : ∀ x, x ≠ 0 → 0 < h x)
    {x : Fin n → ℝ} (hx : x ≠ 0) : ¬ (0 < Q h 0 x) := by
  have := spacelike h hpos hx; linarith

/-- **Indefinite, with the Lorentzian split.** Positive on the time axis, negative
    on (all of) space: the form has exactly one `+` direction (the arc-length
    coordinate) and an n-dimensional negative-definite spatial subspace. With
    Sylvester's law of inertia this fixes the signature at (1, n) — Lorentzian. -/
theorem lorentzian_split (hpos : ∀ x, x ≠ 0 → 0 < h x) (h0 : h 0 = 0)
    {x : Fin n → ℝ} (hx : x ≠ 0) :
    0 < Q h 1 0 ∧ Q h 0 x < 0 :=
  ⟨timelike h h0, spacelike h hpos hx⟩

end DTWilliams.SignatureNecessity
