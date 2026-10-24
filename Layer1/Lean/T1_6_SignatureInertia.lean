/-
  T1.6 (inertia) — the Lorentzian signature is a BASIS-INDEPENDENT invariant.

  Source: P. E. Williams, LA-8370-MS (LASL, 1980), §II, eqs (16)–(19).
  Supplementary Materials: step [L5], Sec. S2 (code index: Sec. S6).
-/
import Mathlib

namespace DTWilliams.SignatureInertia

variable {n : ℕ}

/-- Spacetime carrier: one (arc-length / absolute-velocity) coordinate ⊕ n spatial. -/
abbrev E (n : ℕ) := ℝ × (Fin n → ℝ)

/-- The T1.6 spacetime form: timelike coordinate enters `+`, the positive-definite
    spatial Hessian `h` enters with a single overall minus sign. -/
def Q (h : (Fin n → ℝ) → ℝ) (w : E n) : ℝ := w.1 ^ 2 - h w.2

private lemma sq_pos_of_ne (x : ℝ) (hx : x ≠ 0) : 0 < x ^ 2 := by
  rw [← sq_abs]; exact pow_pos (abs_pos.mpr hx) 2

/-- **Positive (timelike) direction exists.** `Q(c,0) = c² > 0` for `c ≠ 0`
    (using `h 0 = 0`). Hence the positive index is ≥ 1. -/
theorem time_pos (h : (Fin n → ℝ) → ℝ) (h0 : h 0 = 0) {c : ℝ} (hc : c ≠ 0) :
    0 < Q h ((c, 0) : E n) := by
  have hQ : Q h ((c, 0) : E n) = c ^ 2 := by unfold Q; simp [h0]
  rw [hQ]; exact sq_pos_of_ne c hc

/-- **Every spatial direction is negative.** For positive-definite `h`, every
    nonzero `x` gives `Q(0,x) = −h(x) < 0`: the n-dimensional spatial subspace
    `{t = 0}` is negative-definite. Hence the negative index is ≥ n. -/
theorem space_neg (h : (Fin n → ℝ) → ℝ) (hpos : ∀ x, x ≠ 0 → 0 < h x)
    {x : Fin n → ℝ} (hx : x ≠ 0) : Q h ((0, x) : E n) < 0 := by
  have hQ : Q h ((0, x) : E n) = - h x := by unfold Q; simp
  rw [hQ]; linarith [hpos x hx]

/-- **No two-dimensional subspace is positive-definite** (positive index ≤ 1).
    For any linearly independent pair `u, v`, their span contains a nonzero
    vector `w` with `Q w < 0`. Proof: pick a scalar combination whose time
    component vanishes; it lands in the negative spatial cone. Elementary —
    holds for an ARBITRARY positive-definite `h`. -/
theorem no_two_dim_positive (h : (Fin n → ℝ) → ℝ) (hpos : ∀ x, x ≠ 0 → 0 < h x)
    (u v : E n) (hli : LinearIndependent ℝ ![u, v]) :
    ∃ w ∈ Submodule.span ℝ ({u, v} : Set (E n)), w ≠ 0 ∧ Q h w < 0 := by
  have hu0 : u ≠ 0 := by simpa using hli.ne_zero 0
  by_cases hu1 : u.1 = 0
  · -- u itself already has zero time component
    have hu2 : u.2 ≠ 0 := by
      intro h2
      exact hu0 (Prod.ext_iff.mpr ⟨by simpa using hu1, by simpa using h2⟩)
    refine ⟨u, Submodule.subset_span (by simp), hu0, ?_⟩
    have hQ : Q h u = - h u.2 := by unfold Q; rw [hu1]; ring
    rw [hQ]; linarith [hpos u.2 hu2]
  · -- combination w = v.1 • u + (-u.1) • v has zero time component
    set w : E n := v.1 • u + (-u.1) • v with hwdef
    have hw1 : w.1 = 0 := by
      rw [hwdef]
      simp only [Prod.fst_add, Prod.smul_fst, smul_eq_mul]
      ring
    have hwne : w ≠ 0 := by
      intro hzero
      have hcoef := (LinearIndependent.pair_iff.mp hli) v.1 (-u.1) (by rw [hwdef] at hzero; exact hzero)
      exact hu1 (neg_eq_zero.mp hcoef.2)
    have hw2 : w.2 ≠ 0 := by
      intro h2
      exact hwne (Prod.ext_iff.mpr ⟨by simpa using hw1, by simpa using h2⟩)
    refine ⟨w, ?_, hwne, ?_⟩
    · rw [hwdef]
      exact Submodule.add_mem _
        (Submodule.smul_mem _ _ (Submodule.subset_span (by simp)))
        (Submodule.smul_mem _ _ (Submodule.subset_span (by simp)))
    · have hQ : Q h w = - h w.2 := by unfold Q; rw [hw1]; ring
      rw [hQ]; linarith [hpos w.2 hw2]

/-- **T1.6 signature, basis-independent (Sylvester inertia, positive index).**
    For any positive-definite spatial Hessian `h` (with `h 0 = 0`):
      (i)   a positive timelike direction exists                  (pos index ≥ 1);
      (ii)  no 2-dimensional subspace is positive-definite        (pos index ≤ 1);
      (iii) the n-dimensional spatial subspace is negative-definite (neg index ≥ n).
    Hence the positive index is EXACTLY 1 — Euclidean/multi-timelike signatures
    are impossible — and there is an n-dimensional negative-definite subspace.
    With nondegeneracy this is signature (1, n): Lorentzian. -/
theorem signature_one_n (h : (Fin n → ℝ) → ℝ) (h0 : h 0 = 0)
    (hpos : ∀ x, x ≠ 0 → 0 < h x) :
    (0 < Q h ((1, 0) : E n)) ∧
    (∀ u v : E n, LinearIndependent ℝ ![u, v] →
        ∃ w ∈ Submodule.span ℝ ({u, v} : Set (E n)), w ≠ 0 ∧ Q h w < 0) ∧
    (∀ x : Fin n → ℝ, x ≠ 0 → Q h ((0, x) : E n) < 0) :=
  ⟨time_pos h h0 one_ne_zero,
   fun u v hli => no_two_dim_positive h hpos u v hli,
   fun _ hx => space_neg h hpos hx⟩

/-! ## Negative index ≤ n — completing the (1, n) count -/

/-- The spacetime carrier has dimension `n + 1`. -/
theorem finrank_E : Module.finrank ℝ (E n) = n + 1 := by
  show Module.finrank ℝ (ℝ × (Fin n → ℝ)) = n + 1
  rw [Module.finrank_prod, Module.finrank_self, Module.finrank_fin_fun]
  omega

/-- **Negative index ≤ n.** Any negative-definite subspace has dimension ≤ n.
    Were it `(n+1)`-dimensional it would equal all of `E` (dimension count), hence
    contain the timelike vector `(1,0)` with `Q(1,0) = 1 > 0`, contradicting
    negative-definiteness. -/
theorem neg_index_le_n (h : (Fin n → ℝ) → ℝ) (h0 : h 0 = 0)
    (V : Submodule ℝ (E n)) (hV : ∀ w ∈ V, w ≠ 0 → Q h w < 0) :
    Module.finrank ℝ V ≤ n := by
  by_contra hcon
  rw [not_le] at hcon
  have heq : Module.finrank ℝ V = Module.finrank ℝ (E n) := by
    have hle := Submodule.finrank_le V
    rw [finrank_E] at hle ⊢
    omega
  have htop : V = ⊤ := Submodule.eq_top_of_finrank_eq heq
  have hmem : ((1, 0) : E n) ∈ V := by rw [htop]; exact Submodule.mem_top
  have hne : ((1, 0) : E n) ≠ 0 := by simp [Prod.ext_iff]
  have hneg := hV _ hmem hne
  have hpos := time_pos h h0 (one_ne_zero)
  linarith

/-- **Full Lorentzian signature (1, n), basis-independent.**
    Positive index = 1: a timelike direction exists AND no 2-dim subspace is
    positive-definite. Negative index = n: every nonzero spatial vector is
    negative AND no subspace of dimension > n is negative-definite. Hence the
    signature is exactly (1, n) — Lorentzian — for ANY positive-definite stability
    Hessian. Euclidean and all other signatures are impossible. -/
theorem signature_eq_one_n (h : (Fin n → ℝ) → ℝ) (h0 : h 0 = 0)
    (hpos : ∀ x, x ≠ 0 → 0 < h x) :
    (0 < Q h ((1, 0) : E n)) ∧
    (∀ u v : E n, LinearIndependent ℝ ![u, v] →
        ∃ w ∈ Submodule.span ℝ ({u, v} : Set (E n)), w ≠ 0 ∧ Q h w < 0) ∧
    (∀ x : Fin n → ℝ, x ≠ 0 → Q h ((0, x) : E n) < 0) ∧
    (∀ V : Submodule ℝ (E n), (∀ w ∈ V, w ≠ 0 → Q h w < 0) → Module.finrank ℝ V ≤ n) :=
  ⟨time_pos h h0 one_ne_zero,
   fun u v hli => no_two_dim_positive h hpos u v hli,
   fun _ hx => space_neg h hpos hx,
   fun V hV => neg_index_le_n h h0 V hV⟩

end DTWilliams.SignatureInertia
