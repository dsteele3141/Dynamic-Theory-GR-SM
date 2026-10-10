/-
  The 5D Dirac ansatz — first-order factorisation of the 5D mass-shell relation
  ("Axiom X" and the "fourth postulate" are historical internal nicknames for it;
  the object is an ADOPTED ANSATZ, not one of DT's postulates — see status note.)

  Source: P. E. Williams, LA-8370-MS (LASL, 1980), Section III.C, eqs (105)-(112).
  Full derivation, context and provenance:
  Supplementary_Material.pdf — §S2, Appendix A [C6].

  Status (T2.7): Williams derives the quadratic 5D mass-shell relation (eq 105)
  from his three laws, but ADOPTS its first-order (Dirac) factorisation as an
  ansatz ("may be taken to be"; "must be accepted as an independent fundamental
  assumption"), exactly as Dirac did in 1928. He hoped it might one day be derived
  from the three laws; T2.7 shows it cannot be — the anticommuting factorisation
  has no realisation in any commutative ring of characteristic ≠ 2. This is a
  commutative-ring obstruction, NOT a full model-theoretic independence proof from
  the three laws (which are not formalised as axioms). So it is an adopted ansatz,
  genuinely independent of the three laws — NOT a "fourth postulate". DT has three
  postulates, not four.
-/

import Mathlib.Tactic
import Mathlib.Data.Matrix.Basic
import Mathlib.Basic.Complex.Basic
import Mathlib.LinearAlgebra.Matrix.Kronecker
import Mathlib.LinearAlgebra.Matrix.Reindex

-- `FourthPostulate` is a HISTORICAL NICKNAME retained as a stable code identifier;
-- the object is the 5D Dirac ansatz (see header), NOT one of DT's postulates.
namespace DTWilliams.FourthPostulate

/-! ## The 5D metric  (signature +,−,−,−,−)

Williams' eq (112) metric:  g^{jk} = diag(1,−1,−1,−1,−1) on indices j = 0..4,
where x⁰ = ct (time) and x¹..x⁴ = q¹,q²,q³, γ/a₀ (space + mass density). -/

/-- The 5D Dynamic-Theory metric `g^{ji}` valued in an arbitrary ring `R`. -/
def gMetric (R : Type*) [Ring R] (i j : Fin 5) : R :=
  if i = j then (if i = 0 then (1 : R) else -1) else 0

/-! ## Dirac-ansatz data — the first-order factorisation hypothesis (eq 108)

`DiracData R` bundles the four `α̃`-operators and the `β̃`-operator together
with EXACTLY the relations eq (108) demands. Unpacking `(α̃·p)² = p²` for all
momenta `p` gives `α̃_i² = 1` and `α̃_iα̃_j + α̃_jα̃_i = 0` (i ≠ j); together
with `β̃² = 1` and `α̃_iβ̃ + β̃α̃_i = 0` these are the entire content of the
Dirac ansatz. -/
structure DiracData (R : Type*) [Ring R] where
  /-- α̃₁,…,α̃₄  (0-indexed as `α 0`..`α 3`). -/
  α : Fin 4 → R
  /-- β̃. -/
  β : R
  /-- α̃_i² = 1  (from `(α̃·p)² = p²`, diagonal part). -/
  α_sq : ∀ i, α i * α i = 1
  /-- α̃_iα̃_j + α̃_jα̃_i = 0 for i ≠ j  (off-diagonal part of `(α̃·p)² = p²`). -/
  α_anticomm : ∀ i j, i ≠ j → α i * α j + α j * α i = 0
  /-- β̃² = 1. -/
  β_sq : β * β = 1
  /-- α̃_iβ̃ + β̃α̃_i = 0. -/
  αβ_anticomm : ∀ i, α i * β + β * α i = 0

variable {R : Type*} [Ring R]

/-- The 5D gamma operators, eq (111):  γ⁰ ≡ β̃,  γ^μ ≡ −β̃α̃_μ  (μ = 1..4). -/
def DiracData.γ (D : DiracData R) (k : Fin 5) : R :=
  if h : k = 0 then D.β else - (D.β * D.α (k.pred h))

@[simp] lemma DiracData.γ_zero (D : DiracData R) : D.γ 0 = D.β := by
  simp [DiracData.γ]

lemma DiracData.γ_succ (D : DiracData R) (k : Fin 5) (h : k ≠ 0) :
    D.γ k = - (D.β * D.α (k.pred h)) := by
  simp [DiracData.γ, h]

/-! ### Algebraic helper lemmas

These repackage the `DiracData` relations into the forms the case proofs need. -/

/-- The metric `g^{ji}` is symmetric. -/
lemma gMetric_symm (i j : Fin 5) : gMetric R i j = gMetric R j i := by
  unfold gMetric
  rcases eq_or_ne i j with rfl | h
  · rfl
  · simp [h, Ne.symm h]

/-- `g^{00} = 1`. -/
lemma gMetric_zero_zero : gMetric R (0 : Fin 5) 0 = 1 := by simp [gMetric]

/-- `α̃_j β̃ = −β̃ α̃_j` (rearranged `αβ_anticomm`). -/
lemma DiracData.alpha_beta (D : DiracData R) (j : Fin 4) :
    D.α j * D.β = - (D.β * D.α j) :=
  eq_neg_of_add_eq_zero_left (D.αβ_anticomm j)

/-- `β̃ α̃_j β̃ = −α̃_j`  (conjugation by β̃ flips the sign), from `alpha_beta` + `β_sq`. -/
lemma DiracData.beta_conj (D : DiracData R) (j : Fin 4) :
    D.β * D.α j * D.β = - D.α j := by
  calc D.β * D.α j * D.β
      = D.β * (D.α j * D.β) := by rw [mul_assoc]
    _ = D.β * -(D.β * D.α j) := by rw [D.alpha_beta]
    _ = -(D.β * D.β * D.α j) := by rw [mul_neg, ← mul_assoc]
    _ = -(1 * D.α j) := by rw [D.β_sq]
    _ = - D.α j := by rw [one_mul]

/-- Product of two spatial gammas: `γ^μ γ^ν = −(α̃_μ α̃_ν)`. -/
lemma DiracData.gamma_spatial_mul (D : DiracData R) (i k : Fin 5)
    (hi : i ≠ 0) (hk : k ≠ 0) :
    D.γ i * D.γ k = -(D.α (i.pred hi) * D.α (k.pred hk)) := by
  rw [D.γ_succ i hi, D.γ_succ k hk, neg_mul_neg]
  calc D.β * D.α (i.pred hi) * (D.β * D.α (k.pred hk))
      = D.β * D.α (i.pred hi) * D.β * D.α (k.pred hk) := by rw [← mul_assoc]
    _ = - D.α (i.pred hi) * D.α (k.pred hk) := by rw [D.beta_conj]
    _ = -(D.α (i.pred hi) * D.α (k.pred hk)) := by rw [neg_mul]

/-- Square of a spatial gamma: `(γ^μ)² = −1`. -/
lemma DiracData.gamma_sq_spatial (D : DiracData R) (i : Fin 5) (hi : i ≠ 0) :
    D.γ i * D.γ i = -1 := by
  rw [D.gamma_spatial_mul i i hi hi, D.α_sq (i.pred hi)]

/-! ## Theorem 1 — Clifford algebra from the Dirac ansatz

The four cases of `{γ^j,γ^i} = 2g^{ji}`. Each is pure substitution using the
`DiracData` hypotheses; the derivations are recorded in full above each case and
are proved mechanically below. -/

/-- **Case j = i = 0.**  {γ⁰,γ⁰} = 2β̃² = 2 = 2·g⁰⁰.
    Derivation: γ⁰ = β̃, so γ⁰γ⁰ + γ⁰γ⁰ = 2β̃² = 2·1 (by `β_sq`). -/
theorem clifford_time_time (D : DiracData R) :
    D.γ 0 * D.γ 0 + D.γ 0 * D.γ 0 = 2 * gMetric R 0 0 := by
  rw [gMetric_zero_zero, DiracData.γ_zero, ← two_mul, D.β_sq]

/-- **Case one time, one space** (μ ≥ 1).  {γ^μ,γ⁰} = 0 = 2·g^{μ0}.
    Derivation:  γ^μγ⁰ + γ⁰γ^μ = (−β̃α̃_μ)β̃ + β̃(−β̃α̃_μ)
      = −β̃α̃_μβ̃ − β̃²α̃_μ
      = −β̃(−β̃α̃_μ) − α̃_μ            (using α̃_μβ̃ = −β̃α̃_μ, `αβ_anticomm`)
      = β̃²α̃_μ − α̃_μ = α̃_μ − α̃_μ = 0 (using `β_sq`). -/
theorem clifford_time_space (D : DiracData R) (k : Fin 5) (hk : k ≠ 0) :
    D.γ k * D.γ 0 + D.γ 0 * D.γ k = 2 * gMetric R k 0 := by
  have hg : gMetric R k (0 : Fin 5) = 0 := by simp [gMetric, hk]
  rw [hg, mul_zero, DiracData.γ_zero, D.γ_succ k hk,
      neg_mul, mul_neg, D.beta_conj, ← mul_assoc, D.β_sq, one_mul, neg_neg,
      add_neg_cancel]

/-- **Case both space**, μ,ν ≥ 1.  {γ^μ,γ^ν} = −2δ_{μν} = 2g^{μν}.
    Derivation:  γ^μγ^ν + γ^νγ^μ = β̃α̃_μβ̃α̃_ν + β̃α̃_νβ̃α̃_μ.
      Using α̃β̃ = −β̃α̃ and β̃² = 1:  β̃α̃_μβ̃ = −α̃_μ, so
      β̃α̃_μβ̃α̃_ν = −α̃_μα̃_ν  and  β̃α̃_νβ̃α̃_μ = −α̃_να̃_μ, giving
      −(α̃_μα̃_ν + α̃_να̃_μ) = −2δ_{μν}      (by `α_sq` if μ=ν, `α_anticomm` else)
      = 2g^{μν}. -/
theorem clifford_space_space (D : DiracData R) (i k : Fin 5)
    (hi : i ≠ 0) (hk : k ≠ 0) :
    D.γ i * D.γ k + D.γ k * D.γ i = 2 * gMetric R i k := by
  rcases eq_or_ne i k with rfl | hik
  · -- i = k:  (γ^μ)² twice  →  −1 + −1 = −2 = 2·(−1) = 2·g^{μμ}
    have hg : gMetric R i i = -1 := by simp [gMetric, hi]
    rw [hg]
    simp only [D.gamma_sq_spatial i hi]
    noncomm_ring
  · -- i ≠ k:  −(α_m α_n) + −(α_n α_m) = −(α_m α_n + α_n α_m) = 0 = 2·g^{μν}
    have hmn : i.pred hi ≠ k.pred hk := fun h => hik (Fin.pred_inj.mp h)
    have hg : gMetric R i k = 0 := by simp [gMetric, hik]
    rw [hg, mul_zero, D.gamma_spatial_mul i k hi hk, D.gamma_spatial_mul k i hk hi,
        ← neg_add, D.α_anticomm (i.pred hi) (k.pred hk) hmn, neg_zero]

/-- **Theorem 1 (Clifford algebra).**  Given the first-order factorisation data
    (the Dirac ansatz), the gamma operators of eq (111) satisfy the
    5D Clifford relation of eq (112):

        {γ^j, γ^i} = 2 g^{ji}.

    This is the representation-independent core of T2.7 — every realisation of the
    Dirac ansatz (in particular the 8×8 matrices of eq 109) must obey it. -/
theorem gamma_clifford (D : DiracData R) (i j : Fin 5) :
    D.γ i * D.γ j + D.γ j * D.γ i = 2 * gMetric R i j := by
  rcases eq_or_ne i 0 with rfl | hi
  · rcases eq_or_ne j 0 with rfl | hj
    · exact clifford_time_time D
    · -- {γ⁰,γ^j}: symmetric to clifford_time_space (swap order + metric symmetry)
      rw [add_comm, clifford_time_space D j hj, gMetric_symm (R := R) j 0]
  · rcases eq_or_ne j 0 with rfl | hj
    · exact clifford_time_space D i hi
    · exact clifford_space_space D i j hi hj

/-! ## Mass-shell consistency — the link back to postulates 1-3 (eq 105)

Postulates 1-3 give the quadratic mass-shell relation `p_j p^j = γ²c²`
(eq 105, natural units → RHS = 1). The Dirac ansatz demands the first-order
operator `h = -(α̃·p̄ + β̃)` square to it. These two statements record that the
factorisation actually reproduces eq (105):

    (α̃·p̄)² = p²          and          h² = p² + 1.

`p` are the (central) momentum scalars: they commute with the α̃, β̃ and with
each other. -/

/-- `(α̃·p̄)² = Σ p_i²`  — the off-diagonal `α̃`-cross terms cancel by
    `α_anticomm`, the diagonal terms collapse by `α_sq`. (eq 108, first relation) -/
theorem alpha_dot_p_sq (D : DiracData R) (p : Fin 4 → R)
    (hpα : ∀ i j, p i * D.α j = D.α j * p i)
    (hpp : ∀ i j, p i * p j = p j * p i) :
    (∑ i, D.α i * p i) * (∑ i, D.α i * p i) = ∑ i, p i * p i := by
  -- Move each `p_i` past the following `α_j`:  (α_i p_i)(α_j p_j) = α_i α_j (p_i p_j).
  have αp : ∀ i j : Fin 4,
      (D.α i * p i) * (D.α j * p j) = D.α i * D.α j * (p i * p j) := by
    intro i j
    rw [mul_assoc, ← mul_assoc (p i), hpα i j, mul_assoc (D.α j), ← mul_assoc (D.α i)]
  -- Off-diagonal α-products are antisymmetric:  α_i α_j = −(α_j α_i)  for i ≠ j.
  have key : ∀ i j : Fin 4, i ≠ j → D.α i * D.α j = -(D.α j * D.α i) :=
    fun i j hij => eq_neg_of_add_eq_zero_left (D.α_anticomm i j hij)
  -- Expand the product of the two 4-term sums into 16 atoms α_i α_j (p_i p_j).
  simp only [Fin.sum_univ_four, add_mul, mul_add, αp]
  -- Diagonal: α_i² = 1, so those four atoms become p_i p_i.
  simp only [D.α_sq, one_mul]
  -- Lower-triangle atoms → negatives of their upper-triangle partners (α then p reorder).
  rw [key 1 0 (by decide), key 2 0 (by decide), key 3 0 (by decide),
      key 2 1 (by decide), key 3 1 (by decide), key 3 2 (by decide)]
  simp only [neg_mul]
  rw [hpp 1 0, hpp 2 0, hpp 3 0, hpp 2 1, hpp 3 1, hpp 3 2]
  -- Each ±pair cancels; the four diagonal terms remain.
  abel

/-- **Mass-shell recovery (eq 105).**  The first-order Hamiltonian
    `h = -(α̃·p̄ + β̃)` squares to the derived mass-shell operator `p² + 1`,
    so its solutions also solve the quadratic relation postulates 1-3 produced.
    The `α̃β̃`-cross terms vanish by `αβ_anticomm`; `β̃²` collapses by `β_sq`. -/
theorem mass_shell (D : DiracData R) (p : Fin 4 → R)
    (hpα : ∀ i j, p i * D.α j = D.α j * p i)
    (hpβ : ∀ i, p i * D.β = D.β * p i)
    (hpp : ∀ i j, p i * p j = p j * p i) :
    (-(∑ i, D.α i * p i) - D.β) * (-(∑ i, D.α i * p i) - D.β)
      = (∑ i, p i * p i) + 1 := by
  -- The β̃-cross terms vanish:  S β̃ + β̃ S = Σ (α_i β̃ + β̃ α_i) p_i = 0.
  have hcross : (∑ i, D.α i * p i) * D.β + D.β * (∑ i, D.α i * p i) = 0 := by
    rw [Finset.sum_mul, Finset.mul_sum, ← Finset.sum_add_distrib]
    apply Finset.sum_eq_zero
    intro i _
    rw [mul_assoc, hpβ i, ← mul_assoc, ← mul_assoc, ← add_mul, D.αβ_anticomm i, zero_mul]
  -- h = −(S + β̃), so h² = (S + β̃)².
  have hneg : -(∑ i, D.α i * p i) - D.β = -((∑ i, D.α i * p i) + D.β) := by abel
  rw [hneg, neg_mul_neg, add_mul, mul_add, mul_add,
      alpha_dot_p_sq D p hpα hpp, D.β_sq]
  -- Regroup the four expanded terms and kill the β̃-cross via hcross.
  have e : ((∑ i, p i * p i) + (∑ i, D.α i * p i) * D.β)
            + (D.β * (∑ i, D.α i * p i) + 1)
         = (∑ i, p i * p i)
            + ((∑ i, D.α i * p i) * D.β + D.β * (∑ i, D.α i * p i)) + 1 := by abel
  rw [e, hcross, add_zero]

/-! ## Theorem 2 — Independence: the Dirac ansatz requires noncommutativity

Postulates 1-3 deliver the mass-shell relation (eq 105) as a relation among
*commuting* scalar quantities — the momenta `p_i`, the mass density `γ`, and `c`
all commute. The Dirac ansatz instead demands operators `α̃, β̃` that *anticommute*
(eq 108). The next results show these anticommutators cannot be realised in any
commutative ring of characteristic ≠ 2: a `DiracData` over a commutative ring
forces `2 = 0`.

Interpretation. The Dirac ansatz is strictly stronger than the scalar relation it
factors: it injects noncommutativity — hence spin — that is absent from, and not
forced by, the commutative content of postulates 1-3. The second-order
Klein–Gordon operator is the commutative-symbol alternative consistent with the
same dispersion relation but carrying no Clifford/spinor structure. So the
first-order factorisation is genuine extra content — an adopted ansatz genuinely
independent of postulates 1-3, not a theorem of them, unless 1-3 were to inject
noncommutativity by some other means (which their thermodynamic/scalar origin
does not). This is the precise, machine-checked sense in which the "square-root"
choice is independent — a commutative-ring obstruction, NOT a full model-theoretic
independence proof from the three laws (which are not formalised as axioms). -/

/-- **Independence core.**  In any commutative ring, a realisation of the Dirac
    ansatz forces `2 = 0`: anticommuting square roots of `1` collapse the characteristic.
    (From `a² = b² = 1` and `ab + ba = 0` with `a,b` commuting: `2 = 2a²b² =
    (ab+ba)·(ab) = 0`.) -/
theorem diracData_comm_char_two {R : Type*} [CommRing R] (D : DiracData R) :
    (2 : R) = 0 := by
  have h := D.α_anticomm 0 1 (by decide)
  have e1 := D.α_sq 0
  have e2 := D.α_sq 1
  linear_combination (D.α 0 * D.α 1) * h - 2 * (D.α 0 * D.α 0) * e2 - 2 * e1

/-- **No commutative realisation.**  Over a commutative ring with `2 ≠ 0` (e.g.
    ℝ, ℚ, ℂ), the Dirac ansatz is unrealisable: the first-order factorisation
    cannot live in the commutative world where the mass-shell relation sits. -/
theorem no_diracData_of_comm {R : Type*} [CommRing R] (h2 : (2 : R) ≠ 0) :
    IsEmpty (DiracData R) :=
  ⟨fun D => h2 (diracData_comm_char_two D)⟩

/-- Concretely: the Dirac ansatz is unrealisable over the reals, while the mass-shell
    relation (a statement about real scalars) of course holds there — pinpointing
    the extra, independent content of the ansatz. -/
example : IsEmpty (DiracData ℝ) := no_diracData_of_comm (by norm_num)

/-! ## Consistency of the Dirac ansatz — a representation exists

The Dirac ansatz is non-vacuous: the 8×8 matrices of eq (109) form a `DiracData` over the
complex 8×8 matrix ring. Proving it here amounts to
constructing that `DiracData` witness and discharging the six field obligations
by matrix computation (`decide`/`Matrix.ext` + `norm_num`). -/

/-- Transport a `DiracData` along a ring homomorphism `f : M →+* N`.  Each
    defining relation pushes forward since `f` preserves `*`, `+`, `0`, `1`. -/
def DiracData.map {M N : Type*} [Ring M] [Ring N] (f : M →+* N) (D : DiracData M) :
    DiracData N where
  α i := f (D.α i)
  β := f D.β
  α_sq i := by rw [← map_mul, D.α_sq, map_one]
  α_anticomm i j h := by rw [← map_mul, ← map_mul, ← map_add, D.α_anticomm i j h, map_zero]
  β_sq := by rw [← map_mul, D.β_sq, map_one]
  αβ_anticomm i := by rw [← map_mul, ← map_mul, ← map_add, D.αβ_anticomm i, map_zero]

open scoped Kronecker
open Matrix

namespace Pauli

-- These are concrete 2×2 matrix computations; the exact simp argument list is
-- immaterial, so silence the (cosmetic) unused-simp-arg linter here.
set_option linter.unusedSimpArgs false

/-- 2×2 Pauli matrices over ℂ. -/
def sx : Matrix (Fin 2) (Fin 2) ℂ := !![0, 1; 1, 0]
def sy : Matrix (Fin 2) (Fin 2) ℂ := !![0, -Complex.I; Complex.I, 0]
def sz : Matrix (Fin 2) (Fin 2) ℂ := !![1, 0; 0, -1]

lemma sx_sq : sx * sx = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;>
    simp [sx, Matrix.mul_apply, Fin.sum_univ_two, Matrix.one_apply]
lemma sy_sq : sy * sy = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;>
    simp [sy, Matrix.mul_apply, Fin.sum_univ_two, Matrix.one_apply, Complex.I_mul_I]
lemma sz_sq : sz * sz = 1 := by
  ext i j; fin_cases i <;> fin_cases j <;>
    simp [sz, Matrix.mul_apply, Fin.sum_univ_two, Matrix.one_apply]
lemma sxy : sx * sy + sy * sx = 0 := by
  ext i j; fin_cases i <;> fin_cases j <;>
    simp [sx, sy, Matrix.mul_apply]
lemma sxz : sx * sz + sz * sx = 0 := by
  ext i j; fin_cases i <;> fin_cases j <;>
    simp [sx, sz]
lemma syz : sy * sz + sz * sy = 0 := by
  ext i j; fin_cases i <;> fin_cases j <;>
    simp [sy, sz]

end Pauli

/-- A Kronecker product of two involutions is an involution. -/
lemma kron_invol {A B : Matrix (Fin 2) (Fin 2) ℂ} (hA : A * A = 1) (hB : B * B = 1) :
    (A ⊗ₖ B) * (A ⊗ₖ B) = 1 := by
  rw [← mul_kronecker_mul, hA, hB, one_kronecker_one]

/-- Anticommutator when the *first* factors anticommute and the *second* commute. -/
lemma kron_anti1 {A B C D : Matrix (Fin 2) (Fin 2) ℂ}
    (hAC : A * C + C * A = 0) (hBD : D * B = B * D) :
    (A ⊗ₖ B) * (C ⊗ₖ D) + (C ⊗ₖ D) * (A ⊗ₖ B) = 0 := by
  rw [← mul_kronecker_mul, ← mul_kronecker_mul, hBD, ← add_kronecker, hAC, zero_kronecker]

/-- Anticommutator when the *first* factors commute and the *second* anticommute. -/
lemma kron_anti2 {A B C D : Matrix (Fin 2) (Fin 2) ℂ}
    (hAC : A * C = C * A) (hBD : B * D + D * B = 0) :
    (A ⊗ₖ B) * (C ⊗ₖ D) + (C ⊗ₖ D) * (A ⊗ₖ B) = 0 := by
  rw [← mul_kronecker_mul, ← mul_kronecker_mul, hAC, ← kronecker_add, hBD, kronecker_zero]

/-- A 4×4 realisation of the Dirac ansatz (the five anticommuting involutions of eq 108),
    built from Pauli matrices on two qubits:
    α̃₁ = σx⊗1, α̃₂ = σy⊗1, α̃₃ = σz⊗σx, α̃₄ = σz⊗σy, β̃ = σz⊗σz. -/
def D4 : DiracData (Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ) where
  α := ![Pauli.sx ⊗ₖ 1, Pauli.sy ⊗ₖ 1, Pauli.sz ⊗ₖ Pauli.sx, Pauli.sz ⊗ₖ Pauli.sy]
  β := Pauli.sz ⊗ₖ Pauli.sz
  α_sq i := by
    fin_cases i
    · exact kron_invol Pauli.sx_sq (one_mul 1)
    · exact kron_invol Pauli.sy_sq (one_mul 1)
    · exact kron_invol Pauli.sz_sq Pauli.sx_sq
    · exact kron_invol Pauli.sz_sq Pauli.sy_sq
  α_anticomm i j h := by
    fin_cases i <;> fin_cases j <;>
      first
        | exact absurd rfl h
        | exact kron_anti1 Pauli.sxy rfl
        | (rw [add_comm]; exact kron_anti1 Pauli.sxy rfl)
        | exact kron_anti1 Pauli.sxz (by rw [mul_one, one_mul])
        | (rw [add_comm]; exact kron_anti1 Pauli.sxz (by rw [mul_one, one_mul]))
        | exact kron_anti1 Pauli.syz (by rw [mul_one, one_mul])
        | (rw [add_comm]; exact kron_anti1 Pauli.syz (by rw [mul_one, one_mul]))
        | exact kron_anti2 rfl Pauli.sxy
        | (rw [add_comm]; exact kron_anti2 rfl Pauli.sxy)
  β_sq := kron_invol Pauli.sz_sq Pauli.sz_sq
  αβ_anticomm i := by
    fin_cases i <;>
      first
        | exact kron_anti1 Pauli.sxz (by rw [mul_one, one_mul])
        | exact kron_anti1 Pauli.syz (by rw [mul_one, one_mul])
        | exact kron_anti2 rfl Pauli.sxz
        | exact kron_anti2 rfl Pauli.syz

/-- Doubling ring homomorphism `g ↦ 1₂ ⊗ₖ g` (block-diagonal embedding
    M₄(ℂ) ↪ M₈(ℂ)), used to reach Williams' reducible 8×8 form of eq (109). -/
def doubleHom :
    Matrix (Fin 2 × Fin 2) (Fin 2 × Fin 2) ℂ →+*
      Matrix (Fin 2 × (Fin 2 × Fin 2)) (Fin 2 × (Fin 2 × Fin 2)) ℂ where
  toFun g := (1 : Matrix (Fin 2) (Fin 2) ℂ) ⊗ₖ g
  map_one' := one_kronecker_one
  map_mul' a b := by rw [← mul_kronecker_mul, one_mul]
  map_zero' := kronecker_zero _
  map_add' a b := kronecker_add _ _ _

/-- Reindexing equivalence `Fin 2 × (Fin 2 × Fin 2) ≃ Fin 8`. -/
noncomputable def e8 : (Fin 2 × (Fin 2 × Fin 2)) ≃ Fin 8 :=
  Fintype.equivFinOfCardEq (by decide)

/-- **The Dirac ansatz is consistent.**  Williams' eq-(109) construction gives
    a concrete realisation of the Dirac ansatz over the 8×8 complex matrix ring, so
    adopting it introduces no contradiction. Built here as the 4×4 two-qubit
    Pauli realisation `D4`, doubled to 8×8 and reindexed to `Fin 8`.

    (The theorem name `fourth_postulate_consistent` is a historical nickname
    retained as a stable identifier; the object is the 5D Dirac ansatz.) -/
theorem fourth_postulate_consistent :
    Nonempty (DiracData (Matrix (Fin 8) (Fin 8) ℂ)) :=
  ⟨(D4.map doubleHom).map
      (Matrix.reindexAlgEquiv ℂ ℂ e8).toRingEquiv.toRingHom⟩

/-! ## Summary obligation — Theorem 1 (sufficiency), fully assembled

Bundling the pieces: under the Dirac ansatz, the gamma operators
satisfy the 5D Clifford algebra AND the Hamiltonian reproduces the mass-shell
relation — i.e. the 5D Dirac equation (112) is well-posed and consistent with
postulates 1-3. -/
theorem theorem1_sufficiency (D : DiracData R) :
    (∀ i j, D.γ i * D.γ j + D.γ j * D.γ i = 2 * gMetric R i j) :=
  fun i j => gamma_clifford D i j

end DTWilliams.FourthPostulate
