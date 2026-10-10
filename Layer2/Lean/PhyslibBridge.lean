/-
  Physlib bridge — DT's 5D metric and 5D Dirac ansatz against the community
  library Physlib (leanprover-community/physlib, formerly PhysLean / HepLean).

  Source: P. E. Williams, LA-8370-MS (LASL, 1980), Section III, eqs (60), (108)-(112).
  Full derivation, context and provenance:
  Supplementary_Material.pdf — §S2, Appendix A [C1], [C6].

  Purpose. The DT Lean layer builds its 5D objects by hand. This file proves
  that they coincide with Physlib's independent, standard definitions:

  (1) METRIC. DT's 5D Minkowski form `Q5` (eq 60) and metric `gMetric` (eq 112),
      signature (+,−,−,−,−), ARE Physlib's `minkowskiMatrix` at `d = 4`
      (index type `Fin 1 ⊕ Fin 4`). So the 5D space–time–mass geometry is the
      standard (1+4)-dimensional Minkowski geometry, whose symmetry group is
      Physlib's `LorentzGroup 4`.

  (2) DIRAC. Physlib's Dirac-representation gamma endomorphisms γ⁰,γ¹,γ²,γ³
      together with iγ⁵ realise DT's 5D Dirac ansatz (`DiracData`, eq 108), and
      DT's 5D gamma operators (eq 111) come out as exactly those five matrices.
      This is an independent witness to `fourth_postulate_consistent` (which uses
      a hand-built 8×8 realisation): the textbook 4×4 Dirac matrices already
      satisfy Williams' relations, with the fifth (mass-density) generator iγ⁵.

  Nothing here changes any existing DT theorem; it only adds bridges.
-/

import Physlib.Relativity.MinkowskiMatrix
import Physlib.Relativity.CliffordAlgebra
import Physlib.Relativity.PauliMatrices.Basic
import Layer2.Lean.T2_1_Minkowski5D
import Layer2.Lean.T2_7_DiracFactorisation

namespace DTWilliams.PhyslibBridge

open Matrix DTWilliams.FourthPostulate

/-! ## (1) The 5D metric is Physlib's Minkowski matrix at `d = 4` -/

/-- DT's 5D index `Fin 5` ↔ Physlib's `Fin 1 ⊕ Fin 4` (time ⊕ {space, mass}). -/
abbrev idx5 : Fin 1 ⊕ Fin 4 ≃ Fin 5 := finSumFinEquiv

set_option linter.unusedSimpArgs false in
/-- **DT's metric `g^{jk}` (eq 112) equals Physlib's `minkowskiMatrix` (d = 4).** -/
theorem gMetric_eq_minkowskiMatrix (μ ν : Fin 1 ⊕ Fin 4) :
    gMetric ℝ (idx5 μ) (idx5 ν) = minkowskiMatrix μ ν := by
  rcases μ with μ | μ <;> rcases ν with ν | ν <;> fin_cases μ <;> fin_cases ν <;>
    simp [gMetric, minkowskiMatrix.as_diagonal, diagonal_apply] <;> decide

/-- A 5-vector `(x⁰; x¹,x²,x³,x⁴)` in Physlib's index convention. -/
def vec5 (x0 x1 x2 x3 x4 : ℝ) : Fin 1 ⊕ Fin 4 → ℝ := Sum.elim ![x0] ![x1, x2, x3, x4]

/-- **DT's 5D quadratic form `Q5` (eq 60) is the Minkowski form `vᵀ η v` of
    Physlib's `minkowskiMatrix` at `d = 4`.** -/
theorem Q5_eq_minkowski (x0 x1 x2 x3 x4 : ℝ) :
    Minkowski5.Q5 x0 x1 x2 x3 x4 =
      vec5 x0 x1 x2 x3 x4 ⬝ᵥ (minkowskiMatrix *ᵥ vec5 x0 x1 x2 x3 x4) := by
  simp [Minkowski5.Q5, vec5, minkowskiMatrix.as_diagonal, mulVec_diagonal, dotProduct,
    Fintype.sum_sum_type, Fin.sum_univ_four]
  ring

/-! ## (2) Physlib's Dirac matrices realise DT's 5D Dirac ansatz

A generic construction first: any `b` with `b² = 1` and four `G k` with `G k² = −1`,
anticommuting with `b` and with each other, give a `DiracData` whose eq-(111) gamma
operators are exactly `b, G 0, …, G 3`. -/

section Generic

variable {R : Type*} [Ring R] (b : R) (G : Fin 4 → R)
  (hb : b * b = 1) (hG : ∀ k, G k * G k = -1)
  (hbG : ∀ k, G k * b = -(b * G k)) (hGG : ∀ i j, i ≠ j → G i * G j = -(G j * G i))

include hb hbG in
lemma bG_mul_bG (i j : Fin 4) : b * G i * (b * G j) = -(G i * G j) := by
  calc b * G i * (b * G j) = b * (G i * b) * G j := by simp only [mul_assoc]
    _ = b * (-(b * G i)) * G j := by rw [hbG]
    _ = -(b * b * G i * G j) := by noncomm_ring
    _ = -(G i * G j) := by rw [hb, one_mul]

/-- The `DiracData` built from `(b, G)`, with `α̃_k = −b G_k` and `β̃ = b`. -/
def ofClifford : DiracData R where
  α k := -(b * G k)
  β := b
  α_sq k := by rw [neg_mul_neg, bG_mul_bG b G hb hbG, hG, neg_neg]
  α_anticomm i j h := by
    rw [neg_mul_neg, neg_mul_neg, bG_mul_bG b G hb hbG, bG_mul_bG b G hb hbG, hGG i j h,
      neg_neg, add_neg_cancel]
  β_sq := hb
  αβ_anticomm k := by
    rw [neg_mul, mul_assoc, hbG]
    simp only [mul_neg, neg_neg, add_neg_cancel]

/-- The eq-(111) gamma operators of `ofClifford` are exactly `b, G 0, …, G 3`. -/
theorem ofClifford_γ (k : Fin 5) :
    (ofClifford b G hb hG hbG hGG).γ k = (Fin.cons b G : Fin 5 → R) k := by
  refine Fin.cases ?_ (fun j => ?_) k
  · simp [ofClifford]
  · simp only [DiracData.γ, Fin.succ_ne_zero, ↓reduceDIte, Fin.pred_succ, Fin.cons_succ,
      ofClifford]
    rw [mul_neg, neg_neg, ← mul_assoc, hb, one_mul]

end Generic

/-! ### Instantiation with Physlib's gamma endomorphisms -/

open Fermion.Dirac Complex

/-- The four "spacelike" 5D generators: γ¹, γ², γ³ and the mass-density generator iγ⁵. -/
noncomputable def Gphys : Fin 4 → Module.End ℂ Fermion.Dirac :=
  ![gamma (Sum.inr 0), gamma (Sum.inr 1), gamma (Sum.inr 2), I • gamma5]

lemma iγ5_mul_self : (I • gamma5) * (I • gamma5) = -1 := by
  rw [smul_mul_smul_comm, I_mul_I, gamma5_mul_self, neg_one_smul]

lemma iγ5_mul_gamma (μ : Fin 1 ⊕ Fin 3) :
    (I • gamma5) * gamma μ = -(gamma μ * (I • gamma5)) := by
  rw [smul_mul_assoc, gamma5_mul_gamma, mul_smul_comm, smul_neg]

lemma gamma_mul_iγ5 (μ : Fin 1 ⊕ Fin 3) :
    gamma μ * (I • gamma5) = -((I • gamma5) * gamma μ) := by
  rw [iγ5_mul_gamma, neg_neg]

lemma gamma_inr_anti (a b : Fin 3) (h : a ≠ b) :
    gamma (Sum.inr a) * gamma (Sum.inr b) = -(gamma (Sum.inr b) * gamma (Sum.inr a)) :=
  gamma_mul_gamma_of_ne (by simpa using h)

lemma gamma_inr_anti_γ0 (a : Fin 3) :
    gamma (Sum.inr a) * gamma (Sum.inl 0) = -(gamma (Sum.inl 0) * gamma (Sum.inr a)) :=
  gamma_mul_gamma_of_ne (by simp)

-- Case splits below are written out one by one: letting unification search across
-- mismatched gamma indices unfolds Physlib's `gamma` and exhausts the heartbeat limit.
lemma Gphys_sq (k : Fin 4) : Gphys k * Gphys k = -1 := by
  fin_cases k
  · exact gamma_inr_mul_self 0
  · exact gamma_inr_mul_self 1
  · exact gamma_inr_mul_self 2
  · exact iγ5_mul_self

lemma Gphys_anti_γ0 (k : Fin 4) :
    Gphys k * gamma (Sum.inl 0) = -(gamma (Sum.inl 0) * Gphys k) := by
  fin_cases k
  · exact gamma_inr_anti_γ0 0
  · exact gamma_inr_anti_γ0 1
  · exact gamma_inr_anti_γ0 2
  · exact iγ5_mul_gamma _

lemma Gphys_anti (i j : Fin 4) (h : i ≠ j) : Gphys i * Gphys j = -(Gphys j * Gphys i) := by
  fin_cases i <;> fin_cases j
  · exact absurd rfl h
  · exact gamma_inr_anti 0 1 (by decide)
  · exact gamma_inr_anti 0 2 (by decide)
  · exact gamma_mul_iγ5 _
  · exact gamma_inr_anti 1 0 (by decide)
  · exact absurd rfl h
  · exact gamma_inr_anti 1 2 (by decide)
  · exact gamma_mul_iγ5 _
  · exact gamma_inr_anti 2 0 (by decide)
  · exact gamma_inr_anti 2 1 (by decide)
  · exact absurd rfl h
  · exact gamma_mul_iγ5 _
  · exact iγ5_mul_gamma _
  · exact iγ5_mul_gamma _
  · exact iγ5_mul_gamma _
  · exact absurd rfl h

/-- **Physlib's Dirac gammas realise DT's 5D Dirac ansatz** (operator level):
    `β̃ = γ⁰`, `α̃_k = −γ⁰Γ_k` with `Γ = (γ¹, γ², γ³, iγ⁵)`. -/
noncomputable def physlibDirac : DiracData (Module.End ℂ Fermion.Dirac) :=
  ofClifford (gamma (Sum.inl 0)) Gphys gamma_inl_zero_mul_self Gphys_sq Gphys_anti_γ0 Gphys_anti

/-- DT's eq-(111) 5D gamma operators for `physlibDirac` are `γ⁰, γ¹, γ², γ³, iγ⁵`. -/
theorem physlibDirac_γ (k : Fin 5) :
    physlibDirac.γ k = (Fin.cons (gamma (Sum.inl 0)) Gphys : Fin 5 → Module.End ℂ Fermion.Dirac) k :=
  ofClifford_γ _ _ _ _ _ _ k

/-- The same realisation as 4×4 complex matrices (Physlib's Dirac representation,
    via `endEquivMatrix`). Its `β̃` is literally Physlib's matrix `γ0`. -/
noncomputable def physlibDiracMatrix : DiracData (Matrix (Fin 4) (Fin 4) ℂ) :=
  physlibDirac.map endEquivMatrix.toRingEquiv.toRingHom

theorem physlibDiracMatrix_β : physlibDiracMatrix.β = spaceTime.γ0 := by
  simp [physlibDiracMatrix, DiracData.map, physlibDirac, ofClifford, spaceTime.γ0]

/-- **Consequence.** Physlib's matrices satisfy DT's 5D Clifford relation with DT's
    metric — `{γ^i, γ^j} = 2 g^{ij}` on `Fin 5` (Theorem 1, eq 112) — i.e. the 5D
    Dirac equation is realised by the standard 4×4 Dirac matrices plus iγ⁵. -/
theorem physlib_clifford5 (i j : Fin 5) :
    physlibDiracMatrix.γ i * physlibDiracMatrix.γ j +
      physlibDiracMatrix.γ j * physlibDiracMatrix.γ i =
      2 * gMetric (Matrix (Fin 4) (Fin 4) ℂ) i j :=
  theorem1_sufficiency physlibDiracMatrix i j

/-! ## (3) Pauli matrices: DT's hand-built `sx, sy, sz` are Physlib's σ¹, σ², σ³ -/

theorem sx_eq : Pauli.sx = PauliMatrix.pauliMatrix (Sum.inr 0) := by simp [Pauli.sx, PauliMatrix.pauliMatrix]
theorem sy_eq : Pauli.sy = PauliMatrix.pauliMatrix (Sum.inr 1) := by simp [Pauli.sy, PauliMatrix.pauliMatrix]
theorem sz_eq : Pauli.sz = PauliMatrix.pauliMatrix (Sum.inr 2) := by simp [Pauli.sz, PauliMatrix.pauliMatrix]

end DTWilliams.PhyslibBridge
