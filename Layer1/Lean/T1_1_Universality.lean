/-
  T1.1 (necessity) — the universal φ(q̇) is FORCED, and it is GLOBAL

  Source: P. E. Williams, LA-8370-MS (LASL, 1980), §II, eqs (3)–(5).
  Supplementary Materials: step [L2], Sec. S2 (code index: Sec. S6).
-/

import Mathlib.Tactic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Layer1.Lean.T1_1_IntegratingFactor

namespace DTWilliams.Universality

variable {ι : Type*}

/-- **Pairwise composite-consistency** (the family generalisation of
    `IntegratingFactor.ratioDiff_indep`). Forming the composite of any two systems
    `i, j` outputs, via the eqs-3/4 comparison, that the velocity-difference of
    their log integrating factors does not depend on the velocity:
    `ln λ_i(a,σ) − ln λ_j(a,σ') = ln λ_i(b,σ) − ln λ_j(b,σ')` for all velocities
    `a,b` and coordinates `σ,σ'`. (The reflexive case `i = j` is the same output
    of composing a system with a reference copy of itself.) -/
def PairwiseConsistent (L : ι → ℝ → ℝ → ℝ) : Prop :=
  ∀ (i j : ι) (a b s t : ℝ), L i a s - L j a t = L i b s - L j b t

variable {L : ι → ℝ → ℝ → ℝ}

/-- **Universal velocity envelope (zeroth-law transitivity).** The velocity
    increment of `ln λ` is one and the SAME function of q̇ for every system:
    `ln λ_i(a,σ) − ln λ_i(b,σ)` does not depend on the system `i` or its
    coordinate. This is the precise sense in which the velocity dependence φ(q̇)
    is universal — shared, not merely pairwise-matched. -/
theorem velocity_envelope_universal (H : PairwiseConsistent L)
    (i j : ι) (a b s t : ℝ) :
    L i a s - L i b s = L j a t - L j b t := by
  have h := H i j a b s t; linarith

/-- **Global universal product form (eq 5, additive/log statement).** From
    pairwise consistency, fixing ANY reference system `i₀` and coordinate `s₀`,
    there is a SINGLE velocity envelope `Φ(q̇) := ln λ_{i₀}(q̇, s₀)` — the same for
    every system — and system-specific coordinate parts `f i`, with
        ln λ_i(q̇,σ) = Φ(q̇) + f_i(σ)   for ALL systems i.
    The universality of φ is therefore forced, not assumed system-by-system. -/
theorem global_universal_form (H : PairwiseConsistent L) (i₀ : ι) (s₀ : ℝ) :
    ∃ (Φ : ℝ → ℝ) (f : ι → ℝ → ℝ),
      (∀ i q s, L i q s = Φ q + f i s) ∧ (∀ q, Φ q = L i₀ q s₀) := by
  refine ⟨fun q => L i₀ q s₀, fun i s => L i 0 s - L i₀ 0 s₀, ?_, fun _ => rfl⟩
  intro i q s
  show L i q s = L i₀ q s₀ + (L i 0 s - L i₀ 0 s₀)
  have h := H i i₀ q 0 s s₀
  linarith

/-- **Global universal integrating factor (eq 5, multiplicative statement).**
    Exponentiating: there is ONE velocity factor `φ` (`= e^Φ`) shared by the whole
    family, with `λ_i(q̇,σ) = φ(q̇)·f_i(σ)` for every system `i`. So `1/φ(q̇)` is the
    same integrating factor for all systems — Williams' universality, eq (5). -/
theorem global_universal_factor (H : PairwiseConsistent L) (i₀ : ι) (s₀ : ℝ) :
    ∃ (φ : ℝ → ℝ) (f : ι → ℝ → ℝ),
      ∀ i q s, Real.exp (L i q s) = φ q * f i s := by
  obtain ⟨Φ, f, hL, _⟩ := global_universal_form H i₀ s₀
  refine ⟨fun q => Real.exp (Φ q), fun i s => Real.exp (f i s), fun i q s => ?_⟩
  rw [hL i q s, Real.exp_add]

/-- **Velocity-independence IS the state-function property.** The hypothesis
    `ratio1_indep` (the ratio `λ/λ_c` is velocity-independent) holds IFF the
    composite log-factor differs from the subsystem's by a velocity-INDEPENDENT
    function `R(σ,σ̂)`:  `ln λ_c(q̇,σ,σ̂) = ln λ(q̇,σ) − R(σ,σ̂)`.
    The velocity-freeness of `R = ln λ − ln λ_c` is exactly the statement that the
    composite's E-conservative coordinate σ_c is a genuine state function — what
    Carathéodory gives for the composite. Hence the velocity-independence taken as
    a hypothesis in T1.1 is not an extra assumption but Carathéodory-for-composite. -/
theorem ratioIndep_iff_stateFunction (Lsub : ℝ → ℝ → ℝ) (Lc : ℝ → ℝ → ℝ → ℝ) :
    (∀ a b s sh, Lsub a s - Lc a s sh = Lsub b s - Lc b s sh) ↔
      (∃ R : ℝ → ℝ → ℝ, ∀ q s sh, Lc q s sh = Lsub q s - R s sh) := by
  constructor
  · intro h
    refine ⟨fun s sh => Lsub 0 s - Lc 0 s sh, fun q s sh => ?_⟩
    have := h q 0 s sh; linarith
  · rintro ⟨R, hR⟩ a b s sh
    rw [hR a s sh, hR b s sh]; ring

/-- **Bridge to the two-system result.** The original `CompositeIF` data (T1.1)
    already exhibits the universal-velocity-envelope conclusion for its pair: the
    velocity-increment of system 1 equals that of system 2 — they carry the SAME
    φ(q̇). The family notion above generalises this from a pair to any class. -/
theorem compositeIF_velocity_universal
    (C : DTWilliams.IntegratingFactor.CompositeIF) (a b s sh : ℝ) :
    C.L a s - C.L b s = C.Lh a sh - C.Lh b sh := by
  have h := DTWilliams.IntegratingFactor.ratioDiff_indep C a b s sh
  linarith

end DTWilliams.Universality
