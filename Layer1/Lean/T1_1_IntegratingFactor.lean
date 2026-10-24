/-
  T1.1 — Universal Integrating Factor  λ = φ(q̇)·f(σ)

  Source: P. E. Williams, LA-8370-MS (LASL, 1980), §II, eqs (3)–(5).
  Supplementary Materials: step [L2], Sec. S2 (code index: Sec. S6).
-/

import Mathlib.Tactic
import Mathlib.Analysis.SpecialFunctions.Log.Basic

namespace DTWilliams.IntegratingFactor

/-! ## The composite-system data (eqs 3–4 output)

We work with the LOGARITHMS of the integrating factors, so the multiplicative
product form of eq (5) becomes an additive split (and no positivity/division
bookkeeping is needed):

* `L  q s`     = `ln λ(q̇, σ)`        — system 1, velocity `q`, E-coordinate `s = σ`
* `Lh q sh`    = `ln λ̂(q̇, σ̂)`       — system 2, E-coordinate `sh = σ̂`
* `Lc q s sh`  = `ln λ_c(q̇, σ, σ̂)`   — composite

The composite argument (eqs 3–4) outputs exactly that the two log-ratios
`ln λ − ln λ_c` and `ln λ̂ − ln λ_c` are independent of the velocity `q̇`; we take
that as the two hypotheses. "Independent of `q̇`" is stated as: the value at
velocity `a` equals the value at velocity `b`, for all velocities. -/
structure CompositeIF where
  /-- `ln λ(q̇, σ)`. -/
  L  : ℝ → ℝ → ℝ
  /-- `ln λ̂(q̇, σ̂)`. -/
  Lh : ℝ → ℝ → ℝ
  /-- `ln λ_c(q̇, σ, σ̂)`. -/
  Lc : ℝ → ℝ → ℝ → ℝ
  /-- Ratio `λ/λ_c` is independent of velocity (eq 4 vs eq 3, first ratio):
      `ln λ(a,σ) − ln λ_c(a,σ,σ̂)` does not depend on the velocity `a`. -/
  ratio1_indep : ∀ a b s sh, L a s - Lc a s sh = L b s - Lc b s sh
  /-- Ratio `λ̂/λ_c` is independent of velocity (eq 4 vs eq 3, second ratio):
      `ln λ̂(a,σ̂) − ln λ_c(a,σ,σ̂)` does not depend on the velocity `a`. -/
  ratio2_indep : ∀ a b s sh, Lh a sh - Lc a s sh = Lh b sh - Lc b s sh

variable (C : CompositeIF)

/-- **Velocity-difference of the two systems is system-coordinate-free.**
    Subtracting the two ratio hypotheses cancels the composite factor `λ_c`:
    `ln λ(a,σ) − ln λ̂(a,σ̂)` is independent of the velocity `a`. This is the
    crux that lets the *same* `φ(q̇)` serve both systems. -/
theorem ratioDiff_indep (a b s sh : ℝ) :
    C.L a s - C.Lh a sh = C.L b s - C.Lh b sh := by
  have h1 := C.ratio1_indep a b s sh
  have h2 := C.ratio2_indep a b s sh
  linarith

/-- **T1.1 — Universal product form (eq 5), additive (log) statement.**
    From the composite data there exist a single velocity envelope `Φ` (the same
    for both systems — this is the universality) and system-specific
    coordinate parts `f, fh` with

        ln λ(q̇,σ)  = Φ(q̇) + f(σ),
        ln λ̂(q̇,σ̂) = Φ(q̇) + f̂(σ̂).

    Witnesses: `Φ(q̇) = ln λ(q̇, 0)`, `f(σ) = ln λ(0,σ) − ln λ(0,0)`,
    `f̂(σ̂) = ln λ̂(0,σ̂) − ln λ(0,0)` (base velocity `q̇ = 0` as the representative). -/
theorem universal_product_form :
    ∃ Φ f fh : ℝ → ℝ,
      (∀ q s, C.L q s = Φ q + f s) ∧ (∀ q sh, C.Lh q sh = Φ q + fh sh) := by
  refine ⟨fun q => C.L q 0, fun s => C.L 0 s - C.L 0 0,
          fun sh => C.Lh 0 sh - C.L 0 0, ?_, ?_⟩
  · -- product form for system 1
    intro q s
    show C.L q s = C.L q 0 + (C.L 0 s - C.L 0 0)
    have kA := ratioDiff_indep C q 0 s 0
    have kB := ratioDiff_indep C q 0 0 0
    linarith
  · -- same envelope Φ(q) = L q 0 serves system 2 — the universality
    intro q sh
    show C.Lh q sh = C.L q 0 + (C.Lh 0 sh - C.L 0 0)
    have kB := ratioDiff_indep C q 0 0 sh
    linarith

/-- **T1.1 — Universal integrating factor (eq 5), multiplicative statement.**
    Exponentiating the log form: there is ONE velocity factor `φ` (`= e^Φ`,
    the universal `φ(q̇)`) shared by both systems, with

        λ(q̇,σ)  = φ(q̇)·f(σ),
        λ̂(q̇,σ̂) = φ(q̇)·f̂(σ̂),

    where `λ = e^{ln λ}`. This is Williams' eq (5): `1/φ(q̇)` is the same
    integrating factor for every system; the exact, system-specific differential
    is carried by `f(σ)dσ`. -/
theorem universal_integrating_factor :
    ∃ φ f fh : ℝ → ℝ,
      (∀ q s, Real.exp (C.L q s) = φ q * f s) ∧
      (∀ q sh, Real.exp (C.Lh q sh) = φ q * fh sh) := by
  obtain ⟨Φ, f, fh, h1, h2⟩ := universal_product_form C
  exact ⟨fun q => Real.exp (Φ q), fun s => Real.exp (f s), fun sh => Real.exp (fh sh),
    fun q s => by rw [h1 q s, Real.exp_add], fun q sh => by rw [h2 q sh, Real.exp_add]⟩

/-! ## Remark — the integrating factor and exactness (doc p.18, eq after eq 5)

With `λ = φ(q̇)·f(σ)`, the energy 1-form is `đE = φ(q̇) f(σ) dσ`, so
`đE / φ(q̇) = f(σ) dσ` depends on `σ` ALONE and is therefore exact (it is
`dS`, the mechanical entropy of T1.3). In this model that exactness is the
statement that the σ-part `f` returned above is a function of `σ` only — which
it is by construction (`f : ℝ → ℝ`). The velocity factor `φ` is common to all
systems, so `1/φ(q̇)` is the universal integrating factor. The numerical form
`φ(q̇) = √(1 − q̇²/c²)` and the absolute velocity `q̇₀ = c` (where `φ(q̇₀)=0`)
are T1.2; they are empirical/definitional inputs, not part of this theorem. -/

end DTWilliams.IntegratingFactor
