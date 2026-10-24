/-
  T1.4 — Mechanical Maxwell Relations and Energy Capacities

  Source: P. E. Williams, LA-8370-MS (LASL, 1980), §II, eqs (7)–(11).
  Supplementary Materials: step [L3], Sec. S2 (code index: Sec. S6).
-/

import Mathlib.Tactic
import Mathlib.Analysis.Calculus.Deriv.Mul
import Mathlib.Analysis.Calculus.Deriv.Add

namespace DTWilliams.MaxwellRelations

/-! ## (1) Exact differentials of the mechanical potentials (eqs 7-10, corrected)

The state-function values `φ, φ', S, q, F` and the formal differentials
`dS, dq, dF, dq̇` are treated as independent reals (the 1-form algebra). Each
identity records that a Legendre term cancels, leaving the corrected
coefficients (see the companion Sage script, SI Appendix A). -/

/-- **eq (8):** `H = U − F·q ⟹ dH = dU − F·dq − q·dF = φ·dS − q·dF`.
    The `F·dq` term cancels. -/
theorem dH_eq (φ q F dS dq dF : ℝ) :
    (φ * dS + F * dq) - F * dq - q * dF = φ * dS - q * dF := by ring

/-- **eq (9), CORRECTED:** `K = U − φ·S ⟹ dK = dU − φ'·S·dq̇ − φ·dS
    = F·dq − S·φ'·dq̇`. The `φ·dS` term cancels; the coefficient of `dq` is `F`
    (NOT `φ·F` as printed on p.31). -/
theorem dK_eq (φ φ' S F dS dq dqd : ℝ) :
    (φ * dS + F * dq) - φ' * S * dqd - φ * dS = F * dq - S * φ' * dqd := by ring

/-- **eq (10), CORRECTED:** `G = H − φ·S ⟹ dG = dH − φ'·S·dq̇ − φ·dS
    = −q·dF − S·φ'·dq̇`. The coefficient of `dF` is `−q` (the printed `+q·dF`
    has the wrong sign). -/
theorem dG_eq (φ φ' S q dS dF dqd : ℝ) :
    ((φ * dS - q * dF)) - φ' * S * dqd - φ * dS = -(q * dF) - S * φ' * dqd := by ring

/-! ## (2) The Maxwell relations (eq 11) — Clairaut applied to the potentials

For `P : ℝ → ℝ → ℝ`, `px P` / `py P` are the first partials in the two slots. -/

/-- First partial `∂P/∂x` (slot 1). -/
noncomputable def px (P : ℝ → ℝ → ℝ) (x y : ℝ) : ℝ := deriv (fun t => P t y) x
/-- Second partial `∂P/∂y` (slot 2). -/
noncomputable def py (P : ℝ → ℝ → ℝ) (x y : ℝ) : ℝ := deriv (fun t => P x t) y

/-- **Maxwell relations (eq 11) = Clairaut symmetry of the potential.** For a C²
    potential `P` the mixed partials commute (Clairaut/Schwarz — the cited
    standard theorem, here the hypothesis `clairaut`): `∂/∂y(∂P/∂x) = ∂/∂x(∂P/∂y)`.

    The four Maxwell relations are this applied to U, H, K, G with the partials
    identified from eqs (7)-(10):
      • U(S,q):  ∂U/∂S=φ, ∂U/∂q=F   ⟹  ∂φ/∂q|_S = ∂F/∂S|_q              (Maxwell 1)
      • H(S,F):  ∂H/∂S=φ, ∂H/∂F=−q  ⟹  ∂φ/∂F|_S = −∂q/∂S|_F             (Maxwell 2)
      • K(q̇,q): ∂K/∂q̇=−Sφ', ∂K/∂q=F ⟹ φ'·∂S/∂q = −∂F/∂q̇               (Maxwell 3)
      • G(q̇,F): ∂G/∂q̇=−Sφ', ∂G/∂F=−q ⟹ φ'·∂S/∂F = ∂q/∂q̇               (Maxwell 4) -/
theorem maxwell_relation (P : ℝ → ℝ → ℝ)
    (clairaut : ∀ x y, deriv (fun y => px P x y) y = deriv (fun x => py P x y) x)
    (x y : ℝ) :
    deriv (fun y => px P x y) y = deriv (fun x => py P x y) x :=
  clairaut x y

/-- **Product-rule step behind the DT-form of Maxwell 3 and 4.** Since `φ'(q̇)`
    does not depend on `q`, `∂(−S·φ')/∂q = −(∂S/∂q)·φ'`. So Clairaut on `K`,
    `∂(−Sφ')/∂q = ∂F/∂q̇`, becomes the stated DT-form `φ'·(∂S/∂q) = −(∂F/∂q̇)`. -/
theorem coeff_partial_q (S : ℝ → ℝ) (φ' q : ℝ) (hS : DifferentiableAt ℝ S q) :
    deriv (fun q => -(S q * φ')) q = -(deriv S q * φ') := by
  have h : (fun q => -(S q * φ')) = (fun q => S q * (-φ')) := by funext q; ring
  rw [h, deriv_mul_const hS]; ring

/-! ## (3) Energy-capacity relations (page 27)

Scalar partials at a point: `Sqd = (∂S/∂q̇)_q`, `Sq = (∂S/∂q)_q̇`,
`Fqd = (∂F/∂q̇)_q`, `Fq = (∂F/∂q)_q̇`. At constant `F` the chain rule gives
`(∂S/∂q̇)_F = Sqd + Sq·(∂q/∂q̇)_F` with `(∂q/∂q̇)_F = −Fqd/Fq`. The energy
capacities are `Cq = φ·(∂S/∂q̇)_q` and `CF = φ·(∂S/∂q̇)_F`. -/

/-- **Cq − CF (page 27), using Maxwell 3.** With Maxwell 3 in the form
    `S_q = −F_q̇/φ'`, the capacity difference is
    `Cq − CF = (φ/φ')·(∂q/∂q̇)_F·F_q̇`,  `(∂q/∂q̇)_F = −F_q̇/F_q`. -/
theorem Cq_minus_CF (φ φ' Sqd Sq Fqd Fq : ℝ) (hFq : Fq ≠ 0) (hφ' : φ' ≠ 0)
    (maxwell3 : Sq = -Fqd / φ') :
    φ * Sqd - φ * (Sqd + Sq * (-Fqd / Fq)) = (φ / φ') * (-Fqd / Fq) * Fqd := by
  subst maxwell3; field_simp; ring

/-- **CF / Cq (page 27)** — pure chain rule, no Maxwell relation needed:
    `CF/Cq = (∂F/∂q)_S / (∂F/∂q)_q̇`,  `(∂F/∂q)_S = F_q − F_q̇·S_q/S_q̇`. -/
theorem CF_over_Cq (φ Sqd Sq Fqd Fq : ℝ) (hSqd : Sqd ≠ 0) (hFq : Fq ≠ 0) (hφ : φ ≠ 0) :
    (φ * (Sqd + Sq * (-Fqd / Fq))) / (φ * Sqd) = (Fq - Fqd * Sq / Sqd) / Fq := by
  field_simp; ring

end DTWilliams.MaxwellRelations
