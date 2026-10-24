/-
  T1.7 / T1.8 — Riemannian Connection and Weyl Geometry

  Source: P. E. Williams, LA-8370-MS (LASL, 1980), §II, eqs (20)–(37).
  Supplementary Materials: step [L6], Sec. S2 (code index: Sec. S6).
-/

import Mathlib.Tactic

open scoped BigOperators

namespace DTWilliams.Connection

/-- First-kind Christoffel symbol from the metric derivatives `dg i j s = ∂_s g_{ij}`
    (eq 29):  `Γ_{sij} = ½(∂_j g_{si} + ∂_i g_{sj} − ∂_s g_{ij})`. -/
noncomputable def Γ₁ {ι : Type*} (dg : ι → ι → ι → ℝ) (s i j : ι) : ℝ :=
  (dg s i j + dg s j i - dg i j s) / 2

/-! ## T1.7 — symmetric, metric-compatible connection -/

/-- **Gauge metric ⇒ `dg` symmetric in (i,j) (eqs 20-22).** For ĝ_{ij}=a_i·a_j
    (a_i = q⁰_{/i}), the product rule gives `dg i j s = da i s · a j + a i · da j s`,
    which is automatically symmetric under `i ↔ j` (commutativity); this is the
    `dg_symm` hypothesis used below. (Exactness `da i s = da s i`, eq 22, is the
    separate input behind the symmetric connection / Riemann structure.) -/
theorem gauge_dg_symm {ι : Type*} (a : ι → ℝ) (da : ι → ι → ℝ) (i j s : ι) :
    da i s * a j + a i * da j s = da j s * a i + a j * da i s := by ring

variable {ι : Type*} (dg : ι → ι → ι → ℝ)

/-- **Symmetric connection (eq 24, first kind).** `Γ_{sij} = Γ_{sji}` — the
    connection is symmetric in its lower index pair. -/
theorem christoffel_first_symm (dg_symm : ∀ i j s, dg i j s = dg j i s) (s i j : ι) :
    Γ₁ dg s i j = Γ₁ dg s j i := by
  unfold Γ₁; rw [dg_symm i j s]; ring

/-- **Metric compatibility (eq 26).** `∂_s g_{ij} = Γ_{jis} + Γ_{ijs}` — recovering
    the metric derivative from the connection. With cyclic permutation (eqs 27-28)
    this inverts to the Christoffel formula `Γ₁` above. -/
theorem metric_compat (dg_symm : ∀ i j s, dg i j s = dg j i s) (i j s : ι) :
    Γ₁ dg j i s + Γ₁ dg i j s = dg i j s := by
  unfold Γ₁; have h := dg_symm j i s; linarith

/-! ### Second-kind connection -/

section SecondKind
variable [Fintype ι]

/-- Second-kind connection `Γ^k_{ij} = Σ_s g^{ks} Γ_{sij}` (raising with `ginv`). -/
noncomputable def Γ₂ (dg : ι → ι → ι → ℝ) (ginv : ι → ι → ℝ) (k i j : ι) : ℝ :=
  ∑ s, ginv k s * Γ₁ dg s i j

/-- **Symmetric connection (eq 24, second kind).** `Γ^k_{ij} = Γ^k_{ji}`. -/
theorem christoffel_second_symm (ginv : ι → ι → ℝ)
    (dg_symm : ∀ i j s, dg i j s = dg j i s) (k i j : ι) :
    Γ₂ dg ginv k i j = Γ₂ dg ginv k j i := by
  unfold Γ₂
  exact Finset.sum_congr rfl
    (fun s _ => by rw [christoffel_first_symm dg dg_symm s i j])

end SecondKind

/-! ## T1.8 — Weyl geometry from conformal rescaling `g → f·g` -/

/-- Metric derivative of the rescaled metric `f·g` (product rule):
    `∂_s(f·g_{ij}) = f·∂_s g_{ij} + g_{ij}·∂_s f`. -/
def dgf {ι : Type*} (g : ι → ι → ℝ) (dg : ι → ι → ι → ℝ) (f : ℝ) (df : ι → ℝ)
    (i j s : ι) : ℝ := f * dg i j s + g i j * df s

/-- **Conformal Christoffel identity (first kind), the core of eq (33).** With the
    gauge vector `φ_k = ½ ∂_k ln f` (i.e. `∂_k f = 2 f φ_k`), the first-kind
    Christoffel of the rescaled metric `f·g` is

        Γ_{sij}[f·g] = f·(Γ_{sij}[g] + g_{si}φ_j + g_{sj}φ_i − g_{ij}φ_s).

    The product-rule cross-terms assemble into exactly the Weyl correction. -/
theorem weyl_first_kind {ι : Type*} (g : ι → ι → ℝ) (dg : ι → ι → ι → ℝ)
    (f : ℝ) (df φ : ι → ℝ) (hφ : ∀ s, df s = 2 * f * φ s) (s i j : ι) :
    Γ₁ (dgf g dg f df) s i j
      = f * (Γ₁ dg s i j + g s i * φ j + g s j * φ i - g i j * φ s) := by
  unfold Γ₁ dgf
  rw [hφ i, hφ j, hφ s]; ring

/-- **Weyl connection (eq 33).** Raising with the inverse of `f·g` (which is
    `g^{ks}/f`), the conformal factor cancels and the second-kind connection of
    `f·g` is the `g`-connection plus the Weyl term:

        Γ^k_{ij}[f·g] = Γ^k_{ij}[g] + Σ_s g^{ks}(g_{si}φ_j + g_{sj}φ_i − g_{ij}φ_s).

    Weyl geometry and the gauge vector `φ` emerge from the rescaling — not assumed. -/
theorem weyl_connection [Fintype ι] (g : ι → ι → ℝ) (dg : ι → ι → ι → ℝ)
    (ginv : ι → ι → ℝ) (f : ℝ) (df φ : ι → ℝ) (hf : f ≠ 0)
    (hφ : ∀ s, df s = 2 * f * φ s) (k i j : ι) :
    (∑ s, (ginv k s / f) * Γ₁ (dgf g dg f df) s i j)
      = Γ₂ dg ginv k i j
        + ∑ s, ginv k s * (g s i * φ j + g s j * φ i - g i j * φ s) := by
  unfold Γ₂
  rw [← Finset.sum_add_distrib]
  apply Finset.sum_congr rfl
  intro s _
  rw [weyl_first_kind g dg f df φ hφ s i j]
  field_simp
  ring

end DTWilliams.Connection
