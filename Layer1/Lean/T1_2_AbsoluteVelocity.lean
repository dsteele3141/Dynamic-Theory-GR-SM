/-
  T1.2 — The Absolute Velocity

  Source: P. E. Williams, LA-8370-MS (LASL, 1980), §II (energy ratio → φ(q̇₀)=0).
  Supplementary Materials: step [L3], Sec. S2 (code index: Sec. S6).
-/

import Mathlib.Tactic
import Mathlib.Analysis.SpecialFunctions.Sqrt

namespace DTWilliams.AbsoluteVelocity

variable (φ : ℝ → ℝ) (I : ℝ)

/-- Energy transferred at constant velocity `q` between two fixed E-conservative
    curves: `ΔE(q) = φ(q)·I`, with `I = ∫_{σ₁}^{σ₂} f dσ` common to all velocities. -/
def energy (q : ℝ) : ℝ := φ q * I

/-- **Energy ratio law (doc p.19).** `ΔE(q)/ΔE(q₃) = φ(q)/φ(q₃)` — the common
    σ-integral `I` cancels, leaving a pure function of the two velocities. -/
theorem energy_ratio (q q₃ : ℝ) (hI : I ≠ 0) :
    energy φ I q / energy φ I q₃ = φ q / φ q₃ := by
  unfold energy
  rw [mul_div_mul_right _ _ hI]

/-- **Scaling law (doc p.20).** `ΔE(q) = [ΔE₃/φ(q₃)]·φ(q)`. -/
theorem energy_scaled (q q₃ : ℝ) (h₃ : φ q₃ ≠ 0) :
    energy φ I q = (energy φ I q₃ / φ q₃) * φ q := by
  unfold energy
  rw [mul_div_cancel_left₀ I h₃]
  ring

/-- **Energy-free transfer ⟺ φ vanishes (doc p.20).** Between two DISTINCT
    E-conservative curves (`I ≠ 0`), the constant-velocity transfer is zero iff
    the universal factor vanishes at that velocity. -/
theorem no_exchange_iff (q : ℝ) (hI : I ≠ 0) :
    energy φ I q = 0 ↔ φ q = 0 := by
  unfold energy; rw [mul_eq_zero]; tauto

/-- `q₀` is an **absolute velocity** iff the universal factor vanishes there. -/
def IsAbsoluteVelocity (q₀ : ℝ) : Prop := φ q₀ = 0

/-- At the absolute velocity, a constant-velocity process transfers no energy. -/
theorem no_exchange_at_absolute (q₀ : ℝ) (h : IsAbsoluteVelocity φ q₀) :
    energy φ I q₀ = 0 := by
  unfold energy; rw [show φ q₀ = 0 from h, zero_mul]

/-- **Uniqueness.** A strictly monotone universal factor has at most one absolute
    velocity (the empirical φ below is strictly decreasing on the physical range). -/
theorem absolute_velocity_unique (hmono : StrictMono φ ∨ StrictAnti φ)
    {a b : ℝ} (ha : IsAbsoluteVelocity φ a) (hb : IsAbsoluteVelocity φ b) : a = b := by
  have hφ : φ a = φ b := (ha : φ a = 0).trans (hb : φ b = 0).symm
  rcases hmono with h | h
  · exact h.injective hφ
  · exact h.injective hφ

/-! ## The empirical velocity factor (doc p.64) — NOT derived from the axioms -/

/-- Empirical form (mass-spectrometer data, doc p.64): `φ(q̇) = √(1 − q̇²/c²)`. -/
noncomputable def φEmp (c q : ℝ) : ℝ := Real.sqrt (1 - q ^ 2 / c ^ 2)

/-- **The absolute velocity of the empirical factor is q̇₀ = c.** -/
theorem φEmp_absolute (c : ℝ) (hc : c ≠ 0) : IsAbsoluteVelocity (φEmp c) c := by
  show Real.sqrt (1 - c ^ 2 / c ^ 2) = 0
  rw [div_self (pow_ne_zero 2 hc), sub_self, Real.sqrt_zero]

/-- **Below the absolute velocity the factor is strictly positive:** for
    `0 ≤ q̇ < c` (c>0), `φEmp c q̇ > 0`. So no slower constant-velocity process is
    energy-free, pinning q̇₀ = c uniquely on the physical range `[0, c]`. -/
theorem φEmp_pos_below (c q : ℝ) (hc : 0 < c) (hq : 0 ≤ q) (hlt : q < c) :
    0 < φEmp c q := by
  have hq2 : q ^ 2 < c ^ 2 := by
    nlinarith [mul_pos (show (0:ℝ) < c - q by linarith) (show (0:ℝ) < c + q by linarith)]
  have hpos : 0 < 1 - q ^ 2 / c ^ 2 := by
    rw [sub_pos, div_lt_one (by positivity)]; exact hq2
  unfold φEmp
  exact Real.sqrt_pos.mpr hpos

end DTWilliams.AbsoluteVelocity
