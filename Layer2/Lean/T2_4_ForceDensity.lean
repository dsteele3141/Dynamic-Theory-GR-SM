/-
  T2.4 — Force Density Five-Vector  K = Div₅{T}

  Source: P. E. Williams, LA-8370-MS (LASL, 1980), Section III.
  Full derivation, context and provenance:
  Supplementary_Material.pdf — §S2, Appendix A [C5].
-/

import Mathlib.Tactic
import Mathlib.Basic.Complex.Basic
import Mathlib.Data.Matrix.Basic
import Mathlib.Algebra.BigOperators.Fin

namespace DTWilliams.ForceDensity

open Matrix
open Complex (I)

/-- 5×5 antisymmetric Dynamic-Theory field tensor (T2.2/T2.3 convention). Indices
    0..4 = (entropy/time, x, y, z, mass-density); iE,B are the EM field, V,iV₀ the
    new 5D gauge field. -/
noncomputable def Ften (E1 E2 E3 B1 B2 B3 V0 V1 V2 V3 : ℂ) : Matrix (Fin 5) (Fin 5) ℂ :=
  !![  0,      I*E1,  I*E2,  I*E3,  I*V0;
     -I*E1,    0,     B3,   -B2,    V1;
     -I*E2,   -B3,    0,     B1,    V2;
     -I*E3,    B2,   -B1,    0,     V3;
     -I*V0,   -V1,   -V2,   -V3,    0   ]

/-- **F is antisymmetric:** `F_{ij} = −F_{ji}`. -/
theorem Ften_antisymm (E1 E2 E3 B1 B2 B3 V0 V1 V2 V3 : ℂ) (i j : Fin 5) :
    Ften E1 E2 E3 B1 B2 B3 V0 V1 V2 V3 i j
      = - Ften E1 E2 E3 B1 B2 B3 V0 V1 V2 V3 j i := by
  fin_cases i <;> fin_cases j <;> simp [Ften]

/-! ## Force-density components `K_j = (1/c) F_{jk} J^k`

Current five-vector `J^k = ![icρ, J₁, J₂, J₃, J₄]`. Each `K_j = (1/c)·(F·J)_j`. -/

/-- **K₀ (energy-flow component).** `K₀ = (i/c)(J̄·Ē + J₄V₀)`. -/
theorem Kforce_0 (E1 E2 E3 B1 B2 B3 V0 V1 V2 V3 ρ J1 J2 J3 J4 c : ℂ) (hc : c ≠ 0) :
    (1 / c) * (Ften E1 E2 E3 B1 B2 B3 V0 V1 V2 V3).mulVec ![I*c*ρ, J1, J2, J3, J4] 0
      = (I / c) * (J1*E1 + J2*E2 + J3*E3 + J4*V0) := by
  simp [Ften, Matrix.mulVec]
  field_simp
  ring

/-- **K₁ (Lorentz + new V-field force, x-component).**
    `K₁ = ρE₁ + (J₂B₃ − J₃B₂)/c + J₄V₁/c`. -/
theorem Kforce_1 (E1 E2 E3 B1 B2 B3 V0 V1 V2 V3 ρ J1 J2 J3 J4 c : ℂ) (hc : c ≠ 0) :
    (1 / c) * (Ften E1 E2 E3 B1 B2 B3 V0 V1 V2 V3).mulVec ![I*c*ρ, J1, J2, J3, J4] 1
      = ρ*E1 + (J2*B3 - J3*B2)/c + J4*V1/c := by
  simp [Ften, Matrix.mulVec]
  field_simp
  ring_nf
  simp only [Complex.I_sq]
  ring

/-- **K₂ (y-component).** `K₂ = ρE₂ + (J₃B₁ − J₁B₃)/c + J₄V₂/c`. -/
theorem Kforce_2 (E1 E2 E3 B1 B2 B3 V0 V1 V2 V3 ρ J1 J2 J3 J4 c : ℂ) (hc : c ≠ 0) :
    (1 / c) * (Ften E1 E2 E3 B1 B2 B3 V0 V1 V2 V3).mulVec ![I*c*ρ, J1, J2, J3, J4] 2
      = ρ*E2 + (J3*B1 - J1*B3)/c + J4*V2/c := by
  simp [Ften, Matrix.mulVec]
  field_simp
  ring_nf
  simp only [Complex.I_sq]
  ring

/-- **K₃ (z-component).** `K₃ = ρE₃ + (J₁B₂ − J₂B₁)/c + J₄V₃/c`. -/
theorem Kforce_3 (E1 E2 E3 B1 B2 B3 V0 V1 V2 V3 ρ J1 J2 J3 J4 c : ℂ) (hc : c ≠ 0) :
    (1 / c) * (Ften E1 E2 E3 B1 B2 B3 V0 V1 V2 V3).mulVec ![I*c*ρ, J1, J2, J3, J4] 3
      = ρ*E3 + (J1*B2 - J2*B1)/c + J4*V3/c := by
  simp [Ften, Matrix.mulVec]
  field_simp
  ring_nf
  simp only [Complex.I_sq]
  ring

/-- **K₄ (new 5D mass-density force).** `K₄ = ρV₀ − J̄·V̄/c` — depends only on the
    V-field and current (no E, B): a force absent from standard 4D EM. -/
theorem Kforce_4 (E1 E2 E3 B1 B2 B3 V0 V1 V2 V3 ρ J1 J2 J3 J4 c : ℂ) (hc : c ≠ 0) :
    (1 / c) * (Ften E1 E2 E3 B1 B2 B3 V0 V1 V2 V3).mulVec ![I*c*ρ, J1, J2, J3, J4] 4
      = ρ*V0 - (J1*V1 + J2*V2 + J3*V3)/c := by
  simp [Ften, Matrix.mulVec]
  field_simp
  ring_nf
  simp only [Complex.I_sq]
  ring

/-! ## Classical limit `V₀ = V̄ = 0, J₄ = 0`: 4D Lorentz force recovery -/

/-- **Classical limit (spatial): standard 4D Lorentz force density.** With the
    V-field off and no mass current (`J₄=0`), `K₁ = ρE₁ + (J̄×B̄)₁/c`. -/
theorem Kforce_1_classical (E1 E2 E3 B1 B2 B3 ρ J1 J2 J3 c : ℂ) (hc : c ≠ 0) :
    (1 / c) * (Ften E1 E2 E3 B1 B2 B3 0 0 0 0).mulVec ![I*c*ρ, J1, J2, J3, 0] 1
      = ρ*E1 + (J2*B3 - J3*B2)/c := by
  simp [Ften, Matrix.mulVec]
  field_simp
  ring_nf
  simp only [Complex.I_sq]
  ring

/-- **Classical limit: the new mass-density force vanishes.** `K₄ = 0` when the
    V-field is off — the 5th force is purely a Dynamic-Theory effect. -/
theorem Kforce_4_classical (E1 E2 E3 B1 B2 B3 ρ J1 J2 J3 c : ℂ) :
    (1 / c) * (Ften E1 E2 E3 B1 B2 B3 0 0 0 0).mulVec ![I*c*ρ, J1, J2, J3, 0] 4 = 0 := by
  simp [Ften, Matrix.mulVec]

end DTWilliams.ForceDensity
