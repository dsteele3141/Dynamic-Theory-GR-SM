/-
  T2.8 — 5D Dirac Equation with Fields (magnetic-moment coupling)

  Source: P. E. Williams, LA-8370-MS (LASL, 1980), Section III, eqs (116)-(121).
  Full derivation, context and provenance:
  Supplementary_Material.pdf — §S2, Appendix A [C8].
-/

import Mathlib.Tactic
import Mathlib.Basic.Complex.Basic
import Mathlib.Data.Matrix.Basic
import Mathlib.Algebra.BigOperators.Fin

namespace DTWilliams.DiracWithFields

open Matrix
open Complex (I)

/-- The 5D spin tensor `σ^{jk} = ½[γ^j,γ^k]` (Williams, doc p.118): entries are
    the velocity `ẋ^a`, the 3-spin `s_a` (with the 2i factors), and `n^a ≡ σ^{a4}`. -/
noncomputable def sigma (xd1 xd2 xd3 xd4 s1 s2 s3 n1 n2 n3 : ℂ) :
    Matrix (Fin 5) (Fin 5) ℂ :=
  !![   0,        xd1,      xd2,      xd3,     xd4;
      -xd1,        0,   -2*I*s3,   2*I*s2,     n1;
      -xd2,    2*I*s3,       0,  -2*I*s1,      n2;
      -xd3,   -2*I*s2,   2*I*s1,       0,      n3;
      -xd4,      -n1,      -n2,      -n3,       0 ]

/-- The field tensor `F_{jk}` (doc p.118), made antisymmetric (the printed
    `F_{30}=+E₃` is a source typo; antisymmetry requires `−E₃`). -/
def Ffield (E1 E2 E3 B1 B2 B3 V0 V1 V2 V3 : ℂ) : Matrix (Fin 5) (Fin 5) ℂ :=
  !![   0,    E1,    E2,    E3,    V0;
      -E1,     0,   -B3,    B2,    V1;
      -E2,    B3,     0,   -B1,    V2;
      -E3,   -B2,    B1,     0,    V3;
      -V0,   -V1,   -V2,   -V3,     0 ]

/-- **Magnetic-moment coupling (eq 119 → eq 121).** The spin–field term contracts
    to the Pauli coupling `2B·s` plus the new 5D V-field couplings:

        −(i/2) Σ_{j,k} F_{jk} σ^{jk}
            = 2(B₁s₁+B₂s₂+B₃s₃) − i(E₁ẋ¹+E₂ẋ²+E₃ẋ³) − i V₀ẋ⁴ − i(V₁n¹+V₂n²+V₃n³).

    The `2B·s` term is the ordinary magnetic moment; the `V₀ẋ⁴` and `n·V` terms
    are the Dynamic-Theory anomalous-moment couplings (present even for a neutral
    particle — prediction P4). -/
theorem magnetic_coupling
    (xd1 xd2 xd3 xd4 s1 s2 s3 n1 n2 n3
     E1 E2 E3 B1 B2 B3 V0 V1 V2 V3 : ℂ) :
    -(I/2) * (∑ j, ∑ k,
        Ffield E1 E2 E3 B1 B2 B3 V0 V1 V2 V3 j k
          * sigma xd1 xd2 xd3 xd4 s1 s2 s3 n1 n2 n3 j k)
      = 2*(B1*s1 + B2*s2 + B3*s3)
        - I*(E1*xd1 + E2*xd2 + E3*xd3)
        - I*V0*xd4
        - I*(V1*n1 + V2*n2 + V3*n3) := by
  simp [Ffield, sigma, Fin.sum_univ_five]
  ring_nf
  simp only [Complex.I_sq]
  ring

end DTWilliams.DiracWithFields
