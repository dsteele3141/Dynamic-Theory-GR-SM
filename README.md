# Dynamic Theory: the gravitational and particle sectors

Code accompanying the manuscript *Recovery of the gravitational and particle
sectors in Dynamic Theory*. It contains only the derivations
and checks that the manuscript relies on.

## How to read this

- `Layer2/` holds the five-dimensional foundations: the 5D metric, the five-force,
  the 5D Dirac ansatz and its coupling to fields.
- `Layer3/` holds the four-dimensional sectors reached under conservation of mass
  (u4 = 0) and the quantum condition: general relativity, the Standard Model, the
  comparison with data, and the predictions and bounds.
- The three laws themselves are adopted, as in the manuscript, and are not coded here.
- The Lean files check algebraic steps and the consistency of the adopted ansatz
  (for example, that the 5D Dirac factorisation has a matrix realisation and none
  in a commutative ring). They do not test the physics; comparison with experiment
  is done in the Sage scripts.
- Each script header names the step it supports in the Supplementary Materials,
  e.g. `Supplementary Materials: step [G2], Sec. S3`. The same tag is shown after
  each file below (F = Sec. S2, G = S3, S = S4, P = S5; code index in Sec. S6).

## Tools

- **SageMath** (`Layer*/Sage/*.sage`): assertion-driven derivations; every check
  prints PASS or the script stops. Stored output in `Layer*/sage-results/`.
- **Lean 4 + Mathlib + Physlib** (`Layer*/Lean/*.lean`): machine-checked theorems.
  `AxiomCheck.lean` reports the axiom footprint of every `DTWilliams` declaration.

## Running

```
lake build                          # Lean (toolchain pinned in lean-toolchain)
lake env lean AxiomCheck.lean       # axiom footprint: 0 sorry, standard axioms only
cd Layer2 && sage Sage/<name>.sage > sage-results/<name>.out.txt   # Layer2 scripts
cd Layer3 && sage Sage/<name>.sage > sage-results/<name>.out.txt   # Layer3 scripts
python3 linter.py                   # consistency sweep over all stored outputs
```

Assumptions beyond the three laws are tabulated in `ASSUMPTIONS.md`.

## Contents

### The 5D manifold, five-force and 5D Dirac ansatz (manuscript Sec. 2.2, 4.1)

- `Layer2/Sage/T2_1_5D_minkowski.sage` [F1]: 5D metric of signature (+,-,-,-,-) from the arc length; unit 5-velocity; equations of motion
- `Layer2/Lean/T2_1_Minkowski5D.lean` [F1]: the same steps in Lean: signature, determinant, velocity normalisation
- `Layer2/Lean/PhyslibBridge.lean` [F1, F3]: the 5D metric, Dirac matrices and Pauli matrices identified with Physlib's definitions
- `Layer2/Sage/force_density.sage` [F2]: the five-force density as the divergence of the 5D energy-momentum tensor, Div(T) = F.J
- `Layer2/Lean/T2_4_ForceDensity.lean` [F2]: the force-density algebra in Lean
- `Layer2/Sage/dirac_algebra.sage` [F3]: Clifford algebra of the five gamma matrices; the first-order factorisation of the 5D mass shell
- `Layer2/Lean/T2_7_DiracFactorisation.lean` [F3]: the factorisation has an 8x8 (and a 4x4) matrix realisation and none with commuting entries, so it is an adopted ansatz
- `Layer2/Lean/T2_7_PositivityRoute.lean` [F3]: whether positivity selects the first-order (Dirac) equation over the second-order one
- `Layer2/Sage/dirac_with_fields.sage` [F3]: the 5D Dirac equation with electromagnetic fields; the magnetic-moment coupling
- `Layer2/Lean/T2_8_DiracWithFields.lean` [F3]: the field coupling in Lean

### Conservation of mass, quantum condition, 4D surface, weak field (manuscript Sec. 2.3, 3)

- `Layer3/Sage/gateway_london_quantisation.sage` [F4]: the quantum condition (London): trivial holonomy for a quantised loop integral; Bohr radii
- `Layer3/Lean/T2_10_London.lean` [F4]: the same in Lean, with the scope of each step stated
- `Layer3/Sage/T2_11_hypersurface_metric.sage` [F6]: decomposition of the 5D metric into a 4D surface metric plus the mass direction
- `Layer3/Sage/u4_rest_mass_nonconservation.sage` [F5]: size of u4 (rate of failure of mass conservation): about 1e-31 for laboratory matter
- `Layer3/Sage/weak_field_poisson.sage` [F6]: weak-field limit: Poisson's equation with Newton's sign and coefficient, and the extra relation it needs

### Recovery of general relativity (manuscript Sec. 3)

- `Layer3/Sage/gr_on_4d_surface.sage` [G2]: Gauss-Codazzi reduction: with u4 = 0 the 4D surface obeys Einstein's equations
- `Layer3/Sage/einstein_equations_recovery.sage` [G2]: Lovelock's theorem applied: the surface field equations are Einstein's
- `Layer3/Sage/einstein_from_dt_action.sage` [G1]: variation of the Weyl-type action gives Einstein's equations plus a cosmological term
- `Layer3/Sage/weyl_lambda_reading.sage` [G1]: the cosmological constant on Weyl's procedure: |Lambda| = lambda/4, magnitude free
- `Layer3/Sage/weyl_vector_mass.sage` [G1]: the Weyl vector acquires a large mass and decouples, leaving G + Lambda g = kappa T
- `Layer3/Sage/weyl_ghilencea_comparison.sage` [G1]: term-by-term comparison with Weyl quadratic gravity (Ghilencea 2019/2020)
- `Layer3/Sage/kk_radion_g44.sage` [G3]: the g_44 (radion) mode: a free g_44 gives PPN gamma = 1/2; Cassini selects g_44 = -1
- `Layer3/Sage/test_body_motion.sage` [G4]: test bodies move on 4D geodesics; the fifth force vanishes at u4 = 0 (weak equivalence principle)
- `Layer3/Sage/matter_source_tab.sage` [G5]: the 5D matter tensor reduces to the GR source; the surface Bianchi identity requires u4 = 0
- `Layer3/Sage/gravitational_radiation.sage` [G6]: two polarisations at speed c; quadrupole power against the binary pulsar PSR B1913+16
- `Layer3/Sage/exact_ppwave_lift.sage` [G6]: exact pp-waves lift to exact 5D vacuum solutions
- `Layer3/Sage/black_hole_ringdown.sage` [G6]: strong field: ringdown and EHT select the GR (Kerr) solution over the non-singular metric
- `Layer3/Sage/gr_limit_domain.sage` [G7]: where the Newtonian/GR limit holds; the rotation-curve term has the right shape but is far too small

### Recovery of the Standard Model (manuscript Sec. 4)

- `Layer3/Sage/weak_parity.sage` [S1]: a 5D Dirac sector without a boundary is parity-even (the Kaluza-Klein chirality problem)
- `Layer3/Sage/weak_chirality_u4.sage` [S1]: u4 = 0 as a no-flux edge on the half-line leaves a single chiral zero mode
- `Layer3/Sage/noflux_boundary_derivation.sage` [S1]: the no-flux edge condition derived from unitarity and 4D Lorentz invariance
- `Layer3/Sage/weak_va_coupling.sage` [S2]: the V-A form of the charged current follows from the chiral zero modes
- `Layer3/Sage/su2_from_5d.sage` [S3]: the 8x8 representation carries su(2) + u(1), giving the lepton doublet with Y = -1/2
- `Layer3/Sage/weak_neutral_current.sage` [S4]: the neutral current: adopted from the SM, consistent with the data
- `Layer3/Sage/mass_term_su2.sage` [S4]: compatibility of the 5D mass term with SU(2)_L; fermion masses from Yukawa couplings (adopted)
- `Layer3/Sage/su3_larger_rep.sage` [S5]: a field of N blocks carries exactly u(N), which survives u4 = 0; N = 3 for colour is selected
- `Layer3/Sage/strong_sector_su3.sage` [S5]: odd N_c and fractional charges; with q = Ze, confinement acts as a selection rule
- `Layer3/Sage/gauge_fields_halfline.sage` [S6]: gauge fields on the half-line: boundary pairing, no light adjoint scalar, A_4 decouples
- `Layer3/Sage/beta_spectrum_two_body.sage` [S7]: a two-body neutron decay gives a line, not the observed continuous beta spectrum
- `Layer3/Sage/neutrino_helicity.sage` [S7]: the neutrino is the chiral zero mode (spin 1/2, left-handed), as in the SM
- `Layer3/Sage/neutron_lifetime_tunnelling.sage` [S8]: the tunnelling model of the neutron lifetime against the measured lifetime
- `Layer3/Sage/hydrogen_nonsingular_spectrum.sage` [S8]: the non-singular potential in hydrogen and muonic hydrogen bounds its range parameter

### Predictions and bounds from departures from conservation of mass and the quantum condition (manuscript Sec. 6)

- `Layer3/Sage/charge_quantisation_prediction.sage` [P1]: q = Ze fixes the neutron and neutrino charges at exactly zero, also with right-handed neutrinos
- `Layer3/Sage/heavy_massaxis_partners.sage` [P2]: heavy partners along the mass axis: a continuum above a free scale M, bounded at about 1 TeV
- `Layer3/Sage/blackett_coupling_bodies.sage` [P3]: the magnetic moment of a spinning neutral body: excluded as a magnetic field, read as GR frame dragging
- `Layer3/Sage/T2_2_charge_conservation.sage` [P4]: the modified charge-conservation law with the fifth current component J4
- `Layer3/Sage/u4_relaxation_bounds.sage` [P4]: departures from u4 = 0: 4D charge conserved exactly; the breathing mode bounded by Cassini
- `Layer3/Sage/pantheon_dark_energy_fit.sage` [P5]: fit of the secular drift to the Pantheon supernova Hubble diagram
- `Layer3/Sage/dt_lambda_joint_fit.sage` [P5]: joint fit of the drift and a cosmological constant to the supernova data
- `Layer3/Sage/pantheon_sys_DS17f.txt` (data): Pantheon DS17f systematic covariance (Scolnic et al. 2018, https://github.com/dscolnic/Pantheon)
- `Layer3/Sage/hbar_potential_clocks.sage` [P6]: a potential-dependent Planck constant against atomic-clock tests

## Licence

See `LICENSE`.
