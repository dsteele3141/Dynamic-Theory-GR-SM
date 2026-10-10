# Dynamic Theory: General Relativity and the Standard Model

Code accompanying the manuscript *Dynamic Theory: recovery of General Relativity
and the Standard Model from conservation of mass*. It contains only the derivations
and checks that the manuscript relies on.

## Tools

- **SageMath** (`Layer*/Sage/*.sage`): assertion-driven derivations; every check
  prints PASS or the script stops. Stored output in `Layer*/sage-results/`.
- **Cadabra2** (`Layer*/Cadabra/*.cdb`): independent tensor-algebra cross-checks;
  output in `Layer*/cadabra-results/`.
- **Lean 4 + Mathlib + PhysLib** (`Layer*/Lean/*.lean`): machine-checked theorems.
  `AxiomCheck.lean` reports the axiom footprint of every `DTWilliams` declaration.

## Running

```
lake build                          # Lean (toolchain pinned in lean-toolchain)
lake env lean AxiomCheck.lean       # axiom footprint: 0 sorry, standard axioms only
cd Layer3 && sage Sage/<name>.sage > sage-results/<name>.out.txt
python3 linter.py                   # consistency sweep over all stored outputs
```

Assumptions beyond the three laws are tabulated in `ASSUMPTIONS.md`.

## Contents

### The 5D manifold, five-force and 5D Dirac ansatz (manuscript Sec. 2.2, 4.1)

- `Layer2/Sage/T2_1_5D_minkowski.sage`
- `Layer2/Lean/T2_1_Minkowski5D.lean`
- `Layer2/Lean/PhyslibBridge.lean`
- `Layer2/Sage/force_density.sage`
- `Layer2/Lean/T2_4_ForceDensity.lean`
- `Layer2/Sage/dirac_algebra.sage`
- `Layer2/Lean/T2_7_DiracFactorisation.lean`
- `Layer2/Lean/T2_7_PositivityRoute.lean`
- `Layer2/Sage/dirac_with_fields.sage`
- `Layer2/Lean/T2_8_DiracWithFields.lean`

### Layer boundary and Layer 3 inputs: quantum condition, u4 = 0, 4D surface, weak field (manuscript Sec. 2.3, 3)

- `Layer3/Sage/gateway_london_quantisation.sage`
- `Layer3/Lean/T2_10_London.lean`
- `Layer3/Sage/T2_11_hypersurface_metric.sage`
- `Layer3/Sage/u4_rest_mass_nonconservation.sage`
- `Layer3/Sage/weak_field_poisson.sage`

### Recovery of general relativity (manuscript Sec. 3)

- `Layer3/Sage/gr_on_4d_surface.sage`
- `Layer3/Sage/einstein_equations_recovery.sage`
- `Layer3/Sage/einstein_from_dt_action.sage`
- `Layer3/Sage/weyl_lambda_reading.sage`
- `Layer3/Sage/weyl_vector_mass.sage`
- `Layer3/Sage/weyl_ghilencea_comparison.sage`
- `Layer3/Sage/kk_radion_g44.sage`
- `Layer3/Sage/test_body_motion.sage`
- `Layer3/Sage/matter_source_tab.sage`
- `Layer3/Sage/gravitational_radiation.sage`
- `Layer3/Sage/exact_ppwave_lift.sage`
- `Layer3/Sage/black_hole_ringdown.sage`
- `Layer3/Sage/gr_limit_domain.sage`

### Recovery of the Standard Model (manuscript Sec. 4)

- `Layer3/Sage/weak_parity.sage`
- `Layer3/Sage/weak_chirality_u4.sage`
- `Layer3/Sage/noflux_boundary_derivation.sage`
- `Layer3/Sage/weak_va_coupling.sage`
- `Layer3/Sage/su2_from_5d.sage`
- `Layer3/Sage/weak_neutral_current.sage`
- `Layer3/Sage/mass_term_su2.sage`
- `Layer3/Sage/su3_larger_rep.sage`
- `Layer3/Sage/strong_sector_su3.sage`
- `Layer3/Sage/gauge_fields_halfline.sage`
- `Layer3/Sage/beta_spectrum_two_body.sage`
- `Layer3/Sage/neutrino_helicity.sage`
- `Layer3/Sage/neutron_lifetime_tunnelling.sage`
- `Layer3/Sage/hydrogen_nonsingular_spectrum.sage`

### Predictions and bounds from departures from the restrictive assumptions (manuscript Sec. 6)

- `Layer3/Sage/charge_quantisation_prediction.sage`
- `Layer3/Sage/heavy_massaxis_partners.sage`
- `Layer3/Sage/blackett_coupling_bodies.sage`
- `Layer3/Sage/T2_2_charge_conservation.sage`
- `Layer3/Sage/u4_relaxation_bounds.sage`
- `Layer3/Sage/pantheon_dark_energy_fit.sage`
- `Layer3/Sage/dt_lambda_joint_fit.sage`
- `Layer3/Sage/pantheon_sys_DS17f.txt`
- `Layer3/Sage/hbar_potential_clocks.sage`

## Licence

See `LICENSE`.
