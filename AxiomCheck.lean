-- --------------------------------------------------------------------------
-- SOURCE (primary): P. E. Williams, "The Dynamic Theory: A New View of
--   Space, Time, and Matter," LA-8370-MS (LASL, 1980). All equation/page
--   numbers herein refer to this Report. A book ("The Dynamic Theory";
--   "Physics - Against the Odds") is cited ONLY where the Report is
--   insufficient and another source is necessary (e.g. the reversed
--   gauge-sign correction, per Dr. Dan Ross).
-- --------------------------------------------------------------------------

/-
  Axiom-footprint audit for the DTWilliams Lean layer.
  Runs `collectAxioms` (the engine behind `#print axioms`) over every theorem
  and def under the `DTWilliams` namespace, reports the union of axioms used,
  and flags any `sorryAx` or any axiom outside Lean's three standard axioms
  {propext, Classical.choice, Quot.sound}.

  This is the authoritative copy. To run it, copy it into the lake build
  environment alongside the DTWilliams sources and execute:
      lake env lean AxiomCheck.lean
-/
import Basic
import Layer2.Lean.PhyslibBridge
import Layer2.Lean.T2_1_Minkowski5D
import Layer2.Lean.T2_4_ForceDensity
import Layer2.Lean.T2_7_DiracFactorisation
import Layer2.Lean.T2_7_PositivityRoute
import Layer2.Lean.T2_8_DiracWithFields
import Layer3.Lean.T2_10_London

open Lean Elab Command

run_cmd do
  let env ← getEnv
  let pfx : Name := `DTWilliams
  let standard : List Name := [`propext, `Classical.choice, `Quot.sound]
  -- candidate user theorems/defs under `DTWilliams` (skip compiler-internal names)
  let names := env.constants.fold (init := (#[] : Array Name)) (fun acc n ci =>
    let keep := match ci with
      | .thmInfo _  => true
      | .defnInfo _ => true
      | _           => false
    if keep && pfx.isPrefixOf n && !n.isInternal then acc.push n else acc)
  let mut allAx : Array Name := #[]
  let mut sorryNames : Array Name := #[]
  for n in names do
    let ax ← liftCoreM (Lean.collectAxioms n)
    for a in ax do
      if !allAx.contains a then allAx := allAx.push a
    if ax.contains `sorryAx then
      sorryNames := sorryNames.push n
  let nonStd := allAx.filter (fun a => !standard.contains a)
  IO.println "================================================================"
  IO.println "DTWilliams AXIOM FOOTPRINT  (collectAxioms over all theorems/defs)"
  IO.println "================================================================"
  IO.println s!"Theorems + defs scanned under `DTWilliams`: {names.size}"
  IO.println s!"Axioms used (union): {allAx.toList}"
  IO.println s!"Declarations depending on sorryAx: {sorryNames.size}"
  for n in sorryNames do
    IO.println s!"  sorryAx -> {n}"
  if nonStd.isEmpty then
    IO.println "Axioms beyond {propext, Classical.choice, Quot.sound}: NONE"
    IO.println "RESULT: PASS -- every DTWilliams theorem rests only on Lean's"
    IO.println "        standard axioms; no sorryAx anywhere in the footprint."
  else
    IO.println s!"Axioms beyond the standard three: {nonStd.toList}"
    IO.println "RESULT: REVIEW -- non-standard axioms in the footprint (see above)."
  IO.println "================================================================"
