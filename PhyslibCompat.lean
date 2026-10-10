/-
  Compatibility check: Physlib (leanprover-community/physlib, formerly
  PhysLean / HepLean) and every Dynamic-Theory module load together in one
  environment with no name collisions. It contains no declarations; it
  succeeds iff the imports below elaborate.
-/
import Physlib
import DynamicTheory
