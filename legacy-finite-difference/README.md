# Legacy four-body 1D calculation (finite-difference couplings)

Pre-SVD version of the four-body 1D adiabatic hyperspherical calculation, kept for reference.
Recovered 2026-09-22 from `~/Documents/Code/KH-FourBody1D.tar` before that archive was deleted.

## What it does

- B-spline basis in the two hyperangles (`xPoints`, `yPoints`), banded overlap and Hamiltonian
  matrices, eigenvalues via ARPACK `dsband`.
- **Couplings by three-point finite difference in R**: `CalcPMatrix` forms `rPsi - lPsi` over a
  spacing `RDelt`, applies the overlap with `dsbmv` and contracts with `mPsi`, giving
  `P = <psi_i | d psi_j / dR>`. `CalcQMatrix` does the second derivative the same way.
  `FixPhase` keeps eigenvector signs consistent from one R to the next.
- Four hand-tuned hyperangular grid generators (`GridMaker`, `GridMaker111`, `GridMaker222`,
  `GridMakerBetter`).

## Why it is here and not in use

The parent repository supersedes it: `4BodySVD.f90` / `adiabaticSolver.f90` solve the same
four-body 1D problem, and `FourBodyPOT.f90` is the modernized descendant of the `FourBodyPOT.f`
kept here. This directory is retained only because the finite-difference coupling machinery does
not otherwise survive, and because these 2017 files were newer than anything in `~/Developer`
(`FourBody1D.f` by four years, `FourBodyPOT.f` by ten).

For contrast, the 1D three-body delta-function problem computes the same P and Q matrices in
**closed form** — see the `Delta-Function-Recombination` repository.

Library dependencies (`FourBody1D_Bsplines.f`, `FourBody1D_matrix_stuff.f`) were not copied;
equivalents live in `~/Developer/lib` and in this repository.
