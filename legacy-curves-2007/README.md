# Legacy four-body 1D curves and couplings (2007)

The oldest surviving version of the four-body 1D calculation in this repository's lineage.
**Not an SVD calculation** — it computes the adiabatic curves and the couplings directly — but it
is kept here because that is what it produces, and this is where that capability lives.

Recovered 2026-09-23 from `~/Documents/Code/4body200`.

## Why only one copy

`~/Documents/Code` held three directories — `4body100`, `4body200`, `4body300` — whose **source
files are byte-identical** (all 2007-08-23). They differed only in run parameters: `FourBody1D.par`
asks for 10 energies in the `100` variant and 40 in the other two, with identical `adcurve.par`.
The rest of each directory was output. This is the `4body200` copy.

## Relationship to the other versions here

| | |
|---|---|
| `legacy-curves-2007/` (this directory) | 2007-08-23 — earliest |
| `../legacy-finite-difference/` | `FourBody1D.f` 2018-01-09, `FourBodyPOT.f` 2018-01-12 — the same program eleven years on, with the finite-difference coupling routines documented there |
| `../4BodySVD.f90`, `../adiabaticSolver.f90` | the current SVD approach |

## Contents

Sources `FourBody1D.f`, `FourBody1D_Bsplines.f`, `FourBody1D_matrix_stuff.f`, `FourBodyPOT.f`,
`1DPot.f`; build and input files `FourBody1D.mak`, `FourBody1D.par`, `adcurve.par`; data
`Legendre.dat` (Gauss-Legendre nodes/weights) and `index18.dat`; and two representative Grace plots
of the resulting curves, `adcurve_a2_p100.agr` and `chris_4body_odd_fermions.agr`.
