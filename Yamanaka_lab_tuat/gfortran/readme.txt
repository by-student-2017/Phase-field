--------------------------------------------------------------------------
- Installation (fortran compiler and paraview)
1. sudo apt update
2. sudo apt -y install gfortran
3. sudo apt -y install paraview paraview-dev
--------------------------------------------------------------------------
Ref: https://web.tuat.ac.jp/~yamanaka/opensource.html
--------------------------------------------------------------------------
Phase-field model of isothermal solidification in a binary alloy
WBM model
J.A. Warren and W. J. Boettinger, Acta Metall. Mater., 43 (1995), p. 689.

1. gfortran pf_solidification.f90 -o pf_solidification
2. ./pf_solidification
3. paraview
  a1. File -> Open ... -> pf_result_..vtk -> OK
  a2. (click) [Apply] 
  a3. (Coloring) PhaseFields
  a4. (click "play") |>
--------------------------------------------------------------------------
Multi-phase-field model of polycrystalline grain growth

1. gfortran mpf_graingrowth.f90 -o mpf_graingrowth 
2. ./mpf_graingrowth
3. paraview
  a1. File -> Open ... -> mpf_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click "play") |>
--------------------------------------------------------------------------
Phase-field model of martensitic transformation

1. gfortran pf_martensite.f90 -o pf_martensite
2. ./pf_martensite
3. paraview
  a1. File -> Open ... -> mt_result_..vtk -> OK
  a2. (click) [Apply] 
  a3. (Coloring) PhaseFields
  a4. (click "play") |>
--------------------------------------------------------------------------