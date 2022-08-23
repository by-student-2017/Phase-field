-----------------------------------------------------------------------
ubuntu 18.04 LTS

1. sudo apt update
2. sudo apt -y install g++
3. sudo apt -y install paraview paraview-dev
4. g++ fecr3d_mag.cpp -o fecr3d_mag
5. ./fecr3d_mag
6. paraview
  a1. File -> Open ... -> sp_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click) Slice
  a4. (move red line)
  a5. (click) [Apply] 
  a6. (click "play") |>
  a7. (click "Rescale to Data Range")
  a8. (click "play") |>
-----------------------------------------------------------------------
----- dimensionless -----
Energy: 1/RT [1/(J/mol)]
Length: 1/b1 = 1/(L/N) = 1/([m]/(number of mesh))
Time: 1/(b1*b1/D), diffusion coefficient D [m2/s]
Time: Ms*RT, Ms [mol/s/J]=[(1/s)*1/(J/mol)], (mobility Mc=Ms=1)
Mc*RT/D [dimensionless]
Concentration gradient energy, kappa: 1/(b1*b1*RT)

Einstein relations: D = M0*RT, M0 [(m2/s)/(J/mol)]
-----------------------------------------------------------------------