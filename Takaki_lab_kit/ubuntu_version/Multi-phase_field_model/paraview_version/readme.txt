ubuntu 18.04 LTS
-----------------------------------------------------------------------
1. sudo apt update
2. sudo apt -y install gfortran
3. sudo apt -y install paraview paraview-dev
4. gfortran ch14_mpf_paraview.f -o mpf
5. ./mpf
6. paraview
  a1. File -> Open ... -> mpf_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click "play") |>
-----------------------------------------------------------------------