-----------------------------------------------------------------------
ubuntu 18.04 LTS

1. sudo apt update
2. sudo apt -y install g++
3. sudo apt -y install paraview paraview-dev
4. g++ phase_decomp_2D.cpp -o phase_decomp_2D
5. ./phase_decomp_2D
6. paraview
  a1. File -> Open ... -> sp_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click "play") |>
-----------------------------------------------------------------------
ubuntu 18.04 LTS

1. sudo apt update
2. sudo apt -y install g++
3. sudo apt -y install paraview paraview-dev
4. g++ fecr2d_demo.cpp -o fecr2d_demo
5. ./fecr2d_demo
6. paraview
  a1. File -> Open ... -> sp_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click "play") |>
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