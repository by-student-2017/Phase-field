-----------------------------------------------------------------------
ubuntu 20.04 LTS
(underconstructing)

1. sudo apt update
2. sudo apt -y install g++ libfftw3-dev
3. sudo apt -y install paraview paraview-dev
4. g++ NGmicrostructure3D_v3_libfftw3.cpp -lfftw3 -o ngms
5. ./ngms
6. paraview
  a1. File -> Open ... -> NGms_3D_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click) Clip1
  a4. (move red line to origin (0 0 0))
  a5. (click) [Apply] 
  a6. (click "play") |>
  a7. (click "Rescale to Data Range")
  a8. (click "play") |>
-----------------------------------------------------------------------
