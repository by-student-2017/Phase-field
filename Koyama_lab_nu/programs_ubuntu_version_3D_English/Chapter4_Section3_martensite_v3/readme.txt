-----------------------------------------------------------------------
ubuntu 20.04 LTS
(underconstructing)

1. sudo apt update
2. sudo apt -y install g++ libfftw3-dev
3. sudo apt -y install paraview paraview-dev
4. g++ martensite_v3_libfftw3.cpp -lfftw3 -o martensite
5. ./martensite
6. paraview
  a1. File -> Open ... -> sp_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click "play") |>
-----------------------------------------------------------------------
