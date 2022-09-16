-----------------------------------------------------------------------
ubuntu 18.04 LTS
(underconstructing)

1. sudo apt update
2. sudo apt -y install g++ libfftw3-dev
3. sudo apt -y install paraview paraview-dev
4. g++ BaTiO3_dipole_ext_v3_libfftw3.cpp -lfftw3 -o BaTiO3_dipole_ext
5. ./BaTiO3_dipole_ext
6. paraview
  a1. File -> Open ... -> dip_ext_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click "play") |>
-----------------------------------------------------------------------
