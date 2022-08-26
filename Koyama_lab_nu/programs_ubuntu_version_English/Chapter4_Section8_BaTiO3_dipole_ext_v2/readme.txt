-----------------------------------------------------------------------
ubuntu 18.04 LTS

1. sudo apt update
2. sudo apt -y install g++
3. sudo apt -y install paraview paraview-dev
4. sudo apt -y install unzip
5. unzip test.zip
6. g++ BaTiO3_dipole_ext_v2.cpp -o BaTiO3_dipole_ext
7. ./BaTiO3_dipole_ext
8. paraview
  a1. File -> Open ... -> dip_ext_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click "play") |>
-----------------------------------------------------------------------
ubuntu 18.04 LTS

1. sudo apt update
2. sudo apt -y install g++ libfftw3-dev
3. sudo apt -y install paraview paraview-dev
4. sudo apt -y install unzip
5. unzip test.zip
6. g++ BaTiO3_dipole_ext_v2_libfftw3.cpp -lfftw3 -o BaTiO3_dipole_ext
7. ./BaTiO3_dipole_ext
8. paraview
  a1. File -> Open ... -> dip_ext_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click "play") |>
-----------------------------------------------------------------------
