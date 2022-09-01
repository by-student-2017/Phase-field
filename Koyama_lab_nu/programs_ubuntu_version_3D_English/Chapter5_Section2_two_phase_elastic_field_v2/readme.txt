-----------------------------------------------------------------------
ubuntu 18.04 LTS

1. sudo apt update
2. sudo apt -y install g++ libfftw3-dev
3. sudo apt -y install paraview paraview-dev
4. g++ all_elastic_field_v2_libfftw3.cpp -lfftw3 -o all_elastic_field
4. mv *.dat data.dat
  (gedit *.dat (Only the data (ch) of the time you want to calculate.) from output of Chapter5_Section1_two_phase_v2)
6. ./all_elastic_field
7. paraview
  a1. File -> Open ... -> pf_result-000001.vtk -> OK
  a2. (click) [Apply] 
  a3. (click "play") |>
-----------------------------------------------------------------------
