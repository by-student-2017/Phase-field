-----------------------------------------------------------------------
ubuntu 18.04 LTS

1. sudo apt update
2. sudo apt -y install g++
3. sudo apt -y install paraview paraview-dev
4. g++ all_elastic_field_v2.cpp -o all_elastic_field
  (gedit ph.dat (Only the data (ch) of the time you want to calculate.) from output of Chapter5_Section1_two_phase_v2)
5. ./all_elastic_field
6. paraview
  a1. File -> Open ... -> pf_result000000.vtk -> OK
  a2. (click) [Apply] 
  a3. (click) [Glyph]
  a4. Scale Array [No scale array]
  a5. Scale Factor [5]
  a6. (click) [Apply]
-----------------------------------------------------------------------
ubuntu 18.04 LTS

1. sudo apt update
2. sudo apt -y install g++
3. sudo apt -y install paraview paraview-dev
4. g++ all_elastic_field_v2_libfftw3.cpp -lfftw3 -o all_elastic_field
  (gedit ph.dat (Only the data (ch) of the time you want to calculate.) from output of Chapter5_Section1_two_phase_v2)
5. ./all_elastic_field
6. paraview
  a1. File -> Open ... -> pf_result000000.vtk -> OK
  a2. (click) [Apply] 
  a3. (click) [Glyph]
  a4. Scale Array [No scale array]
  a5. Scale Factor [5]
  a6. (click) [Apply]
-----------------------------------------------------------------------
