-----------------------------------------------------------------------
ubuntu 18.04 LTS

1. sudo apt update
2. sudo apt -y install g++
3. sudo apt -y install paraview paraview-dev
4. g++ all_elastic_field_v2.cpp -o all_elastic_field
  (gedit ph.dat (Only the data (ch) of the time you want to calculate.))
5. ./all_elastic_field
6. paraview
  a1. File -> Open ... -> sp_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click "play") |>
-----------------------------------------------------------------------
