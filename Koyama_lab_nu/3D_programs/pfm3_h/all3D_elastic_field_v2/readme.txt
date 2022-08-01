-----------------------------------------------------------------------
ubuntu 18.04 LTS

1. sudo apt update
2. sudo apt -y install g++
3. sudo apt -y install paraview paraview-dev
4. mv *.dat data.dat
  (e.g., *.dat from microstructure3D_v2.cpp)
5. g++ all3D_elastic_field_v2.cpp -o el
6. ./el
7. paraview
  a1. File -> Open ... -> el_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click) Clip1
  a4. (move red line to origin (0 0 0))
  a5. (click) [Apply] 
  a6. (click "play") |>
  a7. (click "Rescale to Data Range")
  a8. (click "play") |>
-----------------------------------------------------------------------
