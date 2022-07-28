-----------------------------------------------------------------------
ubuntu 18.04 LTS

1. sudo apt update
2. sudo apt -y install g++
3. sudo apt -y install paraview paraview-dev
4. g++ lamella_microstructure3D_save.cpp -o lamella
5. ./lamella
6. paraview
  a1. File -> Open ... -> lamella_3D_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click) Slice
  a4. (move red line)
  a5. (click) [Apply] 
  a6. (click "play") |>
  a7. (click "Rescale to Data Range")
  a8. (click "play") |>
-----------------------------------------------------------------------