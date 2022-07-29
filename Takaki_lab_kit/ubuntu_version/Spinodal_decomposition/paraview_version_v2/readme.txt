check ubuntu 18.04 LTS
-----------------------------------------------------------------------
1. sudo apt update
2. sudo apt -y install g++
3. sudo apt -y install paraview paraview-dev
4. g++ main-cpu_v2.cpp -o main-cpu
5. ./main-cpu
6. paraview
  a1. File -> Open ... -> f..vtk -> OK
  a2. (click) [Apply] 
  a3. (click "play") |>
-----------------------------------------------------------------------