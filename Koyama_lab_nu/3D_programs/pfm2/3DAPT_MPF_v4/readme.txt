-----------------------------------------------------------------------
ubuntu 18.04 LTS
(development version)

1. sudo apt update
2. sudo apt -y install g++
3. sudo apt -y install paraview paraview-dev
4. g++ 3DAPT_MPF_v4.cpp -o 3DAPT_MPF
5. ./3DAPT_MPF
6. paraview
  a1. File -> Open ... -> 3DAPT_MPF_N001_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click) Clip1
  a4. (move red line to origin (0 0 0))
  a5. (click) [Apply] 
  a6. (click "play") |>
  a7. (click "Rescale to Data Range")
  a8. (click "play") |>
-----------------------------------------------------------------------