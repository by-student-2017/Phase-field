-----------------------------------------------------------------------
ubuntu 20.04 LTS
(development version)
Spinodal decomposition (diffusive phase separation in A-B-C-D quaternary alloys)

1. sudo apt update
2. sudo apt -y install g++
3. sudo apt -y install paraview paraview-dev
4. g++ four_phase_v2.cpp -o four_phase
5. ./four_phase
6. paraview
  a1. File -> Open ... -> sp_result..vtk -> OK
  a2. [Solid Color] -> [concentration_Total]
  a3. [Outline] -> [Surface]
  a4. (click "play") |>
  a5. (click "Rescale to Data Range")
-----------------------------------------------------------------------
