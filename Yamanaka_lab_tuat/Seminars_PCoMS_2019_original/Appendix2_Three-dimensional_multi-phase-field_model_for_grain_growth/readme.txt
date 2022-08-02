ubuntu 18.04 LTS
-----------------------------------------------------------------------
1. sudo apt update
2. sudo apt -y install python3 python3-numba python3-numpy
3. sudo apt -y install paraview paraview-dev
4. python3 multi_phase_field_for_grain_growth_3d.py
  (wait about 26 min) (you can get from p_3d0.vtk to p_3d2000.vtk)
5. paraview
  a1. File -> Open ... -> lamella_3D_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click) Clip1
  a4. (move red line to origin (0 0 0))
  a5. (click) [Apply] 
  a6. (click "play") |>
  a7. (click "Rescale to Data Range")
  a8. (click "play") |>