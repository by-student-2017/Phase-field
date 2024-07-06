-----------------------------------------------------------------------
Note: Development version
ubuntu 20.04 LTS

1. sudo apt update
2. sudo apt -y install g++ libfftw3-dev
3. sudo apt -y install paraview paraview-dev
4. g++ two_phase_v3_libfftw3.cpp -lfftw3 -o two_phase
5. ./two_phase
6. paraview
  a1. File -> Open ... -> sp_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click "play") |>
7. paraview
  a1. File -> Open ... -> sp_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click) [Glyph]
  a4. Scale Array [No scale array]
  a5. Scale Factor [5]
  a6. (click) [Apply]
-----------------------------------------------------------------------
