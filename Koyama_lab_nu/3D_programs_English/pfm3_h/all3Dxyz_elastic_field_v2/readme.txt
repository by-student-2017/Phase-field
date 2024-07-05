-----------------------------------------------------------------------
ubuntu 18.04 LTS

1. sudo apt update
2. sudo apt -y install g++ libfftw3-dev
3. sudo apt -y install paraview paraview-dev
4. mv *.dat data.dat
  (e.g., *.dat (Only the data of the time you want to calculate.) from 3Dxyz_APT_MPF_v2.cpp, etc)
5. g++ all3Dxyz_elastic_field_v2_libfftw3.cpp -lfftw3 -o el
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
8. paraview
  a1. File -> Open ... -> el_result..vtk -> OK
  a2. (click) [Apply] 
  a3. (click) [Glyph]
  a4. Scale Array [No scale array]
  a5. Scale Factor [5]
  a6. (click) [Apply]
Note: Transparency; Properties -> Opacity
-----------------------------------------------------------------------
