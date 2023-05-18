-----------------------------------------------------------------------
Windows11 and WSL2(Ubuntu-22.04)
(underconstructing)

(On WSL2(Ubuntu-22.04)
1. sudo apt update
2. sudo apt -y install gcc build-essential
3. cmake -S . -B build/ -G"Unix Makefiles"
4. cmake --build build/ --target fem_elast_2d.exe
5. cd ./build
6. cp ../mesh_1.inp ./
7. ./fem_elast_2d.exe

(On windows11)
6. paraview
  a1. File -> Open ... -> time_..vtk -> OK
  a2. (click) [Apply] 
  a3. (click) Clip1
  a4. (move red line to origin (0 0 0))
  a5. (click) [Apply] 
  a6. (click "play") |>
  a7. (click "Rescale to Data Range")
  a8. (click "play") |>
-----------------------------------------------------------------------
