-----------------------------------------------------------------------
Windows11 and WSL2(Ubuntu-22.04)

(On WSL2(Ubuntu-22.04)
1. sudo apt update
2. sudo apt -y install gcc build-essential libfftw3-dev
3. sudo apt -y install libopenmpi-dev libfftw3-mpi-dev
4. cmake -S . -B build/ -G"Unix Makefiles"
5. cmake --build build/ --target pfc_3d_mpi.exe
6. cd ./build
7. export OMP_NUM_THREADS=2
8. mpirun -np 4 ./pfc_3d_mpi.exe

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