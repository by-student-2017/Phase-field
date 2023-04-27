cmake_minimum_required(VERSION 3.22.1 FATAL_ERROR)

# set the project name and version
project(MyProject VERSION 1.0)

find_package(MPI)

add_executable(fft_FeCuMnNi_3d_mpi.exe fft_FeCuMnNi_3d.c init_FeCuMnNi_micro_3d_mpi.c
  prepare_fft_3d.c FeCuMnNi_free_energy_3d.c write_vtk_grid_values_3d.c)

target_compile_options(fft_FeCuMnNi_3d_mpi.exe PUBLIC -O2 -Wall)

#set(CMAKE_CXX_STANDARD 17)
target_compile_features(fft_FeCuMnNi_3d_mpi.exe PUBLIC cxx_std_17)

target_link_libraries(fft_FeCuMnNi_3d_mpi.exe -lm -lfftw3_mpi -lfftw3 ${MPI_LIBRARIES})
target_include_directories(fft_FeCuMnNi_3d_mpi.exe PUBLIC  ${MPI_INCLUDE_PATH})

# usage
# cmake -S . -B build/ -G"Unix Makefiles"
# cmake --build build/ --target fft_FeCuMnNi_3d_mpi.exe
