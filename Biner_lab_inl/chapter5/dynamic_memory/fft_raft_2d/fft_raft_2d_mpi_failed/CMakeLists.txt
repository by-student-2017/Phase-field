cmake_minimum_required(VERSION 3.22.1 FATAL_ERROR)

# set the project name and version
project(MyProject VERSION 1.0)

find_package(MPI)

add_executable(fft_raft_2d_mpi.exe fft_raft_2d_mpi.c
  micro_ch_pre_2d_mpi.c prepare_fft_2d.c green_tensor_2d_mpi.c
  free_energy_ch_2d.c solve_elasticity_2d_mpi.c
  write_vtk_grid_values_2D.c)

target_compile_options(fft_raft_2d_mpi.exe PUBLIC -O2 -Wall)
#target_compile_options(fft_raft_2d_mpi.exe PUBLIC -g -O0 -Q)

#set(CMAKE_CXX_STANDARD 17)
target_compile_features(fft_raft_2d_mpi.exe PUBLIC cxx_std_17)

target_link_libraries(fft_raft_2d_mpi.exe -lm -lfftw3_mpi -lfftw3 ${MPI_LIBRARIES})
target_include_directories(fft_raft_2d_mpi.exe PUBLIC  ${MPI_INCLUDE_PATH})

# usage
# cmake -S . -B build/ -G"Unix Makefiles"
# cmake --build build/ --target fft_raft_2d_mpi.exe
