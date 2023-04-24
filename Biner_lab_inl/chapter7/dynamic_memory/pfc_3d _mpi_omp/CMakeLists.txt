cmake_minimum_required(VERSION 3.22.1 FATAL_ERROR)

# set the project name and version
project(MyProject VERSION 1.0)

find_package(MPI)

add_executable(pfc_3d_mpi.exe pfc_3d_mpi.c prepare_fft_3d.c write_vtk_grid_values_3D.c)

target_compile_options(pfc_3d_mpi.exe PUBLIC -O2 -Wall)

#set(CMAKE_CXX_STANDARD 17)
target_compile_features(pfc_3d_mpi.exe PUBLIC cxx_std_17)

target_link_libraries(pfc_3d_mpi.exe -lfftw3_mpi -lfftw3_omp -lfftw3 -lm -lpthread -fopenmp ${MPI_LIBRARIES})
target_include_directories(pfc_3d_mpi.exe PUBLIC  ${MPI_INCLUDE_PATH})


# usage
# cmake -S . -B build/ -G"Unix Makefiles"
# cmake --build build/ --target pfc_3d_mpi.exe
