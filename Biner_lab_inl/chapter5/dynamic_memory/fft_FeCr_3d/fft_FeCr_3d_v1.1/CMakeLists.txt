cmake_minimum_required(VERSION 3.22.1 FATAL_ERROR)

# set the project name and version
project(MyProject VERSION 1.0)

add_executable(fft_FeCr_3d.exe fft_FeCr_3d.c
  dislo_strain_3d.c FeCr_chem_poten_3d.c
  green_tensor1_3D.c green_tensor2_3D.c
  prepare_fft_3d.c solve_elasticity_3d.c
  micro_ch_pre_3d.c write_vtk_grid_values_3D.c)

target_compile_options(fft_FeCr_3d.exe PUBLIC -O2 -Wall)
#target_compile_options(fft_FeCr_3d.exe PUBLIC -g -O0 -Q)

#set(CMAKE_CXX_STANDARD 17)
target_compile_features(fft_FeCr_3d.exe PUBLIC cxx_std_17)

target_link_libraries(fft_FeCr_3d.exe -lm -lfftw3)

# usage
# cmake -S . -B build/ -G"Unix Makefiles"
# cmake --build build/ --target fft_FeCr_3d.exe
