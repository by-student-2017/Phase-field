cmake_minimum_required(VERSION 3.22.1 FATAL_ERROR)

# set the project name and version
project(MyProject VERSION 1.0)

add_executable(pfc_def_2d.exe pfc_def_2d.c prepare_fft_2d.c write_vtk_grid_values_2D.c)

target_compile_options(pfc_def_2d.exe PUBLIC -O2 -Wall)

#set(CMAKE_CXX_STANDARD 17)
target_compile_features(pfc_def_2d.exe PUBLIC cxx_std_17)

target_link_libraries(pfc_def_2d.exe -lm -lfftw3)

# usage
# cmake -S . -B build/ -G"Unix Makefiles"
# cmake --build build/ --target pfc_def_2d.exe
