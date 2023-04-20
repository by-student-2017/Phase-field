cmake_minimum_required(VERSION 3.22.1 FATAL_ERROR)

# set the project name and version
project(MyProject VERSION 1.0)

add_executable(fft_FeCuMnNi_2d.exe fft_FeCuMnNi_2d.c init_FeCuMnNi_micro_2d.c
  prepare_fft_2d.c FeCuMnNi_free_energy_2d.c write_vtk_grid_values_2D.c)

target_compile_options(fft_FeCuMnNi_2d.exe PUBLIC -O2 -Wall)

#set(CMAKE_CXX_STANDARD 17)
target_compile_features(fft_FeCuMnNi_2d.exe PUBLIC cxx_std_17)

target_link_libraries(fft_FeCuMnNi_2d.exe -lm -lfftw3)

# usage
# cmake -S . -B build/ -G"Unix Makefiles"
# cmake --build build/ --target fft_FeCuMnNi_2d.exe
