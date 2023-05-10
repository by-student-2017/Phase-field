cmake_minimum_required(VERSION 3.22.1 FATAL_ERROR)

# set the project name and version
project(MyProject VERSION 1.0)

add_executable(fem_elast_2d.exe fc_ca_2d.c
	init_grain_micro_2d.c free_energy_fd_ca_2d.c
	write_vtk_grid_values_2D.c)

target_compile_options(fem_elast_2d.exe PUBLIC -O2 -Wall)

#set(CMAKE_CXX_STANDARD 17)
target_compile_features(fem_elast_2d.exe PUBLIC cxx_std_17)

target_link_libraries(fem_elast_2d.exe -lm)

# usage
# cmake -S . -B build/ -G"Unix Makefiles"
# cmake --build build/ --target fem_elast_2d.exe
