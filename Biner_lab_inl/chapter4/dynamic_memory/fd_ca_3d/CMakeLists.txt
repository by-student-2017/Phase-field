cmake_minimum_required(VERSION 3.22.1 FATAL_ERROR)

# set the project name and version
project(MyProject VERSION 1.0)

add_executable(fd_ca_3d.exe fc_ca_3d.c
	init_grain_micro_3d.c free_energy_fd_ca_3d.c
	write_vtk_grid_values_3D.c)

target_compile_options(fd_ca_3d.exe PUBLIC -O2 -Wall)

#set(CMAKE_CXX_STANDARD 17)
target_compile_features(fd_ca_3d.exe PUBLIC cxx_std_17)

target_link_libraries(fd_ca_3d.exe -lm)

# usage
# cmake -S . -B build/ -G"Unix Makefiles"
# cmake --build build/ --target fd_ca_3d.exe
