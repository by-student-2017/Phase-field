cmake_minimum_required(VERSION 3.22.1 FATAL_ERROR)

# set the project name and version
project(MyProject VERSION 1.0)

add_executable(fd_sint_2d.exe fd_sint_2d.c
	micro_sint_pre_2d.c
	free_energy_sint_2d_con.c free_energy_sint_2d_eta.c
	write_vtk_grid_values_2D.c)

target_compile_options(fd_sint_2d.exe PUBLIC -O2 -Wall)

#set(CMAKE_CXX_STANDARD 17)
target_compile_features(fd_sint_2d.exe PUBLIC cxx_std_17)

target_link_libraries(fd_sint_2d.exe -lm)

# usage
# cmake -S . -B build/ -G"Unix Makefiles"
# cmake --build build/ --target fd_sint_2d.exe
