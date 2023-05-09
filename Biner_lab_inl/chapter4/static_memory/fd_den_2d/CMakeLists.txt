cmake_minimum_required(VERSION 3.22.1 FATAL_ERROR)

# set the project name and version
project(MyProject VERSION 1.0)

add_executable(fd_den_2d.exe fd_den_2d.c
	nucleus_2d.c
	write_vtk_grid_values_2D.c)

target_compile_options(fd_den_2d.exe PUBLIC -O2 -Wall)

#set(CMAKE_CXX_STANDARD 17)
target_compile_features(fd_den_2d.exe PUBLIC cxx_std_17)

target_link_libraries(fd_den_2d.exe -lm)

# usage
# cmake -S . -B build/ -G"Unix Makefiles"
# cmake --build build/ --target fd_den_2d.exe
