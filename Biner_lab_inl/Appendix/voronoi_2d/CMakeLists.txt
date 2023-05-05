cmake_minimum_required(VERSION 3.22.1 FATAL_ERROR)

# set the project name and version
project(MyProject VERSION 1.0)

add_executable(voronoi_2d.exe voronoi_2d.c
	Discrete_Voronoi_diagram_2d.c
	write_vtk_grid_values_2D.c )

target_compile_options(voronoi_2d.exe PUBLIC -O2 -Wall)

#set(CMAKE_CXX_STANDARD 17)
target_compile_features(voronoi_2d.exe PUBLIC cxx_std_17)

target_link_libraries(voronoi_2d.exe -lm)

# usage
# cmake -S . -B build/ -G"Unix Makefiles"
# cmake --build build/ --target voronoi_2d.exe
