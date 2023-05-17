cmake_minimum_required(VERSION 3.22.1 FATAL_ERROR)

# set the project name and version
project(MyProject VERSION 1.0)

add_executable(fem_elast_2d.exe
	fem_elast_2d.c
	input_fem_elast_2d.c
	gauss_2d.c
	stiffness_2d.c
		sfr2_2d.c
		jacob2_2d.c
		modps_2d.c
		bmats_2d.c
		dbe_2d.c
	loads_2d.c
	boundary_cond_2d.c
	stress_2d.c
	output_fem_2d.c
		write_vtk_fem_2d.c
	)

target_compile_options(fem_elast_2d.exe PUBLIC -O2 -Wall)

#set(CMAKE_CXX_STANDARD 17)
target_compile_features(fem_elast_2d.exe PUBLIC cxx_std_17)

target_link_libraries(fem_elast_2d.exe -lm -llapack)

# usage
# cmake -S . -B build/ -G"Unix Makefiles"
# cmake --build build/ --target fem_elast_2d.exe
