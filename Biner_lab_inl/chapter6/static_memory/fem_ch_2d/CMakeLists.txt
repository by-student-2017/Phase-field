cmake_minimum_required(VERSION 3.22.1 FATAL_ERROR)

# set the project name and version
project(MyProject VERSION 1.0)

add_executable(fem_ch_2d.exe
	fem_ch_2d.c
	periodic_boundary_2d.c
	gauss_2d.c
	cart_deriv_2d.c
	init_micro_ch_fem_2d.c
	chem_stiff_2d.c
	apply_periodic_bc_2d.c
	recover_slave_dof_2d.c
	write_vtk_fem_2d.c
	)

target_compile_options(fem_ch_2d.exe PUBLIC -O2 -Wall)

#set(CMAKE_CXX_STANDARD 17)
target_compile_features(fem_ch_2d.exe PUBLIC cxx_std_17)

target_link_libraries(fem_ch_2d.exe -lm -llapack)

# usage
# cmake -S . -B build/ -G"Unix Makefiles"
# cmake --build build/ --target fem_ch_2d.exe
