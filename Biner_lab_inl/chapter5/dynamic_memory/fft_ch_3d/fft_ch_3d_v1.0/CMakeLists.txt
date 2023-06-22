cmake_minimum_required(VERSION 3.22.1 FATAL_ERROR)

# set the project name and version
project(MyProject VERSION 1.0)

add_executable(fft_ch_3d.exe
	fft_ch_3d.c
	micro_ch_pre_3d.c
	prepare_fft_3d.c
	free_energy_ch_3d.c
	calculate_energy_3d.c
	write_vtk_grid_values_3D.c
)

target_compile_options(fft_ch_3d.exe PUBLIC -O2 -Wall)
#target_compile_options(fft_ch_3d.exe PUBLIC -g -O0 -Q)

#set(CMAKE_CXX_STANDARD 17)
target_compile_features(fft_ch_3d.exe PUBLIC cxx_std_17)

target_link_libraries(fft_ch_3d.exe -lm -lfftw3)

# usage
# cmake -S . -B build/ -G"Unix Makefiles"
# cmake --build build/ --target fft_ch_3d.exe
