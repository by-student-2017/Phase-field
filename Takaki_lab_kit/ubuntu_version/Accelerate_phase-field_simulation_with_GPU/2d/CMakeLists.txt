cmake_minimum_required(VERSION 3.22.1 FATAL_ERROR)

# set the project name and version
project(MyProject VERSION 1.0)

find_package(CUDA REQUIRED)

cuda_add_executable(main-gpu.exe main-gpu.cu
	write_vtk_grid_values_2D.cu )

target_compile_options(main-gpu.exe PUBLIC -O2 -arch=sm_86 -Wall)

target_compile_features(main-gpu.exe PUBLIC cuda_std_11)

target_link_libraries(main-gpu.exe -lm)

# usage
# cmake -S . -B build/ -G"Unix Makefiles"
# cmake --build build/ --target main-gpu.exe
