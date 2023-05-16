cmake_minimum_required(VERSION 2.8)
project(MyCudaProject)

# For CUDA
find_package( CUDA REQUIRED )
set( CUDA_NVCC_FLAGS "-O2 -arch=native -lm" )
if(CUDA_FOUND)
	message("Cuda version: " ${CUDA_VERSION})
	CUDA_ADD_EXECUTABLE(main-gpu_2d.exe
		main-gpu_2d.cu
		write_vtk_grid_values_2D.cu
	)
	include_directories( ${CUDA_INCLUDE_DIRS} )
endif()

# usage
# cmake -S . -B build/ -G"Unix Makefiles"
# cmake --build build/ --target main-gpu_2d.exe

# help
# cmake --help-module FindCUDA