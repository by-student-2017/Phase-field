cmake_minimum_required(VERSION 2.8)
project(MyCudaProject)

# For CUDA
find_package( CUDA REQUIRED )
set( CUDA_cuda_FLAGS "-O2 -arch=native -lm" )
if(CUDA_FOUND)
	message("Cuda version: " ${CUDA_VERSION})
	CUDA_ADD_EXECUTABLE(fd_ca_2d_cuda.exe
		fd_ca_2d_cuda.cu
		init_grain_micro_2d.cu
		write_vtk_grid_values_2D.cu
	)
	include_directories( ${CUDA_INCLUDE_DIRS} )
endif()

# Usage
# 1. cmake -S . -B build/ -G"Unix Makefiles"
# 2. cmake --build build/ --target fd_ca_2d_cuda.exe
# 3. cd ./build
# 4. cp ../grain_25.inp ./
# 5. ./fd_ca_2d_cuda.exe
# 6. (show time_XX.vtk on ParaView)

# Help
# cmake --help-module FindCUDA