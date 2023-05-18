cmake_minimum_required(VERSION 2.8)
project(MyCudaProject)

# For CUDA
find_package( CUDA REQUIRED )
set( CUDA_NVCC_FLAGS "-O2 -arch=native -lm" )
if(CUDA_FOUND)
	message("Cuda version: " ${CUDA_VERSION})
	CUDA_ADD_EXECUTABLE(main-shared_3d.exe
		main-shared_3d.cu
		write_vtk_grid_values_3D.cu
	)
	include_directories( ${CUDA_INCLUDE_DIRS} )
endif()

# Usage
# 1. cmake -S . -B build/ -G"Unix Makefiles"
# 2. cmake --build build/ --target main-shared_3d.exe
# 3. cd ./build
# 4. ./main-shared_3d.exe
# 5. (show time_XX.vtk on ParaView)

# Help
# cmake --help-module FindCUDA