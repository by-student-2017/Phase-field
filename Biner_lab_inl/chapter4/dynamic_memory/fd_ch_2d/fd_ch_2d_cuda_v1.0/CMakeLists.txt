cmake_minimum_required(VERSION 2.8)
project(MyCudaProject)

# For CUDA
find_package( CUDA REQUIRED )
set( CUDA_NVCC_FLAGS "-O2 -arch=native -lm" )
if(CUDA_FOUND)
	message("Cuda version: " ${CUDA_VERSION})
	CUDA_ADD_EXECUTABLE(fd_ch_2d_cuda.exe
		fd_ch_2d_cuda.cu
		micro_ch_pre_2d.cu
		calculate_energy_2d.cu
		write_vtk_grid_values_2D.cu
	)
	include_directories( ${CUDA_INCLUDE_DIRS} )
	#link_directories(${CUDA_LIBRARY_DIRS})
	target_link_libraries(fd_ch_2d_cuda.exe)
endif()

# Usage
# 1. cmake -S . -B build/ -G"Unix Makefiles"
# 2. cmake --build build/ --target fd_ch_2d_cuda.exe
# 3. cd ./build
# 4. ./fd_ch_2d_cuda.exe
# 5. (show time_XX.vtk on ParaView)

# Help
# cmake --help-module FindCUDA
