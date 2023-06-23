cmake_minimum_required(VERSION 2.8)
project(MyCudaProject)

# For CUDA
find_package( CUDA REQUIRED )
set( CUDA_NVCC_FLAGS "-O2 -arch=native -lm" )
if(CUDA_FOUND)
	message("Cuda version: " ${CUDA_VERSION})
	CUDA_ADD_EXECUTABLE(pfc_def_3d_cufft.exe
		pfc_def_3d_cufft.cu
		prepare_fft_3d.cu
		write_vtk_grid_values_3D.cu
	)
	include_directories( ${CUDA_INCLUDE_DIRS} )
	#link_directories(${CUDA_LIBRARY_DIRS})
	target_link_libraries(pfc_def_3d_cufft.exe -lcufft)
endif()

# Usage
# 1. cmake -S . -B build/ -G"Unix Makefiles"
# 2. cmake --build build/ --target pfc_def_3d_cufft.exe
# 3. cd ./build
# 4. cp ../bi_2r_3d.inp ./
# 5. ./pfc_def_3d_cufft.exe
# 6. (show time_XX.vtk on ParaView)

# Help
# cmake --help-module FindCUDA