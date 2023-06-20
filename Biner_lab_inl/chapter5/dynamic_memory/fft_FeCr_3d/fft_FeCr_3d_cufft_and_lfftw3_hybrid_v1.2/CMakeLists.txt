cmake_minimum_required(VERSION 2.8)
project(MyCudaProject)

# For CUDA
find_package( CUDA REQUIRED )
set( CUDA_NVCC_FLAGS "-O2 -arch=native -lm" )
if(CUDA_FOUND)
	message("Cuda version: " ${CUDA_VERSION})
	CUDA_ADD_EXECUTABLE(fft_FeCr_3d_cufft.exe
		fft_FeCr_3d_cufft.cu
		dislo_strain_3d.cu
		FeCr_chem_poten_3d.cu
		green_tensor1_3D.cu
		green_tensor2_3D.cu
		prepare_fft_3d.cu
		solve_elasticity_3d.cu
		micro_ch_pre_3d.cu
		write_vtk_grid_values_3D.cu
	)
	include_directories( ${CUDA_INCLUDE_DIRS} )
	#link_directories(${CUDA_LIBRARY_DIRS})
	target_link_libraries(fft_FeCr_3d_cufft.exe -lcufft -lfftw3)
endif()

# Usage
# 1. cmake -S . -B build/ -G"Unix Makefiles"
# 2. cmake --build build/ --target fft_FeCr_3d_cufft.exe
# 3. cd ./build
# 4. ./fft_FeCr_3d_cufft.exe
# 5. (show time_XX.vtk on ParaView)

# Help
# cmake --help-module FindCUDA