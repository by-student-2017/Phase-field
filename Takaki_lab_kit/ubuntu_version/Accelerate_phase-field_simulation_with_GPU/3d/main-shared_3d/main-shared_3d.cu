/* Program : 2D Phase-Field Simulation for 
   Spinodal Decomposition in Fe-Cr Alloy by GPU (shared memory version) Computation.
   (e.g., Fe-Cr, Fe-Mo, Al-Zn, etc)
   
   Programmer : Akinori Yamanaka (original version)
   Place : Depertment of Mechanical and Control Engineering Tokyo Institute of Technology
   Date : 7th, July, 2010 
   
   Modified version: By Student
   Place: 2-1 Hirosawa, Wako, Saitama, 351-0198, Japan
   Date: 15th/May/2023
   Test: Ubuntu 22.04 LTS
   
   Compling: nvcc -O2 main-shared_2d.cu write_vtk_grid_values_2D.cu -o main-shared_2d.exe -arch=native -lm --std 'c++17'
   Run: ./main-shared_2d.exe
   ParaView: time_XX.vtk
*/

#include <stdio.h>  //printf()
#include <stdlib.h> //rand() and malloc()
#include <math.h>   //mod() and -lm
//----- ----- -----
#include <cuda.h>   //GPU
/* #include <cuda.h> or
  #include "cuda_runtime.h"
  #include "device_launch_parameters.h" */
//----- ----- -----

#define BSX  8        //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
#define BSY  8        //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
#define BSZ  4        //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
#define TIMES 1
//----- ----- -----
#define NX 256*TIMES //Number of grid points in the x-direction
#define NY 256*TIMES //Number of grid points in the y-direction
#define NZ  16*TIMES //Number of grid points in the z-direction

// Define subroutine "Kernel" for GPU (Device) calculation in detail
__global__ void Kernel
(
	float *f,
	float *fn,
	int    nx,
	int    ny,
	int    nz,
	float  rr,
	float  temp,
	float  L0,
	float  kapa_c,
	float  Da,
	float  Db,
	float  dt,
	float  dx,
	float  dy,
	float  dz
)
{
	int j, jx, jy, jz;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	float  fcc,
		   fce, fcw, fcs, fcn, fcu,  fcd,
		   //----- ----- -----XY
		   fcnw, fcne,
		   fcsw, fcse,
		   //----- ----- -----XZ
		   fcuw, fcue,
		   fcdw, fcde,
		   //----- ----- -----YZ
		   fcun, fcus,
		   fcdn, fcds,
		   //----- ----- -----
		   fcww, fcee, fcnn, fcss, fcuu, fcdd,
		   //----- ----- ----- ----- ----- -----
		   RT = rr*temp,
		   //----- ----- ----- ----- ----- -----
		   mu_chc,
		   mu_chw, mu_che, mu_chn, mu_chs, mu_chu, mu_chd,
		   //----- ----- -----
		   mu_suc,
		   mu_suw, mu_sue, mu_sun, mu_sus, mu_suu, mu_sud,
		   //----- ----- -----
		   mu_c,
		   mu_w, mu_e, mu_n, mu_s, mu_u, mu_d, 
		   //----- ----- -----
		   nab_mu,
		   dfmdx, dfmdy, dfmdz, 
		   //----- ----- -----
		   Dab = Db/Da,
		   mcc, dmc,
		   //----- ----- -----
		   dfdt ;
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// Declare variables used when copying data from global memory to shared memory
	int joff;
	int J0, J1, J2, J3;   // One inner edge
	int J4, J5, J6, J7;   // The most edge
	int J8, J9, J10, J11; // Corner
	//
	int J0z, J1z;   // One inner edge
	int J4z, J5z;   // The most edge
	int J8xz, J9xz, J10xz, J11xz; // XZ
	int J8yz, J9yz, J10yz, J11yz; // XZ
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// 16 kB before GF100 Core, 48 kB after GF100 Core
	const int thread_x = BSX;
	const int thread_y = BSY;
	const int thread_z = BSZ;
	// int=4B, float=4B, double=8B
	// float: 16x16 matrix = 1024 B = 1.024 kB
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* blockDim.x * blockDim.y = 16 * 16. In addition, 
	   add the necessary two adjacent difference grid points to
	   both sides of the x-axis and y-axis, respectively. */
	const int fs_thread_x = (2+thread_x+2);
	const int fs_thread_y = (2+thread_y+2);
	const int fs_thread_z = (2+thread_z+2);
	__shared__ float fs[fs_thread_x][fs_thread_y][fs_thread_z]; //fs is shared memory
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	/* Note ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	   Shared memory, fs = (2+blockDim.x+2) * (2+blockDim.y+2) * (2+blockDim.z+2) matrix
	     jx = threadIdx.x + 2; // +2 when counting from one end
	     jy = threadIdx.y + 2; // +2 when counting from one end
	     jy = threadIdx.z + 2; // +2 when counting from one end
	     ----- ----- ----- ---------- ----- -----
	     threadId+2 (on (2+blockDim+2) * (2+blockDim+2) * (2+blockDim+2) matrix): fine grid
	   ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	   Global memory, f = (gridDim.x*blockDim.x) * (gridDim.y*blockDim.y) * (gridDim.z*blockDim.z) matrix
	     x = (blockDim.x*blockIdx.x + threadIdx.x)
	     y = (blockDim.y*blockIdx.y + threadIdx.y)
	     z = (blockDim.z*blockIdx.z + threadIdx.z)
	     ----- ----- ----- ---------- ----- -----
	     blockId  (on gridDim * gridDim * gridDim  matrix): coarse grid
	       (blockId -(more detail)-> theadId)
	     threadId (on blockDim * blockDim * blockDim matrix): fine grid
	----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----  */
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	jx = threadIdx.x + 2; // +2 when counting from one end
	jy = threadIdx.y + 2; // +2 when counting from one end
	jz = threadIdx.z + 2; // +2 when counting from one end
	joff = nx*ny*(blockDim.z*blockIdx.z) + nx*(blockDim.y*blockIdx.y) + blockDim.x*blockIdx.x; // blockId matrix
	j = joff + nx*ny*threadIdx.z + nx*threadIdx.y + threadIdx.x; // f matrix
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	fcc = f[j]; // Global memory at current (jx,jy) grid point
	fs[jx][jy][jz] = fcc; // Global memory to shared memory at current (jx,jy) grid point
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// Calculating sleeve area in shared memory
	//----- ----- ----- down sleeve area ----- ----- ----- 
	if(blockIdx.z == 0) {J0z =   j +nx*ny*(nz-1),
						 J4z = J0z +nx*ny*(  -1);} // boundary condition at south edge
	else                {J0z =   j +nx*ny*(  -1), 
						 J4z = J0z +nx*ny*(  -1);} // non edge
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.z ==  0)          { fs[jx][jy][ 1] = f[J0z], fs[jx][jy][ 0] = f[J4z];}   // south sleeve area
	//
	//----- ----- ----- up sleeve area ----- ----- ----- 
	if(blockIdx.z == gridDim.z - 1) {J1z = j -nx*ny*(nz-1), 
						 J5z = J1z -nx*ny*(  -1);} // boundary condition at north edge
	else				{J1z =   j -nx*ny*(  -1), 
						 J5z = J1z -nx*ny*(  -1);} // non edge
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.z == (thread_z-1)){ fs[jx][jy][fs_thread_z-2] = f[J1z], fs[jx][jy][fs_thread_z-1] = f[J5z];}   // north sleeve area
	//
	//----- ----- ----- south sleeve area ----- ----- ----- 
	if(blockIdx.y == 0) {J0 =  j +nx*(ny-1),
						 J4 = J0 +nx*(  -1);} // boundary condition at south edge
	else                {J0 =  j +nx*(  -1),
						 J4 = J0 +nx*(  -1);} // non edge (J0=j+nx*(0-1),J4=j+nx(-1-1)=J0+nx*(   -1))
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.y ==  0)          { fs[jx][ 1][jz] = f[J0], fs[jx][ 0][jz] = f[J4];}   // south sleeve area
	//
	//----- ----- ----- north sleeve area ----- ----- ----- 
	if(blockIdx.y == gridDim.y - 1) {J1 = j-nx*(ny-1),
						 J5 = J1 -nx*(  -1);} // boundary condition at north edge
	else				{J1 =  j -nx*(  -1), 
						 J5 = J1 -nx*(  -1);} // non edge
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.y == (thread_y-1)){ fs[jx][fs_thread_y-2][jz] = f[J1], fs[jx][fs_thread_y-1][jz] = f[J5];}   // north sleeve area
	//
	//----- ----- ----- west sleeve area ----- ----- ----- 
	if(blockIdx.x == 0) {J2 =  j +(nx-1),
						 J6 = J2 +(  -1);} // boundary condition at west edge
	else				{J2 =  j +(  -1), 
						 J6 = J2 +(  -1);} // non edge
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.x ==  0)          { fs[ 1][jy][jz] = f[J2], fs[ 0][jy][jz] = f[J6];}   // west  sleeve area
	//
	//----- ----- ----- east sleeve area ----- ----- ----- 
	if(blockIdx.x == gridDim.x - 1) {J3 = j -(nx-1),
						 J7 = J3 -(  -1);} // boundary condition at east edge
	else				{J3 =  j -(  -1),
						 J7 = J3 -(  -1);} // non edge
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.x == (thread_x-1)){ fs[fs_thread_x-2][jy][jz] = f[J3], fs[fs_thread_x-1][jy][jz] = f[J7];}   // east  sleeve area
	//
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- XY = YX
	//----- ----- ----- east and north sleeve area
		 if(blockIdx.y == (gridDim.y - 1) && blockIdx.x == 0) { J8 = j-nx*(ny-1) +(nx-1) ;}
	else if(blockIdx.y == (gridDim.y - 1) && blockIdx.x  > 0) { J8 = j-nx*(ny-1) +(  -1) ;}
	else if(blockIdx.y  < (gridDim.y - 1) && blockIdx.x == 0) { J8 = j-nx*(  -1) +(nx-1) ;}
	else                                                      { J8 = j-nx*(  -1) +(  -1) ;}
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.x == 0            && threadIdx.y == (thread_y-1)) {fs[         1][thread_y+2][jz] = f[J8];}  // east and south sleeve area
	//
	//----- ----- ----- east and north sleeve area
		 if(blockIdx.y == (gridDim.y - 1) && blockIdx.x == (gridDim.x - 1)) { J9 = j-nx*(ny-1) -(nx-1);}
	else if(blockIdx.y == (gridDim.y - 1) && blockIdx.x  < (gridDim.x - 1)) { J9 = j-nx*(ny-1) -(  -1);}
	else if(blockIdx.y  < (gridDim.y - 1) && blockIdx.x == (gridDim.x - 1)) { J9 = j-nx*(  -1) -(nx-1);}
	else                                                                    { J9 = j-nx*(  -1) -(  -1);} // non edge
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.x == (thread_x-1) && threadIdx.y == (thread_y-1)) {fs[thread_x+2][thread_y+2][jz] = f[J9];}  // east and north sleeve area
	//
	//----- ----- ----- west and south sleeve area
		 if(blockIdx.y == 0 && blockIdx.x  > 0) { J10 = j+nx*(ny-1) +(  -1);}
	else if(blockIdx.y  > 0 && blockIdx.x == 0) { J10 = j+nx*(  -1) +(nx-1);}
	else if(blockIdx.y == 0 && blockIdx.x == 0) { J10 = j+nx*(ny-1) +(nx-1);}
	else                                        { J10 = j+nx*(  -1) +(  -1);}
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.x == 0            && threadIdx.y ==  0          ) {fs[         1][         1][jz] = f[J10];} // west and south sleeve area
	//
	//----- ----- ----- west and south sleeve area
		 if(blockIdx.y == 0 && blockIdx.x == (gridDim.x -1)) { J11 = j+nx*(ny-1) -(nx-1);}
	else if(blockIdx.y == 0 && blockIdx.x  < (gridDim.x -1)) { J11 = j+nx*(ny-1) -(  -1);}
	else if(blockIdx.y  > 0 && blockIdx.x == (gridDim.x -1)) { J11 = j+nx*(  -1) -(nx-1);}
	else                                                     { J11 = j+nx*(  -1) -(  -1);}
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.x == (thread_x-1) && threadIdx.y ==  0          ) {fs[thread_x+2][         1][jz] = f[J11];} // west and north sleeve area
	//
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- XZ = ZX
	//----- ----- ----- up and east sleeve area
		 if(blockIdx.z == (gridDim.z - 1) && blockIdx.x == 0) { J8xz = j-nx*ny*(nz-1) +(nx-1) ;}
	else if(blockIdx.z == (gridDim.z - 1) && blockIdx.x  > 0) { J8xz = j-nx*ny*(nz-1) +(  -1) ;}
	else if(blockIdx.z  < (gridDim.z - 1) && blockIdx.x == 0) { J8xz = j-nx*ny*(  -1) +(nx-1) ;}
	else                                                      { J8xz = j-nx*ny*(  -1) +(  -1) ;}
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.x == 0            && threadIdx.z == (thread_z-1)) {fs[         1][jy][thread_z+2] = f[J8xz];}  // east and south sleeve area
	//
	//----- ----- ----- up and east sleeve area
		 if(blockIdx.z == (gridDim.z - 1) && blockIdx.x == (gridDim.x - 1)) { J9xz = j-nx*ny*(nz-1) -(nx-1);}
	else if(blockIdx.z == (gridDim.z - 1) && blockIdx.x  < (gridDim.x - 1)) { J9xz = j-nx*ny*(nz-1) -(  -1);}
	else if(blockIdx.z  < (gridDim.z - 1) && blockIdx.x == (gridDim.x - 1)) { J9xz = j-nx*ny*(  -1) -(nx-1);}
	else                                                                    { J9xz = j-nx*ny*(  -1) -(  -1);} // non edge
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.x == (thread_x-1) && threadIdx.z == (thread_z-1)) {fs[thread_x+2][jy][thread_z+2] = f[J9xz];}  // east and north sleeve area
	//
	//----- ----- ----- down and west sleeve area
		 if(blockIdx.z == 0 && blockIdx.x  > 0) { J10xz = j+nx*ny*(nz-1) +(  -1);}
	else if(blockIdx.z  > 0 && blockIdx.x == 0) { J10xz = j+nx*ny*(  -1) +(nx-1);}
	else if(blockIdx.z == 0 && blockIdx.x == 0) { J10xz = j+nx*ny*(nz-1) +(nx-1);}
	else                                        { J10xz = j+nx*ny*(  -1) +(  -1);}
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.x == 0            && threadIdx.z ==  0          ) {fs[         1][jy][         1] = f[J10xz];} // west and south sleeve area
	//
	//----- ----- ----- down and west sleeve area
		 if(blockIdx.z == 0 && blockIdx.x == (gridDim.x -1)) { J11xz = j+nx*ny*(nz-1) -(nx-1);}
	else if(blockIdx.z == 0 && blockIdx.x  < (gridDim.x -1)) { J11xz = j+nx*ny*(nz-1) -(  -1);}
	else if(blockIdx.z  > 0 && blockIdx.x == (gridDim.x -1)) { J11xz = j+nx*ny*(  -1) -(nx-1);}
	else                                                     { J11xz = j+nx*ny*(  -1) -(  -1);}
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.x == (thread_x-1) && threadIdx.z ==  0          ) {fs[thread_x+2][jy][         1] = f[J11xz];} // west and north sleeve area
	//
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- YZ = ZY
	//----- ----- ----- up and south sleeve area
		 if(blockIdx.z == (gridDim.z - 1) && blockIdx.y == 0) { J8yz = j-nx*ny*(nz-1) +nx*(ny-1) ;}
	else if(blockIdx.z == (gridDim.z - 1) && blockIdx.y  > 0) { J8yz = j-nx*ny*(nz-1) +nx*(  -1) ;}
	else if(blockIdx.z  < (gridDim.z - 1) && blockIdx.y == 0) { J8yz = j-nx*ny*(  -1) +nx*(ny-1) ;}
	else                                                      { J8yz = j-nx*ny*(  -1) +nx*(  -1) ;}
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.y == 0            && threadIdx.z == (thread_z-1)) {fs[jx][         1][thread_z+2] = f[J8yz];}  // east and south sleeve area
	//
	//----- ----- ----- up and north sleeve area
		 if(blockIdx.z == (gridDim.z - 1) && blockIdx.y == (gridDim.y - 1)) { J9yz = j-nx*ny*(nz-1) -nx*(ny-1);}
	else if(blockIdx.z == (gridDim.z - 1) && blockIdx.y  < (gridDim.y - 1)) { J9yz = j-nx*ny*(nz-1) -nx*(  -1);}
	else if(blockIdx.z  < (gridDim.z - 1) && blockIdx.y == (gridDim.y - 1)) { J9yz = j-nx*ny*(  -1) -nx*(ny-1);}
	else                                                                    { J9yz = j-nx*ny*(  -1) -nx*(  -1);} // non edge
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.y == (thread_y-1) && threadIdx.z == (thread_z-1)) {fs[jx][thread_y+2][thread_z+2] = f[J9yz];}  // east and north sleeve area
	//
	//----- ----- ----- down and south sleeve area
		 if(blockIdx.z == 0 && blockIdx.y  > 0) { J10yz = j+nx*ny*(nz-1) +nx*(  -1);}
	else if(blockIdx.z  > 0 && blockIdx.y == 0) { J10yz = j+nx*ny*(  -1) +nx*(ny-1);}
	else if(blockIdx.z == 0 && blockIdx.y == 0) { J10yz = j+nx*ny*(nz-1) +nx*(ny-1);}
	else                                        { J10yz = j+nx*ny*(  -1) +nx*(  -1);}
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.y == 0            && threadIdx.z ==  0          ) {fs[jx][         1][         1] = f[J10yz];} // west and south sleeve area
	//
	//----- ----- ----- down and north sleeve area
		 if(blockIdx.z == 0 && blockIdx.y == (gridDim.y -1)) { J11yz = j+nx*ny*(nz-1) -nx*(ny-1);}
	else if(blockIdx.z == 0 && blockIdx.y  < (gridDim.y -1)) { J11yz = j+nx*ny*(nz-1) -nx*(  -1);}
	else if(blockIdx.z  > 0 && blockIdx.y == (gridDim.y -1)) { J11yz = j+nx*ny*(  -1) -nx*(ny-1);}
	else                                                     { J11yz = j+nx*ny*(  -1) -nx*(  -1);}
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.y == (thread_y-1) && threadIdx.z ==  0          ) {fs[jx][thread_y+2][         1] = f[J11yz];} // west and north sleeve area
	//
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// Wait until all data is secured
	__syncthreads();
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* Consider up to the second nearest neighbor. Therefore, 
	   25 difference grid points are used. (#1 to #25) */
	/* Array data is put in the argument so that 
	   the code on the CPU can be used as it is. */
	fcc  = fs[jx  ][jy  ][jz  ]; // #1 (center: fc)
	fcw  = fs[jx-1][jy  ][jz  ]; // #2 (center: fc)
	fce  = fs[jx+1][jy  ][jz  ]; // #3 (center: fc)
	fcn  = fs[jx  ][jy+1][jz  ]; // #4 (center: fc)
	fcs  = fs[jx  ][jy-1][jz  ]; // #5 (center: fc)
	fcu  = fs[jx  ][jy  ][jz+1]; // #6 (center: fc)
	fcd  = fs[jx  ][jy  ][jz-1]; // #7 (center: fc)
	//
	// XY
	fcnw = fs[jx-1][jy+1][jz  ]; // #8  (center: fcn)
	fcne = fs[jx+1][jy+1][jz  ]; // #9  (center: fcn)
	fcsw = fs[jx-1][jy-1][jz  ]; // #10 (center: fcs)
	fcse = fs[jx+1][jy-1][jz  ]; // #11 (center: fcs)
	//
	// XZ
	fcuw = fs[jx-1][jy  ][jz+1]; // #12 (center: fcu)
	fcue = fs[jx+1][jy  ][jz+1]; // #13 (center: fcu)
	fcdw = fs[jx-1][jy  ][jz-1]; // #14 (center: fcd)
	fcde = fs[jx+1][jy  ][jz-1]; // #15 (center: fcd)
	//
	// YZ
	fcus = fs[jx  ][jy-1][jz+1]; // #16 (center: fcu)
	fcun = fs[jx  ][jy+1][jz+1]; // #17 (center: fcu)
	fcds = fs[jx  ][jy-1][jz-1]; // #18 (center: fcd)
	fcdn = fs[jx  ][jy+1][jz-1]; // #19 (center: fcd)
	//
	fcww = fs[jx-2][jy  ][jz  ]; // #20 (center: fcw)
	fcee = fs[jx+2][jy  ][jz  ]; // #21 (center: fce)
	fcnn = fs[jx  ][jy+2][jz  ]; // #22 (center: fcn)
	fcss = fs[jx  ][jy-2][jz  ]; // #23 (center: fcs)
	fcuu = fs[jx  ][jy  ][jz+2]; // #24 (center: fcu)
	fcdd = fs[jx  ][jy  ][jz-2]; // #25 (center: fcd)
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// term1 = Atomic_interaction*(1-2*f) + RT*{log(f) - log(1-f)}
	mu_chc = L0*(1.0-2.0*fcc) + RT*( log(fcc) - log(1.0-fcc) ); //center: fcc
	//
	mu_chw = L0*(1.0-2.0*fcw) + RT*( log(fcw) - log(1.0-fcw) ); //center: fcw
	mu_che = L0*(1.0-2.0*fce) + RT*( log(fce) - log(1.0-fce) ); //center: fce
	//
	mu_chn = L0*(1.0-2.0*fcn) + RT*( log(fcn) - log(1.0-fcn) ); //center: fcn
	mu_chs = L0*(1.0-2.0*fcs) + RT*( log(fcs) - log(1.0-fcs) ); //center: fcs
	//
	mu_chu = L0*(1.0-2.0*fcu) + RT*( log(fcu) - log(1.0-fcu) ); //center: fcu
	mu_chd = L0*(1.0-2.0*fcd) + RT*( log(fcd) - log(1.0-fcd) ); //center: fcd
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// term2 = -gradient_energy_coefficient * Laplacian(f)
	mu_suc = -kapa_c*( (fce  + fcw  -2.0*fcc)/(dx*dx)
					 + (fcn  + fcs  -2.0*fcc)/(dy*dy)
					 + (fcu  + fcd  -2.0*fcc)/(dz*dz) ); //center: fcc
	//
	mu_suw = -kapa_c*( (fcc  + fcww -2.0*fcw)/(dx*dx)
					 + (fcnw + fcsw -2.0*fcw)/(dy*dy)
					 + (fcuw + fcdw -2.0*fcw)/(dz*dz) ); //fcc=fcwe, fcnw=fcwn, fcsw=fcws, //center: fcw
	mu_sue = -kapa_c*( (fcee + fcc  -2.0*fce)/(dx*dx)
					 + (fcne + fcse -2.0*fce)/(dy*dy)
					 + (fcue + fcde -2.0*fce)/(dz*dz) ); //fcc=fcew, fcne=fcen, fcse=fces, //center: fce
	//
	mu_sun = -kapa_c*( (fcne + fcnw -2.0*fcn)/(dx*dx)
					 + (fcnn + fcc  -2.0*fcn)/(dy*dy)
					 + (fcun + fcdn -2.0*fcn)/(dz*dz) ); //fcc=fcns, fcun=fcnu, fcdn=fcnd, //center: fcn
	mu_sus = -kapa_c*( (fcse + fcsw -2.0*fcs)/(dx*dx)
					 + (fcc  + fcss -2.0*fcs)/(dy*dy)
					 + (fcus + fcds -2.0*fcs)/(dz*dz) ); //fcc=fcsn, fcsu=fcus, fcds=fcsd, //center: fcs
	//
	mu_suu = -kapa_c*( (fcue + fcuw -2.0*fcu)/(dx*dx)
					 + (fcun + fcus -2.0*fcu)/(dy*dy)
					 + (fcuu + fcc  -2.0*fcu)/(dz*dz) ); //fcc=fcud, //center: fcu
	mu_sud = -kapa_c*( (fcde + fcdw -2.0*fcd)/(dx*dx)
					 + (fcdn + fcds -2.0*fcd)/(dy*dy)
					 + (fcc  + fcdd -2.0*fcd)/(dz*dz) ); //fcc=fcdu, //center: fcd
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// mu = dG/df = term1 + term2
	mu_c = mu_chc + mu_suc; // at current (jx,jy,jz) grid point, //center: fcc
	//
	mu_w = mu_chw + mu_suw; // at (jx-1,jy,jz) grid point, //center: fcw
	mu_e = mu_che + mu_sue; // at (jx+1,jy,jz) grid point, //center: fce
	//
	mu_n = mu_chn + mu_sun; // at (jx,jy+1,jz) grid point, //center: fcn
	mu_s = mu_chs + mu_sus; // at (jx,jy-1,jz) grid point, //center: fcs
	//
	mu_u = mu_chu + mu_suu; // at (jx,jy,jz+1) grid point, //center: fcu
	mu_d = mu_chd + mu_sud; // at (jx,jy,jz-1) grid point, //center: fcd
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// Laplacian(mu) = d^2(mu)/dx^2 + d^2(mu)/dy^2
	nab_mu = (mu_e + mu_w -2.0*mu_c)/(dx*dx)  // d^2(mu)/dx^2
		   + (mu_n + mu_s -2.0*mu_c)/(dy*dy)  // d^2(mu)/dy^2
		   + (mu_u + mu_d -2.0*mu_c)/(dz*dz); // d^2(mu)/dz^2
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// (df/dx) * d(mu)/dx, (x is related with e and w), (the center is fc.)
	dfmdx = ( (fce - fcw)/(2.0*dx) * (mu_e - mu_w)/(2.0*dx) );
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// (df/dy) * d(mu)/dy, (y is related with n and s), (the center is fc.)
	dfmdy = ( (fcn - fcs)/(2.0*dy) * (mu_n - mu_s)/(2.0*dy) );
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// (df/dz) * d(mu)/dz, (z is related with u and d), (the center is fc.)
	dfmdz = ( (fcu - fcd)/(2.0*dz) * (mu_u - mu_d)/(2.0*dz) );
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// Mobility, M = { (D_A/RT)*c + (D_B/RT)*(1-c) }*c*(1-c)
	//             = (D_a/RT)*{f + (D_B/D_A)*(1-f)}*f*(1-f)
	mcc = (Da/RT)*(fcc+Dab*(1.0-fcc))*fcc*(1.0-fcc); 
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// dM/df
	dmc = (Da/RT)*((1.0-Dab)*fcc*(1.0-fcc) + (fcc+Dab*(1.0-fcc))*(1.0-2.0*fcc)); 
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// df/dt = M*Laplacian(mu) + (dM/df)*( (df/dx) * d(mu)/dx + (df/dy) * d(mu)/dy + (df/dz) * d(mu)/dz )
	dfdt = mcc*nab_mu + dmc*(dfmdx + dfmdy + dfmdz); 
	fn[j] = f[j] + dfdt*dt;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	// If there are small variations, set the max and min values to the limits
	//if(fn[j]>0.999999){ fn[j]=0.999999; }
	//if(fn[j]<1.0e-6)  { fn[j]=1.0e-6; }
}

void update(float **f, float **fn)
{
	float *tmp = *f;
	*f  = *fn;
	*fn = tmp;
}

void write_vtk_grid_values_3D(int Nx, int Ny, int Nz, float dx, float dy, float dz, int istep, float *data1);

int main(int argc, char** argv)
{
	float *f_d, *fn_d; // name of dynamic memory for GPU, CUDA, device
	float *F_h;        // name of dynamic memory for CPU
	int nx = NX, ny = NY, nz = NZ;
	int times = TIMES;
	
	int nstep=10000;    //Number of time integration steps
	int nprint=1000;    //Output frequency to write the results to file
	//----- ----- ----- -----
	float Lx = 3.0e-07*times, // Simulation length in x-direction [m]
		  Ly = 3.0e-07*times, // Simulation length in y-direction [m]
		  Lz = 3.0e-07*times*nz/sqrt(nx*ny), // Simulation length in z-direction [m]
		  //----- ----- ----- -----
		  dx = Lx/(float)nx, // Grid spacing between two grid pints in x-direction [m]
		  dy = Ly/(float)ny, // Grid spacing between two grid pints in y-direction [m]
		  dz = Lz/(float)nz, // Grid spacing between two grid pints in z-direction [m]
		  //----- ----- ----- -----
		  c_0 = 0.4,    // Initial concentration of Cr (atomic fraction)
		  //----- ----- ----- -----
		  rr  = 8.314,  // Gas constant [J/(mol*K)]
		  temp = 623.0, // Temperature [K]
		  RT = rr*temp, // RT [J/mol]
		  //----- ----- ----- -----
		  L0 = 21020.8-9.31889*temp, // Atomic interaction [J/mol] (Magnetism is omitted)
		  kapa_c = 1.2e-14,  // The value of gradient energy coefficients [J*m^2/mol]
		  //----- ----- ----- -----
		  Da = 1.0e-04*exp(-294000.0/RT), // Self-diffusion coefficient [m^2/s] (Fe)
		  Db = 2.0e-05*exp(-308000.0/RT), // Self-diffusion coefficient [m^2/s] (Cr)
		  //----- ----- ----- -----
		  dt = (dx*dx/Da)*(4.0/6.0)*0.01; // Time increment for the numerical integration [s] (2D -> 3D; n=4 -> n=6)
	
	//----- ----- ----- -----start:(This part is not really needed.)----- ----- ----- ----
	int nDevices;
	cudaGetDeviceCount(&nDevices);
	for (int i = 0; i < nDevices; i++){
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, i);
		printf("--------------------------------------------------\n");
		printf("Device Number: %d\n", i);
		printf("  Device name: %s\n", prop.name);
		printf("  Memory Clock Rate (KHz): %d\n",prop.memoryClockRate);
		printf("  Memory Bus Width (bits): %d\n",prop.memoryBusWidth);
		printf("  Peak Memory Bandwidth (GB/s): %f\n",2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
		printf("--------------------------------------------------\n");
	}
	//----- ----- ----- -----end:(This part is not really needed.)----- ----- ----- ----
	
	f_d  = (float *)malloc(nx*ny*nz*sizeof(float)); //GPU, CUDA, device
	fn_d = (float *)malloc(nx*ny*nz*sizeof(float)); //GPU, CUDA, device
	
	cudaMalloc((void**)&f_d ,nx*ny*nz*sizeof(float)); // define dynamic memory for GPU (device)
	cudaMalloc((void**)&fn_d,nx*ny*nz*sizeof(float)); // define dynamic memory for GPU (device)
	
	F_h  = (float *)malloc(nx*ny*nz*sizeof(float));   // define dynamic memory for CPU (host)
	
	// Initialize the concentration filed f with random modulation
	for(int jz=0; jz<nz ; jz++){
		for(int jy=0; jy<ny ; jy++){
			for(int jx=0; jx<nx ; jx++){
				int j = (jz*ny + jy)*nx + jx; //j = nx*ny*jz + nx*jy + jx;
				float r = (float)rand()/(float)(RAND_MAX)-0.5;
				F_h[j] = c_0 + 0.01*r; //Initial fluctuation is assumed to be 1%.
			}
		}
	}//on CPU calculation
	
	//copy F_h(cpu,host) to f_d(cuda,device)
	cudaMemcpy(f_d,F_h,nx*ny*nz*sizeof(float),cudaMemcpyHostToDevice);
	
	int bsx=BSX, bsy=BSY, bsz=BSZ;     //Number of threads
	dim3 blocks(nx/bsx,ny/bsy,nz/bsz); //nx*ny*nz = blocks * threads
	dim3 threads(bsx,bsy,bsz);         //bsx*bsy*bsz <= 1024
	
	//----- ----- ----- -----start:(This part is not really needed.)----- ----- ----- ----
	//Set recording time
	cudaEvent_t start, stop;
	
	//Initialization
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	//Start recording time
	cudaEventRecord(start);
	//----- ----- ----- -----end:(This part is not really needed.)----- ----- ----- ----
	
	float tc = 0.0; //average concentration
	for(int istep=0; istep<=nstep ; istep++){
		//calculate subroutine "Kernel" on GPU
		Kernel<<<blocks, threads>>>(f_d,fn_d,nx,ny,nz,rr,temp,L0,kapa_c,Da,Db,dt,dx,dy,dz);
		cudaDeviceSynchronize(); //<- new version | old version -> cudaThreadSynchronize();
		
		// replace f_d with new f_d (=fn_d)
		update(&f_d,&fn_d);
		//
		if(istep%(nprint/10) == 0){
			//copy f_d(cuda,device) to F_h(cpu,host)
			cudaMemcpy(F_h,f_d,nx*ny*nz*sizeof(float),cudaMemcpyDeviceToHost);
			
			//check average concentration
			tc = 0.0; //average concentration
			for(int jz=0; jz<nz ; jz++){
				for(int jy=0; jy<ny ; jy++){
					for(int jx=0; jx<nx ; jx++){
						int j = (jz*ny + jy)*nx + jx; //j = nx*ny*jz + nx*jy + jx;
						tc = tc + F_h[j];
					}
				}
			}
			tc = tc/(nx*ny*nz);
			
			//correct concentration
			for(int jz=0; jz<nz ; jz++){
				for(int jy=0; jy<ny ; jy++){
					for(int jx=0; jx<nx ; jx++){
						int j = (jz*ny + jy)*nx + jx; //j = nx*ny*jz + nx*jy + jx;
						F_h[j] = F_h[j] * (c_0/tc);
					}
				}
			}
			
			//copy F_h(cpu,host) to f_d(cuda,device)
			cudaMemcpy(f_d,F_h,nx*ny*nz*sizeof(float),cudaMemcpyHostToDevice);
		}
		//
		if(istep%nprint == 0){
			//copy f_d(cuda,device) to F_h(cpu,host)
			//cudaMemcpy(F_h,f_d,nx*ny*nz*sizeof(float),cudaMemcpyDeviceToHost);
			
			//output vtk format
			write_vtk_grid_values_3D(nx,ny,nz,dx,dy,dz,istep,F_h);
			
			//show current step
			fprintf(stderr,"istep = %5d, average constration = %f, total annealing time= %e [s] at %f [K] \n",istep,tc,(dt*istep),temp);
		}
	}
	
	//----- ----- ----- -----start:(This part is not really needed.)----- ----- ----- ----
	//Stop recording time
	cudaEventRecord(stop);
	
	//Wait all event
	cudaEventSynchronize(stop);
	
	//calculate time. time is [ms] unit.
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);
	
	//Show computing time
	printf("Calculation Time = %9.3f [sec] \n",milliseconds*1.0e-03);
	
	//End processing
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	//----- ----- ----- -----end:(This part is not really needed.)----- ----- ----- ----
	
	cudaFree(f_d);
	cudaFree(fn_d);
	
	free(F_h);
	
	return 0;
}