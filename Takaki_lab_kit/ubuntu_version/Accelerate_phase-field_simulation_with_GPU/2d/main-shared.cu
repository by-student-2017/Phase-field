/* Program : 2D Phase-Field Simulation for 
   Spinodal Decomposition in Fe-Cr Alloy by GPU (shared memory version) Computation.
   
   Programmer : Akinori Yamanaka (original version)
   Place : Depertment of Mechanical and Control Engineering Tokyo Institute of Technology
   Date : 7th, July, 2010 
   
   Modified version: By Student
   Place: 2-1 Hirosawa, Wako, Saitama, 351-0198, Japan
   Date: 15th/May/2023
   Test: Ubuntu 22.04 LTS
   
   Compling: nvcc -O2 main-shared.cu write_vtk_grid_values_2D.cu -o main-shared.exe -arch=native -lm --std 'c++17'
   Run: ./main-shared.exe
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

#define BS 16        //Number of threads, 2^n =< 16 in this code (shared memory = (2+BS+2)^2)
#define TIMES 2
//----- ----- -----
#define NX 256*TIMES //Number of grid points in the x-direction
#define NY 256*TIMES //Number of grid points in the y-direction

// Define subroutine "Kernel" for GPU (Device) calculation in detail
__global__ void Kernel
(
	float *f,
	float *fn,
	int    nx,
	int    ny,
	float  rr,
	float  temp,
	float  L0,
	float  kapa_c,
	float  Da,
	float  Db,
	float  dt,
	float  dx,
	float  dy
)
{
	int j, jx, jy;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	float  fcc,
		   fce, fcw, fcs, fcn,
		   //----- ----- -----
		   fcnw, fcne,
		   fcsw, fcse,
		   //----- ----- -----
		   fcww, fcee, fcnn, fcss,
		   //----- ----- ----- ----- ----- -----
		   RT = rr*temp,
		   //----- ----- ----- ----- ----- -----
		   mu_chc,
		   mu_chw, mu_che, mu_chn, mu_chs,
		   //----- ----- -----
		   mu_suc,
		   mu_suw, mu_sue, mu_sun, mu_sus, 
		   //----- ----- -----
		   mu_c,
		   mu_w, mu_e, mu_n, mu_s, 
		   //----- ----- -----
		   nab_mu,
		   dfmdx, dfmdy,
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
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	const int nthreads = BS; // 16 kB before GF100 Core, 48 kB after GF100 Core
	const int thread_x = nthreads;
	const int thread_y = nthreads;
	// int=4B, float=4B, double=8B
	// float: 16x16 matrix = 1024 B = 1.024 kB
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* blockDim.x * blockDim.y = 16 * 16. In addition, 
	   add the necessary two adjacent difference grid points to
	   both sides of the x-axis and y-axis, respectively. */
	const int fs_thread_x = (2+thread_x+2);
	const int fs_thread_y = (2+thread_y+2);
	__shared__ float fs[fs_thread_x][fs_thread_y]; //fs is shared memory
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	/* Note ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	   Shared memory, fs = (2+blockDim.x+2) * (2+blockDim.y+2) matrix
	     jx = threadIdx.x + 2; // +2 when counting from one end
	     jy = threadIdx.y + 2; // +2 when counting from one end
	     ----- ----- ----- ---------- ----- -----
	     threadId+2 (on (2+blockDim+2) * (2+blockDim+2) matrix): fine grid
	   ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	   Global memory, f = (gridDim.x*blockDim.x) * (gridDim.y*blockDim.y) matrix
	     x = (blockDim.x*blockIdx.x + threadIdx.x)
	     y = (blockDim.y*blockIdx.y + threadIdx.y)
	     ----- ----- ----- ---------- ----- -----
	     blockId  (on gridDim  * gridDim  matrix): coarse grid
	       (blockId -(more detail)-> theadId)
	     threadId (on blockDim * blockDim matrix): fine grid
	----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----  */
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	jx = threadIdx.x + 2; // +2 when counting from one end
	jy = threadIdx.y + 2; // +2 when counting from one end
	joff = nx*(blockDim.y*blockIdx.y) + blockDim.x*blockIdx.x; // blockId matrix
	j = joff + nx*threadIdx.y + threadIdx.x; // f matrix
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	fcc = f[j]; // Global memory at current (jx,jy) grid point
	fs[jx][jy] = fcc; // Global memory to shared memory at current (jx,jy) grid point
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// Calculating sleeve area in shared memory
	//----- ----- ----- south sleeve area ----- ----- ----- 
	if(blockIdx.y == 0) {J0 = nx*(ny-1) + blockDim.x*blockIdx.x + threadIdx.x,
						 J4 = J0 - nx;} // boundary condition at south edge
	else                {J0 =  j - nx, 
						 J4 = J0 - nx;} // non edge
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.y ==  0)          { fs[jx][ 1] = f[J0], fs[jx][ 0] = f[J4];}   // south sleeve area
	//
	//----- ----- ----- north sleeve area ----- ----- ----- 
	if(blockIdx.y == gridDim.y - 1) {J1 = blockDim.x*blockIdx.x + threadIdx.x, 
						 J5 = J1 + nx;} // boundary condition at north edge
	else				{J1 =  j + nx, 
						 J5 = J1 + nx;} // non edge
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.y == (thread_y-1)){ fs[jx][fs_thread_y-2] = f[J1], fs[jx][fs_thread_y-1] = f[J5];}   // north sleeve area
	//
	//----- ----- ----- west sleeve area ----- ----- ----- 
	if(blockIdx.x == 0) {J2 =  j + (nx-1),
						 J6 = J2 - 1;} // boundary condition at west edge
	else				{J2 =  j - 1, 
						 J6 = J2 - 1;} // non edge
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.x ==  0)          { fs[ 1][jy] = f[J2], fs[ 0][jy] = f[J6];}   // west  sleeve area
	//
	//----- ----- ----- east sleeve area ----- ----- ----- 
	if(blockIdx.x == gridDim.x - 1) {J3 = j - (nx-1),
						 J7 = J3 + 1;} // boundary condition at east edge
	else				{J3 =  j + 1,
						 J7 = J3 + 1;} // non edge
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.x == (thread_x-1)){ fs[fs_thread_x-2][jy] = f[J3], fs[fs_thread_x-1][jy] = f[J7];}   // east  sleeve area
	//
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	//----- ----- ----- east and north sleeve area
		 if(blockIdx.x == 0 && blockIdx.y == gridDim.y - 1) { J8 = blockDim.x*gridDim.x - 1 ;} // edge(west and north)
	else if(blockIdx.x  > 0 && blockIdx.y == gridDim.y - 1) { J8 = J1 - 1 ;}            // edge(north)
	else if(blockIdx.x == 0 && blockIdx.y  < gridDim.y - 1) { J8 = j + nx + (nx-1) ;}   // edge(west)
	else                                                    { J8 = j +  0 + (nx-1) ;}   // non edge
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.x == 0            && threadIdx.y == (thread_y-1)) {fs[         1][thread_y+2] = f[J8];}  // east and south sleeve area
	//
	//----- ----- ----- east and north sleeve area
		 if(blockIdx.x == gridDim.x - 1 && blockIdx.y == gridDim.y - 1) { J9 = 0 ;}          // edge(east and north)
	else if(blockIdx.x  < gridDim.x - 1 && blockIdx.y == gridDim.y - 1) { J9 = J1 + 1 ;}     // edge(north)
	else if(blockIdx.x == gridDim.x - 1 && blockIdx.y  < gridDim.y - 1) { J9 = j +  0 + 1 ;}     // edge(east)
	else                                                                { J9 = j + nx + 1 ;} // non edge
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.x == (thread_x-1) && threadIdx.y == (thread_y-1)) {fs[thread_x+2][thread_y+2] = f[J9];}  // east and north sleeve area
	//
	//----- ----- ----- west and south sleeve area
		 if(blockIdx.x  > 0 && blockIdx.y == 0) { J10 = J0 - 1 ;}                       // edge(south)
	else if(blockIdx.x == 0 && blockIdx.y  > 0) { J10 = j - 1  ;}                       // edge(west)
	else if(blockIdx.x == 0 && blockIdx.y == 0) { J10 = ny*blockDim.x*blockDim.x - 1 ;} // edge(west and south)
	else                                        { J10 = j - nx - 1 ;}                   // non edge
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.x == 0            && threadIdx.y ==  0          ) {fs[         1][         1] = f[J10];} // west and south sleeve area
	//
	//----- ----- ----- west and south sleeve area
		 if(blockIdx.x == gridDim.x -1 && blockIdx.y == 0) { J11 = (ny-1)*blockDim.x*blockDim.x;} // edge(east and south)
	else if(blockIdx.x  < gridDim.x -1 && blockIdx.y == 0) { J11 = J0 + 1  ;}         // edge(south)
	else if(blockIdx.x == gridDim.x -1 && blockIdx.y  > 0) { J11 = j - nx - (nx-1) ;} // edge(east)
	else                                                   { J11 = j -  0 - (nx-1) ;} // non edge
	//----- ----- copy Global memory to Shared memory {one inside, edge}
	if(threadIdx.x == (thread_x-1) && threadIdx.y ==  0          ) {fs[thread_x+2][         1] = f[J11];} // west and north sleeve area
	//
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/*
	  copy Global memory to Shared memory {one inside, edge}
	  if(threadIdx.y ==  0)          { fs[jx][ 1] = f[J0], fs[jx][ 0] = f[J4];}   // south sleeve area
	  if(threadIdx.x ==  0)          { fs[ 1][jy] = f[J2], fs[ 0][jy] = f[J6];}   // west  sleeve area
	  if(threadIdx.x == (thread_x-1)){ fs[fs_thread_x-2][jy] = f[J3], fs[fs_thread_x-1][jy] = f[J7];}   // east  sleeve area
	  if(threadIdx.y == (thread_y-1)){ fs[jx][fs_thread_y-2] = f[J1], fs[jx][fs_thread_y-1] = f[J5];}   // north sleeve area
	//----- ----- ----- {one inside}
	  if(threadIdx.x == 0            && threadIdx.y == (thread_y-1)) {fs[         1][thread_y+2] = f[J8];}  // east and south sleeve area
	  if(threadIdx.x == (thread_x-1) && threadIdx.y == (thread_y-1)) {fs[thread_x+2][thread_y+2] = f[J9];}  // east and north sleeve area
	  if(threadIdx.x == 0            && threadIdx.y ==  0          ) {fs[         1][         1] = f[J10];} // west and south sleeve area
	  if(threadIdx.x == (thread_x-1) && threadIdx.y ==  0          ) {fs[thread_x+2][         1] = f[J11];} // west and north sleeve area
	*/
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// Wait until all data is secured
	__syncthreads();
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* Consider up to the second nearest neighbor. Therefore, 
	   13 difference grid points are used. (#1 to #13) */
	/* Array data is put in the argument so that 
	   the code on the CPU can be used as it is. */
	fcc  = fs[jx  ][jy  ]; // #1 (center: fc)
	fcw  = fs[jx-1][jy  ]; // #2 (center: fc)
	fce  = fs[jx+1][jy  ]; // #3 (center: fc)
	fcn  = fs[jx  ][jy+1]; // #4 (center: fc)
	fcs  = fs[jx  ][jy-1]; // #5 (center: fc)
	//
	fcww = fs[jx-2][jy  ]; // #6 (center: fcw)
	fcee = fs[jx+2][jy  ]; // #7 (center: fce)
	fcnn = fs[jx  ][jy+2]; // #8 (center: fcn)
	fcss = fs[jx  ][jy-2]; // #9 (center: fcs)
	//
	fcnw = fs[jx-1][jy+1]; // #10 (center: fcn)
	fcne = fs[jx+1][jy+1]; // #11 (center: fcn)
	fcsw = fs[jx-1][jy-1]; // #12 (center: fcs)
	fcse = fs[jx+1][jy-1]; // #13 (center: fcs)
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// term1 = Atomic_interaction*(1-2*f) + RT*{log(f) - log(1-f)}
	mu_chc = L0*(1.0-2.0*fcc) + RT*( log(fcc) - log(1.0-fcc) ); //center: fcc
	mu_chw = L0*(1.0-2.0*fcw) + RT*( log(fcw) - log(1.0-fcw) ); //center: fcw
	mu_che = L0*(1.0-2.0*fce) + RT*( log(fce) - log(1.0-fce) ); //center: fce
	mu_chn = L0*(1.0-2.0*fcn) + RT*( log(fcn) - log(1.0-fcn) ); //center: fcn
	mu_chs = L0*(1.0-2.0*fcs) + RT*( log(fcs) - log(1.0-fcs) ); //center: fcs
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// term2 = -gradient_energy_coefficient * Laplacian(f)
	mu_suc = -kapa_c*( (fce  + fcw  -2.0*fcc)/(dx*dx) + (fcn  + fcs  -2.0*fcc)/(dy*dy) ); //center: fcc
	mu_suw = -kapa_c*( (fcc  + fcww -2.0*fcw)/(dx*dx) + (fcnw + fcsw -2.0*fcw)/(dy*dy) ); //fcc=fcwe, fcnw=fcwn, fcsw=fcws, //center: fcw
	mu_sue = -kapa_c*( (fcee + fcc  -2.0*fce)/(dx*dx) + (fcne + fcse -2.0*fce)/(dy*dy) ); //fcc=fcew, fcne=fcen, fcse=fces, //center: fce
	mu_sun = -kapa_c*( (fcne + fcnw -2.0*fcn)/(dx*dx) + (fcnn + fcc  -2.0*fcn)/(dy*dy) ); //fcc=fcns, //center: fcn
	mu_sus = -kapa_c*( (fcse + fcsw -2.0*fcs)/(dx*dx) + (fcc  + fcss -2.0*fcs)/(dy*dy) ); //fcc=fcsn, //center: fcs
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// mu = dG/df = term1 + term2
	mu_c = mu_chc + mu_suc; // at current (jx,jy) grid point, //center: fcc
	mu_w = mu_chw + mu_suw; // at (jx-1,jy) grid point, //center: fcw
	mu_e = mu_che + mu_sue; // at (jx+1,jy) grid point, //center: fce
	mu_n = mu_chn + mu_sun; // at (jx,jy+1) grid point, //center: fcn
	mu_s = mu_chs + mu_sus; // at (jx,jy-1) grid point, //center: fcs
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// Laplacian(mu) = d^2(mu)/dx^2 + d^2(mu)/dy^2
	nab_mu = (mu_e + mu_w -2.0*mu_c)/(dx*dx)  // d^2(mu)/dx^2
		   + (mu_n + mu_s -2.0*mu_c)/(dy*dy); // d^2(mu)/dy^2
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// (df/dx) * d(mu)/dx, (x is related with e and w), (the center is fc.)
	dfmdx = ( (fce - fcw)/(2.0*dx) * (mu_e - mu_w)/(2.0*dx) );
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// (df/dy) * d(mu)/dy, (y is related with n and s), (the center is fc.)
	dfmdy = ( (fcn - fcs)/(2.0*dy) * (mu_n - mu_s)/(2.0*dy) );
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
	// df/dt = M*Laplacian(mu) + (dM/df)*( (df/dx) * d(mu)/dx + (df/dy) * d(mu)/dy )
	dfdt = mcc*nab_mu + dmc*(dfmdx + dfmdy); 
	fn[j] = f[j] + dfdt*dt;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
}

void update(float **f, float **fn)
{
	float *tmp = *f;
	*f  = *fn;
	*fn = tmp;
}

void write_vtk_grid_values_2D(int Nx, int Ny, float dx, float dy, int istep, float *data1);

int main(int argc, char** argv)
{
	float *f_d, *fn_d; // name of dynamic memory for GPU, CUDA, device
	float *F_h;        // name of dynamic memory for CPU
	int nx = NX, ny = NY;
	int times = TIMES;
	
	int nstep=10000;    //Number of time integration steps
	int nprint=1000;    //Output frequency to write the results to file
	//----- ----- ----- -----
	float Lx = 3.0e-07*times, // Simulation length in x-direction [micro m]
		  Ly = 3.0e-07*times, // Simulation length in y-direction [micro m]
		  //----- ----- ----- -----
		  dx = Lx/(float)nx, // Grid spacing between two grid pints in x-direction [nm]
		  dy = Ly/(float)ny, // Grid spacing between two grid pints in y-direction [nm]
		  //----- ----- ----- -----
		  c_0 = 0.4,    // Initial concentration of Cr (atomic fraction)
		  //----- ----- ----- -----
		  rr = 8.314,   // Gas constant [J/(mol*K)]
		  temp = 673.0, // Temperature [K]
		  RT = rr*temp,
		  //----- ----- ----- -----
		  L0 = 21020.8-9.31889*temp, // Atomic interaction [J/mol]
		  kapa_c = 1.2e-14,  // The value of gradient energy coefficients [J*m^2/mol]
		  //----- ----- ----- -----
		  Da = 1.0e-04*exp(-294000.0/RT), // Self-diffusion coefficient [m^2/s] (Fe)
		  Db = 2.0e-05*exp(-308000.0/RT), // Self-diffusion coefficient [m^2/s] (Cr)
		  //----- ----- ----- -----
		  dt = (dx*dx/Da)*0.1; // Time increment for the numerical integration [dimensionless]
	
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
	
	f_d  = (float *)malloc(nx*ny*sizeof(float)); //GPU, CUDA, device
	fn_d = (float *)malloc(nx*ny*sizeof(float)); //GPU, CUDA, device
	
	cudaMalloc((void**)&f_d ,nx*ny*sizeof(float)); // define dynamic memory for GPU (device)
	cudaMalloc((void**)&fn_d,nx*ny*sizeof(float)); // define dynamic memory for GPU (device)
	
	F_h  = (float *)malloc(nx*ny*sizeof(float));   // define dynamic memory for CPU (host)
	
	// Initialize the concentration filed F_h with random modulation
	for(int jy=0; jy<ny ; jy++){
		for(int jx=0; jx<nx ; jx++){
			int j = nx*jy + jx;
			float r = (float)rand()/(float)(RAND_MAX)-0.5;
			F_h[j] = c_0 + 0.01*r;
		}
	}//on CPU calculation
	
	//copy F_h(cpu,host) to f_d(cuda,device)
	cudaMemcpy(f_d,F_h,nx*ny*sizeof(float),cudaMemcpyHostToDevice);
	
	int bs=BS; // Number of threads
	dim3 blocks(nx/bs,ny/bs,1); //nx*ny = blocks * threads
	dim3 threads(bs,bs,1);      //bs*bs*1 <= 1024
	
	//----- ----- ----- -----start:(This part is not really needed.)----- ----- ----- ----
	//Set recording time
	cudaEvent_t start, stop;
	
	//Initialization
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	//Start recording time
	cudaEventRecord(start);
	//----- ----- ----- -----end:(This part is not really needed.)----- ----- ----- ----
	
	for(int istep=0; istep<=nstep ; istep++){
		//calculate subroutine "Kernel" on GPU
		Kernel<<<blocks, threads>>>(f_d,fn_d,nx,ny,rr,temp,L0,kapa_c,Da,Db,dt,dx,dy);
		cudaDeviceSynchronize(); //<-new version|old version -> cudaThreadSynchronize();
		
		// replace f_d with new f_d (=fn_d)
		update(&f_d,&fn_d);
		
		if(istep%nprint == 0){
			//copy f_d(cuda,device) to F_h(cpu,host)
			cudaMemcpy(F_h,f_d,nx*ny*sizeof(float),cudaMemcpyDeviceToHost);
			
			//output vtk format
			write_vtk_grid_values_2D(nx,ny,dx,dy,istep,F_h);
			
			//show current step
			fprintf(stderr,"nstep = %5d \n",istep);
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