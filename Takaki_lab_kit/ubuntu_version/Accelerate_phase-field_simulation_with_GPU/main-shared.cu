/* Program : 2D Phase-Field Simulation for 
   Spinodal Decomposition in Fe-Cr Alloy by GPU (shared memory version) Computation.
   
   Programmer : Akinori Yamanaka (original version)
   Place : Depertment of Mechanical and Control Engineering Tokyo Institute of Technology
   Date : 7th, July, 2010 
   
   Modified version: By Student
   Place: 2-1 Hirosawa, Wako, Saitama, 351-0198, Japan
   Date: 11th/May/2023
   Test: Ubuntu 22.04 LTS
   
   Compling: nvcc -O2 main-shared.cu write_vtk_grid_values_2D.cu -o main-shared.exe -arch=native -lm --std 'c++17'
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define TIMES 2
#define NX 256*TIMES //Number of grid points in the x-direction
#define NY 256*TIMES //Number of grid points in the y-direction

// Define subroutine "Kernel" for GPU (Device) calculation in detail
__global__ void Kernel
(
	float *f, float *fn, int nx, int ny,
	float  rr, float temp, float L0,
	float  kapa_c, float da, float db, float dt, float dx, float dy
)
{
	int j, jx, jy;
	float  fcc, fce, fcw, fcs, fcn, fcnw, fcne, fcsw, fcse, fcww, fcee, fcnn, fcss, 
		   mu_chc, mu_chw, mu_che, mu_chn, mu_chs,
		   mu_suc, mu_suw, mu_sue, mu_sun, mu_sus, 
		   mu_c, mu_w, mu_e, mu_n, mu_s, 
		   nab_mu, dfmdx, dfmdy, dab = db/da, mcc, dmc, dfdt ;
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// Declare variables used when copying data from global memory to shared memory
	int joff;
	int J0, J1, J2, J3; // One inner edge
	int J4, J5, J6, J7; // The most edge
	int J8, J9, J10, J11; // Corner
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	//int thread_x = 16, thread_y = 16; // 16=BS (16 kB before GF100 Core, 48 kB after GF100 Core)
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* blockDim.x * blockDim.y = 16 * 16. In addition, 
	   add the necessary two adjacent difference grid points to
	   both sides of the x-axis and y-axis, respectively. */
	__shared__ float fs[2+16+2][2+16+2]; //fs is shared memory
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	/* Note ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	   Shared memory: (2+blockDim.x+2) * (2+blockDim.y+2)
	     jx = threadIdx.x + 2; // +2 when counting from one end
	     jy = threadIdx.y + 2; // +2 when counting from one end
	   ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	   Global memory: (gridDim.x*blockDim.x) * (gridDim.y*blockDim.y)
	     x = (blockIdx.x*blockDim.x + threadIdx.x)
	     y = (blockIdx.y*blockDim.y + threadIdx.y)
	----- ----- ----- ----- ----- ----- ----- ----- ----- -----   */
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	jx = threadIdx.x + 2; // +2 when counting from one end
	jy = threadIdx.y + 2; // +2 when counting from one end
	joff = nx*(blockDim.y*blockIdx.y) + blockDim.x*blockIdx.x;
	j = joff + nx*threadIdx.y + threadIdx.x;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	fcc = f[j]; // Global memory at current (jx,jy) grid point
	fs[jx][jy] = fcc; // Global memory to shared memory at current (jx,jy) grid point
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// Calculating sleeve area in shared memory
	//----- ----- ----- Upper sleeve area
	if(blockIdx.y == 0) {J0 = nx*(ny-1)+blockDim.x*blockIdx.x + threadIdx.x,
						 J4 = J0 - nx;} 
	else                {J0 =  j - nx, 
						 J4 = J0 - nx;}
	//----- ----- ----- Lower sleeve area
	if(blockIdx.y == gridDim.y - 1) {J1 = blockDim.x*blockIdx.x + threadIdx.x, 
						 J5 = J1 + nx;} 
	else				{J1 =  j + nx, 
						 J5 = J1 + nx;}
	//----- ----- ----- Left sleeve area
	if(blockIdx.x == 0) {J2 = joff + nx*threadIdx.x + nx - 1,
						 J6 = J2 - 1;}
	else				{J2 = joff + nx*threadIdx.x - 1, 
						 J6 = J2 - 1;}
	//----- ----- ----- Right sleeve area
	if(blockIdx.x == gridDim.x - 1) {J3 = joff + nx*threadIdx.x + 15 - nx + 1,
						 J7 = J3 + 1;}
	else				{J3 = joff + nx*threadIdx.x + 16,
						 J7 = J3 + 1;}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	//----- ----- ----- Upper and Rgith sleeve area
		 if(blockIdx.x == 0 && blockIdx.y == gridDim.y - 1) { J8 = blockDim.x*16 -1;}
	else if(blockIdx.x  > 0 && blockIdx.y == gridDim.y - 1) { J8 = J1 - 1 ;}
	else if(blockIdx.x == 0 && blockIdx.y  < gridDim.y - 1) { J8 = j + nx + nx -1;}
	else                                                    { J8 = j + nx -1 ;}
	//----- ----- ----- Lower and Rgith sleeve area
		 if(blockIdx.x == gridDim.x - 1 && blockIdx.y == gridDim.y - 1) { J9 = 0 ;}
	else if(blockIdx.x  < gridDim.x - 1 && blockIdx.y == gridDim.y - 1) { J9 = J1 + 1 ;}
	else if(blockIdx.x == gridDim.x - 1 && blockIdx.y  < gridDim.y - 1) { J9 = j  + 1 ;}
	else                                                                { J9 = j + nx +1 ;}
	//----- ----- ----- Upper and Left sleeve area
		 if(blockIdx.x  > 0 && blockIdx.y == 0) { J10 = J0 - 1 ;}
	else if(blockIdx.x == 0 && blockIdx.y  > 0) { J10 =  j -1  ;}
	else if(blockIdx.x == 0 && blockIdx.y == 0) { J10 = nx*blockDim.x*blockDim.y - 1 ;}
	else                                        { J10 = j - nx - 1 ;}
	//----- ----- ----- Lower and Left sleeve area
		 if(blockIdx.x == gridDim.x -1 && blockIdx.y == 0) { J11 = nx*blockDim.x*blockDim.y -1 - nx + 1;}
	else if(blockIdx.x  < gridDim.x -1 && blockIdx.y == 0) { J11 = J0 + 1  ;}
	else if(blockIdx.x == gridDim.x -1 && blockIdx.y  > 0) { J11 =  j - nx - nx + 1 ;}
	else                                                   { J11 = j - nx + 1 ;}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// copy Global memory to Shared memory
	if(threadIdx.y ==  0){ fs[jx][ 1] = f[J0], fs[jx][ 0] = f[J4];}   // Upper sleeve area
	if(threadIdx.y ==  1){ fs[ 1][jx] = f[J2], fs[ 0][jx] = f[J6];}   // Left  sleeve area
	if(threadIdx.y ==  2){ fs[18][jx] = f[J3], fs[19][jx] = f[J7];}   // Right sleeve area
	if(threadIdx.y == 15){ fs[jx][18] = f[J1], fs[jx][19] = f[J5];}   // Lower sleeve area
	//----- ----- ----- 
	if(threadIdx.x ==  0 && threadIdx.y == 15) {fs[ 1][18] = f[J8];}  // Upper and Rgith sleeve area
	if(threadIdx.x == 15 && threadIdx.y == 15) {fs[18][18] = f[J9];}  // Lower and Rgith sleeve area
	if(threadIdx.x ==  0 && threadIdx.y ==  0) {fs[ 1][ 1] = f[J10];} // Upper and Left  sleeve area
	if(threadIdx.x == 15 && threadIdx.y ==  0) {fs[18][ 1] = f[J11];} // Lower and Left  sleeve area
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// Wait until all data is secured
	__syncthreads();
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* Consider up to the second nearest neighbor. Therefore, 
	   13 difference grid points are used. (#1 to #13) */
	fcc  = fs[jx  ][jy  ]; // #1
	fcw  = fs[jx-1][jy  ]; // #2
	fce  = fs[jx+1][jy  ]; // #3
	fcn  = fs[jx  ][jy+1]; // #4
	fcs  = fs[jx  ][jy-1]; // #5
	//
	fcww = fs[jx-2][jy  ]; // #6
	fcee = fs[jx+2][jy  ]; // #7
	fcnn = fs[jx  ][jy+2]; // #8
	fcss = fs[jx  ][jy-2]; // #9
	//
	fcnw = fs[jx-1][jy+1]; // #10
	fcne = fs[jx+1][jy+1]; // #11
	fcsw = fs[jx-1][jy-1]; // #12
	fcse = fs[jx+1][jy-1]; // #13
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// term1 = Atomic_interaction*(1-2*f) + RT*{log(f) - log(1-f)}
	mu_chc = L0*(1.0-2.0*fcc)+rr*temp*(log(fcc)-log(1.0-fcc));
	mu_chw = L0*(1.0-2.0*fcw)+rr*temp*(log(fcw)-log(1.0-fcw));
	mu_che = L0*(1.0-2.0*fce)+rr*temp*(log(fce)-log(1.0-fce));
	mu_chn = L0*(1.0-2.0*fcn)+rr*temp*(log(fcn)-log(1.0-fcn));
	mu_chs = L0*(1.0-2.0*fcs)+rr*temp*(log(fcs)-log(1.0-fcs));
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// term2 = -gradient_energy_coefficient * Laplacian(f)
	mu_suc = -kapa_c*(fce +fcw +fcn +fcs -4.0*fcc)/dx/dx;
	mu_suw = -kapa_c*(fcc +fcww+fcnw+fcsw-4.0*fcw)/dx/dx;
	mu_sue = -kapa_c*(fcee+fcc +fcne+fcse-4.0*fce)/dx/dx;
	mu_sun = -kapa_c*(fcne+fcnw+fcnn+fcc -4.0*fcn)/dx/dx;
	mu_sus = -kapa_c*(fcse+fcsw+fcc +fcss-4.0*fcs)/dx/dx;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// mu = dG/df = term1 + term2
	mu_c = mu_chc + mu_suc; // at current (jx,jy) grid point
	mu_w = mu_chw + mu_suw; // at (jx-1,jy) grid point
	mu_e = mu_che + mu_sue; // at (jx+1,jy) grid point
	mu_n = mu_chn + mu_sun; // at (jx,jy+1) grid point
	mu_s = mu_chs + mu_sus; // at (jx,jy-1) grid point
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// Laplacian(mu) = d^2(mu)/dx^2 + d^2(mu)/dy^2
	nab_mu = (mu_w + mu_e + mu_n + mu_s -4.0*mu_c)/dx/dx;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// (df/dx) * d(mu)/dx
	dfmdx = ((mu_w-mu_e)*(fcw-fce))/(4.0*dx*dx);
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// (df/dy) * d(mu)/dy
	dfmdy = ((mu_n-mu_s)*(fcn-fcs))/(4.0*dx*dx);
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// Mobility, M = { (D_A/RT)*c + (D_B/RT)*(1-c) }*c*(1-c)
	//             = (D_a/RT)*{f + (D_B/D_A)*(1-f)}*f*(1-f)
	mcc = (da/rr/temp)*(fcc+dab*(1.0-fcc))*fcc*(1.0-fcc); 
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// dM/df
	dmc = (da/rr/temp)*((1.0-dab)*fcc*(1.0-fcc)
					  +(fcc+dab*(1.0-fcc))*(1.0-2.0*fcc));
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// df/dt = M*Laplacian(f) + (dM/df)*( (df/dx) * d(mu)/dx + (df/dy) * d(mu)/dy )
	dfdt = mcc*nab_mu + dmc*(dfmdx+dfmdy); 
	fn[j] = f[j]+dfdt*dt;
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
		  c_0 = 0.4,    // Initial concentration (atomic fraction)
		  //----- ----- ----- -----
		  rr = 8.314,   // Gas constant [J/(mol*K)]
		  temp = 673.0, // Temperature [K]
		  //----- ----- ----- -----
		  L0 = 21020.8-9.31889*temp, // Atomic interaction [J/mol]
		  kapa_c = 1.2e-14,  // The value of gradient energy coefficients [J*m^2/mol]
		  //----- ----- ----- -----
		  da = 1.0e-04*exp(-294000.0/rr/temp), // Self-diffusion coefficient [m^2/s] (Fe)
		  db = 2.0e-05*exp(-308000.0/rr/temp), // Self-diffusion coefficient [m^2/s] (Cr)
		  //----- ----- ----- -----
		  dt = (dx*dx/da)*0.1; // Time increment for the numerical integration [dimensionless]
	
	f_d  = (float *)malloc(nx*ny*sizeof(float)); //GPU, CUDA, device
	fn_d = (float *)malloc(nx*ny*sizeof(float)); //GPU, CUDA, device
	
	cudaMalloc((void**)&f_d ,nx*ny*sizeof(float)); // define dynamic memory for GPU (device)
	cudaMalloc((void**)&fn_d,nx*ny*sizeof(float)); // define dynamic memory for GPU (device)
	
	F_h  = (float *)malloc(nx*ny*sizeof(float));   // define dynamic memory for CPU (host)
	
	// Initialize the concentration filed F_h with random modulation
	for(int jy=0; jy<ny ; jy++){
		for(int jx=0; jx<nx ; jx++){
			int j = nx*jy + jx;
			float r = (float)rand()/(float)(RAND_MAX);
			F_h[j] = c_0 + 0.01*r;
		}
	}//on CPU calculation
	
	//copy F_h(cpu,host) to f_d(cuda,device)
	cudaMemcpy(f_d,F_h,nx*ny*sizeof(float),cudaMemcpyHostToDevice);
	
	int bs=16; // Number of threads, 16 only
	dim3 blocks(nx/bs,ny/bs,1); //nx*ny = blocks * threads
	dim3 threads(bs,bs,1);      //bs*bs*1 <= 1024
	
	//----- ----- ----- -----
	//Set recording time
	cudaEvent_t start, stop;
	
	//Initialization
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	//Start recording time
	cudaEventRecord(start);
	//----- ----- ----- -----
	
	for(int istep=0; istep<=nstep ; istep++){
		//calculate subroutine "Kernel" on GPU
		Kernel<<<blocks, threads>>>(f_d,fn_d,nx,ny,rr,temp,L0,kapa_c,da,db,dt,dx,dy);
		//cudaThreadSynchronize();
		
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
	
	//----- ----- ----- -----
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
	//----- ----- ----- -----
	
	cudaFree(f_d);
	cudaFree(fn_d);
	
	free(F_h);
	
	return 0;
}