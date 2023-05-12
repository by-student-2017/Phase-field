/* Program : 2D Phase-Field Simulation for 
   Spinodal Decomposition in Fe-Cr Alloy by GPU Computation.
   
   Programmer : Akinori Yamanaka (original version)
   Place : Depertment of Mechanical and Control Engineering Tokyo Institute of Technology
   Date : 7th, July, 2010 
   
   Modified version: By Student
   Place: 2-1 Hirosawa, Wako, Saitama, 351-0198, Japan
   Date: 12th/May/2023
   Test: Ubuntu 22.04 LTS
   
   Compling: nvcc -O2 main-gpu.cu write_vtk_grid_values_2D.cu -o main-gpu.exe -arch=native -lm --std 'c++17'
*/

#include <stdio.h>  //printf()
#include <stdlib.h> //rand() and malloc()
#include <math.h>   //mod() and -lm
//#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define BS 32        //Number of threads, 2^n=<32, BS*BS*1 <= 1024
#define TIMES 2
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
	float  da,
	float  db,
	float  dt,
	float  dx,
	float  dy
)
{
	int j, jx, jy;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	jx = blockDim.y*blockIdx.y + threadIdx.y; //<-GPU | CPU -> for(jx=0; jx<nx; jx++){
	jy = blockDim.x*blockIdx.x + threadIdx.x; //<-GPU | CPU -> for(jy=0; jy<ny; jy++){
	j  = nx*jy + jx;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	float  fcc,
		   fce,  fcw,  fcs,  fcn,
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
		   //----- ----- ----- ----- ----- -----
		   nab_mu, 
		   dfmdx, dfmdy, 
		   //----- ----- -----
		   dab = db/da, 
		   mcc, dmc,
		   //----- ----- -----
		   dfdt ;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* Consider up to the second nearest neighbor. Therefore, 
	13 difference grid points are used. (#1 to #13)
	   The if statement is used because of periodic boundary conditions. */
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #1 (center: fcc)
	fcc = f[j];
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #2 (center: fcc)
	if(jx == 0)    fcw = f[j+(nx-1)];    //boundary condition at west edge
	else           fcw = f[j+(  -1)];    //non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #3 (center: fcc)
	if(jx == nx-1) fce = f[j-(nx-1)];    //boundary condition at east edge
	else           fce = f[j-(  -1)];    //non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #4 (center: fcc)
	if(jy == 0)    fcs = f[j+nx*(ny-1)]; //boundary condition at south edge
	else           fcs = f[j+nx*(  -1)]; //non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #5 (center: fcc)
	if(jy == ny-1) fcn = f[j-nx*(ny-1)]; //boundary condition at north edge
	else           fcn = f[j-nx*(  -1)]; //non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #6 (center: fcn)
		 if(jx == 0 && jy == ny-1)   { fcnw = f[         nx-1];} // edge(west and north)
	else if(jx == 0 && jy  < ny-1)   { fcnw = f[j+nx    +nx-1];} // edge(west)
	else if(jx  > 0 && jy == ny-1)   { fcnw = f[j-nx*ny +nx-1];} // edge(north)
	else                             { fcnw = f[j       +nx-1];} // non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #7 (center: fcn)
		 if(jx == nx-1 && jy  < ny-1){ fcne = f[j-nx    +nx+1];} // edge(east)
	else if(jx  < nx-1 && jy == ny-1){ fcne = f[j-nx*ny +nx+1];} // edge(north)
	else if(jx == nx-1 && jy == ny-1){ fcne = f[            0];} // edge(east and north)
	else                             { fcne = f[j       +nx+1];} // non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #8 (center: fcs)
		 if(jx == 0 && jy >  0)      { fcsw = f[j+nx    -nx-1];} // edge(west)
	else if(jx  > 0 && jy == 0)      { fcsw = f[j+nx*ny -nx-1];} // edge(south)
	else if(jx == 0 && jy == 0)      { fcsw = f[      nx*ny-1];} // edge(west and south)
	else                             { fcsw = f[j       -nx-1];} // non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #9 (center: fcs)
		 if(jx == nx-1 && jy == 0)   { fcse = f[nx*ny-1 -nx+1];} // edge(east and south)
	else if(jx == nx-1 && jy  > 0)   { fcse = f[j-nx    -nx+1];} // edge(east)
	else if(jx <  nx-1 && jy == 0)   { fcse = f[j+nx*ny -nx+1];} // edge(south)
	else                             { fcse = f[j       -nx+1];} // non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #10 (center: fcw)
		 if(jx == 0)     { fcww = f[j+(nx-2)];}    // edge(west)
	else if(jx == 1)     { fcww = f[j+(nx-2)];}    // edge(west,one inside)
	else                 { fcww = f[j+(  -2)];}    // non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #11 (center: fce)
		 if(jx == nx - 2){ fcee = f[j-(nx-2)];}    // edge(east)
	else if(jx == nx - 1){ fcee = f[j-(nx-2)];}    // edge(east, one inside)
	else                 { fcee = f[j-(  -2)];}    // non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #12 (center: fcn)
		 if(jy == ny - 2){ fcnn = f[j-nx*(ny-2)];} // edge(north)
	else if(jy == ny - 1){ fcnn = f[j-nx*(ny-2)];} // edge(north, one inside)
	else                 { fcnn = f[j-nx*(  -2)];} // non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #13 (center: fcs)
		 if(jy == 0)     { fcss = f[j+nx*(ny-2)];} // edge(south)
	else if(jy == 1)     { fcss = f[j+nx*(ny-2)];} // edge(south, one inside)
	else                 { fcss = f[j+nx*(  -2)];} // non edge
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
	nab_mu = (mu_w + mu_e -2.0*mu_c)/(dx*dx)  // d^2(mu)/dx^2
		   + (mu_n + mu_s -2.0*mu_c)/(dy*dy); // d^2(mu)/dy^2
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// (df/dx) * d(mu)/dx, (x is related with w and e), (the center is fc.)
	dfmdx = ( (fcw - fce)/(2.0*dx) * (mu_w - mu_e)/(2.0*dx) );
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// (df/dy) * d(mu)/dy, (y is related with n and s), (the center is fc.)
	dfmdy = ( (fcn - fcs)/(2.0*dy) * (mu_n - mu_s)/(2.0*dy) );
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// Mobility, M = { (D_A/RT)*c + (D_B/RT)*(1-c) }*c*(1-c)
	//             = (D_a/RT)*{f + (D_B/D_A)*(1-f)}*f*(1-f)
	mcc = (da/RT)*(fcc+dab*(1.0-fcc))*fcc*(1.0-fcc); 
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// dM/df
	dmc = (da/RT)*((1.0-dab)*fcc*(1.0-fcc)
		+ (fcc+dab*(1.0-fcc))*(1.0-2.0*fcc)); 
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// df/dt = M*Laplacian(f) + (dM/df)*( (df/dx) * d(mu)/dx + (df/dy) * d(mu)/dy )
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
	
	int bs=BS; // Number of threads, 16 or 32
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
			fprintf(stderr,"istep = %5d \n",istep);
		}
		//
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
