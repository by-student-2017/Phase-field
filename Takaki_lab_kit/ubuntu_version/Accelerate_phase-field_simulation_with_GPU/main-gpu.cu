/* Program : 2D Phase-Field Simulation for 
   Spinodal Decomposition in Fe-Cr Alloy by GPU Computation.
   
   Programmer : Akinori Yamanaka (original version)
   Place : Depertment of Mechanical and Control Engineering Tokyo Institute of Technology
   Date : 7th, July, 2010 
   
   Modified version: By Student
   Place: 2-1 Hirosawa, Wako, Saitama, 351-0198, Japan
   Date: 11th/May/2023
*/

#include <stdio.h>  //printf()
#include <stdlib.h> //rand() and malloc()
#include <math.h>   //mod() and -lm
//#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define NX 256 //Number of grid points in the x-direction
#define NY 256 //Number of grid points in the y-direction

// Define subroutine "Kernel" for GPU (Device) calculation in detail
__global__ void Kernel
(
	float *f,float *fn, int nx, int ny,
	float  rr, float temp,float L0,
	float kapa_c,float da,float db,float dt,float dx,float dy
)
{
	int j, jx, jy;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	jx = blockDim.y*blockIdx.y + threadIdx.y;
	jy = blockDim.x*blockIdx.x + threadIdx.x;
	j  = nx*jy + jx;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	float fcc, fce, fcw, fcs, fcn, fcnw, fcne, fcsw, fcse, fcww, fcee, fcnn, fcss, 
		  mu_chc, mu_chw, mu_che, mu_chn, mu_chs, mu_suc, mu_suw, mu_sue, mu_sun, 
		  mu_sus, mu_c, mu_w, mu_e, mu_n, mu_s,  
		  nab_mu, dfmdx, dfmdy, dab = db/da, mcc, dmc, dfdt ;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* Consider up to the second nearest neighbor. Therefore, 
	13 difference grid points are used. (#1 to #13)
	   The if statement is used because of periodic boundary conditions.*/
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #1
	fcc = f[j]; // The value of concentration at current (jx,jy) grid point
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #2
	if(jx == 0)    fcw = f[j+nx-1];
	else           fcw = f[j   -1];
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #3
	if(jx == nx-1) fce = f[j-nx+1];
	else           fce = f[j   +1];
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #4
	if(jy == 0)    fcs = f[j+nx*ny-nx];
	else           fcs = f[j      -nx];
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #5
	if(jy == ny-1) fcn = f[j-nx*ny+nx];
	else           fcn = f[j      +nx];
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #6
	if(jx == 0 && jy == ny-1)     { fcnw = f[nx-1];}
	else if(jx == 0 && jy  < ny-1){ fcnw = f[j+nx    +nx-1];}
	else if(jx  > 0 && jy == ny-1){ fcnw = f[j-nx*ny +nx-1];}
	else                          { fcnw = f[j       +nx-1];}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #7
	if(jx == nx-1 && jy  < ny-1)     { fcne = f[j-nx    +nx+1];}
	else if(jx  < nx-1 && jy == ny-1){ fcne = f[j-nx*ny +nx+1];}
	else if(jx == nx-1 && jy == ny-1){ fcne = f[0];}
	else                             { fcne = f[j       +nx+1];}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #8
	if(jx == 0 && jy >  0)     { fcsw = f[j+nx    -nx-1];}
	else if(jx  > 0 && jy == 0){ fcsw = f[j+nx*ny -nx-1];}
	else if(jx == 0 && jy == 0){ fcsw = f[nx*ny-1];} 
	else                       { fcsw = f[j       -nx-1];}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #9
	if(jx == nx-1 && jy == 0)     {fcse = f[nx*ny-1-nx+1];}
	else if(jx == nx-1 && jy  > 0){ fcse = f[j-nx    -nx+1];}
	else if(jx <  nx-1 && jy == 0){ fcse = f[j+nx*ny -nx+1];}
	else                          { fcse = f[j       -nx+1];}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #10
	if(jx == 0) {fcww = f[j+nx-2];}
	else if(jx == 1){ fcww = f[j+nx-2];}
	else            { fcww = f[j   -2];}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #11
	if(jx == nx - 2)     { fcee = f[j-nx+2];}
	else if(jx == nx - 1){ fcee = f[j-nx+2];}
	else                 { fcee = f[j   +2];}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #12
	if(jy == ny - 2)      { fcnn = f[j-nx*ny+nx+nx];}
	else if(jy == ny - 1) { fcnn = f[j-nx*ny+nx+nx];}
	else                  { fcnn = f[j      +nx+nx];}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #13
	if(jy == 0)      { fcss = f[j+nx*ny-nx-nx];}
	else if(jy == 1) { fcss = f[j+nx*ny-nx-nx];}
	else             { fcss = f[j      -nx-nx];}
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
	dfmdy = ((mu_n-mu_s)*(fcn-fcs))/(4.0*dy*dy);
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
	
	int nstep=10000;    //Number of time integration steps
	int nprint=1000;    //Output frequency to write the results to file
	float Lx = 3.0e-07, // Simulation length in x-direction [micro m]
		  Ly = 3.0e-07, // Simulation length in y-direction [micro m]
		  dx = Lx/(float)nx, // Grid spacing between two grid pints in x-direction [nm]
		  dy = Ly/(float)ny, // Grid spacing between two grid pints in y-direction [nm]
		  c_0 = 0.4,    // Initial concentration (atomic fraction)
		  rr = 8.314,   // Gas constant [J/(mol*K)]
		  temp = 673.0, // Temperature [K]
		  L0 = 21020.8-9.31889*temp, // Atomic interaction [J/mol]
		  kapa_c = 1.2e-14,  // The value of gradient energy coefficients [J*m^2/mol]
		  da = 1.0e-04*exp(-294000.0/rr/temp), // Self-diffusion coefficient [m^2/s] (Fe)
		  db = 2.0e-05*exp(-308000.0/rr/temp), // Self-diffusion coefficient [m^2/s] (Cr)
		  dt = (dx*dx/da)*0.1; // Time increment for the numerical integration [dimensionless]
	
	//CUT_DEVICE_INIT(argc, argv);
	
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
	
	int BS=16; // Number of threads, 16 or 32
	dim3 blocks(nx/BS,ny/BS,1); //nx*ny = blocks * threads
	dim3 threads(BS,BS,1);      //BS*BS*1 <= 1024
	
	//unsigned int timer;
	//cutCreateTimer(&timer);
	//cutResetTimer(timer);
	//cutStartTimer(timer);
	
	for(int istep=0; istep<=nstep ; istep++){
		//calculate subroutine "Kernel" on GPU
		Kernel<<<blocks, threads>>>(f_d,fn_d,nx,ny,rr,temp,L0,kapa_c,da,db,dt,dx,dy);
		//cudaThreadSynchronize();
		
		// replace f_d with new f_d (=fn_d)
		update(&f_d,&fn_d);
		//
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
	//cutStopTimer(timer);
	//float calc_time = cutGetTimerValue(timer)*1.0e-03;
	//printf("Calculation Time = %9.3e [sec]\n",calc_time);
	
	cudaFree(f_d);
	cudaFree(fn_d);
	
	free(F_h);
	
	return 0;
}
