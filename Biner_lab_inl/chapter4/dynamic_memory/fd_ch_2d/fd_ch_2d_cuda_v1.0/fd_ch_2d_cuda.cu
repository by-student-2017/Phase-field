/* Finite-difference phase-field code for
   solving Cahn-Hilliard equation */

/* This program solves the Cahn-Hilliard phase-field
   equation with finite difference algorithm using
   five-point stencil. The time integration is
   carried out using explicit Euler scheme. */

#include <stdio.h> //printf()
#include <stdlib.h> //rand()
#include <math.h> //mod() and -lm
#include <time.h>
//----- ----- -----
#include <cuda.h>   //GPU
/* #include <cuda.h> or
  #include "cuda_runtime.h"
  #include "device_launch_parameters.h" */
//----- ----- -----

#define BS 32        //Number of threads, 2^n=<32, BS*BS*1 <= 1024
#define TIMES 1
//----- ----- -----
#define NX 64*TIMES //Number of grid points in the x-direction
#define NY 64*TIMES //Number of grid points in the y-direction

//----- ----- ----- ----- ----- ----- -----
void micro_ch_pre_2d(int Nx, int Ny, float c0, float *con);
//----- ----- -----
float calculate_energy_2d(int Nx, int Ny, float *con, float grad_coef);
//----- ----- -----
void write_vtk_grid_values_2D(int nx, int ny, 
	float dx, float dy,
	int istep, float *data1);
//----- ----- ----- ----- ----- ----- -----

// Define subroutine "Kernel" for GPU (Device) calculation in detail
__global__ void Kernel
(
	float *f,
	float *fn,
	int    nx,
	int    ny,
	float  dx,
	float  dy,
	float  dtime,
	float  mobility,
	float  grad_coef
)
{
	int j, jx, jy;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	jx = blockDim.x*blockIdx.x + threadIdx.x; //<-GPU | CPU -> for(jx=0; jx<nx; jx++){
	jy = blockDim.y*blockIdx.y + threadIdx.y; //<-GPU | CPU -> for(jy=0; jy<ny; jy++){
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
		   mu_chc,
		   mu_chw, mu_che, mu_chn, mu_chs,
		   //----- ----- -----
		   mu_suc,
		   mu_suw, mu_sue, mu_sun, mu_sus, 
		   //----- ----- -----
		   mu_c,
		   mu_w, mu_e, mu_n, mu_s, 
		   //----- ----- ----- ----- ----- -----
		   lap_mu, 
		   //----- ----- -----
		   dfdt ;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* Consider up to the second nearest neighbor. Therefore, 
	13 difference grid points are used. (#1 to #13)
	   The if statement is used because of periodic boundary conditions. */
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #1 (center: fcc)
	fcc = f[j];
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #2 (center: fcc, and + w)
	if(jx == 0)    fcw = f[j+nx-1];    //boundary condition at west edge
	else           fcw = f[j   -1];    //non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #3 (center: fcc, and + e)
	if(jx == nx-1) fce = f[j-nx+1];    //boundary condition at east edge
	else           fce = f[j   +1];    //non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #4 (center: fcc, and + s)
	if(jy == 0)    fcs = f[j+nx*(+ny-1)]; //boundary condition at south edge
	else           fcs = f[j+nx*(   -1)]; //non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #5 (center: fcc, and + n)
	if(jy == ny-1) fcn = f[j+nx*(-ny+1)]; //boundary condition at north edge
	else           fcn = f[j+nx*(   +1)]; //non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #6 (center: fcc, and + n + w)(#5 and #2)
	/* e.g., "if(jx == 0 && jy == ny-1)" is f[ (j+nx*(-ny+1)) + (j+nx-1) - j] = f[j + nx*(-ny+1) +  nx-1] using above #5 and #2 condition.
	         "if(jx == 0 && jy  < ny-1)" is f[ j+nx*(   +1) + (j+nx-1) - j]   = f[j + nx*(   +1) +  nx-1] using above #5 and #2 condition. */
		 if(jx == 0 && jy == ny-1)   { fcnw = f[         nx-1];} // =f[j + nx*(-ny+1) +  nx-1] = f[nx*(ny-1)-nx*(ny-1)+nx-1] = f[nx-1]
	else if(jx == 0 && jy  < ny-1)   { fcnw = f[j+nx    +nx-1];} // =f[j + nx*(   +1) +  nx-1] = f[j+nx    +nx-1]
	else if(jx  > 0 && jy == ny-1)   { fcnw = f[j-nx*ny +nx-1];} // =f[j + nx*(-ny+1) +    -1] = f[j-nx*ny +nx-1]
	else                             { fcnw = f[j       +nx-1];} // =f[j + nx*(   +1) +    -1] = f[j       +nx-1]
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #7 (center: fcc, and + n + e)(#5 and #3)
		 if(jx == nx-1 && jy  < ny-1){ fcne = f[j+nx    -nx+1];} // =f[j + nx*(   +1) + -nx+1] = f[j+nx    -nx+1]
	else if(jx  < nx-1 && jy == ny-1){ fcne = f[j-nx*ny +nx+1];} // =f[j + nx*(-ny+1) +     1] = f[j-nx*ny +nx+1]
	else if(jx == nx-1 && jy == ny-1){ fcne = f[            0];} // =f[j + nx*(-ny+1) + -nx+1] = f[nx*(ny-1)+nx-1+nx*(-ny+1)-(nx-1)]=f[0]
	else                             { fcne = f[j       +nx+1];} // =f[j + nx*(   +1) +     1] = f[j       +nx+1]
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #8 (center: fcc, and + s + w)(#4 and #2)
		 if(jx == 0 && jy >  0)      { fcsw = f[j-nx    +nx-1];} // =f[j + nx*(   -1) +  nx-1] = f[j-nx    +nx-1]
	else if(jx  > 0 && jy == 0)      { fcsw = f[j+nx*ny -nx-1];} // =f[j + nx*(+ny-1) +    -1] = f[j+nx*ny -nx-1]
	else if(jx == 0 && jy == 0)      { fcsw = f[      nx*ny-1];} // =f[j + nx-1       + nx*(+ny-1)] = f[j+nx-1+nx*ny-nx] (and j= nx*jy + jx= nx*0 + 0 = 0)
	else                             { fcsw = f[j       -nx-1];} // =f[j + nx*(   -1) +    -1] = f[j       -nx-1]
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #9 (center: fcc, and + s + e)(#4 and #3)
		 if(jx == nx-1 && jy == 0)   { fcse = f[nx*ny-1 -nx+1];} // =f[j + nx*(+ny-1) + -nx+1] = f[nx-1+nx*ny-nx -nx+1]
	else if(jx == nx-1 && jy  > 0)   { fcse = f[j-nx    -nx+1];} // =f[j + nx*(   -1) + -nx+1] = f[j-nx    -nx+1]
	else if(jx <  nx-1 && jy == 0)   { fcse = f[j+nx*ny -nx+1];} // =f[j + nx*(+ny-1) +     1] = f[j+nx*ny -nx+1]
	else                             { fcse = f[j       -nx+1];} // =f[j + nx*(   -1) +     1] = f[j       -nx+1]
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #10 (center: fcc, and + w + w)
		 if(jx == 0)     { fcww = f[j+nx-2];}    // edge(west)
	else if(jx == 1)     { fcww = f[j+nx-2];}    // edge(west,one inside)
	else                 { fcww = f[j   -2];}    // non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #11 (center: fcc, and + e + e)
		 if(jx == nx - 2){ fcee = f[j-nx+2];}    // edge(east)
	else if(jx == nx - 1){ fcee = f[j-nx+2];}    // edge(east, one inside)
	else                 { fcee = f[j   +2];}    // non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #12 (center: fcc, and + n + n)
		 if(jy == ny - 2){ fcnn = f[j+nx*(-ny+2)];} // edge(north)
	else if(jy == ny - 1){ fcnn = f[j+nx*(-ny+2)];} // edge(north, one inside)
	else                 { fcnn = f[j+nx*(   +2)];} // non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #13 (center: fcc, and + s + s)
		 if(jy == 0)     { fcss = f[j+nx*(+ny-2)];} // edge(south)
	else if(jy == 1)     { fcss = f[j+nx*(+ny-2)];} // edge(south, one inside)
	else                 { fcss = f[j+nx*(   -2)];} // non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// term1 = dfdcon[i][j]
	/* float A=1.0;
	dfdcon[i][j]=A*( 2.0*con[i][j]*(1.0-con[i][j])*(1.0-con[i][j])
				    -2.0*con[i][j]*con[i][j]*(1.0-con[i][j]) );  */
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	mu_chc = 1.0*( 2.0*fcc*(1.0-fcc)*(1.0-fcc) -2.0*fcc*fcc*(1.0-fcc) ); //center: fcc
	mu_chw = 1.0*( 2.0*fcw*(1.0-fcw)*(1.0-fcw) -2.0*fcw*fcw*(1.0-fcw) ); //center: fcw
	mu_che = 1.0*( 2.0*fce*(1.0-fce)*(1.0-fce) -2.0*fce*fce*(1.0-fce) ); //center: fce
	mu_chn = 1.0*( 2.0*fcn*(1.0-fcn)*(1.0-fcn) -2.0*fcn*fcn*(1.0-fcn) ); //center: fcn
	mu_chs = 1.0*( 2.0*fcs*(1.0-fcs)*(1.0-fcs) -2.0*fcs*fcs*(1.0-fcs) ); //center: fcs
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// term2 = -grad_coef * Laplacian(f)
	mu_suc = -grad_coef*( (fce  + fcw  -2.0*fcc)/(dx*dx) + (fcn  + fcs  -2.0*fcc)/(dy*dy) ); //center: fcc
	mu_suw = -grad_coef*( (fcc  + fcww -2.0*fcw)/(dx*dx) + (fcnw + fcsw -2.0*fcw)/(dy*dy) ); //fcc=fcwe, fcnw=fcwn, fcsw=fcws, //center: fcw
	mu_sue = -grad_coef*( (fcee + fcc  -2.0*fce)/(dx*dx) + (fcne + fcse -2.0*fce)/(dy*dy) ); //fcc=fcew, fcne=fcen, fcse=fces, //center: fce
	mu_sun = -grad_coef*( (fcne + fcnw -2.0*fcn)/(dx*dx) + (fcnn + fcc  -2.0*fcn)/(dy*dy) ); //fcc=fcns, //center: fcn
	mu_sus = -grad_coef*( (fcse + fcsw -2.0*fcs)/(dx*dx) + (fcc  + fcss -2.0*fcs)/(dy*dy) ); //fcc=fcsn, //center: fcs
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* Evaluate the terms in Eq.4.17, and
	   accumulate them in array dummy[Nx][Ny] for
	   each grid point in the simulation cell
	   dummy[ij]=dfdcon-grad_coef*lap_con[ij];                  */
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// mu = dG/df = term1 + term2
	mu_c = mu_chc + mu_suc; // at current (jx,jy) grid point, //center: fcc
	mu_w = mu_chw + mu_suw; // at (jx-1,jy) grid point, //center: fcw
	mu_e = mu_che + mu_sue; // at (jx+1,jy) grid point, //center: fce
	mu_n = mu_chn + mu_sun; // at (jx,jy+1) grid point, //center: fcn
	mu_s = mu_chs + mu_sus; // at (jx,jy-1) grid point, //center: fcs
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* Calculate the Laplacing of the terms inside
	   the parenthesis in Eq.4.16
	   lap_dummy[ij] = (hne + hnw -2.0*hnc)/(dx*dx)
	                  +(hns + hnn -2.0*hnc)/(dy*dy);           */
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// Laplacian(mu) = d^2(mu)/dx^2 + d^2(mu)/dy^2
	lap_mu = (mu_e + mu_w -2.0*mu_c)/(dx*dx)  // d^2(mu)/dx^2
		   + (mu_n + mu_s -2.0*mu_c)/(dy*dy); // d^2(mu)/dy^2
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* Explicit Euler time integration of
	   concentration field, Eq.4.16
	   con[ij] = con[ij] + dtime*mobility*lap_dummy[ij];       */
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// df/dt = M*Laplacian(mu)
	dfdt = mobility*lap_mu;
	fn[j] = f[j] + dfdt*dtime;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//}//end for(jy <-GPU | CPU-> }//end for(jy
	//}//end for(jx <-GPU | CPU-> }//end for(jx
}

//f=con_d, fn=new_con_d
void update(float **f, float **fn)
{
	float *tmp = *f;
	*f  = *fn;
	*fn = tmp;
}

int main(){
	//open an output file for writing total bulk energy values
	FILE *out2=fopen("time_energy.out","w");
	
	//simulation cell parameters
	int nx=NX; //Number of grid points in the x-direction
	int ny=NY; //Number of grid points in the y-direction
	
	//The distance between two grid points in x,y-direction
	float dx=1.0; //Grid spacing between two grid pints in x-direction
	float dy=1.0; //Grid spacing between two grid pints in y-direction
	
	//time integration parameters
	int nstep=10000; //Number of time integration steps
	int nprint=50;  //Output frequency to write the results to file
	float dtime=1.e-2; //Time increment for the numerical integration
	float ttime=0.0;    //Total time
	
	//material specific parameters
	
	//Initial concentrations of alloying elements
	float c0=0.40;       //Average composition
	float mobility=1.0;  //The value of mobility coefficient (dimensionless)
	float grad_coef=0.5; //The value of gradient energy coefficients [J(nm)^2/mol]
	
	float energy;
	
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
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	//float       con[Nx][Ny]; //concentraion
	float       *con = (float *)malloc(sizeof(float)*( nx*ny )); // name of dynamic memory for CPU
	//----- ----- ----- -----
	//prepare microstructure
	// Initialize the concentration filed with random modulation
	micro_ch_pre_2d(nx,ny,c0,con);
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	float *con_d, *new_con_d; // name of dynamic memory for GPU, CUDA, device
	//
	con_d     = (float *)malloc(nx*ny*sizeof(float)); //GPU, CUDA, device
	new_con_d = (float *)malloc(nx*ny*sizeof(float)); //GPU, CUDA, device
	//
	cudaMalloc((void**)&con_d,    nx*ny*sizeof(float)); // define dynamic memory for GPU (device)
	cudaMalloc((void**)&new_con_d,nx*ny*sizeof(float)); // define dynamic memory for GPU (device)
	//
	//copy con_h(cpu,host) to con_d(cuda,device)
	cudaMemcpy(con_d,con,nx*ny*sizeof(float),cudaMemcpyHostToDevice); //con = con_h
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	int bs=BS; // Number of threads, 16 or 32
	dim3 blocks(nx/bs,ny/bs,1); //nx*ny = blocks * threads
	dim3 threads(bs,bs,1);      //bs*bs*1 <= 1024
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	
	//----- ----- ----- -----start:(This part is not really needed.)----- ----- ----- ----
	//Set recording time
	cudaEvent_t start, stop;
	
	//Initialization
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	//Start recording time
	cudaEventRecord(start);
	//----- ----- ----- -----end:(This part is not really needed.)----- ----- ----- ----
	
	//Time evolution of concentration filed
	for(int istep=0;istep<=nstep;istep++){
		
		
		//Update the total time
		ttime=ttime+dtime;
		
		//calculate subroutine "Kernel" on GPU
		Kernel<<<blocks, threads>>>(con_d,new_con_d,nx,ny,dx,dy,dtime,mobility,grad_coef);
		cudaDeviceSynchronize(); //<- new version | old version -> cudaThreadSynchronize();
		
		// replace con_d with new con_d
		update(&con_d,&new_con_d);
		
		if(istep%nprint == 0){
			//copy f_d(cuda,device) to F_h(cpu,host)
			cudaMemcpy(con,con_d,nx*ny*sizeof(float),cudaMemcpyDeviceToHost); //con = con_h
			
			/* Calculate the total bulk energy and print the result to
			   time_energy.out file for 2D plots */
			energy = calculate_energy_2d(nx,ny,con,grad_coef);
			
			//print the average free energy density value to file
			fprintf(out2, "%14.6e %14.6e \n",ttime, energy);
			
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_2D(nx,ny,dx,dy,istep,con);
			
			//show current step
			printf("done step: %5d \n",istep);
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
	
	cudaFree(con_d);
	cudaFree(new_con_d);
	
	free(con);
}
