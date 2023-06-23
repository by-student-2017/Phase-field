/* 2D semi-implicit spectral phase-field code
   for solving Cahn-Hilliard equation */

/* This program solves conserved phase-field equation with
   Fourier spactral method by taking into account the effects of
   elastic inhomogeneities and lattice defects, based on
   solution of stress-strain fields ith Green's tensor and
   Fourier transformations. the time integration is
   carried out by using semi-implicit time machining scheme. */

#include <stdio.h>
#include <stdlib.h> //rand() and malloc
#include <math.h> //mod() and -lm
#include <time.h>

//#include <fftw3.h>
//gcc test.c -lfftw3
//#include <mpi.h> //mpi version
//#include <fftw3-mpi.h> //mpi version

/* Memo for complex type
  "float __complex__ " is old version
  "float _Complex " is new version */
//#include <complex.h>
//#include <cuComplex.h>
//#define _Complex_I (1.0iF)
//#define I i
//#undef i
//#undef j
//Ref: http://nalab.mind.meiji.ac.jp/~mk/labo/text/complex-c.pdf (Japanese)

#include <cuda.h>   //GPU
/* #include <cuda.h> or
  #include "cuda_runtime.h"
  #include "device_launch_parameters.h" */
//----- ----- -----
#include <cufft.h> //FFT (GPU)

//typedef float cufftReal;
//typedef cu_Complex cufftComplex;

//----- ----- ----- ----- ----- ----- -----
void micro_ch_pre_3d(int Nx, int Ny, int Nz, float c0, float *con);
//----- ----- -----
void prepare_fft_3d(int Nx, int Ny, int Nz, 
	float dx, float dy, float dz,
	float *kx, float *ky, float *kz, 
	float *k2, float *k4);
//----- ----- -----
float free_energy_ch_3d(float con_ijk);
//----- ----- -----
float calculate_energy_3d(int Nx, int Ny, int Nz, float *con, float grad_coef);
//----- ----- -----
void write_vtk_grid_values_3D(int nx, int ny, int nz, 
	float dx, float dy, float dz,
	int istep, float *data1);
//----- ----- ----- ----- ----- ----- -----

// Define subroutine "Kernel" for GPU (Device) calculation in detail
__global__ void Kernel_semi_implicit_time_integration(
	int   Nx,
	int   Ny,
	int   Nz,
	float dtime,
	float coefA,
	float mobility,
	float grad_coef,
	float *k2_d,
	float *k4_d,
	cufftComplex *conk_d,
	cufftComplex *dfdconk_d
){
	int j, jx, jy, jz;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	jx = blockDim.x*blockIdx.x + threadIdx.x; //<-GPU | CPU -> for(jx=0; jx<nx; jx++){
	jy = blockDim.y*blockIdx.y + threadIdx.y; //<-GPU | CPU -> for(jy=0; jy<ny; jy++){
	jz = blockDim.z*blockIdx.z + threadIdx.z; //<-GPU | CPU -> for(jz=0; jz<nz; jz++){
	j  = (jz*Ny + jy)*Nx + jx; //j = nx*ny*jz + nx*jy + jx;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	float denom;
	//
	denom = 1.0 + dtime*coefA*mobility*grad_coef*k4_d[j];
	conk_d[j].x = ( conk_d[j].x - (dtime*mobility*k2_d[j]*dfdconk_d[j].x) )/denom; //real part
	conk_d[j].y = ( conk_d[j].y - (dtime*mobility*k2_d[j]*dfdconk_d[j].y) )/denom; //imaginary part
	
	/* Note: cufftComplex changed between CUDA 1.0 and 1.1.
	dout[idx].x =  d_signal[idx].y; <- dout[idx][0] = d_signal[idx][1];
	dout[idx].y = -d_signal[idx].x; <- dout[idx][1] = d_signal[idx][0]*(-1.0);
	Ref: https://forums.developer.nvidia.com/t/using-cufftcomplex-type-inside-a-kernel-does-it-work/4039 */
}

int main(){
	clock_t start, end;
	float compute_time;
	
	//get initial wall time
	//(Get initial wall clock time beginning of the execution)
	start = clock();
	
	//simulation cell parameters
	int Nx=64; //Number of grid points in the x-direction
	int Ny=64; //Number of grid points in the y-direction
	int Nz=2; //Number of grid points in the y-direction
	int NxNyNz=Nx*Ny*Nz; //Total number of grid points in the simulation cell
	
	int BSX=8; //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
	int BSY=8; //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
	int BSZ=2; //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
	
	//open an output file for writing total bulk energy values
	FILE *out2=fopen("time_energy.out","w");
	
	//The distance between two grid points in x,y-direction
	float dx=1.0; //Grid spacing between two grid pints in x-direction
	float dy=1.0; //Grid spacing between two grid pints in y-direction
	float dz=1.0; //Grid spacing between two grid pints in z-direction
	
	//time integration parameters
	int nstep=10000; //Number of time integration steps
	int nprint=50;  //Output frequency to write the results to file
	float dtime=1.0e-2; //Time increment for numerical integration
	float ttime=0.0;    //Total time
	float coefA=1.0;
	
	//material specific parameters
	
	//Initial concentrations of alloying elements
	float c0=0.40;       //Initial concentraion (20%Cr-containing Fe-Cr alloy
	float mobility=1.0;  //The value of mobility coefficient (dimensionless)
	float grad_coef=0.5; //The value of gradient energy coefficients [J(nm)^2/mol]
	
	float energy=0.0;
	
	int ijk; //ijk=(i*Ny+j)*Nz+k;
	
	//----- ----- ----- ----- ----- -----
	//const int fftsizex = Nx, fftsizey = Ny;
	//
	cufftComplex *con_d, *dfdcon_d;
	cudaMalloc((void**)&con_d,     sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&dfdcon_d,  sizeof(cufftComplex)*NxNyNz);
	//
	cufftComplex *conk_d, *dfdconk_d;
	cudaMalloc((void**)&conk_d,    sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&dfdconk_d, sizeof(cufftComplex)*NxNyNz);
	//
	cufftHandle plan, iplan;
	cufftPlan3d(&plan,  Nx, Ny, Nz, CUFFT_C2C);
	cufftPlan3d(&iplan, Nx, Ny, Nz, CUFFT_C2C);
	//----- ----- ----- ----- ----- -----
	
	//float con_out[Nx][Ny];
	float *con = (float *)malloc(sizeof(float)*( NxNyNz ));
	
	//prepare microstructure
	// Initialize the concentration filed with random modulation
	micro_ch_pre_3d(Nx,Ny,Nz,c0,con);
	
	//----- ----- ----- ----- ----- ----- ----- ----- -----
	//float kx[Nx];
	float *kx = (float *)malloc(sizeof(float)*( Nx ));
	//float ky[Ny];
	float *ky = (float *)malloc(sizeof(float)*( Ny ));
	//float kz[Ny];
	float *kz = (float *)malloc(sizeof(float)*( Nz ));
	//float k2[Nx][Ny][Nz];
	float *k2 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//float k4[Nx][Ny][Nz];
	float *k4 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//----- ----- ----- ----- ----- ----- ----- ----- -----
	
	//prepare fft (output: kx,ky,kz,k2,k4)
	prepare_fft_3d(Nx,Ny,Nz,dx,dy,dz,kx,ky,kz,k2,k4); //get FFT coefficients
	
	//----- for cufft
	float *k2_d, *k4_d;
	k2_d  = (float *)malloc(NxNyNz*sizeof(float));
	k4_d  = (float *)malloc(NxNyNz*sizeof(float));
	cudaMalloc((void**)&k2_d ,NxNyNz*sizeof(float));
	cudaMalloc((void**)&k4_d ,NxNyNz*sizeof(float));
	cudaMemcpy(k2_d,k2,NxNyNz*sizeof(float),cudaMemcpyHostToDevice); //k2 = k2_h
	cudaMemcpy(k4_d,k4,NxNyNz*sizeof(float),cudaMemcpyHostToDevice); //k4 = k4_h
	//----- ----- ----- -----
	
	float _Complex *dfdconc = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *conc    = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	
	int bsx=BSX, bsy=BSY, bsz=BSZ;     //Number of threads
	dim3 blocks(Nx/bsx,Ny/bsy,Nz/bsz); //nx*ny*nz = blocks * threads
	dim3 threads(bsx,bsy,bsz);         //bsx*bsy*bsz <= 1024
	
	//evolve (Time evolution of microstructure)
	for(int istep=0;istep<=nstep;istep++){
		
		//Update the total time
		ttime=ttime+dtime;
		
		//derivative of free energy
		//Calculate the derivative of elastic energy
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					dfdconc[ijk] = free_energy_ch_3d(con[ijk]);
					//----- ------ ------ ------
					conc[ijk] = con[ijk];
				}
			}
		}
		
		cudaMemcpy(dfdcon_d,dfdconc,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //dfdconc = dfdconc_h
		cudaMemcpy(con_d,conc,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //conc = conc_h
		
		/* Take the values of concentration, derivative of free energy and
		   derivative of elastic energy from real space to Fourier space (forward FFT) */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, con_d, conk_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, dfdcon_d, dfdconk_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		/* Semi-implicit time integration of concentration field at
		   Fourier space (Eq.5.14) */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----  on cuda
		Kernel_semi_implicit_time_integration<<<blocks, threads>>>(Nx,Ny,Nz,
			dtime,coefA,mobility,grad_coef,
			k2_d,k4_d,
			conk_d,dfdconk_d);
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
		
		/* Take concentration field from Fourier space back to
		   real space (inverse FFT) */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(iplan, conk_d, con_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		cudaMemcpy(conc,con_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //conc = conc_h
		
		//for small deviations
		// For small deviations from max and min values, reset the limits
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					//----- ----- ----- -----
					//con[ijk] =  con[ijk]/(NxNyNz);
					con[ijk] = ( __real__ conc[ijk] )/(NxNyNz);
					//con[ijk] =  creal(conc[ijk])/(NxNyNz); //For #include <_Complex.h>
					//----- ----- ----- -----
					if(con[ijk]>=0.9999){
						con[ijk]=0.9999;
					}
					if(con[ijk]<=0.0001){
						con[ijk]=0.0001;
					}
					//----- ----- ----- -----
				}
			}
		}
		
		//print results
		/* If print frequency is reached, output the results to file */
		if(fmod(istep,nprint)==0){
			/* Calculate the total bulk energy and print the result to
			   time_energy.out file for 2D plots */
			energy = calculate_energy_3d(Nx,Ny,Nz,con,grad_coef);
			
			//print the average free energy density value to file
			fprintf(out2, "%14.6e %14.6e \n",ttime, energy);
			
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_3D(Nx,Ny,Nz,dx,dy,dz,istep,con);
			
			printf("done step: %5d, time: %f \n",istep,ttime*(mobility/(dx*dx)));
			/* The quantities having the dimension of distance were normalized with the
			   magnitude of the Burger's vector, the  quantities having the dimension of
			   energy were normalized with RT, and the time t was normalized with M/(dx^2).
			   The initial concentration was modulated by setting the noise term to
			   0.001 in function micro_ch_pre_2d.c */
		}
		
	}//end of time step (evolve,for)
	
	//calculate the execution time and print it
	/* Calculate the compute time and print it to screen */
	end = clock();
	compute_time = ((float) (end - start)) / CLOCKS_PER_SEC;
	printf("Compute Time: %lf \n", compute_time);
	
	//----- ----- ----- ----- ----- -----
	cufftDestroy(plan);
	cufftDestroy(iplan);
	//----- ----- ----- ----- ----- -----
	cudaFree(con_d);
	cudaFree(dfdcon_d);
	//
	cudaFree(conk_d);
	cudaFree(dfdconk_d);
	//----- ----- ----- ----- ----- -----
	free(kx);
	free(ky);
	free(k2);
	free(k4);
	//
	free(con);
	//
	free(conc);
	free(dfdconc);
	//----- ----- ----- ----- ----- -----
}
