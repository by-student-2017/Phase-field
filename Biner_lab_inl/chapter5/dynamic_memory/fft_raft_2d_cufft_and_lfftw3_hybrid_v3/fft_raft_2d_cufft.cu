/* 2D semi-implicit spectral phase-field code 
  for solving precipitation under stress */

/* This program solves Cahn-Hilliard phase-field equation with
   semi-implicit Fourier spectral method by 
   taking into account the effects of elastic inhomogeneities and 
   applied stresses based on solution of stress-strain fields with
   Green's tensor and Fourier transformations.
     The time integration is carried out by using semi-implicit
   time marching scheme. */

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
#include <cuComplex.h>
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
void micro_ch_pre_2d(int Nx, int Ny, float c0, float *con);
//----- ----- -----
void prepare_fft_2d(int Nx, int Ny, 
	float dx, float dy,
	float *kx, float *ky, 
	float *k2, float *k4);
//----- ----- -----
void green_tensor_2d(int Nx, int Ny,
	float *kx, float *ky,
	float cm11, float cm12, float cm44,
	float cp11, float cp12, float cp44,
	float *tmatx);
//----- ----- -----
float free_energy_ch_2d(float con_ij);
//----- ----- -----
void solve_elasticity_2d(int Nx, int Ny,
	float *tmatx,
	float _Complex *s11, float _Complex *s22, float _Complex *s12,
	float _Complex *e11, float _Complex *e22, float _Complex *e12,
	float *ed11, float *ed22, float *ed12,
	float cm11, float cm12, float cm44,
	float cp11, float cp12, float cp44,
	float *ea,
	float ei0,
	float *con, float _Complex *delsdc);
//----- ----- -----
void write_vtk_grid_values_2D(int nx, int ny, 
	float dx, float dy,
	int istep, float *data1);
//----- ----- ----- ----- ----- ----- -----

// Define subroutine "Kernel" for GPU (Device) calculation in detail
__global__ void Kernel_semi_implicit_time_integration(
	int   Nx,
	int   Ny,
	float dtime,
	float coefA,
	float mobility,
	float grad_coef,
	float *k2_d,
	float *k4_d,
	cufftComplex *conk_d,
	cufftComplex *dfdconk_d,
	cufftComplex *delsdck_d
){
	int j, jx, jy;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	jx = blockDim.x*blockIdx.x + threadIdx.x; //<-GPU | CPU -> for(jx=0; jx<nx; jx++){
	jy = blockDim.y*blockIdx.y + threadIdx.y; //<-GPU | CPU -> for(jy=0; jy<ny; jy++){
	j  = Nx*jy + jx;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	float denom;
	//
	denom = 1.0 + dtime*coefA*mobility*grad_coef*k4_d[j];
	conk_d[j].x = ( conk_d[j].x - (dtime*mobility*k2_d[j]*(dfdconk_d[j].x + delsdck_d[j].x)) )/denom;
	conk_d[j].y = ( conk_d[j].y - (dtime*mobility*k2_d[j]*(dfdconk_d[j].y + delsdck_d[j].y)) )/denom;
}

int main(){
	clock_t start, end;
	float compute_time;
	
	//get initial wall time
	//(Get initial wall clock time beginning of the execution)
	start = clock();
	
	//simulation cell parameters
	int Nx=256;
	int Ny=256;
	
	//Number of threads, 2^n=<32, BS*BS*1 <= 1024
	int BS=32;
	
	//Total number of grid points in the simulation cell
	//int NxNy=Nx*Ny;
	
	//The distance between two grid points in x,y-direction
	float dx=1.0; // [nm] unit ?
	float dy=1.0; // [nm] unit ?
	
	//time integration parameters
	int nstep=5000; //Number of time steps
	int nprint=25;  //Print frequency to write the results to file
	float dtime=5.0e-2; //Time increment for numerical integration
	float ttime=0.0;    //Total time
	float coefA=1.0;
	
	//material specific parameters
	
	//Initial concentrations of alloying elements
	float c0=0.40;       //Initial concentraion
	float mobility=1.0;  //The value of mobility coefficient (dimensionless)
	float grad_coef=0.5; //The value of gradient energy coefficients [J(nm)^2/mol]
	
	//elastic constants
	//Elastic constants of matrix phase
	float cm11=1400.0;
	float cm12= 600.0;
	float cm44= 400.0;
	//
	//Elastic constants of second phase
	float cp11=2.0*cm11;
	float cp12=2.0*cm12;
	float cp44=2.0*cm44;
	
	//eigen strains
	float ei0=0.01; //Maginitude of eigenstrains
	
	//Applied strains
	float ea[3]; //Magnitude of applied strains
	ea[0]=0.00;
	ea[1]=0.01;
	ea[2]=0.00;
	
	int ii;
	
	//----- ----- ----- ----- ----- -----
	//const int fftsizex = Nx, fftsizey = Ny;
	//
	cufftComplex *con_d, *dfdcon_d, *delsdc_d;
	cudaMalloc((void**)&con_d,     sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&dfdcon_d,  sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&delsdc_d,  sizeof(cufftComplex)*Nx*Ny);
	//
	cufftComplex *conk_d, *dfdconk_d, *delsdck_d;
	cudaMalloc((void**)&conk_d,    sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&dfdconk_d, sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&delsdck_d, sizeof(cufftComplex)*Nx*Ny);
	//
	cufftHandle plan, iplan;
	cufftPlan2d(&plan,  Nx, Ny, CUFFT_C2C);
	cufftPlan2d(&iplan, Nx, Ny, CUFFT_C2C);
	//----- ----- ----- ----- ----- -----
	
	//----- ----- ----- -----
	//fftw_Complex *s11, *s22, *s12;
	// s11 = (fftw_Complex *)fftw_malloc(sizeof(fftw_Complex) * fftsizex*fftsizey);
	// s22 = (fftw_Complex *)fftw_malloc(sizeof(fftw_Complex) * fftsizex*fftsizey);
	// s12 = (fftw_Complex *)fftw_malloc(sizeof(fftw_Complex) * fftsizex*fftsizey);
	//----- ----- ----- -----
	//fftw_Complex *e11, *e22, *e12;
	// e11 = (fftw_Complex *)fftw_malloc(sizeof(fftw_Complex) * fftsizex*fftsizey);
	// e22 = (fftw_Complex *)fftw_malloc(sizeof(fftw_Complex) * fftsizex*fftsizey);
	// e12 = (fftw_Complex *)fftw_malloc(sizeof(fftw_Complex) * fftsizex*fftsizey);
	//----- ----- ----- ----- ----- -----
	
	//float ed11[Nx][Ny];
	float *ed11 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	//float ed22[Nx][Ny];
	float *ed22 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	//float ed12[Nx][Ny];
	float *ed12 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	
	float _Complex *s11  = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *s22  = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *s12  = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	//
	float _Complex *e11  = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *e22  = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *e12  = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	
	//Initialize stress & strain componentes
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			ii=i*Ny+j;
			//----- ----- ----- -----
			s11[ii] = 0.0;
			s22[ii] = 0.0;
			s12[ii] = 0.0;
			//----- ----- ----- -----
			e11[ii] = 0.0;
			e22[ii] = 0.0;
			e12[ii] = 0.0;
			//----- ----- ----- -----
			//----- ----- ----- -----
			//Strain components due to lattice defects
			ed11[ii] = 0.0;
			ed22[ii] = 0.0;
			ed12[ii] = 0.0;
			//----- ----- ----- -----
		}
	}
	
	//----- prepare microstructure
	float *con = (float *)malloc(sizeof(float)*( Nx*Ny ));
	//
	micro_ch_pre_2d(Nx,Ny,c0,con); //Initialize microstructure
	//----- ----- ----- -----
	
	//----- prepare fft (output: kx,ky,kz,k2,k4)
	//float kx[Nx];
	float *kx = (float *)malloc(sizeof(float)*( Nx ));
	//float ky[Ny];
	float *ky = (float *)malloc(sizeof(float)*( Ny ));
	//float k2[Nx][Ny];
	float *k2 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	//float k4[Nx][Ny];
	float *k4 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	//
	prepare_fft_2d(Nx,Ny,dx,dy,kx,ky,k2,k4); //Calculate coefficients of Fourier transformation
	//----- ----- ----- -----
	
	//----- for cufft
	float *k2_d, *k4_d;
	k2_d  = (float *)malloc(Nx*Ny*sizeof(float));
	k4_d  = (float *)malloc(Nx*Ny*sizeof(float));
	cudaMalloc((void**)&k2_d ,Nx*Ny*sizeof(float));
	cudaMalloc((void**)&k4_d ,Nx*Ny*sizeof(float));
	cudaMemcpy(k2_d,k2,Nx*Ny*sizeof(float),cudaMemcpyHostToDevice); //k2 = k2_h
	cudaMemcpy(k4_d,k4,Nx*Ny*sizeof(float),cudaMemcpyHostToDevice); //k4 = k4_h
	//----- ----- ----- -----
	
	//float tmatx[Nx][Ny][2][2][2][2];
	float *tmatx = (float *)malloc(sizeof(float)*( Nx*Ny*2*2*2*2 )); //real part only
	
	//Greens tensor
	green_tensor_2d(Nx,Ny,kx,ky,cm11,cm12,cm44,cp11,cp12,cp44,tmatx); //Calculate Green's tensor
	
	float _Complex *delsdcc = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *dfdconc = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *conc    = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	
	int bs=BS; // Number of threads, 16 or 32
	dim3 blocks(Nx/bs,Ny/bs,1); //nx*ny = blocks * threads
	dim3 threads(bs,bs,1);      //bs*bs*1 <= 1024
	
	//evolve (Time evolution of microstructure)
	for(int istep=0;istep<=nstep;istep++){
		
		//Update the total time
		ttime=ttime+dtime;
		
		//derivative of elastic energy
		//Calculate the derivative of elastic energy
		solve_elasticity_2d(Nx,Ny,
			tmatx,
			s11,s22,s12,
			e11,e22,e12,
			ed11,ed22,ed12,
			cm11,cm12,cm44,
			cp11,cp12,cp44,
			ea,
			ei0,
			con, delsdcc); // Note: tmatx is real part only
		//----- ----- ----- -----
		
		//derivative of free energy and replacement
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				//Calculate derivative of free energy
				//dfdcon[ii] = free_energy_ch_2d(con[ii]);
				dfdconc[ii] = free_energy_ch_2d(con[ii]);
				//----- ------ ------ ------
				conc[ii] = con[ii];
			}
		}
		cudaMemcpy(dfdcon_d,dfdconc,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //dfdconc = dfdconc_h
		cudaMemcpy(delsdc_d,delsdcc,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //delsdc = delsdc_h
		cudaMemcpy(con_d,conc,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //con = con_h
		
		/* Take the values of concentration, derivative of free energy and
		   derivative of elastic energy from real space to Fourier space (forward FFT) */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, con_d, conk_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, dfdcon_d, dfdconk_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, delsdc_d, delsdck_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		/* Semi-implicit time integration of concentration field at
		   Fourier space (Eq.5.50) */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- on cuda
		Kernel_semi_implicit_time_integration<<<blocks, threads>>>(Nx,Ny,
			dtime,coefA,mobility,grad_coef,
			k2_d,k4_d,
			conk_d,dfdconk_d,delsdck_d);
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
		
		/* Take concentration field from Fourier space back to
		   real space (inverse FFT) */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(iplan, conk_d, con_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		//copy f_d(cuda,device) to F_h(cpu,host)
		//cudaMemcpy(con,con_d,Nx*Ny*sizeof(float),cudaMemcpyDeviceToHost); //con = con_h
		cudaMemcpy(conc,con_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //conc = conc_h
		
		//for small deviations
		// For small deviations from max and min values, reset the limits
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				//----- ----- ----- -----
				//con[ii] =  con[ii]/(Nx*Ny);
				con[ii] = ( __real__ conc[ii] )/(Nx*Ny);
				//con[ii] =  creal(conc[ii])/(Nx*Ny); //For #include <_Complex.h>
				//----- ----- ----- -----
				if(con[ii]>=0.9999){
					con[ii]=0.9999;
				}
				if(con[ii]<=0.0001){
					con[ii]=0.0001;
				}
				//----- ----- ----- -----
			}
		}
		
		//print results
		/* If print frequency is reached, output the results to file */
		if(fmod(istep,nprint)==0){
			printf("done step: %5d \n",istep);
			
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_2D(Nx,Ny,dx,dy,istep,con);
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
	cudaFree(delsdc_d);
	//
	cudaFree(conk_d);
	cudaFree(dfdconk_d);
	cudaFree(delsdck_d);
	//
	free(s11);
	free(s22);
	free(s12);
	//
	free(e11);
	free(e22);
	free(e12);
	//----- ----- ----- ----- ----- -----
	free(ed11);
	free(ed22);
	free(ed12);
	//
	free(kx);
	free(ky);
	free(k2);
	free(k4);
	//
	free(tmatx);
	//
	free(con);
	//
	free(conc);
	free(dfdconc);
	free(delsdcc);
	//----- ----- ----- ----- ----- -----
}
