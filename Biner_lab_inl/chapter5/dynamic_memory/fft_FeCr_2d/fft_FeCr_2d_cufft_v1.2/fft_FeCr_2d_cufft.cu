/* 2D semi-implicit spectral phase-field code 
  for solving precipitation in FeCr alloy */

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
//typedef cuComplex cufftComplex;

//----- ----- ----- ----- ----- ----- -----
void dislo_strain_2d(int Nx, int Ny, int idislo,
	float *ed11, float *ed22, float *ed12);
//----- ----- -----
float FeCr_chem_poten_2d(float cr_ij, float tempr);
//----- ----- -----
void green_tensor_2d(int Nx, int Ny,
	float *kx, float *ky,
	float cm11, float cm12, float cm44,
	float cp11, float cp12, float cp44,
	float *tmatx);
//----- ----- -----
void prepare_fft_2d(int Nx, int Ny, 
	float dx, float dy,
	float *kx, float *ky, 
	float *k2, float *k4);
//----- ----- -----
void solve_elasticity_2d(int Nx, int Ny, int BS,
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
void micro_ch_pre_2d(int Nx, int Ny, float c0, float *con);
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
	conk_d[j].x = ( conk_d[j].x - (dtime*mobility*k2_d[j]*(dfdconk_d[j].x + delsdck_d[j].x)) )/denom; //real part
	conk_d[j].y = ( conk_d[j].y - (dtime*mobility*k2_d[j]*(dfdconk_d[j].y + delsdck_d[j].y)) )/denom; //imaginary part
	
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
	int Nx=128;
	int Ny=128;
	
	//Number of threads, 2^n=<32, BS*BS*1 <= 1024
	int BS=32;
	
	//Total number of grid points in the simulation cell
	//int NxNy=Nx*Ny;
	
	//The distance between two grid points in x,y-direction
	float dx=1.0; // [nm] unit ?
	float dy=1.0; // [nm] unit ?
	
	//time integration parameters
	int nstep=10000; //Number of time steps
	int nprint=50;  //Print frequency to write the results to file
	float dtime=1.0e-2; //Time increment for numerical integration
	float ttime=0.0;    //Total time
	float coefA=2.0;
	
	//material specific parameters
	
	//Initial concentrations of alloying elements
	float c0=0.20;       //Initial concentraion (20%Cr-containing Fe-Cr alloy
	float mobility=1.0;  //The value of mobility coefficient (dimensionless)
	float grad_coef=0.5; //The value of gradient energy coefficients [J(nm)^2/mol]
	
	float tempr=535.0; //Annealing temperature [K]
	float RT=8.314462*tempr; //Gas constant x temperature
	
	//elastic constants
	//Elastic constants of Fe-rich phase [GPa]
	float cm11=233.10e3;
	float cm12=135.44e3;
	float cm44=178.30e3;
	//
	//Elastic constants of Cr-rich phase [GPa]
	float cp11=350.00e3;
	float cp12= 67.80e3;
	float cp44=100.80e3;
	
	//elastic constant of other materials
	//Ref: https://www.jstage.jst.go.jp/article/jsms/69/9/69_657/_pdf
	
	//eigen strains
	//The value of eigenstrains for Cr-rich phase
	float ei0=0.006; //Maginitude of eigenstrains
	
	int ii; //ii=(i*Ny+j);
	
	//----- ----- ----- ----- ----- -----
	//const int fftsizex = Nx, fftsizey = Ny;
	//
	cufftComplex *cr_d, *dfdcr_d, *delsdc_d;
	cudaMalloc((void**)&cr_d,      sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&dfdcr_d,   sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&delsdc_d,  sizeof(cufftComplex)*Nx*Ny);
	//
	cufftComplex *crk_d, *dfdcrk_d, *delsdck_d;
	cudaMalloc((void**)&crk_d,     sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&dfdcrk_d,  sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&delsdck_d, sizeof(cufftComplex)*Nx*Ny);
	//
	cufftHandle plan, iplan;
	//cufftPlan2d(&plan,  Nx, Ny, CUFFT_R2C);
	//cufftPlan2d(&iplan, Nx, Ny, CUFFT_C2R);
	cufftPlan2d(&plan,  Nx, Ny, CUFFT_C2C);
	cufftPlan2d(&iplan, Nx, Ny, CUFFT_C2C);
	//----- ----- ----- ----- ----- -----
	
	//----- ----- ----- -----fftw3
	//fftw_complex *s11, *s22, *s12;
	// s11 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	// s22 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	// s12 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	//----- ----- ----- -----
	//----- ----- ----- -----fftw3
	//fftw_complex *e11, *e22, *e12;
	// e11 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	// e22 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	// e12 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	//----- ----- ----- ----- ----- -----
	
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
		}
	}
	
	//Strain components due to lattice defects
	//float ed11[Nx][Ny];
	float *ed11 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	//float ed22[Nx][Ny];
	float *ed22 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	//float ed12[Nx][Ny];
	float *ed12 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	
	//dislocation eigen strain
	/* idislo=1 for dislocation diploe,
	   idislo=2 for dislocation array */
	int idislo=1;
	dislo_strain_2d(Nx,Ny,idislo,ed11,ed22,ed12);
	
	//Applied strains
	// The components of applied strains
	float ea[3]; //Magnitude of applied strains
	ea[0]=0.0;
	ea[1]=0.0;
	ea[2]=0.0;
	
	//int iflag=1;
	
	//----- prepare microstructure
	float *cr = (float *)malloc(sizeof(float)*( Nx*Ny ));
	//
	micro_ch_pre_2d(Nx,Ny,c0,cr); //Initialize microstructure
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
	
	//prepare fft (output: kx,ky,kz,k2,k4)
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
	
	//float *dfdcr = (float *)malloc(sizeof(float)*( Nx*Ny ));
	float *delsdc = (float *)malloc(sizeof(float)*( Nx*Ny ));
	//
	float _Complex *crc     = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *dfdcrc  = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *delsdcc = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	//
	//float _Complex *crk     = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	//float _Complex *dfdcrk  = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	//float _Complex *delsdck = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	
	//float numer, denom;
	//
	//float denom;
	//float _Complex numer;
	
	int bs=BS; // Number of threads, 16 or 32
	dim3 blocks(Nx/bs,Ny/bs,1); //nx*ny = blocks * threads
	dim3 threads(bs,bs,1);      //bs*bs*1 <= 1024
	
	//evolve (Time evolution of microstructure)
	for(int istep=0;istep<=nstep;istep++){
		
		//Update the total time
		ttime=ttime+dtime;
		
		//derivative of elastic energy
		//Calculate the derivative of elastic energy
		solve_elasticity_2d(Nx,Ny,BS,
			tmatx,
			s11,s22,s12,
			e11,e22,e12,
			ed11,ed22,ed12,
			cm11,cm12,cm44,
			cp11,cp12,cp44,
			ea,
			ei0,
			cr, delsdcc);
		// Note: tmatx is real part only. dslsc is output.
		//----- ----- ----- -----
		
		//derivative of chemical energy
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				//Calculate derivative of chemical energy
				dfdcrc[ii]=FeCr_chem_poten_2d(cr[ii],tempr);
				//----- ------ ------ ------
				//replace cuda array with host array
				delsdcc[ii] = delsdcc[ii]/RT; //And normalize the derivative elastic energy with RT
				crc[ii] = cr[ii];
			}
		}
		cudaMemcpy(dfdcr_d,dfdcrc,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //dfdcrc = dfdcrc_h
		cudaMemcpy(delsdc_d,delsdcc,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //delsdcc = delsdcc_h
		cudaMemcpy(cr_d,crc,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //cr = cr_h
		
		/* Take the values of concentration, derivative of free energy and
		   derivative of elastic energy from real space to Fourier space (forward FFT) */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, cr_d, crk_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, dfdcr_d, dfdcrk_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, delsdc_d, delsdck_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		//cudaMemcpy(crk,crk_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //crk = crk_h
		//cudaMemcpy(dfdcrk,dfdcrk_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //dfdcrk = dfdcrk_h
		//cudaMemcpy(delsdck,delsdck_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //delsdck = delsdck_h
		
		/* Semi-implicit time integration of Cr concentration field at
		   Fourier space (Eq.5.50) */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- on cuda
		Kernel_semi_implicit_time_integration<<<blocks, threads>>>(Nx,Ny,
			dtime,coefA,mobility,grad_coef,
			k2_d,k4_d,
			crk_d,dfdcrk_d,delsdck_d);
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		/* Take concentration field from Fourier space back to
		   real space (inverse FFT) */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(iplan, crk_d, cr_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		//copy f_d(cuda,device) to F_h(cpu,host)
		//cudaMemcpy(cr,cr_d,Nx*Ny*sizeof(float),cudaMemcpyDeviceToHost); //cr = cr_h
		cudaMemcpy(crc,cr_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //crc = crc_h
		
		//for small deviations
		// For small deviations from max and min values, reset the limits
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				//----- ----- ----- -----
				//cr[ii] =  cr[ii]/(Nx*Ny);
				cr[ii] = ( __real__ crc[ii] )/(Nx*Ny);
				//cr[ii] =  creal(crc[ii])/(Nx*Ny); //For #include <complex.h>
				//----- ----- ----- -----
				if(cr[ii]>=0.9999){
					cr[ii]=0.9999;
				}
				if(cr[ii]<=0.0001){
					cr[ii]=0.0001;
				}
				//----- ----- ----- -----
			}
		}
		
		//print results
		/* If print frequency is reached, output the results to file */
		if(fmod(istep,nprint)==0){
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_2D(Nx,Ny,dx,dy,istep,cr);
			
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
	cudaFree(cr_d);
	cudaFree(dfdcr_d);
	cudaFree(delsdc_d);
	//
	cudaFree(crk_d);
	cudaFree(dfdcrk_d);
	cudaFree(delsdck_d);
	//
	cudaFree(k2_d);
	cudaFree(k4_d);
	//----- ----- ----- ----- ----- -----
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
	free(cr);
	//
	free(crc);
	free(dfdcrc);
	free(delsdcc);
	//
	//free(crk);
	//free(dfdcrk);
	//free(delsdck);
	//----- ----- ----- ----- ----- -----
}
