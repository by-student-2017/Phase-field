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

void micro_ch_pre_2d(int Nx, int Ny, float c0, float *con);

void prepare_fft_2d(int Nx, int Ny, 
	float dx, float dy,
	float *kx, float *ky, 
	float *k2, float *k4);

float free_energy_ch_2d(float con_ij);

float calculate_energy_2d(int Nx, int Ny, float *con, float grad_coef);

void write_vtk_grid_values_2D(int nx, int ny, 
	float dx, float dy,
	int istep, float *data1);

int main(){
	clock_t start, end;
	float compute_time;
	
	//get initial wall time
	//(Get initial wall clock time beginning of the execution)
	start = clock();
	
	//simulation cell parameters
	int Nx=64;
	int Ny=64;
	//int NxNy=Nx*Ny; //Total number of grid points in the simulation cell
	
	//Number of threads, 2^n=<32, BS*BS*1 <= 1024
	int BS=32;
	
	//open an output file for writing total bulk energy values
	FILE *out2=fopen("time_energy.out","w");
	
	//The distance between two grid points in x,y-direction
	float dx=1.0; // [nm] unit ?
	float dy=1.0; // [nm] unit ?
	
	//time integration parameters
	int nstep=2000; //Number of time steps
	int nprint=50;  //Print frequency to write the results to file
	float dtime=1.0e-2; //Time increment for numerical integration
	float ttime=0.0;    //Total time
	float coefA=1.0;
	
	//material specific parameters
	
	//Initial concentrations of alloying elements
	float c0=0.40;       //Initial concentraion (20%Cr-containing Fe-Cr alloy
	float mobility=1.0;  //The value of mobility coefficient (dimensionless)
	float grad_coef=0.5; //The value of gradient energy coefficients [J(nm)^2/mol]
	
	float energy=0.0;
	
	int ij; //ij=(i*Ny+j);
	
	//----- ----- ----- ----- ----- -----
	//const int fftsizex = Nx, fftsizey = Ny;
	//
	cufftComplex *con_d, *dfdcon_d;
	cudaMalloc((void**)&con_d,     sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&dfdcon_d,  sizeof(cufftComplex)*Nx*Ny);
	//
	cufftComplex *conk_d, *dfdconk_d;
	cudaMalloc((void**)&conk_d,    sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&dfdconk_d, sizeof(cufftComplex)*Nx*Ny);
	//
	cufftHandle plan, iplan;
	cufftPlan2d(&plan,  Nx, Ny, CUFFT_C2C);
	cufftPlan2d(&iplan, Nx, Ny, CUFFT_C2C);
	//----- ----- ----- ----- ----- -----
	
	//----- ----- ----- ----- ----- -----
	//const int fftsizex = Nx, fftsizey = Ny;
	/* fftw_complex *in, *out; // in[i][0] for real, in[i][1] for imag.
	    in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	   out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	   fftw_plan plan, iplan;
	   plan = fftw_plan_dft_2d(fftsize, fftsize, in, out, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	  iplan = fftw_plan_dft_2d(fftsize, fftsize, in, out, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT */
	//
	//array
	//fftw_complex *con, *conk;
	//con  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	//conk = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	//fftw_plan plan_con, iplan_conk;
	// plan_con  = fftw_plan_dft_2d(fftsizex, fftsizey, con, conk, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//iplan_conk = fftw_plan_dft_2d(fftsizex, fftsizey, conk,con,  FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//array
	//fftw_complex *dfdcon, *dfdconk;
	// dfdcon = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	//dfdconk = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	//fftw_plan plan_dfdcon;
	// plan_dfdcon  = fftw_plan_dft_2d(fftsizex, fftsizey, dfdcon, dfdconk, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//----- ----- ----- ----- ----- -----
	
	//float con_out[Nx][Ny];
	float *con = (float *)malloc(sizeof(float)*( Nx*Ny ));
	
	//prepare microstructure
	// Initialize the concentration filed with random modulation
	micro_ch_pre_2d(Nx,Ny,c0,con);
	
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
	
	float _Complex *dfdconc = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *conc    = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	//
	float _Complex *dfdconk = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *conk    = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	
	float denom;
	float _Complex numer;
	
	int bs=BS; // Number of threads, 16 or 32
	dim3 blocks(Nx/bs,Ny/bs,1); //nx*ny = blocks * threads
	dim3 threads(bs,bs,1);      //bs*bs*1 <= 1024
	
	//evolve (Time evolution of microstructure)
	for(int istep=0;istep<=nstep;istep++){
		
		//Update the total time
		ttime=ttime+dtime;
		
		//derivative of free energy
		//Calculate the derivative of elastic energy
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ij=i*Ny+j;
				dfdconc[ij] = free_energy_ch_2d(con[ij]);
				//----- ------ ------ ------
				conc[ij] = con[ij];
			}
		}
		
		cudaMemcpy(dfdcon_d,dfdconc,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //dfdconc = dfdconc_h
		cudaMemcpy(con_d,conc,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //conc = conc_h
		
		/* Take the values of concentration, derivative of free energy and
		   derivative of elastic energy from real space to Fourier space (forward FFT) */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, con_d, conk_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, dfdcon_d, dfdconk_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		cudaMemcpy(conk,conk_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //conk = conk_h
		cudaMemcpy(dfdconk,dfdconk_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //dfdconk = dfdconk_h
		
		/* Semi-implicit time integration of concentration field at
		   Fourier space (Eq.5.14) */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- from cufft
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ij=i*Ny+j;
				//
				denom=1.0+dtime*coefA*mobility*grad_coef*k4[ij];
				//
				numer=dtime*mobility*k2[ij]*dfdconk[ij];
				conk[ij]=(conk[ij]-numer)/denom;
			}
		}
		cudaMemcpy(conk_d,conk,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //conk = conk_h
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
		
		/* Take concentration field from Fourier space back to
		   real space (inverse FFT) */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(iplan, conk_d, con_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		cudaMemcpy(conc,con_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //conc = conc_h
		
		//for small deviations
		// For small deviations from max and min values, reset the limits
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ij=i*Ny+j;
				//----- ----- ----- -----
				//con[ij] =  con[ij]/(Nx*Ny);
				con[ij] = ( __real__ conc[ij] )/(Nx*Ny);
				//con[ij] =  creal(conc[ij])/(Nx*Ny); //For #include <_Complex.h>
				//----- ----- ----- -----
				if(con[ij]>=0.9999){
					con[ij]=0.9999;
				}
				if(con[ij]<=0.0001){
					con[ij]=0.0001;
				}
				//----- ----- ----- -----
			}
		}
		
		//print results
		/* If print frequency is reached, output the results to file */
		if(fmod(istep,nprint)==0){
			/* Calculate the total bulk energy and print the result to
			   time_energy.out file for 2D plots */
			energy = calculate_energy_2d(Nx,Ny,con,grad_coef);
			
			//print the average free energy density value to file
			fprintf(out2, "%14.6e %14.6e \n",ttime, energy);
			
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_2D(Nx,Ny,dx,dy,istep,con);
			
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
	//
	free(conk);
	free(dfdconk);
	//----- ----- ----- ----- ----- -----
}
