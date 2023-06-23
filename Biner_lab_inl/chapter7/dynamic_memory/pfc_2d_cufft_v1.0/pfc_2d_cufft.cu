/* 2D semi-implicit spectral
  phase-field crystal code */

#include <stdio.h>
#include <stdlib.h> //rand()
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

void prepare_fft_2d(int Nx, int Ny, 
	float dx, float dy,
	float *kx, float *ky, 
	float *k2, float *k4);

void write_vtk_grid_values_2D(int nx, int ny, 
	float dx, float dy,
	int istep, float *data1);

int main(){
	clock_t start, end;
	float compute_time;
	
	//get initial wall time
	start = clock();
	
	FILE *out1=fopen("final_conf.out","w");
	FILE *out2=fopen("energy.out","w");
	
	//simulation cell parameters
	int Nx=64;
	int Ny=64;
	//int NxNy=Nx*Ny; //Total number of grid points in the simulation cell
	
	//Number of threads, 2^n=<32, BS*BS*1 <= 1024
	int BS=32;
	
	//The value of pi
	float pix=4.0*atan(1.0);
	
	//The distance between two grid points in x,y,z-direction
	float dx=pix/4.0;
	float dy=pix/4.0;
	
	//time integration parameters
	int nstep=200000;
	int nprint=2000;
	float dtime=0.05;
	
	//material specific parameters
	//float den0=-0.085; //average density for pfc_3d_v1 (stripe phase)
	float den0=-0.285; //average density for pfc_3d_v2 (triangular phase)
	float tempr=-0.25; //temperature (T-Tm)/Tm, Tm=melting point
	//float tempr0=0.0; //positive value, tempr=tempr+tempr0/isteps;
	float noise=den0*1e-2; //Noise term to modulate the initial density field
	
	int ii;
	
	//----- ----- ----- ----- ----- -----
	//const int fftsizex = Nx, fftsizey = Ny;
	//
	cufftComplex *den_d, *f_den_d;
	cudaMalloc((void**)&den_d,    sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&f_den_d,  sizeof(cufftComplex)*Nx*Ny);
	//
	cufftComplex *den3_d, *f_den3_d;
	cudaMalloc((void**)&den3_d,   sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&f_den3_d, sizeof(cufftComplex)*Nx*Ny);
	//
	cufftComplex *ff_d, *f_ff_d;
	cudaMalloc((void**)&ff_d,     sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&f_ff_d,   sizeof(cufftComplex)*Nx*Ny);
	//
	cufftHandle plan, iplan;
	cufftPlan2d(&plan,  Nx, Ny, CUFFT_C2C);
	cufftPlan2d(&iplan, Nx, Ny, CUFFT_C2C);
	//----- ----- ----- ----- ----- -----
	
	//float kx[Nx];
	float *kx = (float *)malloc(sizeof(float)*( Nx ));
	//float ky[Ny];
	float *ky = (float *)malloc(sizeof(float)*( Ny ));
	//float k2[Nx][Ny];
	float *k2 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	//float k4[Nx][Ny];
	float *k4 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	
	//prepare fft (output: kx,ky,kz,k2,k4)
	prepare_fft_2d(Nx,Ny,dx,dy,kx,ky,k2,k4); //get FFT coefficients
	
	float _Complex *denc   = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *den3c  = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *Nonxc  = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *f_den  = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *f_den3 = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *f_ff   = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *ffc    = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float *den = (float *)malloc(sizeof(float)*( Nx*Ny ));
	float *ff  = (float *)malloc(sizeof(float)*( Nx*Ny ));
	
	//if infile=1 read input from file
	int infile = 0;
	FILE *in1;
	int mx,my;
	if(infile==1){
		//open input file
		in1=fopen("g3_2r.inp","r");
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				fscanf(in1,"%5d %5d %f",&mx,&my,&den[ii]);
				denc[ii] = den[ii];
			}
		}
	}else{
		//initialize density
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				den[ii] = den0 + noise*(0.5-(float)rand()/RAND_MAX);
				denc[ii] = den[ii];
			}
		}
	}
	
	//float Linx[Nx][Ny];
 	float *Linx  = (float *)malloc(sizeof(float)*( Nx*Ny ));
 	//float denom[Nx][Ny];
 	float *denom = (float *)malloc(sizeof(float)*( Nx*Ny ));
 	//
 	float energy;
	//float ss2[Nx][Ny];
	float *ss2 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	//float ss4[Nx][Ny];
	float *ss4 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	
	int bs=BS; // Number of threads, 16 or 32
	dim3 blocks(Nx/bs,Ny/bs,1); //nx*ny = blocks * threads
	dim3 threads(bs,bs,1);      //bs*bs*1 <= 1024
	
	//evolve (evolve microstructure)
	for(int istep=0;istep<=nstep;istep++){
		
		//tempr = tempr + tempr0/istep;
		
		cudaMemcpy(den_d,denc,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //denc = denc_h
		
		//take current density filed from real space to Fourier space (forward FFT transformation)
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, den_d, f_den_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cudaMemcpy(f_den, f_den_d, Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //f_den = f_den_h
		
		//calculate the value of denominator in Eq.7.10 at every grid points
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				Linx[ii]=-k2[ii]*(tempr+1.0-2.0*k2[ii]+k4[ii]);
				denom[ii]=1.0-dtime*Linx[ii];
			}
		}
		
		//calculate the nonlinear term, phi^3, in Eq.7.10
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				den3c[ii]=denc[ii]*denc[ii]*denc[ii];
			}
		}
		cudaMemcpy(den3_d,den3c,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //den3c = den3c_h
		
		//take the value of phi^3 from real space to Fourier space (forward FFT transformation)
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, den3_d, f_den3_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cudaMemcpy(f_den3,f_den3_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //f_den3 = f_den3_h
		
		//calculate the value of phi^(t+1) from Fourier space to real space (inverse FFT transformation)
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				Nonxc[ii]=-k2[ii]*f_den3[ii];
				f_den[ii]=(f_den[ii]+dtime*Nonxc[ii])/denom[ii];
			}
		}
		cudaMemcpy(f_den_d,f_den,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //f_den = f_den_h
		
		//bring back the values of phi^(t+1) from Fourier space to real space (inverse FFT transformation)
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(iplan, f_den_d, den_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cudaMemcpy(denc,den_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //denc = denc_h
		
		//print results
		//if print frequency is reached, output the results to file
		if(fmod(istep,nprint)==0){
		
			printf("done step: %5d \n", istep);
			
			//energy calculation
			//calculate the free energy distribution, Eq.7.6
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					ii=i*Ny+j;
					den[ii] = ( __real__ denc[ii] )/(Nx*Ny);
					ss2[ii]=den[ii]*den[ii];
					ss4[ii]=ss2[ii]*ss2[ii];
				}
			}
			
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					ii=i*Ny+j;
					f_ff[ii]=0.5*f_den[ii]*(1.0-2.0*k2[ii]+k4[ii]);
				}
			}
			cudaMemcpy(f_ff_d,f_ff,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //f_ff = f_ff_h
			
			//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
			cufftExecC2C(iplan, f_ff_d, ff_d, CUFFT_INVERSE); //IFFT
			cudaDeviceSynchronize();
			//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
			cudaMemcpy(ffc,ff_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //ffc = ffc_h
			
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					ii=i*Ny+j;
					ff[ii] = (( __real__ ffc[ii] )/(Nx*Ny)) * den[ii]
					+ 0.5*tempr*ss2[ii]
					+ 0.25*ss4[ii];
				}
			}
			
			//integrate the free energy field
			energy = 0.0;
			
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					ii=i*Ny+j;
					energy = energy + ff[ii];
				}
			}
			
			//average free energy density
			energy = energy/(Nx*Ny);
			
			//print the average free energy density value to file
			fprintf(out2, "%d %14.6e \n",istep, energy);
			
			//output the results in vtk file format for contour plots to be viewed by using paraview
			write_vtk_grid_values_2D(Nx,Ny,dx,dy,istep,den);
			
		}
		
		//if intermediate configuration files are required, print the density field to file
		if(istep==nstep){
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					ii=i*Ny+j;
					fprintf(out1,"%5d %5d %14.6e \n",i,j,den[ii]);
				}
			}
		}
		
		//for recycle
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				denc[ii] = denc[ii]/(Nx*Ny);
			}
		}
		
	}//end of time step (evolve,for)
	
	//calculate the execution time and print it
	end = clock();
	compute_time = ((float) (end - start)) / CLOCKS_PER_SEC;
	printf("Compute Time: %f \n", compute_time);
	
	//----- ----- ----- ----- ----- -----
	cufftDestroy(plan);
	cufftDestroy(iplan);
	//----- ----- ----- ----- ----- -----
	cudaFree(den_d);
	cudaFree(f_den_d);
	//
	cudaFree(den3_d);
	cudaFree(f_den3_d);
	//
	cudaFree(f_ff_d);
	cudaFree(ff_d);
	//----- ----- ----- ----- ----- -----
	free(kx);
	free(ky);
	free(k2);
	free(k4);
	//
	free(denc);
	free(den3c);
	free(Nonxc);
	free(f_den);
	free(f_den3);
	free(f_ff);
	free(ffc);
	free(den);
	free(ff);
	//----- ----- ----- ----- ----- -----
	fclose(out1);
	fclose(out2);
	//----- ----- ----- ----- ----- -----
}
