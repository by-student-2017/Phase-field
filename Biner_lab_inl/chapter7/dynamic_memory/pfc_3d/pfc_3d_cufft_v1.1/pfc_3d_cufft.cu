/* 3D semi-implicit spectral
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

//----- ----- ----- ----- ----- ----- -----
void prepare_fft_3d(int Nx, int Ny, int Nz, 
	float dx, float dy, float dz,
	float *kx, float *ky, float *kz, 
	float *k2, float *k4);
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
	float tempr,
	float *k2_d,
	float *k4_d,
	cufftComplex *f_den_d,
	cufftComplex *f_den3_d
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
	//calculate the value of denominator in Eq.7.10 at every grid points
	denom = 1.0 + dtime*k2_d[j]*(tempr+1.0-2.0*k2_d[j]+k4_d[j]);
	//----- ----- -----
	//calculate the value of phi^(t+1) from Fourier space to real space (inverse FFT transformation)
	f_den_d[j].x = ( f_den_d[j].x - dtime*k2_d[j]*f_den3_d[j].x )/denom; //real part
	f_den_d[j].y = ( f_den_d[j].y - dtime*k2_d[j]*f_den3_d[j].y )/denom; //imaginary part
	
	/* Note: cufftComplex changed between CUDA 1.0 and 1.1.
	dout[idx].x =  d_signal[idx].y; <- dout[idx][0] = d_signal[idx][1];
	dout[idx].y = -d_signal[idx].x; <- dout[idx][1] = d_signal[idx][0]*(-1.0);
	Ref: https://forums.developer.nvidia.com/t/using-cufftcomplex-type-inside-a-kernel-does-it-work/4039 */
}

int main(){
	clock_t start, end;
	float compute_time;
	
	//get initial wall time
	start = clock();
	
	FILE *out1=fopen("final_conf.out","w");
	FILE *out2=fopen("energy.out","w");
	
	//simulation cell parameters
	int Nx=64; //Number of grid points in the x-direction
	int Ny=64; //Number of grid points in the y-direction
	int Nz=2; //Number of grid points in the y-direction
	int NxNyNz=Nx*Ny*Nz; //Total number of grid points in the simulation cell
	
	int BSX=8; //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
	int BSY=8; //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
	int BSZ=2; //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
	
	//The value of pi
	float pix=4.0*atan(1.0);
	
	//The distance between two grid points in x,y,z-direction
	float dx=pix/4.0;
	float dy=pix/4.0;
	float dz=pix/4.0;
	
	//time integration parameters
	//int nstep=20000; //for pfc_3d_v1
	int nstep=200000; //for pfc_3d_v2
	int nprint=5000;
	float dtime=0.05;
	
	//material specific parameters
	//float den0=-0.085; //average density for pfc_3d_v1 (stripe phase)
	float den0=-0.285; //average density for pfc_3d_v2 (triangular phase)
	float tempr=-0.25; //temperature (T-Tm)/Tm, Tm=melting point
	//float tempr0=0.0; //positive value, tempr=tempr+tempr0/isteps;
	float noise=den0*1e-2; //Noise term to modulate the initial density field
	
	int ii; //ii=(i*Ny+j)*Nz+k;
	
	//----- ----- ----- ----- ----- -----
	//const int fftsizex = Nx, fftsizey = Ny;
	//
	cufftComplex *den_d, *f_den_d;
	cudaMalloc((void**)&den_d,    sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&f_den_d,  sizeof(cufftComplex)*NxNyNz);
	//
	cufftComplex *den3_d, *f_den3_d;
	cudaMalloc((void**)&den3_d,   sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&f_den3_d, sizeof(cufftComplex)*NxNyNz);
	//
	cufftComplex *ff_d, *f_ff_d;
	cudaMalloc((void**)&ff_d,     sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&f_ff_d,   sizeof(cufftComplex)*NxNyNz);
	//
	cufftHandle plan, iplan;
	cufftPlan3d(&plan,  Nx, Ny, Nz, CUFFT_C2C);
	cufftPlan3d(&iplan, Nx, Ny, Nz, CUFFT_C2C);
	//----- ----- ----- ----- ----- -----
	
	//----- ----- ----- ----- ----- -----
	//float kx[Nx];
	float *kx = (float *)malloc(sizeof(float)*( Nx ));
	//float ky[Ny];
	float *ky = (float *)malloc(sizeof(float)*( Ny ));
	//float kz[Nz];
	float *kz = (float *)malloc(sizeof(float)*( Nz ));
	//float k2[Nx][Ny][Nz];
	float *k2 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//float k4[Nx][Ny][Nz];
	float *k4 = (float *)malloc(sizeof(float)*( NxNyNz ));
	
	//prepare fft (output: kx,ky,kz,k2,k4)
	prepare_fft_3d(Nx,Ny,Nz,dx,dy,dz,kx,ky,kz,k2,k4); //get FFT coefficients
	//----- ----- ----- ----- ----- -----
	
	//----- for cufft
	float *k2_d, *k4_d;
	k2_d  = (float *)malloc(NxNyNz*sizeof(float));
	k4_d  = (float *)malloc(NxNyNz*sizeof(float));
	cudaMalloc((void**)&k2_d ,NxNyNz*sizeof(float));
	cudaMalloc((void**)&k4_d ,NxNyNz*sizeof(float));
	cudaMemcpy(k2_d,k2,NxNyNz*sizeof(float),cudaMemcpyHostToDevice); //k2 = k2_h
	cudaMemcpy(k4_d,k4,NxNyNz*sizeof(float),cudaMemcpyHostToDevice); //k4 = k4_h
	//----- ----- ----- -----
	
	//----- ----- ----- ----- ----- -----
	float _Complex *denc   = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *den3c  = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *f_den  = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *f_ff   = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *ffc    = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	//----- ----- ----- ----- ----- -----
	float *den = (float *)malloc(sizeof(float)*( NxNyNz ));
	float *ff  = (float *)malloc(sizeof(float)*( NxNyNz ));
	//----- ----- ----- ----- ----- -----
	
	//if infile=1 read input from file
	int infile = 0;
	FILE *in1;
	int mx,my,mz;
	if(infile==1){
		//open input file
		in1=fopen("g3_3r.inp","r");
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ii=(i*Ny+j)*Nz+k;
					fscanf(in1,"%5d %5d %5d %lf",&mx,&my,&mz,&den[ii]);
					denc[ii] = den[ii];
				}
			}
		}
	}else{
		//initialize density
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ii=(i*Ny+j)*Nz+k;
					den[ii] = den0 + noise*(0.5-(float)rand()/RAND_MAX);
					denc[ii] = den[ii];
				}
			}
		}
	}
	
	//----- ----- ----- ----- ----- -----
	float energy;
	//float ss2[Nx][Ny][Nz];
	float *ss2 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//float ss4[Nx][Ny][Nz];
	float *ss4 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//----- ----- ----- ----- ----- -----
	
	int bsx=BSX, bsy=BSY, bsz=BSZ;     //Number of threads
	dim3 blocks(Nx/bsx,Ny/bsy,Nz/bsz); //nx*ny*nz = blocks * threads
	dim3 threads(bsx,bsy,bsz);         //bsx*bsy*bsz <= 1024
	
	//evolve (evolve microstructure)
	for(int istep=0;istep<=nstep;istep++){
		
		//tempr = tempr + tempr0/istep;
		
		cudaMemcpy(den_d,denc,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //denc = denc_h
		
		//take current density filed from real space to Fourier space (forward FFT transformation)
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, den_d, f_den_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		//calculate the nonlinear term, phi^3, in Eq.7.10
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ii=(i*Ny+j)*Nz+k;
					den3c[ii]=denc[ii]*denc[ii]*denc[ii];
				}
			}
		}
		cudaMemcpy(den3_d,den3c,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //den3c = den3c_h
		
		//take the value of phi^3 from real space to Fourier space (forward FFT transformation)
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, den3_d, f_den3_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----  on cuda
		Kernel_semi_implicit_time_integration<<<blocks, threads>>>(Nx,Ny,Nz,
			dtime,tempr,
			k2_d,k4_d,
			f_den_d,f_den3_d);
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
		
		//bring back the values of phi^(t+1) from Fourier space to real space (inverse FFT transformation)
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(iplan, f_den_d, den_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cudaMemcpy(denc,den_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //denc = denc_h
		
		//print results
		//if print frequency is reached, output the results to file
		if(fmod(istep,nprint)==0){
		
			printf("done step: %5d \n", istep);
			
			//energy calculation
			//calculate the free energy distribution, Eq.7.6
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					for(int k=0;k<Nz;k++){
						ii=(i*Ny+j)*Nz+k;
						den[ii] = ( __real__ denc[ii] )/(NxNyNz);
						ss2[ii]=den[ii]*den[ii];
						ss4[ii]=ss2[ii]*ss2[ii];
					}
				}
			}
			
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					for(int k=0;k<Nz;k++){
						ii=(i*Ny+j)*Nz+k;
						f_ff[ii] = 0.5*f_den[ii]*(1.0-2.0*k2[ii]+k4[ii]);
					}
				}
			}
			cudaMemcpy(f_ff_d,f_ff,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //f_ff = f_ff_h
			
			//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
			cufftExecC2C(iplan, f_ff_d, ff_d, CUFFT_INVERSE); //IFFT
			cudaDeviceSynchronize();
			//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
			cudaMemcpy(ffc,ff_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //ffc = ffc_h
			
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					for(int k=0;k<Nz;k++){
						ii=(i*Ny+j)*Nz+k;
						ff[ii] = (( __real__ ffc[ii] )/(NxNyNz)) * den[ii]
						+ 0.5*tempr*ss2[ii]
						+ 0.25*ss4[ii];
					}
				}
			}
			
			//integrate the free energy field
			energy = 0.0;
			
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					for(int k=0;k<Nz;k++){
						ii=(i*Ny+j)*Nz+k;
						energy = energy + ff[ii];
					}
				}
			}
			
			//average free energy density
			energy = energy/(NxNyNz);
			
			//print the average free energy density value to file
			fprintf(out2, "%d %14.6e \n",istep, energy);
			
			//output the results in vtk file format for contour plots to be viewed by using paraview
			write_vtk_grid_values_3D(Nx,Ny,Nz,dx,dy,dz,istep,den);
			
		}
		
		//if intermediate configuration files are required, print the density field to file
		if(istep==nstep){
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					for(int k=0;k<Nz;k++){
						ii=(i*Ny+j)*Nz+k;
						fprintf(out1,"%5d %5d %5d %14.6e \n",i,j,k,den[ii]);
					}
				}
			}
		}
		
		//for recycle
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ii=(i*Ny+j)*Nz+k;
					denc[ii] = denc[ii]/(NxNyNz);
				}
			}
		}
		
	}//end of time step (evolve,for)
	
	//calculate the execution time and print it
	end = clock();
	compute_time = ((float) (end - start)) / CLOCKS_PER_SEC;
	printf("Compute Time: %lf \n", compute_time);
	
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
	//
	cudaFree(k2_d);
	cudaFree(k4_d);
	//----- ----- ----- ----- ----- -----
	free(kx);
	free(ky);
	free(k2);
	free(k4);
	//
	free(denc);
	free(den3c);
	free(f_den);
	free(f_ff);
	free(ffc);
	free(den);
	free(ff);
	//----- ----- ----- ----- ----- -----
	fclose(out1);
	fclose(out2);
	//----- ----- ----- ----- ----- -----
}
