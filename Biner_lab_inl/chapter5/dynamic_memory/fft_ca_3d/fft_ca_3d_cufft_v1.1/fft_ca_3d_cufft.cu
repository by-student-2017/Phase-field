/* 3D semi-implicit spectral phase-field code
   for solving Allen-Cahn equation */

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
void init_grain_micro_3d(int Nx, int Ny, int Nz,
	float dx, float dy, float dz,
	int iflag, int ngrain,
	float *etas, int *glist);
//----- ----- -----
float free_energy_fd_ca_3d(int i, int j, int k,
	int Nx, int Ny, int Nz,
	int ngrain, float *etas, float *eta, int igrain);
//----- ----- -----
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
	float coefA,
	float mobil,
	float grcoef,
	float *k2_d,
	cufftComplex *etak_d,
	cufftComplex *dfdetak_d
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
	denom = 1.0 + dtime*coefA*mobil*grcoef*k2_d[j];
	etak_d[j].x = ( etak_d[j].x - (dtime*mobil*dfdetak_d[j].x) )/denom; //real part
	etak_d[j].y = ( etak_d[j].y - (dtime*mobil*dfdetak_d[j].y) )/denom; //imaginary part
	
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
	
	//open an output file for writing total bulk energy values
	FILE *out2=fopen("area_frac.out","w");
	
	//simulation cell parameters (These values are dummy)
	int Nx=64; //Number of grid points in the x-direction
	int Ny=64; //Number of grid points in the y-direction
	int Nz=2; //Number of grid points in the y-direction
	
	int BSX=8; //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
	int BSY=8; //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
	int BSZ=2; //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
	
	int ngrain=2;
	
	//The distance between two grid points in x,y-direction
	float dx=0.5; //Grid spacing between two grid pints in x-direction
	float dy=0.5; //Grid spacing between two grid pints in y-direction
	float dz=0.5; //Grid spacing between two grid pints in z-direction
	
	//time integration parameters
	int nstep=100000; //Number of time integration steps
	int nprint=100; //Output frequency to write the results to file
	float dtime=0.005; //Time increment for the numerical integration
	float ttime=0.0;   //Total time
	float coefA=1.0;
	
	//material specific parameters
	
	//Initial concentrations of alloying elements
	float mobil=5.0;  //The value of mobility coefficient
	float grcoef=0.1; //The value of gradient energy coefficients [J(nm)^2/mol]
	
	//----- ----- ----- -----
	int ijk; //ijk=(i*Ny+j)*Nz+k;
	//----- ----- ----- -----
	
	/* Generate initial grain microstructure
	   iflag=1 is for bi-crystal and
	   iflag=2 is for polycrystal */
	int iflag=2;
	if(iflag==2){
		//read polycrystal microstructure
		FILE *in=fopen("grain_25.inp","r");
		fscanf(in,"%5d %5d %5d %5d ",&Nx,&Ny,&Nz,&ngrain);
		fclose(in);
	}
	//----- ----- ----- -----
	int NxNyNz=Nx*Ny*Nz; //Total number of grid points in the simulation cell
	float *etas = (float *)malloc(sizeof(float)*( NxNyNz*ngrain ));
	int   *glist = (int *)malloc(sizeof(int)*( ngrain ));
	//----- ----- ----- -----
	if(iflag==2){
		FILE *in=fopen("grain_25.inp","r");
		fscanf(in,"%5d %5d %5d %5d ",&Nx,&Ny,&Nz,&ngrain);
		//
		int nline=1;
		int ri, rj, rk, rigrain;
		float reta;
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					fscanf(in,"%5d %5d %5d %5d %f",&ri,&rj,&rk,&rigrain,&reta);
					//----- ----- ----- -----
					if( i != ri ){ printf("Don't match x data, Line %5d \n",nline); exit(1); }
					if( j != rj ){ printf("Don't match y data, LIne %5d \n",nline); exit(1); }
					if( k != rk ){ printf("Don't match z data, LIne %5d \n",nline); exit(1); }
					nline = nline + 1;
					//----- ----- ----- -----
					ijk=(i*Ny+j)*Nz+k;
					etas[ijk*ngrain+rigrain]=reta;
					//----- ----- ----- -----
				}
			}
		}
		fclose(in);
		//initialize glist
		for(int igrain=0;igrain<ngrain;igrain++){
			glist[igrain]=1.0;
		}
	}
	if(iflag==1){ init_grain_micro_3d(Nx,Ny,Nz,dx,dy,dz,iflag,ngrain,etas,glist); }
	
	//----- ----- ----- ----- ----- -----
	//const int fftsizex = Nx, fftsizey = Ny, fftsizez = Nz;
	//
	cufftComplex *eta_d, *dfdeta_d;
	cudaMalloc((void**)&eta_d,     sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&dfdeta_d,  sizeof(cufftComplex)*NxNyNz);
	//
	cufftComplex *etak_d, *dfdetak_d;
	cudaMalloc((void**)&etak_d,    sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&dfdetak_d, sizeof(cufftComplex)*NxNyNz);
	//
	cufftHandle plan, iplan;
	cufftPlan3d(&plan,  Nx, Ny, Nz, CUFFT_C2C);
	cufftPlan3d(&iplan, Nx, Ny, Nz, CUFFT_C2C);
	//----- ----- ----- ----- ----- -----
	
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
	k2_d  = (float *)malloc(Nx*Ny*sizeof(float));
	k4_d  = (float *)malloc(Nx*Ny*sizeof(float));
	cudaMalloc((void**)&k2_d ,Nx*Ny*sizeof(float));
	cudaMalloc((void**)&k4_d ,Nx*Ny*sizeof(float));
	cudaMemcpy(k2_d,k2,Nx*Ny*sizeof(float),cudaMemcpyHostToDevice); //k2 = k2_h
	cudaMemcpy(k4_d,k4,Nx*Ny*sizeof(float),cudaMemcpyHostToDevice); //k4 = k4_h
	//----- ----- ----- -----
	
	float *eta  = (float *)malloc(sizeof(float)*( NxNyNz ));
	float *eta2 = (float *)malloc(sizeof(float)*( NxNyNz ));
	
	//----- ----- ----- -----
	float denom;
	float _Complex numer;
	//----- ----- ----- -----
	float grain_sum;
	//----- ----- ----- -----
	int ncount;
	//----- ----- ----- -----
	
	float _Complex *dfdetac = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *etac    = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	
	int bsx=BSX, bsy=BSY, bsz=BSZ;     //Number of threads
	dim3 blocks(Nx/bsx,Ny/bsy,Nz/bsz); //nx*ny*nz = blocks * threads
	dim3 threads(bsx,bsy,bsz);         //bsx*bsy*bsz <= 1024
	
	//evolve (Time evolution of microstructure)
	for(int istep=0;istep<=nstep;istep++){
		
		//Update the total time
		ttime=ttime+dtime;
		
		//Loop over each grain
		for(int igrain=0;igrain<ngrain;igrain++){
			
			/* If glist is equal to one, which indicates that
			   the current grain area fraction is greater than 0.001,
			   continue the calculation. Otherwise, the current grain
			   does not exist anymore */
			if(glist[igrain]==1){
				
				/* Assign order parameters to temporary array eta[Nx][Ny] from
				   the common array etas[Nx][Ny][ngrain] for the current grain */
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						for(int k=0;k<Nz;k++){
							ijk=(i*Ny+j)*Nz+k;
							eta[ijk] = etas[ijk*ngrain+igrain];
							//-----
							etac[ijk] = eta[ijk];
						}
					}
				}
				
				//derivative of free energy
				//Calculate the derivative of elastic energy
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						for(int k=0;k<Nz;k++){
							ijk=(i*Ny+j)*Nz+k;
							dfdetac[ijk] = free_energy_fd_ca_3d(i,j,k,Nx,Ny,Nz,ngrain,etas,eta,igrain);
						}
					}
				}
				
				cudaMemcpy(dfdeta_d,dfdetac,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //dfdetac = dfdetac_h
				cudaMemcpy(eta_d,etac,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //etac = etac_h
				
				/* Take the values of concentration, derivative of free energy and
				   derivative of elastic energy from real space to Fourier space (forward FFT) */
				//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
				cufftExecC2C(plan, eta_d, etak_d, CUFFT_FORWARD); //FFT
				cudaDeviceSynchronize();
				//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
				cufftExecC2C(plan, dfdeta_d, dfdetak_d, CUFFT_FORWARD); //FFT
				cudaDeviceSynchronize();
				//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
				
				/* Semi-implicit time integration of concentration field at
				   Fourier space (Eq.5.14) */
				//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- on cuda
				Kernel_semi_implicit_time_integration<<<blocks, threads>>>(Nx,Ny,Nz,
					dtime,coefA,mobil,grcoef,
					k2_d,
					etak_d,dfdetak_d);
				cudaDeviceSynchronize();
				//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
				
				/* Take concentration field from Fourier space back to
				   real space (inverse FFT) */
				//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
				cufftExecC2C(iplan, etak_d, eta_d, CUFFT_INVERSE); //IFFT
				cudaDeviceSynchronize();
				//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
				
				cudaMemcpy(etac,eta_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //etac = etac_h
				
				//for small deviations
				// For small deviations from max and min values, reset the limits
				grain_sum=0.0;
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						for(int k=0;k<Nz;k++){
							ijk=(i*Ny+j)*Nz+k;
							//----- ----- ----- -----
							//eta[ijk] =  eta[ijk]/(NxNyNz);
							eta[ijk] = ( __real__ etac[ijk] )/(NxNyNz);
							//eta[ijk] =  creal(etac[ijk])/(NxNyNz); //For #include <_Complex.h>
							//----- ----- ----- -----
							if(eta[ijk]>=0.9999){
								eta[ijk]=0.9999;
							}
							if(eta[ijk]<=0.0001){
								eta[ijk]=0.0001;
							}
							//----- ----- ----- -----
							/* Calculate the total area of the current grain,
							   also return the order parameter values from
							   the temporary array eta[Nx][Ny] to 
							   common array etas[Nx][Ny][ngrain] */
							grain_sum = grain_sum + eta[ijk];
							etas[ijk*ngrain+igrain] = eta[ijk];
							//----- ----- ----- -----
						}
					}
				}
				
				//Check volume fraction of current grain
				/* Check the area fraction of the current grain.
				   If it is less than 0.001, set the its value in 
				   glist[ngrain] as zero which indicates that it is extinct. 
				   Also print message "grain # is eliminated" to screen. */
				grain_sum = grain_sum/NxNyNz;
				
				if(grain_sum<=0.001){
					glist[igrain]=0;
					printf("grain: No. %3d is eliminated \n",igrain);
				}
				
			}//end if(glist
		}//end igrain
		
		//print results
		// If print frequency reached, print the results to file
		if(fmod(istep,nprint)==0){
			//write vtk file & calculate are function of grains
			/* Prepare the data to be written to vtk file and
			   calculate the area fraction of each grain and
			   print them to file area_fract.out. */
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					for(int k=0;k<Nz;k++){
						ijk=(i*Ny+j)*Nz+k;
						eta2[ijk]=0.0;
					}
				}
			}
			fprintf(out2, "%14.6e ",ttime);
			
			for(int igrain=0;igrain<ngrain;igrain++){
				ncount=0;
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						for(int k=0;k<Nz;k++){
							ijk=(i*Ny+j)*Nz+k;
							//eta2[ijk]=eta2[ijk]+etas[ijk*ngrain+igrain]*etas[ijk*ngrain+igrain];
							eta2[ijk]=eta2[ijk]+etas[ijk*ngrain+igrain]*etas[ijk*ngrain+igrain]*igrain;
							if(etas[ijk*ngrain+igrain]>=0.5){
								ncount=ncount+1;
							}//end if
						}//end for(k
					}//end for(j
				}//end for(i
				ncount=ncount/NxNyNz;
				fprintf(out2, "%5d ",ncount);
			}//end igrain
			fprintf(out2, "\n");
			
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_3D(Nx,Ny,Nz,dx,dy,dz,istep,eta2);
			
			printf("done step: %5d \n",istep);
		}//end if
	}//end for(istep
	
	//calculate the execution time and print it
	/* Calculate the compute time and print it to screen */
	end = clock();
	compute_time = ((float) (end - start)) / CLOCKS_PER_SEC;
	printf("Compute Time: %f \n", compute_time);
	
	//----- ----- ----- ----- ----- -----
	cufftDestroy(plan);
	cufftDestroy(iplan);
	//----- ----- ----- ----- ----- -----
	cudaFree(eta_d);
	cudaFree(dfdeta_d);
	//
	cudaFree(etak_d);
	cudaFree(dfdetak_d);
	//----- ----- ----- ----- ----- -----
	free(kx);
	free(ky);
	free(k2);
	free(k4);
	//
	free(etas);
	free(glist);
	free(eta2);
	free(eta);
	//
	free(etac);
	free(dfdetac);
	//----- ----- ----- ----- ----- -----
}
