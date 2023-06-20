/* This function evaluates the derivative of elastic energy with
   respect to concentration. First, stress and strain values are
   solved with the iterative algorithm described earlier, 
   then derivative of elastic energy is evaluated for all grid points. */

#include <stdlib.h> //rand() and malloc
#include <math.h>
//#include <fftw3.h>

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

/* Variable and array list
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  Nz: Number of grid points in the z-direction
  cm11: C11 component of elasticity matrix for matrix material
  cm12: C12 component of elasticity matrix for matrix material
  cm44: C44 component of elasticity matrix for matrix material
  cp11: C11 component of elasticity matrix for second phase
  cp12: C12 component of elasticity matrix for second phase
  cp44: C44 component of elasticity matrix for second phase
  ed11[Nx][Ny][Nz] to ed13[Nx][Ny][Nz]: Strain component of lattice defects
  ea[6]: Applied strains
  con[Nx][Ny][Nz]: Concentration
  s11[Nx][Ny][Nz] to s13[Nx][Ny][Nz]: Component of stress
  e11[Nx][Ny][Nz] to e13[Nx][Ny][Nz]: Component of strain
  delsdc[Nx][Ny][Nz]: Functional derivative of elastic energy
  //
  omeg11[Nx][Ny][Nz]: Coefficient needed for the Green's tensor
  omeg22[Nx][Ny][Nz]: Coefficient needed for the Green's tensor
  omeg33[Nx][Ny][Nz]: Coefficient needed for the Green's tensor
  omeg12[Nx][Ny][Nz]: Coefficient needed for the Green's tensor
  omeg23[Nx][Ny][Nz]: Coefficient needed for the Green's tensor
  omeg13[Nx][Ny][Nz]: Coefficient needed for the Green's tensor
  tmatx[3][3][3][3]: Green's tensor at i,j,k grid point (real part only)
*/

// Define subroutine "Kernel" for GPU (Device) calculation in detail
__global__ void Kernel_initialization(
	int   Nx, int   Ny, int   Nz,
	float *conr_d,
	float *ei11_d, float *ei22_d, float *ei33_d,
	float *ei12_d, float *ei23_d, float *ei13_d,
	float ei0,
	float *c11_d,  float *c12_d,  float *c44_d,
	float cp11,    float cp12,    float cp44,
	float cm11,    float cm12,    float cm44
){
	int j, jx, jy, jz;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	jx = blockDim.x*blockIdx.x + threadIdx.x; //<-GPU | CPU -> for(jx=0; jx<nx; jx++){
	jy = blockDim.y*blockIdx.y + threadIdx.y; //<-GPU | CPU -> for(jy=0; jy<ny; jy++){
	jz = blockDim.z*blockIdx.z + threadIdx.z; //<-GPU | CPU -> for(jz=0; jz<nz; jz++){
	j  = (jz*Ny + jy)*Nx + jx; //j = nx*ny*jz + nx*jy + jx;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- -----
	// Calculate the eigenstrains (head ei series)
	ei11_d[j] = ei0*conr_d[j];
	ei22_d[j] = ei0*conr_d[j];
	ei33_d[j] = ei0*conr_d[j];
	//
	//ei12_d[j] = 0.0*conr_d[j];
	//ei23_d[j] = 0.0*conr_d[j];
	//ei13_d[j] = 0.0*conr_d[j];
	ei12_d[j] = 0.0;
	ei23_d[j] = 0.0;
	ei13_d[j] = 0.0;
	//----- ----- ----- ----- -----
	
	/* Calculate the effective elastic constants at 
	   the grid points based on the composition and
	   using Vegard's law */
	//----- ----- ----- ----- -----
	c11_d[j] = conr_d[j]*cp11 + (1.0-conr_d[j])*cm11;
	c12_d[j] = conr_d[j]*cp12 + (1.0-conr_d[j])*cm12;
	c44_d[j] = conr_d[j]*cp44 + (1.0-conr_d[j])*cm44;
	//----- ----- ----- ----- -----
}

void green_tensor2_3D(int Nx, int Ny, int Nz,
	float *kx, float *ky, float *kz,
	float *omeg11, float *omeg22, float *omeg33,
	float *omeg12, float *omeg23, float *omeg13,
	int i, int j, int k,
	float *tmatx);

void solve_elasticity_3d(int Nx, int Ny, int Nz,
	int BSX, int BSY, int BSZ,
	float *kx, float *ky, float *kz,
	float *omeg11, float *omeg22, float *omeg33,
	float *omeg12, float *omeg23, float *omeg13,
	float _Complex *s11, float _Complex *s22, float _Complex *s33,
	float _Complex *s12, float _Complex *s23, float _Complex *s13,
	float _Complex *e11, float _Complex *e22, float _Complex *e33,
	float _Complex *e12, float _Complex *e23, float _Complex *e13,
	float *ed11, float *ed22, float *ed33,
	float *ed12, float *ed23, float *ed13,
	float cm11, float cm12, float cm44,
	float cp11, float cp12, float cp44,
	float *ea,
	float ei0,
	float *con,  float _Complex *delsdc){
	
	int NxNyNz=Nx*Ny*Nz;
	
	//----- ----- ----- -----
	cufftComplex *s11_d, *s22_d, *s33_d;
	cudaMalloc((void**)&s11_d, sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&s22_d, sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&s33_d, sizeof(cufftComplex)*NxNyNz);
	//
	cufftComplex *s12_d, *s23_d, *s13_d;
	cudaMalloc((void**)&s12_d, sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&s23_d, sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&s13_d, sizeof(cufftComplex)*NxNyNz);
	//----- ----- ----- -----
	cufftComplex *e11_d, *e22_d, *e33_d;
	cudaMalloc((void**)&e11_d, sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&e22_d, sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&e33_d, sizeof(cufftComplex)*NxNyNz);
	//
	cufftComplex *e12_d, *e23_d, *e13_d;
	cudaMalloc((void**)&e12_d, sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&e23_d, sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&e13_d, sizeof(cufftComplex)*NxNyNz);
	//----- ----- ----- -----
	
	//----- ----- ----- -----
	cufftComplex *s11k_d, *s22k_d, *s33k_d;
	cudaMalloc((void**)&s11k_d, sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&s22k_d, sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&s33k_d, sizeof(cufftComplex)*NxNyNz);
	//
	cufftComplex *s12k_d, *s23k_d, *s13k_d;
	cudaMalloc((void**)&s12k_d, sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&s23k_d, sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&s13k_d, sizeof(cufftComplex)*NxNyNz);
	//----- ----- ----- -----
	cufftComplex *e11k_d, *e22k_d, *e33k_d;
	cudaMalloc((void**)&e11k_d, sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&e22k_d, sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&e33k_d, sizeof(cufftComplex)*NxNyNz);
	//
	cufftComplex *e12k_d, *e23k_d, *e13k_d;
	cudaMalloc((void**)&e12k_d, sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&e23k_d, sizeof(cufftComplex)*NxNyNz);
	cudaMalloc((void**)&e13k_d, sizeof(cufftComplex)*NxNyNz);
	//----- ----- ----- -----
	cufftHandle plan, iplan;
	cufftPlan3d(&plan,  Nx, Ny, Nz, CUFFT_C2C);
	cufftPlan3d(&iplan, Nx, Ny, Nz, CUFFT_C2C);
	//----- ----- ----- -----
	
	float _Complex *s11k = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *s22k = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *s33k = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *s12k = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *s23k = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *s13k = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	//
	float _Complex *e11k = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *e22k = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *e33k = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *e12k = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *e23k = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *e13k = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	
	//----- ----- ----- ----- ----- ----- -----
	// eigenstrains (head ei series) components
	//float ei11[Nx][Ny][Nz];
	float *ei11 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//float ei22[Nx][Ny][Nz];
	float *ei22 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//float ei33[Nx][Ny][Nz];
	float *ei33 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//
	//float ei12[Nx][Ny][Nz];
	float *ei12 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//float ei23[Nx][Ny][Nz];
	float *ei23 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//float ei13[Nx][Ny][Nz];
	float *ei13 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//----- ----- ----- ----- ----- ----- -----
	
	//----- ----- ----- ----- ----- ----- -----
	// elastic modulus components
	//float c11[Nx][Ny][Nz];
	float  *c11 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//float c12[Nx][Ny][Nz];
	float  *c12 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//float c44[Nx][Ny][Nz];
	float  *c44 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//----- ----- ----- ----- ----- ----- -----
	
	//----- ----- -----
	int ijk;
	int klij; //For tmatx
	//----- ----- -----
	
	//----- -----
	float *ei11_d, *ei22_d, *ei33_d; // name of dynamic memory for GPU, CUDA, device
	ei11_d = (float *)malloc(NxNyNz*sizeof(float)); //GPU, CUDA, device
	ei22_d = (float *)malloc(NxNyNz*sizeof(float)); //GPU, CUDA, device
	ei33_d = (float *)malloc(NxNyNz*sizeof(float)); //GPU, CUDA, device
	//
	float *ei12_d, *ei23_d, *ei13_d; // name of dynamic memory for GPU, CUDA, device
	ei12_d = (float *)malloc(NxNyNz*sizeof(float)); //GPU, CUDA, device
	ei23_d = (float *)malloc(NxNyNz*sizeof(float)); //GPU, CUDA, device
	ei13_d = (float *)malloc(NxNyNz*sizeof(float)); //GPU, CUDA, device
	//-----
	cudaMalloc((void**)&ei11_d,NxNyNz*sizeof(float)); // define dynamic memory for GPU (device)
	cudaMalloc((void**)&ei22_d,NxNyNz*sizeof(float)); // define dynamic memory for GPU (device)
	cudaMalloc((void**)&ei33_d,NxNyNz*sizeof(float)); // define dynamic memory for GPU (device)
	//
	cudaMalloc((void**)&ei12_d,NxNyNz*sizeof(float)); // define dynamic memory for GPU (device)
	cudaMalloc((void**)&ei23_d,NxNyNz*sizeof(float)); // define dynamic memory for GPU (device)
	cudaMalloc((void**)&ei13_d,NxNyNz*sizeof(float)); // define dynamic memory for GPU (device)
	//-----
	cudaMemcpy(ei11_d,ei11,NxNyNz*sizeof(float),cudaMemcpyHostToDevice); //ei11 = ei11_h
	cudaMemcpy(ei22_d,ei22,NxNyNz*sizeof(float),cudaMemcpyHostToDevice); //ei22 = ei22_h
	cudaMemcpy(ei33_d,ei33,NxNyNz*sizeof(float),cudaMemcpyHostToDevice); //ei12 = ei12_h
	//
	cudaMemcpy(ei12_d,ei12,NxNyNz*sizeof(float),cudaMemcpyHostToDevice); //ei11 = ei11_h
	cudaMemcpy(ei23_d,ei23,NxNyNz*sizeof(float),cudaMemcpyHostToDevice); //ei22 = ei22_h
	cudaMemcpy(ei13_d,ei13,NxNyNz*sizeof(float),cudaMemcpyHostToDevice); //ei12 = ei12_h
	//----- -----
	float *c11_d, *c12_d, *c44_d; // name of dynamic memory for GPU, CUDA, device
	c11_d = (float *)malloc(NxNyNz*sizeof(float)); //GPU, CUDA, device
	c12_d = (float *)malloc(NxNyNz*sizeof(float)); //GPU, CUDA, device
	c44_d = (float *)malloc(NxNyNz*sizeof(float)); //GPU, CUDA, device
	//-----
	cudaMalloc((void**)&c11_d,NxNyNz*sizeof(float)); // define dynamic memory for GPU (device)
	cudaMalloc((void**)&c12_d,NxNyNz*sizeof(float)); // define dynamic memory for GPU (device)
	cudaMalloc((void**)&c44_d,NxNyNz*sizeof(float)); // define dynamic memory for GPU (device)
	//-----
	cudaMemcpy(c11_d,c11,NxNyNz*sizeof(float),cudaMemcpyHostToDevice); //c11 = c11_h
	cudaMemcpy(c12_d,c12,NxNyNz*sizeof(float),cudaMemcpyHostToDevice); //c12 = c12_h
	cudaMemcpy(c44_d,c44,NxNyNz*sizeof(float),cudaMemcpyHostToDevice); //c44 = c44_h
	//----- -----
	float *conr_d;
	conr_d = (float *)malloc(NxNyNz*sizeof(float)); //GPU, CUDA, device
	cudaMalloc((void**)&conr_d,NxNyNz*sizeof(float)); // define dynamic memory for GPU (device)
	cudaMemcpy(conr_d,con,NxNyNz*sizeof(float),cudaMemcpyHostToDevice); //con = con_h
	//----- -----
	
	int bsx=BSX, bsy=BSY, bsz=BSZ;     //Number of threads
	dim3 blocks(Nx/bsx,Ny/bsy,Nz/bsz); //nx*ny*nz = blocks * threads
	dim3 threads(bsx,bsy,bsz);         //bsx*bsy*bsz <= 1024
	
	Kernel_initialization<<<blocks, threads>>>(Nx,Ny,Nz,
		conr_d,
		ei11_d, ei22_d, ei33_d,
		ei12_d, ei23_d, ei13_d,
		ei0,
		c11_d,  c12_d,  c44_d,
		cp11,   cp12,   cp44,
		cm11,   cm12,   cm44);
	cudaDeviceSynchronize();
	
	//----- -----
	cudaMemcpy(ei11,ei11_d,NxNyNz*sizeof(float),cudaMemcpyDeviceToHost); //ei11 = ei11_h
	cudaMemcpy(ei22,ei22_d,NxNyNz*sizeof(float),cudaMemcpyDeviceToHost); //ei22 = ei22_h
	cudaMemcpy(ei33,ei33_d,NxNyNz*sizeof(float),cudaMemcpyDeviceToHost); //ei33 = ei33_h
	//
	cudaMemcpy(ei12,ei12_d,NxNyNz*sizeof(float),cudaMemcpyDeviceToHost); //ei12 = ei12_h
	cudaMemcpy(ei23,ei23_d,NxNyNz*sizeof(float),cudaMemcpyDeviceToHost); //ei23 = ei23_h
	cudaMemcpy(ei13,ei13_d,NxNyNz*sizeof(float),cudaMemcpyDeviceToHost); //ei13 = ei13_h
	//----- -----
	cudaMemcpy(c11,c11_d,NxNyNz*sizeof(float),cudaMemcpyDeviceToHost); //c11 = c11_h
	cudaMemcpy(c12,c12_d,NxNyNz*sizeof(float),cudaMemcpyDeviceToHost); //c12 = c12_h
	cudaMemcpy(c44,c44_d,NxNyNz*sizeof(float),cudaMemcpyDeviceToHost); //c44 = c44_h
	//----- -----
	
	//for(int i=0;i<Nx;i++){
	//	for(int j=0;j<Ny;j++){
	//		for(int k=0;k<Nz;k++){
	//			ijk=(i*Ny+j)*Nz+k;
	//			//----- ----- ----- ----- -----
	//			// Calculate the eigenstrains (head ei series)
	//			ei11[ijk] = ei0*con[ijk];
	//			ei22[ijk] = ei0*con[ijk];
	//			ei33[ijk] = ei0*con[ijk];
	//			//
	//			ei12[ijk] = 0.0*con[ijk];
	//			ei23[ijk] = 0.0*con[ijk];
	//			ei13[ijk] = 0.0*con[ijk];
	//			//----- ----- ----- ----- -----
	//			
	//			/* Calculate the effective elastic constants at 
	//			   the grid points based on the composition and
	//			   using Vegard's law */
	//			//----- ----- ----- ----- -----
	//			c11[ijk] = con[ijk]*cp11 + (1.0-con[ijk])*cm11;
	//			c12[ijk] = con[ijk]*cp12 + (1.0-con[ijk])*cm12;
	//			c44[ijk] = con[ijk]*cp44 + (1.0-con[ijk])*cm44;
	//			//----- ----- ----- ----- -----
	//		}
	//	}
	//}
	
	/* Note: elastic modulus in this case. */
	//----- ----- ----- -----
	//float c22=c33=c11;
	//float c21=c12;
	//float c31=c13=c12;
	//float c32=c23=c12;
	//float c55=c66=c44;
	//----- ----- ----- -----
	float cm22,cm33;
		cm22=cm33=cm11;
	float cm21;
		cm21=cm12;
	float cm31,cm13;
		cm31=cm13=cm12;
	float cm32,cm23;
		cm32=cm23=cm12;
	float cm55, cm66;
		cm55=cm66=cm44;
	//----- ----- ----- -----
	float cp22, cp33;
		cp22=cp33=cp11;
	float cp21;
		cp21=cp12;
	float cp31, cp13;
		cp31=cp13=cp12;
	float cp32, cp23;
		cp32=cp23=cp12;
	float cp55, cp66;
		cp55=cp66=cp44;
	//----- ----- ----- -----
	float _Complex et21;
	float _Complex et32;
	float _Complex et31;
	//----- ----- ----- -----
	
	//----- ----- ----- ----- ----- ----- -----
	//float smatx_real[Nx][Ny][Nz][3][3];
	//float *smatx_real = (float *)malloc(sizeof(float)*( NxNyNz*3*3 ));
	//
	//float smatx_imag[Nx][Ny][Nz][3][3];
	//float *smatx_imag = (float *)malloc(sizeof(float)*( NxNyNz*3*3 ));
	//----- ----- ----- ----- ----- ----- -----
	//float ematx_real[Nx][Ny][Nz][3][3];
	//float *ematx_real = (float *)malloc(sizeof(float)*( NxNyNz*3*3 ));
	//
	//float ematx_imag[Nx][Ny][Nz][3][3];
	//float *ematx_imag = (float *)malloc(sizeof(float)*( NxNyNz*3*3 ));
	//----- ----- ----- ----- ----- ----- -----
	float _Complex *smatx = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz*3*3 ));
	float _Complex *ematx = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz*3*3 ));
	//----- ----- ----- ----- ----- ----- -----
	
	//float tmatx[3][3][3][3];
	float *tmatx = (float *)malloc(sizeof(float)*( 3*3*3*3 )); //real part only
	
	//----- ----- -----
	float _Complex et11=0.0;
	float _Complex et22=0.0;
	float _Complex et33=0.0;
	//
	float _Complex et12=0.0;
	float _Complex et23=0.0;
	float _Complex et13=0.0;
	//----- ----- -----
	
	//float sum_stress[Nx][Ny][Nz];
	float *sum_stress = (float *)malloc(sizeof(float)*( NxNyNz ));
	
	//----- ----- -----
	//Maximum number of iteration steps
	int niter=10;
	//----- ----- -----
	float old_norm=0.0;
	float normF=0.0;
	//----- ----- -----
	float conver=0.0;
	//----- ----- -----
	//Tolerance value of convergence tests
	float tolerance=0.001;
	//----- ----- -----
	
	/* Solve stress and strain field with 
	   iterative algorithm given in the text */
	for(int iter=0;iter<niter;iter++){
		
		cudaMemcpy(s11_d,s11,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //s11 = s11_h
		cudaMemcpy(s22_d,s22,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //s22 = s22_h
		cudaMemcpy(s33_d,s33,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //s33 = s33_h
		cudaMemcpy(s12_d,s12,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //s12 = s12_h
		cudaMemcpy(s23_d,s23,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //s23 = s23_h
		cudaMemcpy(s13_d,s13,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //s13 = s13_h
		//
		cudaMemcpy(e11_d,e11,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //e11 = e11_h
		cudaMemcpy(e22_d,e22,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //e22 = e22_h
		cudaMemcpy(e33_d,e33,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //e33 = e33_h
		cudaMemcpy(e12_d,e12,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //e12 = e12_h
		cudaMemcpy(e23_d,e23,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //e23 = e23_h
		cudaMemcpy(e13_d,e13,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //e13 = e13_h
		
		/* Take stress and strain components from real space to
		   Fourier space (forward FFT). Step-a */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		// stress (head s series)
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, s11_d, s11k_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, s22_d, s22k_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, s33_d, s33k_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, s12_d, s12k_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, s23_d, s23k_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, s13_d, s13k_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		// strain (head e series)
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, e11_d, e11k_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, e22_d, e22k_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, e33_d, e33k_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, e12_d, e12k_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, e23_d, e23k_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, e13_d, e13k_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		cudaMemcpy(s11k,s11k_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //s11k = s11k_h
		cudaMemcpy(s22k,s22k_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //s22k = s22k_h
		cudaMemcpy(s33k,s33k_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //s33k = s33k_h
		cudaMemcpy(s12k,s12k_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //s12k = s12k_h
		cudaMemcpy(s23k,s23k_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //s23k = s23k_h
		cudaMemcpy(s13k,s13k_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //s13k = s13k_h
		//
		cudaMemcpy(e11k,e11k_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //e11k = e11k_h
		cudaMemcpy(e22k,e22k_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //e22k = e22k_h
		cudaMemcpy(e33k,e33k_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //e33k = e33k_h
		cudaMemcpy(e12k,e12k_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //e12k = e12k_h
		cudaMemcpy(e23k,e23k_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //e23k = e23k_h
		cudaMemcpy(e13k,e13k_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //e13k = e13k_h
		
		//Green operator
		// Calculate strain tensor, Eq.5.46, Step-b
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					
					/* Form stress and strain tensors to be used in 
					   Eq.5.46, Step-b */
					//----- ----- ----- ----- ----- -----
					// stress (smatx, sXY and sXYk)
					//----- ----- ----- ----- ----- -----
					smatx[(ijk*3+0)*3+0]=s11k[ijk];
					smatx[(ijk*3+0)*3+1]=s12k[ijk];
					smatx[(ijk*3+0)*3+2]=s13k[ijk];
					//
					smatx[(ijk*3+1)*3+0]=s12k[ijk];
					smatx[(ijk*3+1)*3+1]=s22k[ijk];
					smatx[(ijk*3+1)*3+2]=s23k[ijk];
					//
					smatx[(ijk*3+2)*3+0]=s13k[ijk];
					smatx[(ijk*3+2)*3+1]=s23k[ijk];
					smatx[(ijk*3+2)*3+2]=s33k[ijk];
					//----- ----- ----- ----- ----- -----
					//
					//----- ----- ----- ----- ----- -----
					// strain (ematx, eXY and eXYk)
					//----- ----- ----- ----- ----- -----
					ematx[(ijk*3+0)*3+0]=e11k[ijk];
					ematx[(ijk*3+0)*3+1]=e12k[ijk];
					ematx[(ijk*3+0)*3+2]=e13k[ijk];
					//
					ematx[(ijk*3+1)*3+0]=e12k[ijk];
					ematx[(ijk*3+1)*3+1]=e22k[ijk];
					ematx[(ijk*3+1)*3+2]=e23k[ijk];
					//
					ematx[(ijk*3+2)*3+0]=e13k[ijk];
					ematx[(ijk*3+2)*3+1]=e23k[ijk];
					ematx[(ijk*3+2)*3+2]=e33k[ijk];
					//----- ----- ----- ----- ----- -----
					
					//----- ----- ----- ----- ----- -----
					green_tensor2_3D(Nx,Ny,Nz,
									kx,ky,kz,
									omeg11,omeg22,omeg33,
									omeg12,omeg23,omeg13,
									i,j,k,
									tmatx);
					//----- ----- ----- ----- ----- -----
					for(int kk=0;kk<3;kk++){
						for(int ll=0;ll<3;ll++){
							for(int ii=0;ii<3;ii++){
								for(int jj=0;jj<3;jj++){
									klij=((kk*3+ll)*3+ii)*3+jj;
									/* Eq.5.46(b): new epsilon(zeta) = epsilon(zeta) - sum( gamma(zeta)*sigma(zeta) )
									   where gamma=tmatx, sigma=smatx
									   Note: tmatx is real part only */
									//----- ----- ----- ----- ----- -----
									//ematx_real[(ijk*3+ii)*3+jj] -= tmatx[klij]*smatx_real[(ijk*3+kk)*3+ll];
									//
									//ematx_imag[(ijk*3+ii)*3+jj] -= tmatx[klij]*smatx_imag[(ijk*3+kk)*3+ll];
									//----- ----- ----- ----- ----- -----
									ematx[(ijk*3+ii)*3+jj] -= tmatx[klij]*smatx[(ijk*3+kk)*3+ll];
									//----- ----- ----- ----- ----- -----
								}//jj
							}//ii
						}//ll
					}//kk
					//----- ----- ----- ----- ----- -----
					
					// Rearrange strain components using symmetry of strain tensor
					//----- ----- ----- ----- ----- -----
					// strain (ematx, eXY and eXYk)
					//----- ----- ----- ----- ----- -----
					e11k[ijk]=ematx[(ijk*3+0)*3+0];
					e12k[ijk]=ematx[(ijk*3+0)*3+1];
					e13k[ijk]=ematx[(ijk*3+0)*3+2];
					//
					//e12k[ijk]=ematx[(ijk*3+1)*3+0];
					e22k[ijk]=ematx[(ijk*3+1)*3+1];
					e23k[ijk]=ematx[(ijk*3+1)*3+2];
					//
					//e13k[ijk]=ematx[(ijk*3+2)*3+0];
					//e23k[ijk]=ematx[(ijk*3+2)*3+1];
					e33k[ijk]=ematx[(ijk*3+2)*3+2];
					//----- ----- ----- ----- ----- -----
					
				}//Nz
			}//Ny
		}//Nx
		
		cudaMemcpy(e11k_d,e11k,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //e11k = e11k_h
		cudaMemcpy(e22k_d,e22k,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //e22k = e22k_h
		cudaMemcpy(e33k_d,e33k,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //e33k = e33k_h
		cudaMemcpy(e12k_d,e12k,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //e12k = e12k_h
		cudaMemcpy(e23k_d,e23k,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //e23k = e23k_h
		cudaMemcpy(e13k_d,e13k,NxNyNz*sizeof(float _Complex),cudaMemcpyHostToDevice); //e13k = e13k_h
		
		//From Fourier space to real space
		/* Take strain components from Fourier space back to
		   real space (inverse FFT), Step-c */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(iplan, e11k_d, e11_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(iplan, e22k_d, e22_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(iplan, e33k_d, e33_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(iplan, e12k_d, e12_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(iplan, e23k_d, e23_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(iplan, e13k_d, e13_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		cudaMemcpy(e11,e11_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //e11 = e11_h
		cudaMemcpy(e22,e22_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //e22 = e22_h
		cudaMemcpy(e33,e33_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //e33 = e33_h
		cudaMemcpy(e12,e12_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //e12 = e12_h
		cudaMemcpy(e23,e23_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //e23 = e23_h
		cudaMemcpy(e13,e13_d,NxNyNz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //e13 = e13_h
		
		//Calculate stresses
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					//----- ----- ----- ----- ----- -----
					// strain (ematx, eXY and eXYk)
					//----- ----- ----- -----
					e11[ijk] /= NxNyNz;
					e22[ijk] /= NxNyNz;
					e33[ijk] /= NxNyNz;
					//
					e12[ijk] /= NxNyNz;
					e23[ijk] /= NxNyNz;
					e13[ijk] /= NxNyNz;
					//----- ----- ----- ----- ----- -----
					
					//----- ----- ----- ----- ----- ----- ----- -----
					/*s11[ijk][0]=c11[ijk]*(ea[0]+e11[ijk][0]-ei11[ijk]-ed11[ijk])
							   +c12[ijk]*(ea[1]+e22[ijk][0]-ei22[ijk]-ed22[ijk])
							   +c13[ijk]*(ea[2]+e33[ijk][0]-ei33[ijk]-ed33[ijk]);
					s22[ijk][0]=c21[ijk]*(ea[0]+e11[ijk][0]-ei11[ijk]-ed11[ijk])
							   +c22[ijk]*(ea[1]+e22[ijk][0]-ei22[ijk]-ed22[ijk])
							   +c23[ijk]*(ea[2]+e33[ijk][0]-ei33[ijk]-ed33[ijk]);
					s22[ijk][0]=c31[ijk]*(ea[0]+e11[ijk][0]-ei11[ijk]-ed11[ijk])
							   +c32[ijk]*(ea[1]+e22[ijk][0]-ei22[ijk]-ed22[ijk])
							   +c33[ijk]*(ea[2]+e33[ijk][0]-ei33[ijk]-ed33[ijk]);*/
					//
					// c13[ijk]=c12[ijk], c23[ijk]=c13[ijk], c22[ijk]=c11[ijk], c33[ijk]=c11[ijk], etc
					/*s11[ijk][0]=c11[ijk]*(ea[0]+e11[ijk][0]-ei11[ijk]-ed11[ijk])
							   +c12[ijk]*(ea[1]+e22[ijk][0]-ei22[ijk]-ed22[ijk])
							   +c12[ijk]*(ea[2]+e33[ijk][0]-ei33[ijk]-ed33[ijk]);
					s22[ijk][0]=c12[ijk]*(ea[0]+e11[ijk][0]-ei11[ijk]-ed11[ijk])
							   +c11[ijk]*(ea[1]+e22[ijk][0]-ei22[ijk]-ed22[ijk])
							   +c12[ijk]*(ea[2]+e33[ijk][0]-ei33[ijk]-ed33[ijk]);
					s33[ijk][0]=c12[ijk]*(ea[0]+e11[ijk][0]-ei11[ijk]-ed11[ijk])
							   +c12[ijk]*(ea[1]+e22[ijk][0]-ei22[ijk]-ed22[ijk])
							   +c11[ijk]*(ea[2]+e33[ijk][0]-ei33[ijk]-ed33[ijk]);*/
					//
					// Calculate strain (head e series) components
					et11=( ea[0]+e11[ijk]-ei11[ijk]-ed11[ijk] );
					et22=( ea[1]+e22[ijk]-ei22[ijk]-ed22[ijk] );
					et33=( ea[2]+e33[ijk]-ei33[ijk]-ed33[ijk] );
					//
					// Calculate stress (head s series) components
					s11[ijk]=( c11[ijk]*et11 + c12[ijk]*et22 + c12[ijk]*et33 );
					s22[ijk]=( c12[ijk]*et11 + c11[ijk]*et22 + c12[ijk]*et33 );
					s33[ijk]=( c12[ijk]*et11 + c12[ijk]*et22 + c11[ijk]*et33 );
					//
					//s11[ijk][1]=0.0;
					//s22[ijk][1]=0.0;
					//s33[ijk][1]=0.0;
					//----- ----- ----- ----- ----- ----- ----- -----
					
					//----- ----- ----- ----- ----- ----- ----- -----
					/*s12[ijk][0]=c44[ijk]*(ea[3]+e12[ijk][0]-ei12[ijk]-ed12[ijk])
							   +c44[ijk]*(ea[3]+e21[ijk][0]-ei21[ijk]-ed21[ijk]);
					s23[ijk][0]=c55[ijk]*(ea[4]+e23[ijk][0]-ei23[ijk]-ed23[ijk])
							   +c55[ijk]*(ea[4]+e32[ijk][0]-ei32[ijk]-ed32[ijk]);
					s13[ijk][0]=c66[ijk]*(ea[5]+e13[ijk][0]-ei13[ijk]-ed13[ijk])
							   +c66[ijk]*(ea[5]+e31[ijk][0]-ei31[ijk]-ed31[ijk]);*/
					//
					// c55[ijk]=c44[ijk], c66[ijk]=c44[ijk], etc
					/*s12[ijk][0]=c44[ijk]*(ea[3]+e12[ijk][0]-ei12[ijk]-ed12[ijk])*2.0;
					  s23[ijk][0]=c44[ijk]*(ea[4]+e23[ijk][0]-ei23[ijk]-ed23[ijk])*2.0;
					  s13[ijk][0]=c44[ijk]*(ea[5]+e13[ijk][0]-ei13[ijk]-ed13[ijk])*2.0;*/
					//
					// Calculate strain (head e series) components
					et12=( ea[3]+e12[ijk]-ei12[ijk]-ed12[ijk] );
					et23=( ea[4]+e23[ijk]-ei23[ijk]-ed23[ijk] );
					et13=( ea[5]+e13[ijk]-ei13[ijk]-ed13[ijk] );
					//
					// Calculate stress (head s series) components
					s12[ijk]=c44[ijk]*et12*2.0;
					s23[ijk]=c44[ijk]*et23*2.0;
					s13[ijk]=c44[ijk]*et13*2.0;
					//
					//s12[ijk][1]=0.0;
					//s23[ijk][1]=0.0;
					//s13[ijk][1]=0.0;
					//----- ----- ----- ----- ----- ----- ----- -----
					
					//check convergence
					sum_stress[ijk] = __real__ ( s11[ijk] + s22[ijk] + s33[ijk]
												+s12[ijk] + s23[ijk] + s13[ijk] );
					//normF=norm(sum_stress,2.0);
					normF = normF + sum_stress[ijk]*sum_stress[ijk];
				}
			}
		}
		
		normF=sqrt(normF);
		
		if(iter==1){
			conver=fabs((normF-old_norm)/(old_norm));
			if(conver<=tolerance){
				break;
			}
		}
		old_norm=normF;
		
	}//end iter
	
	//strain energy
	//Calculate functional derivative of elastic energy
	// sum strain components
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				ijk=(i*Ny+j)*Nz+k;
				//----- ----- ----- ----- ----- ----- ----- -----
				// Calculate strain (head e series) components
				et11=( ea[0]+e11[ijk]-ei11[ijk]-ed11[ijk] );
				et22=( ea[1]+e22[ijk]-ei22[ijk]-ed22[ijk] );
				et33=( ea[2]+e33[ijk]-ei33[ijk]-ed33[ijk] );
				//
				et12=( ea[3]+e12[ijk]-ei12[ijk]-ed12[ijk] );
				et23=( ea[4]+e23[ijk]-ei23[ijk]-ed23[ijk] );
				et13=( ea[5]+e13[ijk]-ei13[ijk]-ed13[ijk] );
				//----- ----- ----- ----- ----- ----- ----- -----
				
				//----- ----- ----- ----- ----- ----- ----- -----
				//Functional derivative of the elastic energy with respect to composition
				/* F=(1/2)*sigma[i][j]*(epsilon[i][j] - epsilon0[i][j])
				   sigma[i][j] = C[i][j][k][l]*(epsilon[k][l] - epsilon0[k][l])
				   epsilon0[i][j] is the position- and composition-dependent eigenstrains */
				/* dF/dc = dF/d(con) = (1/2)*( dCijkl/d(con)*etij*etkl + Cijkl*d(etij)/d(con)*etkl + Cijkl*etij*d(etkl)/d(con) )
				   dCijkl/d(con) = (cpijkl - cmijkl), d(etij)/d(con) = -d(eiij)/d(con) = -ei0 */
				/* delsdc[ijk][0]=0.5*( (cp11-cm11)*et11*et11 -c11[ijk]*ei0*et11 -c11[ijk]*et11*ei0
								    +(cp12-cm12)*et11*et22 -c12[ijk]*ei0*et11 -c12[ijk]*et22*ei0
								    +(cp13-cm13)*et11*et33 -c13[ijk]*ei0*et11 -c13[ijk]*et33*ei0
								    //
								    +(cp21-cm21)*et22*et11 -c21[ijk]*ei0*et22 -c21[ijk]*et11*ei0
								    +(cp22-cm22)*et22*et22 -c22[ijk]*ei0*et22 -c22[ijk]*et22*ei0
								    +(cp23-cm23)*et22*et33 -c23[ijk]*ei0*et22 -c23[ijk]*et33*ei0
								    //
								    +(cp31-cm31)*et33*et11 -c31[ijk]*ei0*et33 -c31[ijk]*et11*ei0
								    +(cp32-cm32)*et33*et22 -c32[ijk]*ei0*et33 -c32[ijk]*et22*ei0
								    +(cp33-cm33)*et33*et33 -c33[ijk]*ei0*et33 -c33[ijk]*et33*ei0
								    //
								    +(cp44-cm44)*et12*et12 -c44[ijk]*ei0*et12 -c44[ijk]*et12*ei0
								    +(cp44-cm44)*et21*et21 -c44[ijk]*ei0*et21 -c44[ijk]*et21*ei0
								    //
								    +(cp55-cm55)*et23*et23 -c55[ijk]*ei0*et23 -c55[ijk]*et23*ei0
								    +(cp55-cm55)*et32*et32 -c55[ijk]*ei0*et32 -c55[ijk]*et32*ei0
								    //
								    +(cp66-cm66)*et13*et13 -c66[ijk]*ei0*et13 -c66[ijk]*et13*ei0
								    +(cp66-cm66)*et31*et31 -c66[ijk]*ei0*et31 -c66[ijk]*et31*ei0
								   ); */
				//
				// c13[ijk]=c12[ijk], c23[ijk]=c13[ijk], c22[ijk]=c11[ijk], c33[ijk]=c11[ijk], etc
				et21=et12;
				et32=et23;
				et31=et13;
				//
				delsdc[ijk]=0.5*( (cp11-cm11)*et11*et11 -c11[ijk]*ei0*( et11 + et11 )
								    +(cp12-cm12)*et11*et22 -c12[ijk]*ei0*( et11 + et22 )
								    +(cp13-cm13)*et11*et33 -c12[ijk]*ei0*( et11 + et33 )
								    //
								    +(cp21-cm21)*et22*et11 -c12[ijk]*ei0*( et22 + et11 )
								    +(cp22-cm22)*et22*et22 -c11[ijk]*ei0*( et22 + et22 )
								    +(cp23-cm23)*et22*et33 -c12[ijk]*ei0*( et22 + et33 )
								    //
								    +(cp31-cm31)*et33*et11 -c12[ijk]*ei0*( et33 + et11 )
								    +(cp32-cm32)*et33*et22 -c12[ijk]*ei0*( et33 + et22 )
								    +(cp33-cm33)*et33*et33 -c11[ijk]*ei0*( et33 + et33 )
								    //
								    +(cp44-cm44)*et12*et12 -c44[ijk]*ei0*( et12 + et12 )
								    +(cp44-cm44)*et21*et21 -c44[ijk]*ei0*( et21 + et21 )
								    //
								    +(cp55-cm55)*et23*et23 -c44[ijk]*ei0*( et23 + et23 )
								    +(cp55-cm55)*et32*et32 -c44[ijk]*ei0*( et32 + et32 )
								    //
								    +(cp66-cm66)*et13*et13 -c44[ijk]*ei0*( et13 + et13 )
								    +(cp66-cm66)*et31*et31 -c44[ijk]*ei0*( et31 + et31 )
								   );
				//delsdc[ijk][1]=0.0;
				//----- ----- ----- ----- ----- ----- ----- -----
			}
		}
	}
	
	//----- ----- ----- ----- -----
	cufftDestroy(plan);
	cufftDestroy(iplan);
	//----- ----- ----- ----- -----
	cudaFree(s11_d);
	cudaFree(s22_d);
	cudaFree(s33_d);
	//
	cudaFree(s12_d);
	cudaFree(s23_d);
	cudaFree(s13_d);
	//----- ----- ----- ----- -----
	cudaFree(e11_d);
	cudaFree(e22_d);
	cudaFree(e33_d);
	//
	cudaFree(e12_d);
	cudaFree(e23_d);
	cudaFree(e13_d);
	//----- ----- ----- ----- -----
	cudaFree(s11k_d);
	cudaFree(s22k_d);
	cudaFree(s33k_d);
	//
	cudaFree(s12k_d);
	cudaFree(s23k_d);
	cudaFree(s13k_d);
	//----- ----- ----- ----- -----
	cudaFree(e11k_d);
	cudaFree(e22k_d);
	cudaFree(e33k_d);
	//
	cudaFree(e12k_d);
	cudaFree(e23k_d);
	cudaFree(e13k_d);
	//----- ----- ----- ----- ----- -----
	cudaFree(ei11_d);
	cudaFree(ei22_d);
	cudaFree(ei33_d);
	//
	cudaFree(ei12_d);
	cudaFree(ei23_d);
	cudaFree(ei13_d);
	//----- ----- ----- ----- ----- -----
	cudaFree(c11_d);
	cudaFree(c12_d);
	cudaFree(c44_d);
	//----- ----- ----- ----- ----- -----
	cudaFree(conr_d);
	//----- ----- ----- ----- ----- -----
	free(s11k);
	free(s22k);
	free(s33k);
	//
	free(s12k);
	free(s23k);
	free(s13k);
	//----- ----- ----- ----- -----
	free(e11k);
	free(e22k);
	free(e33k);
	//
	free(e12k);
	free(e23k);
	free(e13k);
	//----- ----- ----- ----- ----- -----
	//free(smatx_real);
	//free(smatx_imag);
	//free(ematx_real);
	//free(ematx_imag);
	free(smatx);
	free(ematx);
	//----- ----- ----- ----- -----
	free(sum_stress);
	//----- ----- ----- ----- -----
	free(ei11);
	free(ei22);
	free(ei33);
	//
	free(ei12);
	free(ei23);
	free(ei13);
	//----- ----- ----- ----- -----
	free(c11);
	free(c12);
	free(c44);
	//----- ----- ----- ----- -----
	
	return;
}