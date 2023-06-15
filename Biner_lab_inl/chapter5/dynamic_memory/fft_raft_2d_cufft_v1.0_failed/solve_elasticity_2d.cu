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
  cm11: C11 component of elasticity matrix for matrix material
  cm12: C12 component of elasticity matrix for matrix material
  cm44: C44 component of elasticity matrix for matrix material
  cp11: C11 component of elasticity matrix for second phase
  cp12: C12 component of elasticity matrix for second phase
  cp44: C44 component of elasticity matrix for second phase
  ed11: Strain component of lattice defects
  ed22: Strain component of lattice defects
  ed12: Strain component of lattice defects
  ei0: Magnitude of eigenstrains
  ea[3]: Applied strains
  con[Nx][Ny]: Concentration
  s11[Nx][Ny]: Component of stress
  s22[Nx][Ny]: Component of stress
  s12[Nx][Ny]: Component of stress
  e11[Nx][Ny]: Component of strain
  e22[Nx][Ny]: Component of strain
  e12[Nx][Ny]: Component of strain
  delsdc[Nx][Ny]: Functional derivative of elastic energy
  tmatx[Nx][Ny][2][2][2][2]: Values of Green's tensor at all grid points (real part only)
*/

void solve_elasticity_2d(int Nx, int Ny,
	float *tmatx,
	float _Complex *s11, float _Complex *s22, float _Complex *s12,
	float _Complex *e11, float _Complex *e22, float _Complex *e12,
	float *ed11, float *ed22, float *ed12,
	float cm11, float cm12, float cm44,
	float cp11, float cp12, float cp44,
	float *ea,
	float ei0,
	float *con,  float _Complex *delsdc){
	
	//----- ----- ----- -----
	cufftComplex *s11_d, *s22_d, *s12_d;
	cudaMalloc((void**)&s11_d, sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&s22_d, sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&s12_d, sizeof(cufftComplex)*Nx*Ny);
	//----- ----- ----- -----
	cufftComplex *e11_d, *e22_d, *e12_d;
	cudaMalloc((void**)&e11_d, sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&e22_d, sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&e12_d, sizeof(cufftComplex)*Nx*Ny);
	//----- ----- ----- -----
	
	//----- ----- ----- -----
	cufftComplex *s11k_d, *s22k_d, *s12k_d;
	cudaMalloc((void**)&s11k_d, sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&s22k_d, sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&s12k_d, sizeof(cufftComplex)*Nx*Ny);
	//----- ----- ----- -----
	cufftComplex *e11k_d, *e22k_d, *e12k_d;
	cudaMalloc((void**)&e11k_d, sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&e22k_d, sizeof(cufftComplex)*Nx*Ny);
	cudaMalloc((void**)&e12k_d, sizeof(cufftComplex)*Nx*Ny);
	//----- ----- ----- -----
	cufftHandle plan, iplan;
	cufftPlan2d(&plan,  Nx, Ny, CUFFT_C2C);
	cufftPlan2d(&iplan, Nx, Ny, CUFFT_C2C);
	//----- ----- ----- -----
	
	float _Complex *s11k = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *s22k = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *s12k = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	//
	float _Complex *e11k = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *e22k = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	float _Complex *e12k = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny ));
	
	//----- ----- ----- -----
	//const int fftsizex = Nx, fftsizey = Ny;
	
	//fftw_complex *s11k, *s22k, *s12k;
	//s11k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	//s22k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	//s12k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	
	//fftw_plan plan_s11, plan_s22, plan_s12;
	// plan_s11  = fftw_plan_dft_2d(fftsizex, fftsizey, s11, s11k, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	// plan_s22  = fftw_plan_dft_2d(fftsizex, fftsizey, s22, s22k, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	// plan_s12  = fftw_plan_dft_2d(fftsizex, fftsizey, s12, s12k, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//----- ----- ----- -----
	//fftw_complex *e11k, *e22k, *e12k;
	//e11k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	//e22k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	//e12k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	
	//fftw_plan plan_e11, iplan_e11k;
	// plan_e11  = fftw_plan_dft_2d(fftsizex, fftsizey, e11, e11k, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//iplan_e11k = fftw_plan_dft_2d(fftsizex, fftsizey, e11k, e11, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//fftw_plan plan_e22, iplan_e22k;
	// plan_e22  = fftw_plan_dft_2d(fftsizex, fftsizey, e22, e22k, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//iplan_e22k = fftw_plan_dft_2d(fftsizex, fftsizey, e22k, e22, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//fftw_plan plan_e12, iplan_e12k;
	// plan_e12  = fftw_plan_dft_2d(fftsizex, fftsizey, e12, e12k, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//iplan_e12k = fftw_plan_dft_2d(fftsizex, fftsizey, e12k, e12, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//----- ----- ----- -----
	
	//float smatx_real[Nx][Ny][2][2];
	//float *smatx_real = (float *)malloc(sizeof(float)*( Nx*Ny*2*2 ));
	//float ematx_real[Nx][Ny][2][2];
	//float *ematx_real = (float *)malloc(sizeof(float)*( Nx*Ny*2*2 ));
	//
	//float smatx_imag[Nx][Ny][2][2];
	//float *smatx_imag = (float *)malloc(sizeof(float)*( Nx*Ny*2*2 ));
	//float ematx_imag[Nx][Ny][2][2];
	//float *ematx_imag = (float *)malloc(sizeof(float)*( Nx*Ny*2*2 ));
	//
	float _Complex *smatx = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny*2*2 ));
	float _Complex *ematx = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny*2*2 ));
	
	//float sum_stress[Nx][Ny];
	float *sum_stress = (float *)malloc(sizeof(float)*( Nx*Ny ));
	
	float old_norm=0.0;
	float normF=0.0;
	
	float conver=0.0;
	
	float _Complex et11 = 0.0;
	float _Complex et22 = 0.0;
	float _Complex et12 = 0.0;
	
	int ii, ij;
	
	//Maximum number of iteration steps
	int niter=10;
	
	//Tolerance value of convergence tests
	float tolerance=0.001;
	
	//float ei11[Nx][Ny];
	float *ei11 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	//float ei22[Nx][Ny];
	float *ei22 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	//float ei12[Nx][Ny];
	float *ei12 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	
	//float c11[Nx][Ny];
	float  *c11 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	//float c12[Nx][Ny];
	float  *c12 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	//float c44[Nx][Ny];
	float  *c44 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			ii=i*Ny+j;
			
			//Calculate the eigenstrains
			ei11[ii]=ei0*con[ii];
			ei22[ii]=ei0*con[ii];
			ei12[ii]=0.0*con[ii];
			
			/* Calculate the effective elastic constants at 
			   the grid points based on the composition and
			   using Vegard's law */
			c11[ii]=con[ii]*cp11+(1.0-con[ii])*cm11;
			c12[ii]=con[ii]*cp12+(1.0-con[ii])*cm12;
			c44[ii]=con[ii]*cp44+(1.0-con[ii])*cm44;
		}
	}
	
	/* Solve stress and strain field with 
	   iterative algorithm given in the text */
	for(int iter=0;iter<niter;iter++){
		
		cudaMemcpy(s11_d,s11,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //s11 = s11_h
		cudaMemcpy(s22_d,s22,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //s22 = s22_h
		cudaMemcpy(s12_d,s12,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //s12 = s12_h
		//
		cudaMemcpy(e11_d,e11,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //e11 = e11_h
		cudaMemcpy(e22_d,e22,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //e22 = e22_h
		cudaMemcpy(e12_d,e12,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //e12 = e12_h
		
		/* Take stress and strain components from real space to
		   Fourier space (forward FFT). Step-a */
		// stress
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//s11k=fft2(s11);         //Octave or Matlab
		//----- ----- ----- -----
		//fftw_execute(plan_s11); //fftw3
		//----- 
		cufftExecC2C(plan, s11_d, s11k_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//s22k=fft2(s22);         //Octave or Matlab
		//----- ----- ----- -----
		//fftw_execute(plan_s22); //fftw3
		//----- 
		cufftExecC2C(plan, s22_d, s22k_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//s12k=fft2(s12);         //Octave or Matlab
		//----- ----- ----- -----
		//fftw_execute(plan_s12); //fftw3
		//----- 
		cufftExecC2C(plan, s12_d, s12k_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//
		// strain
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//e11k=fft2(e11);         //Octave or Matlab
		//----- ----- ----- -----
		//fftw_execute(plan_e11); //fftw3
		//----- 
		cufftExecC2C(plan, e11_d, e11k_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//e22k=fft2(e22);         //Octave or Matlab
		//----- ----- ----- -----
		//fftw_execute(plan_e22); //fftw3
		//----- 
		cufftExecC2C(plan, e22_d, e22k_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//e12k=fft2(e12);         //Octave or Matlab
		//----- ----- ----- -----
		//fftw_execute(plan_e12); //fftw3
		//----- 
		cufftExecC2C(plan, e12_d, e12k_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		cudaMemcpy(s11k,s11k_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //s11k = s11k_h
		cudaMemcpy(s22k,s22k_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //s22k = s22k_h
		cudaMemcpy(s12k,s12k_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //s12k = s12k_h
		//
		cudaMemcpy(e11k,e11k_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //e11k = e11k_h
		cudaMemcpy(e22k,e22k_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //e22k = e22k_h
		cudaMemcpy(e12k,e12k_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //e12k = e12k_h
		
		/* Form stress and strain tensors to be used in 
		   Eq.5.46, Step-b */
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				//smatx_real[(ii*2+0)*2+0] = __real__ s11k[ii];
				//smatx_real[(ii*2+0)*2+1] = __real__ s12k[ii];
				//smatx_real[(ii*2+1)*2+0] = __real__ s12k[ii];
				//smatx_real[(ii*2+1)*2+1] = __real__ s22k[ii];
				//
				//smatx_imag[(ii*2+0)*2+0] = __imag__ s11k[ii];
				//smatx_imag[(ii*2+0)*2+1] = __imag__ s12k[ii];
				//smatx_imag[(ii*2+1)*2+0] = __imag__ s12k[ii];
				//smatx_imag[(ii*2+1)*2+1] = __imag__ s22k[ii];
				//
				//ematx_real[(ii*2+0)*2+0] = __real__ e11k[ii];
				//ematx_real[(ii*2+0)*2+1] = __real__ e12k[ii];
				//ematx_real[(ii*2+1)*2+0] = __real__ e12k[ii];
				//ematx_real[(ii*2+1)*2+1] = __real__ e22k[ii];
				//
				//ematx_imag[(ii*2+0)*2+0] = __imag__ e11k[ii];
				//ematx_imag[(ii*2+0)*2+1] = __imag__ e12k[ii];
				//ematx_imag[(ii*2+1)*2+0] = __imag__ e12k[ii];
				//ematx_imag[(ii*2+1)*2+1] = __imag__ e22k[ii];
				//
				smatx[(ii*2+0)*2+0] = s11k[ii];
				smatx[(ii*2+0)*2+1] = s12k[ii];
				smatx[(ii*2+1)*2+0] = s12k[ii];
				smatx[(ii*2+1)*2+1] = s22k[ii];
				//
				ematx[(ii*2+0)*2+0] = e11k[ii];
				ematx[(ii*2+0)*2+1] = e12k[ii];
				ematx[(ii*2+1)*2+0] = e12k[ii];
				ematx[(ii*2+1)*2+1] = e22k[ii];
			}
		}
		
		//Green operator
		// Calculate strain tensor, Eq.5.46, Step-b
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ij=i*Ny+j;
				//
				for(int kk=0;kk<2;kk++){
					for(int ll=0;ll<2;ll++){
						for(int ii=0;ii<2;ii++){
							for(int jj=0;jj<2;jj++){
								/* Eq.5.46(b): new epsilon(zeta) = epsilon(zeta) - sum( gamma(zeta)*sigma(zeta) )
								   where gamma=tmatx, sigma=smatx
								   Note: tmatx is real part only */
								//ematx_real[(ij*2+ii)*2+jj]=ematx_real[(ij*2+ii)*2+jj]
								//	-tmatx[(((ij*2+kk)*2+ll)*2+ii)*2+jj]*smatx_real[(ij*2+kk)*2+ll];
								//
								//ematx_imag[(ij*2+ii)*2+jj]=ematx_imag[(ij*2+ii)*2+jj]
								//	-tmatx[(((ij*2+kk)*2+ll)*2+ii)*2+jj]*smatx_imag[(ij*2+kk)*2+ll];
								ematx[(ij*2+ii)*2+jj]=ematx[(ij*2+ii)*2+jj]
									-tmatx[(((ij*2+kk)*2+ll)*2+ii)*2+jj]*smatx[(ij*2+kk)*2+ll];
							}//jj
						}//ii
					}//ll
				}//kk
				//
			}//Ny
		}//Nx
		
		// Rearrange strain components using symmetry of strain tensor
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				
				//e11k[ii][0]=ematx_real[(ii*2+0)*2+0];
				//e12k[ii][0]=ematx_real[(ii*2+0)*2+1];
				//e12k[ii][0]=ematx_real[(ii*2+1)*2+0];
				//e22k[ii][0]=ematx_real[(ii*2+1)*2+1];
				//
				//e11k[ii][1]=ematx_imag[(ii*2+0)*2+0];
				//e12k[ii][1]=ematx_imag[(ii*2+0)*2+1];
				//e12k[ii][1]=ematx_imag[(ii*2+1)*2+0];
				//e22k[ii][1]=ematx_imag[(ii*2+1)*2+1];
				//
				//e11k[ii] = ematx_real[(ii*2+0)*2+0] + ematx_imag[(ii*2+0)*2+0]*I;
				//e12k[ii] = ematx_real[(ii*2+0)*2+1] + ematx_imag[(ii*2+0)*2+1]*I;
				//e12k[ii] = ematx_real[(ii*2+1)*2+0] + ematx_imag[(ii*2+1)*2+0]*I;
				//e22k[ii] = ematx_real[(ii*2+1)*2+1] + ematx_imag[(ii*2+1)*2+1]*I;
				//
				e11k[ii] = ematx[(ii*2+0)*2+0];
				e12k[ii] = ematx[(ii*2+0)*2+1];
				//e12k[ii] = ematx[(ii*2+1)*2+0];
				e22k[ii] = ematx[(ii*2+1)*2+1];
			}
		}
		
		cudaMemcpy(e11k_d,e11k,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //e11k = e11k_h
		cudaMemcpy(e22k_d,e22k,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //e22k = e22k_h
		cudaMemcpy(e12k_d,e12k,Nx*Ny*sizeof(float _Complex),cudaMemcpyHostToDevice); //e12k = e12k_h
		
		//From Fourier space to real space
		/* Take strain components from Fourier space back to
		   real space (inverse FFT), Step-c */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//e11=real(ifft2(e11k));    //Octave or Matlab
		//----- ----- ----- -----
		//fftw_execute(iplan_e11k); //fftw3
		//-----
		cufftExecC2C(iplan, e11k_d, e11_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//e22=real(ifft2(e22k));    //Octave or Matlab
		//----- ----- ----- -----
		//fftw_execute(iplan_e22k); //fftw3
		//-----
		cufftExecC2C(iplan, e22k_d, e22_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//e12=real(ifft2(e12k));    //Octave or Matlab
		//----- ----- ----- -----
		//fftw_execute(iplan_e12k); //fftw3
		//-----
		cufftExecC2C(iplan, e12k_d, e12_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		cudaMemcpy(e11,e11_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //e11 = e11_h
		cudaMemcpy(e22,e22_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //e22 = e22_h
		cudaMemcpy(e12,e12_d,Nx*Ny*sizeof(float _Complex),cudaMemcpyDeviceToHost); //e12 = e12_h
		
		//Calculate stresses
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				//
				//e11[ii][0]=e11[ii][0]/(fftsizex*fftsizey);
				//e22[ii][0]=e22[ii][0]/(fftsizex*fftsizey);
				//e12[ii][0]=e12[ii][0]/(fftsizex*fftsizey);
				//
				//e11[ii][1]=e11[ii][1]/(fftsizex*fftsizey);
				//e22[ii][1]=e22[ii][1]/(fftsizex*fftsizey);
				//e12[ii][1]=e12[ii][1]/(fftsizex*fftsizey);
				//
				e11[ii] = e11[ii]/(Nx*Ny);
				e22[ii] = e22[ii]/(Nx*Ny);
				e12[ii] = e12[ii]/(Nx*Ny);
				//
				/* s11[ii][0]=c11[ii]*(ea[0]+e11[ii][0]-ei11[ii]-ed11[ii])
						  +c12[ii]*(ea[1]+e22[ii][0]-ei22[ii]-ed22[ii]));
				s22[ii][0]=c21[ii]*(ea[0]+e11[ii][0]-ei11[ii]-ed11[ii])
						  +c22[ii]*(ea[1]+e22[ii][0]-ei22[ii]-ed22[ii]); */
				// c21[ii]=c12[ii], c22[ii]=c11[ii], etc
				//-----
				//s11[ii][0]=c11[ii]*(ea[0]+e11[ii][0]-ei11[ii]-ed11[ii])
				//		  +c12[ii]*(ea[1]+e22[ii][0]-ei22[ii]-ed22[ii]);
				//s22[ii][0]=c12[ii]*(ea[0]+e11[ii][0]-ei11[ii]-ed11[ii])
				//		  +c11[ii]*(ea[1]+e22[ii][0]-ei22[ii]-ed22[ii]);
				//s11[ii][1]=0.0;
				//s22[ii][1]=0.0;
				//
				/* s12[ii][0]=c44[ii]*(ea[2]+e12[ii][0]-ei12[ii]-ed12[ii])
						  +c44[ii]*(ea[2]+e21[ii][0]-ei21[ii]-ed21[ii]); */
				// e21[ii]=e12[ii], etc
				//-----
				//s12[ii][0]=c44[ii]*(ea[2]+e12[ii][0]-ei12[ii]-ed12[ii])*2.0;
				//
				//s12[ii][1]=0.0;
				//
				s11[ii] = c11[ii]*(ea[0]+e11[ii]-ei11[ii]-ed11[ii])
						 +c12[ii]*(ea[1]+e22[ii]-ei22[ii]-ed22[ii]);
				s22[ii] = c12[ii]*(ea[0]+e11[ii]-ei11[ii]-ed11[ii])
						 +c11[ii]*(ea[1]+e22[ii]-ei22[ii]-ed22[ii]);
				s12[ii] = c44[ii]*(ea[2]+e12[ii]-ei12[ii]-ed12[ii])*2.0;
			}
		}
		
		//check convergence
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				sum_stress[ii] = ( __real__ s11[ii] + __real__ s22[ii] + __real__ s12[ii] );
			}
		}
		
		//normF=norm(sum_stress,2.0);
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				normF = normF + sum_stress[ii]*sum_stress[ii];
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
			ii=i*Ny+j;
			
			//Calculate strain components
			et11=ea[0]+e11[ii]-ei11[ii]-ed11[ii];
			et22=ea[1]+e22[ii]-ei22[ii]-ed22[ii];
			et12=ea[2]+e12[ii]-ei12[ii]-ed12[ii];
			
			//Functional derivative of the elastic energy with respect to composition
			/* F=(1/2)*sigma[i][j]*(epsilon[i][j] - epsilon0[i][j])
			   sigma[i][j] = C[i][j][k][l]*(epsilon[k][l] - epsilon0[k][l])
			   epsilon0[i][j] is the position- and composition-dependent eigenstrains */
			/* dF/dc = dF/d(con) = (1/2)*( dCijkl/d(con)*etij*etkl + Cijkl*d(etij)/d(con)*etkl + Cijkl*etij*d(etkl)/d(con) )
			   dCijkl/d(con) = (cpijkl - cmijkl), d(etij)/d(con) = -d(eiij)/d(con) = -ei0 */
			// cp21=cp12, cp22=cp11, et21=et12, etc
			/* delsdc[ii][0]=0.5*( (cp11-cm11)*et11*et11 -c11[ii]*et11*ei0 -c11[ii]*et11*ei0
							   +(cp12-cm12)*et11*et22 -c12[ii]*et11*ei0 -c12[ii]*et22*ei0
							   //
							   +(cp21-cm21)*et22*et11 -c21[ii]*et22*ei0 -c21[ii]*et11*ei0
							   +(cp22-cm22)*et22*et22 -c22[ii]*et22*ei0 -c22[ii]*et22*ei0
							   //
							   +(cp44-cm44)*et12*et12 -c44[ii]*et12*ei0 -c44[ii]*et21*ei0
							   +(cp44-cm44)*et21*et21 -c44[ii]*et21*ei0 -c44[ii]*et12*ei0
							  ); */
			delsdc[ii]=0.5*(et11*( (cp12-cm12)*et22 + (cp11-cm11)*et11 - c12[ii]*ei0 - c11[ii]*ei0 )
							   -ei0*(     c12[ii]*et22 +     c11[ii]*et11 )
						      +et22*( (cp11-cm11)*et22 + (cp12-cm12)*et11 - c12[ii]*ei0 - c11[ii]*ei0 )
							   -ei0*(     c11[ii]*et22 +     c12[ii]*et11 )
							   +2.0*(cp44-cm44)*et12*et12 -4.0*ei0*c44[ii]*et12
							  );
			//delsdc[ii][1]=0.0;
		}
	}
	
	//----- ----- ----- ----- ----- -----
	cufftDestroy(plan);
	cufftDestroy(iplan);
	//----- ----- ----- ----- ----- -----
	cudaFree(s11k_d);
	cudaFree(s22k_d);
	cudaFree(s12k_d);
	//
	cudaFree(e11k_d);
	cudaFree(e22k_d);
	cudaFree(e12k_d);
	//----- ----- ----- ----- ----- -----
	free(s11k);
	free(s22k);
	free(s12k);
	//
	free(e11k);
	free(e22k);
	free(e12k);
	//----- ----- ----- ----- ----- -----
	//free(smatx_real);
	//free(smatx_imag);
	//free(ematx_real);
	//free(ematx_imag);
	free(smatx);
	free(ematx);
	//
	free(sum_stress);
	//
	free(ei11);
	free(ei22);
	free(ei12);
	//
	free(c11);
	free(c12);
	free(c44);
	//----- ----- ----- ----- ----- -----
	
	return;
}