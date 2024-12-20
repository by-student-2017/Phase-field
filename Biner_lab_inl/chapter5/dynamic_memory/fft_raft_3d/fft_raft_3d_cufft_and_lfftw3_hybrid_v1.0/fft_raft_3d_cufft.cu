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

#include <fftw3.h>
//gcc test.c -lfftw3
//#include <mpi.h> //mpi version
//#include <fftw3-mpi.h> //mpi version

/* Memo for complex type
  "float __complex__ " is old version
  "float _Complex " is new version */
//#include <complex.h>
//#include <cuComplex.h>
#define _Complex_I (1.0iF)
#define I i
#undef i
#undef j
//Ref: http://nalab.mind.meiji.ac.jp/~mk/labo/text/complex-c.pdf (Japanese)

#include <cuda.h>   //GPU
/* #include <cuda.h> or
  #include "cuda_runtime.h"
  #include "device_launch_parameters.h" */
//----- ----- -----
#include <cufft.h> //FFT (GPU)

//typedef float cufftReal;
//typedef cuComplex cufftComplex;

void micro_ch_pre_3d(int Nx, int Ny, int Nz, float c0, float *con);

void prepare_fft_3d(int Nx, int Ny, int Nz, 
	float dx, float dy, float dz,
	float *kx, float *ky, float *kz, 
	float *k2, float *k4);

void green_tensor1_3D(int Nx, int Ny, int Nz,
	float *kx, float *ky, float *kz,
	float cm11, float cm12, float cm44,
	float cp11, float cp12, float cp44,
	float *omeg11, float *omeg22, float *omeg33,
	float *omeg12, float *omeg23, float *omeg13);
//void green_tensor2_3D();

float free_energy_ch_3d(float con_ijk);

void solve_elasticity_3d(int Nx, int Ny, int Nz,
	float *kx, float *ky, float *kz,
	float *omeg11, float *omeg22, float *omeg33,
	float *omeg12, float *omeg23, float *omeg13,
	fftw_complex *s11, fftw_complex *s22, fftw_complex *s33,
	fftw_complex *s12, fftw_complex *s23, fftw_complex *s13,
	fftw_complex *e11, fftw_complex *e22, fftw_complex *e33,
	fftw_complex *e12, fftw_complex *e23, fftw_complex *e13,
	float *ed11, float *ed22, float *ed33,
	float *ed12, float *ed23, float *ed13,
	float cm11, float cm12, float cm44,
	float cp11, float cp12, float cp44,
	float *ea, float ei0,
	float *con, float *delsdc);

void write_vtk_grid_values_3D(int nx, int ny, int nz, 
	float dx, float dy, float dz,
	int istep, float *data1);

int main(){
	clock_t start, end;
	float compute_time;
	
	//get initial wall time
	//(Get initial wall clock time beginning of the execution)
	start = clock();
	
	int times=1;
	
	//simulation cell parameters
	int Nx=256*times; //Number of grid points in the x-direction
	int Ny=256*times; //Number of grid points in the y-direction
	int Nz=  2*times; //Number of grid points in the z-direction
	
	int BSX=8; //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
	int BSY=8; //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
	int BSZ=2; //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
	
	//Total number of grid points in the simulation cell
	int NxNyNz=Nx*Ny*Nz;
	
	//The distance between two grid points in x,y-direction
	float dx=1.0; // [nm] unit ?
	float dy=1.0; // [nm] unit ?
	float dz=1.0; // [nm] unit ?
	
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
	
	/* Note: elastic modulus in this case.
	float c22=c33=c11;
	float c21=c12;
	float c31=c13=c12;
	float c32=c23=c12;
	float c55=c66=c44;
	//
	float cm22=cm33=cm11;
	float cm21=cm12;
	float cm31=cm13=cm12;
	float cm32=cm23=cm12;
	float cm55=cm66=cm44;
	//
	float cp22=cp33=cp11;
	float cp21=cp12;
	float cp31=cp13=cp12;
	float cp32=cp23=cp12;
	float cp55=cp66=cp44; */
	
	//eigen strains
	float ei0=0.01; //Maginitude of eigenstrains
	
	//Applied strains
	float ea[6]; //Magnitude of applied strains
	ea[0]=0.00; //e11
	ea[1]=0.01; //e22
	ea[2]=0.00; //e33
	ea[3]=0.00; //e12
	ea[4]=0.00; //e23
	ea[5]=0.00; //e13
	
	int ijk;
	
	//----- ----- ----- ----- ----- -----
	const int fftsizex = Nx, fftsizey = Ny, fftsizez = Nz;
	const int scale=fftsizex*fftsizey*fftsizez;
	float fftw3d_scale = (float)scale;
	//
	cufftComplex *con_d, *dfdcon_d, *delsdc_d;
	cudaMalloc((void**)&con_d,     sizeof(cufftComplex)*Nx*Ny*Nz);
	cudaMalloc((void**)&dfdcon_d,  sizeof(cufftComplex)*Nx*Ny*Nz);
	cudaMalloc((void**)&delsdc_d,  sizeof(cufftComplex)*Nx*Ny*Nz);
	//
	cufftComplex *conk_d, *dfdconk_d, *delsdck_d;
	cudaMalloc((void**)&conk_d,    sizeof(cufftComplex)*Nx*Ny*Nz);
	cudaMalloc((void**)&dfdconk_d, sizeof(cufftComplex)*Nx*Ny*Nz);
	cudaMalloc((void**)&delsdck_d, sizeof(cufftComplex)*Nx*Ny*Nz);
	//
	//cufftComplex *conc_d;
	//cudaMalloc((void**)&conc_d,    sizeof(cufftComplex)*Nx*Ny*Nz);
	//
	cufftHandle plan, iplan;
	//cufftPlan3d(&plan,  Nx, Ny, Nz, CUFFT_R2C);
	//cufftPlan3d(&iplan, Nx, Ny, Nz, CUFFT_C2R);
	cufftPlan3d(&plan,  Nx, Ny, Nz, CUFFT_C2C);
	cufftPlan3d(&iplan, Nx, Ny, Nz, CUFFT_C2C);
	//----- ----- ----- ----- ----- -----
	
	//----- ----- ----- ----- ----- -----
	//const int fftsizex = Nx, fftsizey = Ny, fftsizez = Nz;
	//const int scale=fftsizex*fftsizey*fftsizez;
	//float fftw3d_scale = (float)scale;
	/* fftw_complex *in, *out; // in[i][0] for real, in[i][1] for imag.
	    in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	   out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	   fftw_plan plan, iplan;
	   plan = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);	//For forward FFT
	  iplan = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);	//For inverse FFT */
	//
	//array
	//fftw_complex *con, *conk;
	// con = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	//conk = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	//fftw_plan plan_con, iplan_conk;
	// plan_con  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, con, conk, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//iplan_conk = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, conk,con,  FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//array
	//fftw_complex *dfdcon, *dfdconk;
	// dfdcon = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	//dfdconk = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	//fftw_plan plan_dfdcon, iplan_dfdconk;
	// plan_dfdcon  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, dfdcon, dfdconk, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//iplan_dfdconk = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, dfdconk, dfdcon, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//array
	//fftw_complex *delsdc, *delsdck;
	// delsdc = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	//delsdck = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	//fftw_plan plan_delsdc, iplan_delsdck;
	// plan_delsdc  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, delsdc, delsdck, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//iplan_delsdck = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, delsdck,delsdc,  FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//----- ----- ----- ----- ----- -----
	
	//----- ----- ----- -----
	//stress (head s series) components
	fftw_complex *s11, *s22, *s33;
	 s11 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 s22 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 s33 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 //
	fftw_complex *s12, *s23, *s13;
	 s12 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 s23 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 s13 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	//----- ----- ----- -----
	//strain (head e series) components
	fftw_complex *e11, *e22, *e33;
	 e11 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 e22 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 e33 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 //
	fftw_complex *e12, *e23, *e13;
	 e12 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 e23 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 e13 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	//----- ----- ----- ----- ----- -----
	
	//----- ----- ----- ----- ----- -----
	//float ed11[Nx][Ny][Nz];
	float *ed11 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//float ed22[Nx][Ny][Nz];
	float *ed22 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//float ed33[Nx][Ny][Nz];
	float *ed33 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//
	//float ed12[Nx][Ny][Nz];
	float *ed12 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//float ed23[Nx][Ny][Nz];
	float *ed23 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//float ed13[Nx][Ny][Nz];
	float *ed13 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//----- ----- ----- ----- ----- -----
	
	//Initialize stress & strain componentes
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				ijk=(i*Ny+j)*Nz+k;
				//----- ----- ----- ----- -----
				// stress (smatx, sXYk)
				s11[ijk][0] = s11[ijk][1] = 0.0;
				s22[ijk][0] = s22[ijk][1] = 0.0;
				s33[ijk][0] = s33[ijk][1] = 0.0;
				//
				s12[ijk][0] = s12[ijk][1] = 0.0;
				s23[ijk][0] = s23[ijk][1] = 0.0;
				s13[ijk][0] = s13[ijk][1] = 0.0;
				//----- ----- ----- ----- -----
				// strain (ematx, eXYk)
				e11[ijk][0] = e11[ijk][1] = 0.0;
				e22[ijk][0] = e22[ijk][1] = 0.0;
				e33[ijk][0] = e33[ijk][1] = 0.0;
				//
				e12[ijk][0] = e12[ijk][1] = 0.0;
				e23[ijk][0] = e23[ijk][1] = 0.0;
				e13[ijk][0] = e13[ijk][1] = 0.0;
				//----- ----- ----- ----- -----
				//----- ----- ----- ----- -----
				//Strain components due to lattice defects
				ed11[ijk] = 0.0;
				ed22[ijk] = 0.0;
				ed33[ijk] = 0.0;
				//
				ed12[ijk] = 0.0;
				ed23[ijk] = 0.0;
				ed13[ijk] = 0.0;
				//----- ----- ----- ----- -----
			}
		}
	}
	
	//----- prepare microstructure
	float *con = (float *)malloc(sizeof(float)*( NxNyNz ));
	//
	micro_ch_pre_3d(Nx,Ny,Nz,c0,con); //Initialize microstructure
	//----- ----- ----- -----
	
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
	
	//float omeg11[Nx][Ny][Nz], Coefficient needed for the Green's tensor
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	float *omeg11 = (float *)malloc(sizeof(float)*( NxNyNz ));
	float *omeg22 = (float *)malloc(sizeof(float)*( NxNyNz ));
	float *omeg33 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//
	float *omeg12 = (float *)malloc(sizeof(float)*( NxNyNz ));
	float *omeg23 = (float *)malloc(sizeof(float)*( NxNyNz ));
	float *omeg13 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	
	//Greens tensor
	green_tensor1_3D(Nx,Ny,Nz,
					 kx,ky,kz,
					 cm11,cm12,cm44,
					 cp11,cp12,cp44,
					 omeg11,omeg22,omeg33,
					 omeg12,omeg23,omeg13);
	
	//----- -----
	//float *dfdcon = (float *)malloc(sizeof(float)*( NxNyNz ));
	float *delsdc = (float *)malloc(sizeof(float)*( NxNyNz ));
	//
	float _Complex *conc    = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *dfdconc = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *delsdcc = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	//
	float _Complex *conk    = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *dfdconk = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	float _Complex *delsdck = (float _Complex *)malloc(sizeof(float _Complex)*( NxNyNz ));
	//----- -----
	
	//----- -----
	// semi-implicit scheme
	float denom;
	float _Complex numer;
	//----- -----
	
	int bsx=BSX, bsy=BSY, bsz=BSZ;     //Number of threads
	dim3 blocks(Nx/bsx,Ny/bsy,Nz/bsz); //nx*ny*nz = blocks * threads
	dim3 threads(bsx,bsy,bsz);         //bsx*bsy*bsz <= 1024
	
	//evolve (Time evolution of microstructure)
	for(int istep=0;istep<=nstep;istep++){
		
		//Update the total time
		ttime=ttime+dtime;
		
		//derivative of free energy
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					//Calculate derivative of free energy
					dfdconc[ijk]=free_energy_ch_3d(con[ijk]) + 0.0*I;
				}
			}
		}
		cudaMemcpy(dfdcon_d,dfdconc,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyHostToDevice); //dfdconc = dfdconc_h
		
		//derivative of elastic energy
		//Calculate the derivative of elastic energy
		solve_elasticity_3d(Nx,Ny,Nz,
			kx,ky,kz,
			omeg11,omeg22,omeg33,
			omeg12,omeg23,omeg13,
			s11,s22,s33,
			s12,s23,s13,
			e11,e22,e33,
			e12,e23,e13,
			ed11,ed22,ed33,
			ed12,ed23,ed13,
			cm11,cm12,cm44,
			cp11,cp12,cp44,
			ea,ei0,
			con, delsdc);
		// Note: tmatx is real part only
		//----- ----- ----- -----
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					//replace cuda array with host array
					delsdcc[ijk] = delsdc[ijk];
					conc[ijk] = con[ijk];
				}
			}
		}
		cudaMemcpy(delsdc_d,delsdcc,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyHostToDevice); //delsdc = delsdc_h
		cudaMemcpy(con_d,conc,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyHostToDevice); //con = con_h
		
		/* Take the values of concentration, derivative of free energy and
		   derivative of elastic energy from real space to Fourier space (forward FFT) */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//conk=fft3(con);           //Octave or Matlab
		//----- ----- ----- -----
		//fftw_execute(plan_con);   //fftw3
		//----- ----- ----- -----
		cufftExecC2C(plan, con_d, conk_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//dfdconk=fft3(dfdcon);     //Octave or Matlab
		//----- ----- ----- -----
		//fftw_execute(plan_dfdcon);//fftw3
		//----- ----- ----- -----
		cufftExecC2C(plan, dfdcon_d, dfdconk_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//delsdck=fft3(delsdc);     //Octave or Matlab
		//----- ----- ----- -----
		//fftw_execute(plan_delsdc);//fftw3
		//----- ----- ----- -----
		cufftExecC2C(plan, delsdc_d, delsdck_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		cudaMemcpy(conk,conk_d,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //conk = conk_h
		cudaMemcpy(dfdconk,dfdconk_d,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //dfdconk = dfdconk_h
		cudaMemcpy(delsdck,delsdck_d,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //delsdck = delsdck_h
		
		/* Semi-implicit time integration of concentration field at
		   Fourier space (Eq.5.50) */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- from fftw3
		//for(int i=0;i<Nx;i++){
		//	for(int j=0;j<Ny;j++){
		//		for(int k=0;k<Nz;k++){
		//			ijk=(i*Ny+j)*Nz+k;
		//			//
		//			denom=1.0+dtime*coefA*mobility*grad_coef*k4[ijk]; // real part only
		//			//
		//			numer=dtime*mobility*k2[ijk]*(dfdconk[ijk][0]+delsdck[ijk][0]);
		//			conk[ijk][0]=(conk[ijk][0]-numer)/denom;
		//			//
		//			numer=dtime*mobility*k2[ijk]*(dfdconk[ijk][1]+delsdck[ijk][1]);
		//			conk[ijk][1]=(conk[ijk][1]-numer)/denom;
		//		}
		//	}
		//}
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- from cufft
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					//
					denom=1.0+dtime*coefA*mobility*grad_coef*k4[ijk];
					numer=dtime*mobility*k2[ijk]*(dfdconk[ijk]+delsdck[ijk]);
					conk[ijk]=(conk[ijk]-numer)/denom;
				}
			}
		}
		cudaMemcpy(conk_d,conk,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyHostToDevice); //conk = conk_h
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		/* Take concentration field from Fourier space back to
		   real space (inverse FFT) */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//con=real(ifft3(conk));    //Octave or Matlab
		//----- ----- ----- -----
		//fftw_execute(iplan_conk); //fftw3
		//----- ----- ----- -----
		cufftExecC2C(iplan, conk_d, con_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		cudaMemcpy(conc,con_d,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //conc = conc_h
		
		//for small deviations
		// For small deviations from max and min values, reset the limits
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					//----- ----- ----- -----
					// con[ijk][0] /= fftw3d_scale; //For fftw3
					// con[ijk][1] /= fftw3d_scale; //For fftw3
					con[ijk] = ( __real__ conc[ijk] )/(Nx*Ny*Nz);
					//con[ijk] =  creal(conc[ijk])/(Nx*Ny*Nz); //For #include <complex.h>
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
			printf("done step: %5d \n",istep);
			
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_3D(Nx,Ny,Nz,dx,dy,dz,istep,con);
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
	//----- ----- ----- ----- -----
	fftw_free(s11);
	fftw_free(s22);
	fftw_free(s33);
	//
	fftw_free(s12);
	fftw_free(s23);
	fftw_free(s13);
	//----- ----- ----- ----- -----
	fftw_free(e11);
	fftw_free(e22);
	fftw_free(e33);
	//
	fftw_free(e12);
	fftw_free(e23);
	fftw_free(e13);
	//----- ----- ----- ----- ----- -----
	free(ed11);
	free(ed22);
	free(ed33);
	//
	free(ed12);
	free(ed23);
	free(ed13);
	//----- ----- ----- ----- -----
	free(kx);
	free(ky);
	free(k2);
	free(k4);
	//----- ----- ----- ----- -----
	free(con);
	//----- ----- ----- ----- -----
	free(conc);
	free(dfdconc);
	free(delsdcc);
	//
	free(conk);
	free(dfdconk);
	free(delsdck);
	//----- ----- ----- ----- ----- -----
}
