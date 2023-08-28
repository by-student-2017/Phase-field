/* 3D semi-implicit spectral phase-field code 
  for solving precipitation in Fe_Cu_Ni_Mn alloy */
  
// The dimension of energy were normalized with RT.
/* The time t was normalized with dx^2/Dcua(T),
where Dacu(T) is the diffusion constant of Cu in
the alpha phase at temperature T. */

/* The phase-field variable eta(r,t) characterize
  the phases distribution of alpha(bcc) and gamma(fcc) phase in
  the Cu precipitates and takes the values of 0 < h(eta) < 1, 
  corresponding h(eta)=0 for alpha(bcc) and
  h(eta)=1 for gamma(fcc) phase, respectively. */

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

void init_FeCuMnNi_micro_3d(int Nx, int Ny, int Nz, 
	float cu0, float mn0, float ni0,
	float *cu, float *mn, float *ni, float *orp);

void prepare_fft_3d(int Nx, int Ny, int Nz, 
	float dx, float dy, float dz,
	float *kx, float *ky, float *kz, 
	float *k2, float *k4);

void FeCuMnNi_free_energy_3d(int Nx, int Ny, int Nz, 
	float *cu, float *mn, float *ni, float *orp, float tempr,
	float *dgdcu, float *dgdmn, float *dgdni, float *dgdor);

void write_vtk_grid_values_3D(int nx, int ny, int nz, 
	float dx, float dy, float dz, int istep, 
	float *data1, float *data2, float *data3, float *data4);

int main(){
	clock_t start, end;
	float compute_time;
	
	//get initial wall time
	start = clock();
	
	//simulation cell parameters
	int Nx=128;
	int Ny=128;
	int Nz=2;
	
	int BSX=8; //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
	int BSY=8; //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
	int BSZ=2; //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
	
	//Total number of grid points in the simulation cell
	//int NxNy=Nx*Ny;
	
	//The distance between two grid points in x,y-direction
	float dx=0.5; //[nm] unit ?
	float dy=0.5; //[nm] unit ?
	float dz=0.5; //[nm] unit ?
	
	//time integration parameters
	int nstep=5000;      //Number of time steps
	int nprint=50;       //Output frequency to write the results to file
	float dtime=1.0e-2; //Time increment for numerical integration
	float ttime=0.0;    //Total time
	//float coefA=1.0;
	
	//material specific parameters
	
	//Initial concentrations of alloying elements
	// e.g., 15 at.% Cu, 1 at.% Mn, 1 at.% Ni, (100-15-1-1) at.% Fe
	float cu0=0.15;
	float mn0=0.01;
	float ni0=0.01;
	
	//diffusivity parameters for mobility
	float gconst=8.314472; //The value of gas constant [J/mol/K]
	float tempr=823.0;   //Temperature [K]
	float RT=gconst*tempr;
	
	//gradient energy coefficients [J(nm)^2/mol]
	float Kc=5.0e3;   //[J(nm)^2/mol], 5.0e-15 [Jm^2/mol], Cahn-Hilliard eq. (conservation)
	float Keta=1.0e3; //[J(nm)^2/mol], 1.0e-15 [Jm^2/mol], Allen-Cahn eq. (non-conservation)
	
	//gradient energy of the composition
	//Cahn-Hilliard eq. (conservation)
	float grcoef_cu=Kc/RT; //0.68884=5000.0/(8.314472*873) for 873 [K]
	float grcoef_ni=Kc/RT; //0.68884=5000.0/(8.314472*873) for 873 [K]
	float grcoef_mn=Kc/RT; //0.68884=5000.0/(8.314472*873) for 873 [K]
	
	//gradient energy of the phase field
	//Allen-Cahn eq. (non-conservation)
	float grcoef_or=Keta/RT; //0.13729=1000/(8.314472*873) for 873 [K]
	
	// Diffusion constant [m^2/s]
	// A=alpha phase, G=gamma phase
	// D = D0*exp(-Q/RT)
	// D0[m^2/s], Q[J/mol]
	//
	float D0ACu=4.7e-5;
	float D0GCu=4.3e-5;
	//
	float QACu=2.44e5;
	float QGCu=2.80e5;
	//
	float D0ANi=1.4e-4;
	float D0GNi=1.08e-5;
	//
	//QANi=2.46e5;
	float QANi=2.56e5;
	float QGNi=2.74e5;
	//
	float D0AMn=1.49e-4;
	float D0GMn=2.78e-5;
	//
	//QAMn=2.33e5;
	float QAMn=2.64e5;
	float QGMn=2.64e5;
	//
	float DCuA=(D0ACu*exp(-QACu/RT));
	float DCuG=(D0GCu*exp(-QGCu/RT))/DCuA;
	//
	float DNiA=(D0ANi*exp(-QANi/RT))/DCuA;
	float DNiG=(D0GNi*exp(-QGNi/RT))/DCuA;
	//
	float DMnA=(D0AMn*exp(-QAMn/RT))/DCuA;
	float DMnG=(D0GMn*exp(-QGMn/RT))/DCuA;
	//
	DCuA=1.0;
	
	int ii;
	
	//----- ----- ----- ----- ----- -----
	//const int fftsizex = Nx, fftsizey = Ny;
	//----- ----- ----- ----- 
	cufftComplex *cu_d, *cuk_d; // Cu
	cudaMalloc((void**)&cu_d,    sizeof(cufftComplex)*Nx*Ny*Nz);
	cudaMalloc((void**)&cuk_d,   sizeof(cufftComplex)*Nx*Ny*Nz);
	//
	cufftComplex *mn_d, *mnk_d; // Mn
	cudaMalloc((void**)&mn_d,    sizeof(cufftComplex)*Nx*Ny*Nz);
	cudaMalloc((void**)&mnk_d,   sizeof(cufftComplex)*Nx*Ny*Nz);
	//
	cufftComplex *ni_d, *nik_d; // Ni
	cudaMalloc((void**)&ni_d,    sizeof(cufftComplex)*Nx*Ny*Nz);
	cudaMalloc((void**)&nik_d,   sizeof(cufftComplex)*Nx*Ny*Nz);
	//
	cufftComplex *orp_d, *orpk_d; // order parameter
	cudaMalloc((void**)&orp_d,   sizeof(cufftComplex)*Nx*Ny*Nz);
	cudaMalloc((void**)&orpk_d,  sizeof(cufftComplex)*Nx*Ny*Nz);
	//----- ----- ----- ----- 
	cufftComplex *dgdcu_d, *dgdcuk_d; // dG/dCu
	cudaMalloc((void**)&dgdcu_d,  sizeof(cufftComplex)*Nx*Ny*Nz);
	cudaMalloc((void**)&dgdcuk_d, sizeof(cufftComplex)*Nx*Ny*Nz);
	//
	cufftComplex *dgdmn_d, *dgdmnk_d; // dG/dMn
	cudaMalloc((void**)&dgdmn_d,  sizeof(cufftComplex)*Nx*Ny*Nz);
	cudaMalloc((void**)&dgdmnk_d, sizeof(cufftComplex)*Nx*Ny*Nz);
	//
	cufftComplex *dgdni_d, *dgdnik_d; // dG/dNi
	cudaMalloc((void**)&dgdni_d,  sizeof(cufftComplex)*Nx*Ny*Nz);
	cudaMalloc((void**)&dgdnik_d, sizeof(cufftComplex)*Nx*Ny*Nz);
	//
	cufftComplex *dgdor_d, *dgdork_d; // dG/dor
	cudaMalloc((void**)&dgdor_d,  sizeof(cufftComplex)*Nx*Ny*Nz);
	cudaMalloc((void**)&dgdork_d, sizeof(cufftComplex)*Nx*Ny*Nz);
	//----- ----- ----- ----- 
	cufftHandle plan, iplan;
	cufftPlan3d(&plan,  Nx, Ny, Nz, CUFFT_C2C);
	cufftPlan3d(&iplan, Nx, Ny, Nz, CUFFT_C2C);
	//----- ----- ----- ----- ----- -----
	
	//float kx[Nx];
	float *kx = (float *)malloc(sizeof(float)*( Nx ));
	//float ky[Ny];
	float *ky = (float *)malloc(sizeof(float)*( Ny ));
	//float kz[Ny];
	float *kz = (float *)malloc(sizeof(float)*( Nz ));
	//float k2[Nx][Ny][Nz];
	float *k2 = (float *)malloc(sizeof(float)*( Nx*Ny*Nz ));
	//float k4[Nx][Ny][Nz];
	float *k4 = (float *)malloc(sizeof(float)*( Nx*Ny*Nz ));
	
	//prepare fft (output: kx,ky,kz,k2,k4)
	prepare_fft_3d(Nx,Ny,Nz,dx,dy,dz,kx,ky,kz,k2,k4); //get FFT coefficients
	
	//float  cu[Nx][Ny];
	 float *cu = (float *)malloc(sizeof(float)*( Nx*Ny*Nz ));
	//float  mn[Nx][Ny];
	 float *mn = (float *)malloc(sizeof(float)*( Nx*Ny*Nz ));
	//float  ni[Nx][Ny];
	 float *ni = (float *)malloc(sizeof(float)*( Nx*Ny*Nz ));
	//float orp[Nx][Ny];
	float *orp = (float *)malloc(sizeof(float)*( Nx*Ny*Nz ));
	
	//prepare microstructure
	init_FeCuMnNi_micro_3d(Nx,Ny,Nz,cu0,mn0,ni0,cu,mn,ni,orp);
	
	float _Complex  *cuc   = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny*Nz )); // Cu
	float _Complex  *mnc   = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny*Nz )); // Mn
	float _Complex  *nic   = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny*Nz )); // Ni
	float _Complex *orpc   = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny*Nz )); // Order parameter
	
	float _Complex  *cuk   = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny*Nz )); // Cu
	float _Complex  *mnk   = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny*Nz )); // Mn
	float _Complex  *nik   = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny*Nz )); // Ni
	float _Complex *orpk   = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny*Nz )); // Order parameter
	
	//float dgdcu[Nx][Ny][Nz];
	 float *dgdcu = (float *)malloc(sizeof(float)*( Nx*Ny*Nz ));
	//float dgdmn[Nx][Ny][Nz];
	 float *dgdmn = (float *)malloc(sizeof(float)*( Nx*Ny*Nz ));
	//float dgdni[Nx][Ny][Nz];
	 float *dgdni = (float *)malloc(sizeof(float)*( Nx*Ny*Nz ));
	//float dgdor[Nx][Ny][Nz];
	 float *dgdor = (float *)malloc(sizeof(float)*( Nx*Ny*Nz ));
	
	float _Complex *dgdcuc = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny*Nz ));
	float _Complex *dgdmnc = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny*Nz ));
	float _Complex *dgdnic = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny*Nz ));
	float _Complex *dgdorc = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny*Nz ));
	
	float _Complex *dgdcuk = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny*Nz )); // dG/dCu
	float _Complex *dgdmnk = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny*Nz )); // dG/dMn
	float _Complex *dgdnik = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny*Nz )); // dG/dNi
	float _Complex *dgdork = (float _Complex *)malloc(sizeof(float _Complex)*( Nx*Ny*Nz )); // dG/d(Order parameter)
	
	float mcoef_cu;
	float mcoef_mn;
	float mcoef_ni;
	float mcoef_orp;
	
	int bsx=BSX, bsy=BSY, bsz=BSZ;     //Number of threads
	dim3 blocks(Nx/bsx,Ny/bsy,Nz/bsz); //nx*ny*nz = blocks * threads
	dim3 threads(bsx,bsy,bsz);         //bsx*bsy*bsz <= 1024
	
	//evolve (Evolve the microstructure)
	for(int istep=0;istep<=nstep;istep++){
		
		//Update the total time
		ttime=ttime+dtime;
		
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					//ii=i*Ny*Nz+j*Nz+k;
					ii=(i*Ny+j)*Nz+k;
					cuc[ii] = cu[ii];
					mnc[ii] = mn[ii];
					nic[ii] = ni[ii];
					orpc[ii] = orp[ii];
				}
			}
		}
		//
		cudaMemcpy(cu_d,cuc,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyHostToDevice); //cuc = cuc_h, Cu
		cudaMemcpy(mn_d,mnc,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyHostToDevice); //mnc = mnc_h, Mn
		cudaMemcpy(ni_d,nic,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyHostToDevice); //nic = nic_h, Ni
		cudaMemcpy(orp_d,orpc,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyHostToDevice); //orpc = orpc_h, Order parameter
		//
		/* Transform the alloying elements and
		   the order parameter values from
		   real space to Fourier space (forward FFT transformation) */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, cu_d, cuk_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, mn_d, mnk_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, ni_d, nik_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, orp_d, orpk_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//
		cudaMemcpy(cuk,cuk_d,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //cu = cu_h, Cu
		cudaMemcpy(mnk,mnk_d,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //mn = mn_h, Mn
		cudaMemcpy(nik,nik_d,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //ni = ni_h, Ni
		cudaMemcpy(orpk,orpk_d,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //orp = orp_h, Order
		
		//derivative of free energy
		/* Calculate the derivatives of free energy for
		   each alloying elements and the order parameter */
		//Note: input array (real): cu, mn, ni, orp
		FeCuMnNi_free_energy_3d(Nx,Ny,Nz,cu,mn,ni,orp,tempr,dgdcu,dgdmn,dgdni,dgdor);
		//Note: output array (real): dgdcu, dgdmn, dgdni, dgdor
		
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					//ii=i*Ny*Nz+j*Nz+k;
					ii=(i*Ny+j)*Nz+k;
					//----- ----- ----- -----
					dgdcuc[ii] = dgdcu[ii];
					//----- ----- ----- -----
					dgdmnc[ii] = dgdmn[ii];
					//----- ----- ----- -----
					dgdnic[ii] = dgdni[ii];
					//----- ----- ----- -----
					dgdorc[ii] = dgdor[ii];
					//----- ----- ----- -----
				}
			}
		}
		//
		cudaMemcpy(dgdcu_d,dgdcuc,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyHostToDevice); //dgdcuc = dgdcuc_h, dG/dCu
		cudaMemcpy(dgdmn_d,dgdmnc,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyHostToDevice); //dgdmnc = dgdmnc_h, dG/dMn
		cudaMemcpy(dgdni_d,dgdnic,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyHostToDevice); //dgdnic = dgdnic_h, dG/dNi
		cudaMemcpy(dgdor_d,dgdorc,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyHostToDevice); //dgdorc = dgdorc_h, dG/d(Order parameter)
		//
		/* Transform the derivative values from
		   real space to Fourier space (forward FFT transformation) */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, dgdcu_d, dgdcuk_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, dgdmn_d, dgdmnk_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, dgdni_d, dgdnik_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(plan, dgdor_d, dgdork_d, CUFFT_FORWARD); //FFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//
		cudaMemcpy(dgdcuk,dgdcuk_d,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //dgdcu = dgdcu_h, dG/dCu
		cudaMemcpy(dgdmnk,dgdmnk_d,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //dgdmn = dgdmn_h, dG/dMn
		cudaMemcpy(dgdnik,dgdnik_d,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //dgdni = dgdni_h, dG/dNi
		cudaMemcpy(dgdork,dgdork_d,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //dgdor = dgdor_h, dG/d(Order parameter)
		
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					//ii=i*Ny*Nz+j*Nz+k;
					ii=(i*Ny+j)*Nz+k;
					
					//mobilities
					/* Calculate the mobility of each alloying elements,
					   Eq.5.24, based on the current values of
					   concentration and order parameter */
					/* Mi(eta,T) = ci0*(1.0-ci0)*{(1.0-eta)*Dia(T)/RT + eta*Dig(T)/RT} (Eq.5.24)
					   eta=orp[i][j] */
					mcoef_cu=cu0*(1.0-cu0)*( (1.0-orp[ii])*DCuA + orp[ii]*DCuG );
					mcoef_mn=mn0*(1.0-mn0)*( (1.0-orp[ii])*DMnA + orp[ii]*DMnG );
					mcoef_ni=ni0*(1.0-ni0)*( (1.0-orp[ii])*DNiA + orp[ii]*DNiG );
					//
					mcoef_orp=0.1;
					
					//time integration
					/* Semi-implicit time integration of each alloying elements
					   and the order parameter */
					// denominator is a real number
					// These are related by Eq.5.14 (Cahn-Hilliard for ci) or Eq.5.21 (Allen-Cahn for eta)
					//----- ----- ----- -----
					 cuk[ii]= (cuk[ii]-dtime*k2[ii]*mcoef_cu*dgdcuk[ii])/
						(1.0+dtime*k4[ii]*mcoef_cu*grcoef_cu);
					 mnk[ii]= (mnk[ii]-dtime*k2[ii]*mcoef_mn*dgdmnk[ii])/
						(1.0+dtime*k4[ii]*mcoef_mn*grcoef_mn);
					 nik[ii]= (nik[ii]-dtime*k2[ii]*mcoef_ni*dgdnik[ii])/
						(1.0+dtime*k4[ii]*mcoef_ni*grcoef_ni);
					//
					orpk[ii]=(orpk[ii]-dtime*mcoef_orp*dgdork[ii])/
						(1.0+dtime*k2[ii]*mcoef_orp*grcoef_or);
					//----- ----- ----- -----
				}
			}
		}
		
		cudaMemcpy(cuk_d,cuk,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyHostToDevice); //cuk = cuk_h, Cu
		cudaMemcpy(mnk_d,mnk,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyHostToDevice); //mnk = mnk_h, Mn
		cudaMemcpy(nik_d,nik,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyHostToDevice); //nik = nik_h, Ni
		cudaMemcpy(orpk_d,orpk,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyHostToDevice); //orpk = orpk_h, Order
		
		/* Bring back the integrated values from
		Fourier space to real space (inverse FFT transformations) */
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(iplan, cuk_d, cu_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(iplan, mnk_d, mn_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(iplan, nik_d, ni_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		cufftExecC2C(iplan, orpk_d, orp_d, CUFFT_INVERSE); //IFFT
		cudaDeviceSynchronize();
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		cudaMemcpy(cuc,cu_d,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //cuc = cuc_h, Cu
		cudaMemcpy(mnc,mn_d,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //mnc = mnc_h, Mn
		cudaMemcpy(nic,ni_d,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //nic = nic_h, Ni
		cudaMemcpy(orpc,orp_d,Nx*Ny*Nz*sizeof(float _Complex),cudaMemcpyDeviceToHost); //orpc = orpc_h, Order parameter
		
		//for small deviations
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					//ii=i*Ny*Nz+j*Nz+k;
					ii=(i*Ny+j)*Nz+k;
					//----- ----- ----- -----
					cu[ii] = ( __real__ cuc[ii] )/(Nx*Ny);
					//cu[ii] =  creal(cuc[ii])/(Nx*Ny); //For #include <_Complex.h>
					//----- -----
					mn[ii] = ( __real__ mnc[ii] )/(Nx*Ny);
					//mn[ii] =  creal(mnc[ii])/(Nx*Ny); //For #include <_Complex.h>
					//----- -----
					ni[ii] = ( __real__ nic[ii] )/(Nx*Ny);
					//ni[ii] =  creal(nic[ii])/(Nx*Ny); //For #include <_Complex.h>
					//----- -----
					orp[ii] = ( __real__ orpc[ii] )/(Nx*Ny);
					//orp[ii] =  creal(orpc[ii])/(Nx*Ny); //For #include <_Complex.h>
					//----- ----- ----- -----
					//Cu
					if(cu[ii]>=0.9999){
					cu[ii]=0.9999;
					}
					if(cu[ii]<=0.0001){
						cu[ii]=0.0001;
					}
					//----- ----- ----- -----
					//Mn
					if(mn[ii]>=0.9999){
						mn[ii]=0.9999;
					}
					if(mn[ii]<=0.0001){
						mn[ii]=0.0001;
					}
					//----- ----- ----- -----
					//Ni
					if(ni[ii]>=0.9999){
						ni[ii]=0.9999;
					}
					if(ni[ii]<=0.0001){
						ni[ii]=0.0001;
					}
					//----- ----- ----- -----
				}
			}
		}
		
		//print results
		/* If print frequency is reached, output the results to file */
		if(fmod(istep,nprint)==0){
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_3D(Nx,Ny,Nz,dx,dy,dz,istep,cu,mn,ni,orp);
			
			printf("done step: %5d, time %e [s], Temp: %8.3lf [K] \n",istep,ttime,tempr);
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
	cudaFree(cu_d);
	cudaFree(cuk_d);
	//
	cudaFree(mn_d);
	cudaFree(mnk_d);
	//
	cudaFree(ni_d);
	cudaFree(nik_d);
	//
	cudaFree(orp_d);
	cudaFree(orpk_d);
	//----- ----- ----- ----- 
	cudaFree(dgdcu_d);
	cudaFree(dgdcuk_d);
	//
	cudaFree(dgdmn_d);
	cudaFree(dgdmnk_d);
	//
	cudaFree(dgdni_d);
	cudaFree(dgdnik_d);
	//
	cudaFree(dgdor_d);
	cudaFree(dgdork_d);
	//----- ----- ----- ----- ----- -----
	free(kx);
	free(ky);
	free(k2);
	free(k4);
	//----- ----- ----- ----- 
	free(cu);
	free(cuk);
	free(cuc);
	//
	free(mn);
	free(mnk);
	free(mnc);
	//
	free(ni);
	free(nik);
	free(nic);
	//
	free(orp);
	free(orpk);
	free(orpc);
	//----- ----- ----- ----- 
	free(dgdcu);
	free(dgdcuk);
	//
	free(dgdmn);
	free(dgdmnk);
	//
	free(dgdni);
	free(dgdnik);
	//
	free(dgdor);
	free(dgdork);
	//----- ----- ----- ----- ----- -----
}
