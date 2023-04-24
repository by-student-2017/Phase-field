#include <math.h> //M_PI
#include <stdlib.h> //rand()

#include <mpi.h> //mpi version
#include <fftw3-mpi.h> //mpi version

void prepare_fft_2d_mpi(int Nx, int Ny, 
	double dx, double dy,
	double *kx, double *ky, 
	double *k2, double *k4,
	ptrdiff_t local_n0, ptrdiff_t local_0_start){
	
	int Nx2=Nx/2;
	int Ny2=Ny/2;
	
	//int Nxm1=Nx-1;
	int Nym1=Ny-1;
	
	double delkx=(2.0*M_PI)/(Nx*dx);
	double delky=(2.0*M_PI)/(Ny*dy);
	
	double fk1;
	double fk2;
	
	int ii;
	
	for(int i=0;i<local_n0;i++){
		if((local_0_start+i)<Nx2){
			fk1=delkx*(local_0_start+i);
		}else{
			fk1=delkx*((local_0_start+i)-Nx);
		}
		kx[local_0_start+i]=fk1;
	}
	for(int j=0;j<Ny2;j++){
		fk2=delky*j;
		ky[j]=fk2;
		ky[Nym1-j]=-fk2-delky;
	}
	
	//for(int i=0;i<Nx;i++){
	for(int i=0;i<local_n0;i++){
		for(int j=0;j<Ny;j++){
			//k2[i][j] = kx[i]*kx[i] + ky[j]*ky[j];
			//k4[i][j] = k2[i][j]*k2[i][j];
			ii=(local_0_start+i)*Ny+j;
			k2[ii] = kx[local_0_start+i]*kx[local_0_start+i] + ky[j]*ky[j];
			k4[ii] = k2[ii]*k2[ii];
		}
	}
	
	return;
}