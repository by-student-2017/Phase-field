#include <math.h> //M_PI
#include <stdlib.h> //rand()

void prepare_fft_2d(int Nx, int Ny, 
	float dx, float dy,
	float *kx, float *ky, 
	float *k2, float *k4){
	
	int Nx2=Nx/2;
	int Ny2=Ny/2;
	
	int Nxm1=Nx-1;
	int Nym1=Ny-1;
	
	float delkx=(2.0*M_PI)/(Nx*dx);
	float delky=(2.0*M_PI)/(Ny*dy);
	
	float fk1;
	float fk2;
	
	int ii;
	
	for(int i=0;i<Nx2;i++){
		fk1=delkx*i;
		kx[i]=fk1;
		kx[Nxm1-i]=-fk1-delkx;
	}
	for(int j=0;j<Ny2;j++){
		fk2=delky*j;
		ky[j]=fk2;
		ky[Nym1-j]=-fk2-delky;
	}
	
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			//k2[i][j] = kx[i]*kx[i] + ky[j]*ky[j];
			//k4[i][j] = k2[i][j]*k2[i][j];
			ii=i*Ny+j;
			k2[ii] = kx[i]*kx[i] + ky[j]*ky[j];
			k4[ii] = k2[ii]*k2[ii];
		}
	}
	
	return;
}