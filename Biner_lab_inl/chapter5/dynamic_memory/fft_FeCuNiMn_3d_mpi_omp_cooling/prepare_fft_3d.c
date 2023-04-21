#include <math.h> //M_PI
#include <stdlib.h> //rand()

void prepare_fft_3d(int Nx, int Ny, int Nz, 
	double dx, double dy, double dz,
	double *kx, double *ky, double *kz, 
	double *k2, double *k4){
	
	int Nx2=Nx/2;
	int Ny2=Ny/2;
	int Nz2=Nz/2;
	
	int Nxm1=Nx-1;
	int Nym1=Ny-1;
	int Nzm1=Nz-1;
	
	double delkx=(2.0*M_PI)/(Nx*dx);
	double delky=(2.0*M_PI)/(Ny*dy);
	double delkz=(2.0*M_PI)/(Nz*dz);
	
	double fk1;
	double fk2;
	double fk3;
	
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
	for(int k=0;k<Nz2;k++){
		fk3=delkz*k;
		kz[k]=fk3;
		kz[Nzm1-k]=-fk3-delkz;
	}
	
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				//k2[i][j][k] = kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k];
				//k4[i][j][k] = k2[i][j][k]*k2[i][j][k];
				ii=i*Ny*Nz+j*Nz+k;
				k2[ii] = kx[i]*kx[i] + ky[j]*ky[j] + kz[k]*kz[k];
				k4[ii] = k2[ii]*k2[ii];
			}
		}
	}
	
	return;
}