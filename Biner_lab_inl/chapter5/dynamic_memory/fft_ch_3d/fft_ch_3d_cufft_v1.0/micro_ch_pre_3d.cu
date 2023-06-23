/* This function initializes the microstructure for
   given average composition modulated with
   a noise term to account the thermal fluctuations in
   Cahn-Hilliard equation. */

#include <stdlib.h> //rand()
//#include <fftw3.h>

/* Variable and array list
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  Nz: Number of grid points in the z-direction
  c0: Average alloy composition
  con[Nx][Ny][Nz]: Concentration field for un-optimized mode (iflag=1)
*/

void micro_ch_pre_3d(int Nx, int Ny, int Nz, float c0, float *con){
	
	// Total number of grid points in the simulation cell
	//int NxNy=Nx*Ny;
	
	//Set the magnitude of the noise term for fluctuations
	float noise=0.02;
	
	int ijk;
	
	//Introduce random flucturation to concentration, con[Nx][Ny]
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				ijk=(i*Ny+j)*Nz+k;
				con[ijk] = c0 + noise*(0.5-(float)rand()/RAND_MAX);
			}
		}
	}
	
	return;
}