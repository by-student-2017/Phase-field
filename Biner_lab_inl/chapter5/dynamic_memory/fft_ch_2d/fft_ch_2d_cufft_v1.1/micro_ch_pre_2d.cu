/* This function initializes the microstructure for
   given average composition modulated with
   a noise term to account the thermal fluctuations in
   Cahn-Hilliard equation. */

#include <stdlib.h> //rand()
//#include <fftw3.h>

/* Variable and array list
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  c0: Average alloy composition
  con[Nx][Ny]: Concentration field for un-optimized mode (iflag=1)
*/

void micro_ch_pre_2d(int Nx, int Ny, float c0, float *con){
	
	// Total number of grid points in the simulation cell
	//int NxNy=Nx*Ny;
	
	//Set the magnitude of the noise term for fluctuations
	float noise=0.02;
	
	int ii;
	
	//Introduce random flucturation to concentration, con[Nx][Ny]
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			ii=i*Ny+j;
			con[ii] = c0 + noise*(0.5-(float)rand()/RAND_MAX);
		}
	}
	
	return;
}