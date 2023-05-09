/* This function introduces initial solid nuclei in
   the center of the simulation cell. The size of
   the nuclei is in terms of grid numbers. */

#include <math.h> //mod() and -lm
//#include <stdlib.h> //rand()
//#include <fftw3.h>

/* Variable and array list
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  seed: Initial nuclei size (in terms of grid numbers)
  phi[Nx][Ny]: Phase-field parametr
  tempr[Nx][Ny]: temperature
*/

void nucleus_2d(int Nx, int Ny,
	double seed,
	double *phi, double *tempr){
	
	int ij; //ij=i*Ny+j;
	
	// Initialize phi and tempr arrays with longhand format
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			ij=i*Ny+j;
			phi[ij]=0.0;
			tempr[ij]=0.0;
		}
	}
	
	/* Compare the vector length starting from
	   the center of the simulation cell. If it is smaller or
	   equal to the value of seed, assign value of one to phi. */
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			ij=i*Ny+j;
			if( (i-Nx/2)*(i-Nx/2)+(j-Ny/2)*(j-Ny/2) < seed ){
				phi[ij]=1.0;
			}
		}
	}
	
	return;
}