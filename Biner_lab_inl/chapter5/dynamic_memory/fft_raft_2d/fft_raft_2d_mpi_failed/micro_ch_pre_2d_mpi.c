/* This function initializes the microstructure for
   given average composition modulated with
   a noise term to account the thermal fluctuations in
   Cahn-Hilliard equation. */

#include <stdlib.h> //rand()
//#include <fftw3.h>
#include <mpi.h> //mpi version
#include <fftw3-mpi.h> //mpi version

/* Variable and array list
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  c0: Average alloy composition
  con[Nx][Ny]: Concentration field for un-optimized mode (iflag=1)
*/

void micro_ch_pre_2d_mpi(int Nx, int Ny, double c0,
	fftw_complex *con, ptrdiff_t local_n0, ptrdiff_t local_0_start){
	
	// Total number of grid points in the simulation cell
	//int NxNy=Nx*Ny;
	
	//Set the magnitude of the noise term for fluctuations
	double noise=0.02;
	
	int ii;
	
	srand((int)local_0_start);
	
	//Introduce random flucturation to concentration, con[Nx][Ny]
	for(int i=0;i<local_n0;i++){
		for(int j=0;j<Ny;j++){
			ii=i*Ny+j;
			con[ii][0] = c0 + noise*(0.5-(double)rand()/RAND_MAX);
			con[ii][1] = 0.0;
		}
	}
	
	return;
}