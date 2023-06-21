#include <math.h> //M_PI
#include <stdlib.h> //rand()

#include <mpi.h> //mpi version
#include <fftw3-mpi.h> //mpi version

/* Variable and array list (Note: imaginary part = 0.0)
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  Nz: Number of grid points in the z-direction
  cu0: Initial Cu concentration
  mn0: Initial Mn concentration
  ni0: Initial Ni concentration
  cu(Nx,Ny,Nz): Modulated Cu values at the grid points
  mn(Nx,Ny,Nz): Modulated Mn values at the grid points
  ni(Nx,Ny,Nz): Modulated Ni values at the grid points
  orp(Nx,Ny,Nz): Modulated order parameter values at the grid points */

void init_FeCuMnNi_micro_3d_mpi(int Nx, int Ny, int Nz, 
	double cu0, double mn0, double ni0,
	double *cu, double *mn, double *ni, double *orp,
	ptrdiff_t local_n0, ptrdiff_t local_0_start){
	
	//The value of the noise term
	double noise=0.001;
	
	int iimpi;
	
	srand((unsigned int)local_0_start);
	
	/* Assign the modulated alloy concentrations and 
	   the order parameter values to the grid points */
	for(int i=0;i<local_n0;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				iimpi=((local_0_start+i)*Ny+j)*Nz+k;
				//
				 cu[iimpi] =   cu0 + noise*(0.5-(double)rand()/RAND_MAX);
				 mn[iimpi] =   mn0 + noise*(0.5-(double)rand()/RAND_MAX);
				 ni[iimpi] =   ni0 + noise*(0.5-(double)rand()/RAND_MAX);
				//
				orp[iimpi] = 0.001 + noise*(0.5-(double)rand()/RAND_MAX);
			}
		}
	}
	
	return;
}