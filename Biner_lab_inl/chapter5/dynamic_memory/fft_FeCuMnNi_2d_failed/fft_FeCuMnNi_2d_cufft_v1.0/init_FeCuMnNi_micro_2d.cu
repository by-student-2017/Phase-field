#include <math.h> //M_PI
#include <stdlib.h> //rand()

/* Variable and array list (Note: imaginary part = 0.0)
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  cu0: Initial Cu concentration
  mn0: Initial Mn concentration
  ni0: Initial Ni concentration
  cu(Nx,Ny): Modulated Cu values at the grid points
  mn(Nx,Ny): Modulated Mn values at the grid points
  ni(Nx,Ny): Modulated Ni values at the grid points
  orp(Nx,Ny): Modulated order parameter values at the grid points */

void init_FeCuMnNi_micro_2d(int Nx, int Ny, 
	float cu0, float mn0, float ni0,
	float *cu, float *mn, float *ni, float *orp){
	
	//The value of the noise term
	float noise=0.001;
	
	int ii;
	
	/* Assign the modulated alloy concentrations and 
	   the order parameter values to the grid points */
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			ii=i*Ny+j;
			//
			 cu[ii] =   cu0 + noise*(0.5-(float)rand()/RAND_MAX);
			 mn[ii] =   mn0 + noise*(0.5-(float)rand()/RAND_MAX);
			 ni[ii] =   ni0 + noise*(0.5-(float)rand()/RAND_MAX);
			//
			orp[ii] = 0.001 + noise*(0.5-(float)rand()/RAND_MAX);
		}
	}
	
	return;
}