/* This function evaluates the functional derivative of 
   free energy function, Eq.4.59, for the current cell under
   consideration at current grid locations given by the main program. */

/* Variable and array list
  i: Current grid id in the x-direction
  j: Current grid id int he y-direction
  icell: Cell id of the current cell being considered in the main program
  ncell: Total number of cell in the simulation
  gamma: Parameter in the free energy functional, Eq.4.59
  kappa: Parameter in the free energy functional, Eq.4.59
  dfdphi: Functional derivative of the free energy (out)
  phi[Nx][Ny]: The value of the order parameter at grid point (i,j)
  phis[Nx][Ny][ncell]: Common array containing the order parameters of the cells
*/

#include <math.h> //M_PI
#include <stdlib.h> //rand()

double free_energy_3d(int i, int j, int k, 
	int Nx, int Ny, int Nz,
	int icell, int ncell, int used_ncell,
	double gamma, double kappa, 
	double *phi, double *phis){
	
	double sum_phi=0.0;
	
	int ijk=(i*Ny+j)*Nz+k;
	int ijkcj;
	
	/* Go over total number of cells in the simulation
	   if jcell counter is not icell add their value to sum_phi */
	for(int jcell=0;jcell<used_ncell;jcell++){
		if( icell != jcell ){
			ijkcj=(ijk*ncell+jcell);
			//sum_phi = sum_phi + phis[i][j][k][jcell]*phis[i][j][k][jcell];
			sum_phi = sum_phi + phis[ijkcj]*phis[ijkcj];
		}
	}
	
	/* Calculate the functional derivative of the free energy for 
	   the current cell, the terms in the first bracket, Eq.4.59, at
	   the current grid point location. */
	//double dfdphi = gamma*phi[i][j][k]*(1.0-phi[i][j][k])*(1.0-2.0*phi[i][j][k])
	//		   +2.0*kappa*phi[i][j][k]*sum_phi;
	double dfdphi = gamma*phi[ijk]*(1.0-phi[ijk])*(1.0-2.0*phi[ijk])
			   +2.0*kappa*phi[ijk]*sum_phi;
	
	return dfdphi;
}