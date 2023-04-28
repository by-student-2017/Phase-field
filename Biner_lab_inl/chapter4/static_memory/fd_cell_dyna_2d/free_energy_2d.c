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

double free_energy_2d(int i, int j, int Nx, int Ny,
	int icell, int ncell,
	double gamma, double kappa, 
	double *phi, double *phis){
	
	double sum_phi=0.0;
	
	/* Go over total number of cells in the simulation
	   if jcell counter is not icell add their value to sum_phi */
	int ijc;
	for(int jcell=0;jcell<ncell;jcell++){
		if( icell != jcell ){
			ijc=((i*Ny+j)*ncell+jcell);
			//sum_phi = sum_phi + phis[i][j][jcell]*phis[i][j][jcell];
			sum_phi = sum_phi + phis[ijc]*phis[ijc];
		}
	}
	
	/* Calculate the functional derivative of the free energy for 
	   the current cell, the terms in the first bracket, Eq.4.59, at
	   the current grid point location. */
	//double dfdphi = gamma*phi[i][j]*(1.0-phi[i][j])*(1.0-2.0*phi[i][j])
	//		   +2.0*kappa*phi[i][j]*sum_phi;
	int ij=i*Ny+j;
	double dfdphi = gamma*phi[ij]*(1.0-phi[ij])*(1.0-2.0*phi[ij])
			   +2.0*kappa*phi[ij]*sum_phi;
	
	return dfdphi;
}