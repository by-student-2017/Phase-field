/* This function calculates the derivative of free energy, Eq. 4.18 at
   the current gid point in the main program */

/* Variable and array list
  i: Current grid point in the x-direction
  j: Current grid point in the y-direction
  ngrain: Number of grains in the simulation
  igrain: Current grain number under consideration in the main program
  dfdeta: The derivative of the free energy for the current grid point (output)
  etas[Nx][Ny][ngrain]: common array containing order parameters of grains, isolve=1
  eta[Nx][Ny]: Order parameters of current grain under consideration in the main program
*/

double free_energy_fd_ca_2d(int i, int j, int Nx, int Ny,
	int ngrain, double *etas, double *eta, int igrain){
	
	// Constant in free energy function, Eq.4.33
	double A=1.0;
	double B=1.0;
	
	double sum=0.0;
	
	int ij=i*Nx+j;
	
	// Summation in the free energy function (Eq.4.33)
	for(int jgrain=0;jgrain<ngrain;jgrain++){
		if(jgrain != igrain){
			sum = sum + etas[ij*ngrain+jgrain]*etas[ij*ngrain+jgrain];
		}
	}
	
	/* The derivative of the free energy for the grid point that is under
	   consideration in the main program */
	double dfdeta_ij=A*(2.0*B*eta[ij]*sum + eta[ij]*eta[ij]*eta[ij] - eta[ij]);
	
	return dfdeta_ij;
}