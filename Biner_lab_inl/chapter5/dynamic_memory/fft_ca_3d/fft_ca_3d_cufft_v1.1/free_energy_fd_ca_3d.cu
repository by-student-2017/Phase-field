/* This function calculates the derivative of free energy, Eq. 4.18 at
   the current gid point in the main program */

/* Variable and array list
  i: Current grid point in the x-direction
  j: Current grid point in the y-direction
  k: Current grid point in the z-direction
  ngrain: Number of grains in the simulation
  igrain: Current grain number under consideration in the main program
  dfdeta: The derivative of the free energy for the current grid point (output)
  etas[Nx][Ny][Nz][ngrain]: common array containing order parameters of grains, isolve=1
  eta[Nx][Ny][Nz]: Order parameters of current grain under consideration in the main program
*/

float free_energy_fd_ca_3d(int i, int j, int k,
	int Nx, int Ny, int Nz,
	int ngrain, float *etas, float *eta, int igrain){
	
	// Constant in free energy function, Eq.4.33
	float A=1.0;
	float B=1.0;
	
	float sum=0.0;
	
	int ijk=(i*Ny+j)*Nz+k;
	
	// Summation in the free energy function (Eq.4.33)
	for(int jgrain=0;jgrain<ngrain;jgrain++){
		if(jgrain != igrain){
			sum = sum + etas[ijk*ngrain+jgrain]*etas[ijk*ngrain+jgrain];
		}
	}
	
	/* The derivative of the free energy for the grid point that is under
	   consideration in the main program */
	float dfdeta_ijk=A*(2.0*B*eta[ijk]*sum + eta[ijk]*eta[ijk]*eta[ijk] - eta[ijk]);
	
	return dfdeta_ijk;
}