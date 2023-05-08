/* This function calculates the derivative for
   concentration field and for order parameters at
   the current grid point under consideration in
   the main program */

/* Variable and array list
  i: Current grid point in the x-direction
  j: Current grid point in the y-direction
  iflag: iflag=1 derivative is provided for the concentration filed,
    iflag=2 the derivative is provided for the order parameter
  npart: Number of particles in the simulation
  dfdcon: The value of derivative of free energy function respect to
    concentration at the current grid point
  dfdeta: The value of derivative of free energy function respect to
    current order parameter at the current grid point
  con[Nx][Ny]: Concentration field
  eta[Nx][Ny]: Order parameter that is under consideration in the main program
  etas[Nx][Ny][ngrain]: Common array containing order parameters
*/

double free_energy_sint_2d_eta(int i, int j, int Nx, int Ny,
	double *con, double *eta, double *etas,
	int npart){
	
	// Constant in free energy function, Eq.4.39
	//double A=16.0;
	double B=1.0;
	
	//derivative of free energy for etas
	double sum2=0.0;
	
	int ij=i*Nx+j;
		
	double etas_ijp;
	
	//Summations in Eq.4.39 for ther current point
	for(int ipart=0;ipart<npart;ipart++){
		etas_ijp = etas[ij*npart+ipart];
		sum2 = sum2 + ( etas_ijp * etas_ijp );
	}
	
	/* Functional derivative of the free energy for
	   the current order parameter at current grid point (i,j) */
	double dfdeta_ij=B*(-12.0*eta[ij]*eta[ij]*(2.0-con[ij])
						+12.0*eta[ij]*(1.0-con[ij])
						+12.0*eta[ij]*sum2);
	
	return dfdeta_ij;
}