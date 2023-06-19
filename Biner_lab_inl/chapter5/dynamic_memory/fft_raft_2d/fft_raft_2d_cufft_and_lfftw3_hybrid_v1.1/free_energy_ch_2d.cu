/* This function calculates the derivative of free energy, Eq. 4.18 at
   the current gid point in the main program */

/* Variable and array list
  i: Current grid point in the x-direction
  j: Current grid point in the y-direction
  dfdcon: Derivative of free energy (dummy)
  con[Nx][Ny]: Concentration filed for un-optimized mode (iflag=1)
*/
float free_energy_ch_2d(float con_ij){
	
	// Constant in free energy function, Eq.4.18
	float A=1.0;
	
	float dfdcon_ij;
	
	/* Derivative of free energy, Eq.4.18 at the current grid point in the main program.
	   
	   mu = df/dc = A*( 2*c[i][j]*(1-c[i][j])^2 + 2*(c[i][j])^2*(1-c[i][j]) ) (Eq.4.18)
	   where mu is The functional derivative of chemical/bulk energy,
	          c is concentration.
	          chemical/bulk energy, f(c) = A*c^2*(1-c)^2 (Eq.4.13)
	*/
	dfdcon_ij=A*(2.0*con_ij*(1.0-con_ij)*(1.0-con_ij)
				-2.0*con_ij*con_ij*(1.0-con_ij)
				);
	
	return dfdcon_ij;
}