/* This function calculates the derivative of free energy, Eq. 4.18 at
   the current gid point in the main program */

/* Variable and array list
  i: Current grid point in the x-direction
  j: Current grid point in the y-direction
  k: Current grid point in the z-direction
  dfdcon: Derivative of free energy (dummy)
  con[Nx][Ny]: Concentration filed for un-optimized mode (iflag=1)
*/

float free_energy_ch_3d(float con_ijk){
	
	// Constant in free energy function, Eq.4.18
	float A=1.0;
	
	float dfdcon_ijk;
	
	/* Derivative of free energy, Eq.4.18 at the current grid point in the main program.
	   
	   mu = df/dc = A*( 2*c[i][j][k]*(1-c[i][j][k])^2 + 2*(c[i][j][k])^2*(1-c[i][j][k]) ) (Eq.4.18)
	   where mu is The functional derivative of chemical/bulk energy,
	          c is concentration.
	          chemical/bulk energy, f(c) = A*c^2*(1-c)^2 (Eq.4.13)
	*/
	dfdcon_ijk=A*(2.0*con_ijk*(1.0-con_ijk)*(1.0-con_ijk)
				 -2.0*con_ijk*con_ijk*(1.0-con_ijk)
				 );
	
	return dfdcon_ijk;
}