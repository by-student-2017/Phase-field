/* This function calculates the derivatives of free energy,
   based on the composition value at the integration points. */

/* Variable and array list
  c: The value of the concentration at the current integration point
  dfdc: First derivative of free energy
  df2dc: Second derivative of free energy
*/

void free_energy_fem_2d(double c, double dfdc, double df2dc){
	
	// Value of constant A in free energy function (Eq.6.53)
	double constA = 1.0;
	
	/* First derivative of free energy with respect to
	   concentration at current integration point. */
	dfdc = constA * (2.0*c - 6.0*c*c + 4.0*c*c*c);
	
	/* Second derivative of free energy with respect to
	   concentration at current integration point. */
	df2dc = constA * (2.0 - 12.0*c + 12.0*c*c);
	
	return;
}