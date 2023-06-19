/* This function calculates the derivative of
   the chemical energy with respect to Cr concentration at
   the current grid point under consideration in the main program */

/* Variable and array list
  i: Current grid point in the x-direction
  j: Current grid point in the y-direction
  tempr: Temperature
  dfdcon: Derivative of free energy (dummy)
  cr[Nx][Ny][Nz]: Cr concentration in the grid points
*/

#include <math.h> //mod() and -lm

double FeCr_chem_poten_3d(double cr_ijk, double tempr){
	
	// Value of gas constant x temeprature
	double RT=8.314462*tempr;
	
	double dfdcr_ijk;
	
	/* Derivative of chemical energy with respect to Cr at current grid point.
	   Note that it is normallized with the value of RT */
	dfdcr_ijk=(-cr_ijk*(20500.0-9.68*tempr)+(10.0-cr_ijk)*(20500.0-9.68*tempr)
			 +(log(cr_ijk)-log(1.0-cr_ijk))*RT
			  )/RT;
	
	return dfdcr_ijk;
}