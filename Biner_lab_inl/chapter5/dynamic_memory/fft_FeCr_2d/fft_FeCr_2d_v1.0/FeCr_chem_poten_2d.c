/* This function calculates the derivative of
   the chemical energy with respect to Cr concentration at
   the current grid point under consideration in the main program */

/* Variable and array list
  i: Current grid point in the x-direction
  j: Current grid point in the y-direction
  tempr: Temperature
  dfdcon: Derivative of free energy (dummy)
  cr[Nx][Ny]: Cr concentration in the grid points
*/

#include <math.h> //mod() and -lm

double FeCr_chem_poten_2d(double cr_ij, double tempr){
	
	// Value of gas constant x temeprature
	double RT=8.314462*tempr;
	
	double dfdcr_ij;
	
	/* Derivative of chemical energy with respect to Cr at current grid point.
	   Note that it is normallized with the value of RT */
	dfdcr_ij=(-cr_ij*(20500.0-9.68*tempr)+(1.0-cr_ij)*(20500.0-9.68*tempr)
			  +(log(cr_ij)-log(1.0-cr_ij))*RT
			 )/RT;
	
	return dfdcr_ij;
}