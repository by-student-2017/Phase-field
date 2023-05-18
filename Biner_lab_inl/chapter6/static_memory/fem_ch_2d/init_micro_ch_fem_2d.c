/* This function modulates the initial microstructure and
   initializes nodal variable list, according to Eq.6.66. */

#include <stdlib.h> //rand()

/* Variable and array list
  npoin: Total number of nodes in the solution (int)
  ndofn: Number of DOF per node (int)
  conc0: Alloy concentration (double)
  nodco[ntotv]: Nodal variable list, ntotv = npoin * ndofn.
*/

void init_micro_ch_fem_2d(int npoin, int ndofn, double conc0, 
	double *nodco){
	
	//Total number of variables in the solution
	int npoin2 = 2*npoin;
	
	double noise = 0.02;
	
	double nodco[npoin2];
	for(int ipoin=0;ipoin<npoin;ipoin++){
		// Modulate the composition field
		nodco[ipoin] = conc0 + noise*(0.5-(double)rand()/RAND_MAX);
		nodco[ipoin + npoin] = 0.0;
	}
	
	return;
}