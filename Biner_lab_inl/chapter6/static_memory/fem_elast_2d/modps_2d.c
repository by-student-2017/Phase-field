/* This function calculates the elasticity materix for
   plane-stress or plane-strain. */

/* Variable and array list
  mtype: Material type, for the current element (int)
  ntype: solution type, ntype=1 for plane-stress,
  ntype=2 for plane-strain (int)
  nstre: Number of stress components (int)
  props[ngaus]: Position of sampling points (double)
  dmatx[nstre][nstre]: Elasticity matrix (double)
*/

void modps_2d(int mtype, int ntype, int nstre,
	double *props, double *dmatx){
	
	//Material Parameters
	/* Determine the values of the Young's modulus and Poisson's ratio */
	young=props[mtype][0];
	poiss=props[mtype][1];
	
	for(int istre=0;istre<3;istre++){
		for(int jstre=0;jstre<3;jstre++){
			dmatx[istre][jstre]=0.0;
		}
	}
	
	double dconst;
	
	/* Evaluate the elasticity matrix for plane-stress (Eq.6.39) */
	if(ntype == 1){
		//Plane Stress
		dconst = young/(1.0 - poiss * poiss);
		dmatx[0][0] = dconst;
		dmatx[1][1] = dconst;
		dmatx[0][1] = dconst * poiss;
		dmatx[1][0] = dconst * poiss;
		dmatx[2][2] = (1.0 - 2.0*poiss)*dconst/2.0;
	}
	
	/* Evaluate the elasticity materix for plane-strain (Eq.6.40). */
	if(ntype == 2){
		//Plane strain
		dconst = young * (1.0 - poiss)/( (1+poiss)*(1.0 - 2.0*poiss) );
		dmatx[0][0] = dconst;
		dmatx[1][1] = dconst;
		dmatx[0][1] = dconst * poiss / (1.0 - poiss);
		dmatx[1][0] = dconst * poiss / (1.0 - poiss);
		dmatx[2][2] = (1.0 - 2.0*poiss)*dconst/(2.0*(1.0 - poiss));
	}
	
	return;
}