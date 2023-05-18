/* This function calculates the elasticity materix for
   plane-stress or plane-strain. */

/* Variable and array list
  mtype: Material type, for the current element (int)
  ntype: solution type, ntype=1 for plane-stress,
    ntype=2 for plane-strain (int)
  nstre: Number of stress components (int)
  props[nelem][ndime]: Properties of material type
  dmatx[nstre][nstre]: Elasticity matrix (double)
*/

void modps_2d(int mtype, int ntype, int nstre,
	double *props, double *dmatx,
	int nelem, int ndime){
	
	//Material Parameters
	/* Determine the values of the Young's modulus and Poisson's ratio */
	double young=props[mtype*ndime+0];
	double poiss=props[mtype*ndime+1];
	
	for(int istre=0;istre<3;istre++){
		for(int jstre=0;jstre<3;jstre++){
			dmatx[istre*nstre+jstre]=0.0;
		}
	}
	
	double dconst;
	
	/* Evaluate the elasticity matrix for plane-stress (Eq.6.39) */
	if(ntype == 1){
		//Plane Stress
		dconst = young/(1.0 - poiss * poiss);
		dmatx[0*nstre+0] = dconst;
		dmatx[1*nstre+1] = dconst;
		dmatx[0*nstre+1] = dconst * poiss;
		dmatx[1*nstre+0] = dconst * poiss;
		dmatx[2*nstre+2] = (1.0 - 2.0*poiss)*dconst/2.0;
	}
	
	/* Evaluate the elasticity materix for plane-strain (Eq.6.40). */
	if(ntype == 2){
		//Plane strain
		dconst = young * (1.0 - poiss)/( (1+poiss)*(1.0 - 2.0*poiss) );
		dmatx[0*nstre+0] = dconst;
		dmatx[1*nstre+1] = dconst;
		dmatx[0*nstre+1] = dconst * poiss / (1.0 - poiss);
		dmatx[1*nstre+0] = dconst * poiss / (1.0 - poiss);
		dmatx[2*nstre+2] = (1.0 - 2.0*poiss)*dconst/(2.0*(1.0 - poiss));
	}
	
	return;
}