/* This function multiplies the elasticity matrix and the strain matrix */

/* Variable and array list
  nevab: Total number of element variables (nevab = nnode * ndofn).
  nstre: Number of stress components
  bmatx[nstre][nevab]: Strain matrix at the current integration point
    (nevab = nnode * ndofn ).
  dmatx[nstre][nstre]: Elasticity matrix
  dbmat[nstre][nevab]: Multiplication results of bmatx and dmatx.
*/

void dbe_2d(int nevab, int nstre,
	double *bmatx, double *dmatx, double *dbmat){
	
	/* Matrix multiplication of elasticity matrix (dmatx) and
		strain matrix (bmatx). */
	for(int istre=0;istre<nstre;istre++){
		for(int ievab=0;ievab<nevab;ievab++){
			dbmat[istre][ievab]=0.0;
			for(int jstre=0;jstre<nstre;jstre++){
				dbmat[istre][ievab] = dbmat[istre][ievab]
									+ dmatx[istre][jstre] * bmatx[jstre][ievab];
			}
		}
	}
	
	return;
}