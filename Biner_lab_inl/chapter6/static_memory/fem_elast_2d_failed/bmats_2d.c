/* This function froms the strain matrix by using
   the Cartesian derivatives of the shape functions. */

/* Variable and array list
  nnode: Number of node per element. (int)
  shape[nnode]: Shape function values. (double)
  cartd[ndime][nnode]: cartesian derivatives of
    shape functions at the current integration point. (double)
  bmatx[nstre][nevab]: Strain matrix at the current integration point
    (nevab = nnode * ndofn ). (double)
*/

void bmats_2d(double *cartd, double *shape,
	int nnode, double *bmatx,
	int nevab){
	
	int ngash=0;
	
	for(int inode=0;inode<nnode;inode++){
		
		// Counter for the elements of the strain matrix
		int mgash = ngash + 1;
		int ngash = mgash + 1;
		
		// Evalutate the strain matrix (Eq.6.37)
		bmatx[0*nevab+mgash] = cartd[0*nnode+inode];
		bmatx[0*nevab+ngash] = 0.0;
		//
		bmatx[1*nevab+mgash] = 0.0;
		bmatx[1*nevab+ngash] = cartd[1*nnode+inode];
		//
		bmatx[2*nevab+mgash] = cartd[1*nnode+inode];
		bmatx[2*nevab+ngash] = cartd[0*nnode+inode];
	}
	
	return;
}