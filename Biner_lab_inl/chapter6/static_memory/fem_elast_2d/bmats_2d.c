/* This function froms the strain matrix by using
   the Cartesian derivatives of the shape functions. */

/* Variable and array list
  nnode: Number of node per element.
  shape[nnode]: Shape function values.
  cartd[ndime][nnode]: cartesian derivatives of
    shape functions at the current integration point.
  bmatx[nstre][nevab]: Strain matrix at the current integration point
    (nevab = nnode * ndofn ).
*/

void modps_2d(double *cartd, double *shape,
	int nnode, double *bmatx){
	
	int ngash=0;
	
	for(int inode=0;inode<nnode;inode++){
		
		// Counter for the elements of the strain matrix
		mgash = ngash + 1;
		ngash = mgash + 1;
		
		// Evalutate the strain matrix (Eq.6.37)
		bmatx[0][mgash] = cartd[0][inode];
		bmatx[0][ngash] = 0.0;
		//
		bmatx[1][mgash] = 0.0;
		bmatx[1][ngash] = cartd[1][inode];
		//
		bmatx[2][mgash] = cartd[1][inode];
		bmatx[2][ngash] = cartd[0][inode];
	}
	
	return;
}