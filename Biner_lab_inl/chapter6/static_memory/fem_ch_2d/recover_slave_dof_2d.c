/* This function assigns the solution values to the slave nodes. */

/* Variable and array list
  npoin: Total number of nodes in the solution.
  ndofn: Number of DOF per node
  ncountm: Total number of master nodes
  ncounts: Total number of slave node
  master[ncountm]: Listing of master nodes
  slave[ncounts]: Listing of slave nodes
  asdis[ntotv]: Nodal values, ntotv = npoin * ndofn.
*/

void recover_slave_dof_2d(double *asdis,
	int ncountm, int ncounts,
	int *master, int *slave,
	int npoin, int ndofn){
	
	/* Loop over the number of master nodes and
	   assign the solution obtained for the master nodes to
	   corresponding slave nodes. */
	int im, is;
	for(int ipbc=0;ipbc<ncountm;ipbc++){
		im = master[ipbc];
		is = slave[ipbc];
		asdis[is] = asdis[im];
	}
	
	return;
}