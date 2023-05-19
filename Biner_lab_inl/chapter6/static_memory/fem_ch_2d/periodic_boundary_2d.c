/* This function finds the nodes that are at the boundaries of 
   the simulation cell. For each pair of edge (left-right and bottom-top) assigns 
   master and slave node designation. */

/* Variable and array list
  npoin: Total number of nodes in the solution.
  ndofn: Number of DOF per node
  ncountm: Total number of master nodes
  ncounts: Total number of slave node
  master[ncountm]: Listing of master nodes
  slave[ncounts]: Listing of slave nodes
  asdis[ntotv]: Nodal values, ntotv = npoin * ndofn.
  coord[npoin][ndime]: Cartesian coordinates of each node (double)
*/

void periodic_boundary_2d(int npoin, double *coord,
	int ncountm, int ncounts,
	int *master, int *slave){
	
	// Initialize counters for master and slave nodes
	int ncountm = 0;
	int ncounts = 0;
	
	//determine xmax, xmin, ymax, ymin
	// Determine the maximum and minimum x- and y0 coordinates of the FEM mesh.
	double xmax, xmin;
	double ymax, ymin;
		xmax = xmin = coord[ipoin][0];
		ymax = ymin = coord[ipoin][1];
	for(int ipoin=0;ipoin<npoin;ipoin++){
		if(coord[ipoin][0] > xmax){ xmax=coord[ipoin][0]; }
		if(coord[ipoin][0] < xmin){ xmin=coord[ipoin][0]; }
		//
		if(coord[ipoin][1] > ymax){ ymax=coord[ipoin][1]; }
		if(coord[ipoin][1] < ymin){ ymin=coord[ipoin][1]; }
	}
	
	/* Find the node numbers on the left and right boundaries of the FEM mesh.
	   Assign the nodes residing on the left boundary as master nodes and
	   on the right ones as slave nodes. */
	for(int ipoin=0;ipoin<npoin;ipoin++){
		//
		//left mater nodes
		diff = xmin - coord[ipoin][0];
		if(fabs(diff) <= 1.0e-6){
			ncountm = ncountm + 1;
			master[ncountm] = ipoin;
		}
		
		//right slave nodes
		diff = xmax - coord[ipoin][0];
		if(fabs(diff) <= 1.0e-6){
			ncounts = ncounts + 1;
			slave[ncounts] = ipoin;
		}
		//
	}
	
	/* Find the node numbers on the top and bottom boundaries of the FEM mesh.
	   Assign the nodes residing on the bottom boundary as master nodes and
	   on the top ones as slave nodes. */
	for(int ipoin=0;ipoin<npoin;ipoin++){
		//
		//bottom master nodes
		diff = ymin - coord[ipoin][1];
		if(fabs(diff) <= 1.0e-6){
			ncountm = ncountm + 1;
			master[ncountm] = ipoin;
		}
		
		//top slave nodes
		diff = ymax - coord[ipoin][1];
		if(fabs(diff) <= 1.0e-6){
			ncounts = ncounts + 1;
			slave[ncounts] = ipoin;
		}
		//
	}
	
	if(ncountm != ncounts){
		error("ncountm should be equal ncounts");
	}
	
	return;
}