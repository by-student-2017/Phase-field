/* This function modifies the global stiffness matrix and rhs,
   load, vector to impose the periodic boundary conditions to FEM mesh.
   the nodal infromation is determined in function periodic_boundary_2d.c.
   Note this algorithm only works if mesh is uniform and the Cartesian coordinates of 
   the boundary nodes also have symmetry properties. The algorithms to impose
   periodic boundary conditions for arbitrary FEM meshes are given in [5-7]. */

/* Variable and array list
  ncountm: Number of master nodes (int)
  ncounts: Number of slave nodess (int)
  ndofn: Number DOF per node (int)
  iter: Iteration number in Newton-Raphson algorithm (int)
  master[ncountm]: Array containing the node numbers of master nodes (int)
  slave[ncounts]: Array containing the node numbers of slave nodes (int)
  gforce[ntotv]: global rhs, load, vector (double)
  gstif[ntotv][ntotv]: global stiffness matrix (ntotv = npoin * ndofn). (double)
*/

void apply_periodic_bc_2d(int ncountm, int ncounts,
	int *master, int *slave,
	int ndofn, int npoin,
	double *gstif, double *gforce,
	int iter){
	
	int im,is;
	for(int ipbc=0;ipbc<ncountm;ipbc++){
		im = master[ipbc];
		is = slave[ipbc];
		
		if(iter == 0){
			for(int itotv=0;itotv<ntotv;itotv++){
				//Add rows
				gstif[im*ntotv + itotv] += gstif[is*ntotv + itotv];
				
				//Add columns
				gstif[itotv*ntotv + im] += gstif[itotv*ntotv + is];
				
				//zero slave dofs
				gstif[is*ntotv + itotv] = 0.0;
				gstif[itotv*ntotv + is] = 0.0;
				
				gstif[is*ntotv + is] = 1.0;
			}
		}
		
		//add rhs
		gforce[im] += gforce[is];
		gforce[is] = 0.0;
	}
	
	return;
}