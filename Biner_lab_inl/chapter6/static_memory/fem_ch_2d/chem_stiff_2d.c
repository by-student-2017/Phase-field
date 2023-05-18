/* This function forms the global stiffness and global rhs, load, vector for solution of
   Cahn-Hilliard equation for the modified Newton-Raphson solultion algorithm. */

/* Variable and array list
  npoin: Total number of nodes in the solution
  nelem: Total number of elements in the solution
  nnode: Number of nodes per element
  nstre: Number of stress compoents
  ndime: Number of Cartesian coordinate dimension
  ndofn: Number of DOF per node
  ngaus: The order of numerical integration
  
*/

void sfr2_2d();
void jacob2_2d();
void free_energy_fem_2d();

void chem_stiff_2d(int ncountm, int ncounts,
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