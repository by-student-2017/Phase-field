/* This function precalculates the Cartesian derivatives of
   shape functions at all integration of points of all elements in the simulation. 
   In addition, function also calculates the area/volume contributions of the integration points.*/

/* Variable and array list
  npoin: Total number of nodes in the solution (int)
  nelem: Total number of elements in the solution (int)
  nnode: Number of nodes per elements (int)
  nstre: Number of stress components (int)
  ndime: Number of Cartesian coordinate dimension (int)
  ndofn: Number of DOF per node (int)
  ngaus: The order of numerical integration (int)
  ntype: Solution type (ntype=1 for plane-stress and,
    ntype=2 for plane-strain.) (int)
  posgp[ngaus]: Position of sampling points (double)
  weigp[ngaus]: Weighting factors at the sampling points (double)
  lnods[nelem][nnode]: Element nodal connectivity list (int)
  coord[npoin][ndime]: Cartesian coordinates of nodes (double)
  dvolum[nelem][mgaus]: Contributions of integration points to element area/volume (double)
  dgdx[nelem][mgaus][ndime][nnode]: Cartesian derivatives of shape functions at integration points (double)
*/

void sfr2_2d();
void jacob3_2d();

void cart_deriv_2d(int ncountm, int ncounts,
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