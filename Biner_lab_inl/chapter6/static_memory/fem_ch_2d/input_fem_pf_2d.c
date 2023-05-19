/* This function is the shortened version function input_fem_elast_2d.c and
  reads the FEM mesh. Therefore, just its listing is given in here. */

#include <stdio.h> //printf()

/* Variable and array list
  npoin: Total number of nodes in the solution (int)
  nelem: Total number of elements (int)
  nvfix: Total number of constrained nodes (int)
  ntype: Solution type (ntyep=1, plane-stress and ntype=2 plane-strain) (int)
  nnode: Number of nodes per element (int)
  ndofn: Number of degree of freedom (DOF) per node (int)
  ndime: Number of spatial dimensions (int)
  ngaus: The order of numerical integration (int)
  nmats: Total number of different materials in the solution (int)
  nstre: Number of stress components (int)
  nprop: Number of material properties (int)
  matno[nelem]: Material types for the elements (int)
  nofix[nvfix]: Node numbers at which one or more DOFs are constrained (int)
  lnods[nelem][nnode]: Element nodal connectivity list (int)
  coord[npoin][ndime]: Cartesian coordinates of each node (double)
  iffix[nvfix][ndofn]: List of constrained DOFs (int)
  fixed[nvfix][ndofn]: Prescribed value of any constrained DOFs (double)
  props[nmats][nprops]: For each different material, the properties of that material (double)
*/

void input_fem_pf_2d(int npoin, int nelem, int nvfix, int ntype,
	int nnode, int ndofn, int ndime, int ngaus, int nstre, int nmats, int nprop,
	int *lnods,int *matno, double *coord,
	FILE *in, FILE *out){
	
	// Read the control data (npoin...nprop)
	fscanf(in,"%5d",&npoin);
	fscanf(in,"%5d",&nelem);
	fscanf(in,"%5d",&nvfix);
	fscanf(in,"%5d",&ntype);
	fscanf(in,"%5d",&nnode);
	fscanf(in,"%5d",&ndofn);
	fscanf(in,"%5d",&ndime);
	fscanf(in,"%5d",&ngaus);
	fscanf(in,"%5d",&nstre);
	fscanf(in,"%5d",&nmats);
	fscanf(in,"%5d",&nprop);
	
	//Element node numbers & material property number
	/* Read element connectivity list (lnode) and 
	   material number list (matno). */
	int jelem;
	for(int ielem=0;ielem<nelem;ielem++){
		fscanf(in,"%5d",&jelem);
		for(int inode=0;inode<nnode;inode++){
			fscanf(in,"%5d",&lnods[jelem*nnode+inode]);
		}
		fscanf(in,"%5d",&matno[jelem]);
	}
	
	//Nodal coordinates
	// Read nodal coordinates (coord)
	int jpoin;
	for(int ipoin=0;ipoin<npoin;ipoin++){
		fscanf(in,"%5d",&jpoin);
		for(int idime=0;idime<ndime;idime++){
			fscanf(in,"%lf",&coord[ipoin*ndime+idime]);
		}
	}
	
	for(int ipoin=0;ipoin<npoin;ipoin++){
		if(coord[ipoin][0] < 0.0){
			coord[ipoin][0] = 0.0;
		}
		if(coord[ipoin][1] < 0.0){
			coord[ipoin][1] = 0.0;
		}
	}
	
	return;
}