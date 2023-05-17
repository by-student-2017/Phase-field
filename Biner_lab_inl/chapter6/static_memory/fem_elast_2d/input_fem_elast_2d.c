/* This function reads the input file for the FEM analysis of
   two-dimensional elasticity using three-,four-,and eight-node
   isoparametric elements. */

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

void input_fem_elast_2d(int npoin, int nelem, int nvfix, int ntype,
	int nnode, int ndofn, int ndime, int ngaus, int nmats, int nstre, int nprop,
	int *matno, int *nofix, int *lnods, double *coord, int *iffix, double *fixed, double *props,
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
			fscanf(in,"%5d",&lnods[jelem][inods]);
		}
		fscanf(in,"%5d",&natno[jelem]);
	}
	
	//Nodal coordinates
	// Read nodal coordinates (coord)
	int jpoin;
	for(int ipoin=0;ipoin<npoin;ipoin++){
		fscanf(in,"%5d",&jpoin);
		for(int idime=0;idime<ndime;idime++){
			fscanf(in,"%lf",&coord[ipoin][idime]);
		}
	}
	
	//Constraint nodes and their values
	/* Read the data on constrained nodes
	   (nofix, iffix, and fixed). */
	for(int ivfix=0;ivfix<nvfix;ivfix++){
		fscanf(in,"%5d",&nofix[ivfix]);
		for(int idofn=0;idofn<ndofn;idofn++){
			fscanf(in,"%5d",&iffix[ivfix][idofn]);
			fscanf(in,"%lf",&fixed[ivfix][idofn]);
		}
	}
	
	//Material properties
	// Read each material properties
	int jmats;
	for(int imats=0;imats<nmats;imats++){
		fscanf(in,"%5d",&jmats);
		for(int iprop=0;iprop<nprop;iprop++){
			fscanf(in,"%lf",&props[jmats][iprop]);
		}
	}
	
	//print out
	// Print out solution input parameters to the output file
	fprintf(out, "******************\n");
	fprintf(out, "* FEM input data *\n");
	fprintf(out, "******************\n");
	fprintf(out, "\n");
	//
	fprintf(out, "Number of Elements            : %5d \n",nelem);
	fprintf(out, "Number of Node                : %5d \n",npoin);
	fprintf(out, "Number of Fixed nodes         : %5d \n",nvfix);
	fprintf(out, "Number of Nodes per element   : %5d \n",nnode);
	fprintf(out, "Number of Integreation points : %5d \n",ngaus);
	fprintf(out, "Number of Materials           : %5d \n",nmats);
	fprintf(out, "Number of Propertiers         : %5d \n",nprop);
	fprintf(out, "\n");
	
	if(ntype == 1){
		fprintf(out, "Plane-stress elasticity solution \n");
	}
	if(ntype == 2){
		fprintf(out, "Plane-strain elasticity solution \n");
	}
	
	
	
	return;
}