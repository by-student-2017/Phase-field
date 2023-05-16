/* This function writes the FEM mesh information and
   the nodal values in vtk file format to be viewed by
   using Paraview */

#include <stdio.h> //printf()

/* Variable and array list
npoin: Total number of nodes in the solution (int)
  nelem: Total number of elements in the solution (int)
  nnode: Number of nodes per element
  cont1[npin]: Number of nodes per element (int)
  lnods[nelem][nnode]: Nodal connectivity list of elements (int)
  coord[npoin][ndime]: Cartesian coordinates of the nodes (double)
*/

void write_vtk_fem_2d(int npion, int nelem, int nnode,
	int *lnods, double *coord, int istep, int cont1){
	
	//open output file
	char fname[256];
	sprintf(fname,"time_%d.vtk",istep);
	FILE *out=fopen(fname,"w");
	
	//start writing ASCII VTK file
	
	//header of VTK file
	// Header of the vtk file, these lines should not be modified
	fprintf(out,"# vtk DataFile Version 2.0 \n");
	fprintf(out,"time_10.vtk \n");
	fprintf(out,"ASCII \n");
	fprintf(out,"DATASET UNSTRUCTURED_GRID \n");
	
	//write nodal coordinates
	// Write the Cartesian coordinates of the nodes
	fprintf(out,"DIMENSIONS %5d float \n",npoin);
	
	double dummy = 0.0;
	
	int ndime = 2;
	for(int ipoin=0;ipoin<npoin;ipoin++){
		fprintf(out,"%14.6e %14.6e %14.6e \n",coord[ipoin*ndime + 0],coord[ipoin*ndime + 1],dummy);
	}
	
	//write element connectivity
	// Write the element connectivity list in vtk format
	int iconst1 = nelem * nnode;
	fprintf(out,"CELLS %5d %5d \n",nelem,iconst1);
	
	for(int ielem=0;ielem<nelem;ielem++){
		fprintf(out,"%5d \n",nnode);
		for(int inode=0;inode<nnode;inode++){
			fprintf(out,"%5d \n",(lnods[ielem][inode]-1));
		}
		fprintf(out,"%5d \n");
	}
	
	//write cell types
	// For each element write its cell type
	if(nnode == 8){
		ntype = 23;
	}
	if(nnode == 4){
		ntype = 9;
	}
	if(nnode == 3){
		ntype = 5;
	}
		
	fprintf(out,"CELL_TYPES %5d \n",nelem);
	for(int i=0;i<nelem;i++){
		fprintf(out,"%2d \n",ntype);
	}
	
	//write Nodal scaler & vector values
	/* Write the values of the nodal variables.
	   In this case, they are the stress components. */
	fprintf(out,"POINT_DATA %5d float \n",npoin);
	
	//write stress values as scalar
	fprintf(out,"SCALARS sigma_xx float 1 \n");
	fprintf(out,"LOOKUP_TABLE default \n");
	for(int ipoin=0;ipoin<npoin;ipoin++){
		fprintf(out,"%14.6e \n",cont1[ipoin][0]);
	}
	
	fprintf(out,"SCALARS sigma_yy float 1 \n");
	fprintf(out,"LOOKUP_TABLE default \n");
	for(int ipoin=0;ipoin<npoin;ipoin++){
		fprintf(out,"%14.6e \n",cont1[ipoin][1]);
	}
		
	fprintf(out,"SCALARS sigma_xy float 1 \n");
	fprintf(out,"LOOKUP_TABLE default \n");
	for(int ipoin=0;ipoin<npoin;ipoin++){
		fprintf(out,"%14.6e \n",cont1[ipoin][2]);
	}
	
	fclose(out);
	return;
}