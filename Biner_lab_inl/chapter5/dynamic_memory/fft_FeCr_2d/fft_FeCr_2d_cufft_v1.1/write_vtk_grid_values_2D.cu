#include <stdio.h>
#include <stdlib.h> //rand()
#include <math.h> //M_PI

void write_vtk_grid_values_2D(int nx, int ny, 
	float dx, float dy,
	int istep, float *data1){
	
	//open output file
	char fname[256];
	sprintf(fname,"time_%d.vtk",istep);
	FILE *out=fopen(fname,"w");
	
	int nz=1;
	int npoin=nx*ny*nz;
	float x,y,z;
	
	//start writing ASCII VTK file
	
	//header of VTK file
	fprintf(out,"# vtk DataFile Version 2.0 \n");
	fprintf(out,"time_10.vtk \n");
	fprintf(out,"ASCII \n");
	fprintf(out,"DATASET STRUCTURED_GRID \n");
	
	//coords of grid points
	fprintf(out,"DIMENSIONS %5d %5d %5d \n",nx,ny,nz);
	fprintf(out,"POINTS %7d float \n",npoin);
	
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			x=dx*i;
			y=dy*j;
			z=0.0;
			
			fprintf(out,"%14.6e %14.6e %14.6e \n",x,y,z);
		}
	}
	
	//write grid point values
	fprintf(out,"POINT_DATA %10d \n",npoin);
	fprintf(out,"SCALARS Concentration float 1 \n");
	fprintf(out,"LOOKUP_TABLE default \n");
	
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			//ii=i*ny+j;
			//fprintf(out,"%14.6e \n",data1[i][j]);
			fprintf(out,"%14.6e \n",data1[i*ny+j]);
		}
	}
	
	fclose(out);
	return;
}