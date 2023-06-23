#include <stdio.h>
#include <stdlib.h> //rand()
#include <math.h> //M_PI

void write_vtk_grid_values_2D(int nx, int ny, 
	double dx, double dy, double dx0,
	int istep, double *data1, double *data2){
	
	//open output file
	char fname[256];
	sprintf(fname,"time_%d.vtk",istep);
	FILE *out=fopen(fname,"w");
	
	int nz=1;
	int npoin=nx*ny*nz;
	
	//start writing ASCII VTK file
	
	//header of VTK file
	fprintf(out,"# vtk DataFile Version 3.0 \n");
	fprintf(out,"time_10.vtk \n");
	fprintf(out,"ASCII \n");
	fprintf(out,"DATASET STRUCTURED_POINTS \n");
	
	//coords of grid points
	fprintf(out,"DIMENSIONS %5d %5d %5d \n",nx,ny,nz);
	fprintf(out,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(out,"ASPECT_RATIO %f %f %f \n",(float)(dx*nx)/(float)(dx0*nx),(float)(dy*ny)/(float)(dx0*nx),1.0);
	
	//write grid point values
	fprintf(out,"POINT_DATA %15d \n",npoin);
	fprintf(out,"SCALARS DEN float 1 \n");
	fprintf(out,"LOOKUP_TABLE default \n");
	
	for(int j=0;j<ny;j++){
		for(int i=0;i<nx;i++){
			//ii=i*ny+j;
			//fprintf(out,"%14.6e \n",data1[i][j]);
			fprintf(out,"%14.6e \n",data1[i*ny+j]);
		}
	}
	
	fprintf(out,"SCALARS ENEG float 1 \n");
	fprintf(out,"LOOKUP_TABLE default \n");
	
	for(int j=0;j<ny;j++){
		for(int i=0;i<nx;i++){
			//ii=i*ny+j;
			//fprintf(out,"%14.6e \n",data2[i][j]);
			fprintf(out,"%14.6e \n",data2[i*ny+j]);
		}
	}
	
	fclose(out);
	return;
}