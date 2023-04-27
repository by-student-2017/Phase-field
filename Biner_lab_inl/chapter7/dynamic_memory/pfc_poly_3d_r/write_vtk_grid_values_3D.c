#include <stdio.h>
#include <stdlib.h> //rand()
#include <math.h> //M_PI

void write_vtk_grid_values_3D(int nx, int ny, int nz, 
	double dx, double dy, double dz,
	int istep, double *data1, double *data2){
	
	//open output file
	char fname[256];
	sprintf(fname,"time_%d.vtk",istep);
	FILE *out=fopen(fname,"w");
	
	int npoin=nx*ny*nz;
	double x,y,z;
	
	//start writing ASCII VTK file
	
	//header of VTK file
	fprintf(out,"# vtk DataFile Version 2.0 \n");
	fprintf(out,"time_10.vtk \n");
	fprintf(out,"ASCII \n");
	fprintf(out,"DATASET STRUCTURED_GRID \n");
	
	//coords of grid points
	fprintf(out,"DIMENSIONS %5d %5d %5d \n",nx,ny,nz);
	fprintf(out,"POINTS %15d float \n",npoin);
	
	for(int k=0;k<nz;k++){
		for(int j=0;j<ny;j++){
			for(int i=0;i<nx;i++){
				x=dx*i;
				y=dy*j;
				z=dz*k;
				
				fprintf(out,"%14.6e %14.6e %14.6e \n",x,y,z);
			}
		}
	}
	
	//write grid point values
	fprintf(out,"POINT_DATA %15d \n",npoin);
	fprintf(out,"SCALARS DEN float 1 \n");
	fprintf(out,"LOOKUP_TABLE default \n");
	
	for(int k=0;k<nz;k++){
		for(int j=0;j<ny;j++){
			for(int i=0;i<nx;i++){
				//ii=(i*ny+j)*nz+k;
				//fprintf(out,"%14.6e \n",data1[i][j][k]);
				fprintf(out,"%14.6e \n",data1[(i*ny+j)*nz+k]);
			}
		}
	}
	
	fprintf(out,"SCALARS ENEG float 1 \n");
	fprintf(out,"LOOKUP_TABLE default \n");
	
	for(int k=0;k<nz;k++){
		for(int j=0;j<ny;j++){
			for(int i=0;i<nx;i++){
				//ii=i*ny*nz+j*nz+k;
				//fprintf(out,"%14.6e \n",data2[i][j][k]);
				fprintf(out,"%14.6e \n",data2[(i*ny+j)*nz+k]);
			}
		}
	}
	
	fclose(out);
	return;
}