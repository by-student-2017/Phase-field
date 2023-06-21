#include <stdio.h>
#include <stdlib.h> //rand()
#include <math.h> //M_PI

void write_vtk_grid_values_2D(int nx, int ny, 
	float dx, float dy, int istep, 
	float *data1, float *data2, float *data3, float *data4){
	
	//open output file
	char fname[256];
	sprintf(fname,"time_%d.vtk",istep);
	FILE *out=fopen(fname,"w");
	
	int nz=1;
	int npoin=nx*ny*nz;
	float x,y,z;
		
	float data0;
	
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
	fprintf(out,"POINT_DATA %5d \n",npoin);
	fprintf(out,"SCALARS C_Cu float 1 \n");
	fprintf(out,"LOOKUP_TABLE default \n");
	
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			//ii=i*ny+j;
			//fprintf(out,"%14.6e \n",data1[i][j]);
			fprintf(out,"%14.6e \n",data1[i*ny+j]);
		}
	}
	
	fprintf(out,"SCALARS C_Mn float 1 \n");
	fprintf(out,"LOOKUP_TABLE default \n");
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			fprintf(out,"%14.6e \n",data2[i*ny+j]);
		}
	}
	
	fprintf(out,"SCALARS C_Ni float 1 \n");
	fprintf(out,"LOOKUP_TABLE default \n");
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			fprintf(out,"%14.6e \n",data3[i*ny+j]);
		}
	}
	
	fprintf(out,"SCALARS Order_param float 1 \n");
	fprintf(out,"LOOKUP_TABLE default \n");
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			fprintf(out,"%14.6e \n",data4[i*ny+j]);
		}
	}
	
	fprintf(out,"SCALARS C_Fe float 1 \n");
	fprintf(out,"LOOKUP_TABLE default \n");
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			data0=(1.0-data2[i*ny+j]-data2[i*ny+j]-data3[i*ny+j]);
			if(data0<0){
				data0=0.0;
			}
			fprintf(out,"%14.6e \n",data0);
		}
	}
	
	fclose(out);
	return;
}