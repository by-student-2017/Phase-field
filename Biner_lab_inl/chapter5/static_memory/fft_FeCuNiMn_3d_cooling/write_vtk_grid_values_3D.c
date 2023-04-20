#include <stdio.h>
#include <stdlib.h> //rand()
#include <math.h> //M_PI

void write_vtk_grid_values_3D(int nx, int ny, int nz, 
	double dx, double dy, double dz, int istep, 
	double *data1, double *data2, double *data3, double *data4){
	
	//open output file
	char fname[256];
	sprintf(fname,"time_%d.vtk",istep);
	FILE *out=fopen(fname,"w");
	
	int npoin=nx*ny*nz;
	double x,y,z;
	
	int ii;
	double data0;
	
	//start writing ASCII VTK file
	
	//header of VTK file
	fprintf(out,"# vtk DataFile Version 2.0 \n");
	fprintf(out,"time_10.vtk \n");
	fprintf(out,"ASCII \n");
	fprintf(out,"DATASET STRUCTURED_GRID \n");
	
	//coords of grid points
	fprintf(out,"DIMENSIONS %5d %5d %5d \n",nx,ny,nz);
	fprintf(out,"POINTS %15d float \n",npoin);
	
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			for(int k=0;k<nz;k++){
				x=dx*i;
				y=dy*j;
				z=dz*k;
				
				fprintf(out,"%14.6e %14.6e %14.6e \n",x,y,z);
			}
		}
	}
	
	//write grid point values
	fprintf(out,"POINT_DATA %15d \n",npoin);
	fprintf(out,"SCALARS Concentration_Cu float 1 \n");
	fprintf(out,"LOOKUP_TABLE default \n");
	
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			for(int k=0;k<nz;k++){
				//ii=i*ny*nz+j*nz+k;
				//fprintf(out,"%14.6e \n",data1[i][j][k]);
				fprintf(out,"%14.6e \n",data1[i*ny*nz+j*nz+k]);
			}
		}
	}
	
	fprintf(out,"SCALARS Concentration_Mn float 1 \n");
	fprintf(out,"LOOKUP_TABLE default \n");
	
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			for(int k=0;k<nz;k++){
				fprintf(out,"%14.6e \n",data2[i*ny*nz+j*nz+k]);
			}
		}
	}
	
	
	fprintf(out,"SCALARS Concentration_Ni float 1 \n");
	fprintf(out,"LOOKUP_TABLE default \n");
	
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			for(int k=0;k<nz;k++){
				fprintf(out,"%14.6e \n",data3[i*ny*nz+j*nz+k]);
			}
		}
	}
	
	fprintf(out,"SCALARS Order_parameter float 1 \n");
	fprintf(out,"LOOKUP_TABLE default \n");
	
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			for(int k=0;k<nz;k++){
				fprintf(out,"%14.6e \n",data4[i*ny*nz+j*nz+k]);
			}
		}
	}
	
	fprintf(out,"SCALARS Concentration_Fe float 1 \n");
	fprintf(out,"LOOKUP_TABLE default \n");
	for(int i=0;i<nx;i++){
		for(int j=0;j<ny;j++){
			for(int k=0;k<nz;k++){
				ii=i*ny*nz+j*nz+k;
				data0=(1.0-data2[ii]-data2[ii]-data3[ii]);
				if(data0<0){
					data0=0.0;
				}
				fprintf(out,"%14.6e \n",data0);
			}
		}
	}
	
	fclose(out);
	return;
}