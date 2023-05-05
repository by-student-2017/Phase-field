#include <stdio.h> //printf()
#include <stdlib.h> //rand()
#include <math.h> //mod() and -lm
#include <time.h>

/* This program generates random polycrystalline microstructure by using
   Voronoi tessellation for the phase-field models */

void Discrete_Voronoi_diagram_3d();
void write_vtk_grid_values_3D();

int main(){
	
	//open output file
	//----- Assign unit names for the output files
	
	//Output file name for Voronoi tessellation
	//FILE  *out=fopen("Voronoi_vertices.out","w");
	
	//Intermediate plot file name
	//FILE *out1=fopen("plot_1.out","w");
	
	//Final graphical output file name
	//FILE *out2=fopen("final_plot.p","w");
	
	/* File name that contains the coordinates of
	   the randomly generated points */
	//FILE *out3=fopen("original_points.out","w");
	
	/* File name that conatins the Cartesian coordinates of
	   the corners of the simulation cell */
	//FILE *out4=fopen("cell_1.out","w");
	
	/* Tabulated output file that will be used in
	   phase-field simulations */
	FILE *out5=fopen("grain_25.inp","w");
	//----- ----- ----- ---- -----
	
	//----- ----- ----- ---- -----
	int ngrains=25;   //Number of grains in the simulation cell
	//
	double xmax=64.0; //The length of the simulation cell in the x-direction
	double ymax=64.0; //The length of the simulation cell in the y-direction
	double zmax=2.0;  //The length of the simulation cell in the z-direction
	//
	double x0=0.0;    //The origin of the coordinate system (x-direction)
	double y0=0.0;    //The origin of the coordinate system (y-direction)
	double z0=0.0;    //The origin of the coordinate system (y-direction)
	//
	double dx=1.0;
	double dy=1.0;
	double dz=1.0;
	//
	int ntimes=3;     //Number of duplicate (x-,y-,and,z-direction)
	//----- ----- ----- ---- -----
	
	int ijk;
	
	//----- ----- ----- ---- ----- ----- ----- ----- ----
	//  generate random points for voronoi tessellation
	//----- ----- ----- ---- ----- ----- ----- ----- ----
	// Randomly generate the x- and y-coordinates of the points
	
	srand(time(NULL));
	
	double *x = (double *)malloc(sizeof(double)*( ngrains*(ntimes * ntimes * ntimes) ));
	double *y = (double *)malloc(sizeof(double)*( ngrains*(ntimes * ntimes * ntimes) ));
	double *z = (double *)malloc(sizeof(double)*( ngrains*(ntimes * ntimes * ntimes) ));
	
	//generate random value in the unit cell
	for(int igrain=0;igrain<ngrains;igrain++){
		x[igrain] = xmax*( (double)rand()/RAND_MAX );
		y[igrain] = ymax*( (double)rand()/RAND_MAX );
		z[igrain] = zmax*( (double)rand()/RAND_MAX );
	}
	
	int npoints=ngrains; // unit cell
	
	//----- ----- ----- ---- ----- -----
	// generate voronoi diagram in the unit cell
	//----- ----- ----- ---- ----- -----
	
	int Nx=(int)(xmax-x0)/dx;
	int Ny=(int)(ymax-y0)/dy;
	int Nz=(int)(zmax-z0)/dz;
	//
	double *cu = (double *)malloc(sizeof(double)*( Nx*Ny*Nz ));
	int *fu = (int *)malloc(sizeof(int)*( Nx*Ny*Nz ));
	//
	//calculate voronoi diagram of the unit cell
	Discrete_Voronoi_diagram_3d(npoints, x, y, z,
		Nx, Ny, Nz, cu, fu,
		dx, dy, dz);
	
	//----- ----- ----- ---- ----- -----
	// Add vertex markers (vu[Nx][Ny][Nz])
	//----- ----- ----- ---- ----- -----
	
	int *vu = (int *)malloc(sizeof(int)*( Nx*Ny*Nz ));
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				ijk=(i*Ny+j)*Nz+k;
				vu[ijk]=fu[ijk];
				for(int igrain=0;igrain<ngrains;igrain++){
					if( (fabs(dx*i - x[igrain]) < dx)
					 && (fabs(dy*j - y[igrain]) < dy)
					 && (fabs(dz*k - z[igrain]) < dz) ){
					 	vu[ijk]=ngrains+(int)(ngrains/5);
					}
				}//for (igrain
			}//for(k
		}//for(j
	}//for(i
	
	//----- ----- ----- ---- ----- -----
	// Duplicate the voronoi diagram
	//----- ----- ----- ---- ----- -----
	
	int Nxn=Nx*ntimes;
	int Nyn=Ny*ntimes;
	int Nzn=Nz*ntimes;
	double *c = (double *)malloc(sizeof(double)*( Nxn*Nyn*Nzn ));
	int *f = (int *)malloc(sizeof(int)*( Nxn*Nyn*Nzn ));
	int *v = (int *)malloc(sizeof(int)*( Nxn*Nyn*Nzn ));
	//
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				ijk=(i*Ny+j)*Nz+k;
				//Duplicate the c and the f
				for(int tx=0;tx<ntimes;tx++){
					for(int ty=0;ty<ntimes;ty++){
						for(int tz=0;tz<ntimes;tz++){
							c[((i+Nx*tx)*Nyn+(j+Ny*ty))*Nzn+(k+Nz*tz)]=cu[ijk];
							f[((i+Nx*tx)*Nyn+(j+Ny*ty))*Nzn+(k+Nz*tz)]=fu[ijk];
							v[((i+Nx*tx)*Nyn+(j+Ny*ty))*Nzn+(k+Nz*tz)]=vu[ijk];
						}//tz
					}//ty
				}//tx
				//
			}//for(k
		}//for(j
	}//for(i
	
	//----- ----- ----- ---- ----- -----
	// print out as input file
	//----- ----- ----- ---- ----- -----
	
	fprintf(out5, "%5d %5d %5d %5d \n",Nxn,Nyn,Nzn,ngrains);
	
	for(int i=0;i<Nxn;i++){
		for(int j=0;j<Nyn;j++){
			for(int k=0;k<Nzn;k++){
				ijk=(i*Nyn+j)*Nzn+k;
				fprintf(out5, "%5d %5d %5d %5d %14.6e \n",i,j,k,f[ijk],(1.0-c[ijk]));
			}
		}
	}
	
	//----- ----- ----- ---- ----- -----
	//write vtk file
	/* Write the results in vtk format for contour plots
	   to be viewed by using Paraview */
	int istep=0;
	//----- ----- ----- ---- ----- -----
	//write_vtk_grid_values_3D(Nx,Ny,Nz,dx,dy,dz,istep,fu);
	write_vtk_grid_values_3D(Nx,Ny,Nz,dx,dy,dz,istep,vu);
	//----- ----- ----- ---- ----- -----
	//write_vtk_grid_values_3D(Nxn,Nyn,Nzn,dx,dy,dz,istep,f);
	//write_vtk_grid_values_3D(Nxn,Nyn,Nzn,dx,dy,dz,istep,v);
	//----- ----- ----- ---- ----- -----
}