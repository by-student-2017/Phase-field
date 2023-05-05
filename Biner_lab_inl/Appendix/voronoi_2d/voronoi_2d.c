#include <stdio.h> //printf()
#include <stdlib.h> //rand()
#include <math.h> //mod() and -lm
#include <time.h>

/* This program generates random polycrystalline microstructure by using
   Voronoi tessellation for the phase-field models */

void Discrete_Voronoi_diagram_2d();
void write_vtk_grid_values_2D();

int main(){
	
	//open output file
	//----- Assign unit names for the output files
	
	//Output file name for Voronoi tessellation
	FILE  *out=fopen("Voronoi_vertices.out","w");
	
	//Intermediate plot file name
	//FILE *out1=fopen("plot_1.out","w");
	
	//Final graphical output file name
	//FILE *out2=fopen("final_plot.p","w");
	
	/* File name that contains the coordinates of
	   the randomly generated points */
	FILE *out3=fopen("original_points.out","w");
	
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
	//
	double x0=0.0;    //The origin of the coordinate system (x-direction)
	double y0=0.0;    //The origin of the coordinate system (y-direction)
	//
	double dx=1.0;
	double dy=1.0;
	//----- ----- ----- ---- -----
	
	int ij;
	
	//----- ----- ----- ---- ----- ----- ----- ----- ----
	//  generate random points for voronoi tessellation
	//----- ----- ----- ---- ----- ----- ----- ----- ----
	// Randomly generate the x- and y-coordinates of the points
	
	srand(time(NULL));
	
	double *x = (double *)malloc(sizeof(double)*( ngrains*9 ));
	double *y = (double *)malloc(sizeof(double)*( ngrains*9 ));
	
	//generate random value in the unit cell
	for(int igrain=0;igrain<ngrains;igrain++){
		x[igrain] = xmax*( (double)rand()/RAND_MAX );
		y[igrain] = ymax*( (double)rand()/RAND_MAX );
	}
	
	//----- ----- ----- ---- ----- -----
	//  Duplicate the points for symmetry
	//----- ----- ----- ---- ----- -----
	/* Replicate the points for their nine periodic images,
	   for generation of periodic simulation cell in phase-field simulations */
	
	int npoints=ngrains; // unit cell
	int jpoint; // all cells
	
	for(int ipoint=0;ipoint<npoints;ipoint++){
		jpoint = npoints*1 + ipoint;
		x[jpoint]=x[ipoint];
		y[jpoint]=y[ipoint]-ymax;
	}
	//
	for(int ipoint=0;ipoint<npoints;ipoint++){
		jpoint = npoints*2 + ipoint;
		x[jpoint]=x[ipoint]+xmax;
		y[jpoint]=y[ipoint]-ymax;
	}
	//
	for(int ipoint=0;ipoint<npoints;ipoint++){
		jpoint = npoints*3 + ipoint;
		x[jpoint]=x[ipoint]+xmax;
		y[jpoint]=y[ipoint];
	}
	//
	for(int ipoint=0;ipoint<npoints;ipoint++){
		jpoint = npoints*4 + ipoint;
		x[jpoint]=x[ipoint]+xmax;
		y[jpoint]=y[ipoint]+ymax;
	}
	//
	for(int ipoint=0;ipoint<npoints;ipoint++){
		jpoint = npoints*5 + ipoint;
		x[jpoint]=x[ipoint];
		y[jpoint]=y[ipoint]+ymax;
	}
	//
	for(int ipoint=0;ipoint<npoints;ipoint++){
		jpoint = npoints*6 + ipoint;
		x[jpoint]=x[ipoint]-xmax;
		y[jpoint]=y[ipoint]+ymax;
	}
	//
	for(int ipoint=0;ipoint<npoints;ipoint++){
		jpoint = npoints*7 + ipoint;
		x[jpoint]=x[ipoint]-xmax;
		y[jpoint]=y[ipoint];
	}
	//
	for(int ipoint=0;ipoint<npoints;ipoint++){
		jpoint = npoints*8 + ipoint;
		x[jpoint]=x[ipoint]-xmax;
		y[jpoint]=y[ipoint]-ymax;
	}
	
	//----- ----- ----- ---- ----- -----
	// Print-out original random points
	//----- ----- ----- ---- ----- -----
	//Output the coordinates of the points to file
	
	for(int jpoint=0;jpoint<npoints*9;jpoint++){
		fprintf(out3,"%4.6e %14.6e \n",x[jpoint],y[jpoint]);
	}
	
	//----- ----- ----- ---- ----- -----
	// generate voronoi diagram in the unit cell
	//----- ----- ----- ---- ----- -----
	
	int Nx=(int)(xmax-x0)/dx;
	int Ny=(int)(ymax-y0)/dy;
	//
	double *cu = (double *)malloc(sizeof(double)*( Nx*Ny ));
	double *fu = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//
	//calculate voronoi diagram of the unit cell
	Discrete_Voronoi_diagram_2d(npoints, x, y,
		Nx, Ny, cu, fu,
		dx, dy);
	
	//----- ----- ----- ---- ----- -----
	// Add vertex markers (vu[Nx][Ny])
	//----- ----- ----- ---- ----- -----
	
	double *vu = (double *)malloc(sizeof(double)*( Nx*Ny ));
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			ij=i*Ny+j;
			vu[ij]=fu[ij];
			for(int igrain=0;igrain<ngrains;igrain++){
				if( (fabs(dx*i - x[igrain]) < dx)
				 && (fabs(dy*j - y[igrain]) < dy) ){
				 	vu[ij]=ngrains*1.2;
				}
			}
		}
	}
	
	//----- ----- ----- ---- ----- -----
	// Duplicate the voronoi diagram
	//----- ----- ----- ---- ----- -----
	
	int Nx3=Nx*3;
	int Ny3=Ny*3;
	double *c = (double *)malloc(sizeof(double)*( Nx3*Ny3 ));
	double *f = (double *)malloc(sizeof(double)*( Nx3*Ny3 ));
	double *v = (double *)malloc(sizeof(double)*( Nx3*Ny3 ));
	//
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			ij=i*Ny+j;
			//----- ----- ----- ---- ----- -----
			//Duplicate the c and the f
			//-:0, non:1, +:2
			//----- ----- ----- ---- ----- -----#0
			c[(i+Nx*1)*Ny3+(j+Ny*1)]=cu[ij]; //non,non
			f[(i+Nx*1)*Ny3+(j+Ny*1)]=fu[ij]; //non,non
			v[(i+Nx*1)*Ny3+(j+Ny*1)]=vu[ij]; //non,non
			//----- ----- ----- ---- ----- -----
			//----- ----- ----- ---- ----- -----#1
			c[(i+Nx*1)*Ny3+(j+Ny*0)]=cu[ij]; //non,-ymax
			f[(i+Nx*1)*Ny3+(j+Ny*0)]=fu[ij]; //non,-ymax
			v[(i+Nx*1)*Ny3+(j+Ny*0)]=vu[ij]; //non,-ymax
			//----- ----- ----- ---- ----- -----#2
			c[(i+Nx*2)*Ny3+(j+Ny*0)]=cu[ij]; //+xmax,-ymax
			f[(i+Nx*2)*Ny3+(j+Ny*0)]=fu[ij]; //+xmax,-ymax
			v[(i+Nx*2)*Ny3+(j+Ny*0)]=vu[ij]; //+xmax,-ymax
			//----- ----- ----- ---- ----- -----#3
			c[(i+Nx*2)*Ny3+(j+Ny*1)]=cu[ij]; //+xmax,non
			f[(i+Nx*2)*Ny3+(j+Ny*1)]=fu[ij]; //+xmax,non
			v[(i+Nx*2)*Ny3+(j+Ny*1)]=vu[ij]; //+xmax,non
			//----- ----- ----- ---- ----- -----#4
			c[(i+Nx*2)*Ny3+(j+Ny*2)]=cu[ij]; //+xmax,+ymax
			f[(i+Nx*2)*Ny3+(j+Ny*2)]=fu[ij]; //+xmax,+ymax
			v[(i+Nx*2)*Ny3+(j+Ny*2)]=vu[ij]; //+xmax,+ymax
			//----- ----- ----- ---- ----- -----#5
			c[(i+Nx*1)*Ny3+(j+Ny*2)]=cu[ij]; //non,+ymax
			f[(i+Nx*1)*Ny3+(j+Ny*2)]=fu[ij]; //non,+ymax
			v[(i+Nx*1)*Ny3+(j+Ny*2)]=vu[ij]; //non,+ymax
			//----- ----- ----- ---- ----- -----#6
			c[(i+Nx*0)*Ny3+(j+Ny*2)]=cu[ij]; //-xmax,+ymax
			f[(i+Nx*0)*Ny3+(j+Ny*2)]=fu[ij]; //-xmax,+ymax
			v[(i+Nx*0)*Ny3+(j+Ny*2)]=vu[ij]; //-xmax,+ymax
			//----- ----- ----- ---- ----- -----#7
			c[(i+Nx*0)*Ny3+(j+Ny*1)]=cu[ij]; //-xmax,non
			f[(i+Nx*0)*Ny3+(j+Ny*1)]=fu[ij]; //-xmax,non
			v[(i+Nx*0)*Ny3+(j+Ny*1)]=vu[ij]; //-xmax,non
			//----- ----- ----- ---- ----- -----#8
			c[(i+Nx*0)*Ny3+(j+Ny*0)]=cu[ij]; //-xmax,-ymax
			f[(i+Nx*0)*Ny3+(j+Ny*0)]=fu[ij]; //-xmax,-ymax
			v[(i+Nx*0)*Ny3+(j+Ny*0)]=vu[ij]; //-xmax,-ymax
			//----- ----- ----- ---- ----- -----
		}
	}
	
	//----- ----- ----- ---- ----- ----- ----- ----- ----- ----- ----- ----- -----
	// print voronoi results to file to be viewed later in Voroni_vertices.out
	//----- ----- ----- ---- ----- ----- ----- ----- ----- ----- ----- ----- -----
	// double lnods[ncount][nnode];
	// int rows = sizeof(lnods)/sizeof(lnods[0]);
	// int colums  = sizeof(lnods[0)/sizeof(lnods[0][0]);
	
	for(int igrain=0;igrain<ngrains;igrain++){
		fprintf(out, "# i %d \n",igrain);
		//
		for(int ipoint=0;ipoint<npoints;ipoint++){
			for(int ni=0;ni<9;ni++){
				jpoint = npoints*ni + ipoint;
				fprintf(out, "%14.6e %14.6e \n",x[jpoint],y[jpoint]);
			}
		}
		//
	}
	
	//----- ----- ----- ---- ----- -----
	// print out as input file
	//----- ----- ----- ---- ----- -----
	
	fprintf(out5, "%5d %5d %5d \n",Nx3,Ny3,ngrains);
	
	for(int i=0;i<Nx3;i++){
		for(int j=0;j<Ny3;j++){
			ij=i*Ny3+j;
			fprintf(out5, "%5d %5d %14.6e %14.6e \n",i,j,f[ij],(1.0-c[ij])); //for etas[i][j][ngrain]=(1.0-f[ij])
		}
	}
	
	//----- ----- ----- ---- ----- -----
	//write vtk file
	/* Write the results in vtk format for contour plots
	   to be viewed by using Paraview */
	int istep=0;
	//----- ----- ----- ---- ----- -----
	//write_vtk_grid_values_2D(Nx,Ny,dx,dy,istep,fu);
	write_vtk_grid_values_2D(Nx,Ny,dx,dy,istep,vu);
	//----- ----- ----- ---- ----- -----
	//write_vtk_grid_values_2D(Nx3,Ny3,dx,dy,istep,f);
	//write_vtk_grid_values_2D(Nx3,Ny3,dx,dy,istep,v);
	//----- ----- ----- ---- ----- -----
}