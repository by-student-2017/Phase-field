#include <math.h> //mod() and -lm

/* Variable and array list
  npoints: Number of grains in the simulation cell
  x[npoints]: x position of igrain
  y[npoints]: y position of igrain
  Nx: Number of grid points in the x-direction (Nx=(int)(xmax-x0)/distance_x);)
  Ny: Number of grid points in the y-direction (Ny=(int)(ymax-y0)/distance_y);)
  f[Nx][Ny]: The indices of the Voronoi points
  c[Nx][Ny]: The points of the Voronoi facets (0:non-boundary, 1:boundary)
*/

void Discrete_Voronoi_diagram_2d(
	int npoints, double *x, double *y,
	int Nx, int Ny, double *c, double *f,
	double distance_x, double distance_y){
	
	int ij; //ij=i*Ny+j;
	
	//Set the indices of the Voronoi points
	double dx;
	double dy;
	double distance;
	double min;
	//
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			ij=i*Ny+j;
			//
			min = Nx + Ny; //dummy
			f[ij]=99999;
			for(int ipoint=0;ipoint<npoints;ipoint++){
				//
				//----- ----- ----- ---- ----- -----#0
				// main area
				dx = distance_x*i - x[ipoint];
				dy = distance_y*j - y[ipoint];
				distance = sqrt(dx*dx + dy*dy);
				if(distance<=min){
					min=distance;
					f[ij]=ipoint;
				}
				//----- ----- ----- ---- ----- -----#1
				// periodic boundary conditions
				dx = distance_x*i - (x[ipoint]);
				dy = distance_y*j - (y[ipoint]-Ny);
				distance = sqrt(dx*dx + dy*dy);
				if(distance<=min){
					min=distance;
					f[ij]=ipoint;
				}
				//----- ----- ----- ---- ----- -----#2
				// periodic boundary conditions
				dx = distance_x*i - (x[ipoint]+Nx);
				dy = distance_y*j - (y[ipoint]-Ny);
				distance = sqrt(dx*dx + dy*dy);
				if(distance<=min){
					min=distance;
					f[ij]=ipoint;
				}
				//----- ----- ----- ---- ----- -----#3
				// periodic boundary conditions
				dx = distance_x*i - (x[ipoint]+Nx);
				dy = distance_y*j - (y[ipoint]);
				distance = sqrt(dx*dx + dy*dy);
				if(distance<=min){
					min=distance;
					f[ij]=ipoint;
				}
				//----- ----- ----- ---- ----- -----#4
				// periodic boundary conditions
				dx = distance_x*i - (x[ipoint]+Nx);
				dy = distance_y*j - (y[ipoint]+Ny);
				distance = sqrt(dx*dx + dy*dy);
				if(distance<=min){
					min=distance;
					f[ij]=ipoint;
				}
				//----- ----- ----- ---- ----- -----#5
				// periodic boundary conditions
				dx = distance_x*i - (x[ipoint]);
				dy = distance_y*j - (y[ipoint]+Ny);
				distance = sqrt(dx*dx + dy*dy);
				if(distance<=min){
					min=distance;
					f[ij]=ipoint;
				}
				//----- ----- ----- ---- ----- -----#6
				// periodic boundary conditions
				dx = distance_x*i - (x[ipoint]-Nx);
				dy = distance_y*j - (y[ipoint]+Ny);
				distance = sqrt(dx*dx + dy*dy);
				if(distance<=min){
					min=distance;
					f[ij]=ipoint;
				}
				//----- ----- ----- ---- ----- -----#7
				// periodic boundary conditions
				dx = distance_x*i - (x[ipoint]-Nx);
				dy = distance_y*j - (y[ipoint]);
				distance = sqrt(dx*dx + dy*dy);
				if(distance<=min){
					min=distance;
					f[ij]=ipoint;
				}
				//----- ----- ----- ---- ----- -----#8
				// periodic boundary conditions
				dx = distance_x*i - (x[ipoint]-Nx);
				dy = distance_y*j - (y[ipoint]-Ny);
				distance = sqrt(dx*dx + dy*dy);
				if(distance<=min){
					min=distance;
					f[ij]=ipoint;
				}
				//
			}//end for(ipoint
		}//end for(j
	}//end for(i
	
	//Set boundaries (set the points of the Voronoi facets)
	int ip,im;
	int jp,jm;
	//
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			ij=i*Ny+j;
			//
			ip=i+1; if(ip==Nx){ip=0;}
			im=i-1; if(im==-1){im=(Nx-1);}
			jp=j+1; if(jp==Ny){jp=0;}
			jm=j-1; if(jm==-1){jm=(Ny-1);}
			//
			c[ij]=0;
			if(  fabs(f[ip*Ny+j]-f[ij]) >= 0.0
			  || fabs(f[im*Ny+j]-f[ij]) >= 0.0
			  || fabs(f[i*Ny+jp]-f[ij]) >= 0.0
			  || fabs(f[i*Ny+jm]-f[ij]) >= 0.0
			){
				c[ij]=1;
			}
		}//end for(j
	}//end for(i
	
	return;
}