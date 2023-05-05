#include <math.h> //mod() and -lm

/* Variable and array list
  npoints: Number of grains in the simulation cell
  x[npoints]: x position of igrain
  y[npoints]: y position of igrain
  z[npoints]: z position of igrain
  Nx: Number of grid points in the x-direction (Nx=(int)(xmax-x0)/distance_x);)
  Ny: Number of grid points in the y-direction (Ny=(int)(ymax-y0)/distance_y);)
  Nz: Number of grid points in the z-direction (Nz=(int)(zmax-z0)/distance_z);)
  f[Nx][Ny][Nz]: The indices of the Voronoi points
  c[Nx][Ny][Nz]: The points of the Voronoi facets (0:non-boundary, 1:boundary)
*/

void Discrete_Voronoi_diagram_3d(
	int npoints, double *x, double *y, double *z,
	int Nx, int Ny, int Nz, double *c, double *f,
	double distance_x, double distance_y, double distance_z){
	
	int ijk; //ijk=(i*Ny+j)*Nz+k;
	
	//Set the indices of the Voronoi points
	double dx;
	double dy;
	double dz;
	//
	double distance;
	double min;
	double min_start=distance_x*Nx + distance_y*Ny + distance_z*Nz;
	//
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				ijk=(i*Ny+j)*Nz+k;
				//
				min = min_start; //dummy
				f[ijk]=99999;
				for(int ipoint=0;ipoint<npoints;ipoint++){
					//
					for(int tx=-1;tx<=1;tx++){
						for(int ty=-1;ty<=1;ty++){
							for(int tz=-1;tz<=1;tz++){
								dx = distance_x*i - (x[ipoint] + distance_x*Nx*tx);
								dy = distance_y*j - (y[ipoint] + distance_y*Ny*ty);
								dz = distance_z*k - (z[ipoint] + distance_z*Nz*tz);
								distance = sqrt(dx*dx + dy*dy + dz*dz);
								if(distance<=min){
									min=distance;
									f[ijk]=ipoint;
								}
							}//tz
						}//ty
					}//tx
				//
				}//end for(ipoint
			}//end for(k
		}//end for(j
	}//end for(i
	
	//Set boundaries (set the points of the Voronoi facets)
	int ip,im;
	int jp,jm;
	int kp,km;
	//
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				ijk=(i*Ny+j)*Nz+k;
				//
				ip=i+1; if(ip==Nx){ip=0;}
				im=i-1; if(im==-1){im=(Nx-1);}
				jp=j+1; if(jp==Ny){jp=0;}
				jm=j-1; if(jm==-1){jm=(Ny-1);}
				kp=k+1; if(kp==Nz){kp=0;}
				km=k-1; if(km==-1){km=(Nz-1);}
				//
				c[ijk]=0;
				if(  fabs(f[(ip*Ny+j)*Nz+k]-f[ijk]) >= 0.0
				  || fabs(f[(im*Ny+j)*Nz+k]-f[ijk]) >= 0.0
				  || fabs(f[(i*Ny+jp)*Nz+k]-f[ijk]) >= 0.0
				  || fabs(f[(i*Ny+jm)*Nz+k]-f[ijk]) >= 0.0
				  || fabs(f[(i*Ny+j)*Nz+kp]-f[ijk]) >= 0.0
				  || fabs(f[(i*Ny+j)*Nz+km]-f[ijk]) >= 0.0
				){
					c[ijk]=1;
				}
			}//end for(k
		}//end for(j
	}//end for(i
	
	return;
}