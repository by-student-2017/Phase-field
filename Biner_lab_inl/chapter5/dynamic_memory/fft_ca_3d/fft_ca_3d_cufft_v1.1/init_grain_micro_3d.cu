/* This function grenerates the order parameters for
   a circular bicrystal or a polycrystalline micro-
   structure consist of 25 grains which is generated
   with the source code given in Appendix A. */

#include <stdio.h> //printf()
#include <stdlib.h> //rand()
#include <math.h> //mod() and -lm
//#include <fftw3.h>

/* Variable and array list
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  Nz: Number of grid points in the z-direction
  dx: Grid spacing between two grid points in the x-direction
  dy: Grid spacing between two grid points in the y-direction
  dz: Grid spacing between two grid points in the z-direction
  iflag: If iflag=1 bicrystal simulation,
            iflag=2 polycrystal simulation
        (in this case, iflag=1 only)
  isolve:(in this case, isolve=1 only)
  ngrain: Number of grains in the simulation
  glist[ngrain]: Flag for existing grains
  etas[Nx][Ny][Nz][ngrain]: common array containing order parameters of grains, isolve=1
*/

void init_grain_micro_3d(int Nx, int Ny, int Nz,
	float dx, float dy, float dz,
	int iflag, int ngrain,
	float *etas, int *glist){
	
	int x0,y0,z0;
	float radius;
	float xlength;
	int ijk;
	
	//----- ----- ----- -----
	//generate two garins
	//----- ----- ----- -----
	if(iflag==1){
		
		//ngrain=2;
		
		//etas[][0] 1st grain
		//etas[][1] 2nd grain
		
		x0=(int)Nx/2;
		y0=(int)Ny/2;
		z0=(int)Nz/2;
		
		radius=14.0; //radius of second grain
		
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					
					etas[ijk*ngrain+0]=1.0;
					etas[ijk*ngrain+1]=0.0;
					
					xlength=sqrt( (i-x0)*(i-x0) + (j-y0)*(j-y0) + (k-z0)*(k-z0) );
					if(xlength <= radius){
						etas[ijk*ngrain+0]=0.0;
						etas[ijk*ngrain+1]=1.0;
					}
				}//end for(k
			}//end for(j
		}//end for(i
		//
	}//end if(iflag
	
	//initialize glist
	for(int igrain=0;igrain<ngrain;igrain++){
		glist[igrain]=1.0;
	}
	
	return;
}