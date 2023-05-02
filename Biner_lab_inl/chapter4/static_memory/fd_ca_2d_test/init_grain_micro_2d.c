/* This function grenerates the order parameters for
   a circular bicrystal or a polycrystalline micro-
   structure consist of 25 grains which is generated
   with the source code given in Appendix A. */

#include <stdlib.h> //rand()
#include <math.h> //mod() and -lm
//#include <fftw3.h>

/* Variable and array list
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  dx: Grid spacing between two grid points in the x-direction
  dy: Grid spacing between two grid points in the y-direction
  iflag: If iflag=1 bicrystal simulation,
            iflag=2 polycrystal simulation
        (in this case, iflag=1 only)
  isolve:(in this case, isolve=1 only)
  ngrain: Number of grains in the simulation
  glist[ngrain]: Flag for existing grains
  etas[Nx][Ny][ngrain]: common array containing order parameters of grains, isolve=1
*/

void init_grain_micro_2d(int Nx, int Ny,
	double dx, double dy, 
	int iflag, int ngrain,
	double *etas, int *glist){
	
	int x0,y0;
	double radius;
	double xlength;
	int ij;
	
	//----- ----- ----- -----
	//generate two garins
	//----- ----- ----- -----
	if(iflag==1){
		
		//ngrain=2;
		
		//etas[][0] 1st grain
		//etas[][1] 2nd grain
		
		x0=(int)Nx/2;
		y0=(int)Ny/2;
		
		radius=14.0; //radius of second grain
		
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ij=i*Ny+j;
				
				etas[ij*ngrain+0]=1.0;
				etas[ij*ngrain+1]=0.0;
				
				xlength=sqrt( (i-x0)*(i-x0) + (j-y0)*(j-y0) );
				if(xlength <= radius){
					etas[ij*ngrain+0]=0.0;
					etas[ij*ngrain+1]=1.0;
				}
			}//end for(j
		}//end for(i
		//
	}//end if(iflag
	
	//----- ----- ----- -----
	//generate polycrystal microstructure
	//----- ----- ----- -----
	//if(iflag==2){
	//	FILE *in=fopen("grain_25.inp","r");
	//}
	
	//initialize glist
	for(int igrain=0;igrain<ngrain;igrain++){
		glist[igrain]=1.0;
	}
	
	return;
}