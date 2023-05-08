/* This function initializes the concentration and
   the order parameters in the simulation cell and
   defines the particle positions for either
   two particles or nine particle simulations.
   It is used in both longhand(iflag=1). */

#include <math.h> //mod() and -lm
#include <stdlib.h> //rand() and malloc
//#include <fftw3.h>

/* Variable and array list
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  Nz: Number of grid points in the z-direction
  npart: Number of grid points of the simulation
  iflag: iflag=1 for longhand code (iflag=1 only)
  etas[Nx][Ny][Nz][npart]: Common array of non-conserved order
    parameters for the particles. In the optimized code,
    this array is in the form of etas[NxNyNz][npart] where,
    NxNyNz is total number of grid points in the simulation cell.
  con[Nx][Ny][Nz]: Concentration field in takes values of 
    one in the particles and zero elsewhere. In the optimized code,
    this is one-dimensional array in the from of con[NxNyNz].
*/

void micro_sint_pre_3d(int Nx, int Ny, int Nz,
	int npart, int iflag, 
	double *etas, double *con){
	
	int ijk; //ijk=(i*Ny+j)*Nz+k;
	
	//initialize
	if(iflag==1){
		for(int ipart=0;ipart<npart;ipart++){
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					for(int k=0;k<Nz;k++){
						ijk=(i*Ny+j)*Nz+k;
						etas[ijk*npart+ipart]=0.0;
					}
				}
			}
		}
	}//end if
	
	double R;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	//double xc[9];
	double      *xc = (double *)malloc(sizeof(double)*( npart ));
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	//double yc[9];
	double      *yc = (double *)malloc(sizeof(double)*( npart ));
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	double zc;
	
	double Rx;
	double xx1;
	
	// Place nine spherical particles into the simulation cell
	if(npart != 2){
		//
		//The radius of the first five large particles
		R=10.0;
		
		//Center coordinates of the particles
		xc[0]=29.0;
		yc[0]=50.0;
		
		xc[1]=50.0;
		yc[1]=50.0;
		
		xc[2]=71.0;
		yc[2]=50.0;
		
		xc[3]=50.0;
		yc[3]=29.0;
		
		xc[4]=50.0;
		yc[4]=71.0;
		
		xc[5]=39.0;
		yc[5]=39.0;
		
		xc[6]=61.0;
		yc[6]=39.0;
		
		xc[7]=39.0;
		yc[7]=61.0;
		
		xc[8]=61.0;
		yc[8]=61.0;
		
		zc = (int)Nz/2;
		
		//Set the particle radius (first firve are the large ones)
		//for(int ipart=0;ipart<npart;ipart++){
		for(int ipart=0;ipart<9;ipart++){
			//
			Rx=R;
			
			if(ipart>=5){
				Rx=0.5*R;
			}
			
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					for(int k=0;k<Nz;k++){
						ijk=(i*Ny+j)*Nz+k;
						xx1=sqrt( (i-xc[ipart])*(i-xc[ipart])
								 +(j-yc[ipart])*(j-yc[ipart])
							     +(k-zc       )*(k-zc       )
								);
						//
						if(xx1 <= Rx){
							con[ijk]=0.999;
							etas[ijk*npart+ipart]=0.999;
						}
					}
				}
			}
			//
		}//end for(ipart
	}//end if
	
	double R1,R2;
	double x1,y1,y2,z1;
	double xx2;
	
	//Place two particles into the simulation cell
	if(npart==2){
		//
		//Radius of the particles
		R1=20.0;
		R2=0.5*R1;
		
		//Center corrdinates of the particles
		x1=(int)Nx/2;
		y1=40.0;
		y2=70.0;
		z1=(int)Nz/2;
		
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					//
					xx1=sqrt( (i-x1)*(i-x1) + (j-y1)*(j-y1) + (k-z1)*(j-z1) );
					xx2=sqrt( (i-x1)*(i-x1) + (j-y2)*(j-y2) + (k-z1)*(j-z1) );
					
					if(xx1 <= R1){
						con[ijk]=0.999;
						etas[ijk*npart+0]=0.999;
					}
					
					//Repeat the same procedure for the second particle
					if(xx2 <= R2){
						con[ijk]=0.999;
						etas[ijk*npart+0]=0.0;
						etas[ijk*npart+1]=0.999;
					}
				}
			}
		}
		//
	}//end if
	
	free(xc);
	free(yc);
	
	return;
}