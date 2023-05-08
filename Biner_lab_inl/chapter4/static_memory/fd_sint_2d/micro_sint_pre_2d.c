/* This function initializes the concentration and
   the order parameters in the simulation cell and
   defines the particle positions for either
   two particles or nine particle simulations.
   It is used in both longhand(iflag=1). */

#include <math.h> //mod() and -lm
//#include <stdlib.h> //rand()
//#include <fftw3.h>

/* Variable and array list
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  npart: Number of grid points of the simulation
  iflag: iflag=1 for longhand code (iflag=1 only)
  etas[Nx][Ny][npart]: Common array of non-conserved order
    parameters for the particles. In the optimized code,
    this array is in the form of etas[NxNy][npart] where,
    NxNy is total number of grid points in
    the simulation cell.
  con[Nx][Ny]: Concentration field in takes values of 
    one in the particles and zero elsewhere. In the optimized code,
    this is one-dimensional array in the from of con[NxNy].
*/

void micro_sint_pre_2d(int Nx, int Ny, 
	int npart, int iflag, 
	double *etas, double *con){
	
	int ij; //ij=i*Ny+j;
	
	//initialize
	if(iflag==1){
		for(int ipart=0;ipart<npart;ipart++){
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					ij=i*Ny+j;
					etas[ij*npart+ipart]=0.0;
				}
			}
		}
	}
	
	double R;
	double xc[9],yc[9];
	
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
					ij=i*Ny+j;
					xx1=sqrt( (i-xc[ipart])*(i-xc[ipart])
							 +(j-yc[ipart])*(j-yc[ipart]));
					//
					if(xx1 <= Rx){
						con[ij]=0.999;
						etas[ij*npart+ipart]=0.999;
					}
				}
			}
			//
		}//end for(ipart
	}//end if
	
	double R1,R2;
	double x1,y1,y2;
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
		
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ij=i*Ny+j;
				//
				xx1=sqrt( (i-x1)*(i-x1) + (j-y1)*(j-y1) );
				xx2=sqrt( (i-x1)*(i-x1) + (j-y2)*(j-y2) );
				
				if(xx1 <= R1){
					con[ij]=0.999;
					etas[ij*npart+0]=0.999;
				}
				
				//Repeat the same procedure for the second particle
				if(xx2 <= R2){
					con[ij]=0.999;
					etas[ij*npart+0]=0.0;
					etas[ij*npart+1]=0.999;
				}
			}
		}
		//
	}//end if
	
	return;
}