/* This function introduces randomly distributed and
   slightly overlapping cells having radium of R into 
   the simulation cell. */

/* Variable and array list
  Nx: Number of grid points of the simulation cell in the x-direction
  Ny: Number of grid points of the simulation cell in the y-direction
  R: The radius of the cells in terms of number of grid points
  ncell: Input as desired cell number, output is the actual number of cell generated
  phis[Nx][Ny][ncell]: Common array containing the order parameters of the cells
  vac[ncell]: The self-propulsion values of the cells
  nccell: Number of soft cells
  ccell[nccell]: Cell number of soft cells
*/

#include <math.h> //M_PI
#include <stdlib.h> //rand()
#include <stdio.h> //printf()

void micro_poly_cell_2d(int Nx, int Ny, double R,
	int ncell, double *phis, double *vac,
	int *nccell, int *ccell){
	
	int max_ncell=ncell;
	int used_ncell=ncell;
	
	int icell;
	int ij; //ij=(i*Ny+j);
	int ijci; //ijci=(ij*max_ncell+icell);
	
	//initialize (Longhand initialization)
	// Initialize the order parameter of the cells
	for(icell=0;icell<max_ncell;icell++){
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				//phis[i][j][icell]=0.0;
				ij=(i*Ny+j);
				ijci=(ij*max_ncell+icell);
				phis[ijci]=0.0;
			}
		}
	}
	
	//----- ----- ----- -----
	int iter;
	//----- ----- ----- -----
	double R2, Rsq;
	double xmin, ymin;
	double xmax, ymax;
	//----- ----- ----- -----
	//max_ncell=ncell=80
	//double xc_1d[80];
	double *xc_1d = (double *)malloc(sizeof(double)*( max_ncell ));
	//double yc_1d[80];
	double *yc_1d = (double *)malloc(sizeof(double)*( max_ncell ));
	//----- ----- ----- -----
	int ncol=1;
	int nrow=2;
	//double xc_2d[ncol][nrow];
	double *xc_2d = (double *)malloc(sizeof(double)*( ncol*nrow ));
	//double yc_2d[ncol][nrow];
	double *yc_2d = (double *)malloc(sizeof(double)*( ncol*nrow ));
	//----- ----- ----- -----
	double xnc;
	double ync;
	//----- ----- ----- -----
	int iflag;
	//----- ----- ----- -----
	double xdist;
	//----- ----- ----- -----
	double ix;
	//----- ----- ----- -----
	double dx;
	double dy;
	//----- ----- ----- -----
	double velc;
	//----- ----- ----- -----
	
	/* Generate slightly overlapping cell
	   microstructure, if ncell>2 */
	if(max_ncell>2){
		
		// Calculate the diameter and the square of cell radius, R
		R2 = 2.0*R;
		Rsq = R*R;
		
		// Initialize the total number cells in the simulation
		used_ncell = 0;
		
		/* Determine the minimum and maximum cell dimensions in
		   the x and y directions */
		xmin = 0.0;
		ymin = 0.0;
		
		xmax=(double)Nx;
		ymax=(double)Ny;
		
		// Initialize the coordinates of the cell centeres
		for(int i=0;i<max_ncell;i++){
			xc_1d[i]=0.0;
			yc_1d[i]=0.0;
		}
		
		// Start randomly introducing cells into the simulation cell
		for(iter=0;iter<500000;iter++){
			
			// Randomly assign cell center coordinates
			xnc=(double)Nx*( (double)rand()/RAND_MAX );
			ync=(double)Ny*( (double)rand()/RAND_MAX );
			
			/* Set iflag=1. In the following tests, if its value
			   still remains as iflag=1, a new cell will be generated. */
			iflag=1;
			
			/* Check, if new cell as a whole is in the bounding box of 
			   the simulation cell. if not set flag=0. */
			if( ((xnc-R) < xmin) || ((xnc+R) > xmax) ){
				iflag=0;
			}
			
			/* Check, if new cell as a whole is in the bounding box of 
			   the simulation cell. if not set flag=0. */
			if( ((ync-R) < ymin) || ((ync+R) > ymax) ){
				iflag=0;
			}
			
			// If new cell fits into the simulation cell continue
			if(iflag==1){
				/* Calculate the distance to the previously generated
				   cell centers. If distance<1.6R set iflag=0. */
				for(int i=0;i<used_ncell;i++){
					xdist=sqrt((xc_1d[i]-xnc)*(xc_1d[i]-xnc)
							  +(yc_1d[i]-ync)*(yc_1d[i]-ync)
							  );
					if(xdist <= R*1.6){
						iflag=0;
					}
				}
			}//end if
			
			/* After two tests above, if iflag=1, set the 
			   order parameter of the new cell to 0.999. */
			if(iflag==1){
				
				used_ncell = used_ncell + 1;
				
				xc_1d[(used_ncell-1)] = xnc;
				yc_1d[(used_ncell-1)] = ync;
				
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						if( ((i-xnc)*(i-xnc) + (j-ync)*(j-ync))<Rsq ){
							ij=i*Ny+j;
							phis[ij*max_ncell+(used_ncell-1)]=0.999;
						}
					}
				}
			}//end iflag
			
			/* If the total number of cell equals to 80 exit
			   from the loop irand. */
			if(used_ncell==max_ncell){
				break;
			}
		}//end for(iter)
		
		//randomize cell self propulsion velocities
		// Randomly change the sign of self-propulsion values of the cells
		//Absolute value of the self-propulsion value of the cells
		velc=0.2;
		
		for(int i=0;i<used_ncell;i++){
			ix = (double)rand()/RAND_MAX;
			if(ix<=0.5){
				vac[i]=-velc;
			}else{
				vac[i]=+velc;
			}
		}
		
		/* Print how many cells are generated after
		   500000 random steps to screen. This number can be
		   set to larger value if more cell needs to be generated. */
		printf("iteration done: %6d \n",iter);
		printf("number of cell created: %5d \n",used_ncell);
		nccell[1]=used_ncell;
		
		//soft cells
		// Assign the following cell numbers as soft cell
		nccell[0]=5;
		ccell[0]=32;
		ccell[1]=11;
		ccell[2]=16;
		ccell[3]=21;
		ccell[4]=46;
	}//end if(ncell==80
	
	//For two cell simulation
	/* If ncell=2, place two cells into the center of
	   the simulation cell with having radius R and
	   assign their self-propulsion values. */
	if(ncell==2){
		
		nccell[1]=2;
		nccell[0]=1;
		ccell[0]=1;
		
		R2=R*R;
		
		xc_2d[0*nrow+0]=(double)Nx/2.0-1.25*R;
		yc_2d[0*nrow+0]=(double)Ny/2.0;
		
		xc_2d[0*nrow+1]=(double)Nx/2.0+1.25*R;
		yc_2d[0*nrow+1]=(double)Ny/2.0;
		
		//ncol=1;
		//nrow=2;
		
		icell=0;
		
		for(int icol=0;icol<ncol;icol++){
			for(int irow=0;irow<nrow;irow++){
				//
				icell=icell+1;
				dx=xc_2d[icol*nrow+irow];
				dy=yc_2d[icol*nrow+irow];
				//
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						if( ((i-dx)*(i-dx) + (j-dy)*(j-dy))<R2 ){
							ij=i*Ny+j;
							phis[ij*max_ncell+(icell-1)]=0.999;
						}
					}
				}
				//
			}
		}
		
		//cell self propulsion velocity
		vac[0]=+0.5;
		vac[1]=-0.5;
	}//end if(ncell==2
	
	return;
	
	free(xc_1d);
	free(yc_1d);
	free(xc_2d);
	free(yc_2d);
}