/* Finite-difference phase-field code for
   solving Cahn-Hilliard equation */

/* This program solves the Cahn-Hilliard phase-field
   equation with finite difference algorithm using
   five-point stencil. The time integration is
   carried out using explicit Euler scheme. */

#include <stdio.h> //printf()
#include <stdlib.h> //rand()
#include <math.h> //mod() and -lm
#include <time.h>

void micro_ch_pre_3d();
double calculate_energy_3d();
double free_energy_ch_3d();
void write_vtk_grid_values_3D();

int main(){
	// Get initial wall clock time beginning of the execution
	clock_t start, end;
	double compute_time;
	//get initial wall time
	start = clock();
	
	//open an output file for writing total bulk energy values
	FILE *out2=fopen("time_energy.out","w");
	
	//simulation cell parameters
	int Nx=64; //Number of grid points in the x-direction
	int Ny=64; //Number of grid points in the y-direction
	int Nz=2; //Number of grid points in the y-direction
	int NxNyNz=Nx*Ny*Nz; //Total number of grid points in the simulation cell
	
	//The distance between two grid points in x,y-direction
	double dx=1.0; //Grid spacing between two grid pints in x-direction
	double dy=1.0; //Grid spacing between two grid pints in y-direction
	double dz=1.0; //Grid spacing between two grid pints in z-direction
	
	//time integration parameters
	int nstep=10000; //Number of time integration steps
	int nprint=50;  //Output frequency to write the results to file
	double dtime=1.e-2; //Time increment for the numerical integration
	double ttime=0.0;    //Total time
	
	//material specific parameters
	
	//Initial concentrations of alloying elements
	double c0=0.40;       //Average composition
	double mobility=1.0;  //The value of mobility coefficient (dimensionless)
	double grad_coef=0.5; //The value of gradient energy coefficients [J(nm)^2/mol]
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	//double       con[Nx][Ny][Nz]; //concentraion
	double       *con = (double *)malloc(sizeof(double)*( NxNyNz ));
	//double   lap_con[Nx][Ny][Nz]; //Laplacian
	double   *lap_con = (double *)malloc(sizeof(double)*( NxNyNz ));
	//double     dummy[Nx][Ny][Nz];
	double     *dummy = (double *)malloc(sizeof(double)*( NxNyNz ));
	//double lap_dummy[Nx][Ny][Nz]; //Laplacian
	double *lap_dummy = (double *)malloc(sizeof(double)*( NxNyNz ));
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	int ijk; //ijk=(i*Ny+j)*Nz+k;
	//----- ----- ----- -----
	
	//prepare microstructure
	// Initialize the concentration filed with random modulation
	micro_ch_pre_3d(Nx,Ny,Nz,c0,con);
	
	//----- ----- ----- -----
	int ip,im;
	int jp,jm;
	int kp,km;
	//
	double hne,hnw;
	double hns,hnn;
	double hnu,hnd;
	double hnc;
	//----- ----- ----- -----
	double dfdcon;
	//----- ----- ----- -----
	double energy;
	//----- ----- ----- -----
	
	//Time evolution of concentration filed
	for(int istep=0;istep<=nstep;istep++){
		
		//Update the total time
		ttime=ttime+dtime;
		
		/* Calculate the Laplacian of concentration with
		   five-point finite difference stencil */
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					
					ip=i+1;
					im=i-1;
					
					jp=j+1;
					jm=j-1;
					
					kp=k+1;
					km=k-1;
					
					if(ip==Nx){
						ip=0;
					}
					if(im==-1){
						im=(Nx-1);
					}
					
					if(jp==Ny){
						jp=0;
					}
					if(jm==-1){
						jm=(Ny-1);
					}
					
					if(kp==Nz){
						kp=0;
					}
					if(km==-1){
						km=(Nz-1);
					}
					
					// x
					hne=con[(ip*Ny+j)*Nz+k];
					hnw=con[(im*Ny+j)*Nz+k];
					// y
					hns=con[(i*Ny+jm)*Nz+k];
					hnn=con[(i*Ny+jp)*Nz+k];
					// z
					hnu=con[(i*Ny+j)*Nz+kp];
					hnd=con[(i*Ny+j)*Nz+km];
					// current
					hnc=con[ijk];
					
					lap_con[ijk] = (hne + hnw -2.0*hnc)/(dx*dx)
								  +(hns + hnn -2.0*hnc)/(dy*dy)
								  +(hnu + hnd -2.0*hnc)/(dz*dz);
					
					//Calculate the derivative of free energy
					dfdcon=free_energy_ch_3d((double)con[ijk]);
					
					/* Evaluate the terms in Eq.4.17, and
					   accumulate them in array dummy[Nx][Ny][Nz] for
					   each grid point in the simulation cell */
					dummy[ijk]=dfdcon-grad_coef*lap_con[ijk];
					//
				}//end for(k
			}//end for(j
		}//end for(i
		
		
		//Time integration of Eq.4.16
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					
					ip=i+1;
					im=i-1;
					
					jp=j+1;
					jm=j-1;
					
					kp=k+1;
					km=k-1;
					
					if(ip==Nx){
						ip=0;
					}
					if(im==-1){
						im=(Nx-1);
					}
					
					if(jp==Ny){
						jp=0;
					}
					if(jm==-1){
						jm=(Ny-1);
					}
					
					if(kp==Nz){
						kp=0;
					}
					if(km==-1){
						km=(Nz-1);
					}
					
					// x
					hne=dummy[(ip*Ny+j)*Nz+k];
					hnw=dummy[(im*Ny+j)*Nz+k];
					// y
					hns=dummy[(i*Ny+jm)*Nz+k];
					hnn=dummy[(i*Ny+jp)*Nz+k];
					// z
					hnu=dummy[(i*Ny+j)*Nz+kp];
					hnd=dummy[(i*Ny+j)*Nz+km];
					// current
					hnc=dummy[ijk];
					
					/* Calculate the Laplacing of the terms inside
					   the parenthesis in Eq.4.16 */
					lap_dummy[ijk] = (hne + hnw -2.0*hnc)/(dx*dx)
									+(hns + hnn -2.0*hnc)/(dy*dy)
									+(hnu + hnd -2.0*hnc)/(dz*dz);
					
					/* Explicit Euler time integration of
					   concentration field, Eq.4.16 */
					con[ijk] = con[ijk] + dtime*mobility*lap_dummy[ijk];
					
					/* If there are small variations, 
					   set the max and min values to file */
					if(con[ijk]>=0.9999){
					   con[ijk]=0.9999;
					}
					if(con[ijk]<0.0000){
					   con[ijk]=0.0000;
					}
				//
				}//end for(k
			}//end for(j
		}//end for(i
		
		//print results
		// If print frequency reached, print the results to file
		if(fmod(istep,nprint)==0){
			/* Calculate the total bulk energy and print the result to
			   time_energy.out file for 2D plots */
			energy = calculate_energy_3d(Nx,Ny,Nz,con,grad_coef);
			
			//print the average free energy density value to file
			fprintf(out2, "%14.6e %14.6e \n",ttime, energy);
			
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_3D(Nx,Ny,Nz,dx,dy,dz,istep,con);
			
			printf("done step: %5d \n",istep);
		}//end if
	}//end for(istep
	
	//calculate the execution time and print it
	// Calculate the compute time and print it to screen
	end = clock();
	compute_time = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Compute Time: %lf \n", compute_time);
}
