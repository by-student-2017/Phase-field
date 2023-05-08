/* Finite-difference phase-field code for
   solving Allen-Cahn equation */

/* This program solves the non-concerved multicomponent
   Allen-Cahn equation with finite difference algorithm using
   five-point stencil. The time integration is carried out with
   explicit Euler scheme. */

#include <stdio.h> //printf()
#include <stdlib.h> //rand()
#include <math.h> //mod() and -lm
#include <time.h>

#define Nx 100 //Number of grid points in the x-direction
#define Ny 100 //Number of grid points in the y-direction
#define npart 2 //Number of particles in the simulation, set to either 2 or 9

	double       eta[Nx][Ny];
	double   lap_eta[Nx][Ny]; //Laplacian
	double       con[Nx][Ny];
	double   lap_con[Nx][Ny]; //Laplacian
	double     dummy[Nx][Ny];
	double lap_dummy[Nx][Ny];
	double      etas[Nx][Ny][npart];
	double      phi2[Nx][Ny];

void micro_sint_pre_2d();
double free_energy_sint_2d_con();
double free_energy_sint_2d_eta();
void write_vtk_grid_values_2D();

int main(){
	// Get initial wall clock time beginning of the execution
	clock_t start, end;
	double compute_time;
	//get initial wall time
	start = clock();
	
	//simulation cell parameters
	//int Nx=64; //Number of grid points in the x-direction
	//int Ny=64; //Number of grid points in the y-direction
	//int NxNy=Nx*Ny; //Total number of grid points in the simulation cell
	
	//The distance between two grid points in x,y-direction
	double dx=0.5; //Grid spacing between two grid pints in x-direction
	double dy=0.5; //Grid spacing between two grid pints in y-direction
	
	//time integration parameters
	int nstep=12500; //Number of time integration steps
	int nprint=100; //Output frequency to write the results to file
	double dtime=1.0e-4; //Time increment for the numerical integration
	double ttime=0.0;   //Total time
	
	//int npart=2; //Number of particles in the simulation, set to either 2 or 9
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* Model-specific parameters (see Eqs.4.40,4.43,and Table 4.1) 
	   non-dimensionalized parameters */
	double coefm=5.0; //Gradient coefficient for the concentration field, Eq.4.38, kappa_rho
	double coefk=2.0; //Gradient coefficient for the order parameters, Eq.4.38, kappa_eta
	double coefl=10.0; //Mobility of the order parameters, Eq.4.43, (L is the grain boundary mobility)
	//
	double dvol=0.040; //Bulk diffusion coefficient, Dvol. Dvol is the bulk diffusivity in the lattice.
	double dvap=0.002; //Vapor diffusion coefficient, Dvap. Dvap is the diffusivity of the vapor phase
	double dsur=16.0;  //surface diffusion coefficient, Dsurf. Dsurf is the surface diffusivity.
	double dgrb=1.6;   //Grain goundary diffusion coefficient, DGB. DGB is the grain boundary diffusibity.
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* Note
	   where kappa_rho and kappa_eta are the gradient energy coefficients for
	   concentration and grain boundary energies, respectively. */
	
	// Generate initial grain microstructure
	int iflag=1;
	micro_sint_pre_2d(Nx,Ny,npart,iflag,etas,con); // iflag=1 only
	
	/* initialize temporary array eta[Nx][Ny] for
	   the time integration of the order parameters */
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			eta[i][j]=0.0;
		}
	}
	
	//----- ----- ----- -----
	int ip,im;
	int jp,jm;
	//
	double hne,hnw,hns,hnn;
	double hnc;
	//----- ----- ----- -----
	double dfdcon;
	double dfdeta;
	//----- ----- ----- -----
	double phi;
	double sum;
	double mobil;
	//----- ----- ----- -----
	
	//Time integration
	// Time evolution of concentration and order parameter
	for(int istep=0;istep<=nstep;istep++){
		
		//Update the total time
		ttime=ttime+dtime;
		
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//Time evolution of concentration field, rho (Cahn-Hilliard equation)
		iflag=1; //Set iflag=1 for the functional derivative of free energy
		//
		/* Calculate the Laplacian of the term inside the parenthesis in Eq.4.40 */
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				
				ip=i+1;
				im=i-1;
				
				jp=j+1;
				jm=j-1;
				
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
				
				hne=con[ip][j];
				hnw=con[im][j];
				hns=con[i][jm];
				hnn=con[i][jp];
				hnc=con[i][j];
				
				// Calculate the Laplacian of concentration
				lap_con[i][j] = (hne + hnw -2.0*hnc)/(dx*dx)
							   +(hns + hnn -2.0*hnc)/(dy*dy);
				
				//Calculate the derivative of free energy
				dfdcon = free_energy_sint_2d_con(i,j,Nx,Ny,con,eta,etas,npart);
				
				// Calculate the terms inside parenthesis in Eq.4.40
				//  coefm = kappa_rho, lap_con = Laplacian of the rho
				dummy[i][j] = dfdcon - ( 0.5 * coefm * lap_con[i][j] );
			}//end for(j
		}//end for(i
		
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				
				ip=i+1;
				im=i-1;
				
				jp=j+1;
				jm=j-1;
				
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
				
				hne=dummy[ip][j];
				hnw=dummy[im][j];
				hns=dummy[i][jm];
				hnn=dummy[i][jp];
				hnc=dummy[i][j];
				
				// Calculate the laplacian of the terms inside parenthesis in Eq.4.40
				lap_dummy[i][j] = (hne + hnw -2.0*hnc)/(dx*dx)
							 	 +(hns + hnn -2.0*hnc)/(dy*dy);
				
				//Mobility
				/* Calculate the diffusivity/mobility parameter for
				   the current grid point */
				
				//Calculate the value of interpolation function, Eq.4.42
				// phi = rho^3 * (10 - 15*rho + 6*rho^2) (Eq.4.42)
				phi = con[i][j]*con[i][j]*con[i][j]*(
						10.0 - 15.0*con[i][j] + 6.0*con[i][j]*con[i][j]
						);
				
				//Calculate the summations in Eq.4.41
				sum=0.0;
				for(int ipart=0;ipart<npart;ipart++){
					for(int jpart=0;jpart<npart;jpart++){
						if(ipart != jpart){
							//Not eta^2, but eta_i*eta_j
							sum = sum + etas[i][j][ipart]*etas[i][j][jpart];
						}
					}
				}
				
				/* The value of the diffusivity/mobility parameter at
				   the grid point */
				mobil = dvol*phi
					  + dvap*(1.0-phi)
					  + dsur*con[i][j]*(1.0-con[i][j])
					  + dgrb*sum;
				
				/* Explicit Euler time integration of
				   concentration field, Eq.4.40 for
				   the current grid point*/
				con[i][j] = con[i][j] + ( dtime * mobil * lap_dummy[i][j] );
				
				/* If there are small variations, 
				   set the max and min values to the limits */
				if(con[i][j]>=0.9999){
				   con[i][j]=0.9999;
				}
				if(con[i][j]<0.0001){
				   con[i][j]=0.0001;
				}//Not 0, but 0.00001
			}//end for(j
		}//end for(i
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//evolve etas
		// Evolute the non-conserved order parameters, eta (Allen-Cahn equation)
		iflag=2; //set iflag=2 for functional derivative of free energy
		//
		//Loop over the number of particles
		for(int ipart=0;ipart<npart;ipart++){
			//
			/* Assign the current particle order parameter values from
			   common array etas[Nx][Ny][npart] to eta[Nx][Ny] */
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					eta[i][j]=etas[i][j][ipart];
				}
			}
			//
			/* Calculate the Laplacian of concentration with
			   five-point finite difference stencil */
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					
					ip=i+1;
					im=i-1;
					
					jp=j+1;
					jm=j-1;
					
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
					
					hne=eta[ip][j];
					hnw=eta[im][j];
					hns=eta[i][jm];
					hnn=eta[i][jp];
					hnc=eta[i][j];
					
					// The Laplacian of the order parameter, Eq.4.43
					lap_eta[i][j] = (hne + hnw -2.0*hnc)/(dx*dx)
								   +(hns + hnn -2.0*hnc)/(dy*dy);
					
					//Calculate the derivative of free energy
					/* Functional derivative of the free energy for
					   the current order parameter at the current grid point */
					dfdeta = free_energy_sint_2d_eta(i,j,Nx,Ny,con,eta,etas,npart);
					
					/* Explicit Euler time integration of the order parameter,
					   Eq.4.43 for the current grid point */
					//  coefk=kappa_eta, lap_eta=Laplacian of the eta
					eta[i][j] = eta[i][j] - (dtime * coefl)*(dfdeta - (0.5 * coefk * lap_eta[i][j]));
					
					/* If there are small variations, 
					   set the max and min values to the limits */
					if(eta[i][j]>=0.9999){
					   eta[i][j]=0.9999;
					}
					if(eta[i][j]<0.0001){
					   eta[i][j]=0.0001;
					}//Not 0, but 0.00001
				}//end for(j
			}//end for(i
			
			/* Return the current values of the order parameter, eta[Nx][Ny] to
			   common array etas[Nx][Ny][npart] */
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					etas[i][j][ipart]=eta[i][j];
				}
			}
			//
		}//end for(ipart
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
		//print results
		// If print frequency reached, print the results to file
		if(fmod(istep,nprint)==0){
			
			
			//initialization
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					phi2[i][j]=0.0;
				}
			}
			
			for(int ipart=0;ipart<npart;ipart++){
				//
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						phi2[i][j] = phi2[i][j] + (etas[i][j][ipart] * etas[i][j][ipart]);
						//phi2[i][j] = phi2[i][j] + (etas[i][j][ipart] * etas[i][j][ipart])*(ipart+1.0);
					}//end for(j
				}//end for(i
				//
			}//end igrain
			
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_2D(Nx,Ny,dx,dy,istep,con,phi2);
			
			//printf("done step: %5d, %14.6e \n",istep,ttime);
			printf("done step: %5d \n",istep);
		}//end if
		//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	}//end for(istep
	
	//calculate the execution time and print it
	// Calculate the compute time and print it to screen
	end = clock();
	compute_time = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Compute Time: %lf \n", compute_time);
}
