/* Finite-difference phase-field code for
   solving Allen-Cahn equation */

/* This program solves the non-concerved multicomponent
   Allen-Cahn equation with finite difference algorithm using
   five-point stencil. The time integration is carried out with
   explicit Euler scheme. */

#include <stdio.h> //printf()
#include <stdlib.h> //rand() and malloc
#include <math.h> //mod() and -lm
#include <time.h>

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
	int Nx=100; //Number of grid points in the x-direction
	int Ny=100; //Number of grid points in the y-direction
	int NxNy=Nx*Ny; //Total number of grid points in the simulation cell
	
	//The distance between two grid points in x,y-direction
	double dx=0.5; //Grid spacing between two grid pints in x-direction
	double dy=0.5; //Grid spacing between two grid pints in y-direction
	
	//time integration parameters
	int nstep=12500; //Number of time integration steps
	int nprint=100; //Output frequency to write the results to file
	double dtime=1.0e-4; //Time increment for the numerical integration
	double ttime=0.0;   //Total time
	
	//Number of particles in the simulation, set to either 2 or 10 (bug fix: why not 9)
	int npart=10; //"int npart=9" doesn't work for some reason.
	
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
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ---------- ----- ----- ----- ----- ----- ----- 
	//double       con[Nx][Ny]; //concentraion
	double       *con = (double *)malloc(sizeof(double)*( NxNy ));
	//double   lap_con[Nx][Ny]; //Laplacian
	double   *lap_con = (double *)malloc(sizeof(double)*( NxNy ));
	//----- ----- ----- ---------- ----- ----- ----- ----- ----- ----- 
	//double     dummy[Nx][Ny];
	double     *dummy = (double *)malloc(sizeof(double)*( NxNy ));
	//double lap_dummy[Nx][Ny]; //Laplacian
	double *lap_dummy = (double *)malloc(sizeof(double)*( NxNy ));
	//----- ----- ----- ---------- ----- ----- ----- ----- ----- ----- 
	//double       eta[Nx][Ny]; //order parameter
	double       *eta = (double *)malloc(sizeof(double)*( NxNy ));
	//double   lap_eta[Nx][Ny]; //Laplacian
	double   *lap_eta = (double *)malloc(sizeof(double)*( NxNy ));
	//----- ----- ----- ---------- ----- ----- ----- ----- ----- ----- 
	//double      etas[Nx][Ny][npart];
	double      *etas = (double *)malloc(sizeof(double)*( NxNy*npart ));
	//----- ----- ----- ---------- ----- ----- ----- ----- ----- ----- 
	//double      phi2[Nx][Ny];
	double      *phi2 = (double *)malloc(sizeof(double)*( NxNy ));
	//----- ----- ----- ---------- ----- ----- ----- ----- ----- ----- 
	
	int ij; //ij=(i*Ny+j);
	
	/* initialize temporary array eta[Nx][Ny] for
	   the time integration of the order parameters */
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			ij=(i*Ny+j);
			//----- ----- -----
			con[ij]=0.0;
			dummy[ij]=0.0;
			eta[ij]=0.0;
			//----- ----- -----
			lap_con[ij]=0.0;
			lap_dummy[ij]=0.0;
			lap_eta[ij]=0.0;
			//----- ----- -----
			for(int ipart=0;ipart<npart;ipart++){
				etas[ij*npart+ipart]=0.0;
			}
			//----- ----- -----
		}
	}
	
	// Generate initial grain microstructure
	int iflag=1;
	micro_sint_pre_2d(Nx,Ny,npart,iflag,etas,con); // iflag=1 only
	
	//----- ----- ----- -----
	int ip,im;
	int jp,jm;
	//
	double hne,hnw;
	double hns,hnn;
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
				ij=(i*Ny+j);
				
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
				
				hne=con[ip*Ny+j];
				hnw=con[im*Ny+j];
				hns=con[i*Ny+jm];
				hnn=con[i*Ny+jp];
				hnc=con[ij];
				
				// Calculate the Laplacian of concentration
				lap_con[ij] = (hne + hnw -2.0*hnc)/(dx*dx)
							 +(hns + hnn -2.0*hnc)/(dy*dy);
				
				//Calculate the derivative of free energy
				dfdcon = free_energy_sint_2d_con(i,j,Nx,Ny,con,eta,etas,npart);
				
				// Calculate the terms inside parenthesis in Eq.4.40
				//  coefm = kappa_rho, lap_con = Laplacian of the rho
				dummy[ij] = dfdcon - ( 0.5 * coefm * lap_con[ij] );
			}//end for(j
		}//end for(i
		
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ij=(i*Ny+j);
				
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
				
				hne=dummy[ip*Ny+j];
				hnw=dummy[im*Ny+j];
				hns=dummy[i*Ny+jm];
				hnn=dummy[i*Ny+jp];
				hnc=dummy[ij];
				
				// Calculate the laplacian of the terms inside parenthesis in Eq.4.40
				lap_dummy[ij] = (hne + hnw -2.0*hnc)/(dx*dx)
							   +(hns + hnn -2.0*hnc)/(dy*dy);
				
				//Mobility
				/* Calculate the diffusivity/mobility parameter for
				   the current grid point */
				
				//Calculate the value of interpolation function, Eq.4.42
				// phi = rho^3 * (10 - 15*rho + 6*rho^2) (Eq.4.42)
				phi = con[ij]*con[ij]*con[ij]*(
						10.0 - 15.0*con[ij] + 6.0*con[ij]*con[ij]
						);
				
				//Calculate the summations in Eq.4.41
				sum=0.0;
				for(int ipart=0;ipart<npart;ipart++){
					for(int jpart=0;jpart<npart;jpart++){
						if(ipart != jpart){
							//Not eta^2, but eta_i*eta_j
							sum = sum + etas[ij*npart+ipart]*etas[ij*npart+jpart];
						}
					}
				}
				
				/* The value of the diffusivity/mobility parameter at
				   the grid point */
				mobil = dvol*phi
					  + dvap*(1.0-phi)
					  + dsur*con[ij]*(1.0-con[ij])
					  + dgrb*sum;
				
				/* Explicit Euler time integration of
				   concentration field, Eq.4.40 for
				   the current grid point*/
				con[ij] = con[ij] + ( dtime * mobil * lap_dummy[ij] );
				
				/* If there are small variations, 
				   set the max and min values to the limits */
				if(con[ij]>=0.9999){
				   con[ij]=0.9999;
				}
				if(con[ij]<0.0001){
				   con[ij]=0.0001;
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
					ij=(i*Ny+j);
					eta[ij]=etas[ij*npart+ipart];
				}
			}
			//
			/* Calculate the Laplacian of concentration with
			   five-point finite difference stencil */
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					ij=(i*Ny+j);
					
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
					
					hne=eta[ip*Ny+j];
					hnw=eta[im*Ny+j];
					hns=eta[i*Ny+jm];
					hnn=eta[i*Ny+jp];
					hnc=eta[ij];
					
					// The Laplacian of the order parameter, Eq.4.43
					lap_eta[ij] = (hne + hnw -2.0*hnc)/(dx*dx)
								 +(hns + hnn -2.0*hnc)/(dy*dy);
					
					//Calculate the derivative of free energy
					/* Functional derivative of the free energy for
					   the current order parameter at the current grid point */
					dfdeta = free_energy_sint_2d_eta(i,j,Nx,Ny,con,eta,etas,npart);
					
					/* Explicit Euler time integration of the order parameter,
					   Eq.4.43 for the current grid point */
					//  coefk=kappa_eta, lap_eta=Laplacian of the eta
					eta[ij] = eta[ij] - (dtime * coefl)*(dfdeta - (0.5 * coefk * lap_eta[ij]));
					
					/* If there are small variations, 
					   set the max and min values to the limits */
					if(eta[ij]>=0.9999){
					   eta[ij]=0.9999;
					}
					if(eta[ij]<0.0001){
					   eta[ij]=0.0001;
					}//Not 0, but 0.00001
				}//end for(j
			}//end for(i
			
			/* Return the current values of the order parameter, eta[Nx][Ny] to
			   common array etas[Nx][Ny][npart] */
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					ij=(i*Ny+j);
					etas[ij*npart+ipart]=eta[ij];
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
					ij=(i*Ny+j);
					phi2[ij]=0.0;
				}
			}
			
			for(int ipart=0;ipart<npart;ipart++){
				//
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						ij=(i*Ny+j);
						phi2[ij] = phi2[ij] + (etas[ij*npart+ipart] * etas[ij*npart+ipart]);
						//phi2[ij] = phi2[ij] + (etas[ij*npart+ipart] * etas[ij*npart+ipart])*(ipart+1.0);
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
	
	//----- ----- ----- 
	free(con);
	free(lap_con);
	free(dummy);
	free(lap_dummy);
	free(eta);
	free(lap_eta);
	free(etas);
	free(phi2);
	//----- ----- ----- 
}
