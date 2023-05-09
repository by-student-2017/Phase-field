/* Finite-difference phase-field code for
   dendritic solidification */

/* This program solves the non-concerved multicomponent
   Allen-Cahn equation with finite difference algorithm using
   five-point stencil. The time integration is carried out with
   explicit Euler scheme. */

#include <stdio.h> //printf()
#include <stdlib.h> //rand() and malloc
#include <math.h> //mod() and -lm
#include <time.h>

void nucleus_2d();
void write_vtk_grid_values_2D();

int main(){
	// Get initial wall clock time beginning of the execution
	clock_t start, end;
	double compute_time;
	//get initial wall time
	start = clock();
	
	//simulation cell parameters
	int Nx=300; //Number of grid points in the x-direction
	int Ny=300; //Number of grid points in the y-direction
	int NxNy=Nx*Ny; //Total number of grid points in the simulation cell
	
	//The distance between two grid points in x,y-direction
	double dx=0.03; //Grid spacing between two grid pints in x-direction
	double dy=0.03; //Grid spacing between two grid pints in y-direction
	
	//time integration parameters
	int nstep=4000; //Number of time integration steps
	int nprint=50; //Output frequency to write the results to file
	double dtime=1.0e-4; //Time increment for the numerical integration
	double ttime=0.0;   //Total time
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	//Material specific parameters (see Eqs.4.51 and 4.52, also Table 4.2)
	double tau=0.0003; //tau=1/M_phi
	double epsilonb=0.01; //The mean value of epsilon
	//double mu=1.0;     //Kinetic growth coefficient [m/(K*s)]
	double kappa=1.8;  //Dimensionless latent heat
	double delta=0.02; //strength of the anisotropy
	double aniso=6.0;  //4:cubit lattices, 6:hexagonal lattices
	double alpha=0.9;  //a positive constant
	double gamma=10.0; //a constant
	double teq=1.0;   //Equilibrium temperature, Teq
	double theta0=0.2; //Initial offset angle and taken as a constant
	double seed=5.0;   //The size of the initial seed (In grid numbers, see function nucleus)
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	/* tempr is nondimensionalized so that the characteristic cooling tmperature is zero and the equilibirum temperature is one.
	   kappa is proportional to the latent heat and inversely proportional to the strength of cooling. */
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	/* The interfacial mobility, M_phi = sqrt(2*Wsl)/(6*epsilon0)*(Tm/L)*mu = (b*Tm*mu)/(3*(3*dx)*L)
	   The energy barrier, Wsl = (6*gamma0*b)/(3*dx), where dx is [nm] units,
	   b = 2*tanh^-1(1-2*lambda), lambda=0.1 (usually),
	   The mean value of epsilon, epsilonb= sqrt((3*(3*dx)*gamma0)/b),
	   Tm = melting point [K], gamma0 is interfacial energy density [J/m^2], L is latent heat of solidification [J/m^3]
	   The thermal diffusivity, a = K/Cp = (Thermal conductivity [W/(m*K)])/(specific heat per unit volume [J/(K*m^3)]) */
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	double pi=4.0*atan(1.0); // The value of pi
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	//double           phi[Nx][Ny]; //non-conserved phase-field parameter
	double           *phi = (double *)malloc(sizeof(double)*( NxNy ));
	//double       lap_phi[Nx][Ny]; //Laplacian of phi
	double       *lap_phi = (double *)malloc(sizeof(double)*( NxNy ));
	//double         tempr[Nx][Ny]; //Temperature field which is derived from the conservation law of enthalpy
	double         *tempr = (double *)malloc(sizeof(double)*( NxNy ));
	//double     lap_tempr[Nx][Ny]; //Laplacian of tempr
	double     *lap_tempr = (double *)malloc(sizeof(double)*( NxNy ));
	//double         phidx[Nx][Ny]; //d(phi)/dx which is related to theta
	double         *phidx = (double *)malloc(sizeof(double)*( NxNy ));
	//double         phidy[Nx][Ny]; //d(phi)/dy which is related to theta
	double         *phidy = (double *)malloc(sizeof(double)*( NxNy ));
	//double       epsilon[Nx][Ny]; //the anisotrophic gradient energy coefficient which determines the thickness of the interface layer
	double       *epsilon = (double *)malloc(sizeof(double)*( NxNy ));
	//double epsilon_deriv[Nx][Ny]; //d(epsilon)/d(theta)
	double *epsilon_deriv = (double *)malloc(sizeof(double)*( NxNy ));
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	
	int ij;  //ij=(i*Ny+j);
	int ipj; //ipj=ip*Ny+j;
	int imj; //imj=im*Ny+j;
	int ijp; //ijp=i*Ny+jp;
	int ijm; //ijm=i*Ny+jm;
	
	//Initialize and introduce initial nuclei
	/* Initialize the phi and tempr arrays and
	   introduce the seed in the center of the simulation cell */
	nucleus_2d(Nx,Ny,seed,phi,tempr);
	
	//----- ----- ----- -----
	int ip,im;
	int jp,jm;
	//
	double hne,hnw,hns,hnn;
	double hnc;
	//----- ----- ----- -----
	double theta;
	double phiold;
	//----- ----- ----- -----
	double term1;
	double term2;
	double m; //driving force for the interface motion proportional to supercooling
	//----- ----- ----- -----
	
	//Time integration
	// Time evolution of Eqs.4.51 and 4.52
	for(int istep=0;istep<=nstep;istep++){
		
		//Update the total time
		ttime=ttime+dtime;
		
		// Calculate the Laplacian and epsilon
		/* Calculate the Laplacian of phi and the Laplacian of temperature and
		   the epsilon and its derivatives depsilon/dtheta in Eqs.4.51 and 4.52 */
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
				
				ipj=ip*Ny+j;
				imj=im*Ny+j;
				ijm=i*Ny+jm;
				ijp=i*Ny+jp;
				
				hne=phi[ipj];
				hnw=phi[imj];
				hns=phi[ijm];
				hnn=phi[ijp];
				hnc=phi[ij];
				
				// Calculate the Laplacian of phi in Eq.4.51
				lap_phi[ij] = (hne + hnw -2.0*hnc)/(dx*dx)
							 +(hns + hnn -2.0*hnc)/(dy*dy);
				
				hne=tempr[ipj];
				hnw=tempr[imj];
				hns=tempr[ijm];
				hnn=tempr[ijp];
				hnc=tempr[ij];
				
				// Calculate the Laplacian of temperature in Eq.4.52
				lap_tempr[ij] = (hne + hnw -2.0*hnc)/(dx*dx)
							   +(hns + hnn -2.0*hnc)/(dy*dy);
				
				//gradients of phi
				// Calculate the Cartesian derivatives of phi
				phidx[ij]=(phi[ipj]-phi[imj])/(2.0*dx);
				phidy[ij]=(phi[ijp]-phi[ijm])/(2.0*dy);
				
				//calculate angle, Eq.4.48
				// theta = arctan( (d(phi)/dy)/(d(phi)/dx) )
				//  theta=atan2(y,x) [radian] (range:[-pi,pi])
				theta = atan2(phidy[ij],phidx[ij]);
				
				//epsilon and its derivative
				// Calculate epsilon, Eqs4.46 and 4.47
				epsilon[ij] = epsilonb*( 1.0 + delta*cos(aniso*(theta-theta0)) );
				
				// Calculate derivative of epsilon, d(epsilon)/d(theta)
				epsilon_deriv[ij] = -epsilonb * aniso * delta * sin(aniso*(theta-theta0));
				//
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
				
				// Assign the values of phi from previous time step to phiold
				phiold=phi[ij];
				
				ijp=i*Ny+jp;
				ijm=i*Ny+jm;
				
				//first term, d(...)/dy
				// Calculate the first term of Eq.4.51
				term1 = (epsilon[ijp]*epsilon_deriv[ijp]*phidx[ijp]
					    -epsilon[ijm]*epsilon_deriv[ijm]*phidx[ijm]
						)/(2.0*dy);
				
				ipj=ip*Ny+j;
				imj=im*Ny+j;
				
				//second term, d(...)/dy
				// Calculate the second term of Eq.4.51
				term2 =-(epsilon[ipj]*epsilon_deriv[ipj]*phidy[ipj]
					    -epsilon[imj]*epsilon_deriv[imj]*phidy[imj]
						)/(2.0*dx);
				
				//facter m
				// Calculate value m of Eq.4.49
				//  m(T) = (alpha/pi)*arctan( gamma*(Teq - T) )
				//   theta = atan(x) [radian] (range:[-pi/2,pi/2])
				//   tan(atan(x))=x, atan(tan(x))=x
				m = (alpha/pi)*atan( gamma*(teq - tempr[ij]) );
				
				//time integration
				// Evolve phi with Euler time integration, Eq.4.50 and 4.51
				// The time-dependent Ginzburg-Landau or Allen-Cahn equation
				phi[ij] = phi[ij] + (dtime/tau)*( term1 + term2
						+epsilon[ij]*epsilon[ij]*lap_phi[ij]
						+phiold*(1.0 - phiold)*(phiold - 0.5 + m)
						);
				
				//evolve temperature
				// Evolve temperature with Euler time integration, Eq.4.52
				tempr[ij] = tempr[ij] + dtime * lap_tempr[ij]
							+ kappa*(phi[ij] - phiold);
				
				/* If there are small variations, 
				   set the max and min values to the limits */
				if(phi[ij]<1e-6){
				   phi[ij]=1e-6;
				}
				if(tempr[ij]<1e-6){
				   tempr[ij]=1e-6;
				}
				//
			}//end for(j
		}//end for(i
		
		//print results
		// If print frequency reached, print the results to file
		if(fmod(istep,nprint)==0){
			
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_2D(Nx,Ny,dx,dy,istep,phi,tempr);
			
			printf("done step: %5d, %14.6e \n",istep,ttime);
			//printf("done step: %5d \n",istep);
		}//end if
		
	}//end for(istep
	
	//calculate the execution time and print it
	// Calculate the compute time and print it to screen
	end = clock();
	compute_time = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Compute Time: %lf \n", compute_time);
}
