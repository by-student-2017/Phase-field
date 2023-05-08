/* Finite-difference phase-field code for cell-dynamics */

/* This program solves the phase-field modeling of 
   the cell dynamics with finite difference algorithm. */

#include <stdio.h> //printf()
#include <stdlib.h> //rand()
#include <math.h> //mod() and -lm
#include <time.h>

void micro_poly_cell_2d();
double free_energy_2d();
void write_vtk_grid_values_2D();

int main(){
	
	// Get initial wall clock time beginning of the execution
	clock_t start, end;
	double compute_time;
	//get initial wall time
	start = clock();
	
	//simulation cell parameters
	int Nx=200; //Number of grid points in the x-direction
	int Ny=200; //Number of grid points in the y-direction
	//int NxNy=Nx*Ny; //Total number of grid points in the simulation cell
	
	//The distance between two grid points in x,y-direction
	double dx=1.0; //Grid spacing between two grid pints in x-direction
	double dy=1.0; //Grid spacing between two grid pints in y-direction
	
	//time integration parameters
	int nstep=3000; //Number of time integration steps
	int nprint=50;  //Output frequency to write the results to file
	double dtime=5.e-3; //Time increment for the numerical integration
	
	//material specific parameters
	/* For poly-cell simulations set the desired number of cells and
	   radius of the cells (in grid points) */
	int ncell=80;
	double R=12.0;
	
	/* For two-cell simulations set the ncell=2 and 
	   radius of the cells (in grid points). */
	//int ncell=2;
	//double R=25.0;
	
	// Set the values of the parameters (Eq.4.53 and 4.59, and Table 4.3)
	//The non-dimensionalized values of the parameters used in the model
	//double gamma=5.0; // For ncell 2
	double gamma=2.5; // For ncell 80
	double lambda=7.0;
	double kappa=60.0;
	double mu=40.0;
	double kisa=1.5e3;
	
	//The value of pi
	/* Pre-calculate the some of the expressions as constants
	   (see Eqs. 4.55, 4.58 and 4.59). */
	double pix=4.0*atan(1.0);
	double const1=30.0/(lambda*lambda);
	double const2=2.0*mu/(pix*R*R);
	double const3=60.0*kappa/(lambda*lambda*kisa);
	
	//set array
	//double     phi[Nx][Ny];
	double *phi = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double    phis[Nx][Ny][ncell];
	double *phis = (double *)malloc(sizeof(double)*( Nx*Ny*ncell ));
	//double lap_phi[Nx][Ny];
	double *lap_phi = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double   phidx[Nx][Ny];
	double *phidx = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double   phidy[Nx][Ny];
	double *phidy = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double    phi1[Nx][Ny]; //write vtk file
	double *phi1 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double  term4x[Nx][Ny];
	//double             vac[ncell];
	double *vac = (double *)malloc(sizeof(double)*( ncell ));
	//int              ccell[ncell];
	int *ccell = (int *)malloc(sizeof(int)*( ncell ));
	
	//int nccell=0;
	//int used_ncell=0;
	int nccell[2]; //0:nccell, 1:used_ncell;
	
	//prepare initial microstructure
	micro_poly_cell_2d(Nx,Ny,R,
					ncell,phis,vac,
					nccell,ccell);
	printf("check nccell[0]=nccell: %5d \n",nccell[0]);
	printf("check nccell[1]=used_ncell: %5d \n",nccell[1]);
	
	//----- ----- ----- -----
	double gamma_cell;
	//----- ----- ----- -----
	int ip,im;
	int jp,jm;
	//
	double hne,hnw,hns,hnn;
	double hnc;
	//----- ----- ----- -----
	double vinteg;
	double vintegx;
	double vintegy;
	//----- ----- ----- -----
	double sum_phi;
	//----- ----- ----- -----
	double dfdphi;
	//----- ----- ----- -----
	double term2;
	double term3;
	//----- ----- ----- -----
	double vac_cell;
	double vnx;
	double vny;
	double vnphi; //term4x[i][j]=vnphi;
	//----- ----- ----- -----
	double add;
	//----- ----- ----- -----
	int ij; //ij=(i*Ny+j);
	int ijci; //ijci=(ij*ncell)+icell;
	int ijcj; //ijcj=(ij*ncell)+jcell;
	//----- ----- ----- -----
	
	//Evolution (Time evolution)
	for(int istep=0;istep<=nstep;istep++){
		
		//After 500 time steps, increase the time increment to 0.01
		if(istep>=500){
			dtime=1.0e-2;
		}
		
		//Loop over each cell in the simulation
		for(int icell=0;icell<nccell[1];icell++){
			
			//cell elasticity parameter
			/* Assign the cell elasticity parameter depending on
			   its character, determined by ccell. */
			gamma_cell = gamma;
			for(int iccell=0; iccell<nccell[0];iccell++){ //Note: "ccell", not "cell"
				if(icell == ccell[iccell]){ //Note: "[iccell]", not "[icell]"
					gamma_cell=0.5*gamma;
				}
			}
			
			//Assign cells
			/* Assign current cell order parameter values from
			   common array phis[Nx][Ny][ncell] to phi[Nx][Ny]. */
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					ij=(i*Ny+j);
					ijci=(ij*ncell)+icell;
					phi[ij]=phis[ijci];
				}
			}
			
			//calculate the laplacian and gradient terms
			/* Calculate the Laplacian and Cartesian derivative of
			   the current cell with five-point stencil. */
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
					
					hne=phi[ip*Ny+j];
					hnw=phi[im*Ny+j];
					hns=phi[i*Ny+jm];
					hnn=phi[i*Ny+jp];
					hnc=phi[i*Ny+j];
					
					lap_phi[ij] = (hne + hnw -2.0*hnc)/(dx*dx)
								 +(hns + hnn -2.0*hnc)/(dy*dy);
					
					//gradients of phi
					phidx[ij]=(phi[ip*Ny+j]-phi[im*Ny+j])/(2.0*dx);
					phidy[ij]=(phi[i*Ny+jp]-phi[i*Ny+jm])/(2.0*dy);
				}//end for(j
			}//end for(i
			
			//volume integrations
			// Calculate the volume integrals that appears in Eqs.4.58 and 4.60
			vinteg =0.0;
			vintegx=0.0;
			vintegy=0.0;
			
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					ij=(i*Ny+j);
					
					// Calculate the integral value in Eq.4.58
					vinteg = vinteg + phi[ij]*phi[ij];
					
					sum_phi=0.0;
					
					for(int jcell=0;jcell<nccell[1];jcell++){
						if(icell != jcell){
							ijcj=(ij*ncell)+jcell;
							sum_phi = sum_phi + phis[ijcj]*phis[ijcj];
						}
					}
					
					// Calculate the volume integral for x-component of the velocity vector in Eq.4.60
					vintegx = vintegx + phi[ij] * phidx[ij] * sum_phi;
					
					// Calculate the volume integral for y-component of the velocity vector in Eq.4.60
					vintegy = vintegy + phi[ij] * phidy[ij] * sum_phi;
				}
			}
			
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					ij=(i*Ny+j);
					
					//Second term
					// For the current grid point, calculate the second term of Eq.4.59
					//derivative of free energy
					dfdphi = free_energy_2d(i,j,Nx,Ny,icell,ncell,nccell[1],gamma,kappa,phi,phis);
					term2 = -const1 * dfdphi;
					
					//Third term
					// For the current grid point, calculate the third term of Eq.4.59
					term3 = -const2 * (vinteg - pix*R*R)*phi[ij];
					
					//cell velocity
					// If time steps <= 200 do not include the self-propulsion
					if(istep<=200){
						vac_cell = 0.0;
					}else{
						/* Otherwise take the randomly assigned value of 
						   self-propulsion value in function micro_poly_cell. */
						vac_cell = vac[icell];
					}
					
					//cell velocity vector
					// Determine the values of the velocity vector of current cell.
					vnx = vac_cell + const3 * vintegx; //Cell velocity value in the x-direction
					vny = vac_cell + const3 * vintegy; //Cell velocity value in the y-direction
					
					// If ncell=2 take the self-propulsion only in the x-direction
					if(ncell==2){
						vnx = vac_cell;
						vny = 0.0;
					}
					
					//Fourth term
					// Calculate the fourth term in Eq.4.59
					vnphi = vnx * phidx[ij] + vny * phidy[ij];
					//term4x[i][j] = vnphi;
					
					//time integration
					// Evolve phi with explicit Euler time integration scheme (Eq.4.59)
					phi[ij] = phi[ij] + dtime * ( gamma_cell * lap_phi[ij]
												+ term2 + term3 - vnphi );
				}//end for(j
			}//end for(i
			
			/* Return the current values of cell under
			   consideration to common array phis[Nx][Ny][ncell]. */
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					ij=(i*Ny+j);
					//for small deviations
					// For small deviations, set the order parameters values to 0.9999 and 0
					if(phi[ij]>=0.9999){
					   phi[ij]=0.9999;
					}
					if(phi[ij]<0.0000){
					   phi[ij]=0.0000;
					}
					
					ijci=(ij*ncell)+icell;
					phis[ijci]=phi[ij];
				}//end for(j
			}//end for(i
		}//end for(icell
		
		//print results
		// If print frequency reached, print the results to file
		if(fmod(istep,nprint)==0){
			//write vtk file
			/* Write the results for contour plot in vtk format
			   to be viewed by Paraview. */
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					ij=(i*Ny+j);
					phi1[ij]=0.0;
				}
			}
			
			/* Adds a small value to the order parameters of the soft cells,
			   in order to distinguish them from the other cells during
			   the contour plotting in Paraview. */
			for(int icell=0;icell<nccell[1];icell++){
				add=1.0;
				for(int iccell=0;iccell<nccell[0];iccell++){ //Note: "ccell", not "cell"
					if(icell == ccell[iccell]){ //Note: "[iccell]", not "[icell]"
						add=1.2;
					}
				}
				
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						ij=(i*Ny+j);
						ijci=(ij*ncell)+icell;
						phi1[ij] = phi1[ij] + ( phis[ijci]*phis[ijci]*add );
					}
				}
				write_vtk_grid_values_2D(Nx,Ny,dx,dy,istep,phi1);
			}
			printf("done step: %5d \n",istep);
		}//end if
	}//end for(istep
	
	//calculate the execution time and print it
	// Calculate the compute time and print it to screen
	end = clock();
	compute_time = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Compute Time: %lf \n", compute_time);
	
	free(phi);
	free(phis);
	free(lap_phi);
	free(phidx);
	free(phidy);
	free(phi1);
	free(vac);
	free(ccell);
}