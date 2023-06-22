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

void init_grain_micro_2d();
double free_energy_fd_ca_2d();
void write_vtk_grid_values_2D();

int main(){
	// Get initial wall clock time beginning of the execution
	clock_t start, end;
	double compute_time;
	//get initial wall time
	start = clock();
	
	//open an output file for writing total bulk energy values
	FILE *out2=fopen("area_frac.out","w");
	
	//simulation cell parameters (These values are dummy)
	int Nx=64; //Number of grid points in the x-direction
	int Ny=64; //Number of grid points in the y-direction
	
	int ngrain=2;
	
	//The distance between two grid points in x,y-direction
	double dx=0.5; //Grid spacing between two grid pints in x-direction
	double dy=0.5; //Grid spacing between two grid pints in y-direction
	
	//time integration parameters
	int nstep=100000; //Number of time integration steps
	int nprint=100; //Output frequency to write the results to file
	double dtime=0.005; //Time increment for the numerical integration
	double ttime=0.0;   //Total time
	
	//material specific parameters
	
	//Initial concentrations of alloying elements
	double mobil=5.0;  //The value of mobility coefficient
	double grcoef=0.1; //The value of gradient energy coefficients [J(nm)^2/mol]
	
	//----- ----- ----- -----
	int ij; //ij=i*Ny+j;
	//----- ----- ----- -----
	
	/* Generate initial grain microstructure
	   iflag=1 is for bi-crystal and
	   iflag=2 is for polycrystal */
	int iflag=2;
	if(iflag==2){
		//read polycrystal microstructure
		FILE *in=fopen("grain_25.inp","r");
		fscanf(in,"%5d %5d %5d ",&Nx,&Ny,&ngrain);
		fclose(in);
	}
	//----- ----- ----- -----
	int NxNy=Nx*Ny; //Total number of grid points in the simulation cell
	double     *eta = (double *)malloc(sizeof(double)*( NxNy ));
	double *lap_eta = (double *)malloc(sizeof(double)*( NxNy ));
	double    *etas = (double *)malloc(sizeof(double)*( NxNy*ngrain ));
	int *glist = (int *)malloc(sizeof(int)*( ngrain ));
	double    *eta2 = (double *)malloc(sizeof(double)*( NxNy ));
	//----- ----- ----- -----
	if(iflag==2){
		FILE *in=fopen("grain_25.inp","r");
		fscanf(in,"%5d %5d %5d ",&Nx,&Ny,&ngrain);
		//
		int nline=1;
		int ri, rj, rigrain;
		double reta;
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				fscanf(in,"%5d %5d %5d %lf",&ri,&rj,&rigrain,&reta);
				//----- ----- ----- -----
				if( i != ri ){ printf("Don't match x data, Line %5d \n",nline); exit(1); }
				if( j != rj ){ printf("Don't match y data, LIne %5d \n",nline); exit(1); }
				nline = nline + 1;
				//----- ----- ----- -----
				ij=i*Ny+j;
				etas[ij*ngrain+rigrain]=reta;
				//----- ----- ----- -----
			}
		}
		fclose(in);
		//initialize glist
		for(int igrain=0;igrain<ngrain;igrain++){
			glist[igrain]=1.0;
		}
	}
	if(iflag==1){ init_grain_micro_2d(Nx,Ny,dx,dy,iflag,ngrain,etas,glist); }
	
	//----- ----- ----- -----
	int ip,im;
	int jp,jm;
	//
	double hne,hnw,hns,hnn;
	double hnc;
	//----- ----- ----- -----
	double dfdeta;
	//----- ----- ----- -----
	double grain_sum;
	//----- ----- ----- -----
	int ncount;
	//----- ----- ----- -----
	
	//Time integration
	for(int istep=0;istep<=nstep;istep++){
		
		//Update the total time
		ttime=ttime+dtime;
		
		//Loop over each grain
		for(int igrain=0;igrain<ngrain;igrain++){
			
			/* If glist is equal to one, which indicates that
			   the current grain area fraction is greater than 0.001,
			   continue the calculation. Otherwise, the current grain
			   does not exist anymore */
			if(glist[igrain]==1){
				
				/* Assign order parameters to temporary array eta[Nx][Ny] from
				   the common array etas[Nx][Ny][ngrain] for the current grain */
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						ij=i*Ny+j;
						eta[ij]=etas[ij*ngrain+igrain];
					}
				}
				
				// Calculate the Laplacian term in Eq.4.35
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						ij=i*Ny+j;
						
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
						
						lap_eta[ij] = (hne + hnw -2.0*hnc)/(dx*dx)
									 +(hns + hnn -2.0*hnc)/(dy*dy);
						
						/* Determine the derivative of the free energy for
						   the current grain at the current grain point under
						   consideration, first three terms inside the bracket in Eq.4.35 */
						dfdeta=free_energy_fd_ca_2d(i,j,Nx,Ny,ngrain,etas,eta,igrain);
						
						/* With explicit Euler integration, time integration of order parameter of
						   current point at current grid point, Eq.4.35 */
						eta[ij] = eta[ij] - dtime*mobil*(dfdeta-grcoef*lap_eta[ij]);
						
						/* If there are small variations, 
						   set the max and min values to the limits */
						if(eta[ij]>=0.9999){
						   eta[ij]=0.9999;
						}
						if(eta[ij]<1.0e-6){
						   eta[ij]=1.0e-6;
						}
						//
					}//end for(j
				}//end for(i
				
				/* Calculate the total area of the current grain,
				   also return the order parameter values from
				   the temporary array eta[Nx][Ny] to 
				   common array etas[Nx][Ny][ngrain] */
				grain_sum=0.0;
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						ij=i*Ny+j;
						etas[ij*ngrain+igrain]=eta[ij];
						grain_sum=grain_sum+eta[ij];
					}//end for(j
				}//end for(i
				
				//Check volume fraction of current grain
				/* Check the area fraction of the current grain.
				   If it is less than 0.001, set the its value in 
				   glist[ngrain] as zero which indicates that it is extinct. 
				   Also print message "grain # is eliminated" to screen. */
				grain_sum = grain_sum/NxNy;
				
				if(grain_sum<=0.001){
					glist[igrain]=0;
					printf("grain: No. %3d is eliminated \n",igrain);
				}
				
			}//end if(glist
		}//end igrain
		
		//print results
		// If print frequency reached, print the results to file
		if(fmod(istep,nprint)==0){
			//write vtk file & calculate are function of grains
			/* Prepare the data to be written to vtk file and
			   calculate the area fraction of each grain and
			   print them to file area_fract.out. */
			
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					ij=i*Ny+j;
					eta2[ij]=0.0;
				}
			}
			fprintf(out2, "%14.6e ",ttime);
			
			for(int igrain=0;igrain<ngrain;igrain++){
				ncount=0;
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						ij=i*Ny+j;
						//eta2[ij]=eta2[ij]+etas[ij*ngrain+igrain]*etas[ij*ngrain+igrain];
						eta2[ij]=eta2[ij]+etas[ij*ngrain+igrain]*etas[ij*ngrain+igrain]*(igrain+1.0);
						if(etas[ij*ngrain+igrain]>=0.5){
							ncount=ncount+1;
						}//end if
					}//end for(j
				}//end for(i
				ncount=ncount/NxNy;
				fprintf(out2, "%5d ",ncount);
			}//end igrain
			fprintf(out2, "\n");
			
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_2D(Nx,Ny,dx,dy,istep,eta2);
			
			printf("done step: %5d \n",istep);
		}//end if
	}//end for(istep
	
	//calculate the execution time and print it
	// Calculate the compute time and print it to screen
	end = clock();
	compute_time = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Compute Time: %lf \n", compute_time);
}
