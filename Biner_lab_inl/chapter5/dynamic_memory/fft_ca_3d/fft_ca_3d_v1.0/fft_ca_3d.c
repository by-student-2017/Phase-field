/* 3D semi-implicit spectral phase-field code
   for solving Allen-Cahn equation */

/* This program solves conserved phase-field equation with
   Fourier spactral method by taking into account the effects of
   elastic inhomogeneities and lattice defects, based on
   solution of stress-strain fields ith Green's tensor and
   Fourier transformations. the time integration is
   carried out by using semi-implicit time machining scheme. */

#include <stdio.h>
#include <stdlib.h> //rand() and malloc
#include <math.h> //mod() and -lm
#include <time.h>

#include <fftw3.h>
//gcc test.c -lfftw3
//#include <mpi.h> //mpi version
//#include <fftw3-mpi.h> //mpi version

void init_grain_micro_3d();
double free_energy_fd_ca_3d();
void prepare_fft_3d();
void write_vtk_grid_values_3D();

int main(){
	clock_t start, end;
	double compute_time;
	
	//get initial wall time
	//(Get initial wall clock time beginning of the execution)
	start = clock();
	
	//open an output file for writing total bulk energy values
	FILE *out2=fopen("area_frac.out","w");
	
	//simulation cell parameters (These values are dummy)
	int Nx=64; //Number of grid points in the x-direction
	int Ny=64; //Number of grid points in the y-direction
	int Nz=2; //Number of grid points in the y-direction
	
	int ngrain=2;
	
	//The distance between two grid points in x,y-direction
	double dx=0.5; //Grid spacing between two grid pints in x-direction
	double dy=0.5; //Grid spacing between two grid pints in y-direction
	double dz=0.5; //Grid spacing between two grid pints in z-direction
	
	//time integration parameters
	int nstep=100000; //Number of time integration steps
	int nprint=100; //Output frequency to write the results to file
	double dtime=0.005; //Time increment for the numerical integration
	double ttime=0.0;   //Total time
	double coefA=1.0;
	
	//material specific parameters
	
	//Initial concentrations of alloying elements
	double mobil=5.0;  //The value of mobility coefficient
	double grcoef=0.1; //The value of gradient energy coefficients [J(nm)^2/mol]
	
	//----- ----- ----- -----
	int ijk; //ijk=(i*Ny+j)*Nz+k;
	//----- ----- ----- -----
	
	/* Generate initial grain microstructure
	   iflag=1 is for bi-crystal and
	   iflag=2 is for polycrystal */
	int iflag=2;
	if(iflag==2){
		//read polycrystal microstructure
		FILE *in=fopen("grain_25.inp","r");
		fscanf(in,"%5d %5d %5d %5d ",&Nx,&Ny,&Nz,&ngrain);
		fclose(in);
	}
	//----- ----- ----- -----
	int NxNyNz=Nx*Ny*Nz; //Total number of grid points in the simulation cell
	double *etas = (double *)malloc(sizeof(double)*( NxNyNz*ngrain ));
	int   *glist = (int *)malloc(sizeof(int)*( ngrain ));
	//----- ----- ----- -----
	if(iflag==2){
		FILE *in=fopen("grain_25.inp","r");
		fscanf(in,"%5d %5d %5d %5d ",&Nx,&Ny,&Nz,&ngrain);
		//
		int nline=1;
		int ri, rj, rk, rigrain;
		double reta;
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					fscanf(in,"%5d %5d %5d %5d %lf",&ri,&rj,&rk,&rigrain,&reta);
					//----- ----- ----- -----
					if( i != ri ){ printf("Don't match x data, Line %5d \n",nline); exit(1); }
					if( j != rj ){ printf("Don't match y data, LIne %5d \n",nline); exit(1); }
					if( k != rk ){ printf("Don't match z data, LIne %5d \n",nline); exit(1); }
					nline = nline + 1;
					//----- ----- ----- -----
					ijk=(i*Ny+j)*Nz+k;
					etas[ijk*ngrain+rigrain]=reta;
					//----- ----- ----- -----
				}
			}
		}
		fclose(in);
		//initialize glist
		for(int igrain=0;igrain<ngrain;igrain++){
			glist[igrain]=1.0;
		}
	}
	if(iflag==1){ init_grain_micro_3d(Nx,Ny,Nz,dx,dy,dz,iflag,ngrain,etas,glist); }
	
	//----- ----- ----- ----- ----- -----
	const int fftsizex = Nx, fftsizey = Ny, fftsizez = Nz;
	const int scale=fftsizex*fftsizey*fftsizez;
	double fftw3d_scale = (double)scale;
	/* fftw_complex *in, *out; // in[i][0] for real, in[i][1] for imag.
	    in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	   out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	   fftw_plan plan, iplan;
	   plan = fftw_plan_dft_2d(fftsize, fftsize, in, out, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	  iplan = fftw_plan_dft_2d(fftsize, fftsize, in, out, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT */
	//
	//array
	fftw_complex *eta, *etak;
	eta  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	etak = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	fftw_plan plan_eta, iplan_etak;
	 plan_eta  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez,  eta, etak, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_etak = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez,  etak,eta,  FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//array
	fftw_complex *dfdeta, *dfdetak;
	 dfdeta = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	dfdetak = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	fftw_plan plan_dfdeta;
	 plan_dfdeta  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, dfdeta, dfdetak, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//----- ----- ----- ----- ----- -----
	
	//----- ----- ----- ----- ----- ----- ----- ----- -----
	//double kx[Nx];
	double *kx = (double *)malloc(sizeof(double)*( Nx ));
	//double ky[Ny];
	double *ky = (double *)malloc(sizeof(double)*( Ny ));
	//double kz[Ny];
	double *kz = (double *)malloc(sizeof(double)*( Nz ));
	//double k2[Nx][Ny][Nz];
	double *k2 = (double *)malloc(sizeof(double)*( NxNyNz ));
	//double k4[Nx][Ny][Nz];
	double *k4 = (double *)malloc(sizeof(double)*( NxNyNz ));
	//----- ----- ----- ----- ----- ----- ----- ----- -----
	
	//prepare fft (output: kx,ky,kz,k2,k4)
	prepare_fft_3d(Nx,Ny,Nz,dx,dy,dz,kx,ky,kz,k2,k4); //get FFT coefficients
	
	double *eta_in = (double *)malloc(sizeof(double)*( NxNyNz ));
	double   *eta2 = (double *)malloc(sizeof(double)*( NxNyNz ));
	
	//----- ----- ----- -----
	double numer, denom;
	//----- ----- ----- -----
	double grain_sum;
	//----- ----- ----- -----
	int ncount;
	//----- ----- ----- -----
	
	//evolve (Time evolution of microstructure)
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
						for(int k=0;k<Nz;k++){
							ijk=(i*Ny+j)*Nz+k;
							eta[ijk][0] = etas[ijk*ngrain+igrain];
							eta[ijk][1] = 0.0;
							//-----
							eta_in[ijk] = eta[ijk][0];
						}
					}
				}
				
				//derivative of free energy
				//Calculate the derivative of elastic energy
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						for(int k=0;k<Nz;k++){
							ijk=(i*Ny+j)*Nz+k;
							dfdeta[ijk][0] = free_energy_fd_ca_3d(i,j,k,Nx,Ny,Nz,ngrain,etas,eta_in,igrain);
							dfdeta[ijk][1] = 0.0;
						}
					}
				}
				
				/* Take the values of concentration, derivative of free energy and
				   derivative of elastic energy from real space to Fourier space (forward FFT) */
				//etak=fft3(eta);
				fftw_execute(plan_eta);
				//dfdetak=fft3(dfdeta);
				fftw_execute(plan_dfdeta);
				
				/* Semi-implicit time integration of concentration field at
				   Fourier space (Eq.5.14) */
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						for(int k=0;k<Nz;k++){
							ijk=(i*Ny+j)*Nz+k;
							//
							denom=1.0+dtime*coefA*mobil*grcoef*k2[ijk];
							//
							numer=dtime*mobil*dfdetak[ijk][0];
							etak[ijk][0]=(etak[ijk][0]-numer)/denom;
							//
							numer=dtime*mobil*dfdetak[ijk][1];
							etak[ijk][1]=(etak[ijk][1]-numer)/denom;
						}
					}
				}
				
				/* Take concentration field from Fourier space back to
				   real space (inverse FFT) */
				//eta=real(ifft3(etak));
				fftw_execute(iplan_etak);
				
				//for small deviations
				// For small deviations from max and min values, reset the limits
				grain_sum=0.0;
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						for(int k=0;k<Nz;k++){
							ijk=(i*Ny+j)*Nz+k;
							//----- ----- ----- -----
							 eta[ijk][0] =  eta[ijk][0]/fftw3d_scale;
							 eta[ijk][1] =  eta[ijk][1]/fftw3d_scale;
							//----- ----- ----- -----
							if(eta[ijk][0]>=0.9999){
								eta[ijk][0]=0.9999;
							}
							if(eta[ijk][0]<=0.0001){
								eta[ijk][0]=0.0001;
							}
							eta[ijk][1]=0.0;
							//----- ----- ----- -----
							/* Calculate the total area of the current grain,
							   also return the order parameter values from
							   the temporary array eta[Nx][Ny] to 
							   common array etas[Nx][Ny][ngrain] */
							grain_sum = grain_sum + eta[ijk][0];
							etas[ijk*ngrain+igrain] = eta[ijk][0];
							//----- ----- ----- -----
						}
					}
				}
				
				//Check volume fraction of current grain
				/* Check the area fraction of the current grain.
				   If it is less than 0.001, set the its value in 
				   glist[ngrain] as zero which indicates that it is extinct. 
				   Also print message "grain # is eliminated" to screen. */
				grain_sum = grain_sum/NxNyNz;
				
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
					for(int k=0;k<Nz;k++){
						ijk=(i*Ny+j)*Nz+k;
						eta2[ijk]=0.0;
					}
				}
			}
			fprintf(out2, "%14.6e ",ttime);
			
			for(int igrain=0;igrain<ngrain;igrain++){
				ncount=0;
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						for(int k=0;k<Nz;k++){
							ijk=(i*Ny+j)*Nz+k;
							//eta2[ijk]=eta2[ijk]+etas[ijk*ngrain+igrain]*etas[ijk*ngrain+igrain];
							eta2[ijk]=eta2[ijk]+etas[ijk*ngrain+igrain]*etas[ijk*ngrain+igrain]*igrain;
							if(etas[ijk*ngrain+igrain]>=0.5){
								ncount=ncount+1;
							}//end if
						}//end for(k
					}//end for(j
				}//end for(i
				ncount=ncount/NxNyNz;
				fprintf(out2, "%5d ",ncount);
			}//end igrain
			fprintf(out2, "\n");
			
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_3D(Nx,Ny,Nz,dx,dy,dz,istep,eta2);
			
			printf("done step: %5d \n",istep);
		}//end if
	}//end for(istep
	
	//calculate the execution time and print it
	/* Calculate the compute time and print it to screen */
	end = clock();
	compute_time = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Compute Time: %lf \n", compute_time);
	
	//----- ----- ----- ----- ----- -----
	fftw_destroy_plan(plan_eta);
	fftw_destroy_plan(iplan_etak);
	fftw_destroy_plan(plan_dfdeta);
	//----- ----- ----- ----- ----- -----
	fftw_free(eta);
	fftw_free(etak);
	fftw_free(dfdeta);
	fftw_free(dfdetak);
	//
	free(kx);
	free(ky);
	free(k2);
	free(k4);
	//
	free(etas);
	free(glist);
	free(eta2);
	free(eta_in);
	//----- ----- ----- ----- ----- -----
}
