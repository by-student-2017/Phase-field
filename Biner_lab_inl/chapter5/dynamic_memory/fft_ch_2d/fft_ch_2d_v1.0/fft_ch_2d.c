/* 2D semi-implicit spectral phase-field code
   for solving Cahn-Hilliard equation */

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

void micro_ch_pre_2d();
void prepare_fft_2d();
double free_energy_ch_2d();
double calculate_energy_2d();
void write_vtk_grid_values_2D();

int main(){
	clock_t start, end;
	double compute_time;
	
	//get initial wall time
	//(Get initial wall clock time beginning of the execution)
	start = clock();
	
	//simulation cell parameters
	int Nx=64;
	int Ny=64;
	//int NxNy=Nx*Ny; //Total number of grid points in the simulation cell
	
	//open an output file for writing total bulk energy values
	FILE *out2=fopen("time_energy.out","w");
	
	//The distance between two grid points in x,y-direction
	double dx=1.0; // [nm] unit ?
	double dy=1.0; // [nm] unit ?
	
	//time integration parameters
	int nstep=2000; //Number of time steps
	int nprint=50;  //Print frequency to write the results to file
	double dtime=1.0e-2; //Time increment for numerical integration
	double ttime=0.0;    //Total time
	double coefA=1.0;
	
	//material specific parameters
	
	//Initial concentrations of alloying elements
	double c0=0.40;       //Initial concentraion (20%Cr-containing Fe-Cr alloy
	double mobility=1.0;  //The value of mobility coefficient (dimensionless)
	double grad_coef=0.5; //The value of gradient energy coefficients [J(nm)^2/mol]
	
	double energy=0.0;
	
	int ij; //ij=(i*Ny+j);
	
	//----- ----- ----- ----- ----- -----
	const int fftsizex = Nx, fftsizey = Ny;
	/* fftw_complex *in, *out; // in[i][0] for real, in[i][1] for imag.
	    in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	   out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	   fftw_plan plan, iplan;
	   plan = fftw_plan_dft_2d(fftsize, fftsize, in, out, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	  iplan = fftw_plan_dft_2d(fftsize, fftsize, in, out, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT */
	//
	//array
	fftw_complex *con, *conk;
	con  = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	conk = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	fftw_plan plan_con, iplan_conk;
	 plan_con  = fftw_plan_dft_2d(fftsizex, fftsizey, con, conk, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_conk = fftw_plan_dft_2d(fftsizex, fftsizey, conk,con,  FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//array
	fftw_complex *dfdcon, *dfdconk;
	 dfdcon = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	dfdconk = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	fftw_plan plan_dfdcon;
	 plan_dfdcon  = fftw_plan_dft_2d(fftsizex, fftsizey, dfdcon, dfdconk, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//----- ----- ----- ----- ----- -----
	
	//double con_out[Nx][Ny];
	double *con_out = (double *)malloc(sizeof(double)*( Nx*Ny ));
	
	//prepare microstructure
	// Initialize the concentration filed with random modulation
	micro_ch_pre_2d(Nx,Ny,c0,con_out);
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			ij=i*Ny+j;
			con[ij][0] = con_out[ij];
			con[ij][1] = 0.0;
		}
	}
	
	//double kx[Nx];
	double *kx = (double *)malloc(sizeof(double)*( Nx ));
	//double ky[Ny];
	double *ky = (double *)malloc(sizeof(double)*( Ny ));
	//double k2[Nx][Ny];
	double *k2 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double k4[Nx][Ny];
	double *k4 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	
	//prepare fft (output: kx,ky,kz,k2,k4)
	prepare_fft_2d(Nx,Ny,dx,dy,kx,ky,k2,k4); //Calculate coefficients of Fourier transformation
	
	double numer, denom;
	
	//evolve (Time evolution of microstructure)
	for(int istep=0;istep<=nstep;istep++){
		
		//Update the total time
		ttime=ttime+dtime;
		
		//derivative of free energy
		//Calculate the derivative of elastic energy
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ij=i*Ny+j;
				dfdcon[ij][0] = free_energy_ch_2d((double)con_out[ij]);
				dfdcon[ij][1] = 0.0;
			}
		}
		
		/* Take the values of concentration, derivative of free energy and
		   derivative of elastic energy from real space to Fourier space (forward FFT) */
		//conk=fft2(con);
		fftw_execute(plan_con);
		//dfdconk=fft2(dfdcon);
		fftw_execute(plan_dfdcon);
		
		/* Semi-implicit time integration of concentration field at
		   Fourier space (Eq.5.14) */
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ij=i*Ny+j;
				//
				denom=1.0+dtime*coefA*mobility*grad_coef*k4[ij];
				//
				numer=dtime*mobility*k2[ij]*dfdconk[ij][0];
				conk[ij][0]=(conk[ij][0]-numer)/denom;
				//
				numer=dtime*mobility*k2[ij]*dfdconk[ij][1];
				conk[ij][1]=(conk[ij][1]-numer)/denom;
			}
		}
		
		/* Take concentration field from Fourier space back to
		   real space (inverse FFT) */
		//con=real(ifft2(conk));
		fftw_execute(iplan_conk);
		
		//for small deviations
		// For small deviations from max and min values, reset the limits
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ij=i*Ny+j;
				//----- ----- ----- -----
				 con[ij][0] =  con[ij][0]/(fftsizex*fftsizey);
				 con[ij][1] =  con[ij][1]/(fftsizex*fftsizey);
				//----- ----- ----- -----
				if(con[ij][0]>=0.9999){
					con[ij][0]=0.9999;
				}
				if(con[ij][0]<=0.0001){
					con[ij][0]=0.0001;
				}
				con[ij][1]=0.0;
				//----- ----- ----- -----
				con_out[ij]  = con[ij][0];
				//----- ----- ----- -----
			}
		}
		
		//print results
		/* If print frequency is reached, output the results to file */
		if(fmod(istep,nprint)==0){
			/* Calculate the total bulk energy and print the result to
			   time_energy.out file for 2D plots */
			energy = calculate_energy_2d(Nx,Ny,con_out,grad_coef);
			
			//print the average free energy density value to file
			fprintf(out2, "%14.6e %14.6e \n",ttime, energy);
			
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_2D(Nx,Ny,dx,dy,istep,con_out);
			
			printf("done step: %5d, time: %f \n",istep,ttime*(mobility/(dx*dx)));
			/* The quantities having the dimension of distance were normalized with the
			   magnitude of the Burger's vector, the  quantities having the dimension of
			   energy were normalized with RT, and the time t was normalized with M/(dx^2).
			   The initial concentration was modulated by setting the noise term to
			   0.001 in function micro_ch_pre_2d.c */
		}
		
	}//end of time step (evolve,for)
	
	//calculate the execution time and print it
	/* Calculate the compute time and print it to screen */
	end = clock();
	compute_time = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Compute Time: %lf \n", compute_time);
	
	//----- ----- ----- ----- ----- -----
	fftw_destroy_plan(plan_con);
	fftw_destroy_plan(iplan_conk);
	fftw_destroy_plan(plan_dfdcon);
	//----- ----- ----- ----- ----- -----
	fftw_free(con);
	fftw_free(conk);
	fftw_free(dfdcon);
	fftw_free(dfdconk);
	//
	free(kx);
	free(ky);
	free(k2);
	free(k4);
	//
	free(con_out);
	//----- ----- ----- ----- ----- -----
}
