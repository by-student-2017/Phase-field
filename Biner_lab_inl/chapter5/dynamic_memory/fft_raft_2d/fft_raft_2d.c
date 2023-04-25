/* 2D semi-implicit spectral phase-field code 
  for solving precipitation under stress */

/* This program solves Cahn-Hilliard phase-field equation with
   semi-implicit Fourier spectral method by 
   taking into account the effects of elastic inhomogeneities and 
   applied stresses based on solution of stress-strain fields with
   Green's tensor and Fourier transformations.
     The time integration is carried out by using semi-implicit
   time marching scheme. */

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
void green_tensor_2d();
double free_energy_ch_2d();
void solve_elasticity_2d();
void write_vtk_grid_values_2D();

int main(){
	clock_t start, end;
	double compute_time;
	
	//get initial wall time
	//(Get initial wall clock time beginning of the execution)
	start = clock();
	
	//simulation cell parameters
	int Nx=256;
	int Ny=256;
	
	//Total number of grid points in the simulation cell
	//int NxNy=Nx*Ny;
	
	//The distance between two grid points in x,y-direction
	double dx=1.0; // [nm] unit ?
	double dy=1.0; // [nm] unit ?
	
	//time integration parameters
	int nstep=5000; //Number of time steps
	int nprint=25;  //Print frequency to write the results to file
	double dtime=5.0e-2; //Time increment for numerical integration
	double ttime=0.0;    //Total time
	double coefA=1.0;
	
	//material specific parameters
	
	//Initial concentrations of alloying elements
	double c0=0.40;       //Initial concentraion
	double mobility=1.0;  //The value of mobility coefficient (dimensionless)
	double grad_coef=0.5; //The value of gradient energy coefficients [J(nm)^2/mol]
	
	//elastic constants
	//Elastic constants of matrix phase
	double cm11=1400.0;
	double cm12= 600.0;
	double cm44= 400.0;
	//
	//Elastic constants of second phase
	double cp11=2.0*cm11;
	double cp12=2.0*cm12;
	double cp44=2.0*cm44;
	
	//eigen strains
	double ei0=0.01; //Maginitude of eigenstrains
	
	//Applied strains
	double ea[3]; //Magnitude of applied strains
	ea[0]=0.00;
	ea[1]=0.01;
	ea[2]=0.00;
	
	int ii;
	
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
	 con = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	conk = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	fftw_plan plan_con, iplan_conk;
	 plan_con  = fftw_plan_dft_2d(fftsizex, fftsizey, con, conk, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_conk = fftw_plan_dft_2d(fftsizex, fftsizey, conk, con, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//array
	fftw_complex *dfdcon, *dfdconk;
	 dfdcon = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	dfdconk = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	fftw_plan plan_dfdcon, iplan_dfdconk;
	 plan_dfdcon  = fftw_plan_dft_2d(fftsizex, fftsizey, dfdcon, dfdconk, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_dfdconk = fftw_plan_dft_2d(fftsizex, fftsizey, dfdconk, dfdcon, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//array
	fftw_complex *delsdc, *delsdck;
	 delsdc = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	delsdck = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	fftw_plan plan_delsdc, iplan_delsdck;
	 plan_delsdc  = fftw_plan_dft_2d(fftsizex, fftsizey, delsdc, delsdck, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_delsdck = fftw_plan_dft_2d(fftsizex, fftsizey, delsdck, delsdc, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//----- ----- ----- ----- ----- -----
	
	//----- ----- ----- -----
	fftw_complex *s11, *s22, *s12;
	 s11 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	 s22 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	 s12 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	//----- ----- ----- -----
	
	fftw_complex *e11, *e22, *e12;
	 e11 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	 e22 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	 e12 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	//----- ----- ----- ----- ----- -----
	
	//double ed11[Nx][Ny];
	double *ed11 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double ed22[Nx][Ny];
	double *ed22 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double ed12[Nx][Ny];
	double *ed12 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	
	//Initialize stress & strain componentes
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			ii=i*Ny+j;
			//----- ----- ----- -----
			s11[ii][0] = 0.0;
			s22[ii][0] = 0.0;
			s12[ii][0] = 0.0;
			//
			s11[ii][1] = 0.0;
			s22[ii][1] = 0.0;
			s12[ii][1] = 0.0;
			//----- ----- ----- -----
			e11[ii][0] = 0.0;
			e22[ii][0] = 0.0;
			e12[ii][0] = 0.0;
			//
			e11[ii][1] = 0.0;
			e22[ii][1] = 0.0;
			e12[ii][1] = 0.0;
			//----- ----- ----- -----
			//----- ----- ----- -----
			//Strain components due to lattice defects
			ed11[ii] = 0.0;
			ed22[ii] = 0.0;
			ed12[ii] = 0.0;
			//----- ----- ----- -----
		}
	}
	
	//prepare microstructure
	micro_ch_pre_2d(Nx,Ny,c0,con); //Initialize microstructure
	
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
	
	//double tmatx[Nx][Ny][2][2][2][2];
	double *tmatx = (double *)malloc(sizeof(double)*( Nx*Ny*2*2*2*2 )); //real part only
	
	//Greens tensor
	green_tensor_2d(Nx,Ny,kx,ky,cm11,cm12,cm44,cp11,cp12,cp44,tmatx); //Calculate Green's tensor
	
	double numer, denom;
	
	//double con_out[Nx][Ny];
	double *con_out = (double *)malloc(sizeof(double)*( Nx*Ny ));
	
	//evolve (Time evolution of microstructure)
	for(int istep=0;istep<=nstep;istep++){
		
		//Update the total time
		ttime=ttime+dtime;
		
		//derivative of free energy
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				//Calculate derivative of free energy
				dfdcon[ii][0]=free_energy_ch_2d((double)con[ii][0]);
				dfdcon[ii][1]=0.0;
			}
		}
		
		//derivative of elastic energy
		//Calculate the derivative of elastic energy
		solve_elasticity_2d(Nx,Ny,
			tmatx,
			s11,s22,s12,
			e11,e22,e12,
			ed11,ed22,ed12,
			cm11,cm12,cm44,
			cp11,cp12,cp44,
			ea,
			ei0,
			con, delsdc);
		// Note: tmatx is real part only
		
		/* Take the values of concentration, derivative of free energy and
		   derivative of elastic energy from real space to Fourier space (forward FFT) */
		//conk=fft2(con);
		fftw_execute(plan_con);
		//dfdconk=fft2(dfdcon);
		fftw_execute(plan_dfdcon);
		//delsdck=fft2(delsdc);
		fftw_execute(plan_delsdc);
		
		/* Semi-implicit time integration of concentration field at
		   Fourier space (Eq.5.50) */
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				//
				denom=1.0+dtime*coefA*mobility*grad_coef*k4[ii];
				//
				numer=dtime*mobility*k2[ii]*(dfdconk[ii][0]+delsdck[ii][0]);
				conk[ii][0]=(conk[ii][0]-numer)/denom;
				//
				numer=dtime*mobility*k2[ii]*(dfdconk[ii][1]+delsdck[ii][1]);
				conk[ii][1]=(conk[ii][1]-numer)/denom;
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
				ii=i*Ny+j;
				//----- ----- ----- -----
				 con[ii][0] =  con[ii][0]/(fftsizex*fftsizey);
				 con[ii][1] =  con[ii][1]/(fftsizex*fftsizey);
				//----- ----- ----- -----
				if(con[ii][0]>=0.9999){
					con[ii][0]=0.9999;
				}
				if(con[ii][0]<=0.0001){
					con[ii][0]=0.0001;
				}
				con[ii][1]=0.0;
				//----- ----- ----- -----
				con_out[ii]  = con[ii][0];
				//----- ----- ----- -----
			}
		}
		
		//print results
		/* If print frequency is reached, output the results to file */
		if(fmod(istep,nprint)==0){
			printf("done step: %5d \n",istep);
			
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_2D(Nx,Ny,dx,dy,istep,con_out);
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
	fftw_destroy_plan(iplan_dfdconk);
	fftw_destroy_plan(plan_delsdc);
	fftw_destroy_plan(iplan_delsdck);
	//----- ----- ----- ----- ----- -----
	fftw_free(con);
	fftw_free(conk);
	fftw_free(dfdcon);
	fftw_free(dfdconk);
	fftw_free(delsdc);
	fftw_free(delsdck);
	//
	fftw_free(s11);
	fftw_free(s22);
	fftw_free(s12);
	//
	fftw_free(e11);
	fftw_free(e22);
	fftw_free(e12);
	//----- ----- ----- ----- ----- -----
	free(ed11);
	free(ed22);
	free(ed12);
	//
	free(kx);
	free(ky);
	free(k2);
	free(k4);
	//
	free(tmatx);
	//
	free(con_out);
	//----- ----- ----- ----- ----- -----
}
