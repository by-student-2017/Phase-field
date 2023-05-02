/* 2D semi-implicit spectral phase-field code 
  for solving precipitation in FeCr alloy */

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

void dislo_strain_2d();
double FeCr_chem_poten_2d();
void green_tensor_2d();
void prepare_fft_2d();
void solve_elasticity_2d();
void micro_ch_pre_2d();
void write_vtk_grid_values_2D();

int main(){
	clock_t start, end;
	double compute_time;
	
	//get initial wall time
	//(Get initial wall clock time beginning of the execution)
	start = clock();
	
	//simulation cell parameters
	int Nx=128;
	int Ny=128;
	
	//Total number of grid points in the simulation cell
	//int NxNy=Nx*Ny;
	
	//The distance between two grid points in x,y-direction
	double dx=1.0; // [nm] unit ?
	double dy=1.0; // [nm] unit ?
	
	//time integration parameters
	int nstep=10000; //Number of time steps
	int nprint=50;  //Print frequency to write the results to file
	double dtime=1.0e-2; //Time increment for numerical integration
	double ttime=0.0;    //Total time
	double coefA=2.0;
	
	//material specific parameters
	
	//Initial concentrations of alloying elements
	double c0=0.20;       //Initial concentraion
	double mobility=1.0;  //The value of mobility coefficient (dimensionless)
	double grad_coef=0.5; //The value of gradient energy coefficients [J(nm)^2/mol]
	
	double tempr=535.0;
	double RT=8.314462*tempr;
	
	//elastic constants
	//Elastic constants of Fe-rich phase
	double cm11=233.10e3;
	double cm12=135.44e3;
	double cm44=178.30e3;
	//
	//Elastic constants of Cr-rich phase
	double cp11=350.00e3;
	double cp12= 67.80e3;
	double cp44=100.80e3;
	
	//eigen strains
	//The value of eigenstrains for Cr-rich phase
	double ei0=0.006; //Maginitude of eigenstrains
	
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
	fftw_complex *cr, *crk;
	 cr = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	crk = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	fftw_plan plan_cr, iplan_crk;
	 plan_cr  = fftw_plan_dft_2d(fftsizex, fftsizey, cr, crk, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_crk = fftw_plan_dft_2d(fftsizex, fftsizey, crk,cr,  FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//array
	fftw_complex *dfdcr, *dfdcrk;
	 dfdcr = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	dfdcrk = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	//fftw_plan plan_dfdcr, iplan_dfdcrk;
	fftw_plan plan_dfdcr;
	 plan_dfdcr  = fftw_plan_dft_2d(fftsizex, fftsizey, dfdcr, dfdcrk, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//iplan_dfdcrk = fftw_plan_dft_2d(fftsizex, fftsizey, dfdcrk, dfdcr, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//array
	fftw_complex *delsdc, *delsdck;
	 delsdc = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	delsdck = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	//fftw_plan plan_delsdc, iplan_delsdck;
	fftw_plan plan_delsdc;
	 plan_delsdc  = fftw_plan_dft_2d(fftsizex, fftsizey, delsdc, delsdck, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//iplan_delsdck = fftw_plan_dft_2d(fftsizex, fftsizey, delsdck, delsdc, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
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
	
	int ii;
	
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
	
	/* idislo=1 for dislocation diploe,
	   idislo=2 for dislocation array */
	int idislo=1;
	dislo_strain_2d(Nx,Ny,idislo,ed11,ed22,ed12);
	
	//Applied strains
	// The components of applied strains
	double ea[3]; //Magnitude of applied strains
	ea[0]=0.0;
	ea[1]=0.0;
	ea[2]=0.0;
	
	//int iflag=1;
	
	//prepare microstructure
	micro_ch_pre_2d(Nx,Ny,c0,cr); //Initialize microstructure
	
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
	
	//double cr_out[Nx][Ny];
	double *cr_out = (double *)malloc(sizeof(double)*( Nx*Ny ));
	
	//evolve (Time evolution of microstructure)
	for(int istep=0;istep<=nstep;istep++){
		
		//Update the total time
		ttime=ttime+dtime;
		
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
			cr, delsdc);
		// Note: tmatx is real part only
		
		//derivative of chemical energy
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				//Normalize the derivative elastic energy with RT
				delsdc[ii][0]=delsdc[ii][0]/RT;
				delsdc[ii][1]=0.0;
				//Calculate derivative of chemical energy
				dfdcr[ii][0]=FeCr_chem_poten_2d((double)cr[ii][0],tempr);
				dfdcr[ii][1]=0.0;
			}
		}
		
		/* Take the values of concentration, derivative of free energy and
		   derivative of elastic energy from real space to Fourier space (forward FFT) */
		//crk=fft2(cr);
		fftw_execute(plan_cr);
		//dfdcrk=fft2(dfdcr);
		fftw_execute(plan_dfdcr);
		//delsdck=fft2(delsdc);
		fftw_execute(plan_delsdc);
		
		/* Semi-implicit time integration of Cr concentration field at
		   Fourier space (Eq.5.50) */
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				//
				denom=1.0+dtime*coefA*mobility*grad_coef*k4[ii];
				//
				numer=dtime*mobility*k2[ii]*(dfdcrk[ii][0]+delsdck[ii][0]);
				crk[ii][0]=(crk[ii][0]-numer)/denom;
				//
				numer=dtime*mobility*k2[ii]*(dfdcrk[ii][1]+delsdck[ii][1]);
				crk[ii][1]=(crk[ii][1]-numer)/denom;
			}
		}
		
		/* Take concentration field from Fourier space back to
		   real space (inverse FFT) */
		//cr=real(ifft2(crk));
		fftw_execute(iplan_crk);
		
		//for small deviations
		// For small deviations from max and min values, reset the limits
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				//----- ----- ----- -----
				 cr[ii][0] =  cr[ii][0]/(fftsizex*fftsizey);
				 cr[ii][1] =  cr[ii][1]/(fftsizex*fftsizey);
				//----- ----- ----- -----
				if(cr[ii][0]>=0.9999){
					cr[ii][0]=0.9999;
				}
				if(cr[ii][0]<=0.0001){
					cr[ii][0]=0.0001;
				}
				cr[ii][1]=0.0;
				//----- ----- ----- -----
				cr_out[ii]  = cr[ii][0];
				//----- ----- ----- -----
			}
		}
		
		//print results
		/* If print frequency is reached, output the results to file */
		if(fmod(istep,nprint)==0){
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_2D(Nx,Ny,dx,dy,istep,cr_out);
			
			printf("done step: %5d \n",istep);
		}
		
	}//end of time step (evolve,for)
	
	//calculate the execution time and print it
	/* Calculate the compute time and print it to screen */
	end = clock();
	compute_time = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Compute Time: %lf \n", compute_time);
	
	//----- ----- ----- ----- ----- -----
	fftw_destroy_plan(plan_cr);
	fftw_destroy_plan(iplan_crk);
	fftw_destroy_plan(plan_dfdcr);
	//fftw_destroy_plan(iplan_dfdcrk);
	fftw_destroy_plan(plan_delsdc);
	//fftw_destroy_plan(iplan_delsdck);
	//----- ----- ----- ----- ----- -----
	fftw_free(cr);
	fftw_free(crk);
	fftw_free(dfdcr);
	fftw_free(dfdcrk);
	fftw_free(delsdc);
	fftw_free(delsdck);
	//
	fftw_free(s11);
	fftw_free(s22);
	fftw_free(s12);
	//
	//fftw_free(e11); //why ? free(): invalid pointer
	//fftw_free(e22); //why ? free(): invalid pointer
	//fftw_free(e12); //why ? free(): invalid pointer
	//----- ----- ----- ----- ----- -----
	//free(ed11); //why ? free(): invalid pointer
	//free(ed22); //why ? free(): invalid pointer
	//free(ed12); //why ? free(): invalid pointer
	//
	free(kx);
	free(ky);
	free(k2);
	free(k4);
	//
	free(tmatx);
	//
	free(cr_out);
	//----- ----- ----- ----- ----- -----
}
