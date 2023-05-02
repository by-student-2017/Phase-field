/* 3D semi-implicit spectral phase-field code 
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

void dislo_strain_3d();
double FeCr_chem_poten_3d();
void green_tensor1_3D();
//void green_tensor2_3D();
void prepare_fft_3d();
void solve_elasticity_3d();
void micro_ch_pre_3d();
void write_vtk_grid_values_3D();

int main(){
	clock_t start, end;
	double compute_time;
	
	//get initial wall time
	//(Get initial wall clock time beginning of the execution)
	start = clock();
	
	//simulation cell parameters
	int Nx=128;
	int Ny=128;
	int Nz=2;
	
	//Total number of grid points in the simulation cell
	int NxNyNz=Nx*Ny*Nz;
	
	//The distance between two grid points in x,y-direction
	double dx=1.0; // [nm] unit ?
	double dy=1.0; // [nm] unit ?
	double dz=1.0; // [nm] unit ?
	
	//time integration parameters
	int nstep=10000; //Number of time steps
	int nprint=50;  //Print frequency to write the results to file
	double dtime=1.0e-2; //Time increment for numerical integration
	double ttime=0.0;    //Total time
	double coefA=2.0;
	
	//material specific parameters
	
	//Initial concentrations of alloying elements
	double c0=0.20;       //Initial concentraion (20%Cr-containing Fe-Cr alloy
	double mobility=1.0;  //The value of mobility coefficient (dimensionless)
	double grad_coef=0.5; //The value of gradient energy coefficients [J(nm)^2/mol]
	
	double tempr=535.0; //Annealing temperature [K]
	double RT=8.314462*tempr; //Gas constant x temperature
	
	//elastic constants
	//Elastic constants of Fe-rich phase [GPa]
	double cm11=233.10e3;
	double cm12=135.44e3;
	double cm44=178.30e3;
	//
	//Elastic constants of Cr-rich phase [GPa]
	double cp11=350.00e3;
	double cp12= 67.80e3;
	double cp44=100.80e3;
	
	//elastic constant of other materials
	//Ref: https://www.jstage.jst.go.jp/article/jsms/69/9/69_657/_pdf
	
	//eigen strains
	//The value of eigenstrains for Cr-rich phase
	double ei0=0.006; //Maginitude of eigenstrains
	
	//Applied strains
	// The components of applied strains
	double ea[6]; //Magnitude of applied strains
	ea[0]=0.0; //e11
	ea[1]=0.0; //e22
	ea[2]=0.0; //e33
	ea[3]=0.0; //e12
	ea[4]=0.0; //e23
	ea[5]=0.0; //e13
	
	int ijk; //ijk=(i*Ny+j)*Nz+k;
	
	//----- ----- ----- ----- ----- -----
	const int fftsizex = Nx, fftsizey = Ny, fftsizez = Nz;
	const int scale=fftsizex*fftsizey*fftsizez;
	double fftw3d_scale = (double)scale;
	/* fftw_complex *in, *out; // in[i][0] for real, in[i][1] for imag.
	    in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	   out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	   fftw_plan plan, iplan;
	   plan = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);	//For forward FFT
	  iplan = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);	//For inverse FFT */
	//
	//array
	fftw_complex *cr, *crk;
	 cr = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	crk = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	fftw_plan plan_cr, iplan_crk;
	 plan_cr  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, cr, crk, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_crk = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, crk,cr,  FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//array
	fftw_complex *dfdcr, *dfdcrk;
	 dfdcr = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	dfdcrk = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	//fftw_plan plan_dfdcr, iplan_dfdcrk;
	fftw_plan plan_dfdcr;
	 plan_dfdcr  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, dfdcr, dfdcrk, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//iplan_dfdcrk = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, dfdcrk, dfdcr, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//array
	fftw_complex *delsdc, *delsdck;
	 delsdc = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	delsdck = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	//fftw_plan plan_delsdc, iplan_delsdck;
	fftw_plan plan_delsdc;
	 plan_delsdc  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, delsdc, delsdck, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//iplan_delsdck = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, delsdck,delsdc,  FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//----- ----- ----- ----- ----- -----
	
	//----- ----- ----- -----
	//stress (head s series) components
	fftw_complex *s11, *s22, *s33;
	 s11 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 s22 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 s33 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 //
	fftw_complex *s12, *s23, *s13;
	 s12 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 s23 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 s13 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	//----- ----- ----- -----
	//strain (head e series) components
	fftw_complex *e11, *e22, *e33;
	 e11 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 e22 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 e33 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 //
	fftw_complex *e12, *e23, *e13;
	 e12 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 e23 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 e13 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	//----- ----- ----- ----- ----- -----
	
	//Initialize stress & strain componentes
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				ijk=(i*Ny+j)*Nz+k;
				//----- ----- ----- ----- -----
				// stress (smatx, sXYk)
				s11[ijk][0] = s11[ijk][1] = 0.0;
				s22[ijk][0] = s22[ijk][1] = 0.0;
				s33[ijk][0] = s33[ijk][1] = 0.0;
				//
				s12[ijk][0] = s12[ijk][1] = 0.0;
				s23[ijk][0] = s23[ijk][1] = 0.0;
				s13[ijk][0] = s13[ijk][1] = 0.0;
				//----- ----- ----- ----- -----
				// strain (ematx, eXYk)
				e11[ijk][0] = e11[ijk][1] = 0.0;
				e22[ijk][0] = e22[ijk][1] = 0.0;
				e33[ijk][0] = e33[ijk][1] = 0.0;
				//
				e12[ijk][0] = e12[ijk][1] = 0.0;
				e23[ijk][0] = e23[ijk][1] = 0.0;
				e13[ijk][0] = e13[ijk][1] = 0.0;
				//----- ----- ----- ----- -----
			}
		}
	}
	
	//Strain components due to lattice defects
	//----- ----- ----- ----- ----- -----
	//double ed11[Nx][Ny][Nz];
	double *ed11 = (double *)malloc(sizeof(double)*( NxNyNz ));
	//double ed22[Nx][Ny][Nz];
	double *ed22 = (double *)malloc(sizeof(double)*( NxNyNz ));
	//double ed33[Nx][Ny][Nz];
	double *ed33 = (double *)malloc(sizeof(double)*( NxNyNz ));
	//
	//double ed12[Nx][Ny][Nz];
	double *ed12 = (double *)malloc(sizeof(double)*( NxNyNz ));
	//double ed23[Nx][Ny][Nz];
	double *ed23 = (double *)malloc(sizeof(double)*( NxNyNz ));
	//double ed13[Nx][Ny][Nz];
	double *ed13 = (double *)malloc(sizeof(double)*( NxNyNz ));
	//----- ----- ----- ----- ----- -----
	
	//dislocation eigen strain
	/* idislo=1 for dislocation diploe,
	   idislo=2 for dislocation array */
	int idislo=1;
	dislo_strain_3d(Nx,Ny,Nz,idislo,
					ed11,ed22,ed33,
					ed12,ed23,ed13);
	
	//int iflag=1;
	
	//prepare microstructure
	micro_ch_pre_3d(Nx,Ny,Nz,c0,cr); //Initialize microstructure
	
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
	
	//double omeg11[Nx][Ny][Nz], Coefficient needed for the Green's tensor
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	double *omeg11 = (double *)malloc(sizeof(double)*( NxNyNz ));
	double *omeg22 = (double *)malloc(sizeof(double)*( NxNyNz ));
	double *omeg33 = (double *)malloc(sizeof(double)*( NxNyNz ));
	//
	double *omeg12 = (double *)malloc(sizeof(double)*( NxNyNz ));
	double *omeg23 = (double *)malloc(sizeof(double)*( NxNyNz ));
	double *omeg13 = (double *)malloc(sizeof(double)*( NxNyNz ));
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	
	//Greens tensor
	green_tensor1_3D(Nx,Ny,Nz,
					 kx,ky,kz,
					 cm11,cm12,cm44,
					 cp11,cp12,cp44,
					 omeg11,omeg22,omeg33,
					 omeg12,omeg23,omeg13);
	
	//----- -----
	// semi-implicit scheme
	double numer;
	double denom;
	//----- -----
	
	//double cr_out[Nx][Ny][Nz];
	double *cr_out = (double *)malloc(sizeof(double)*( NxNyNz ));
	
	//evolve (Time evolution of microstructure)
	for(int istep=0;istep<=nstep;istep++){
		
		//Update the total time
		ttime=ttime+dtime;
		
		//derivative of elastic energy
		//Calculate the derivative of elastic energy
		solve_elasticity_3d(Nx,Ny,Nz,
			kx,ky,kz,
			omeg11,omeg22,omeg33,
			omeg12,omeg23,omeg13,
			s11,s22,s33,
			s12,s23,s13,
			e11,e22,e33,
			e12,e23,e13,
			ed11,ed22,ed33,
			ed12,ed23,ed13,
			cm11,cm12,cm44,
			cp11,cp12,cp44,
			ea,ei0,
			cr, delsdc);
		// Note: tmatx is real part only
		
		//derivative of chemical energy
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					//Normalize the derivative elastic energy with RT
					delsdc[ijk][0]=delsdc[ijk][0]/RT;
					delsdc[ijk][1]=0.0;
					//Calculate derivative of chemical energy
					dfdcr[ijk][0]=FeCr_chem_poten_3d((double)cr[ijk][0],tempr);
					dfdcr[ijk][1]=0.0;
				}
			}
		}
		
		/* Take the values of concentration, derivative of free energy and
		   derivative of elastic energy from real space to Fourier space (forward FFT) */
		//crk=fft3(cr);
		fftw_execute(plan_cr);
		//dfdcrk=fft3(dfdcr);
		fftw_execute(plan_dfdcr);
		//delsdck=fft3(delsdc);
		fftw_execute(plan_delsdc);
		
		/* Semi-implicit time integration of Cr concentration field at
		   Fourier space (Eq.5.50) */
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					//
					denom=1.0+dtime*coefA*mobility*grad_coef*k4[ijk];
					//
					numer=dtime*mobility*k2[ijk]*(dfdcrk[ijk][0]+delsdck[ijk][0]);
					crk[ijk][0]=(crk[ijk][0]-numer)/denom;
					//
					numer=dtime*mobility*k2[ijk]*(dfdcrk[ijk][1]+delsdck[ijk][1]);
					crk[ijk][1]=(crk[ijk][1]-numer)/denom;
				}
			}
		}
		
		/* Take concentration field from Fourier space back to
		   real space (inverse FFT) */
		//cr=real(ifft3(crk));
		fftw_execute(iplan_crk);
		
		//for small deviations
		// For small deviations from max and min values, reset the limits
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					//----- ----- ----- -----
					 cr[ijk][0] /= fftw3d_scale;
					 cr[ijk][1] /= fftw3d_scale;
					//----- ----- ----- -----
					if(cr[ijk][0]>=0.9999){
						cr[ijk][0]=0.9999;
					}
					if(cr[ijk][0]<=0.0001){
						cr[ijk][0]=0.0001;
					}
					cr[ijk][1]=0.0;
					//----- ----- ----- -----
					cr_out[ijk]  = cr[ijk][0];
					//----- ----- ----- -----
				}
			}
		}
		
		//print results
		/* If print frequency is reached, output the results to file */
		if(fmod(istep,nprint)==0){
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_3D(Nx,Ny,Nz,dx,dy,dz,istep,cr_out);
			
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
	//----- ----- ----- ----- -----
	fftw_free(s11);
	fftw_free(s22);
	fftw_free(s33);
	//
	fftw_free(s12);
	fftw_free(s23);
	fftw_free(s13);
	//----- ----- ----- ----- -----
	//fftw_free(e11);
	//fftw_free(e22);
	//fftw_free(e33);
	//
	//fftw_free(e12);
	//fftw_free(e23);
	//fftw_free(e13);
	//----- ----- ----- ----- ----- -----
	//free(ed11);
	//free(ed22);
	//free(ed33);
	//
	//free(ed12);
	//free(ed23);
	//free(ed13);
	//----- ----- ----- ----- -----
	free(kx);
	free(ky);
	free(k2);
	free(k4);
	//----- ----- ----- ----- -----
	free(cr_out);
	//----- ----- ----- ----- ----- -----
}
