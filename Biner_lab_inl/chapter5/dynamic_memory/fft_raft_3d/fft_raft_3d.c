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

void micro_ch_pre_3d();
void prepare_fft_3d();
void green_tensor1_3D();
//void green_tensor2_3D();
double free_energy_ch_3d();
void solve_elasticity_3d();
void write_vtk_grid_values_3D();

int main(){
	clock_t start, end;
	double compute_time;
	
	//get initial wall time
	//(Get initial wall clock time beginning of the execution)
	start = clock();
	
	//simulation cell parameters
	int Nx=64;
	int Ny=64;
	int Nz=64;
	
	//Total number of grid points in the simulation cell
	int NxNyNz=Nx*Ny*Nz;
	
	//The distance between two grid points in x,y-direction
	double dx=1.0; // [nm] unit ?
	double dy=1.0; // [nm] unit ?
	double dz=1.0; // [nm] unit ?
	
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
	
	/* Note: elastic modulus in this case.
	double c22=c33=c11;
	double c21=c12;
	double c31=c13=c12;
	double c32=c23=c12;
	double c55=c66=c44;
	//
	double cm22=cm33=cm11;
	double cm21=cm12;
	double cm31=cm13=cm12;
	double cm32=cm23=cm12;
	double cm55=cm66=cm44;
	//
	double cp22=cp33=cp11;
	double cp21=cp12;
	double cp31=cp13=cp12;
	double cp32=cp23=cp12;
	double cp55=cp66=cp44; */
	
	//eigen strains
	double ei0=0.01; //Maginitude of eigenstrains
	
	//Applied strains
	double ea[6]; //Magnitude of applied strains
	ea[0]=0.00; //e11
	ea[1]=0.01; //e22
	ea[2]=0.00; //e33
	ea[3]=0.00; //e12
	ea[4]=0.00; //e23
	ea[5]=0.00; //e13
	
	int ijk;
	
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
	fftw_complex *con, *conk;
	 con = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	conk = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	fftw_plan plan_con, iplan_conk;
	 plan_con  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, con, conk, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_conk = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, conk,con,  FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//array
	fftw_complex *dfdcon, *dfdconk;
	 dfdcon = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	dfdconk = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	fftw_plan plan_dfdcon, iplan_dfdconk;
	 plan_dfdcon  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, dfdcon, dfdconk, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_dfdconk = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, dfdconk, dfdcon, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//array
	fftw_complex *delsdc, *delsdck;
	 delsdc = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	delsdck = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	fftw_plan plan_delsdc, iplan_delsdck;
	 plan_delsdc  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, delsdc, delsdck, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_delsdck = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, delsdck,delsdc,  FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
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
				//----- ----- ----- ----- -----
				//Strain components due to lattice defects
				ed11[ijk] = 0.0;
				ed22[ijk] = 0.0;
				ed33[ijk] = 0.0;
				//
				ed12[ijk] = 0.0;
				ed23[ijk] = 0.0;
				ed13[ijk] = 0.0;
				//----- ----- ----- ----- -----
			}
		}
	}
	
	//prepare microstructure
	micro_ch_pre_3d(Nx,Ny,Nz,c0,con); //Initialize microstructure
	
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
	
	//double con_out[Nx][Ny][Nz];
	double *con_out = (double *)malloc(sizeof(double)*( NxNyNz ));
	
	//evolve (Time evolution of microstructure)
	for(int istep=0;istep<=nstep;istep++){
		
		//Update the total time
		ttime=ttime+dtime;
		
		//derivative of free energy
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					//Calculate derivative of free energy
					dfdcon[ijk][0]=free_energy_ch_3d((double)con[ijk][0]);
					dfdcon[ijk][1]=0.0;
				}
			}
		}
		
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
			con, delsdc);
		// Note: tmatx is real part only
		
		/* Take the values of concentration, derivative of free energy and
		   derivative of elastic energy from real space to Fourier space (forward FFT) */
		//----- ----- ----- -----
		//conk=fft3(con);
		fftw_execute(plan_con);
		//dfdconk=fft3(dfdcon);
		fftw_execute(plan_dfdcon);
		//delsdck=fft3(delsdc);
		fftw_execute(plan_delsdc);
		//----- ----- ----- -----
		
		/* Semi-implicit time integration of concentration field at
		   Fourier space (Eq.5.50) */
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					//
					denom=1.0+dtime*coefA*mobility*grad_coef*k4[ijk]; // real part only
					//
					numer=dtime*mobility*k2[ijk]*(dfdconk[ijk][0]+delsdck[ijk][0]);
					conk[ijk][0]=(conk[ijk][0]-numer)/denom;
					//
					numer=dtime*mobility*k2[ijk]*(dfdconk[ijk][1]+delsdck[ijk][1]);
					conk[ijk][1]=(conk[ijk][1]-numer)/denom;
				}
			}
		}
		
		/* Take concentration field from Fourier space back to
		   real space (inverse FFT) */
		//----- ----- ----- -----
		//con=real(ifft3(conk));
		fftw_execute(iplan_conk);
		//----- ----- ----- -----
		
		//for small deviations
		// For small deviations from max and min values, reset the limits
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					//----- ----- ----- -----
					 con[ijk][0] /= fftw3d_scale;
					 con[ijk][1] /= fftw3d_scale;
					//----- ----- ----- -----
					if(con[ijk][0]>=0.9999){
						con[ijk][0]=0.9999;
					}
					if(con[ijk][0]<=0.0001){
						con[ijk][0]=0.0001;
					}
					con[ijk][1]=0.0;
					//----- ----- ----- -----
					con_out[ijk]  = con[ijk][0];
					//----- ----- ----- -----
				}
			}
		}
		
		//print results
		/* If print frequency is reached, output the results to file */
		if(fmod(istep,nprint)==0){
			printf("done step: %5d \n",istep);
			
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_3D(Nx,Ny,Nz,dx,dy,dz,istep,con_out);
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
	//----- ----- ----- ----- -----
	fftw_free(s11);
	fftw_free(s22);
	fftw_free(s33);
	//
	fftw_free(s12);
	fftw_free(s23);
	fftw_free(s13);
	//----- ----- ----- ----- -----
	fftw_free(e11);
	fftw_free(e22);
	fftw_free(e33);
	//
	fftw_free(e12);
	fftw_free(e23);
	fftw_free(e13);
	//----- ----- ----- ----- ----- -----
	free(ed11);
	free(ed22);
	free(ed33);
	//
	free(ed12);
	free(ed23);
	free(ed13);
	//----- ----- ----- ----- -----
	free(kx);
	free(ky);
	free(k2);
	free(k4);
	//----- ----- ----- ----- -----
	free(con_out);
	//----- ----- ----- ----- ----- -----
}
