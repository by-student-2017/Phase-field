/* 2D semi-implicit spectral
  phase-field crystal code */

#include <stdio.h>
#include <stdlib.h> //rand()
#include <math.h> //mod() and -lm
#include <time.h>

#include <fftw3.h>
//gcc test.c -lfftw3
//#include <mpi.h> //mpi version
//#include <fftw3-mpi.h> //mpi version

void prepare_fft_2d(int Nx, int Ny, 
	double dx, double dy,
	double *kx, double *ky, 
	double *k2, double *k4);

void write_vtk_grid_values_2D(int nx, int ny, 
	double dx, double dy,
	int istep, double *data1);

int main(){
	clock_t start, end;
	double compute_time;
	
	//get initial wall time
	start = clock();
	
	FILE *out1=fopen("final_conf.out","w");
	FILE *out2=fopen("energy.out","w");
	
	//simulation cell parameters
	int Nx=64;
	int Ny=64;
	
	//The value of pi
	double pix=4.0*atan(1.0);
	
	//The distance between two grid points in x,y,z-direction
	double dx=pix/4.0;
	double dy=pix/4.0;
	
	//time integration parameters
	int nstep=200000;
	int nprint=2000;
	double dtime=0.05;
	
	//material specific parameters
	//double den0=-0.085; //average density for pfc_3d_v1 (stripe phase)
	double den0=-0.285; //average density for pfc_3d_v2 (triangular phase)
	double tempr=-0.25; //temperature (T-Tm)/Tm, Tm=melting point
	//double tempr0=0.0; //positive value, tempr=tempr+tempr0/isteps;
	double noise=den0*1e-2; //Noise term to modulate the initial density field
	
	int ii;
	
	const int fftsizex = Nx, fftsizey = Ny;
	/* fftw_complex *in, *out; // in[i][0] for real, in[i][1] for imag.
	   in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	   out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	   fftw_plan plan, iplan;
	   plan = fftw_plan_dft_2d(fftsize, fftsize, in, out, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	  iplan = fftw_plan_dft_2d(fftsize, fftsize, in, out, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT */
	
	//array
	fftw_complex *den, *f_den;
	den = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	f_den = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	fftw_plan plan_den, iplan_f_den;
	plan_den = fftw_plan_dft_2d(fftsizex, fftsizey, den, f_den, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_f_den = fftw_plan_dft_2d(fftsizex, fftsizey, f_den, den, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	fftw_complex *den3, *f_den3;
	den3 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	f_den3 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	fftw_plan plan_den3;
	plan_den3 = fftw_plan_dft_2d(fftsizex, fftsizey, den3, f_den3, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//
	//double den_out[Nx][Ny];
	double *den_out = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//
	//double kx[Nx];
	double *kx = (double *)malloc(sizeof(double)*( Nx ));
	//double ky[Ny];
	double *ky = (double *)malloc(sizeof(double)*( Ny ));
	//double k2[Nx][Ny];
	double *k2 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double k4[Nx][Ny];
	double *k4 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//
	//double Linx[Nx][Ny];
	double *Linx = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double denom[Nx][Ny];
	double *denom = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//
	fftw_complex *Nonx;
	Nonx = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	//
	double energy;
	//double ss2[Nx][Ny];
	double *ss2 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double ss4[Nx][Ny];
	double *ss4 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//
	fftw_complex *ff, *f_ff;
	ff = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	f_ff = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	fftw_plan iplan_f_ff;
	iplan_f_ff = fftw_plan_dft_2d(fftsizex, fftsizey, f_ff, ff, FFTW_BACKWARD, FFTW_ESTIMATE);	//For inverse FFT
	
	//if infile=1 read input from file
	int infile = 0;
	FILE *in1;
	int mx,my;
	if(infile==1){
		//open input file
		in1=fopen("g3_2r.inp","r");
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				fscanf(in1,"%5d %5d %lf",&mx,&my,&den[ii][0]);
				den[ii][1] = 0.0;
			}
		}
	}else{
		//initialize density
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				//den[i][j]=den0+noise*(0.5-(double)rand()/RAND_MAX);
				ii=i*Ny+j;
				den[ii][0] = den0 + noise*(0.5-(double)rand()/RAND_MAX);
				den[ii][1] = 0.0;
			}
		}
	}
	
	//prepare fft (output: kx,ky,kz,k2,k4)
	prepare_fft_2d(Nx,Ny,dx,dy,kx,ky,k2,k4); //get FFT coefficients
	
	//evolve (evolve microstructure)
	for(int istep=0;istep<=nstep;istep++){
		
		//tempr = tempr + tempr0/istep;
		
		//take current density filed from real space to Fourier space (forward FFT transformation)
		//f_den = fftn(den);
		fftw_execute(plan_den);
		
		//calculate the value of denominator in Eq.7.10 at every grid points
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				Linx[ii]=-k2[ii]*(tempr+1.0-2.0*k2[ii]+k4[ii]);
				denom[ii]=1.0-dtime*Linx[ii];
			}
		}
		
		//calculate the nonlinear term, phi^3, in Eq.7.10
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				//den3[i][j]=den[i][j]*den[i][j]*den[i][j];
				ii=i*Ny+j;
				den3[ii][0] =   den[ii][0]*den[ii][0]*den[ii][0] -3*den[ii][0]*den[ii][1]*den[ii][1];
				den3[ii][1] = 3*den[ii][0]*den[ii][0]*den[ii][1]   -den[ii][1]*den[ii][1]*den[ii][1];
				//(a+bj)*(a+bj)*(a+bj)=(a*a+2*a*bj-b*b)*(a+bj)
				//= a*a*a -3*a*b*b + 3*a*a*bj -b*b*bj
			}
		}
		
		//take the value of phi^3 from real space to Fourier space (forward FFT transformation)
		//f_den3=fftn(den3);
		fftw_execute(plan_den3);
		
		//calculate the value of phi^(t+1) from Fourier space to real space (inverse FFT transformation)
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				//Nonx[i][j]=-k2[i][j]*f_den3[i][j];
				//f_den[i][j]=(f_den[i][j]+dtime*Nonx[i][j])/denom[i][j];
				ii=i*Ny+j;
				Nonx[ii][0] = -k2[ii]*f_den3[ii][0];
				Nonx[ii][1] = -k2[ii]*f_den3[ii][1];
				f_den[ii][0] = (f_den[ii][0]+dtime*Nonx[ii][0])/denom[ii];
				f_den[ii][1] = (f_den[ii][1]+dtime*Nonx[ii][1])/denom[ii];
			}
		}
		
		//bring back the values of phi^(t+1) from Fourier space to real space (inverse FFT transformation)
		//den=real(ifftn(f_den));
		fftw_execute(iplan_f_den);
		
		//print results
		//if print frequency is reached, output the results to file
		if(fmod(istep,nprint)==0){
			
			//energy calculation
			//calculate the free energy distribution, Eq.7.6
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					//ss2[i][j]=creal(den[i][j])*creal(den[i][j]);
					ii=i*Ny+j;
					den_out[ii]=den[ii][0]/(fftsizex*fftsizey);
					ss2[ii]=den_out[ii]*den_out[ii];
					ss4[ii]=ss2[ii]*ss2[ii];
				}
			}
			
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					//f_ff[i][j]=0.5*f_den[i][j]*(1.0-2.0*k2[i][j]+k4[i][j]);
					ii=i*Ny+j;
					f_ff[ii][0] = 0.5*f_den[ii][0]*(1.0-2.0*k2[ii]+k4[ii]);
					f_ff[ii][1] = 0.5*f_den[ii][1]*(1.0-2.0*k2[ii]+k4[ii]);
				}
			}
			
			//ff=real(ifftn(f_ff));
			fftw_execute(iplan_f_ff);
			
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					//f_ff[i][j]=creal(ff[i][j])*creal(den[i][j])
					ii=i*Ny+j;
					ff[ii][0] = (ff[ii][0]/(fftsizex*fftsizey))
							  *(den[ii][0]/(fftsizex*fftsizey))
					+ 0.5*tempr*ss2[ii]
					+ 0.25*ss4[ii];
				}
			}
			
			//integrate the free energy field
			energy = 0.0;
			
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					//energy = energy + creal(ff[i][j]);
					ii=i*Ny+j;
					energy = energy + ff[ii][0];
				}
			}
			
			//average free energy density
			energy = energy/(Nx*Ny);
			
			//print the average free energy density value to file
			fprintf(out2, "%d %14.6e \n",istep, energy);
			
			//output the results in vtk file format for contour plots to be viewed by using paraview
			write_vtk_grid_values_2D(Nx,Ny,dx,dy,istep,den_out);
			
			printf("done step: %5d \n", istep);
		}
		
		//if intermediate configuration files are required, print the density field to file
		if(istep==nstep){
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					ii=i*Ny+j;
					fprintf(out1,"%5d %5d %14.6e \n",i,j,den_out[ii]);
				}
			}
		}
		
		//for recycle
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				den[ii][0] = den[ii][0]/(fftsizex*fftsizey);
				den[ii][1] = den[ii][1]/(fftsizex*fftsizey);
			}
		}
		
	}//end of time step (evolve,for)
	
	//calculate the execution time and print it
	end = clock();
	compute_time = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Compute Time: %lf \n", compute_time);
	
	fftw_destroy_plan(plan_den);
	fftw_destroy_plan(plan_den3);
	fftw_destroy_plan(iplan_f_den);
	fftw_destroy_plan(iplan_f_ff);
	//
	fftw_free(den);
	fftw_free(f_den);
	fftw_free(den3);
	fftw_free(f_den3);
	fftw_free(ff);
	fftw_free(f_ff);
	
	fclose(out1);
	fclose(out2);
}
