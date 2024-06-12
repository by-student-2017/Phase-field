/* 3D semi-implicit spectral
  phase-field crystal code */

#include <stdio.h>
#include <stdlib.h> //rand()
#include <math.h> //mod() and -lm
#include <time.h>

#include <fftw3.h>
//gcc test.c -lfftw3
//#include <mpi.h> //mpi version
//#include <fftw3-mpi.h> //mpi version

void prepare_fft_3d(int Nx, int Ny, int Nz, 
	double dx, double dy, double dz,
	double *kx, double *ky, double *kz, 
	double *k2, double *k4);

void write_vtk_grid_values_3D(int nx, int ny, int nz, 
	double dx, double dy, double dz, double dx0,
	int istep, double *data1, double *data2);

int main(){
	clock_t start, end;
	double compute_time;
	
	//get initial wall time
	start = clock();
	
	FILE *out1=fopen("final_conf.out","w");
	FILE *out2=fopen("energy.out","w");
	
	//simulation cell parameters
	int Nx=512;
	int Ny=512;
	int Nz=2;
	
	//The value of pi
	double pix=4.0*atan(1.0);
	
	//The distance between two grid points in x,y,z-direction
	double dx=pix/4.0;
	double dy=pix/4.0;
	double dz=pix/4.0;
	
	double dx0=dx;
	//double dy0=dy;
	
	//time integration parameters
	//int nstep=20000; //for pfc_3d_v1
	int nstep=200000; //for pfc_3d_v2
	int nprint=5000;
	double dtime=0.05;
	
	int ndefor=1000;
	double stran=0.0;
	
	//material specific parameters
	//double den0=-0.085; //average density for pfc_3d_v1 (stripe phase)
	double den0=-0.285; //average density for pfc_3d_v2 (triangular phase)
	double tempr=-0.25; //temperature (T-Tm)/Tm, Tm=melting point
	//double tempr0=0.0; //positive value, tempr=tempr+tempr0/isteps;
	double noise=den0*1e-2; //Noise term to modulate the initial density field
	
	int ii;
	
	const int fftsizex = Nx, fftsizey = Ny, fftsizez = Nz;
	/* fftw_complex *in, *out; // in[i][0] for real, in[i][1] for imag.
	   in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey*fftsizez);
	   out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey*fftsizez);
	   fftw_plan plan, iplan;
	   plan = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);	//For forward FFT
	  iplan = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);	//For inverse FFT */
	
	//array
	fftw_complex *den, *f_den;
	den = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey*fftsizez);
	f_den = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey*fftsizez);
	fftw_plan plan_den, iplan_f_den;
	plan_den = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, den, f_den, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_f_den = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, f_den, den, FFTW_BACKWARD, FFTW_ESTIMATE);	//For inverse FFT
	//
	fftw_complex *den3, *f_den3;
	den3 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey*fftsizez);
	f_den3 = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey*fftsizez);
	fftw_plan plan_den3;
	plan_den3 = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, den3, f_den3, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//
	//double den_out[Nx][Ny][Nz];
	double *den_out = (double *)malloc(sizeof(double)*( Nx*Ny*Nz ));
	//
	//double kx[Nx];
	double *kx = (double *)malloc(sizeof(double)*( Nx ));
	//double ky[Ny];
	double *ky = (double *)malloc(sizeof(double)*( Ny ));
	//double kz[Nz];
	double *kz = (double *)malloc(sizeof(double)*( Nz ));
	//double k2[Nx][Ny][Nz];
	double *k2 = (double *)malloc(sizeof(double)*( Nx*Ny*Nz ));
	//double k4[Nx][Ny][Nz];
	double *k4 = (double *)malloc(sizeof(double)*( Nx*Ny*Nz ));
	//
	//double Linx[Nx][Ny][Nz];
	double *Linx = (double *)malloc(sizeof(double)*( Nx*Ny*Nz ));
	//double denom[Nx][Ny][Nz];
	double *denom = (double *)malloc(sizeof(double)*( Nx*Ny*Nz ));
	//
	fftw_complex *Nonx;
	Nonx = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey*fftsizez);
	//
	double energy=0.0;
	double energy0=0.0;
	//double ss2[Nx][Ny][Nz];
	double *ss2 = (double *)malloc(sizeof(double)*( Nx*Ny*Nz ));
	//double ss4[Nx][Ny][Nz];
	double *ss4 = (double *)malloc(sizeof(double)*( Nx*Ny*Nz ));
	//
	fftw_complex *ff, *f_ff;
	ff = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey*fftsizez);
	f_ff = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey*fftsizez);
	fftw_plan iplan_f_ff;
	iplan_f_ff = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, f_ff, ff, FFTW_BACKWARD, FFTW_ESTIMATE);	//For inverse FFT
	//
	//double ff_out[Nx][Ny][Nz];
	double *ff_out = (double *)malloc(sizeof(double)*( Nx*Ny*Nz ));
	
	//if infile=1 read input from file
	int infile = 1;
	FILE *in1;
	int mx,my,mz;
	if(infile==1){
		//open input file
		in1=fopen("bi_2r_3d.inp","r");
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ii=i*Ny*Nz+j*Nz+k;
					fscanf(in1,"%5d %5d %5d %lf",&mx,&my,&mz,&den[ii][0]);
					den[ii][1] = 0.0;
				}
			}
		}
	}else{
		//initialize density
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					//den[i][j][k]=den0+noise*(0.5-(double)rand()/RAND_MAX);
					ii=i*Ny*Nz+j*Nz+k;
					den[ii][0] = den0 + noise*(0.5-(double)rand()/RAND_MAX);
					den[ii][1] = 0.0;
					//modulate the density field with given noise term
				}
			}
		}
	}
	
	//prepare fft (output: kx,ky,kz,k2,k4)
	prepare_fft_3d(Nx,Ny,Nz,dx,dy,dz,kx,ky,kz,k2,k4); //get FFT coefficients
	
	//evolve (evolve microstructure)
	for(int istep=0;istep<=nstep;istep++){
		
		//tempr = tempr + tempr0/istep;
		
		//loding
		if(fmod(istep,ndefor)==0){
			dx=dx+5.0e-4*dx0;
			dy=dy-5.0e-4*dx0;
			
			stran=stran+5.0e-4*dx0;
			
			prepare_fft_3d(Nx,Ny,Nz,dx,dy,dz,kx,ky,kz,k2,k4); //get FFT coefficients
		}
		
		//take current density filed from real space to Fourier space (forward FFT transformation)
		//f_den = fftn(den);
		fftw_execute(plan_den);
		
		//calculate the value of denominator in Eq.7.10 at every grid points
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ii=i*Ny*Nz+j*Nz+k;
					Linx[ii]=-k2[ii]*(tempr+1.0-2.0*k2[ii]+k4[ii]);
					denom[ii]=1.0-dtime*Linx[ii];
				}
			}
		}
		
		//calculate the nonlinear term, phi^3, in Eq.7.10
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					//den3[i][j][k]=den[i][j][k]*den[i][j][k]*den[i][j][k];
					ii=i*Ny*Nz+j*Nz+k;
					den3[ii][0] =   den[ii][0]*den[ii][0]*den[ii][0] -3*den[ii][0]*den[ii][1]*den[ii][1];
					den3[ii][1] = 3*den[ii][0]*den[ii][0]*den[ii][1]   -den[ii][1]*den[ii][1]*den[ii][1];
					//(a+bj)*(a+bj)*(a+bj)=(a*a+2*a*bj-b*b)*(a+bj)
					//= a*a*a -3*a*b*b + 3*a*a*bj -b*b*bj
				}
			}
		}
		
		//take the value of phi^3 from real space to Fourier space (forward FFT transformation)
		//f_den3=fftn(den3);
		fftw_execute(plan_den3);
		
		//calculate the value of phi^(t+1) from Fourier space to real space (inverse FFT transformation)
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					//Nonx[i][j][k]=-k2[i][j][k]*f_den3[i][j][k];
					//f_den[i][j][k]=(f_den[i][j][k]+dtime*Nonx[i][j][k])/denom[i][j][k];
					ii=i*Ny*Nz+j*Nz+k;
					Nonx[ii][0] = -k2[ii]*f_den3[ii][0];
					Nonx[ii][1] = -k2[ii]*f_den3[ii][1];
					f_den[ii][0] = (f_den[ii][0]+dtime*Nonx[ii][0])/denom[ii];
					f_den[ii][1] = (f_den[ii][1]+dtime*Nonx[ii][1])/denom[ii];
				}
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
					for(int k=0;k<Nz;k++){
						//ss2[i][j][k]=creal(den[i][j][k])*creal(den[i][j][k]);
						ii=i*Ny*Nz+j*Nz+k;
						den_out[ii]=den[ii][0]/(fftsizex*fftsizey*fftsizez);
						ss2[ii]=den_out[ii]*den_out[ii];
						ss4[ii]=ss2[ii]*ss2[ii];
					}
				}
			}
			
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					for(int k=0;k<Nz;k++){
						//f_ff[i][j][k]=0.5*f_den[i][j][k]*(1.0-2.0*k2[i][j][k]+k4[i][j][k]);
						ii=i*Ny*Nz+j*Nz+k;
						f_ff[ii][0] = 0.5*f_den[ii][0]*(1.0-2.0*k2[ii]+k4[ii]);
						f_ff[ii][1] = 0.5*f_den[ii][1]*(1.0-2.0*k2[ii]+k4[ii]);
					}
				}
			}
			
			//ff=real(ifftn(f_ff));
			fftw_execute(iplan_f_ff);
			
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					for(int k=0;k<Nz;k++){
						//f_ff[i][j][k]=creal(ff[i][j][k])*creal(den[i][j][k])
						ii=i*Ny*Nz+j*Nz+k;
						ff[ii][0] = (ff[ii][0]/(fftsizex*fftsizey*fftsizez))
								  *(den[ii][0]/(fftsizex*fftsizey*fftsizez))
						+ 0.5*tempr*ss2[ii]
						+ 0.25*ss4[ii];
						ff_out[ii]=ff[ii][0];
					}
				}
			}
			
			//integrate the free energy field
			energy = 0.0;
			
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					for(int k=0;k<Nz;k++){
						//energy = energy + creal(ff[i][j][k]);
						ii=i*Ny*Nz+j*Nz+k;
						energy = energy + ff[ii][0];
					}
				}
			}
			
			//average free energy density
			energy = energy/(Nx*Ny*Nz);
			
			if(istep==nprint){
				energy0=energy;
			}
			
			energy=energy-energy0;
			
			//print the average free energy density value to file
			fprintf(out2, "%14.6e %14.6e \n",stran, energy);
			
			//output the results in vtk file format for contour plots to be viewed by using paraview
			write_vtk_grid_values_3D(Nx,Ny,Nz,dx,dy,dz,dx0,istep,den_out,ff_out);
			
			printf("done step: %5d \n", istep);
		}
		
		//if intermediate configuration files are required, print the density field to file
		if(istep==nstep){
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					for(int k=0;k<Nz;k++){
						ii=i*Ny*Nz+j*Nz+k;
						fprintf(out1,"%5d %5d %5d %14.6e \n",i,j,k,den_out[ii]);
					}
				}
			}
		}
		
		//for recycle
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ii=i*Ny*Nz+j*Nz+k;
					den[ii][0] = den[ii][0]/(fftsizex*fftsizey*fftsizez);
					den[ii][1] = den[ii][1]/(fftsizex*fftsizey*fftsizez);
				}
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
