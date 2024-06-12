/* 2D semi-implicit spectral
  phase-field crystal code */

#include <stdio.h>
#include <stdlib.h> //rand()
#include <math.h> //mod() and -lm
#include <time.h>

//#include <fftw3.h>
//gcc test.c -lfftw3 -lm
#include <mpi.h> //mpi version
#include <fftw3-mpi.h> //mpi version
//mpicc test.c -lfftw3_mpi -lfftw3 -lm
//Note: need -lfftw3 for " undefined reference to symbol 'fftw_malloc'"

void prepare_fft_2d(int Nx, int Ny, 
	double dx, double dy,
	double *kx, double *ky, 
	double *k2, double *k4);

void write_vtk_grid_values_2D(int nx, int ny, 
	double dx, double dy, double dx0,
	int istep, double *data1, double *data2);

int main(int argc, char **argv){
	clock_t start, end;
	double compute_time;
	
	//get initial wall time
	start = clock();
	
	FILE *out1=fopen("final_conf.out","w");
	FILE *out2=fopen("energy.out","w");
	
	//simulation cell parameters
	int Nx=512;
	int Ny=512;
	
	//The value of pi
	double pix=4.0*atan(1.0);
	
	//The distance between two grid points in x,y,z-direction
	double dx=pix/4.0;
	double dy=pix/4.0;
	
	double dx0=dx;
	//double dy0=dy;
	
	//time integration parameters
	int nstep=350000;
	int nprint=5000;
	double dtime=0.01;
	
	int ndefor=1000;
	double stran=0.0;
	
	//material specific parameters
	//double den0=-0.085; //average density for pfc_3d_v1 (stripe phase)
	double den0=-0.285; //average density for pfc_3d_v2 (triangular phase)
	double tempr=-0.5; //temperature (T-Tm)/Tm, Tm=melting point
	//double tempr0=tempr; //positive value, tempr=tempr+tempr0/isteps;
	double noise=den0*1e-2; //Noise term to modulate the initial density field
	
	int ii;
	int iimpi;
	
	//double den_out[Nx][Ny];
	double *den_out = (double *)malloc(sizeof(double)*( Nx*Ny ));
	
	//if infile=1 read input from file
	int infile = 1;
	FILE *in1;
	int mx,my;
	if(infile==1){
		//open input file
		in1=fopen("bi_2r_2d.inp","r");
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				fscanf(in1,"%5d %5d %lf",&mx,&my,&den_out[ii]);
			}
		}
	}else{
		//initialize density
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				den_out[ii] = den0 + noise*(0.5-(double)rand()/RAND_MAX);
			}
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
	prepare_fft_2d(Nx,Ny,dx,dy,kx,ky,k2,k4); //get FFT coefficients
	
	//double Linx[Nx][Ny];
	double *Linx = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double denom[Nx][Ny];
	double *denom = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//
	double energy=0.0;
	double energy0=0.0;
	//double ss2[Nx][Ny];
	double *ss2 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double ss4[Nx][Ny];
	double *ss4 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//
	//double ff_out[Nx][Ny];
	double *ff_out = (double *)malloc(sizeof(double)*( Nx*Ny ));
	
	//----- ----- ----- ----- ----- -----
	const ptrdiff_t fftsizex = Nx, fftsizey = Ny;
	ptrdiff_t alloc_local, local_n0, local_0_start;
	
	int rank;
	//int num_proc;
	
	MPI_Init(&argc, &argv);
	//MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	fftw_mpi_init();
	
	alloc_local = fftw_mpi_local_size_2d(fftsizex, fftsizey, MPI_COMM_WORLD, &local_n0, &local_0_start);
	
	//array
	fftw_complex *den, *f_den;
	  den = fftw_alloc_complex(alloc_local);
	f_den = fftw_alloc_complex(alloc_local);
	fftw_plan plan_den, iplan_f_den;
	   plan_den = fftw_mpi_plan_dft_2d(fftsizex, fftsizey, den, f_den, MPI_COMM_WORLD, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_f_den = fftw_mpi_plan_dft_2d(fftsizex, fftsizey, f_den, den, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	fftw_complex *den3, *f_den3;
	  den3 = fftw_alloc_complex(alloc_local);
	f_den3 = fftw_alloc_complex(alloc_local);
	fftw_plan plan_den3;
	plan_den3 = fftw_mpi_plan_dft_2d(fftsizex, fftsizey, den3, f_den3, MPI_COMM_WORLD, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//
	fftw_complex *Nonx;
	  Nonx = fftw_alloc_complex(alloc_local);
	//
	fftw_complex *ff, *f_ff;
	  ff = fftw_alloc_complex(alloc_local);
	f_ff = fftw_alloc_complex(alloc_local);
	fftw_plan iplan_f_ff;
	iplan_f_ff = fftw_mpi_plan_dft_2d(fftsizex, fftsizey, f_ff, ff, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE);	//For inverse FFT
	//----- ----- ----- ----- ----- -----
	
	for(int i=0;i<local_n0;i++){
		for(int j=0;j<Ny;j++){
			ii=i*Ny+j;
			iimpi=(local_0_start+i)*Ny+j;
			den[ii][0] = den_out[iimpi];
			den[ii][1] = 0.0;
		}
	}
	
	//evolve (evolve microstructure)
	for(int istep=0;istep<=nstep;istep++){
		
		//tempr = tempr + tempr0/istep;
		
		//loding
		if(fmod(istep,ndefor)==0){
			dx=dx+5.0e-4*dx0;
			dy=dy-5.0e-4*dx0;
			
			stran=stran+5.0e-4*dx0;
			
			prepare_fft_2d(Nx,Ny,dx,dy,kx,ky,k2,k4); //get FFT coefficients
		}
		
		//take current density filed from real space to Fourier space (forward FFT transformation)
		//f_den = fftn(den);
		fftw_execute(plan_den);
		
		//calculate the value of denominator in Eq.7.10 at every grid points
		for(int i=0;i<local_n0;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				iimpi=(local_0_start+i)*Ny+j;
				Linx[iimpi]=-k2[iimpi]*(tempr+1.0-2.0*k2[iimpi]+k4[iimpi]);
				denom[iimpi]=1.0-dtime*Linx[iimpi];
			}
		}
		
		//calculate the nonlinear term, phi^3, in Eq.7.10
		for(int i=0;i<local_n0;i++){
			for(int j=0;j<Ny;j++){
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
		for(int i=0;i<local_n0;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				iimpi=(local_0_start+i)*Ny+j;
				Nonx[ii][0] = -k2[iimpi]*f_den3[ii][0];
				Nonx[ii][1] = -k2[iimpi]*f_den3[ii][1];
				f_den[ii][0] = (f_den[ii][0]+dtime*Nonx[ii][0])/denom[iimpi];
				f_den[ii][1] = (f_den[ii][1]+dtime*Nonx[ii][1])/denom[iimpi];
			}
		}
		
		//bring back the values of phi^(t+1) from Fourier space to real space (inverse FFT transformation)
		//den=real(ifftn(f_den));
		fftw_execute(iplan_f_den);
		
		//print results
		//if print frequency is reached, output the results to file
		if(fmod(istep,nprint)==0){
		
			//printf("done step: %5d \n", istep);
			
			//energy calculation
			//calculate the free energy distribution, Eq.7.6
			for(int i=0;i<local_n0;i++){
				for(int j=0;j<Ny;j++){
					ii=i*Ny+j;
					iimpi=(local_0_start+i)*Ny+j;
					den_out[iimpi]=den[ii][0]/(fftsizex*fftsizey);
					ss2[iimpi]=den_out[iimpi]*den_out[iimpi];
					ss4[iimpi]=ss2[iimpi]*ss2[iimpi];
				}
			}
			
			for(int i=0;i<local_n0;i++){
				for(int j=0;j<Ny;j++){
					ii=i*Ny+j;
					f_ff[ii][0] = 0.5*f_den[ii][0]*(1.0-2.0*k2[iimpi]+k4[iimpi]);
					f_ff[ii][1] = 0.5*f_den[ii][1]*(1.0-2.0*k2[iimpi]+k4[iimpi]);
				}
			}
			
			//ff=real(ifftn(f_ff));
			fftw_execute(iplan_f_ff);
			
			for(int i=0;i<local_n0;i++){
				for(int j=0;j<Ny;j++){
					ii=i*Ny+j;
					iimpi=(local_0_start+i)*Ny+j;
					ff[ii][0] = (ff[ii][0]/(fftsizex*fftsizey))
							  *(den[ii][0]/(fftsizex*fftsizey))
					+ 0.5*tempr*ss2[iimpi]
					+ 0.25*ss4[iimpi];
					ff_out[ii]=ff[ii][0];
				}
			}
			
			//integrate the free energy field
			energy = 0.0;
			
			for(int i=0;i<local_n0;i++){
				for(int j=0;j<Ny;j++){
					//energy = energy + creal(ff[i][j]);
					ii=i*Ny+j;
					energy = energy + ff[ii][0];
				}
			}
			
			//average free energy density
			energy = energy/(Nx*Ny);
			
			if(istep==1){
				energy0=energy;
			}
			
			energy=energy-energy0;
			
			MPI_Gather(&den_out[local_0_start*Ny], local_n0*Ny, MPI_DOUBLE, &den_out[local_0_start*Ny], local_n0*Ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			//
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			if (rank == 0){
				//print the average free energy density value to file
				fprintf(out2, "%14.6e %14.6e \n",stran, energy);
				
				//output the results in vtk file format for contour plots to be viewed by using paraview
				write_vtk_grid_values_2D(Nx,Ny,dx,dy,dx0,istep,den_out,ff_out);
				
				printf("done step: %5d \n", istep);
			}
		}
		
		//if intermediate configuration files are required, print the density field to file
		if(istep==nstep){
			MPI_Gather(&den_out[local_0_start*Ny], local_n0*Ny, MPI_DOUBLE, &den_out[local_0_start*Ny], local_n0*Ny, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			//
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			if (rank == 0){
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						ii=i*Ny+j;
						fprintf(out1,"%5d %5d %14.6e \n",i,j,den_out[ii]);
					}
				}
			}
		}
		
		//for recycle
		for(int i=0;i<local_n0;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				den[ii][0] = den[ii][0]/(fftsizex*fftsizey);
				den[ii][1] = den[ii][1]/(fftsizex*fftsizey);
			}
		}
		
	}//end of time step (evolve,for)
	
	if (rank == 0){
		//calculate the execution time and print it
		end = clock();
		compute_time = ((double) (end - start)) / CLOCKS_PER_SEC;
		printf("Compute Time: %lf \n", compute_time);
	}
	
	//----- ----- ----- ----- ----- -----
	fftw_destroy_plan(plan_den);
	fftw_destroy_plan(plan_den3);
	fftw_destroy_plan(iplan_f_den);
	fftw_destroy_plan(iplan_f_ff);
	//----- ----- ----- ----- ----- -----
	fftw_free(den);
	fftw_free(f_den);
	fftw_free(den3);
	fftw_free(f_den3);
	fftw_free(ff);
	fftw_free(f_ff);
	//----- ----- ----- ----- ----- -----
	fclose(out1);
	fclose(out2);
	//----- ----- ----- ----- ----- -----
	MPI_Finalize(); //for fftw3_mpi
}
