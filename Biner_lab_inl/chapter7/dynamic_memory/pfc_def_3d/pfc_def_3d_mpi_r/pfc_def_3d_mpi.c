/* 3D semi-implicit spectral
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

void prepare_fft_3d(int Nx, int Ny, int Nz, 
	double dx, double dy, double dz,
	double *kx, double *ky, double *kz, 
	double *k2, double *k4);

void write_vtk_grid_values_3D(int nx, int ny, int nz, 
	double dx, double dy, double dz, double dx0,
	int istep, double *data1, double *data2);

int main(int argc, char **argv){
	clock_t start, end;
	double compute_time;
	
	//get initial wall time
	start = clock();
	
	FILE *out1=fopen("final_conf.out","w");
	FILE *out2=fopen("energy.out","w");
	
	//simulation cell parameters
	int Nx, Ny, Nz;
	//
	//if infile=1 read input from file
	int infile = 1;
	FILE *in1;
	if(infile==1){
		//open input file
		in1=fopen("bi_2r_3d.inp","r");
		fscanf(in1,"%5d %5d %5d ",&Nx,&Ny,&Nz);
	}else{
		Nx=64; Ny=64; Nz=1;
	}
	
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
	int iimpi;
	
	//double den_out[Nx][Ny][Nz];
	double *den_out = (double *)malloc(sizeof(double)*( Nx*Ny*Nz ));
	
	//if infile=1 read input from file
	//int infile = 1;
	//FILE *in1;
	int mx,my,mz;
	if(infile==1){
		//open input file
		//in1=fopen("bi_2r_3d.inp","r");
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ii=i*Ny*Nz+j*Nz+k;
					fscanf(in1,"%5d %5d %5d %lf",&mx,&my,&mz,&den_out[ii]);
				}
			}
		}
	}else{
		//initialize density
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ii=i*Ny*Nz+j*Nz+k;
					den_out[ii] = den0 + noise*(0.5-(double)rand()/RAND_MAX);
					//modulate the density field with given noise term
				}
			}
		}
	}
	
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
	
	//prepare fft (output: kx,ky,kz,k2,k4)
	prepare_fft_3d(Nx,Ny,Nz,dx,dy,dz,kx,ky,kz,k2,k4); //get FFT coefficients
	
	//double Linx[Nx][Ny][Nz];
	double *Linx = (double *)malloc(sizeof(double)*( Nx*Ny*Nz ));
	//double denom[Nx][Ny][Nz];
	double *denom = (double *)malloc(sizeof(double)*( Nx*Ny*Nz ));
	//
	double energy=0.0;
	double energy0=0.0;
	//double ss2[Nx][Ny][Nz];
	double *ss2 = (double *)malloc(sizeof(double)*( Nx*Ny*Nz ));
	//double ss4[Nx][Ny][Nz];
	double *ss4 = (double *)malloc(sizeof(double)*( Nx*Ny*Nz ));
	//
	//double ff_out[Nx][Ny][Nz];
	double *ff_out = (double *)malloc(sizeof(double)*( Nx*Ny*Nz ));
	
	//----- ----- ----- ----- ----- -----
	int rank;
	//int num_proc;
	
	MPI_Init(&argc, &argv);
	//MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	fftw_mpi_init();
	
	// Ref: https://www.fftw.org/fftw2_doc/fftw_4.html
	const ptrdiff_t fftsizex = Nx, fftsizey = Ny, fftsizez = Nz;
	ptrdiff_t size, local_n0, local_0_start;
	
	size = fftw_mpi_local_size_3d(fftsizex, fftsizey, fftsizez, MPI_COMM_WORLD, &local_n0, &local_0_start);
	
	//array
	fftw_complex *den, *f_den;
	  den = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
	f_den = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
	fftw_plan plan_den, iplan_f_den;
	   plan_den = fftw_mpi_plan_dft_3d(fftsizex, fftsizey, fftsizez, den, f_den, MPI_COMM_WORLD, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_f_den = fftw_mpi_plan_dft_3d(fftsizex, fftsizey, fftsizez, f_den, den, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE);	//For inverse FFT
	//
	fftw_complex *den3, *f_den3;
	  den3 = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
	f_den3 = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
	fftw_plan plan_den3;
	plan_den3 = fftw_mpi_plan_dft_3d(fftsizex, fftsizey, fftsizez, den3, f_den3, MPI_COMM_WORLD, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//
	fftw_complex *Nonx;
	  Nonx = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
	//
	fftw_complex *ff, *f_ff;
	  ff = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
	f_ff = (fftw_complex *)fftw_malloc(size*sizeof(fftw_complex));
	fftw_plan iplan_f_ff;
	iplan_f_ff = fftw_mpi_plan_dft_3d(fftsizex, fftsizey, fftsizez, f_ff, ff, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE);	//For inverse FFT
	//----- ----- ----- ----- ----- -----
	
	for(int i=0;i<local_n0;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				ii=(i*Ny+j)*Nz+k;
				iimpi=((local_0_start+i)*Ny+j)*Nz+k;
				den[ii][0] = den_out[iimpi];
				den[ii][1] = 0.0;
				//modulate the density field with given noise term
			}
		}
	}
	
	//evolve (evolve microstructure)
	for(int istep=0;istep<=nstep;istep++){
		
		//tempr = tempr + tempr0/istep;
		
		//loding
		if(fmod(istep,ndefor)==0){
			dx=dx+5.0e-4*dx0;
			dy=dy-5.0e-4/2.0*dx0;
			dz=dz-5.0e-4/2.0*dx0;
			
			stran=stran+5.0e-4*dx0;
			
			prepare_fft_3d(Nx,Ny,Nz,dx,dy,dz,kx,ky,kz,k2,k4); //get FFT coefficients
		}
		
		//take current density filed from real space to Fourier space (forward FFT transformation)
		//f_den = fftn(den);
		fftw_execute(plan_den);
		
		//calculate the value of denominator in Eq.7.10 at every grid points
		for(int i=0;i<local_n0;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ii=(i*Ny+j)*Nz+k;
					iimpi=((local_0_start+i)*Ny+j)*Nz+k;
					Linx[iimpi]=-k2[iimpi]*(tempr+1.0-2.0*k2[iimpi]+k4[iimpi]);
					denom[iimpi]=1.0-dtime*Linx[iimpi];
				}
			}
		}
		
		//calculate the nonlinear term, phi^3, in Eq.7.10
		for(int i=0;i<local_n0;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ii=(i*Ny+j)*Nz+k;
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
		for(int i=0;i<local_n0;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ii=(i*Ny+j)*Nz+k;
					iimpi=((local_0_start+i)*Ny+j)*Nz+k;
					Nonx[ii][0] = -k2[iimpi]*f_den3[ii][0];
					Nonx[ii][1] = -k2[iimpi]*f_den3[ii][1];
					f_den[ii][0] = (f_den[ii][0]+dtime*Nonx[ii][0])/denom[iimpi];
					f_den[ii][1] = (f_den[ii][1]+dtime*Nonx[ii][1])/denom[iimpi];
				}
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
					for(int k=0;k<Nz;k++){
						ii=(i*Ny+j)*Nz+k;
						iimpi=((local_0_start+i)*Ny+j)*Nz+k;
						den_out[iimpi]=den[ii][0]/(fftsizex*fftsizey*fftsizez);
						ss2[iimpi]=den_out[iimpi]*den_out[iimpi];
						ss4[iimpi]=ss2[iimpi]*ss2[iimpi];
					}
				}
			}
			
			for(int i=0;i<local_n0;i++){
				for(int j=0;j<Ny;j++){
					for(int k=0;k<Nz;k++){
						ii=(i*Ny+j)*Nz+k;
						iimpi=((local_0_start+i)*Ny+j)*Nz+k;
						f_ff[ii][0] = 0.5*f_den[ii][0]*(1.0-2.0*k2[iimpi]+k4[iimpi]);
						f_ff[ii][1] = 0.5*f_den[ii][1]*(1.0-2.0*k2[iimpi]+k4[iimpi]);
					}
				}
			}
			
			//ff=real(ifftn(f_ff));
			fftw_execute(iplan_f_ff);
			
			for(int i=0;i<local_n0;i++){
				for(int j=0;j<Ny;j++){
					for(int k=0;k<Nz;k++){
						ii=(i*Ny+j)*Nz+k;
						iimpi=((local_0_start+i)*Ny+j)*Nz+k;
						ff[ii][0] = (ff[ii][0]/(fftsizex*fftsizey*fftsizez))
								  *(den[ii][0]/(fftsizex*fftsizey*fftsizez))
						+ 0.5*tempr*ss2[iimpi]
						+ 0.25*ss4[iimpi];
						ff_out[ii]=ff[ii][0];
					}
				}
			}
			
			//integrate the free energy field
			energy = 0.0;
			
			for(int i=0;i<local_n0;i++){
				for(int j=0;j<Ny;j++){
					for(int k=0;k<Nz;k++){
						ii=(i*Ny+j)*Nz+k;
						energy = energy + ff[ii][0];
					}
				}
			}
			
			//average free energy density
			energy = energy/(Nx*Ny*Nz);
			
			if(istep==0){
				energy0=energy;
			}
			
			energy=energy-energy0;
			
			MPI_Gather(&den_out[local_0_start*Ny*Nz], local_n0*Ny*Nz, MPI_DOUBLE, &den_out[local_0_start*Ny*Nz], local_n0*Ny*Nz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			//
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			if (rank == 0){
				//print the average free energy density value to file
				fprintf(out2, "%14.6e %14.6e \n",stran, energy);
				
				//output the results in vtk file format for contour plots to be viewed by using paraview
				write_vtk_grid_values_3D(Nx,Ny,Nz,dx,dy,dz,dx0,istep,den_out,ff_out);
				
				printf("done step: %8d \n", istep);
			}
		}
		
		//if intermediate configuration files are required, print the density field to file
		if(istep==nstep){
			MPI_Gather(&den_out[local_0_start*Ny*Nz], local_n0*Ny*Nz, MPI_DOUBLE, &den_out[local_0_start*Ny*Nz], local_n0*Ny*Nz, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			//
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			if (rank == 0){
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						for(int k=0;k<Nz;k++){
							ii=(i*Ny+j)*Nz+k;
							fprintf(out1,"%5d %5d %5d %14.6e \n",i,j,k,den_out[ii]);
						}
					}
				}
			}
		}
		
		//for recycle
		for(int i=0;i<local_n0;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ii=(i*Ny+j)*Nz+k;
					den[ii][0] = den[ii][0]/(fftsizex*fftsizey*fftsizez);
					den[ii][1] = den[ii][1]/(fftsizex*fftsizey*fftsizez);
				}
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
