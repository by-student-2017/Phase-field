/* This function evaluates the derivative of elastic energy with
   respect to concentration. First, stress and strain values are
   solved with the iterative algorithm described earlier, 
   then derivative of elastic energy is evaluated for all grid points. */

#include <stdlib.h> //rand() and malloc
#include <math.h>
//#include <fftw3.h>
#include <mpi.h> //mpi version
#include <fftw3-mpi.h> //mpi version

/* Variable and array list
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  cm11: C11 component of elasticity matrix for matrix material
  cm12: C12 component of elasticity matrix for matrix material
  cm44: C44 component of elasticity matrix for matrix material
  cp11: C11 component of elasticity matrix for second phase
  cp12: C12 component of elasticity matrix for second phase
  cp44: C44 component of elasticity matrix for second phase
  ed11: Strain component of lattice defects
  ed22: Strain component of lattice defects
  ed12: Strain component of lattice defects
  ei0: Magnitude of eigenstrains
  ea[3]: Applied strains
  con[Nx][Ny]: Concentration
  s11[Nx][Ny]: Component of stress
  s22[Nx][Ny]: Component of stress
  s12[Nx][Ny]: Component of stress
  e11[Nx][Ny]: Component of strain
  e22[Nx][Ny]: Component of strain
  e12[Nx][Ny]: Component of strain
  delsdc[Nx][Ny]: Functional derivative of elastic energy
  tmatx[Nx][Ny][2][2][2][2]: Values of Green's tensor at all grid points (real part only)
*/

void solve_elasticity_2d_mpi(int Nx, int Ny,
	double *tmatx,
	fftw_complex *s11, fftw_complex *s22, fftw_complex *s12,
	fftw_complex *e11, fftw_complex *e22, fftw_complex *e12,
	double *ed11, double *ed22, double *ed12,
	double cm11, double cm12, double cm44,
	double cp11, double cp12, double cp44,
	double *ea,
	double ei0,
	fftw_complex *con, fftw_complex *delsdc){
	
	//----- ----- ----- -----
	//----- ----- ----- ----- ----- -----
	const ptrdiff_t fftsizex = Nx, fftsizey = Ny;
	ptrdiff_t alloc_local, local_n0, local_0_start;
	
	//int rank;
	//int num_proc;
	
	//MPI_Init(&argc, &argv);
	//MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
	//fftw_mpi_init();
	
	alloc_local = fftw_mpi_local_size_2d(fftsizex, fftsizey, MPI_COMM_WORLD, &local_n0, &local_0_start);
	
	fftw_complex *s11k, *s22k, *s12k;
	s11k = fftw_alloc_complex(alloc_local);
	s22k = fftw_alloc_complex(alloc_local);
	s12k = fftw_alloc_complex(alloc_local);
	
	fftw_plan plan_s11, plan_s22, plan_s12;
	 plan_s11  = fftw_mpi_plan_dft_2d(fftsizex, fftsizey, s11, s11k, MPI_COMM_WORLD, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	 plan_s22  = fftw_mpi_plan_dft_2d(fftsizex, fftsizey, s22, s22k, MPI_COMM_WORLD, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	 plan_s12  = fftw_mpi_plan_dft_2d(fftsizex, fftsizey, s12, s12k, MPI_COMM_WORLD, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//----- ----- ----- -----
	fftw_complex *e11k, *e22k, *e12k;
	e11k = fftw_alloc_complex(alloc_local);
	e22k = fftw_alloc_complex(alloc_local);
	e12k = fftw_alloc_complex(alloc_local);
	
	fftw_plan plan_e11, iplan_e11k;
	 plan_e11  = fftw_mpi_plan_dft_2d(fftsizex, fftsizey, e11, e11k, MPI_COMM_WORLD, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_e11k = fftw_mpi_plan_dft_2d(fftsizex, fftsizey, e11k, e11, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	fftw_plan plan_e22, iplan_e22k;
	 plan_e22  = fftw_mpi_plan_dft_2d(fftsizex, fftsizey, e22, e22k, MPI_COMM_WORLD, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_e22k = fftw_mpi_plan_dft_2d(fftsizex, fftsizey, e22k, e22, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	fftw_plan plan_e12, iplan_e12k;
	 plan_e12  = fftw_mpi_plan_dft_2d(fftsizex, fftsizey, e12, e12k, MPI_COMM_WORLD, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_e12k = fftw_mpi_plan_dft_2d(fftsizex, fftsizey, e12k, e12, MPI_COMM_WORLD, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//----- ----- ----- -----
	
	//double smatx_real[Nx][Ny][2][2];
	double *smatx_real = (double *)malloc(sizeof(double)*( Nx*Ny*2*2 ));
	//double ematx_real[Nx][Ny][2][2];
	double *ematx_real = (double *)malloc(sizeof(double)*( Nx*Ny*2*2 ));
	//
	//double smatx_imag[Nx][Ny][2][2];
	double *smatx_imag = (double *)malloc(sizeof(double)*( Nx*Ny*2*2 ));
	//double ematx_imag[Nx][Ny][2][2];
	double *ematx_imag = (double *)malloc(sizeof(double)*( Nx*Ny*2*2 ));
	
	//double sum_stress[Nx][Ny];
	double *sum_stress = (double *)malloc(sizeof(double)*( Nx*Ny ));
	
	double old_norm=0.0;
	double normF=0.0;
	
	double conver=0.0;
	
	double et11=0.0;
	double et22=0.0;
	double et12=0.0;
	
	int ii;
	int iimpi, ijmpi;
	
	//Maximum number of iteration steps
	int niter=10;
	
	//Tolerance value of convergence tests
	double tolerance=0.001;
	
	//double ei11[Nx][Ny];
	double *ei11 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double ei22[Nx][Ny];
	double *ei22 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double ei12[Nx][Ny];
	double *ei12 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	
	//double c11[Nx][Ny];
	double  *c11 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double c12[Nx][Ny];
	double  *c12 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double c44[Nx][Ny];
	double  *c44 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	
	for(int i=0;i<local_n0;i++){
		for(int j=0;j<Ny;j++){
			ii=i*Ny+j;
			iimpi=(local_0_start+i)*Ny+j;
			
			//Calculate the eigenstrains
			ei11[iimpi]=ei0*con[ii][0];
			ei22[iimpi]=ei0*con[ii][0];
			ei12[iimpi]=0.0*con[ii][0];
			
			/* Calculate the effective elastic constants at 
			   the grid points based on the composition and
			   using Vegard's law */
			c11[iimpi]=con[ii][0]*cp11+(1.0-con[ii][0])*cm11;
			c12[iimpi]=con[ii][0]*cp12+(1.0-con[ii][0])*cm12;
			c44[iimpi]=con[ii][0]*cp44+(1.0-con[ii][0])*cm44;
		}
	}
	
	/* Solve stress and strain field with 
	   iterative algorithm given in the text */
	for(int iter=0;iter<niter;iter++){
		
		/* Take stress and strain components from real space to
		   Fourier space (forward FFT). Step-a */
		// stress
		//s11k=fft2(s11);
		fftw_execute(plan_s11);
		//s22k=fft2(s22);
		fftw_execute(plan_s22);
		//s12k=fft2(s12);
		fftw_execute(plan_s12);
		//
		// strain
		//e11k=fft2(e11);
		fftw_execute(plan_e11);
		//e22k=fft2(e22);
		fftw_execute(plan_e22);
		//e12k=fft2(e12);
		fftw_execute(plan_e12);
		
		/* Form stress and strain tensors to be used in 
		   Eq.5.46, Step-b */
		for(int i=0;i<local_n0;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				iimpi=(local_0_start+i)*Ny+j;
				
				smatx_real[(iimpi*2+0)*2+0]=s11k[ii][0];
				smatx_real[(iimpi*2+0)*2+1]=s12k[ii][0];
				smatx_real[(iimpi*2+1)*2+0]=s12k[ii][0];
				smatx_real[(iimpi*2+1)*2+1]=s22k[ii][0];
				//
				smatx_imag[(iimpi*2+0)*2+0]=s11k[ii][1];
				smatx_imag[(iimpi*2+0)*2+1]=s12k[ii][1];
				smatx_imag[(iimpi*2+1)*2+0]=s12k[ii][1];
				smatx_imag[(iimpi*2+1)*2+1]=s22k[ii][1];
				//
				ematx_real[(iimpi*2+0)*2+0]=e11k[ii][0];
				ematx_real[(iimpi*2+0)*2+1]=e12k[ii][0];
				ematx_real[(iimpi*2+1)*2+0]=e12k[ii][0];
				ematx_real[(iimpi*2+1)*2+1]=e22k[ii][0];
				//
				ematx_imag[(iimpi*2+0)*2+0]=e11k[ii][1];
				ematx_imag[(iimpi*2+0)*2+1]=e12k[ii][1];
				ematx_imag[(iimpi*2+1)*2+0]=e12k[ii][1];
				ematx_imag[(iimpi*2+1)*2+1]=e22k[ii][1];
			}
		}
		
		//Green operator
		// Calculate strain tensor, Eq.5.46, Step-b
		for(int i=0;i<local_n0;i++){
			for(int j=0;j<Ny;j++){
				//ij=i*Ny+j;
				ijmpi=(local_0_start+i)*Ny+j;
				//
				for(int kk=0;kk<2;kk++){
					for(int ll=0;ll<2;ll++){
						for(int ii=0;ii<2;ii++){
							for(int jj=0;jj<2;jj++){
								/* Eq.5.46(b): new epsilon(zeta) = epsilon(zeta) - sum( gamma(zeta)*sigma(zeta) )
								   where gamma=tmatx, sigma=smatx
								   Note: tmatx is real part only */
								ematx_real[(ijmpi*2+ii)*2+jj]=ematx_real[(ijmpi*2+ii)*2+jj]
									-tmatx[(((ijmpi*2+kk)*2+ll)*2+ii)*2+jj]*smatx_real[(ijmpi*2+kk)*2+ll];
								//
								ematx_imag[(ijmpi*2+ii)*2+jj]=ematx_imag[(ijmpi*2+ii)*2+jj]
									-tmatx[(((ijmpi*2+kk)*2+ll)*2+ii)*2+jj]*smatx_imag[(ijmpi*2+kk)*2+ll];
							}//jj
						}//ii
					}//ll
				}//kk
				//
			}//Ny
		}//Nx
		
		// Rearrange strain components using symmetry of strain tensor
		for(int i=0;i<local_n0;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				iimpi=(local_0_start+i)*Ny+j;
				
				e11k[ii][0]=ematx_real[(iimpi*2+0)*2+0];
				e12k[ii][0]=ematx_real[(iimpi*2+0)*2+1];
				e12k[ii][0]=ematx_real[(iimpi*2+1)*2+0];
				e22k[ii][0]=ematx_real[(iimpi*2+1)*2+1];
				//
				e11k[ii][1]=ematx_imag[(iimpi*2+0)*2+0];
				e12k[ii][1]=ematx_imag[(iimpi*2+0)*2+1];
				e12k[ii][1]=ematx_imag[(iimpi*2+1)*2+0];
				e22k[ii][1]=ematx_imag[(iimpi*2+1)*2+1];
			}
		}
		
		//From Fourier space to real space
		/* Take strain components from Fourier space back to
		   real space (inverse FFT), Step-c */
		//e11=real(ifft2(e11k));
		fftw_execute(iplan_e11k);
		//e22=real(ifft2(e22k));
		fftw_execute(iplan_e22k);
		//e12=real(ifft2(e12k));
		fftw_execute(iplan_e12k);
		
		//Calculate stresses
		for(int i=0;i<local_n0;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				iimpi=(local_0_start+i)*Ny+j;
				//
				e11[ii][0]=e11[ii][0]/(fftsizex*fftsizey);
				e22[ii][0]=e22[ii][0]/(fftsizex*fftsizey);
				e12[ii][0]=e12[ii][0]/(fftsizex*fftsizey);
				//
				e11[ii][1]=e11[ii][1]/(fftsizex*fftsizey);
				e22[ii][1]=e22[ii][1]/(fftsizex*fftsizey);
				e12[ii][1]=e12[ii][1]/(fftsizex*fftsizey);
				//
				s11[ii][0]=c11[iimpi]*(ea[0]+e11[ii][0]-ei11[iimpi]-ed11[iimpi])
						  +c12[iimpi]*(ea[1]+e22[ii][0]-ei22[iimpi]-ed22[iimpi]);
				s11[ii][1]=0.0;
				s22[ii][0]=c11[iimpi]*(ea[1]+e22[ii][0]-ei22[iimpi]-ed22[iimpi])
						  +c12[iimpi]*(ea[0]+e11[ii][0]-ei11[iimpi]-ed11[iimpi]);
				s22[ii][1]=0.0;
				s12[ii][0]=2.0
						  *c44[iimpi]*(ea[2]+e12[ii][0]-ei12[iimpi]-ed12[iimpi]);
				s12[ii][1]=0.0;
			}
		}
		
		//check convergence
		for(int i=0;i<local_n0;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				iimpi=(local_0_start+i)*Ny+j;
				sum_stress[iimpi] = ( s11[ii][0] + s22[ii][0] + s12[ii][0] );
			}
		}
		
		//normF=norm(sum_stress,2.0);
		for(int i=0;i<local_n0;i++){
			for(int j=0;j<Ny;j++){
				//ii=i*Ny+j;
				iimpi=(local_0_start+i)*Ny+j;
				normF = normF + sum_stress[iimpi]*sum_stress[iimpi];
			}
		}
		normF=sqrt(normF);
		
		if(iter==1){
			conver=fabs((normF-old_norm)/(old_norm));
			if(conver<=tolerance){
				break;
			}
		}
		old_norm=normF;
		
	}//end iter
	
	//strain energy
	//Calculate functional derivative of elastic energy
	// sum strain components
	for(int i=0;i<local_n0;i++){
		for(int j=0;j<Ny;j++){
			ii=i*Ny+j;
			iimpi=(local_0_start+i)*Ny+j;
			
			//Calculate strain components
			et11=ea[0]+e11[ii][0]-ei11[iimpi]-ed11[iimpi];
			et22=ea[1]+e22[ii][0]-ei22[iimpi]-ed22[iimpi];
			et12=ea[2]+e12[ii][0]-ei12[iimpi]-ed12[iimpi];
			
			//Functional derivative of the elastic energy with respect to composition
			delsdc[ii][0]=0.5*(et11*( (cp12-cm12)*et22 + (cp11-cm11)*et11 - c12[iimpi]*ei0 - c11[iimpi]*ei0 )
							   -ei0*(  c12[iimpi]*et22 +  c11[iimpi]*et11 )
						      +et22*( (cp11-cm11)*et22 + (cp12-cm12)*et11 - c12[iimpi]*ei0 - c11[iimpi]*ei0 )
							   -ei0*(  c11[iimpi]*et22 +  c12[iimpi]*et11 )
							   +2.0*(cp44-cm44)*et12*et12
							   -4.0*ei0*c44[iimpi]*et12
							  );
			delsdc[ii][1]=0.0;
		}
	}
	
	//----- ----- ----- ----- ----- -----
	fftw_destroy_plan(plan_s11);
	fftw_destroy_plan(plan_s22);
	fftw_destroy_plan(plan_s12);
	//
	fftw_destroy_plan(plan_e11);
	fftw_destroy_plan(iplan_e11k);
	fftw_destroy_plan(plan_e22);
	fftw_destroy_plan(iplan_e22k);
	fftw_destroy_plan(plan_e12);
	fftw_destroy_plan(iplan_e12k);
	//----- ----- ----- ----- ----- -----
	fftw_free(s11k);
	fftw_free(s22k);
	fftw_free(s12k);
	//
	fftw_free(e11k);
	fftw_free(e22k);
	fftw_free(e12k);
	//----- ----- ----- ----- ----- -----
	free(smatx_real);
	free(smatx_imag);
	free(ematx_real);
	free(ematx_imag);
	//
	free(sum_stress);
	//
	free(ei11);
	free(ei22);
	free(ei12);
	//
	free(c11);
	free(c12);
	free(c44);
	//----- ----- ----- ----- ----- -----
	
	return;
}