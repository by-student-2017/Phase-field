/* This function evaluates the derivative of elastic energy with
   respect to concentration. First, stress and strain values are
   solved with the iterative algorithm described earlier, 
   then derivative of elastic energy is evaluated for all grid points. */

#include <stdlib.h> //rand() and malloc
#include <math.h>
#include <fftw3.h>

/* Variable and array list
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  Nz: Number of grid points in the z-direction
  cm11: C11 component of elasticity matrix for matrix material
  cm12: C12 component of elasticity matrix for matrix material
  cm44: C44 component of elasticity matrix for matrix material
  cp11: C11 component of elasticity matrix for second phase
  cp12: C12 component of elasticity matrix for second phase
  cp44: C44 component of elasticity matrix for second phase
  ed11[Nx][Ny][Nz] to ed13[Nx][Ny][Nz]: Strain component of lattice defects
  ea[6]: Applied strains
  con[Nx][Ny][Nz]: Concentration
  s11[Nx][Ny][Nz] to s13[Nx][Ny][Nz]: Component of stress
  e11[Nx][Ny][Nz] to e13[Nx][Ny][Nz]: Component of strain
  delsdc[Nx][Ny][Nz]: Functional derivative of elastic energy
  //
  omeg11[Nx][Ny][Nz]: Coefficient needed for the Green's tensor
  omeg22[Nx][Ny][Nz]: Coefficient needed for the Green's tensor
  omeg33[Nx][Ny][Nz]: Coefficient needed for the Green's tensor
  omeg12[Nx][Ny][Nz]: Coefficient needed for the Green's tensor
  omeg23[Nx][Ny][Nz]: Coefficient needed for the Green's tensor
  omeg13[Nx][Ny][Nz]: Coefficient needed for the Green's tensor
  tmatx[3][3][3][3]: Green's tensor at i,j,k grid point (real part only)
*/

void green_tensor2_3D(int Nx, int Ny, int Nz,
	float *kx, float *ky, float *kz,
	float *omeg11, float *omeg22, float *omeg33,
	float *omeg12, float *omeg23, float *omeg13,
	int i, int j, int k,
	float *tmatx);

void solve_elasticity_3d(int Nx, int Ny, int Nz,
	float *kx, float *ky, float *kz,
	float *omeg11, float *omeg22, float *omeg33,
	float *omeg12, float *omeg23, float *omeg13,
	fftw_complex *s11, fftw_complex *s22, fftw_complex *s33,
	fftw_complex *s12, fftw_complex *s23, fftw_complex *s13,
	fftw_complex *e11, fftw_complex *e22, fftw_complex *e33,
	fftw_complex *e12, fftw_complex *e23, fftw_complex *e13,
	float *ed11, float *ed22, float *ed33,
	float *ed12, float *ed23, float *ed13,
	float cm11, float cm12, float cm44,
	float cp11, float cp12, float cp44,
	float *ea, float ei0,
	float *con, float *delsdc){
	
	//----- ----- ----- -----
	const int fftsizex = Nx, fftsizey = Ny, fftsizez = Nz;
	const int scale=fftsizex*fftsizey*fftsizez;
	float fftw3d_scale = (float)scale;
	
	//stress (head s series) components
	//----- ----- ----- -----
	fftw_complex *s11k, *s22k, *s33k;
	 s11k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 s22k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 s33k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	//
	fftw_complex *s12k, *s23k, *s13k;
	 s12k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 s23k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 s13k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	//----- ----- ----- -----
	fftw_plan plan_s11, plan_s22, plan_s33;
	 plan_s11  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, s11, s11k, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	 plan_s22  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, s22, s22k, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	 plan_s33  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, s33, s33k, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//
	fftw_plan plan_s12, plan_s23, plan_s13;
	 plan_s12  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, s12, s12k, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	 plan_s23  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, s23, s23k, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	 plan_s13  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, s13, s13k, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//----- ----- ----- -----
	
	//strain (head e series) components
	//----- ----- ----- -----
	fftw_complex *e11k, *e22k, *e33k;
	 e11k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 e22k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 e33k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	//
	fftw_complex *e12k, *e23k, *e13k;
	 e12k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 e23k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	 e13k = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * scale);
	//----- ----- ----- -----
	fftw_plan plan_e11, iplan_e11k;
	 plan_e11  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, e11, e11k, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_e11k = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, e11k,e11,  FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//----- ----- ----- -----
	fftw_plan plan_e22, iplan_e22k;
	 plan_e22  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, e22, e22k, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_e22k = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, e22k,e22,  FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//----- ----- ----- -----
	fftw_plan plan_e33, iplan_e33k;
	 plan_e33  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, e33, e33k, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_e33k = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, e33k,e33,  FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//----- ----- ----- -----
	//----- ----- ----- -----
	fftw_plan plan_e12, iplan_e12k;
	 plan_e12  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, e12, e12k, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_e12k = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, e12k,e12,  FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//----- ----- ----- -----
	fftw_plan plan_e23, iplan_e23k;
	 plan_e23  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, e23, e23k, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_e23k = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, e23k,e23,  FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//----- ----- ----- -----
	fftw_plan plan_e13, iplan_e13k;
	 plan_e13  = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, e13, e13k, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_e13k = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, e13k,e13,  FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//----- ----- ----- -----
	
	int NxNyNz=Nx*Ny*Nz;
	
	//----- ----- ----- ----- ----- ----- -----
	// eigenstrains (head ei series) components
	//float ei11[Nx][Ny][Nz];
	float *ei11 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//float ei22[Nx][Ny][Nz];
	float *ei22 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//float ei33[Nx][Ny][Nz];
	float *ei33 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//
	//float ei12[Nx][Ny][Nz];
	float *ei12 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//float ei23[Nx][Ny][Nz];
	float *ei23 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//float ei13[Nx][Ny][Nz];
	float *ei13 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//----- ----- ----- ----- ----- ----- -----
	
	//----- ----- ----- ----- ----- ----- -----
	// elastic modulus components
	//float c11[Nx][Ny][Nz];
	float  *c11 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//float c12[Nx][Ny][Nz];
	float  *c12 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//float c44[Nx][Ny][Nz];
	float  *c44 = (float *)malloc(sizeof(float)*( NxNyNz ));
	//----- ----- ----- ----- ----- ----- -----
	
	//----- ----- -----
	int ijk;
	int klij; //For tmatx
	//----- ----- -----
	
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				ijk=(i*Ny+j)*Nz+k;
				//----- ----- ----- ----- -----
				// Calculate the eigenstrains (head ei series)
				ei11[ijk] = ei0*con[ijk];
				ei22[ijk] = ei0*con[ijk];
				ei33[ijk] = ei0*con[ijk];
				//
				ei12[ijk] = 0.0*con[ijk];
				ei23[ijk] = 0.0*con[ijk];
				ei13[ijk] = 0.0*con[ijk];
				//----- ----- ----- ----- -----
				
				/* Calculate the effective elastic constants at 
				   the grid points based on the composition and
				   using Vegard's law */
				//----- ----- ----- ----- -----
				c11[ijk] = con[ijk]*cp11 + (1.0-con[ijk])*cm11;
				c12[ijk] = con[ijk]*cp12 + (1.0-con[ijk])*cm12;
				c44[ijk] = con[ijk]*cp44 + (1.0-con[ijk])*cm44;
				//----- ----- ----- ----- -----
			}
		}
	}
	
	/* Note: elastic modulus in this case. */
	//----- ----- ----- -----
	//float c22=c33=c11;
	//float c21=c12;
	//float c31=c13=c12;
	//float c32=c23=c12;
	//float c55=c66=c44;
	//----- ----- ----- -----
	float cm22,cm33;
		cm22=cm33=cm11;
	float cm21;
		cm21=cm12;
	float cm31,cm13;
		cm31=cm13=cm12;
	float cm32,cm23;
		cm32=cm23=cm12;
	float cm55, cm66;
		cm55=cm66=cm44;
	//----- ----- ----- -----
	float cp22, cp33;
		cp22=cp33=cp11;
	float cp21;
		cp21=cp12;
	float cp31, cp13;
		cp31=cp13=cp12;
	float cp32, cp23;
		cp32=cp23=cp12;
	float cp55, cp66;
		cp55=cp66=cp44;
	//----- ----- ----- -----
	float et21;
	float et32;
	float et31;
	//----- ----- ----- -----
	
	//----- ----- ----- ----- ----- ----- -----
	//float smatx_real[Nx][Ny][Nz][3][3];
	float *smatx_real = (float *)malloc(sizeof(float)*( NxNyNz*3*3 ));
	//
	//float smatx_imag[Nx][Ny][Nz][3][3];
	float *smatx_imag = (float *)malloc(sizeof(float)*( NxNyNz*3*3 ));
	//----- ----- ----- ----- ----- ----- -----
	//float ematx_real[Nx][Ny][Nz][3][3];
	float *ematx_real = (float *)malloc(sizeof(float)*( NxNyNz*3*3 ));
	//
	//float ematx_imag[Nx][Ny][Nz][3][3];
	float *ematx_imag = (float *)malloc(sizeof(float)*( NxNyNz*3*3 ));
	//----- ----- ----- ----- ----- ----- -----
	
	//float tmatx[3][3][3][3];
	float *tmatx = (float *)malloc(sizeof(float)*( 3*3*3*3 )); //real part only
	
	//----- ----- -----
	float et11=0.0;
	float et22=0.0;
	float et33=0.0;
	//
	float et12=0.0;
	float et23=0.0;
	float et13=0.0;
	//----- ----- -----
	
	//float sum_stress[Nx][Ny][Nz];
	float *sum_stress = (float *)malloc(sizeof(float)*( NxNyNz ));
	
	//----- ----- -----
	//Maximum number of iteration steps
	int niter=10;
	//----- ----- -----
	float old_norm=0.0;
	float normF=0.0;
	//----- ----- -----
	float conver=0.0;
	//----- ----- -----
	//Tolerance value of convergence tests
	float tolerance=0.001;
	//----- ----- -----
	
	/* Solve stress and strain field with 
	   iterative algorithm given in the text */
	for(int iter=0;iter<niter;iter++){
		
		/* Take stress and strain components from real space to
		   Fourier space (forward FFT). Step-a */
		//----- ----- ----- -----
		// stress (head s series)
		//s11k=fft3(s11);
		fftw_execute(plan_s11);
		//s22k=fft3(s22);
		fftw_execute(plan_s22);
		//s33k=fft3(s33);
		fftw_execute(plan_s33);
		//
		//s12k=fft3(s12);
		fftw_execute(plan_s12);
		//s23k=fft3(s23);
		fftw_execute(plan_s23);
		//s13k=fft3(s13);
		fftw_execute(plan_s13);
		//----- ----- ----- -----
		// strain (head e series)
		//e11k=fft3(e11);
		fftw_execute(plan_e11);
		//e22k=fft3(e22);
		fftw_execute(plan_e22);
		//e33k=fft3(e33);
		fftw_execute(plan_e33);
		//
		//e12k=fft3(e12);
		fftw_execute(plan_e12);
		//e23k=fft3(e23);
		fftw_execute(plan_e23);
		//e13k=fft3(e13);
		fftw_execute(plan_e13);
		//----- ----- ----- -----
		
		/* Form stress and strain tensors to be used in 
		   Eq.5.46, Step-b */
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					//----- ----- ----- ----- ----- -----
					// stress (smatx, sXY and sXYk)
					smatx_real[(ijk*3+0)*3+0]=s11k[ijk][0];
					smatx_real[(ijk*3+0)*3+1]=s12k[ijk][0];
					smatx_real[(ijk*3+0)*3+2]=s13k[ijk][0];
					//
					smatx_real[(ijk*3+1)*3+0]=s12k[ijk][0];
					smatx_real[(ijk*3+1)*3+1]=s22k[ijk][0];
					smatx_real[(ijk*3+1)*3+2]=s23k[ijk][0];
					//
					smatx_real[(ijk*3+2)*3+0]=s13k[ijk][0];
					smatx_real[(ijk*3+2)*3+1]=s23k[ijk][0];
					smatx_real[(ijk*3+2)*3+2]=s33k[ijk][0];
					//----- ----- ----- -----
					smatx_imag[(ijk*3+0)*3+0]=s11k[ijk][1];
					smatx_imag[(ijk*3+0)*3+1]=s12k[ijk][1];
					smatx_imag[(ijk*3+0)*3+2]=s13k[ijk][1];
					//
					smatx_imag[(ijk*3+1)*3+0]=s12k[ijk][1];
					smatx_imag[(ijk*3+1)*3+1]=s22k[ijk][1];
					smatx_imag[(ijk*3+1)*3+2]=s23k[ijk][1];
					//
					smatx_imag[(ijk*3+2)*3+0]=s13k[ijk][1];
					smatx_imag[(ijk*3+2)*3+1]=s23k[ijk][1];
					smatx_imag[(ijk*3+2)*3+2]=s33k[ijk][1];
					//----- ----- ----- ----- ----- -----
					
					//----- ----- ----- ----- ----- -----
					// strain (ematx, eXY and eXYk)
					ematx_real[(ijk*3+0)*3+0]=e11k[ijk][0];
					ematx_real[(ijk*3+0)*3+1]=e12k[ijk][0];
					ematx_real[(ijk*3+0)*3+2]=e13k[ijk][0];
					//
					ematx_real[(ijk*3+1)*3+0]=e12k[ijk][0];
					ematx_real[(ijk*3+1)*3+1]=e22k[ijk][0];
					ematx_real[(ijk*3+1)*3+2]=e23k[ijk][0];
					//
					ematx_real[(ijk*3+2)*3+0]=e13k[ijk][0];
					ematx_real[(ijk*3+2)*3+1]=e23k[ijk][0];
					ematx_real[(ijk*3+2)*3+2]=e33k[ijk][0];
					//----- ----- ----- -----
					ematx_imag[(ijk*3+0)*3+0]=e11k[ijk][1];
					ematx_imag[(ijk*3+0)*3+1]=e12k[ijk][1];
					ematx_imag[(ijk*3+0)*3+2]=e13k[ijk][1];
					//
					ematx_imag[(ijk*3+1)*3+0]=e12k[ijk][1];
					ematx_imag[(ijk*3+1)*3+1]=e22k[ijk][1];
					ematx_imag[(ijk*3+1)*3+2]=e23k[ijk][1];
					//
					ematx_imag[(ijk*3+2)*3+0]=e13k[ijk][1];
					ematx_imag[(ijk*3+2)*3+1]=e23k[ijk][1];
					ematx_imag[(ijk*3+2)*3+2]=e33k[ijk][1];
					//----- ----- ----- ----- ----- -----
				}
			}
		}
		
		//Green operator
		// Calculate strain tensor, Eq.5.46, Step-b
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					//----- ----- ----- ----- ----- -----
					green_tensor2_3D(Nx,Ny,Nz,
									kx,ky,kz,
									omeg11,omeg22,omeg33,
									omeg12,omeg23,omeg13,
									i,j,k,
									tmatx);
					//----- ----- ----- ----- ----- -----
					for(int kk=0;kk<3;kk++){
						for(int ll=0;ll<3;ll++){
							for(int ii=0;ii<3;ii++){
								for(int jj=0;jj<3;jj++){
									klij=((kk*3+ll)*3+ii)*3+jj;
									/* Eq.5.46(b): new epsilon(zeta) = epsilon(zeta) - sum( gamma(zeta)*sigma(zeta) )
									   where gamma=tmatx, sigma=smatx
									   Note: tmatx is real part only */
									//----- ----- ----- ----- ----- -----
									ematx_real[(ijk*3+ii)*3+jj] -= tmatx[klij]*smatx_real[(ijk*3+kk)*3+ll];
									//
									ematx_imag[(ijk*3+ii)*3+jj] -= tmatx[klij]*smatx_imag[(ijk*3+kk)*3+ll];
									//----- ----- ----- ----- ----- -----
								}//jj
							}//ii
						}//ll
					}//kk
					//----- ----- ----- ----- ----- -----
				}//Nz
			}//Ny
		}//Nx
		
		// Rearrange strain components using symmetry of strain tensor
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					//----- ----- ----- ----- ----- -----
					// strain (ematx, eXY and eXYk)
					e11k[ijk][0]=ematx_real[(ijk*3+0)*3+0];
					e12k[ijk][0]=ematx_real[(ijk*3+0)*3+1];
					e13k[ijk][0]=ematx_real[(ijk*3+0)*3+2];
					//
					e12k[ijk][0]=ematx_real[(ijk*3+1)*3+0];
					e22k[ijk][0]=ematx_real[(ijk*3+1)*3+1];
					e23k[ijk][0]=ematx_real[(ijk*3+1)*3+2];
					//
					e13k[ijk][0]=ematx_real[(ijk*3+2)*3+0];
					e23k[ijk][0]=ematx_real[(ijk*3+2)*3+1];
					e33k[ijk][0]=ematx_real[(ijk*3+2)*3+2];
					//----- ----- ----- -----
					e11k[ijk][1]=ematx_imag[(ijk*3+0)*3+0];
					e12k[ijk][1]=ematx_imag[(ijk*3+0)*3+1];
					e13k[ijk][1]=ematx_imag[(ijk*3+0)*3+2];
					//
					e12k[ijk][1]=ematx_imag[(ijk*3+1)*3+0];
					e22k[ijk][1]=ematx_imag[(ijk*3+1)*3+1];
					e23k[ijk][1]=ematx_imag[(ijk*3+1)*3+2];
					//
					e13k[ijk][1]=ematx_imag[(ijk*3+2)*3+0];
					e23k[ijk][1]=ematx_imag[(ijk*3+2)*3+1];
					e33k[ijk][1]=ematx_imag[(ijk*3+2)*3+2];
					//----- ----- ----- ----- ----- -----
				}
			}
		}
		
		//From Fourier space to real space
		/* Take strain components from Fourier space back to
		   real space (inverse FFT), Step-c */
		//----- ----- ----- -----
		//e11=real(ifft3(e11k));
		fftw_execute(iplan_e11k);
		//e22=real(ifft3(e22k));
		fftw_execute(iplan_e22k);
		//e33=real(ifft3(e33k));
		fftw_execute(iplan_e33k);
		//
		//e12=real(ifft3(e12k));
		fftw_execute(iplan_e12k);
		//e23=real(ifft3(e23k));
		fftw_execute(iplan_e23k);
		//e13=real(ifft3(e13k));
		fftw_execute(iplan_e13k);
		//----- ----- ----- -----
		
		//Calculate stresses
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					//----- ----- ----- ----- ----- -----
					// strain (ematx, eXY and eXYk)
					//----- ----- ----- -----
					e11[ijk][0] /= fftw3d_scale;
					e22[ijk][0] /= fftw3d_scale;
					e33[ijk][0] /= fftw3d_scale;
					//
					e12[ijk][0] /= fftw3d_scale;
					e23[ijk][0] /= fftw3d_scale;
					e13[ijk][0] /= fftw3d_scale;
					//----- ----- ----- -----
					e11[ijk][1] /= fftw3d_scale;
					e22[ijk][1] /= fftw3d_scale;
					e33[ijk][1] /= fftw3d_scale;
					//
					e12[ijk][1] /= fftw3d_scale;
					e23[ijk][1] /= fftw3d_scale;
					e13[ijk][1] /= fftw3d_scale;
					//----- ----- ----- ----- ----- -----
					
					//----- ----- ----- ----- ----- ----- ----- -----
					/*s11[ijk][0]=c11[ijk]*(ea[0]+e11[ijk][0]-ei11[ijk]-ed11[ijk])
							   +c12[ijk]*(ea[1]+e22[ijk][0]-ei22[ijk]-ed22[ijk])
							   +c13[ijk]*(ea[2]+e33[ijk][0]-ei33[ijk]-ed33[ijk]);
					s22[ijk][0]=c21[ijk]*(ea[0]+e11[ijk][0]-ei11[ijk]-ed11[ijk])
							   +c22[ijk]*(ea[1]+e22[ijk][0]-ei22[ijk]-ed22[ijk])
							   +c23[ijk]*(ea[2]+e33[ijk][0]-ei33[ijk]-ed33[ijk]);
					s22[ijk][0]=c31[ijk]*(ea[0]+e11[ijk][0]-ei11[ijk]-ed11[ijk])
							   +c32[ijk]*(ea[1]+e22[ijk][0]-ei22[ijk]-ed22[ijk])
							   +c33[ijk]*(ea[2]+e33[ijk][0]-ei33[ijk]-ed33[ijk]);*/
					//
					// c13[ijk]=c12[ijk], c23[ijk]=c13[ijk], c22[ijk]=c11[ijk], c33[ijk]=c11[ijk], etc
					/*s11[ijk][0]=c11[ijk]*(ea[0]+e11[ijk][0]-ei11[ijk]-ed11[ijk])
							   +c12[ijk]*(ea[1]+e22[ijk][0]-ei22[ijk]-ed22[ijk])
							   +c12[ijk]*(ea[2]+e33[ijk][0]-ei33[ijk]-ed33[ijk]);
					s22[ijk][0]=c12[ijk]*(ea[0]+e11[ijk][0]-ei11[ijk]-ed11[ijk])
							   +c11[ijk]*(ea[1]+e22[ijk][0]-ei22[ijk]-ed22[ijk])
							   +c12[ijk]*(ea[2]+e33[ijk][0]-ei33[ijk]-ed33[ijk]);
					s33[ijk][0]=c12[ijk]*(ea[0]+e11[ijk][0]-ei11[ijk]-ed11[ijk])
							   +c12[ijk]*(ea[1]+e22[ijk][0]-ei22[ijk]-ed22[ijk])
							   +c11[ijk]*(ea[2]+e33[ijk][0]-ei33[ijk]-ed33[ijk]);*/
					//
					// Calculate strain (head e series) components
					et11=( ea[0]+e11[ijk][0]-ei11[ijk]-ed11[ijk] );
					et22=( ea[1]+e22[ijk][0]-ei22[ijk]-ed22[ijk] );
					et33=( ea[2]+e33[ijk][0]-ei33[ijk]-ed33[ijk] );
					//
					// Calculate stress (head s series) components
					s11[ijk][0]=( c11[ijk]*et11 + c12[ijk]*et22 + c12[ijk]*et33 );
					s22[ijk][0]=( c12[ijk]*et11 + c11[ijk]*et22 + c12[ijk]*et33 );
					s33[ijk][0]=( c12[ijk]*et11 + c12[ijk]*et22 + c11[ijk]*et33 );
					//
					s11[ijk][1]=0.0;
					s22[ijk][1]=0.0;
					s33[ijk][1]=0.0;
					//----- ----- ----- ----- ----- ----- ----- -----
					
					//----- ----- ----- ----- ----- ----- ----- -----
					/*s12[ijk][0]=c44[ijk]*(ea[3]+e12[ijk][0]-ei12[ijk]-ed12[ijk])
							   +c44[ijk]*(ea[3]+e21[ijk][0]-ei21[ijk]-ed21[ijk]);
					s23[ijk][0]=c55[ijk]*(ea[4]+e23[ijk][0]-ei23[ijk]-ed23[ijk])
							   +c55[ijk]*(ea[4]+e32[ijk][0]-ei32[ijk]-ed32[ijk]);
					s13[ijk][0]=c66[ijk]*(ea[5]+e13[ijk][0]-ei13[ijk]-ed13[ijk])
							   +c66[ijk]*(ea[5]+e31[ijk][0]-ei31[ijk]-ed31[ijk]);*/
					//
					// c55[ijk]=c44[ijk], c66[ijk]=c44[ijk], etc
					/*s12[ijk][0]=c44[ijk]*(ea[3]+e12[ijk][0]-ei12[ijk]-ed12[ijk])*2.0;
					  s23[ijk][0]=c44[ijk]*(ea[4]+e23[ijk][0]-ei23[ijk]-ed23[ijk])*2.0;
					  s13[ijk][0]=c44[ijk]*(ea[5]+e13[ijk][0]-ei13[ijk]-ed13[ijk])*2.0;*/
					//
					// Calculate strain (head e series) components
					et12=( ea[3]+e12[ijk][0]-ei12[ijk]-ed12[ijk] );
					et23=( ea[4]+e23[ijk][0]-ei23[ijk]-ed23[ijk] );
					et13=( ea[5]+e13[ijk][0]-ei13[ijk]-ed13[ijk] );
					//
					// Calculate stress (head s series) components
					s12[ijk][0]=c44[ijk]*et12*2.0;
					s23[ijk][0]=c44[ijk]*et23*2.0;
					s13[ijk][0]=c44[ijk]*et13*2.0;
					//
					s12[ijk][1]=0.0;
					s23[ijk][1]=0.0;
					s13[ijk][1]=0.0;
					//----- ----- ----- ----- ----- ----- ----- -----
				}
			}
		}
		
		//check convergence
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					sum_stress[ijk] = ( s11[ijk][0] + s22[ijk][0] + s33[ijk][0]
									   +s12[ijk][0] + s23[ijk][0] + s13[ijk][0] );
				}
			}
		}
		
		//normF=norm(sum_stress,2.0);
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				for(int k=0;k<Nz;k++){
					ijk=(i*Ny+j)*Nz+k;
					normF = normF + sum_stress[ijk]*sum_stress[ijk];
				}
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
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				ijk=(i*Ny+j)*Nz+k;
				//----- ----- ----- ----- ----- ----- ----- -----
				// Calculate strain (head e series) components
				et11=( ea[0]+e11[ijk][0]-ei11[ijk]-ed11[ijk] );
				et22=( ea[1]+e22[ijk][0]-ei22[ijk]-ed22[ijk] );
				et33=( ea[2]+e33[ijk][0]-ei33[ijk]-ed33[ijk] );
				//
				et12=( ea[3]+e12[ijk][0]-ei12[ijk]-ed12[ijk] );
				et23=( ea[4]+e23[ijk][0]-ei23[ijk]-ed23[ijk] );
				et13=( ea[5]+e13[ijk][0]-ei13[ijk]-ed13[ijk] );
				//----- ----- ----- ----- ----- ----- ----- -----
				
				//----- ----- ----- ----- ----- ----- ----- -----
				//Functional derivative of the elastic energy with respect to composition
				/* F=(1/2)*sigma[i][j]*(epsilon[i][j] - epsilon0[i][j])
				   sigma[i][j] = C[i][j][k][l]*(epsilon[k][l] - epsilon0[k][l])
				   epsilon0[i][j] is the position- and composition-dependent eigenstrains */
				/* dF/dc = dF/d(con) = (1/2)*( dCijkl/d(con)*etij*etkl + Cijkl*d(etij)/d(con)*etkl + Cijkl*etij*d(etkl)/d(con) )
				   dCijkl/d(con) = (cpijkl - cmijkl), d(etij)/d(con) = -d(eiij)/d(con) = -ei0 */
				/* delsdc[ijk][0]=0.5*( (cp11-cm11)*et11*et11 -c11[ijk]*ei0*et11 -c11[ijk]*et11*ei0
								    +(cp12-cm12)*et11*et22 -c12[ijk]*ei0*et11 -c12[ijk]*et22*ei0
								    +(cp13-cm13)*et11*et33 -c13[ijk]*ei0*et11 -c13[ijk]*et33*ei0
								    //
								    +(cp21-cm21)*et22*et11 -c21[ijk]*ei0*et22 -c21[ijk]*et11*ei0
								    +(cp22-cm22)*et22*et22 -c22[ijk]*ei0*et22 -c22[ijk]*et22*ei0
								    +(cp23-cm23)*et22*et33 -c23[ijk]*ei0*et22 -c23[ijk]*et33*ei0
								    //
								    +(cp31-cm31)*et33*et11 -c31[ijk]*ei0*et33 -c31[ijk]*et11*ei0
								    +(cp32-cm32)*et33*et22 -c32[ijk]*ei0*et33 -c32[ijk]*et22*ei0
								    +(cp33-cm33)*et33*et33 -c33[ijk]*ei0*et33 -c33[ijk]*et33*ei0
								    //
								    +(cp44-cm44)*et12*et12 -c44[ijk]*ei0*et12 -c44[ijk]*et12*ei0
								    +(cp44-cm44)*et21*et21 -c44[ijk]*ei0*et21 -c44[ijk]*et21*ei0
								    //
								    +(cp55-cm55)*et23*et23 -c55[ijk]*ei0*et23 -c55[ijk]*et23*ei0
								    +(cp55-cm55)*et32*et32 -c55[ijk]*ei0*et32 -c55[ijk]*et32*ei0
								    //
								    +(cp66-cm66)*et13*et13 -c66[ijk]*ei0*et13 -c66[ijk]*et13*ei0
								    +(cp66-cm66)*et31*et31 -c66[ijk]*ei0*et31 -c66[ijk]*et31*ei0
								   ); */
				//
				// c13[ijk]=c12[ijk], c23[ijk]=c13[ijk], c22[ijk]=c11[ijk], c33[ijk]=c11[ijk], etc
				et21=et12;
				et32=et23;
				et31=et13;
				//
				delsdc[ijk]=0.5*( (cp11-cm11)*et11*et11 -c11[ijk]*ei0*( et11 + et11 )
								    +(cp12-cm12)*et11*et22 -c12[ijk]*ei0*( et11 + et22 )
								    +(cp13-cm13)*et11*et33 -c12[ijk]*ei0*( et11 + et33 )
								    //
								    +(cp21-cm21)*et22*et11 -c12[ijk]*ei0*( et22 + et11 )
								    +(cp22-cm22)*et22*et22 -c11[ijk]*ei0*( et22 + et22 )
								    +(cp23-cm23)*et22*et33 -c12[ijk]*ei0*( et22 + et33 )
								    //
								    +(cp31-cm31)*et33*et11 -c12[ijk]*ei0*( et33 + et11 )
								    +(cp32-cm32)*et33*et22 -c12[ijk]*ei0*( et33 + et22 )
								    +(cp33-cm33)*et33*et33 -c11[ijk]*ei0*( et33 + et33 )
								    //
								    +(cp44-cm44)*et12*et12 -c44[ijk]*ei0*( et12 + et12 )
								    +(cp44-cm44)*et21*et21 -c44[ijk]*ei0*( et21 + et21 )
								    //
								    +(cp55-cm55)*et23*et23 -c44[ijk]*ei0*( et23 + et23 )
								    +(cp55-cm55)*et32*et32 -c44[ijk]*ei0*( et32 + et32 )
								    //
								    +(cp66-cm66)*et13*et13 -c44[ijk]*ei0*( et13 + et13 )
								    +(cp66-cm66)*et31*et31 -c44[ijk]*ei0*( et31 + et31 )
								   );
				//delsdc[ijk][1]=0.0;
				//----- ----- ----- ----- ----- ----- ----- -----
			}
		}
	}
	
	//----- ----- ----- ----- -----
	fftw_destroy_plan(plan_s11);
	fftw_destroy_plan(plan_s22);
	fftw_destroy_plan(plan_s33);
	//
	fftw_destroy_plan(plan_s12);
	fftw_destroy_plan(plan_s23);
	fftw_destroy_plan(plan_s13);
	//----- ----- ----- ----- -----
	 fftw_destroy_plan(plan_e11);
	fftw_destroy_plan(iplan_e11k);
	 fftw_destroy_plan(plan_e22);
	fftw_destroy_plan(iplan_e22k);
	 fftw_destroy_plan(plan_e33);
	fftw_destroy_plan(iplan_e33k);
	//
	 fftw_destroy_plan(plan_e12);
	fftw_destroy_plan(iplan_e12k);
	 fftw_destroy_plan(plan_e23);
	fftw_destroy_plan(iplan_e23k);
	 fftw_destroy_plan(plan_e13);
	fftw_destroy_plan(iplan_e13k);
	//----- ----- ----- ----- -----
	fftw_free(s11k);
	fftw_free(s22k);
	fftw_free(s33k);
	//
	fftw_free(s12k);
	fftw_free(s23k);
	fftw_free(s13k);
	//----- ----- ----- ----- -----
	fftw_free(e11k);
	fftw_free(e22k);
	fftw_free(e33k);
	//
	fftw_free(e12k);
	fftw_free(e23k);
	fftw_free(e13k);
	//----- ----- ----- ----- ----- -----
	free(smatx_real);
	free(smatx_imag);
	free(ematx_real);
	free(ematx_imag);
	//----- ----- ----- ----- -----
	free(sum_stress);
	//----- ----- ----- ----- -----
	free(ei11);
	free(ei22);
	free(ei33);
	//
	free(ei12);
	free(ei23);
	free(ei13);
	//----- ----- ----- ----- -----
	free(c11);
	free(c12);
	free(c44);
	//----- ----- ----- ----- -----
	
	return;
}