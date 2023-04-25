#include <stdlib.h> //rand() and malloc
#include <mpi.h> //mpi version
#include <fftw3-mpi.h> //mpi version

/* In here, it is given for 2D case. The 3D version which
   takes into account the six components of the strains and
   stresses is given in Appendix-B for those who wish to expand
   their analysis into the 3D. */

/* Variable and array list
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  cm11: C11 component of elasticity matrix for matrix material
  cm12: C12 component of elasticity matrix for matrix material
  cm44: C44 component of elasticity matrix for matrix material
  cp11: C11 component of elasticity matrix for second phase
  cp12: C12 component of elasticity matrix for second phase
  cp44: C44 component of elasticity matrix for second phase
  kx[Nx]: Fourier vector component in x-direction
  ky[Ny]: Fourier vector component in y-direction
  tmatx[Nx][Ny][2][2][2][2]: Values of Green's tensor at all grid points (real part only)
*/

void green_tensor_2d_mpi(int Nx, int Ny,
	double *kx, double *ky,
	double cm11, double cm12, double cm44,
	double cp11, double cp12, double cp44,
	double *tmatx,
	ptrdiff_t local_n0, ptrdiff_t local_0_start){
	
	//double omeg11[Nx][Ny];
	double *omeg11 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double omeg22[Nx][Ny];
	double *omeg22 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double omeg12[Nx][Ny];
	double *omeg12 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//
	double gmatx[2][2];
	//
	double dvect[2];
	
	/* Elastic constants of the homogeneous material, C0, which
	   are arithmetic averages of the elastic constants of the phases */
	double c11=0.5*(cm11+cp11);
	double c12=0.5*(cm12+cp12);
	double c44=0.5*(cm44+cp44);
	
	double chi=(c11-c12-2.0*c44)/c44;
	
	double rr;
	double d0;
	
	//int ij;
	int ijmpi;
	
	/* K0[i][k][zeta] = C0[i][j][k][l]*zeta[l]*zeta[k] (Eq.5.40)
	   where K0(zeta) is the acoustic tensor of the homogenous material.
	   N0(eta) = [K0(eta)]^-1 (Eq.5.42) */
	for(int i=0;i<local_n0;i++){
		for(int j=0;j<Ny;j++){
			//ij=i*Ny+j;
			ijmpi=(local_0_start+i)*Ny+j;
			
			rr = kx[local_0_start+i]*kx[local_0_start+i] + ky[j]*ky[j];
			
			d0 = c11*(rr*rr*rr) + chi*(c11+c12)*rr*(kx[local_0_start+i]*kx[local_0_start+i] * ky[j]*ky[j]);
			
			if(rr<1.0e-8){
				d0=1.0;
			}
			
			omeg11[ijmpi]=(c44*(rr*rr)+(c11-c44)*rr*ky[j]*ky[j])/(c44*d0);
			omeg22[ijmpi]=(c44*(rr*rr)+(c11-c44)*rr*kx[local_0_start+i]*kx[local_0_start+i])/(c44*d0);
			omeg12[ijmpi]=-(c12+c44)*kx[local_0_start+i]*ky[j]*rr/(c44*d0);
		}
	}
	
	for(int i=0;i<local_n0;i++){
		for(int j=0;j<Ny;j++){
			//ij=i*Ny+j;
			ijmpi=(local_0_start+i)*Ny+j;
			
			//Green's tensor
			gmatx[0][0]=omeg11[ijmpi];
			gmatx[0][1]=omeg12[ijmpi];
			gmatx[1][0]=omeg12[ijmpi];
			gmatx[1][1]=omeg22[ijmpi];
			
			//position vector
			dvect[0]=kx[local_0_start+i];
			dvect[1]=ky[j];
			
			//Green operator
			for(int kk=0;kk<2;kk++){
				for(int ll=0;ll<2;ll++){
					for(int ii=0;ii<2;ii++){
						for(int jj=0;jj<2;jj++){
							/* Eq.5.44
							G0[k][h][i][j] = (1/4)*(N0[h][i](zeta)*zeta[j]*zeta[k]
												   +N0[k][i](zeta)*zeta[j]*zeta[h]
												   +N0[h][i](zeta)*zeta[i]*zeta[k]
												   +N0[k][j](zeta)*zeta[i]*zeta[j])
							where G0[k][h][i][j] = tmatx[(((ij*2+k)*2+h)*2+i)*2+j]
							      N0[h][i](zeta) = gmatx[h][i]
							      zeta[j] = dvect[j] */
							tmatx[(((ijmpi*2+kk)*2+ll)*2+ii)*2+jj]=0.25*
								(gmatx[ll][ii]*dvect[jj]*dvect[kk]
								+gmatx[kk][ii]*dvect[jj]*dvect[ll]
								+gmatx[ll][jj]*dvect[ii]*dvect[kk]
								+gmatx[kk][jj]*dvect[ii]*dvect[ll]
								);
						}//jj
					}//ii
				}//ll
			}//kk
			
		}//j
	}//i
	
	return;
}