#include <stdlib.h> //rand() and malloc

/* The functions ( green_tensor1_3D() and green_tensor2_3D() ) given in 
   this section generate the Green's tensor in three-dimension. In order to
   avoid the exhaustion of the available memory, in three-dimension,
   it is generated with two separate functions. In the main program,
   the function green_tensor1_3D.c should be called only once after
   establishing FFT coefficients. Then, the function green_tensor2_3D.c should
   be called at every grid points, whenever, the Green's tensor is needed. */

/* This function evaluates the Green's tensor. It needs to be called at
   every grid points, whenever the value of Green's tensor is needed in
   the main program. */

/* Variable and array list
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  Nz: Number of grid points in the z-direction
  i: Current grid number in the x-direction
  j: Current grid number in the y-direction
  k: Current grid number in the z-direction
  kx[Nx]: Fourier vector component in x-direction
  ky[Ny]: Fourier vector component in y-direction
  kz[Nz]: Fourier vector component in z-direction
  cm11,cm12,cm44: Elasticity parameters of the matrix phase
  cp11,cp12,cp44: Elasticity parameters of the second phase
  omeg11[Nx][Ny][Nz]: Coefficient needed for the Green's tensor
  omeg22[Nx][Ny][Nz]: Coefficient needed for the Green's tensor
  omeg33[Nx][Ny][Nz]: Coefficient needed for the Green's tensor
  omeg12[Nx][Ny][Nz]: Coefficient needed for the Green's tensor
  omeg23[Nx][Ny][Nz]: Coefficient needed for the Green's tensor
  omeg13[Nx][Ny][Nz]: Coefficient needed for the Green's tensor
  tmatx[3][3][3][3]: Green's tensor
*/

void green_tensor2_3D(int Nx, int Ny, int Nz,
	double *kx, double *ky, double *kz,
	double *omeg11, double *omeg22, double *omeg33,
	double *omeg12, double *omeg23, double *omeg13,
	int i, int j, int k,
	double *tmatx){
	
	int ijk=(i*Ny+j)*Nz+k;
	
	double gmatx[3][3];
	//
	gmatx[0][0] = omeg11[ijk];
	gmatx[0][1] = omeg12[ijk];
	gmatx[0][2] = omeg13[ijk];
	//
	gmatx[1][0] = omeg12[ijk];
	gmatx[1][1] = omeg22[ijk];
	gmatx[1][2] = omeg23[ijk];
	//
	gmatx[2][0] = omeg13[ijk];
	gmatx[2][1] = omeg23[ijk];
	gmatx[2][2] = omeg33[ijk];
	
	//position vector
	double dvect[3];
	//
	dvect[0] = kx[i];
	dvect[1] = ky[j];
	dvect[2] = kz[k];
	
	int klij;
	
	//Green operator
	for(int kk=0;kk<3;kk++){
		for(int ll=0;ll<3;ll++){
			for(int ii=0;ii<3;ii++){
				for(int jj=0;jj<3;jj++){
					klij=((kk*3+ll)*3+ii)*3+jj;
					/* Eq.5.44
						G0[k][h][i][j] = (1/4)*(N0[h][i](zeta)*zeta[j]*zeta[k]
											   +N0[k][i](zeta)*zeta[j]*zeta[h]
											   +N0[h][i](zeta)*zeta[i]*zeta[k]
											   +N0[k][j](zeta)*zeta[i]*zeta[j])
						where G0[k][h][i][j] = tmatx[((k*3+h)*3+i)*3+j]
						      N0[h][i](zeta) = gmatx[h][i]
						      zeta[j] = dvect[j] */
					//
					tmatx[klij]=0.25*(gmatx[ll][ii]*dvect[jj]*dvect[kk]
									 +gmatx[kk][ii]*dvect[jj]*dvect[ll]
									 +gmatx[ll][jj]*dvect[ii]*dvect[kk]
									 +gmatx[kk][jj]*dvect[ii]*dvect[ll]
									 );
					//
				}
			}
		}
	}
	
	return;
}