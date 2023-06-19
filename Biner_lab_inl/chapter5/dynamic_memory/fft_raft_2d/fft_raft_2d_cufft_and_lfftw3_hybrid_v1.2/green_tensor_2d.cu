#include <stdlib.h> //rand() and malloc

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

void green_tensor_2d(int Nx, int Ny,
	float *kx, float *ky,
	float cm11, float cm12, float cm44,
	float cp11, float cp12, float cp44,
	float *tmatx){
	
	//float omeg11[Nx][Ny];
	float *omeg11 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	//float omeg22[Nx][Ny];
	float *omeg22 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	//float omeg12[Nx][Ny];
	float *omeg12 = (float *)malloc(sizeof(float)*( Nx*Ny ));
	//
	float gmatx[2][2];
	//
	float dvect[2];
	
	/* Elastic constants of the homogeneous material, C0, which
	   are arithmetic averages of the elastic constants of the phases */
	float c11=0.5*(cm11+cp11);
	float c12=0.5*(cm12+cp12);
	float c44=0.5*(cm44+cp44);
	
	float chi=(c11-c12-2.0*c44)/c44;
	
	float rr;
	float d0;
	
	int ij;
	
	/* K0[i][j][zeta] = C0[i][j][k][l]*zeta[l]*zeta[k] (Eq.5.40)
	   where K0(zeta) is the acoustic tensor of the homogenous material.
	   zeta is the Fourier frequencies or vectors and i=sqrt(-1).
	   
	   N0(zeta) = [K0(zeta)]^-1 (Eq.5.42)
	   K0(zeta) = C[i][j][k][l]*zeta[l]*zeta[k] (Eq.5.40) 
	   
	   sigma[i][j] = C[i][j][k][l] * epsilon[k][l] (Eq.5.35)
	   where C[i][j][k][l] is position-dependent (also composition-dependent) elastic modulus tensor. */
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			ij=i*Ny+j;
			
			rr = kx[i]*kx[i] + ky[j]*ky[j];
			//
			d0 = c11*(rr*rr*rr) + chi*(c11+c12)*rr*(kx[i]*kx[i] * ky[j]*ky[j]);
			
			if(rr<1.0e-8){
				d0=1.0;
			}
			
			omeg11[ij]=(c44*(rr*rr)+(c11-c44)*rr*ky[j]*ky[j])/(c44*d0);
			omeg22[ij]=(c44*(rr*rr)+(c11-c44)*rr*kx[i]*kx[i])/(c44*d0);
			omeg12[ij]=-(c12+c44)*kx[i]*ky[j]*rr/(c44*d0);
		}
	}
	
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			ij=i*Ny+j;
			
			//Green's tensor
			gmatx[0][0]=omeg11[ij];
			gmatx[0][1]=omeg12[ij];
			gmatx[1][0]=omeg12[ij];
			gmatx[1][1]=omeg22[ij];
			
			//position vector
			dvect[0]=kx[i];
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
							tmatx[(((ij*2+kk)*2+ll)*2+ii)*2+jj]=0.25*
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

	/*----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	Note
	
	//alnn=sqrt((float)ii*(float)ii+(float)jj*(float)jj); // 2D
	alnn=sqrt((float)ii*(float)ii+(float)jj*(float)jj+(float)kk*(float)kk);
	if(alnn==0.0){alnn=1.0;}
	nec[1]=nx=(float)ii/alnn;	// n is the unit vector in the k direction, n=k/|k|
	nec[2]=ny=(float)jj/alnn;	// n is the unit vector in the k direction, n=k/|k|
	nec[3]=nz=(float)kk/alnn;	// n is the unit vector in the k direction, n=k/|k|
	//nec[3]=nz=0.0; // 2D
	
	// C[i][k][j][l]*n[j]*n[l], C: elastic modulus, n: unit vector
	C0[1][1][1][1]=c11;
	C0[2][2][2][2]=c22; -> c11
	C0[3][3][3][3]=c33; -> c11
	C0[1][2][1][2]=C0[1][2][2][1]=C0[2][1][1][2]=C0[2][1][2][1]=c66; -> c44
	C0[2][3][2][3]=C0[2][3][3][2]=C0[3][2][2][3]=C0[3][2][3][2]=c44;
	C0[1][3][1][3]=C0[1][3][3][1]=C0[3][1][1][3]=C0[3][1][3][1]=c55; -> c44
	C0[1][1][2][2]=C0[2][2][1][1]=c12=c21;
	C0[1][1][3][3]=C0[3][3][1][1]=c13=c31; -> c12
	C0[2][2][3][3]=C0[3][3][2][2]=c23=c32; -> c12
	
	a11=C0[1][1][1][1]*nx*nx+C0[1][2][1][2]*ny*ny+C0[1][3][1][3]*nz*nz;
	a22=C0[1][2][1][2]*nx*nx+C0[2][2][2][2]*ny*ny+C0[2][3][2][3]*nz*nz;
	a33=C0[3][1][3][1]*nx*nx+C0[2][3][2][3]*ny*ny+C0[3][3][3][3]*nz*nz;
	a12=(C0[1][1][2][2]+C0[1][2][1][2])*nx*ny;
	a23=(C0[2][2][3][3]+C0[2][3][2][3])*ny*nz;
	a31=(C0[3][3][1][1]+C0[3][1][3][1])*nx*nz;
	a21=a12;
	a32=a23;
	a13=a31;
	
	// cofactor
	b11=a22*a33-a23*a32;
	b22=a11*a33-a13*a31;
	b33=a11*a22-a12*a21;
	b12=-(a21*a33-a23*a31);
	b23=-(a11*a32-a12*a31);
	b31=-(a22*a31-a21*a32);
	b21=b12;
	b32=b23;
	b13=b31;
	
	// det (C[i][k][j][l]*n[j]*n[l])
	det1=a11*a22*a33+a12*a23*a31+a13*a32*a21-a13*a31*a22-a11*a23*a32-a33*a12*a21;
	if(det1==0.0){det1=1.0;}
	
	// inverse matrix
	om[1][1]=b11/det1;
	om[2][2]=b22/det1;
	om[3][3]=b33/det1;
	om[1][2]=om[2][1]=b12/det1;
	om[2][3]=om[3][2]=b23/det1;
	om[3][1]=om[1][3]=b31/det1;
	----- ----- ----- ----- ----- ----- ----- ----- ----- -----*/