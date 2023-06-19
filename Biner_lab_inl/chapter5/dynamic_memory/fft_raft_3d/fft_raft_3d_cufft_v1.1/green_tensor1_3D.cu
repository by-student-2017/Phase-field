#include <stdlib.h> //rand() and malloc

/* The functions ( green_tensor1_3D() and green_tensor2_3D() ) given in 
   this section generate the Green's tensor in three-dimension. In order to
   avoid the exhaustion of the available memory, in three-dimension,
   it is generated with two separate functions. In the main program,
   the function green_tensor1_3D.c should be called only once after
   establishing FFT coefficients. Then, the function green_tensor2_3D.c should
   be called at every grid points, whenever, the Green's tensor is needed. */

/* This function evaluates the coefficients that
   will be used in the green_tensor2_3D.c. It needs to be called 
   only once after establishing the FFT coefficients in the main program. */

/* Variable and array list
  Nx: Number of grid points in the x-direction
  Ny: Number of grid points in the y-direction
  Nz: Number of grid points in the z-direction
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
*/

void green_tensor1_3D(int Nx, int Ny, int Nz,
	float *kx, float *ky, float *kz,
	float cm11, float cm12, float cm44,
	float cp11, float cp12, float cp44,
	float *omeg11, float *omeg22, float *omeg33,
	float *omeg12, float *omeg23, float *omeg13){
	
	/* Elastic constants of the homogeneous material, C0, which
	   are arithmetic averages of the elastic constants of the phases */
	float c11=0.5*(cm11+cp11);
	float c12=0.5*(cm12+cp12);
	float c44=0.5*(cm44+cp44); // "c44=0.5*(cm44+cm44) in text" is a mistake ?
	
	float chi=(c11-c12-2.0*c44)/c44;
	float chi2 = chi * chi;
	
	float kx1,ky1,kz1;
	float kx2,ky2,kz2;
	
	//float rr;
	float r2,r4;
	
	float d0;
	
	int ijk;
	
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			for(int k=0;k<Nz;k++){
				ijk=(i*Ny+j)*Nz+k;
				
				kx1 = kx[i];
				ky1 = ky[j];
				kz1 = kz[k];
				//
				kx2 = kx1 * kx1;
				ky2 = ky1 * ky1;
				kz2 = kz1 * kz1;
				
				r2 = kx2 + ky2 + kz2; //rr=r2
				//
				r4 = r2 * r2;
				
				d0 = c11*(r2 * r2 * r2) 
				   + chi*(c11+c12)*r2*( (kx2 * ky2)
									   +(ky2 * kz2)
									   +(kx2 * kz2)
									  )
				   +chi2*(c11+2.0*c12+c44)*(kx2 * ky2 * kz2);
				
				if(r2<1.0e-8){
					d0=1.0;
				}
				
				omeg11[ijk]=( c44*r4 + (c11+c12)*chi*(ky2*kz2) + (c11-c44)*r2*(ky2 + kz2) )/(c44*d0); // numerator(y and z, except x)
				omeg22[ijk]=( c44*r4 + (c11+c12)*chi*(kz2*kx2) + (c11-c44)*r2*(kz2 + kx2) )/(c44*d0); // numerator(z and x, except y)
				omeg33[ijk]=( c44*r4 + (c11+c12)*chi*(kx2*ky2) + (c11-c44)*r2*(kx2 + ky2) )/(c44*d0); // numerator(x and y, except z)
				//
				omeg12[ijk]=-(c12+c44)*(kx1*ky1)*(r2 + chi*kz2)/(c44*d0); // check: x*y, z
				omeg23[ijk]=-(c12+c44)*(ky1*kz1)*(r2 + chi*kx2)/(c44*d0); // check: y*z, x
				omeg13[ijk]=-(c12+c44)*(kx1*kz1)*(r2 + chi*ky2)/(c44*d0); // check: x*z, y
			}
		}
	}
	
	return;
}

	/*----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	Note
	
	alnn=sqrt((float)ii*(float)ii+(float)jj*(float)jj+(float)kk*(float)kk);
	if(alnn==0.0){alnn=1.0;}
	nec[1]=nx=(float)ii/alnn;	// n is the unit vector in the k direction, n=k/|k|
	nec[2]=ny=(float)jj/alnn;	// n is the unit vector in the k direction, n=k/|k|
	nec[3]=nz=(float)kk/alnn;	// n is the unit vector in the k direction, n=k/|k|
	
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