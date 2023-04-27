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
	double *kx, double *ky, double *kz,
	double cm11, double cm12, double cm44,
	double cp11, double cp12, double cp44,
	double *omeg11, double *omeg22, double *omeg33,
	double *omeg12, double *omeg23, double *omeg13){
	
	/* Elastic constants of the homogeneous material, C0, which
	   are arithmetic averages of the elastic constants of the phases */
	double c11=0.5*(cm11+cp11);
	double c12=0.5*(cm12+cp12);
	double c44=0.5*(cm44+cp44); // "c44=0.5*(cm44+cm44) in text" is a mistake ?
	
	double chi=(c11-c12-2.0*c44)/c44;
	double chi2 = chi * chi;
	
	double kx1,ky1,kz1;
	double kx2,ky2,kz2;
	
	//double rr;
	double r2,r4;
	
	double d0;
	
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