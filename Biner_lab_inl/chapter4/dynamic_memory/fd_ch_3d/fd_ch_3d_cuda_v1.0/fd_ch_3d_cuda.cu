/* Finite-difference phase-field code for
   solving Cahn-Hilliard equation */

/* This program solves the Cahn-Hilliard phase-field
   equation with finite difference algorithm using
   five-point stencil. The time integration is
   carried out using explicit Euler scheme. */

#include <stdio.h> //printf()
#include <stdlib.h> //rand()
#include <math.h> //mod() and -lm
#include <time.h>
//----- ----- -----
#include <cuda.h>   //GPU
/* #include <cuda.h> or
  #include "cuda_runtime.h"
  #include "device_launch_parameters.h" */
//----- ----- -----

#define BSX  8        //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
#define BSY  8        //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
#define BSZ  2        //Number of threads, 2^n=<32, BSX*BSY*BSZ <= 1024
#define TIMES 1
//----- ----- -----
#define NX 64*TIMES //Number of grid points in the x-direction
#define NY 64*TIMES //Number of grid points in the y-direction
#define NZ  2*TIMES //Number of grid points in the z-direction

//----- ----- ----- ----- ----- ----- -----
void micro_ch_pre_3d(int Nx, int Ny, int Nz, float c0, float *con);
//----- ----- -----
float calculate_energy_3d(int Nx, int Ny, int Nz, float *con, float grad_coef);
//----- ----- -----
void write_vtk_grid_values_3D(int nx, int ny, int nz, 
	float dx, float dy, float dz,
	int istep, float *data1);
//----- ----- ----- ----- ----- ----- -----

// Define subroutine "Kernel" for GPU (Device) calculation in detail
__global__ void Kernel
(
	float *f, 
	float *fn,
	int    nx,
	int    ny,
	int    nz,
	float  dx,
	float  dy,
	float  dz,
	float  dtime,
	float  mobility,
	float  grad_coef
)
{
	int j, jx, jy, jz;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	jx = blockDim.x*blockIdx.x + threadIdx.x; //<-GPU | CPU -> for(jx=0; jx<nx; jx++){
	jy = blockDim.y*blockIdx.y + threadIdx.y; //<-GPU | CPU -> for(jy=0; jy<ny; jy++){
	jz = blockDim.z*blockIdx.z + threadIdx.z; //<-GPU | CPU -> for(jz=0; jz<nz; jz++){
	j  = (jz*ny + jy)*nx + jx; //j = nx*ny*jz + nx*jy + jx;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	float  fcc,
		   fce,  fcw,  fcs,  fcn,  fcu,  fcd,
		   //----- ----- -----XY
		   fcnw, fcne,
		   fcsw, fcse,
		   //----- ----- -----XZ
		   fcuw, fcue,
		   fcdw, fcde,
		   //----- ----- -----YZ
		   fcun, fcus,
		   fcdn, fcds,
		   //----- ----- -----
		   fcww, fcee, fcnn, fcss, fcuu, fcdd,
		   //----- ----- ----- ----- ----- -----
		   mu_chc,
		   mu_chw, mu_che, mu_chn, mu_chs, mu_chu, mu_chd,
		   //----- ----- -----
		   mu_suc,
		   mu_suw, mu_sue, mu_sun, mu_sus, mu_suu, mu_sud,
		   //----- ----- -----
		   mu_c,
		   mu_w, mu_e, mu_n, mu_s, mu_u, mu_d, 
		   //----- ----- ----- ----- ----- -----
		   lap_mu, 
		   //----- ----- -----
		   dfdt ;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* Consider up to the second nearest neighbor. Therefore, 
	13 difference grid points are used. (#1 to #13)
	   The if statement is used because of periodic boundary conditions. */
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #1 (center: fcc)
	fcc = f[j];
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #2 (center: fcc)
	if(jx == 0)    fcw = f[j+(nx-1)];    //boundary condition at west edge
	else           fcw = f[j+(  -1)];    //non edge    [if(jx  > 0)]
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #3 (center: fcc)
	if(jx == nx-1) fce = f[j-(nx-1)];    //boundary condition at east edge
	else           fce = f[j-(  -1)];    //non edge    [if(jx  < nx-1)]
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #4 (center: fcc)
	if(jy == 0)    fcs = f[j+nx*(ny-1)]; //boundary condition at south edge
	else           fcs = f[j+nx*(  -1)]; //non edge    [if(jy  > 0)]
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #5 (center: fcc)
	if(jy == ny-1) fcn = f[j-nx*(ny-1)]; //boundary condition at north edge
	else           fcn = f[j-nx*(  -1)]; //non edge    [if(jy  < ny-1)]
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #6 (center: fcc)
	if(jz == 0)    fcu = f[j+nx*ny*(nz-1)]; //boundary condition at edge (down)
	else           fcu = f[j+nx*ny*(  -1)]; //non edge [if(jz  > 0)]
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #7 (center: fcc)
	if(jz == nz-1) fcd = f[j-nx*ny*(nz-1)]; //boundary condition at edge (up)
	else           fcd = f[j-nx*ny*(  -1)]; //non edge [if(jz  < nz-1)]
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----YX = XY
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #8 (center: fcc, and + n + w)(#5 and #2)
	/* e.g., "if(jy == ny-1 && jx == 0)" is f[ (j+nx*(-ny+1)) + (j+nx-1) - j] = f[j + nx*(-ny+1) +  nx-1] using above #5 and #2 condition.
	         "if(jy  < ny-1 && jx == 0)" is f[ j+nx*(   +1) + (j+nx-1) - j]   = f[j + nx*(   +1) +  nx-1] using above #5 and #2 condition. */
		 if(jy == ny-1 && jx == 0)   { fcnw = f[j-nx*(ny-1) +(nx-1)];}
	else if(jy  < ny-1 && jx == 0)   { fcnw = f[j-nx*(  -1) +(nx-1)];}
	else if(jy == ny-1 && jx  > 0)   { fcnw = f[j-nx*(ny-1) +(  -1)];}
	else                             { fcnw = f[j-nx*(  -1) +(  -1)];}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #9 (center: fcc, and + n + e)(#5 and #3)
		 if(jy  < ny-1 && jx == nx-1){ fcne = f[j-nx*(  -1) -(nx-1)];}
	else if(jy == ny-1 && jx  < nx-1){ fcne = f[j-nx*(ny-1) -(  -1)];}
	else if(jy == ny-1 && jx == nx-1){ fcne = f[j-nx*(ny-1) -(nx-1)];}
	else                             { fcne = f[j-nx*(  -1) -(  -1)];}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #10 (center: fcc, and + s + w)(#4 and #2)
		 if(jy >  0 && jx == 0)      { fcsw = f[j+nx*(  -1) +(nx-1)];}
	else if(jy == 0 && jx  > 0)      { fcsw = f[j+nx*(ny-1) +(  -1)];}
	else if(jy == 0 && jx == 0)      { fcsw = f[j+nx*(ny-1) +(nx-1)];}
	else                             { fcsw = f[j+nx*(  -1) +(  -1)];}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #11 (center: fcc, and + s + e)(#4 and #3)
		 if(jy == 0 && jx == nx-1)   { fcse = f[j+nx*(ny-1) -(nx-1)];}
	else if(jy  > 0 && jx == nx-1)   { fcse = f[j+nx*(  -1) -(nx-1)];}
	else if(jy == 0 && jx <  nx-1)   { fcse = f[j+nx*(ny-1) -(  -1)];}
	else                             { fcse = f[j+nx*(  -1) -(  -1)];}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----ZX = XZ
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #12 (center: fcc, and + u + w)(#7 and #2)
		 if(jz == nz-1 && jx == 0)   { fcuw = f[j-nx*ny*(nz-1) +(nx-1)];}
	else if(jz  < nz-1 && jx == 0)   { fcuw = f[j-nx*ny*(  -1) +(nx-1)];}
	else if(jz == nz-1 && jx  > 0)   { fcuw = f[j-nx*ny*(nz-1) +(  -1)];}
	else                             { fcuw = f[j-nx*ny*(  -1) +(  -1)];}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #9 (center: fcc, and + u + e)(#7 and #3)
		 if(jz  < nz-1 && jx == nx-1){ fcue = f[j-nx*ny*(  -1) -(nx-1)];}
	else if(jz == nz-1 && jx  < nx-1){ fcue = f[j-nx*ny*(nz-1) -(  -1)];}
	else if(jz == nz-1 && jx == nx-1){ fcue = f[j-nx*ny*(nz-1) -(nx-1)];}
	else                             { fcue = f[j-nx*ny*(  -1) -(  -1)];}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #10 (center: fcc, and + d + w)(#6 and #2)
		 if(jz >  0 && jx == 0)      { fcdw = f[j+nx*ny*(  -1) +(nx-1)];}
	else if(jz == 0 && jx  > 0)      { fcdw = f[j+nx*ny*(nz-1) +(  -1)];}
	else if(jz == 0 && jx == 0)      { fcdw = f[j+nx*ny*(nz-1) +(nx-1)];}
	else                             { fcdw = f[j+nx*ny*(  -1) +(  -1)];}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #11 (center: fcc, and + d + e)(#6 and #3)
		 if(jz == 0 && jx == nx-1)   { fcde = f[j+nx*ny*(nz-1) -(nx-1)];}
	else if(jz  > 0 && jx == nx-1)   { fcde = f[j+nx*ny*(  -1) -(nx-1)];}
	else if(jz == 0 && jx <  nx-1)   { fcde = f[j+nx*ny*(nz-1) -(  -1)];}
	else                             { fcde = f[j+nx*ny*(  -1) -(  -1)];}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----ZY = YZ
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #12 (center: fcc, and + u + s)(#7 and #4)
		 if(jz == nz-1 && jy == 0)   { fcus = f[j-nx*ny*(nz-1) +nx*(ny-1)];}
	else if(jz  < nz-1 && jy == 0)   { fcus = f[j-nx*ny*(  -1) +nx*(ny-1)];}
	else if(jz == nz-1 && jy  > 0)   { fcus = f[j-nx*ny*(nz-1) +nx*(  -1)];}
	else                             { fcus = f[j-nx*ny*(  -1) +nx*(  -1)];}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #13 (center: fcc, and + u + n)(#7 and #5)
		 if(jz  < nz-1 && jy == ny-1){ fcun = f[j-nx*ny*(  -1) -nx*(ny-1)];}
	else if(jz == nz-1 && jy  < ny-1){ fcun = f[j-nx*ny*(nz-1) -nx*(  -1)];}
	else if(jz == nz-1 && jy == ny-1){ fcun = f[j-nx*ny*(nz-1) -nx*(ny-1)];}
	else                             { fcun = f[j-nx*ny*(  -1) -nx*(  -1)];}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #14 (center: fcc, and + d + s)(#6 and #4)
		 if(jz >  0 && jy == 0)      { fcds = f[j+nx*ny*(  -1) +nx*(ny-1)];}
	else if(jz == 0 && jy  > 0)      { fcds = f[j+nx*ny*(nz-1) +nx*(  -1)];}
	else if(jz == 0 && jy == 0)      { fcds = f[j+nx*ny*(nz-1) +nx*(ny-1)];}
	else                             { fcds = f[j+nx*ny*(  -1) +nx*(  -1)];}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #15 (center: fcc, and + d + n)(#6 and #5)
		 if(jz == 0 && jy == ny-1)   { fcdn = f[j+nx*ny*(nz-1) -nx*(ny-1)];}
	else if(jz  > 0 && jy == ny-1)   { fcdn = f[j+nx*ny*(  -1) -nx*(ny-1)];}
	else if(jz == 0 && jy <  ny-1)   { fcdn = f[j+nx*ny*(nz-1) -nx*(  -1)];}
	else                             { fcdn = f[j+nx*ny*(  -1) -nx*(  -1)];}
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #20 (center: fcw)
		 if(jx == 0)     { fcww = f[j+(nx-2)];}    // edge(west)
	else if(jx == 1)     { fcww = f[j+(nx-2)];}    // edge(west,one inside)
	else                 { fcww = f[j+(  -2)];}    // non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #21 (center: fce)
		 if(jx == nx - 2){ fcee = f[j-(nx-2)];}    // edge(east)
	else if(jx == nx - 1){ fcee = f[j-(nx-2)];}    // edge(east, one inside)
	else                 { fcee = f[j-(  -2)];}    // non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #22 (center: fcn)
		 if(jy == ny - 2){ fcnn = f[j-nx*(ny-2)];} // edge(north)
	else if(jy == ny - 1){ fcnn = f[j-nx*(ny-2)];} // edge(north, one inside)
	else                 { fcnn = f[j-nx*(  -2)];} // non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #23 (center: fcs)
		 if(jy == 0)     { fcss = f[j+nx*(ny-2)];} // edge(south)
	else if(jy == 1)     { fcss = f[j+nx*(ny-2)];} // edge(south, one inside)
	else                 { fcss = f[j+nx*(  -2)];} // non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #24 (center: fcu)
		 if(jz == nz - 2){ fcuu = f[j-nx*ny*(nz-2)];} // edge(up)
	else if(jz == nz - 1){ fcuu = f[j-nx*ny*(nz-2)];} // edge(down, one inside)
	else                 { fcuu = f[j-nx*ny*(  -2)];} // non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #25 (center: fcd)
		 if(jz == 0)     { fcdd = f[j+nx*ny*(nz-2)];} // edge(up)
	else if(jz == 1)     { fcdd = f[j+nx*ny*(nz-2)];} // edge(down, one inside)
	else                 { fcdd = f[j+nx*ny*(  -2)];} // non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// term1 = dfdcon[i][j][k]
	/* float A=1.0;
	dfdcon[i][j][k]=A*( 2.0*con[i][j][k]*(1.0-con[i][j][k])*(1.0-con[i][j][k])
					   -2.0*con[i][j][k]*con[i][j][k]*(1.0-con[i][j][k]) );  */
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	mu_chc = 1.0*( 2.0*fcc*(1.0-fcc)*(1.0-fcc) -2.0*fcc*fcc*(1.0-fcc) ); //center: fcc
	//
	mu_chw = 1.0*( 2.0*fcw*(1.0-fcw)*(1.0-fcw) -2.0*fcw*fcw*(1.0-fcw) ); //center: fcw
	mu_che = 1.0*( 2.0*fce*(1.0-fce)*(1.0-fce) -2.0*fce*fce*(1.0-fce) ); //center: fce
	//
	mu_chn = 1.0*( 2.0*fcn*(1.0-fcn)*(1.0-fcn) -2.0*fcn*fcn*(1.0-fcn) ); //center: fcn
	mu_chs = 1.0*( 2.0*fcs*(1.0-fcs)*(1.0-fcs) -2.0*fcs*fcs*(1.0-fcs) ); //center: fcs
	//
	mu_chu = 1.0*( 2.0*fcu*(1.0-fcu)*(1.0-fcu) -2.0*fcu*fcu*(1.0-fcu) ); //center: fcu
	mu_chd = 1.0*( 2.0*fcd*(1.0-fcd)*(1.0-fcd) -2.0*fcd*fcd*(1.0-fcd) ); //center: fcd
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// term2 = -grad_coef * Laplacian(f)
	mu_suc = -grad_coef*( (fce  + fcw  -2.0*fcc)/(dx*dx)
						+ (fcn  + fcs  -2.0*fcc)/(dy*dy)
						+ (fcu  + fcd  -2.0*fcc)/(dz*dz) ); //center: fcc
	//
	mu_suw = -grad_coef*( (fcc  + fcww -2.0*fcw)/(dx*dx)
						+ (fcnw + fcsw -2.0*fcw)/(dy*dy)
						+ (fcuw + fcdw -2.0*fcw)/(dz*dz) ); //fcc=fcwe, fcnw=fcwn, fcsw=fcws, //center: fcw
	mu_sue = -grad_coef*( (fcee + fcc  -2.0*fce)/(dx*dx)
						+ (fcne + fcse -2.0*fce)/(dy*dy)
						+ (fcue + fcde -2.0*fce)/(dz*dz) ); //fcc=fcew, fcne=fcen, fcse=fces, //center: fce
	//
	mu_sun = -grad_coef*( (fcne + fcnw -2.0*fcn)/(dx*dx)
						+ (fcnn + fcc  -2.0*fcn)/(dy*dy)
						+ (fcun + fcdn -2.0*fcn)/(dz*dz) ); //fcc=fcns, fcun=fcnu, fcdn=fcnd, //center: fcn
	mu_sus = -grad_coef*( (fcse + fcsw -2.0*fcs)/(dx*dx)
						+ (fcc  + fcss -2.0*fcs)/(dy*dy)
						+ (fcus + fcds -2.0*fcs)/(dz*dz) ); //fcc=fcsn, fcsu=fcus, fcds=fcsd, //center: fcs
	//
	mu_suu = -grad_coef*( (fcue + fcuw -2.0*fcu)/(dx*dx)
						+ (fcun + fcus -2.0*fcu)/(dy*dy)
						+ (fcuu + fcc  -2.0*fcu)/(dz*dz) ); //fcc=fcud, //center: fcu
	mu_sud = -grad_coef*( (fcde + fcdw -2.0*fcd)/(dx*dx)
						+ (fcdn + fcds -2.0*fcd)/(dy*dy)
						+ (fcc  + fcdd -2.0*fcd)/(dz*dz) ); //fcc=fcdu, //center: fcd
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* Evaluate the terms in Eq.4.17, and
	   accumulate them in array dummy[Nx][Ny][Nz] for
	   each grid point in the simulation cell
	   dummy[ijk]=dfdcon-grad_coef*lap_con[ijk];               */
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// mu = dG/df = term1 + term2
	mu_c = mu_chc + mu_suc; // at current (jx,jy,jz) grid point, //center: fcc
	//
	mu_w = mu_chw + mu_suw; // at (jx-1,jy,jz) grid point, //center: fcw
	mu_e = mu_che + mu_sue; // at (jx+1,jy,jz) grid point, //center: fce
	//
	mu_n = mu_chn + mu_sun; // at (jx,jy+1,jz) grid point, //center: fcn
	mu_s = mu_chs + mu_sus; // at (jx,jy-1,jz) grid point, //center: fcs
	//
	mu_u = mu_chu + mu_suu; // at (jx,jy,jz+1) grid point, //center: fcu
	mu_d = mu_chd + mu_sud; // at (jx,jy,jz-1) grid point, //center: fcd
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* Calculate the Laplacing of the terms inside
	   the parenthesis in Eq.4.16
	   lap_dummy[ijk] = (hne + hnw -2.0*hnc)/(dx*dx)
					   +(hns + hnn -2.0*hnc)/(dy*dy)
					   +(hnu + hnd -2.0*hnc)/(dz*dz);          */
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// Laplacian(mu) = d^2(mu)/dx^2 + d^2(mu)/dy^2 + d^2(mu)/dz^2
	lap_mu = (mu_w + mu_e -2.0*mu_c)/(dx*dx)  // d^2(mu)/dx^2
		   + (mu_n + mu_s -2.0*mu_c)/(dy*dy)  // d^2(mu)/dy^2
		   + (mu_u + mu_d -2.0*mu_c)/(dz*dz); // d^2(mu)/dz^2
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* Explicit Euler time integration of
	   concentration field, Eq.4.16
	   con[ijk] = con[ijk] + dtime*mobility*lap_dummy[ijk];     */
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// df/dt = M*Laplacian(mu)
	dfdt = mobility*lap_mu;
	fn[j] = f[j] + dfdt*dtime;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	// If there are small variations, set the max and min values to the limits
	//if(fn[j]>0.999999){ fn[j]=0.999999; }
	//if(fn[j]<1.0e-6)  { fn[j]=1.0e-6; }
	
	//}//end for(jz <-GPU | CPU -> }//end for(jz
	//}//end for(jy <-GPU | CPU -> }//end for(jy
	//}//end for(jx <-GPU | CPU -> }//end for(jx
}

void update(float **f, float **fn)
{
	float *tmp = *f;
	*f  = *fn;
	*fn = tmp;
}

int main(){
	//open an output file for writing total bulk energy values
	FILE *out2=fopen("time_energy.out","w");
	
	//simulation cell parameters
	int nx=NX; //Number of grid points in the x-direction
	int ny=NY; //Number of grid points in the y-direction
	int nz=NZ; //Number of grid points in the y-direction
	
	//The distance between two grid points in x,y-direction
	float dx=1.0; //Grid spacing between two grid pints in x-direction
	float dy=1.0; //Grid spacing between two grid pints in y-direction
	float dz=1.0; //Grid spacing between two grid pints in z-direction
	
	//time integration parameters
	int nstep=10000; //Number of time integration steps
	int nprint=50;  //Output frequency to write the results to file
	float dtime=1.e-2; //Time increment for the numerical integration
	float ttime=0.0;    //Total time
	
	//material specific parameters
	
	//Initial concentrations of alloying elements
	float c0=0.40;       //Average composition
	float mobility=1.0;  //The value of mobility coefficient (dimensionless)
	float grad_coef=0.5; //The value of gradient energy coefficients [J(nm)^2/mol]
	
	float energy;
	
	//----- ----- ----- -----start:(This part is not really needed.)----- ----- ----- ----
	int nDevices;
	cudaGetDeviceCount(&nDevices);
	for (int i = 0; i < nDevices; i++){
		cudaDeviceProp prop;
		cudaGetDeviceProperties(&prop, i);
		printf("--------------------------------------------------\n");
		printf("Device Number: %d\n", i);
		printf("  Device name: %s\n", prop.name);
		printf("  Memory Clock Rate (KHz): %d\n",prop.memoryClockRate);
		printf("  Memory Bus Width (bits): %d\n",prop.memoryBusWidth);
		printf("  Peak Memory Bandwidth (GB/s): %f\n",2.0*prop.memoryClockRate*(prop.memoryBusWidth/8)/1.0e6);
		printf("--------------------------------------------------\n");
	}
	//----- ----- ----- -----end:(This part is not really needed.)----- ----- ----- ----
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	//float       con[Nx][Ny][Nz]; //concentraion
	float       *con = (float *)malloc(sizeof(float)*( nx*ny*nz ));
	//----- ----- ----- -----
	//prepare microstructure
	// Initialize the concentration filed with random modulation
	micro_ch_pre_3d(nx,ny,nz,c0,con);
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	float *con_d, *new_con_d; // name of dynamic memory for GPU, CUDA, device
	//
	con_d     = (float *)malloc(nx*ny*nz*sizeof(float)); //GPU, CUDA, device
	new_con_d = (float *)malloc(nx*ny*nz*sizeof(float)); //GPU, CUDA, device
	//
	cudaMalloc((void**)&con_d ,   nx*ny*nz*sizeof(float)); // define dynamic memory for GPU (device)
	cudaMalloc((void**)&new_con_d,nx*ny*nz*sizeof(float)); // define dynamic memory for GPU (device)
	//
	//copy F_h(cpu,host) to f_d(cuda,device)
	cudaMemcpy(con_d,con,nx*ny*nz*sizeof(float),cudaMemcpyHostToDevice); //con = con_h
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	int bsx=BSX, bsy=BSY, bsz=BSZ;     //Number of threads
	dim3 blocks(nx/bsx,ny/bsy,nz/bsz); //nx*ny*nz = blocks * threads
	dim3 threads(bsx,bsy,bsz);         //bsx*bsy*bsz <= 1024
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	
	//----- ----- ----- -----start:(This part is not really needed.)----- ----- ----- ----
	//Set recording time
	cudaEvent_t start, stop;
	
	//Initialization
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	//Start recording time
	cudaEventRecord(start);
	//----- ----- ----- -----end:(This part is not really needed.)----- ----- ----- ----
	
	//Time evolution of concentration filed
	for(int istep=0;istep<=nstep;istep++){
		
		//Update the total time
		ttime=ttime+dtime;
		
		//calculate subroutine "Kernel" on GPU
		Kernel<<<blocks, threads>>>(con_d,new_con_d,nx,ny,nz,dx,dy,dz,dtime,mobility,grad_coef);
		cudaDeviceSynchronize(); //<- new version | old version -> cudaThreadSynchronize();
		
		// replace con_d with new con_d
		update(&con_d,&new_con_d);
		
		//print results
		// If print frequency reached, print the results to file
		if(istep%nprint == 0){
			//copy f_d(cuda,device) to F_h(cpu,host)
			cudaMemcpy(con,con_d,nx*ny*nz*sizeof(float),cudaMemcpyDeviceToHost); //con = con_h
			
			/* Calculate the total bulk energy and print the result to
			   time_energy.out file for 2D plots */
			energy = calculate_energy_3d(nx,ny,nz,con,grad_coef);
			
			//print the average free energy density value to file
			fprintf(out2, "%14.6e %14.6e \n",ttime, energy);
			
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_3D(nx,ny,nz,dx,dy,dz,istep,con);
			
			printf("done step: %5d \n",istep);
		}//end if
	}//end for(istep
	
	//----- ----- ----- -----start:(This part is not really needed.)----- ----- ----- ----
	//Stop recording time
	cudaEventRecord(stop);
	
	//Wait all event
	cudaEventSynchronize(stop);
	
	//calculate time. time is [ms] unit.
	float milliseconds = 0;
	cudaEventElapsedTime(&milliseconds, start, stop);
	
	//Show computing time
	printf("Calculation Time = %9.3f [sec] \n",milliseconds*1.0e-03);
	
	//End processing
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	//----- ----- ----- -----end:(This part is not really needed.)----- ----- ----- ----
	
	cudaFree(con_d);
	cudaFree(new_con_d);
	
	free(con);
}
