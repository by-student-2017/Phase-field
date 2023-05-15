/* Program : 2D Phase-Field Simulation for 
   Spinodal Decomposition in Fe-Cr Alloy by GPU Computation.
   
   Programmer : Akinori Yamanaka (original version)
   Place : Depertment of Mechanical and Control Engineering Tokyo Institute of Technology
   Date : 7th, July, 2010 
   
   Modified version: By Student
   Place: 2-1 Hirosawa, Wako, Saitama, 351-0198, Japan
   Date: 15th/May/2023
   Test: Ubuntu 22.04 LTS
   
   Compiling: nvcc -O2 main-cpu.cu write_vtk_grid_values_3D.cu -o main-cpu.exe -lm -Wall
   Run: ./main-cpu.exe
   ParaView: time_XX.vtk
*/

#include <stdio.h>  //printf()
#include <stdlib.h> //rand() and malloc()
#include <math.h>   //mod() and -lm
//----- ----- -----
#include <cuda.h>   //GPU
/* #include <cuda.h> or
  #include "cuda_runtime.h"
  #include "device_launch_parameters.h" */
//----- ----- -----

#define TIMES 1
//----- ----- -----
#define NX 128*TIMES //Number of grid points in the x-direction
#define NY 128*TIMES //Number of grid points in the y-direction
#define NZ  16*TIMES //Number of grid points in the z-direction

void Kernel
(
	float *f, 
	float *fn,
	int    nx,
	int    ny,
	int    nz,
	float  rr,
	float  temp,
	float  L0,
	float  kapa_c,
	float  Da,
	float  Db,
	float  dt,
	float  dx,
	float  dy,
	float  dz
)
{
	int j, jx, jy, jz;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	for(jx=0; jx<nx; jx++){ //<-CPU | GPU-> jx = blockDim.z*blockIdx.z + threadIdx.z;
	for(jy=0; jy<ny; jy++){ //<-CPU | GPU-> jy = blockDim.y*blockIdx.y + threadIdx.y;
	for(jz=0; jz<nz; jz++){ //<-CPU | GPU-> jz = blockDim.x*blockIdx.x + threadIdx.x;
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
		   RT = rr*temp,
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
		   nab_mu, 
		   dfmdx, dfmdy, dfmdz, 
		   //----- ----- -----
		   Dab = Db/Da, 
		   mcc, dmc,
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
	// term1 = Atomic_interaction*(1-2*f) + RT*{log(f) - log(1-f)}
	mu_chc = L0*(1.0-2.0*fcc) + RT*( log(fcc) - log(1.0-fcc) ); //center: fcc
	//
	mu_chw = L0*(1.0-2.0*fcw) + RT*( log(fcw) - log(1.0-fcw) ); //center: fcw
	mu_che = L0*(1.0-2.0*fce) + RT*( log(fce) - log(1.0-fce) ); //center: fce
	//
	mu_chn = L0*(1.0-2.0*fcn) + RT*( log(fcn) - log(1.0-fcn) ); //center: fcn
	mu_chs = L0*(1.0-2.0*fcs) + RT*( log(fcs) - log(1.0-fcs) ); //center: fcs
	//
	mu_chu = L0*(1.0-2.0*fcu) + RT*( log(fcu) - log(1.0-fcu) ); //center: fcu
	mu_chd = L0*(1.0-2.0*fcd) + RT*( log(fcd) - log(1.0-fcd) ); //center: fcd
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// term2 = -gradient_energy_coefficient * Laplacian(f)
	mu_suc = -kapa_c*( (fce  + fcw  -2.0*fcc)/(dx*dx)
					 + (fcn  + fcs  -2.0*fcc)/(dy*dy)
					 + (fcu  + fcd  -2.0*fcc)/(dz*dz) ); //center: fcc
	//
	mu_suw = -kapa_c*( (fcc  + fcww -2.0*fcw)/(dx*dx)
					 + (fcnw + fcsw -2.0*fcw)/(dy*dy)
					 + (fcuw + fcdw -2.0*fcw)/(dz*dz) ); //fcc=fcwe, fcnw=fcwn, fcsw=fcws, //center: fcw
	mu_sue = -kapa_c*( (fcee + fcc  -2.0*fce)/(dx*dx)
					 + (fcne + fcse -2.0*fce)/(dy*dy)
					 + (fcue + fcde -2.0*fce)/(dz*dz) ); //fcc=fcew, fcne=fcen, fcse=fces, //center: fce
	//
	mu_sun = -kapa_c*( (fcne + fcnw -2.0*fcn)/(dx*dx)
					 + (fcnn + fcc  -2.0*fcn)/(dy*dy)
					 + (fcun + fcdn -2.0*fcn)/(dz*dz) ); //fcc=fcns, fcun=fcnu, fcdn=fcnd, //center: fcn
	mu_sus = -kapa_c*( (fcse + fcsw -2.0*fcs)/(dx*dx)
					 + (fcc  + fcss -2.0*fcs)/(dy*dy)
					 + (fcus + fcds -2.0*fcs)/(dz*dz) ); //fcc=fcsn, fcsu=fcus, fcds=fcsd, //center: fcs
	//
	mu_suu = -kapa_c*( (fcue + fcuw -2.0*fcu)/(dx*dx)
					 + (fcun + fcus -2.0*fcu)/(dy*dy)
					 + (fcuu + fcc  -2.0*fcu)/(dz*dz) ); //fcc=fcud, //center: fcu
	mu_sud = -kapa_c*( (fcde + fcdw -2.0*fcd)/(dx*dx)
					 + (fcdn + fcds -2.0*fcd)/(dy*dy)
					 + (fcc  + fcdd -2.0*fcd)/(dz*dz) ); //fcc=fcdu, //center: fcd
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
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
	// Laplacian(mu) = d^2(mu)/dx^2 + d^2(mu)/dy^2 + d^2(mu)/dz^2
	nab_mu = (mu_e + mu_w -2.0*mu_c)/(dx*dx)  // d^2(mu)/dx^2
		   + (mu_n + mu_s -2.0*mu_c)/(dy*dy)  // d^2(mu)/dy^2
		   + (mu_u + mu_d -2.0*mu_c)/(dz*dz); // d^2(mu)/dz^2
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// (df/dx) * d(mu)/dx, (x is related with e and w), (the center is fc.)
	dfmdx = ( (fce - fcw)/(2.0*dx) * (mu_e - mu_w)/(2.0*dx) );
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// (df/dy) * d(mu)/dy, (y is related with n and s), (the center is fc.)
	dfmdy = ( (fcn - fcs)/(2.0*dy) * (mu_n - mu_s)/(2.0*dy) );
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// (df/dz) * d(mu)/dz, (z is related with u and d), (the center is fc.)
	dfmdz = ( (fcu - fcd)/(2.0*dz) * (mu_u - mu_d)/(2.0*dz) );
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// Mobility, M = { (D_A/RT)*c + (D_B/RT)*(1-c) }*c*(1-c)
	//             = (D_a/RT)*{f + (D_B/D_A)*(1-f)}*f*(1-f)
	mcc = (Da/RT)*(fcc+Dab*(1.0-fcc))*fcc*(1.0-fcc); 
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// dM/df
	dmc = (Da/RT)*((1.0-Dab)*fcc*(1.0-fcc) + (fcc+Dab*(1.0-fcc))*(1.0-2.0*fcc)); 
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// df/dt = M*Laplacian(mu) + (dM/df)*( (df/dx) * d(mu)/dx + (df/dy) * d(mu)/dy + (df/dz) * d(mu)/dz )
	dfdt = mcc*nab_mu + dmc*(dfmdx + dfmdy + dfmdz); 
	fn[j] = f[j] + dfdt*dt;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	// If there are small variations, set the max and min values to the limits
	//if(fn[j]>0.999999){ fn[j]=0.999999; }
	//if(fn[j]<1.0e-6)  { fn[j]=1.0e-6; }
	
	}//end for(jz
	}//end for(jy
	}//end for(jx
}

void update(float **f, float **fn)
{
	float *tmp = *f;
	*f  = *fn;
	*fn = tmp;
}

void write_vtk_grid_values_3D(int Nx, int Ny, int Nz, float dx, float dy, float dz, int istep, float *data1);

int main(int argc, char** argv)
{ 
	float *f, *fn;
	int nx = NX, ny = NY, nz = NZ;
	int times = TIMES;
	char filename[] = "f000";
	
	int nstep=10000;    //Number of time integration steps
	int nprint=1000;    //Output frequency to write the results to file
	//----- ----- ----- -----
	float Lx = 3.0e-07*times, // Simulation length in x-direction [micro m]
		  Ly = 3.0e-07*times, // Simulation length in y-direction [micro m]
		  Lz = 3.0e-07*times*nz/sqrt(nx*ny), // Simulation length in z-direction [micro m]
		  //----- ----- ----- -----
		  dx = Lx/(float)nx, // Grid spacing between two grid pints in x-direction [nm]
		  dy = Ly/(float)ny, // Grid spacing between two grid pints in y-direction [nm]
		  dz = Lz/(float)nz, // Grid spacing between two grid pints in z-direction [nm]
		  //----- ----- ----- -----
		  c_0 = 0.4,    // Initial concentration (atomic fraction)
		  //----- ----- ----- -----
		  rr = 8.314,   // Gas constant [J/(mol*K)]
		  temp = 673.0, // Temperature [K]
		  RT = rr*temp,
		  //----- ----- ----- -----
		  L0 = 21020.8-9.31889*temp, // Atomic interaction [J/mol]
		  kapa_c = 1.2e-14,  // The value of gradient energy coefficients [J*m^2/mol]
		  //----- ----- ----- -----
		  Da = 1.0e-04*exp(-294000.0/RT), // Self-diffusion coefficient [m^2/s] (Fe)
		  Db = 2.0e-05*exp(-308000.0/RT), // Self-diffusion coefficient [m^2/s] (Cr)
		  //----- ----- ----- -----
		  dt = (dx*dx/Da)*0.1; // Time increment for the numerical integration [dimensionless]
	
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
	
	f  = (float *)malloc(nx*ny*nz*sizeof(float));
	fn = (float *)malloc(nx*ny*nz*sizeof(float));
	
	// Initialize the concentration filed f with random modulation
	for(int jz=0; jz<nz ; jz++){
		for(int jy=0; jy<ny ; jy++){
			for(int jx=0; jx<nx ; jx++){
				int j = (jz*ny + jy)*nx + jx; //j = nx*ny*jz + nx*jy + jx;
				float r = (float)rand()/(float)(RAND_MAX);
				f[j] = c_0 + 0.01*r;
			}
		}
	}
	
	//----- ----- ----- -----start:(This part is not really needed.)----- ----- ----- ----
	//Set recording time
	cudaEvent_t start, stop;
	
	//Initialization
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	//Start recording time
	cudaEventRecord(start);
	//----- ----- ----- -----end:(This part is not really needed.)----- ----- ----- ----
	
	for(int istep=0; istep<=nstep ; istep++){
		//calculate subroutine "Kernel" on CPU
		Kernel(f,fn,nx,ny,nz,rr,temp,L0,kapa_c,Da,Db,dt,dx,dy,dz);
		
		// replace f with new f (=fn)
		update(&f,&fn);
		//
		if(istep%nprint == 0){
			sprintf(filename,"f%03d",istep/nprint);
			
			//output vtk format
			write_vtk_grid_values_3D(nx,ny,nz,dx,dy,dz,istep,f);
			
			//show current step
			fprintf(stderr,"istep = %5d \n",istep);
		}
	}
	
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
	
	free(f);
	free(fn);
	
	return 0;
}