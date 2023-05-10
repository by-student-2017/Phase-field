#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#define NX 256
#define NY 256

__global__ void Kernel
(
 float *f,float *fn, int nx, int ny,
 float  rr, float temp,float L0,
 float kapa_c,float da,float db,float dt,float dx,float dy
)
{
	int j, jx, jy;
	float fcc, fce, fcw, fcs, fcn, fcnw, fcne, fcsw, fcse, fcww, fcee, fcnn, fcss, 
		  mu_chc, mu_chw, mu_che, mu_chn, mu_chs, mu_suc, mu_suw, mu_sue, mu_sun, 
		  mu_sus, mu_c, mu_w, mu_e, mu_n, mu_s, 
		  nab_mu, dfmdx, dfmdy, dab = db/da, mcc, dmc, dfdt ;
	
	jx = blockDim.y*blockIdx.y + threadIdx.y;
	jy = blockDim.x*blockIdx.x + threadIdx.x;
	j  = nx*jy + jx;
	
	fcc = f[j];
	
	if(jx == 0)    fcw = f[j+nx-1];
	else           fcw = f[j   -1];
	if(jx == nx-1) fce = f[j-nx+1];
	else           fce = f[j   +1];
	if(jy == 0)    fcs = f[j+nx*ny-nx];
	else           fcs = f[j      -nx];
	if(jy == ny-1) fcn = f[j-nx*ny+nx];
	else           fcn = f[j      +nx];
	
	if(jx == 0 && jy == ny-1)     { fcnw = f[nx-1];}
	else if(jx == 0 && jy  < ny-1){ fcnw = f[j+nx    +nx-1];}
	else if(jx  > 0 && jy == ny-1){ fcnw = f[j-nx*ny +nx-1];}
	else                          { fcnw = f[j       +nx-1];}
	
	if(jx == nx-1 && jy  < ny-1)     { fcne = f[j-nx    +nx+1];} 
	else if(jx  < nx-1 && jy == ny-1){ fcne = f[j-nx*ny +nx+1];}
	else if(jx == nx-1 && jy == ny-1){ fcne = f[0];}
	else                             { fcne = f[j       +nx+1];} 
	
	if(jx == 0 && jy >  0)     { fcsw = f[j+nx    -nx-1];}
	else if(jx  > 0 && jy == 0){ fcsw = f[j+nx*ny -nx-1];}
	else if(jx == 0 && jy == 0){ fcsw = f[nx*ny-1];} 
	else                       { fcsw = f[j       -nx-1];} 
	
	if(jx == nx-1 && jy == 0)     {fcse = f[nx*ny-1-nx+1];}             
	else if(jx == nx-1 && jy  > 0){ fcse = f[j-nx    -nx+1];}
	else if(jx <  nx-1 && jy == 0){ fcse = f[j+nx*ny -nx+1];}
	else                          { fcse = f[j       -nx+1];} 
	
	if(jx == 0) {fcww = f[j+nx-2];} 
	else if(jx == 1){ fcww = f[j+nx-2];}
	else            { fcww = f[j   -2];} 
	
	if(jx == nx - 2)     { fcee = f[j-nx+2];}
	else if(jx == nx - 1){ fcee = f[j-nx+2];}
	else                 { fcee = f[j   +2];} 
	
	if(jy == ny - 2)      { fcnn = f[j-nx*ny+nx+nx];}
	else if(jy == ny - 1) { fcnn = f[j-nx*ny+nx+nx];}
	else                  { fcnn = f[j      +nx+nx];} 
	
	if(jy == 0)      { fcss = f[j+nx*ny-nx-nx];}
	else if(jy == 1) { fcss = f[j+nx*ny-nx-nx];}
	else             { fcss = f[j      -nx-nx];} 
	
	mu_chc = L0*(1.0-2.0*fcc)+rr*temp*(log(fcc)-log(1.0-fcc));
	mu_chw = L0*(1.0-2.0*fcw)+rr*temp*(log(fcw)-log(1.0-fcw));
	mu_che = L0*(1.0-2.0*fce)+rr*temp*(log(fce)-log(1.0-fce));
	mu_chn = L0*(1.0-2.0*fcn)+rr*temp*(log(fcn)-log(1.0-fcn));
	mu_chs = L0*(1.0-2.0*fcs)+rr*temp*(log(fcs)-log(1.0-fcs));
	
	mu_suc = -kapa_c*(fce +fcw +fcn +fcs -4.0*fcc)/dx/dx;  
	mu_suw = -kapa_c*(fcc +fcww+fcnw+fcsw-4.0*fcw)/dx/dx;  
	mu_sue = -kapa_c*(fcee+fcc +fcne+fcse-4.0*fce)/dx/dx;  
	mu_sun = -kapa_c*(fcne+fcnw+fcnn+fcc -4.0*fcn)/dx/dx; 
	mu_sus = -kapa_c*(fcse+fcsw+fcc +fcss-4.0*fcs)/dx/dx;  
	
	mu_c = mu_chc + mu_suc; 
	mu_w = mu_chw + mu_suw; 
	mu_e = mu_che + mu_sue; 
	mu_n = mu_chn + mu_sun; 
	mu_s = mu_chs + mu_sus; 
	
	nab_mu = (mu_w + mu_e + mu_n + mu_s -4.0*mu_c)/dx/dx;  
	dfmdx = ((mu_w-mu_e)*(fcw-fce))/(4.0*dx*dx); 
	dfmdy = ((mu_n-mu_s)*(fcn-fcs))/(4.0*dy*dy); 
	
	mcc = (da/rr/temp)*(fcc+dab*(1.0-fcc))*fcc*(1.0-fcc); 
	dmc = (da/rr/temp)*((1.0-dab)*fcc*(1.0-fcc)
		+(fcc+dab*(1.0-fcc))*(1.0-2.0*fcc)); 
	
	dfdt = mcc*nab_mu + dmc*(dfmdx+dfmdy); 
	fn[j] = f[j]+dfdt*dt;
	
}

void update(float **f, float **fn)
{
	float *tmp = *f;
	*f  = *fn;
	*fn = tmp;
}
void write_vtk_grid_values_2D(int Nx, int Ny, float dx, float dy, int istep, float *data1);
int main(int argc, char** argv)
{
	float *f, *fn, *F;
	int nx = NX, ny = NY;
	int j; //int j, jx, jy, nstep;
	
	int nend = 10000, nout = 1000;
	float Lx = 3.0e-07, 
		  Ly = 3.0e-07,
		  dx = Lx/(float)nx,
		  dy = Ly/(float)ny,
		  c_0 = 0.4,
		  rr = 8.314, 
		  temp = 673.0,
		  L0 = 21020.8-9.31889*temp, 
		  kapa_c = 1.2e-14, 
		  da = 1.0e-04*exp(-294000.0/rr/temp), 
		  db = 2.0e-05*exp(-308000.0/rr/temp),
		  dt = (dx*dx/da)*0.1;  
	
	//CUT_DEVICE_INIT(argc, argv);
	
	f  = (float *)malloc(nx*ny*sizeof(float)); //GPU, CUDA, device
	fn = (float *)malloc(nx*ny*sizeof(float)); //GPU, CUDA, device
	
	cudaMalloc((void**)&f ,nx*ny*sizeof(float));
	cudaMalloc((void**)&fn,nx*ny*sizeof(float));
	
	F  = (float *)malloc(nx*ny*sizeof(float)); //CPU, host
	for(int jy=0; jy<ny ; jy++){
		for(int jx=0; jx<nx ; jx++){
			j = nx*jy + jx;
			float r = (float)rand()/(float)(RAND_MAX);
			F[j] = c_0 + 0.01*r;
		}
	}//CPU
	
	cudaMemcpy(f,F,nx*ny*sizeof(float),cudaMemcpyHostToDevice); //copy F(cpu) to f(cuda)
	//free(F);
	
	dim3 blocks(nx/16,ny/16,1);
	dim3 threads(16,16,1);
	
	//unsigned int timer;
	//cutCreateTimer(&timer);
	//cutResetTimer(timer);
	//cutStartTimer(timer);
	
	for(int nstep=0; nstep<=nend ; nstep++){
		
		Kernel<<<blocks, threads>>>(f,fn,nx,ny,rr,temp,L0,kapa_c,da,db,dt,dx,dy);
		//cudaThreadSynchronize();
		update(&f,&fn);
		
		if(nstep%nout == 0){
			fprintf(stderr,"nstep = %5d \n",nstep);
			cudaMemcpy(F,f,nx*ny*sizeof(float),cudaMemcpyDeviceToHost); //copy f(cuda) to F(cpu)
			write_vtk_grid_values_2D(nx,ny,dx,dy,nstep,F);
		}
	}
	//cutStopTimer(timer);       
	//float calc_time = cutGetTimerValue(timer)*1.0e-03;
	//printf("Calculation Time = %9.3e [sec]\n",calc_time);
	
	cudaFree(f);
	cudaFree(fn);
	
	free(F);
	
	return 0;
}
