#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cuda.h>
#include <cutil.h>

#define NX 256
#define NY 256

__global__ void Kernel
(
 float *f, float *fn, int nx, int ny,
 float  rr, float temp, float L0,
 float  kapa_c, float da, float db, float dt, float dx, float dy
)
{
  int j, jx, jy;
  float  fcc, fce, fcw, fcs, fcn, fcnw, fcne, fcsw, fcse, fcww, fcee, fcnn, fcss, 
         mu_chc, mu_chw, mu_che, mu_chn, mu_chs,
         mu_suc, mu_suw, mu_sue, mu_sun, mu_sus, 
         mu_c, mu_w, mu_e, mu_n, mu_s, 
         nab_mu, dfmdx, dfmdy, dab = db/da, mcc, dmc, dfdt ;

  int joff;
  int J0,J1,J2,J3,J4,J5,J6,J7,J8,J9,J10,J11;

  __shared__ float fs[20][20];

  jx = threadIdx.x + 2;
  jy = threadIdx.y + 2;
  joff = nx*(blockDim.y*blockIdx.y) + blockDim.x*blockIdx.x;
  j = joff + nx*threadIdx.y + threadIdx.x;

  fcc = f[j];
  fs[jx][jy] = fcc;

  if(blockIdx.y == 0) {J0 = nx*(ny-1)+blockDim.x*blockIdx.x + threadIdx.x,
                       J4 = J0 - nx;} 
  else                {J0 =  j - nx, 
                       J4 = J0 - nx;}

  if(blockIdx.y == gridDim.y - 1) {J1 = blockDim.x*blockIdx.x + threadIdx.x, 
                                   J5 = J1 + nx;} 
  else                            {J1 =  j + nx, 
                                   J5 = J1 + nx;}

  if(blockIdx.x == 0) {J2 = joff + nx*threadIdx.x + nx - 1,
                       J6 = J2 - 1;}
  else                {J2 = joff + nx*threadIdx.x - 1, 
                       J6 = J2 - 1;}

  if(blockIdx.x == gridDim.x - 1) {J3 = joff + nx*threadIdx.x + 15 - nx + 1,
                                   J7 = J3 + 1;}
  else                            {J3 = joff + nx*threadIdx.x + 16,
                                   J7 = J3 + 1;}

       if(blockIdx.x == 0 && blockIdx.y == gridDim.y - 1) { J8 = blockDim.x*16 -1;}
  else if(blockIdx.x  > 0 && blockIdx.y == gridDim.y - 1) { J8 = J1 - 1 ;}
  else if(blockIdx.x == 0 && blockIdx.y  < gridDim.y - 1) { J8 = j + nx + nx -1;}
  else                                                    { J8 = j + nx -1 ;}

       if(blockIdx.x == gridDim.x - 1 && blockIdx.y == gridDim.y - 1) { J9 = 0 ;}
  else if(blockIdx.x  < gridDim.x - 1 && blockIdx.y == gridDim.y - 1) { J9 = J1 + 1 ;}
  else if(blockIdx.x == gridDim.x - 1 && blockIdx.y  < gridDim.y - 1) { J9 = j  + 1 ;}
  else                                                                { J9 = j + nx +1 ;}

       if(blockIdx.x  > 0 && blockIdx.y == 0) { J10 = J0 - 1 ;}
  else if(blockIdx.x == 0 && blockIdx.y  > 0) { J10 =  j -1  ;}
  else if(blockIdx.x == 0 && blockIdx.y == 0) { J10 = nx*blockDim.x*blockDim.y - 1 ;}
  else                                        { J10 = j - nx - 1 ;}

       if(blockIdx.x == gridDim.x -1 && blockIdx.y == 0) { J11 = nx*blockDim.x*blockDim.y -1 - nx + 1;}
  else if(blockIdx.x  < gridDim.x -1 && blockIdx.y == 0) { J11 = J0 + 1  ;}
  else if(blockIdx.x == gridDim.x -1 && blockIdx.y  > 0) { J11 =  j - nx - nx + 1 ;}
  else                                                   { J11 = j - nx + 1 ;}

  if(threadIdx.y ==  0){ fs[jx][ 1] = f[J0], fs[jx][ 0] = f[J4];}
  if(threadIdx.y ==  1){ fs[ 1][jx] = f[J2], fs[ 0][jx] = f[J6];}
  if(threadIdx.y ==  2){ fs[18][jx] = f[J3], fs[19][jx] = f[J7];}
  if(threadIdx.y == 15){ fs[jx][18] = f[J1], fs[jx][19] = f[J5];}
  if(threadIdx.x ==  0 && threadIdx.y == 15) {fs[ 1][18] = f[J8];}
  if(threadIdx.x == 15 && threadIdx.y == 15) {fs[18][18] = f[J9];}
  if(threadIdx.x ==  0 && threadIdx.y ==  0) {fs[ 1][ 1] = f[J10];}
  if(threadIdx.x == 15 && threadIdx.y ==  0) {fs[18][ 1] = f[J11];}

  __syncthreads();

  fcc  = fs[jx  ][jy  ];
  fcw  = fs[jx-1][jy  ];
  fce  = fs[jx+1][jy  ];
  fcn  = fs[jx  ][jy+1];
  fcs  = fs[jx  ][jy-1];

  fcww = fs[jx-2][jy  ];
  fcee = fs[jx+2][jy  ];
  fcnn = fs[jx  ][jy+2];
  fcss = fs[jx  ][jy-2];

  fcnw = fs[jx-1][jy+1];
  fcne = fs[jx+1][jy+1];
  fcsw = fs[jx-1][jy-1];
  fcse = fs[jx+1][jy-1];

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
  dfmdy = ((mu_n-mu_s)*(fcn-fcs))/(4.0*dx*dx); 

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

int main(int argc, char** argv)
{ 
      float *f, *fn, *F;
      int nx = NX, ny = NY, 
          j, jx, jy, nstep;

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
      
      CUT_DEVICE_INIT(argc, argv);
      
      f  = (float *)malloc(nx*ny*sizeof(float));
      fn = (float *)malloc(nx*ny*sizeof(float));

      CUDA_SAFE_CALL(cudaMalloc((void**)&f ,nx*ny*sizeof(float)));
      CUDA_SAFE_CALL(cudaMalloc((void**)&fn,nx*ny*sizeof(float)));

      F  = (float *)malloc(nx*ny*sizeof(float));
      for(jy=0; jy<ny ; jy++){
       for(jx=0; jx<nx ; jx++){
        j = nx*jy + jx;
        float r = (float)rand()/(float)(RAND_MAX);
        F[j] = c_0 + 0.01*r;
         }
       }   
   
      CUDA_SAFE_CALL(cudaMemcpy(f,F,nx*ny*sizeof(float),
                     cudaMemcpyHostToDevice));
      free(F);
      
      dim3 blocks(nx/16,ny/16,1);
      dim3 threads(16,16,1);

     unsigned int timer;
     CUT_SAFE_CALL(cutCreateTimer(&timer));       
     CUT_SAFE_CALL(cutResetTimer(timer));       
     CUT_SAFE_CALL(cutStartTimer(timer));       

     for(nstep=0; nstep<=nend ; nstep++){

      Kernel<<<blocks, threads>>>(f,fn,nx,ny,rr,
                                  temp,L0,kapa_c,da,db,dt,dx,dy);
      cudaThreadSynchronize();
      update(&f,&fn);

      if(nstep%nout == 0) fprintf(stderr,"nstep = %4d\n",nstep);
      }

     CUT_SAFE_CALL(cutStopTimer(timer));       
     float calc_time = cutGetTimerValue(timer)*1.0e-03;       
     printf("Calculation Time = %9.3e [sec]\n",calc_time);

     CUDA_SAFE_CALL(cudaFree(f));
     CUDA_SAFE_CALL(cudaFree(fn));
      
     return 0;

}


