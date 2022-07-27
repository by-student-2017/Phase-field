//
//
//  Program : 2D Phase-Field Simulation for 
//            Spinodal Decomposition in Fe-Cr Alloy
//            by CPU Computation
//
//  Programmer : Akinori Yamanaka
//   
//  Place : Depertment of Mechanical and Control Engineering
//          Tokyo Institute of Technology
//
//   Date : 7th, July, 2010

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <cuda.h>
//#include <cutil.h>
//#include <time.h> // for clock()
#define NX 128
#define NY 128

void Kernel
(
 float *f, 
 float *fn, 
 int    nx,
 int    ny,
 float  rr,
 float  temp,
 float  L0,
 float  kapa_c,
 float  da,
 float  db,
 float  dt,
 float  dx
)
{
  int j, jx, jy;
  float  fcc,  fce,  fcw,  fcs,  fcn, 
         fcnw, fcne, fcsw, fcse, fcww, 
         fcee, fcnn, fcss, 
         mu_chc, mu_chw, mu_che, mu_chn, mu_chs,
         mu_suc, mu_suw, mu_sue, mu_sun, mu_sus, 
         mu_c, mu_w, mu_e, mu_n, mu_s, 
         nab_mu, dfmdx, dfmdy, dab = db/da, mcc, dmc, dfdt ;

  for(jx=0; jx<nx; jx++){
  for(jy=0; jy<ny; jy++){
  j  = nx*jy + jx;

  fcc = f[j];

  if(jx == 0) fcw = f[j+nx-1];
  else        fcw = f[j   -1];
  if(jx == nx-1) fce = f[j-nx+1];
  else           fce = f[j   +1];
  if(jy == 0) fcs = f[j+nx*ny-nx];
  else        fcs = f[j      -nx];
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
  dfmdy = ((mu_n-mu_s)*(fcn-fcs))/(4.0*dx*dx); 

  mcc = (da/rr/temp)*(fcc+dab*(1.0-fcc))*fcc*(1.0-fcc); 
  dmc = (da/rr/temp)*((1.0-dab)*fcc*(1.0-fcc)
                     +(fcc+dab*(1.0-fcc))*(1.0-2.0*fcc)); 

  dfdt = mcc*nab_mu + dmc*(dfmdx+dfmdy); 
  fn[j] = f[j]+dfdt*dt;
   
   }
   }
} 

void update(float **f, float **fn)
{
  float *tmp = *f;
  *f  = *fn;
  *fn = tmp;
}

//void micro_avs
//(int nx, 
// int ny,
// float dx,
// float dy, 
// float *f, 
// char *fileBodyName)
//{
//  FILE *fp;
//  char fName[256];
//  int j,k;

//  sprintf(fName,"%s.inp",fileBodyName);
//  fp = fopen(fName,"w");

//  fprintf(fp,"%5d %5d %5d %5d %5d\n",nx*ny,(nx-1)*(ny-1),1,0,0);

//  k=0;
//  for(int jy=1; jy<=ny; jy++){
//   for(int jx=1; jx<=nx; jx++){
//    k=k+1;
//    fprintf(fp,"%5d %5d %5d %5d\n",k,jx,jy,0);
//   }
//  }
//  k=0;
//  for(int jy=1; jy<=ny-1; jy++){
//   for(int jx=1; jx<=nx-1; jx++){
//    k=k+1;
//    fprintf(fp,"%5d %5d quad %5d %5d %5d %5d\n",k,0,
//            nx*(jy-1)+jx,nx*(jy-1)+jx+1,nx*jy+jx+1,nx*jy+jx);
//   }
//  }

//  fprintf(fp,"1 1\n");
//  fprintf(fp,"phase\n");

//  k=0;
//  for(int jy=0; jy<ny; jy++){
//   for(int jx=0; jx<nx; jx++){
//    j = nx*jy + jx;
//    k=k+1;
//    fprintf(fp,"%5d %16.8f\n",k,f[j]);
//   }
//  }
//  fclose(fp);;
//}

void gnuplot
(int nx, 
 int ny,
 float dx,
 float dy, 
 float *f, 
 char *fileBodyName)
{
  FILE *fp;
  char fName[256];
  int   j;

  sprintf(fName,"%s.dat",fileBodyName);
  fp = fopen(fName,"w");
  
  for(int jx=0; jx<nx; jx++){
   for(int jy=0; jy<ny; jy++){
    j = nx*jy + jx;
    fprintf(fp,"%5d %5d %16.8f\n",jx,jy,f[j]);
   }
    fprintf(fp,"\n");
  }
  fclose(fp);;
}

int main(int argc, char** argv)
{
      //clock_t start = clock();	
      float *f, *fn;
      int nx = NX, ny = NY, 
          j, jx, jy, nstep;
      char filename[] = "f000";

      int nend = 1000, nout = 100;
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
      
      f  = (float *)malloc(nx*ny*sizeof(float));
      fn = (float *)malloc(nx*ny*sizeof(float));

      for(jy=0; jy<ny ; jy++){
       for(jx=0; jx<nx ; jx++){
        j = nx*jy + jx;
        float r = (float)rand()/(float)(RAND_MAX);
        f[j] = c_0 + 0.01*r;
         }
       }   

     //unsigned int timer;
     //CUT_SAFE_CALL(cutCreateTimer(&timer));       
     //CUT_SAFE_CALL(cutResetTimer(timer));       
     //CUT_SAFE_CALL(cutStartTimer(timer));       

     for(nstep=0; nstep<=nend ; nstep++){

      Kernel(f,fn,nx,ny,rr,temp,L0,kapa_c,da,db,dt,dx);
      update(&f,&fn);

      if(nstep%nout == 0){
       fprintf(stderr,"nstep = %4d\n",nstep);

       sprintf(filename,"f%03d",nstep/nout);
       //micro_avs(nx,ny,dx,dy,f,filename);
       gnuplot(nx,ny,dx,dy,f,filename);

        }        
      }

     //CUT_SAFE_CALL(cutStopTimer(timer));
     //float calc_time = cutGetTimerValue(timer)*1.0e-03;
     //
     //clock_t end = clock();
     //float calc_time = (end - start)/CLOCKS_PER_SEC;
     //printf("Calculation Time = %9.3e [sec]\n",calc_time);

     return 0;

}


