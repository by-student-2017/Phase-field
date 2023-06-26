/* Finite-difference phase-field code for
   solving Allen-Cahn equation */

/* This program solves the non-concerved multicomponent
   Allen-Cahn equation with finite difference algorithm using
   five-point stencil. The time integration is carried out with
   explicit Euler scheme. */

#include <stdio.h>  //printf()
#include <stdlib.h> //rand() and malloc()
#include <math.h>   //mod() and -lm
#include <time.h>
//----- ----- -----
#include <cuda.h>   //GPU
/* #include <cuda.h> or
  #include "cuda_runtime.h"
  #include "device_launch_parameters.h" */
//----- ----- -----

#define BS 32        //Number of threads, 2^n=<32, BS*BS*1 <= 1024
#define TIMES 2
//----- ----- -----
#define NX 64*TIMES //Number of grid points in the x-direction
#define NY 64*TIMES //Number of grid points in the y-direction

//----- ----- ----- ----- ----- ----- -----
void init_grain_micro_2d(int Nx, int Ny,
	float dx, float dy, 
	int iflag, int ngrain,
	float *etas, int *glist);
//----- ----- -----
void write_vtk_grid_values_2D(int Nx, int Ny, 
	float dx, float dy,
	int istep, float *data1);
//----- ----- ----- ----- ----- ----- -----

// Define subroutine "Kernel" for GPU (Device) calculation in detail
__global__ void Kernel
(
	float *f, 
	float *fn,
	float *s,
	int    nx,
	int    ny,
	float  dx,
	float  dy,
	float  dtime,
	float  mobil,
	float  grcoef
)
{
	int j, jx, jy;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	jx = blockDim.x*blockIdx.x + threadIdx.x; //<-GPU | CPU -> for(jx=0; jx<nx; jx++){
	jy = blockDim.y*blockIdx.y + threadIdx.y; //<-GPU | CPU -> for(jy=0; jy<ny; jy++){
	j  = nx*jy + jx;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	float  fcc,
		   fce,  fcw,  fcs,  fcn,
		   //----- ----- -----
		   mu_chc,
		   //----- ----- -----
		   mu_suc,
		   //----- ----- -----
		   mu_c,
		   //----- ----- -----
		   dfdt,
		   //----- ----- -----
		   scc ;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* Consider up to the second nearest neighbor. Therefore, 
	5 difference grid points are used. (#1 to #5)
	   The if statement is used because of periodic boundary conditions. */
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #1 (center: fcc)
	fcc = f[j];
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #2 (center: fcc, and + w)
	if(jx == 0)    fcw = f[j+nx-1];    //boundary condition at west edge
	else           fcw = f[j   -1];    //non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #3 (center: fcc, and + e)
	if(jx == nx-1) fce = f[j-nx+1];    //boundary condition at east edge
	else           fce = f[j   +1];    //non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #4 (center: fcc, and + s)
	if(jy == 0)    fcs = f[j+nx*(+ny-1)]; //boundary condition at south edge
	else           fcs = f[j+nx*(   -1)]; //non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- #5 (center: fcc, and + n)
	if(jy == ny-1) fcn = f[j+nx*(-ny+1)]; //boundary condition at north edge
	else           fcn = f[j+nx*(   +1)]; //non edge
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	
	scc = s[j];
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	/* The derivative of the free energy for the grid point that is under
	   consideration in the main program
	   double A=1.0; double B=1.0;
	   double dfdeta_ij=A*(2.0*B*eta[ij]*sum + eta[ij]*eta[ij]*eta[ij] - eta[ij]); */
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	mu_chc = 1.0*( 2.0*1.0*fcc*scc + fcc*fcc*fcc - fcc ); //center: fcc
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// term2 = -grcoef * Laplacian(f)
	mu_suc = -grcoef*( (fce  + fcw  -2.0*fcc)/(dx*dx) + (fcn  + fcs  -2.0*fcc)/(dy*dy) ); //center: fcc
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// mu = dG/df = term1 + term2
	mu_c = mu_chc + mu_suc; // at current (jx,jy) grid point, //center: fcc
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	// df/dt = M*Laplacian(mu)
	dfdt = mobil*mu_c;
	fn[j] = f[j] - dfdt*dtime;
	//----- ----- ----- ----- ----- ----- ----- ----- ----- ----- 
	
	//}//end for(jy <-GPU | CPU-> }//end for(jy
	//}//end for(jx <-GPU | CPU-> }//end for(jx
}

//void update(float **f, float **fn)
//{
//	float *tmp = *f;
//	*f  = *fn;
//	*fn = tmp;
//}

int main(){
	//open an output file for writing total bulk energy values
	FILE *out2=fopen("area_frac.out","w");
	
	//simulation cell parameters
	int Nx; //Number of grid points in the x-direction
	int Ny; //Number of grid points in the y-direction
	
	int ngrain;
	
	//The distance between two grid points in x,y-direction
	float dx=0.5; //Grid spacing between two grid pints in x-direction
	float dy=0.5; //Grid spacing between two grid pints in y-direction
	
	//time integration parameters
	int nstep=100000; //Number of time integration steps
	int nprint=100; //Output frequency to write the results to file
	float dtime=0.005; //Time increment for the numerical integration
	float ttime=0.0;   //Total time
	
	//material specific parameters
	
	//Initial concentrations of alloying elements
	float mobil=5.0;  //The value of mobility coefficient
	float grcoef=0.1; //The value of gradient energy coefficients [J(nm)^2/mol]
	
	//----- ----- ----- -----
	int ij; //ij=i*Ny+j;
	//----- ----- ----- -----
	float grain_sum;
	//----- ----- ----- -----
	int ncount;
	//----- ----- ----- -----
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	/* Generate initial grain microstructure
	   iflag=1 is for bi-crystal and
	   iflag=2 is for polycrystal */
	int iflag=2;
	if(iflag==1){ Nx=NX; Ny=NY; ngrain=2; }
	if(iflag==2){
		//read polycrystal microstructure
		FILE *in=fopen("grain_25.inp","r");
		fscanf(in,"%5d %5d %5d ",&Nx,&Ny,&ngrain);
		fclose(in);
	}
	//----- ----- ----- -----
	int NxNy=Nx*Ny; //Total number of grid points in the simulation cell
	float     *eta = (float *)malloc(sizeof(float)*( NxNy ));
	float    *etas = (float *)malloc(sizeof(float)*( NxNy*ngrain ));
	int *glist = (int *)malloc(sizeof(int)*( ngrain ));
	float    *eta2 = (float *)malloc(sizeof(float)*( NxNy ));
	//
	// Summation in the free energy function (Eq.4.33)
	float     *sum = (float *)malloc(sizeof(float)*( NxNy ));
	//----- ----- ----- -----
	if(iflag==2){
		FILE *in=fopen("grain_25.inp","r");
		fscanf(in,"%5d %5d %5d ",&Nx,&Ny,&ngrain);
		//
		int nline=1;
		int ri, rj, rigrain;
		float reta;
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				fscanf(in,"%5d %5d %5d %f",&ri,&rj,&rigrain,&reta);
				//----- ----- ----- -----
				if( i != ri ){ printf("Don't match x data, Line %5d \n",nline); exit(1); }
				if( j != rj ){ printf("Don't match y data, LIne %5d \n",nline); exit(1); }
				nline = nline + 1;
				//----- ----- ----- -----
				ij=i*Ny+j;
				etas[ij*ngrain+rigrain]=reta;
				//----- ----- ----- -----
			}
		}
		fclose(in);
		//initialize glist
		for(int igrain=0;igrain<ngrain;igrain++){
			glist[igrain]=1.0;
		}
	}
	if(iflag==1){ init_grain_micro_2d(Nx,Ny,dx,dy,iflag,ngrain,etas,glist); }
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	
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
	float *eta_d, *new_eta_d, *sum_d; // name of dynamic memory for GPU, CUDA, device
	//----- ----- ----- -----
	eta_d      = (float *)malloc(NxNy*sizeof(float)); //GPU, CUDA, device
	new_eta_d  = (float *)malloc(NxNy*sizeof(float)); //GPU, CUDA, device
	sum_d      = (float *)malloc(NxNy*sizeof(float)); //GPU, CUDA, device
	//----- ----- ----- -----
	cudaMalloc((void**)&eta_d,     NxNy*sizeof(float)); // define dynamic memory for GPU (device)
	cudaMalloc((void**)&new_eta_d, NxNy*sizeof(float)); // define dynamic memory for GPU (device)
	cudaMalloc((void**)&sum_d,     NxNy*sizeof(float)); // define dynamic memory for GPU (device)
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	
	//----- ----- ----- ----- ----- ----- ----- ----- ----- -----
	int bs=BS; // Number of threads, 16 or 32
	dim3 blocks(Nx/bs,Ny/bs,1); //nx*ny = blocks * threads
	dim3 threads(bs,bs,1);      //bs*bs*1 <= 1024
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
	
	//Time integration
	for(int istep=0;istep<=nstep;istep++){
		
		//Update the total time
		ttime=ttime+dtime;
		
		//Loop over each grain
		for(int igrain=0;igrain<ngrain;igrain++){
			
			/* If glist is equal to one, which indicates that
			   the current grain area fraction is greater than 0.001,
			   continue the calculation. Otherwise, the current grain
			   does not exist anymore */
			if(glist[igrain]==1){
				
				/* Assign order parameters to temporary array eta[Nx][Ny] from
				   the common array etas[Nx][Ny][ngrain] for the current grain */
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						ij=i*Ny+j;
						eta[ij]=etas[ij*ngrain+igrain];
						//----- ----- ----- -----
						/* Determine the derivative of the free energy for
						   the current grain at the current grain point under
						   consideration, first three terms inside the bracket in Eq.4.35 */
						sum[ij] = 0.0;
						for(int jgrain=0;jgrain<ngrain;jgrain++){
							if(jgrain != igrain){
								sum[ij] = sum[ij] + etas[ij*ngrain+jgrain]*etas[ij*ngrain+jgrain];
							}
						}
						//----- ----- ----- -----
					}
				}
				//copy eta_d(cuda,device) to eta_h(cpu,host)
				cudaMemcpy(eta_d,eta,NxNy*sizeof(float),cudaMemcpyHostToDevice); //eta = eta_h
				cudaMemcpy(sum_d,sum,NxNy*sizeof(float),cudaMemcpyHostToDevice); //sum = sum_h
				
				// Calculate the Laplacian term in Eq.4.35
				//calculate subroutine "Kernel" on GPU
				Kernel<<<blocks, threads>>>(eta_d,new_eta_d,sum_d,Nx,Ny,dx,dy,dtime,mobil,grcoef);
				cudaDeviceSynchronize(); //<- new version | old version -> cudaThreadSynchronize();
				
				// replace eta_d with new eta_d
				//update(&eta_d,&new_eta_d);
				//cudaMemcpy(eta_d,new_eta_d,NxNy*sizeof(float),cudaMemcpyDeviceToDevice); //eta = eta_h
				
				//copy new_eta_d(cuda,device) to eta_h(cpu,host)
				cudaMemcpy(eta,new_eta_d,NxNy*sizeof(float),cudaMemcpyDeviceToHost); //eta = eta_h
				
				/* Calculate the total area of the current grain,
				   also return the order parameter values from
				   the temporary array eta[Nx][Ny] to 
				   common array etas[Nx][Ny][ngrain] */
				grain_sum=0.0;
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						ij=i*Ny+j;
						//----- ----- -----
						/* If there are small variations, 
						   set the max and min values to file */
						if(eta[ij]>=0.9999){
						   eta[ij]=0.9999;
						}
						if(eta[ij]<1.0e-6){
						   eta[ij]=1.0e-6;
						}
						//----- ----- -----
						etas[ij*ngrain+igrain]=eta[ij];
						grain_sum=grain_sum+eta[ij];
					}//end for(j
				}//end for(i
				
				//Check volume fraction of current grain
				/* Check the area fraction of the current grain.
				   If it is less than 0.001, set the its value in 
				   glist[ngrain] as zero which indicates that it is extinct. 
				   Also print message "grain # is eliminated" to screen. */
				grain_sum = grain_sum/NxNy;
				
				if(grain_sum<=0.001){
					glist[igrain]=0;
					printf("grain: No. %3d is eliminated \n",igrain);
				}
				
			}//end if(glist
		}//end igrain
		
		//print results
		// If print frequency reached, print the results to file
		if(fmod(istep,nprint)==0){
			//write vtk file & calculate are function of grains
			/* Prepare the data to be written to vtk file and
			   calculate the area fraction of each grain and
			   print them to file area_fract.out. */
			for(int i=0;i<Nx;i++){
				for(int j=0;j<Ny;j++){
					ij=i*Ny+j;
					eta2[ij]=0.0;
				}
			}
			fprintf(out2, "%14.6e ",ttime);
			
			for(int igrain=0;igrain<ngrain;igrain++){
				ncount=0;
				for(int i=0;i<Nx;i++){
					for(int j=0;j<Ny;j++){
						ij=i*Ny+j;
						//eta2[ij]=eta2[ij]+etas[ij*ngrain+igrain]*etas[ij*ngrain+igrain];
						eta2[ij]=eta2[ij]+etas[ij*ngrain+igrain]*etas[ij*ngrain+igrain]*igrain;
						if(etas[ij*ngrain+igrain]>=0.5){
							ncount=ncount+1;
						}//end if
					}//end for(j
				}//end for(i
				ncount=ncount/NxNy;
				fprintf(out2, "%5d ",ncount);
			}//end igrain
			fprintf(out2, "\n");
			
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_2D(Nx,Ny,dx,dy,istep,eta2);
			
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
	
	cudaFree(eta_d);
	cudaFree(new_eta_d);
	cudaFree(sum_d);
	
	free(eta);
	free(glist);
	free(etas);
	free(eta2);
	free(sum);
	
	return 0;
}