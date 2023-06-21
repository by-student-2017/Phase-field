/* 2D semi-implicit spectral phase-field code 
  for solving precipitation in Fe_Cu_Ni_Mn alloy */
  
// The dimension of energy were normalized with RT.
/* The time t was normalized with dx^2/Dcua(T),
where Dacu(T) is the diffusion constant of Cu in
the alpha phase at temperature T. */

/* The phase-field variable eta(r,t) characterize
  the phases distribution of alpha(bcc) and gamma(fcc) phase in
  the Cu precipitates and takes the values of 0 < h(eta) < 1, 
  corresponding h(eta)=0 for alpha(bcc) and
  h(eta)=1 for gamma(fcc) phase, respectively. */

#include <stdio.h>
#include <stdlib.h> //rand()
#include <math.h> //mod() and -lm
#include <time.h>

#include <fftw3.h>
//gcc test.c -lfftw3
//#include <mpi.h> //mpi version
//#include <fftw3-mpi.h> //mpi version

void init_FeCuMnNi_micro_2d(int Nx, int Ny, 
	double cu0, double mn0, double ni0,
	double *cu, double *mn, double *ni, double *orp);

void prepare_fft_2d(int Nx, int Ny, 
	double dx, double dy,
	double *kx, double *ky, 
	double *k2, double *k4);

void FeCuMnNi_free_energy_2d(int Nx, int Ny, 
	double *cu, double *mn, double *ni, double *orp, double tempr,
	double *dgdcu, double *dgdmn, double *dgdni, double *dgdor);

void write_vtk_grid_values_2D(int nx, int ny, 
	double dx, double dy, int istep, 
	double *data1, double *data2, double *data3, double *data4);

int main(){
	clock_t start, end;
	double compute_time;
	
	//get initial wall time
	start = clock();
	
	//simulation cell parameters
	int Nx=128;
	int Ny=128;
	
	//Total number of grid points in the simulation cell
	//int NxNy=Nx*Ny;
	
	//The distance between two grid points in x,y-direction
	double dx=0.5; //[nm] unit ?
	double dy=0.5; //[nm] unit ?
	
	//time integration parameters
	int nstep=5000;      //Number of time steps
	int nprint=50;       //Output frequency to write the results to file
	double dtime=1.0e-2; //Time increment for numerical integration
	double ttime=0.0;    //Total time
	//double coefA=1.0;
	
	//material specific parameters
	
	//Initial concentrations of alloying elements
	// e.g., 15 at.% Cu, 1 at.% Mn, 1 at.% Ni, (100-15-1-1) at.% Fe
	double cu0=0.15;
	double mn0=0.01;
	double ni0=0.01;
	
	//diffusivity parameters for mobility
	double gconst=8.314472; //The value of gas constant [J/mol/K]
	double tempr=823.0;   //Temperature [K]
	double dtempr=1000.0/(3600.0*6.0); //cooling rate [K/s]
	double RT=gconst*tempr;
	
	//gradient energy coefficients [J(nm)^2/mol]
	double Kc=5.0e3;   //[J(nm)^2/mol], 5.0e-15 [Jm^2/mol], Cahn-Hilliard eq. (conservation)
	double Keta=1.0e3; //[J(nm)^2/mol], 1.0e-15 [Jm^2/mol], Allen-Cahn eq. (non-conservation)
	
	//gradient energy of the composition
	//Cahn-Hilliard eq. (conservation)
	double grcoef_cu=Kc/RT; //0.68884=5000.0/(8.314472*873) for 873 [K]
	double grcoef_ni=Kc/RT; //0.68884=5000.0/(8.314472*873) for 873 [K]
	double grcoef_mn=Kc/RT; //0.68884=5000.0/(8.314472*873) for 873 [K]
	
	//gradient energy of the phase field
	//Allen-Cahn eq. (non-conservation)
	double grcoef_or=Keta/RT; //0.13729=1000/(8.314472*873) for 873 [K]
	
	// Diffusion constant [m^2/s]
	// A=alpha phase, G=gamma phase
	// D = D0*exp(-Q/RT)
	// D0[m^2/s], Q[J/mol]
	//
	double D0ACu=4.7e-5;
	double D0GCu=4.3e-5;
	//
	double QACu=2.44e5;
	double QGCu=2.80e5;
	//
	double D0ANi=1.4e-4;
	double D0GNi=1.08e-5;
	//
	//QANi=2.46e5;
	double QANi=2.56e5;
	double QGNi=2.74e5;
	//
	double D0AMn=1.49e-4;
	double D0GMn=2.78e-5;
	//
	//QAMn=2.33e5;
	double QAMn=2.64e5;
	double QGMn=2.64e5;
	//
	double DCuA=(D0ACu*exp(-QACu/RT));
	double DCuG=(D0GCu*exp(-QGCu/RT))/DCuA;
	//
	double DNiA=(D0ANi*exp(-QANi/RT))/DCuA;
	double DNiG=(D0GNi*exp(-QGNi/RT))/DCuA;
	//
	double DMnA=(D0AMn*exp(-QAMn/RT))/DCuA;
	double DMnG=(D0GMn*exp(-QGMn/RT))/DCuA;
	//
	DCuA=1.0;
	
	int ii;
	
	//double  cu[Nx][Ny];
	 double *cu = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double  mn[Nx][Ny];
	 double *mn = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double  ni[Nx][Ny];
	 double *ni = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double orp[Nx][Ny];
	double *orp = (double *)malloc(sizeof(double)*( Nx*Ny ));
	
	//prepare microstructure
	init_FeCuMnNi_micro_2d(Nx,Ny,cu0,mn0,ni0,cu,mn,ni,orp);
	
	//----- ----- ----- ----- ----- -----
	const int fftsizex = Nx, fftsizey = Ny;
	/* fftw_complex *in, *out; // in[i][0] for real, in[i][1] for imag.
	    in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	   out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	   fftw_plan plan, iplan;
	   plan = fftw_plan_dft_2d(fftsize, fftsize, in, out, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	  iplan = fftw_plan_dft_2d(fftsize, fftsize, in, out, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT */
	//
	//array, Cu
	fftw_complex *cuc, *cuck;
	 cuc = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	cuck = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	fftw_plan plan_cuc, iplan_cuck;
	plan_cuc = fftw_plan_dft_2d(fftsizex, fftsizey, cuc, cuck, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_cuck = fftw_plan_dft_2d(fftsizex, fftsizey, cuck, cuc, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//array, Mn
	fftw_complex *mnc, *mnck;
	 mnc = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	mnck = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	fftw_plan plan_mnc, iplan_mnck;
	plan_mnc = fftw_plan_dft_2d(fftsizex, fftsizey, mnc, mnck, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_mnck = fftw_plan_dft_2d(fftsizex, fftsizey, mnck, mnc, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//array, Ni
	fftw_complex *nic, *nick;
	 nic = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	nick = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	fftw_plan plan_nic, iplan_nick;
	plan_nic = fftw_plan_dft_2d(fftsizex, fftsizey, nic, nick, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_nick = fftw_plan_dft_2d(fftsizex, fftsizey, nick, nic, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//array, orp
	fftw_complex *orpc, *orpck;
	 orpc = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	orpck = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	fftw_plan plan_orpc, iplan_orpck;
	plan_orpc = fftw_plan_dft_2d(fftsizex, fftsizey, orpc, orpck, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	iplan_orpck = fftw_plan_dft_2d(fftsizex, fftsizey, orpck, orpc, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//
	//array, dG/dCu
	fftw_complex *dgdcuc, *dgdcuck;
	 dgdcuc = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	dgdcuck = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	fftw_plan plan_dgdcuc;
	plan_dgdcuc = fftw_plan_dft_2d(fftsizex, fftsizey, dgdcuc, dgdcuck, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//fftw_plan iplan_dgdcuck;
	//iplan_dgdcuck = fftw_plan_dft_2d(fftsizex, fftsizey, dgdcuck, dgdcuc, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//array, dG/dMn
	fftw_complex *dgdmnc, *dgdmnck;
	 dgdmnc = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	dgdmnck = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	fftw_plan plan_dgdmnc;
	plan_dgdmnc = fftw_plan_dft_2d(fftsizex, fftsizey, dgdmnc, dgdmnck, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//fftw_plan iplan_dgdmnck;
	//iplan_dgdmnck = fftw_plan_dft_2d(fftsizex, fftsizey, dgdmnck, dgdmnc, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//array, dG/dNi
	fftw_complex *dgdnic, *dgdnick;
	 dgdnic = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	dgdnick = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	fftw_plan plan_dgdnic;
	plan_dgdnic = fftw_plan_dft_2d(fftsizex, fftsizey, dgdnic, dgdnick, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//fftw_plan iplan_dgdnick;
	//iplan_dgdnick = fftw_plan_dft_2d(fftsizex, fftsizey, dgdnick, dgdnic, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//
	//array, dG/dor
	fftw_complex *dgdorc, *dgdorck;
	 dgdorc = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	dgdorck = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey);
	fftw_plan plan_dgdorc;
	plan_dgdorc = fftw_plan_dft_2d(fftsizex, fftsizey, dgdorc, dgdorck, FFTW_FORWARD,  FFTW_ESTIMATE); //For forward FFT
	//fftw_plan iplan_dgdorck;
	//iplan_dgdorck = fftw_plan_dft_2d(fftsizex, fftsizey, dgdorck, dgdorc, FFTW_BACKWARD, FFTW_ESTIMATE); //For inverse FFT
	//----- ----- ----- ----- ----- -----
	//
	//initialization
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			ii=i*Ny+j;
			//----- ----- ----- -----
			 cuc[ii][0] =  cu[ii];
			 cuc[ii][1] =  0.0;
			//----- ----- ----- -----
			 mnc[ii][0] =  mn[ii];
			 mnc[ii][1] =  0.0;
			//----- ----- ----- -----
			 nic[ii][0] =  ni[ii];
			 nic[ii][1] =  0.0;
			//----- ----- ----- -----
			orpc[ii][0] =  orp[ii];
			orpc[ii][1] =  0.0;
			//----- ----- ----- -----
		}
	}
	
	//double kx[Nx];
	double *kx = (double *)malloc(sizeof(double)*( Nx ));
	//double ky[Ny];
	double *ky = (double *)malloc(sizeof(double)*( Ny ));
	//double k2[Nx][Ny];
	double *k2 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double k4[Nx][Ny];
	double *k4 = (double *)malloc(sizeof(double)*( Nx*Ny ));
	
	//prepare fft (output: kx,ky,k2,k4)
	prepare_fft_2d(Nx,Ny,dx,dy,kx,ky,k2,k4); //get FFT coefficients
	
	//double dgdcu[Nx][Ny];
	double *dgdcu = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double dgdmn[Nx][Ny];
	double *dgdmn = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double dgdni[Nx][Ny];
	double *dgdni = (double *)malloc(sizeof(double)*( Nx*Ny ));
	//double dgdor[Nx][Ny];
	double *dgdor = (double *)malloc(sizeof(double)*( Nx*Ny ));
	
	double mcoef_cu;
	double mcoef_mn;
	double mcoef_ni;
	double mcoef_orp;
	
	double rtcoef=0.0;
	
	//evolve (Evolve the microstructure)
	for(int istep=0;istep<=nstep;istep++){
		
		//Update the total time
		ttime=ttime+(dtime*rtcoef);
		
		//Update the cooling temperature
		tempr=tempr-dtempr*(dtime*rtcoef);
		
		//Update the diffusion constant [m^2/s]
		DCuA=(D0ACu*exp(-QACu/RT));
		DCuG=(D0GCu*exp(-QGCu/RT))/DCuA;
		//
		DNiA=(D0ANi*exp(-QANi/RT))/DCuA;
		DNiG=(D0GNi*exp(-QGNi/RT))/DCuA;
		//
		DMnA=(D0AMn*exp(-QAMn/RT))/DCuA;
		DMnG=(D0GMn*exp(-QGMn/RT))/DCuA;
		//
		rtcoef=(dx*1e-9*dx*1e-9)/DCuA;
		//
		DCuA=1.0;
		
		/* Transform the alloying elements and
		   the order parameter values from
		   real space to Fourier space (forward FFT transformation) */
		//cuck=fft2(cuc);
		fftw_execute(plan_cuc);
		//mnck=fft2(mnc);
		fftw_execute(plan_mnc);
		//nick=fft2(nic);
		fftw_execute(plan_nic);
		//orpck=fft2(orpc);
		fftw_execute(plan_orpc);
		
		//derivative of free energy
		/* Calculate the derivatives of free energy for
		   each alloying elements and the order parameter */
		//Note: input array (real): cu, mn, ni, orp
		FeCuMnNi_free_energy_2d(Nx,Ny,cu,mn,ni,orp,tempr,dgdcu,dgdmn,dgdni,dgdor);
		//Note: output array (real): dgdcu, dgdmn, dgdni, dgdor
		
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				//----- ----- ----- -----
				dgdcuc[ii][0] = dgdcu[ii];
				dgdcuc[ii][1] = 0.0;
				//----- ----- ----- -----
				dgdmnc[ii][0] = dgdmn[ii];
				dgdmnc[ii][1] = 0.0;
				//----- ----- ----- -----
				dgdnic[ii][0] = dgdni[ii];
				dgdnic[ii][1] = 0.0;
				//----- ----- ----- -----
				dgdorc[ii][0] = dgdor[ii];
				dgdorc[ii][1] = 0.0;
				//----- ----- ----- -----
			}
		}
		
		/* Transform the derivative values from
		   real space to Fourier space (forward FFT transformation) */
		//dgdcuck=fft2(dgdcuc);
		fftw_execute(plan_dgdcuc);
		//dgdmnck=fft2(dgdmnc);
		fftw_execute(plan_dgdmnc);
		//dgdnick=fft2(dgdnic);
		fftw_execute(plan_dgdnic);
		//dgdorck=fft2(dgdorc);
		fftw_execute(plan_dgdorc);
		
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				
				//mobilities
				/* Calculate the mobility of each alloying elements,
				   Eq.5.24, based on the current values of
				   concentration and order parameter */
				/* Mi(eta,T) = ci0*(1.0-ci0)*{(1.0-eta)*Dia(T)/RT + eta*Dig(T)/RT} (Eq.5.24)
				   eta=orp[ii] */
				mcoef_cu=cu0*(1.0-cu0)*( (1.0-orp[ii])*DCuA + orp[ii]*DCuG );
				mcoef_mn=mn0*(1.0-mn0)*( (1.0-orp[ii])*DMnA + orp[ii]*DMnG );
				mcoef_ni=ni0*(1.0-ni0)*( (1.0-orp[ii])*DNiA + orp[ii]*DNiG );
				mcoef_orp=0.1;
				
				//time integration
				/* Semi-implicit time integration of each alloying elements
				   and the order parameter */
				// denominator is a real number
				// These are related by Eq.5.14 (Cahn-Hilliard for ci) or Eq.5.21 (Allen-Cahn for eta)
				//----- ----- ----- -----
				// real part
				 cuck[ii][0]= (cuck[ii][0]-dtime*k2[ii]*mcoef_cu*dgdcuck[ii][0])/
					(1.0+dtime*k4[ii]*mcoef_cu*grcoef_cu);
				 mnck[ii][0]= (mnck[ii][0]-dtime*k2[ii]*mcoef_mn*dgdmnck[ii][0])/
					(1.0+dtime*k4[ii]*mcoef_mn*grcoef_mn);
				 nick[ii][0]= (nick[ii][0]-dtime*k2[ii]*mcoef_ni*dgdnick[ii][0])/
					(1.0+dtime*k4[ii]*mcoef_ni*grcoef_ni);
				orpck[ii][0]=(orpck[ii][0]-dtime*k2[ii]*mcoef_orp*dgdorck[ii][0])/
					(1.0+dtime*k4[ii]*mcoef_orp*grcoef_or);
				//----- ----- ----- -----
				// imaginary part
				 cuck[ii][1]= (cuck[ii][1]-dtime*k2[ii]*mcoef_cu*dgdcuck[ii][1])/
					(1.0+dtime*k4[ii]*mcoef_cu*grcoef_cu);
				 mnck[ii][1]= (mnck[ii][1]-dtime*k2[ii]*mcoef_mn*dgdmnck[ii][1])/
					(1.0+dtime*k4[ii]*mcoef_mn*grcoef_mn);
				 nick[ii][1]= (nick[ii][1]-dtime*k2[ii]*mcoef_ni*dgdnick[ii][1])/
					(1.0+dtime*k4[ii]*mcoef_ni*grcoef_ni);
				orpck[ii][1]=(orpck[ii][1]-dtime*k2[ii]*mcoef_orp*dgdorck[ii][1])/
					(1.0+dtime*k4[ii]*mcoef_orp*grcoef_or);
				//----- ----- ----- -----
			}
		}
		
		/* Bring back the integrated values from
		Fourier space to real space (inverse FFT transformations) */
		//cuc=real(ifft2(cuck));
		fftw_execute(iplan_cuck);
		//mnc=real(ifft2(mnck));
		fftw_execute(iplan_mnck);
		//nic=real(ifft2(nic);
		fftw_execute(iplan_nick);
		//orpc=real(ifft2(orpck));
		fftw_execute(iplan_orpck);
		
		//for small deviations
		for(int i=0;i<Nx;i++){
			for(int j=0;j<Ny;j++){
				ii=i*Ny+j;
				//----- ----- ----- -----
				 cuc[ii][0] =  cuc[ii][0]/(fftsizex*fftsizey);
				 cuc[ii][1] =  cuc[ii][1]/(fftsizex*fftsizey);
				//----- -----
				 mnc[ii][0] =  mnc[ii][0]/(fftsizex*fftsizey);
				 mnc[ii][1] =  mnc[ii][1]/(fftsizex*fftsizey);
				//----- -----
				 nic[ii][0] =  nic[ii][0]/(fftsizex*fftsizey);
				 nic[ii][1] =  nic[ii][1]/(fftsizex*fftsizey);
				//----- -----
				orpc[ii][0] = orpc[ii][0]/(fftsizex*fftsizey);
				orpc[ii][1] = orpc[ii][1]/(fftsizex*fftsizey);
				//----- ----- ----- -----
				//Cu
				if(cuc[ii][0]>=0.9999){
					cuc[ii][0]=0.9999;
				}
				if(cuc[ii][0]<=0.0001){
					cuc[ii][0]=0.0001;
				}
				cuc[ii][1]=0.0;
				//----- ----- ----- -----
				//Mn
				if(mnc[ii][0]>=0.9999){
					mnc[ii][0]=0.9999;
				}
				if(mnc[ii][0]<=0.0001){
					mnc[ii][0]=0.0001;
				}
				mnc[ii][1]=0.0;
				//----- ----- ----- -----
				//Ni
				if(nic[ii][0]>=0.9999){
					nic[ii][0]=0.9999;
				}
				if(nic[ii][0]<=0.0001){
					nic[ii][0]=0.0001;
				}
				nic[ii][1]=0.0;
				//----- ----- ----- -----
				//order parameter
				if(orpc[ii][0]>=0.9999){
					orpc[ii][0]=0.9999;
				}
				if(orpc[ii][0]<=0.0001){
					orpc[ii][0]=0.0001;
				}
				orpc[ii][1]=0.0;
				//----- ----- ----- -----
				 cu[ii]  = cuc[ii][0];
				 mn[ii]  = mnc[ii][0];
				 ni[ii]  = nic[ii][0];
				orp[ii] = orpc[ii][0];
				//----- ----- ----- -----
			}
		}
		
		//print results
		/* If print frequency is reached, output the results to file */
		if(fmod(istep,nprint)==0){
			//write vtk file
			/* Write the results in vtk format for contour plots
			   to be viewed by using Paraview */
			write_vtk_grid_values_2D(Nx,Ny,dx,dy,istep,cu,mn,ni,orp);
			
			printf("done step: %5d, time %e [s], Temp: %8.3lf [K] \n",istep,ttime,tempr);
		}
		
	}//end of time step (evolve,for)
	
	//calculate the execution time and print it
	/* Calculate the compute time and print it to screen */
	end = clock();
	compute_time = ((double) (end - start)) / CLOCKS_PER_SEC;
	printf("Compute Time: %lf \n", compute_time);
	
	//----- ----- ----- ----- ----- -----
	fftw_destroy_plan(plan_cuc);
	fftw_destroy_plan(plan_mnc);
	fftw_destroy_plan(plan_nic);
	fftw_destroy_plan(plan_orpc);
	fftw_destroy_plan(iplan_cuck);
	fftw_destroy_plan(iplan_mnck);
	fftw_destroy_plan(iplan_nick);
	fftw_destroy_plan(iplan_orpck);
	//
	fftw_destroy_plan(plan_dgdcuc);
	fftw_destroy_plan(plan_dgdmnc);
	fftw_destroy_plan(plan_dgdnic);
	fftw_destroy_plan(plan_dgdorc);
	//----- ----- ----- ----- ----- -----
	fftw_free(cuc);
	fftw_free(cuck);
	fftw_free(mnc);
	fftw_free(mnck);
	fftw_free(nic);
	fftw_free(nick);
	fftw_free(orpc);
	fftw_free(orpck);
	//
	fftw_free(dgdcuc);
	fftw_free(dgdcuck);
	fftw_free(dgdmnc);
	fftw_free(dgdmnck);
	fftw_free(dgdnic);
	fftw_free(dgdnick);
	fftw_free(dgdorc);
	fftw_free(dgdorck);
	//----- ----- ----- ----- ----- -----
}
