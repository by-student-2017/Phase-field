#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#include <fftw3.h>
//g++ martensite_v2_libfftw3.cpp -lfftw3
//#include <mpi.h> //mpi version
//#include <fftw3-mpi.h> //mpi version

#define DRND(x) ((double)(x)/RAND_MAX*rand())//Random number function setting
//#define ND 256					//The number of divisions per side of the computational domain in 
									//difference computation (a power of 2 due to the use of fast Fourier transform)
//#define IG 8					//2^IG=ND

	//int nd=ND, ndm=ND-1; 		//Define the number of difference divisions
								// (number of difference blocks) on one side of the calculation area, ND-1
	//int nd2=ND/2;				//Define ND/2: used in Fast Fourier Transform
	//int ig=IG;				//2^ig=ND
	double PI=3.14159;			//Pi
	double rr=8.3145;			//gas constant
	double time1;				//Number of calculation counts (proportional to time)
	int Nstep, iout;

	//double s1h[ND][ND], s2h[ND][ND];	//Phase field of martensite

	void ini000(double *s1h, double *s2h, int ND);	//Initial field setup subroutine
	void datsave(double *s1h, double *s2h, int ND);			//data save subroutine
	void datsave_paraview(double *s1h, double *s2h, int ND);//data save subroutine

//******* main program ******************************************
int main(void)
{
    int ND, IG;
	int nd, ndm, nd2, ig;
	
	double s1, s2;									//Phase field of martensite
	//double ep11h0[ND][ND], ep22h0[ND][ND];		//Transformational strain in tissue
	//double ep11qrh0[ND][ND],	ep11qih0[ND][ND];	//Fourier transform of constraint strain variation
	//double ep22qrh0[ND][ND],	ep22qih0[ND][ND];	//Fourier transform of constraint strain variation
	double s1k_chem, s1k_str;						//potential
	//double s1k_su[ND][ND];						//potential
	double s2k_chem, s2k_str;						//potential
	//double s2k_su[ND][ND];						//potential
	double c11, c12, c44, lam0, mu0, nu0; 			//elastic constant
	double c11c, c12c, c44c;			 			//elastic constant
	double eta_s1[4][4], eta_s2[4][4];				//Eigen distortion component
	//double ec11[ND][ND], ec22[ND][ND];			//Constraint strain variation (real space)
	double ep11T, ep22T;
	double ep11_0, ep22_0;							//Mean value of transformation strain in the structure
	double ep11_a, ep22_a, ep12_a, ep21_a;			//Strain caused by external force
	double sig11_a, sig22_a;						//external force
	double Z11ep, Z12ep, Z21ep, Z22ep;				//Coefficients for inverse Fourier transform
	double sum11, sum22;							//spatial integral of s1 and s2
	double s1ddtt, s2ddtt;							//Time variation of s1 and s2 (left side of evolution equation)
	double el_fac;									//Normalization constant of elastic constant

	int   i, j, k, l, ii, jj, kk, iii, jjj;			//integer
	int   p, q, m, n;								//integer
	int   ip, im, jp, jm, Nstep;					//integer
	double al, temp, delt;							//Computational domain, temperature, time step
	double time1max;								//Maximum calculation count (calculation end count)
	double b1, vm0, atom_n;							//Side length of difference block, molar volume, number of atoms in unit cell
	double smob;									//Mobility (relaxation coefficient of crystal transformation)
	double nxx, nyy, nxy, alnn;						//product of basis vectors in Fourier space, norm

	double AA0, AA1, AA2, AA3;						//factor in Gibbs energy
	double AA0e;									//factor in Gibbs energy
	double a1_c, b1_c, c1_c;						//lattice constant
	double a1_t, b1_t, c1_t;						//lattice constant
	double kappa_s1, kappa_s2;						//gradient energy factor
	double kappa_s1c, kappa_s2c;					//gradient energy factor
	double ds_fac;									//Fluctuation coefficient of crystal transformation
	double eta1, eta2;

//****** Setting calculation conditions and material constants ****************************************
	printf("---------------------------------\n");
	printf("read parameters from parameters.txt\n");
	FILE *fp;
	char name[40], comment[72];
	float param;
	float data[40];
	i = 0;
	fp = fopen("parameters.txt","r");
	while(fscanf(fp,"%s %e %[^\n] ",name, &param, comment) != EOF)
	{
		data[i] = param;
		printf("%d, %s %e \n",i,name,data[i]);
		i = i + 1;
	}
	fclose(fp);
	//
	ND      = int(data[0]);
	delt    = data[1];
	temp    = data[2];
	al      = data[3];	// [nm]
	smob    = data[4];
	ds_fac  = data[5];
	AA0e    = data[6];
	AA1     = data[7];
	AA2     = data[8];
	AA3     = data[9];
	kappa_s1c=data[10];
	kappa_s2c=data[11];
	vm0     = data[12];
	time1max= int(data[13]);
	Nstep   = int(data[14]);
	c11c    = data[15];
	c12c    = data[16];
	c44c    = data[17];
	lam0    = data[18];
	mu0     = data[19];
	eta1    = data[20];
	eta2    = data[21];
	sig22_a = data[22];
	printf("---------------------------------\n");
	//
	IG=int(log2(ND));
	nd=ND; 					
	ndm=ND-1;				//define ND-1
	nd2=ND/2;				//define ND/2 for FFT
	ig=IG;					//2^ig=ND
	//
	double *s1h      = (double *)malloc(sizeof(double)*( ND*ND + ND ));//Phase field of martensite
	double *s2h      = (double *)malloc(sizeof(double)*( ND*ND + ND ));//Phase field of martensite
	//
	double *ep11h0   = (double *)malloc(sizeof(double)*( ND*ND + ND ));//Transformational strain in tissue
	double *ep22h0   = (double *)malloc(sizeof(double)*( ND*ND + ND ));//Transformational strain in tissue
	double *ep11qrh0 = (double *)malloc(sizeof(double)*( ND*ND + ND ));//Fourier transform of constraint strain variation
	double *ep11qih0 = (double *)malloc(sizeof(double)*( ND*ND + ND ));//Fourier transform of constraint strain variation
	double *ep22qrh0 = (double *)malloc(sizeof(double)*( ND*ND + ND ));//Fourier transform of constraint strain variation
	double *ep22qih0 = (double *)malloc(sizeof(double)*( ND*ND + ND ));//Fourier transform of constraint strain variation
	//
	double *s1k_su   = (double *)malloc(sizeof(double)*( ND*ND + ND ));//potential
	double *s2k_su   = (double *)malloc(sizeof(double)*( ND*ND + ND ));//potential
	//
	double *ec11     = (double *)malloc(sizeof(double)*( ND*ND + ND ));//Constraint strain variation (real space)
	double *ec22     = (double *)malloc(sizeof(double)*( ND*ND + ND ));//Constraint strain variation (real space)
	//
	double *xi       = (double *)malloc(sizeof(double)*( ND*ND + ND ));//array of real and imaginary parts of the Fourier transform
	double *xr       = (double *)malloc(sizeof(double)*( ND*ND + ND ));//array of real and imaginary parts of the Fourier transform
	//
	const int fftsize = ND;
	fftw_complex *in, *out; // in[i][0] for real, in[i][1] for imag.
	fftw_plan plan, iplan;
	in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsize*fftsize);
	out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsize*fftsize);
	//plan = fftw_plan_dft_1d(fftsize, in, out, FFTW_FORWARD, FFTW_ESTIMATE); //FFT
	//fftw_execute(plan); //FFT
	//iplan = fftw_plan_dft_1d(fftsize, in, out, FFTW_BACKWARD, FFTW_ESTIMATE); //IFFT
	//fftw_execute(iplan); //IFFT
	//
	//After calc.: fftw_destroy_plan(plan); fftw_free(in); fftw_free(out);
	//
	//printf("DELT(0.2)=  ");	scanf(" %lf",&delt);	//Time step input	//delt=0.2;

	//temp=500.0;					//temperature (K)
	//al=500.0*1.0E-09;			//Computational area (m)
	al=al*1.0E-09;				//Computational area (m)
	b1=al/nd;					//Diff block length

	time1=0.0;					//Initial calculation count setting
	//time1max=1.0+1.0e+07;		//Setting the maximum calculation count

	//smob=1.0;					//Mobility (relaxation coefficient of crystal transformation)
	//ds_fac=0.01;				//Fluctuation coefficient of crystal transformation

	//AA0=1000.0/rr/temp;			//Chemical driving force for martensite transformation
	AA0=AA0e/rr/temp;			//Chemical driving force for martensite transformation
	//AA1=1.0;  AA2=3.0*AA1+12.0;  AA3=2.0*AA1+12.0;	//factor in Gibbs energy

	//kappa_s1=kappa_s2=5.0e-15/rr/temp/b1/b1;	//gradient energy factor
	kappa_s1=kappa_s1c/rr/temp/b1/b1;	//gradient energy factor
	kappa_s2=kappa_s2c/rr/temp/b1/b1;	//gradient energy factor

	//a1_c=b1_c=c1_c=3.5E-10;						//Lattice constant (m)
	//atom_n=4.0;  vm0=6.02E23*a1_c*b1_c*c1_c/atom_n;	//Calculation of molar volume (assuming fcc)

//*** Setting the Eigen distortion of the s1 field ***************
	//eta_s1[1][1]=0.05; eta_s1[2][2]=-0.05;
	eta_s1[1][1]=eta1; eta_s1[2][2]=eta2;
	eta_s1[3][3]=0.;
	eta_s1[1][2]=eta_s1[2][1]=eta_s1[1][3]=eta_s1[3][1]=eta_s1[2][3]=eta_s1[3][2]=0.;

//*** Setting the Eigen distortion of the s2 field ***************
	eta_s2[1][1]=eta_s1[2][2];
	eta_s2[2][2]=eta_s1[1][1];
	eta_s2[3][3]=0.;
	eta_s2[1][2]=eta_s2[2][1]=eta_s2[1][3]=eta_s2[3][1]=eta_s2[2][3]=eta_s2[3][2]=0.;

//***** elastic constant of Ni ****************************
	el_fac=1.0E+11*vm0/rr/temp;
	//c11=2.508*el_fac;
	//c44=1.235*el_fac;
	//c12=1.500*el_fac;
	c11=c11c*el_fac;
	c12=c12c*el_fac;
	c44=c44c*el_fac;
	//c12=c11-2.0*c44;
	if(lam0==0.0){
		lam0=c12;//Set Lahme constants from cij data
	}
	if(mu0==0.0){
		mu0=c44;//Set Lahme constants from cij data
	}
	nu0=lam0/2.0/(lam0+mu0);//Poisson's ratio
	printf("nu0= %f  \n", nu0);//Display the value of Poisson's ratio (to check the validity of the Lame constant)

//*** External force setting *******************************
 	//sig22_a=0.0;//In this calculation, set 0 because external force is not considered.
	ep11_a=-0.5*lam0/mu0/(3.0*lam0+2.0*mu0)*sig22_a;
	ep22_a=(lam0+mu0)/mu0/(3.0*lam0+2.0*mu0)*sig22_a;
	ep12_a=ep21_a=0.0;

//*** Setting up sin and cos tables, bit reversal tables, and initial fields ***************

	ini000(s1h, s2h, ND);		//Initial field setting

//**** Simulation start ******************************
//Nstep = 10;
iout = -1;
 plan = fftw_plan_dft_2d(fftsize, fftsize, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);	//For FFT
iplan = fftw_plan_dft_2d(fftsize, fftsize, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);	//For IFFT
start:;

	//if(time1<=100.){Nstep=10;} else{Nstep=200;}		//Changing the time interval for saving data
	//if((((int)(time1) % Nstep)==0)) {datsave(s1h, s2h, ND);} 	//Save tissue data every fixed repeat count
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(s1h, s2h, ND);} 	//Save tissue data every fixed repeat count

//***** gradient potential ***********************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}	//periodic boundary conditions
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
			//s1k_su[i][j]=-kappa_s1*(s1h[ip][j]+s1h[im][j]+s1h[i][jp]+s1h[i][jm]-4.0*s1h[i][j]);//Ž®(4.2.4)
			//s2k_su[i][j]=-kappa_s2*(s2h[ip][j]+s2h[im][j]+s2h[i][jp]+s2h[i][jm]-4.0*s2h[i][j]);
			s1k_su[i*ND+j]=-kappa_s1*(s1h[ip*ND+j]+s1h[im*ND+j]+s1h[i*ND+jp]+s1h[i*ND+jm]-4.0*s1h[i*ND+j]);//Ž®(4.2.4)
			s2k_su[i*ND+j]=-kappa_s2*(s2h[ip*ND+j]+s2h[im*ND+j]+s2h[i*ND+jp]+s2h[i*ND+jm]-4.0*s2h[i*ND+j]);
		}
	}

//**** Fourier transform of the Eigen strain field [equation (4.7)] ep11 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=ep11h0[i][j]=eta_s1[1][1]*s1h[i][j]+eta_s2[1][1]*s2h[i][j];
			//xi[i][j]=0.0;
			xr[i*ND+j]=ep11h0[i*ND+j]=eta_s1[1][1]*s1h[i*ND+j]+eta_s2[1][1]*s2h[i*ND+j];
			xi[i*ND+j]=0.0;
			//
			in[i*ND+j][0] = xr[i*ND+j]; //For FFT (real)
			in[i*ND+j][1] = xi[i*ND+j]; //For FFT (imag.)
		}
	}
	//qs=-1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	fftw_execute(plan); //For FFT
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			xr[i*ND+j] = out[i*ND+j][0]; //For FFT (real)
			xi[i*ND+j] = out[i*ND+j][1]; //For FFT (imag.)
			//
			//ep11qrh0[i][j]=xr[i][j];
			//ep11qih0[i][j]=xi[i][j];
			ep11qrh0[i*ND+j]=xr[i*ND+j];
			ep11qih0[i*ND+j]=xi[i*ND+j];
		}
	}
	//ep11qrh0[0][0]=ep11qih0[0][0]=0.;
	ep11qrh0[0]=ep11qih0[0]=0.0;

//**** Fourier transform of the Eigen strain field [equation (4.7)] ep22 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=ep22h0[i][j]=eta_s1[2][2]*s1h[i][j]+eta_s2[2][2]*s2h[i][j];
			//xi[i][j]=0.0;
			xr[i*ND+j]=ep22h0[i*ND+j]=eta_s1[2][2]*s1h[i*ND+j]+eta_s2[2][2]*s2h[i*ND+j];
			xi[i*ND+j]=0.0;
			//
			in[i*ND+j][0] = xr[i*ND+j]; //For FFT (real)
			in[i*ND+j][1] = xi[i*ND+j]; //For FFT (imag.)
		}
	}
	//qs=-1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	fftw_execute(plan); //For FFT
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){ 
			xr[i*ND+j] = out[i*ND+j][0]; //For FFT (real)
			xi[i*ND+j] = out[i*ND+j][1]; //For FFT (imag.)
			//
			//ep22qrh0[i][j]=xr[i][j];
			//ep22qih0[i][j]=xi[i][j];
			ep22qrh0[i*ND+j]=xr[i*ND+j];
			ep22qih0[i*ND+j]=xi[i*ND+j];
		}
	}
	//ep22qrh0[0][0]=ep22qih0[0][0]=0.;
	ep22qrh0[0]=ep22qih0[0]=0.0;

//*** Calculation of the average value of the Eigen strain field ***
	sum11=sum22=0.;
	for(i=0;i<=ndm;i++){
		//for(j=0;j<=ndm;j++){ sum11+=ep11h0[i][j];  sum22+=ep22h0[i][j]; }
		for(j=0;j<=ndm;j++){ sum11+=ep11h0[i*ND+j];  sum22+=ep22h0[i*ND+j]; }
	}
	ep11_0=sum11/nd/nd;  ep22_0=sum22/nd/nd;

//***** Calculation of constraint strain variation ec11 [formula (4.9)] *************************************
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj);
			if(alnn==0.){alnn=1.;}
			nxx=(double)ii/alnn*(double)ii/alnn;
			nyy=(double)jj/alnn*(double)jj/alnn;
			Z11ep=nxx*(2.0*(1.0-nu0)-nxx-nu0/(1.0-nu0)*nyy)/(1.0-2.0*nu0);
			Z12ep=nxx*(2.0*nu0      -nyy-nu0/(1.0-nu0)*nxx)/(1.0-2.0*nu0);
			//xr[i][j]=Z11ep*ep11qrh0[i][j]+Z12ep*ep22qrh0[i][j];
			//xi[i][j]=Z11ep*ep11qih0[i][j]+Z12ep*ep22qih0[i][j];
			xr[i*ND+j]=Z11ep*ep11qrh0[i*ND+j]+Z12ep*ep22qrh0[i*ND+j];
			xi[i*ND+j]=Z11ep*ep11qih0[i*ND+j]+Z12ep*ep22qih0[i*ND+j];
			//
			in[i*ND+j][0] = xr[i*ND+j]; //For IFFT (real)
			in[i*ND+j][1] = xi[i*ND+j]; //For IFFT (imag.)
		}
	}
	//qs=1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			xr[i*ND+j] = out[i*ND+j][0]/(fftsize*fftsize); //For IFFT (real)
			//
			//ec11[i][j]=xr[i][j];
			ec11[i*ND+j]=xr[i*ND+j];
		}
	}

//***** Calculation of constraint strain variation ec22 [Equation (4.9)] *****************************
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj);
			if(alnn==0.){alnn=1.;}
			nxx=(double)ii/alnn*(double)ii/alnn;
			nyy=(double)jj/alnn*(double)jj/alnn;
			Z21ep=nyy*(2.0*nu0      -nxx-nu0/(1.0-nu0)*nyy)/(1.0-2.0*nu0);
			Z22ep=nyy*(2.0*(1.0-nu0)-nyy-nu0/(1.0-nu0)*nxx)/(1.0-2.0*nu0);
			//xr[i][j]=Z21ep*ep11qrh0[i][j]+Z22ep*ep22qrh0[i][j];
			//xi[i][j]=Z21ep*ep11qih0[i][j]+Z22ep*ep22qih0[i][j];
			xr[i*ND+j]=Z21ep*ep11qrh0[i*ND+j]+Z22ep*ep22qrh0[i*ND+j];
			xi[i*ND+j]=Z21ep*ep11qih0[i*ND+j]+Z22ep*ep22qih0[i*ND+j];
			//
			in[i*ND+j][0] = xr[i*ND+j]; //For IFFT (real)
			in[i*ND+j][1] = xi[i*ND+j]; //For IFFT (imag.)
		}
	}
	//qs=1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			xr[i*ND+j] = out[i*ND+j][0]/(fftsize*fftsize); //For IFFT (real)
			//ec22[i][j]=xr[i][j];
			ec22[i*ND+j]=xr[i*ND+j];
		}
	}

//****** Potential calculation ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){

			//s1=s1h[i][j];  	s2=s2h[i][j];
			s1=s1h[i*ND+j];  	s2=s2h[i*ND+j];

//****** Calculation of chemical potential [equation (4.4)] ********************************
			s1k_chem=AA0*s1*(AA1-AA2*s1+AA3*(s1*s1+s2*s2));
			s2k_chem=AA0*s2*(AA1-AA2*s2+AA3*(s1*s1+s2*s2));

//****** Calculation of elastic potential [equation (4.8)] ********************************

			//ep11T=ep11h0[i][j]-ep11_0-ec11[i][j]-ep11_a;
			//ep22T=ep22h0[i][j]-ep22_0-ec22[i][j]-ep22_a;
			ep11T=ep11h0[i*ND+j]-ep11_0-ec11[i*ND+j]-ep11_a;
			ep22T=ep22h0[i*ND+j]-ep22_0-ec22[i*ND+j]-ep22_a;

			s1k_str=ep11T*((lam0+2.0*mu0)*eta_s1[1][1]+lam0*eta_s1[2][2])
				   +ep22T*((lam0+2.0*mu0)*eta_s1[2][2]+lam0*eta_s1[1][1]);
			s2k_str=ep11T*((lam0+2.0*mu0)*eta_s2[1][1]+lam0*eta_s2[2][2])
				   +ep22T*((lam0+2.0*mu0)*eta_s2[2][2]+lam0*eta_s2[1][1]);

//****** Calculate the time evolution of the phase field [equation (4.10)] ********************************
			//s1ddtt=-smob*(s1k_chem+s1k_su[i][j]+s1k_str);
			//s2ddtt=-smob*(s2k_chem+s2k_su[i][j]+s2k_str);
			//s1h[i][j]=s1h[i][j]+( s1ddtt+ds_fac*(2.0*DRND(1.)-1.0) )*delt;//Explicit method
			//s2h[i][j]=s2h[i][j]+( s2ddtt+ds_fac*(2.0*DRND(1.)-1.0) )*delt;
			s1ddtt=-smob*(s1k_chem+s1k_su[i*ND+j]+s1k_str);
			s2ddtt=-smob*(s2k_chem+s2k_su[i*ND+j]+s2k_str);
			s1h[i*ND+j]=s1h[i*ND+j]+( s1ddtt+ds_fac*(2.0*DRND(1.)-1.0) )*delt;//Explicit method
			s2h[i*ND+j]=s2h[i*ND+j]+( s2ddtt+ds_fac*(2.0*DRND(1.)-1.0) )*delt;

//*** Correction for s domain (0<=s<=1) ***
			//if(s1h[i][j]>=1.0){s1h[i][j]=1.0;}  if(s1h[i][j]<=0.0){s1h[i][j]=0.0;}
			//if(s2h[i][j]>=1.0){s2h[i][j]=1.0;}  if(s2h[i][j]<=0.0){s2h[i][j]=0.0;}
			if(s1h[i*ND+j]>=1.0){s1h[i*ND+j]=1.0;}
			if(s1h[i*ND+j]<=0.0){s1h[i*ND+j]=0.0;}
			if(s2h[i*ND+j]>=1.0){s2h[i*ND+j]=1.0;}
			if(s2h[i*ND+j]<=0.0){s2h[i*ND+j]=0.0;}
		}
	}

	//if(keypress()){return 0;}//Waiting for key

	time1=time1+1.0;								//Add calculation count
	if(time1<time1max){goto start;}	//Determining if the maximum count has been reached

end:;
  if(plan) fftw_destroy_plan(plan);		//For FFT
  if(iplan) fftw_destroy_plan(iplan);	//For IFFT
  fftw_free(in); fftw_free(out);		// For FFT and IFFT
  return 0;
}

//************ Initial field setup subroutine *************
void ini000(double *s1h, double *s2h, int ND)
{
	int i, j;
	//srand(time(NULL)); // random number initialization
	int nd=ND, ndm=ND-1, nd2=ND/2;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s1h[i][j]=0.2*DRND(1.0); s2h[i][j]=0.2*DRND(1.0);//Set the field with random numbers up to 20%
			//if(abs(j-nd2)<(nd/40)){s1h[i][j]=DRND(1.0); s2h[i][j]=DRND(1.0);}
			s1h[i*ND+j]=0.2*DRND(1.0);//Set the field with random numbers up to 20%
			s2h[i*ND+j]=0.2*DRND(1.0);//Set the field with random numbers up to 20%
		}
	}
}

//************ data save subroutine *******************************
void datsave(double *s1h, double *s2h, int ND)
{
	FILE *stream;		//Stream pointer setting
	int i, j;			//integer
	int nd=ND, ndm=ND-1, nd2=ND/2;

	stream = fopen("test.dat", "a");	//Open the write destination file for appending
	fprintf(stream, "%e\n", time1);		//Saving calculation counts
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//fprintf(stream, "%e  %e ", s1h[i][j], s2h[i][j]);//Field data storage
			fprintf(stream, "%e  %e ", s1h[i*ND+j], s2h[i*ND+j]);
		}
	}
	fprintf(stream, "\n");	//writing a newline
	fclose(stream);			//close file
}

void datsave_paraview(double *s1h, double *s2h, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j;
	int nd=ND, ndm=ND-1, nd2=ND/2;
	
	iout = iout + 1;
	printf("pf_result%06d.vtk \n",iout);
	sprintf(fName,"pf_result%06d.vtk",iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),(1));
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*1));
	fprintf(fp,"SCALARS Phase_field_1 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e ", s1h[i][j], s2h[i][j]);//Field data storage
			//fprintf(fp,"%10.6f\n", s1h[i][j]);
			fprintf(fp,"%10.6f\n", s1h[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Phase_field_2 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e ", s1h[i][j], s2h[i][j]);//Field data storage
			//fprintf(fp,"%10.6f\n", s2h[i][j]);
			fprintf(fp,"%10.6f\n", s2h[i*ND+j]);
		}
	}
	fclose(fp);
}
