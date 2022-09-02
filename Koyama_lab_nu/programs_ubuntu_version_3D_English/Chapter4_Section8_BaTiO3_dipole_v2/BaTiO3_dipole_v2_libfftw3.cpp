#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#include <fftw3.h>
//g++ martensite_v2_libfftw3.cpp -lfftw3
//#include <mpi.h> //mpi version
//#include <fftw3-mpi.h> //mpi version

#define DRND(x) ((double)(x)/RAND_MAX*rand())//Random number function setting

//#define ND 128					//The number of divisions per side of the computational domain in difference computation
									// (a power of 2 due to the use of fast Fourier transform)
//#define IG 7						// 2^IG=ND

	//int nd=ND, ndm=ND-1; 			//Define the number of difference divisions
									// (number of difference blocks) on one side of the calculation area, ND-1
	//int nd2=ND/2;				 	//Define ND/2: used in Fast Fourier Transform
	//int ig=IG;					//2^ig=ND
	double PI=3.14159;				//Pi
	double RR=8.3145;				//gas constant
	double time1;					//Number of calculation counts (proportional to time)
	double ep0=8.8541878e-12; 		//Vacuum permittivity (F/m)
	double Peq;						//equilibrium value of the moment of polarization
	int iout;
	
	//double s1h[ND][ND][ND], s2h[ND][ND][ND];//polarization moment in x direction, polarization moment in y direction

	void ini000(double *s1h, double *s2h, double *s3h, int ND);			//initial profile of the polarization moment at time 0
	void datsave(double *s1h, double *s2h, double *s3h, int ND);	//data save subroutine
	void datsave_paraview(double *s1h, double *s2h, double *s3h, int ND);//data save subroutine
	void datin(double *s1h, double *s2h, double *s3h, int ND);	//Initial data loading

//******* main program ******************************************
int main(void)
{
    int ND;
	int nd, ndm, nd2, ig;
	
	int   i, j, k, l, ii, jj, kk, Nstep;	//integer
	int   ip, im, jp, jm, kp, km;			//integer

	//double s1qrh[ND][ND][ND], s1qih[ND][ND][ND];	//Fourier transform of s1 (real part, imaginary part)
	//double s2qrh[ND][ND][ND], s2qih[ND][ND][ND];	//Fourier transform of s2 (real part, imaginary part)
	//double s1h2[ND][ND][ND], s2h2[ND][ND][ND];		//Auxiliary arrays for s1 and s2

	//double ss1qrh[ND][ND][ND], ss1qih[ND][ND][ND];	//Fourier transform of s1*s1 (real part, imaginary part)
	//double ss2qrh[ND][ND][ND], ss2qih[ND][ND][ND];	//Fourier transform of s2*s2 (real part, imaginary part)
	//double s1s2qrh[ND][ND][ND], s1s2qih[ND][ND][ND];//Fourier transform of s1*s2 (real part, imaginary part)
	double ss1ss2;							//Work variables for correction of numerical errors in domain
	double ss1ss3;							//Work variables for correction of numerical errors in domain
	double ss2ss3;							//Work variables for correction of numerical errors in domain

	double a0_aa, a0_a, a0_c;				//Lattice constant of BaTiO3 (tetragonal)
	double al, temp, delt, vm0;				//Side length of calculation domain, temperature, time step, molar volume
	double time1max;						//Maximum calculation count (calculation end count)
	double b1;								//Length of one side of difference block

	double s1, s2, s3;						//x and y components of the moment of polarization
	double s1ip, s1im, s1jp, s1jm, s1kp, s1km;//left, right, top, bottom values of s1
	double s2ip, s2im, s2jp, s2jm, s2kp, s2km;//left, right, top, bottom values of s2
	double s3ip, s3im, s3jp, s3jm, s3kp, s3km;//left, right, top, bottom values of s3

	double s1k, s1k_chem, s1k_surf, s1k_str, s1k_ddi;	//potential for s1
	//double s1k_dd[ND][ND][ND];					//dipole-dipole interaction potential
	double s2k, s2k_chem, s2k_surf, s2k_str, s2k_ddi;	//potential for s2
	//double s2k_dd[ND][ND][ND];					//dipole-dipole interaction potential
	double s3k, s3k_chem, s3k_surf, s3k_str, s3k_ddi;	//potential for s2
	//double s3k_dd[ND][ND][ND];					//dipole-dipole interaction potential
	double smob1, smob2, smob3;					//Domain interface mobility
	double s1ddtt, s2ddtt, s3ddtt;				//left-hand side of the evolution equation

	double A1, A11, A12;					//parameters in chemical free energy
	double A1e, A11e, A12e;
	double A111, A112, A123;
	double A111e, A112e, A123e;
	double A1111, A1112, A1122, A1123;
	double A1111e, A1112e, A1122e, A1123e;
	double Tc0;

	double nx, ny, nz, alnn;				//Unit vector components of reciprocal space and reciprocal lattice vector length
	double kapaP, kapaPc;					//gradient energy factor

	double E1_ex;							//external electric field
	double E1_ex_0x, E1_ex_x; 				//external electric field (x)
	double E1_ex_0y, E1_ex_y; 				//external electric field (y)
	double E1_ex_0z, E1_ex_z; 				//external electric field (z)
	double ep00;							//dielectric constant
	double Add0, Add0c;						//Coefficients in dipole-dipole interaction calculations

	double t1, t2, t3, tQ, tR, tS, tT;		//Working variables for calculation of equilibrium moment of polarization
	int readff;

//****** Setting calculation conditions and material constants ****************************************
	printf("---------------------------------\n");
	printf("read parameters from parameters.txt\n");
	FILE *fp;
	char name[40], comment[172];
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
	time1max= int(data[2]);
	Nstep   = int(data[3]);
	temp    = data[4];
	al      = data[5];
	vm0     = data[6];
	smob1   = data[7];
	smob2   = data[8];
	smob3   = data[9];
	A1e     = data[10]; //Landau expansion form
	A11e    = data[11]; //Landau expansion form
	A12e    = data[12]; //Landau expansion form
	A111e   = data[13]; //Landau expansion form
	A112e   = data[14]; //Landau expansion form
	A123e   = data[15]; //Landau expansion form
	A1111e  = data[16]; //Landau expansion form
	A1112e  = data[17]; //Landau expansion form
	A1122e  = data[18]; //Landau expansion form
	A1123e  = data[19]; //Landau expansion form
	Tc0     = data[20];
	kapaPc  = data[21];
	E1_ex_x = data[22];
	E1_ex_y = data[23];
	E1_ex_z = data[24];
	ep00    = data[25];
	Add0c   = data[26];
	readff  = int(data[27]);
	printf("---------------------------------\n");
	//
	ig = int(log2(ND));
	nd=ND;
	ndm=ND-1;
	nd2=ND/2;
	//
	double *s1h     = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//polarization moment in x direction
	double *s2h     = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//polarization moment in y direction
	double *s3h     = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//polarization moment in y direction
	//
	double *xi      = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//array of real and imaginary parts of the Fourier transform
	double *xr      = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//array of real and imaginary parts of the Fourier transform
	//
	const int fftsize = ND;
	fftw_complex *in, *out; // in[i][0] for real, in[i][1] for imag.
	fftw_plan plan, iplan;
	in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsize*fftsize*fftsize);
	out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsize*fftsize*fftsize);
	//plan = fftw_plan_dft_1d(fftsize, in, out, FFTW_FORWARD, FFTW_ESTIMATE); //FFT
	//fftw_execute(plan); //FFT
	//iplan = fftw_plan_dft_1d(fftsize, in, out, FFTW_BACKWARD, FFTW_ESTIMATE); //IFFT
	//fftw_execute(iplan); //IFFT
	//
	//After calc.: fftw_destroy_plan(plan); fftw_free(in); fftw_free(out);
	//
	double *s1qrh   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s1 (real part)
	double *s1qih   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s1 (imaginary part)
	double *s2qrh   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s2 (real part)
	double *s2qih   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s2 (imaginary part)
	double *s3qrh   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s2 (real part)
	double *s3qih   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s2 (imaginary part)
	double *s1h2    = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//auxiliary array of s1
	double *s2h2    = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//auxiliary array of s2
	double *s3h2    = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//auxiliary array of s3
	//double *ss1qrh  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s1*s1 (real part)
	//double *ss1qih  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s1*s1 (imaginary part)
	//double *ss2qrh  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s2*s2 (real part)
	//double *ss2qih  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s2*s2 (imaginary part)
	//double *ss3qrh  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s2*s2 (real part)
	//double *ss3qih  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s2*s2 (imaginary part)
	//double *s1s2qrh = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s1*s2 (real part)
	//double *s1s2qih = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s1*s2 (imaginary part)
	//
	double *s1k_dd  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//dipole-dipole interaction potential
	double *s2k_dd  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//dipole-dipole interaction potential
	double *s3k_dd  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//dipole-dipole interaction potential
	//
	//printf("DELT(0.1)=  "); scanf(" %lf",&delt);//time step setting	//delt=0.1;

	time1=0.0;								//Initial calculation count setting
	//time1max=1.0+1.0e+07;					//Setting the maximum calculation count
	//Note that all time is dimensionless.

	//temp=298.0;							//temperature (K)

	//al=250.;								//Length of one side of calculation area (micro meter)
	b1=al*1.0E-06/nd;						//Length of one side of difference block (m)

	//a0_aa=0.4;  a0_a=0.3992;  a0_c=0.4036;//Lattice constant (nm)
	//vm0=a0_a*a0_a*a0_c*1.0e-27*6.02e+23;	//molar volume (molecule 1 mole)

	//smob1=1.; smob2=1.;					//Mobility in structural phase transitions (normalized and set to 1)

//--- parameter value in chemical free energy [see Table 4.7]----------------------
	//A1=4.124e+05*vm0/RR/temp;
	A1=A1e*vm0/RR/temp;
	//A11=-2.097e+08*vm0/RR/temp;
	A11=A11e*vm0/RR/temp;
	//A12=7.974e+08*vm0/RR/temp;
	A12=A12e*vm0/RR/temp;
	//A111=1.294e+09*vm0/RR/temp;
	A111=A111e*vm0/RR/temp;
	//A112=-1.950e+09*vm0/RR/temp;
	A112=A112e*vm0/RR/temp;
	//A123=-2.500e+09*vm0/RR/temp;
	A123=A123e*vm0/RR/temp;
	//A1111=3.863e+10*vm0/RR/temp;
	A1111=A1111e*vm0/RR/temp;
	//A1112=2.529e+10*vm0/RR/temp;
	A1112=A1112e*vm0/RR/temp;
	//A1122=1.637e+10*vm0/RR/temp;
	A1122=A1122e*vm0/RR/temp;
	//A1123=1.367e+10*vm0/RR/temp;
	A1123=A1123e*vm0/RR/temp;

	//Tc0=115.0+273.0;  // [K]

//--- Calculation of equilibrium moment of polarization ----------------------
	t1=3.0*A111/4.0/A1111;
	t2=A11/2.0/A1111;
	t3=A1*(temp-Tc0)/4.0/A1111;
	tQ=(3.0*t2-t1*t1)/9.0;
	tR=(9.0*t1*t2-27.0*t3-2.0*t1*t1*t1)/54.0;
	tS=pow((tR+sqrt(tQ*tQ*tQ+tR*tR)),(1.0/3.0));
	tT=pow((tR-sqrt(tQ*tQ*tQ+tR*tR)),(1.0/3.0));
	Peq=sqrt(fabs(tS+tR-t1/3.0));//equilibrium moment of polarization
	printf("Peq=  %f  \n", Peq);//Display of equilibrium moment of polarization

	//kapaP=1.0e-15/ep0*vm0/RR/temp/b1/b1;//gradient energy factor
	kapaP=kapaPc/ep0*vm0/RR/temp/b1/b1;//gradient energy factor

	//E1_ex_0=0.0;//External electric field (x direction)
	//E1_ex_0=2.0e+6*vm0/RR/temp;//External electric field (x direction)
	E1_ex_0x=E1_ex_x*vm0/RR/temp;//External electric field (x direction)
	E1_ex_0y=E1_ex_y*vm0/RR/temp;//External electric field (y direction)
	E1_ex_0z=E1_ex_z*vm0/RR/temp;//External electric field (z direction)

	//ep00=200.0;//dielectric constant

	//Add0=1.0/ep0/ep00*vm0/RR/temp;//Coefficient in dipole-dipole interaction calculation [see formula (4.52)]
	Add0=Add0c/ep0/ep00*vm0/RR/temp;//Coefficients in dipole-dipole interaction calculations

//*** Setting up sin and cos tables, bit reversal tables, and initial fields ***************

	if(readff == 0){
		ini000(s1h, s2h, s3h, ND);	//initial profile of the polarization moment at time 0
	} else {
		datin(s1h, s2h, s3h, ND);	//Input initial tissue field
	}

//**** Simulation start ******************************
//Nstep = 10;
 plan = fftw_plan_dft_3d(fftsize, fftsize, fftsize, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);	//For FFT
iplan = fftw_plan_dft_3d(fftsize, fftsize, fftsize, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);	//For IFFT
start: ;

	if((((int)(time1) % Nstep)==0)){datsave(s1h, s2h, s3h, ND);} //Save tissue data every fixed repeat count
	if((((int)(time1) % Nstep)==0)){datsave_paraview(s1h, s2h, s3h, ND);} //Save tissue data every fixed repeat count

//**** Fourier transform of s1 [equation (3.27)] ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=s1h[i][j][k]; xi[i][j][k]=0.;
				xr[i*ND*ND+j*ND+k]=s1h[i*ND*ND+j*ND+k];
				xi[i*ND*ND+j*ND+k]=0.;
				//
				in[i*ND*ND+j*ND+k][0] = xr[i*ND*ND+j*ND+k]; //For FFT (real)
				in[i*ND*ND+j*ND+k][1] = xi[i*ND*ND+j*ND+k]; //For FFT (imag.)
			}
		}
	}
	fftw_execute(plan); //For FFT
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k] = out[i*ND*ND+j*ND+k][0]; //For FFT (real)
				xi[i*ND*ND+j*ND+k] = out[i*ND*ND+j*ND+k][1]; //For FFT (imag.)
				//
				//s1qrh[i][j][k]=xr[i][j][k];  s1qih[i][j][k]=xi[i][j][k];
				s1qrh[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
				s1qih[i*ND*ND+j*ND+k]=xi[i*ND*ND+j*ND+k];
			}
		}
	}
	//s1qrh[0][0]=s1qih[0][0]=0.;
	s1qrh[0]=s1qih[0]=0.0;

//**** Fourier transform of s2 [equation (3.27)] ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=s2h[i][j][k]; xi[i][j][k]=0.;
				xr[i*ND*ND+j*ND+k]=s2h[i*ND*ND+j*ND+k];
				xi[i*ND*ND+j*ND+k]=0.0;
				//
				in[i*ND*ND+j*ND+k][0] = xr[i*ND*ND+j*ND+k]; //For FFT (real)
				in[i*ND*ND+j*ND+k][1] = xi[i*ND*ND+j*ND+k]; //For FFT (imag.)
			}
		}
	}
	fftw_execute(plan); //For FFT
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k] = out[i*ND*ND+j*ND+k][0]; //For FFT (real)
				xi[i*ND*ND+j*ND+k] = out[i*ND*ND+j*ND+k][1]; //For FFT (imag.)
				//
				//s2qrh[i][j][k]=xr[i][j][k];  s2qih[i][j][k]=xi[i][j][k];
				s2qrh[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
				s2qih[i*ND*ND+j*ND+k]=xi[i*ND*ND+j*ND+k];
			}
		}
	}
	//s2qrh[0][0]=s2qih[0][0]=0.;
	s2qrh[0]=s2qih[0]=0.0;

//**** Fourier transform of s3 [equation (3.27)] ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=s3h[i][j][k]; xi[i][j][k]=0.;
				xr[i*ND*ND+j*ND+k]=s3h[i*ND*ND+j*ND+k];
				xi[i*ND*ND+j*ND+k]=0.0;
				//
				in[i*ND*ND+j*ND+k][0] = xr[i*ND*ND+j*ND+k]; //For FFT (real)
				in[i*ND*ND+j*ND+k][1] = xi[i*ND*ND+j*ND+k]; //For FFT (imag.)
			}
		}
	}
	fftw_execute(plan); //For FFT
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k] = out[i*ND*ND+j*ND+k][0]; //For FFT (real)
				xi[i*ND*ND+j*ND+k] = out[i*ND*ND+j*ND+k][1]; //For FFT (imag.)
				//
				//s3qrh[i][j][k]=xr[i][j][k];  s3qih[i][j][k]=xi[i][j][k];
				s3qrh[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
				s3qih[i*ND*ND+j*ND+k]=xi[i*ND*ND+j*ND+k];
			}
		}
	}
	//s3qrh[0][0]=s3qih[0][0]=0.;
	s3qrh[0]=s3qih[0]=0.0;

//***** Computation of the dipole-dipole interaction of s1 [convolution integral in Eq. (3.33)] ***********************
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			for(k=0;k<=ndm;k++){
				if(k<=nd2-1){kk=k;}  if(k>=nd2){kk=k-nd;}
				//
				alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj+(double)kk*(double)kk);
				if(alnn==0.){alnn=1.;}
				nx=(double)ii/alnn;  ny=(double)jj/alnn;  nz=(double)kk/alnn;
				//xr[i][j][k]=s1qrh[i][j][k]*nx*nx+s2qrh[i][j][k]*nx*ny+s3qrh[i][j][k]*nx*nz;
				//xi[i][j][k]=s1qih[i][j][k]*nx*nx+s2qih[i][j][k]*nx*ny+s3qih[i][j][k]*nx*nz;
				xr[i*ND*ND+j*ND+k]=s1qrh[i*ND*ND+j*ND+k]*nx*nx+s2qrh[i*ND*ND+j*ND+k]*nx*ny+s3qrh[i*ND*ND+j*ND+k]*nx*nz;
				xi[i*ND*ND+j*ND+k]=s1qih[i*ND*ND+j*ND+k]*nx*nx+s2qih[i*ND*ND+j*ND+k]*nx*ny+s3qih[i*ND*ND+j*ND+k]*nx*nz;
				//
				in[i*ND*ND+j*ND+k][0] = xr[i*ND*ND+j*ND+k]; //For IFFT (real)
				in[i*ND*ND+j*ND+k][1] = xi[i*ND*ND+j*ND+k]; //For IFFT (imag.)
			}
		}
	}
	fftw_execute(iplan); //For IFFT
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ s1k_dd[i][j][k]=xr[i][j][k]; } }
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k] = out[i*ND*ND+j*ND+k][0]/(fftsize*fftsize*fftsize); //For IFFT (real)
				//
				s1k_dd[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//***** Computation of the dipole-dipole interaction of s2 [convolution integral in Eq. (3.33)] ************************
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			for(k=0;k<=ndm;k++){
				if(k<=nd2-1){kk=k;}  if(k>=nd2){kk=k-nd;}
				//
				alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj+(double)kk*(double)kk);
				if(alnn==0.){alnn=1.;}
				nx=(double)ii/alnn;  ny=(double)jj/alnn;  nz=(double)kk/alnn;
				//xr[i][j][k]=s1qrh[i][j][k]*ny*nx+s2qrh[i][j][k]*ny*ny+s3qrh[i][j][k]*ny*nz;
				//xi[i][j][k]=s1qih[i][j][k]*ny*nx+s2qih[i][j][k]*ny*ny+s3qih[i][j][k]*ny*nz;
				xr[i*ND*ND+j*ND+k]=s1qrh[i*ND*ND+j*ND+k]*ny*nx+s2qrh[i*ND*ND+j*ND+k]*ny*ny+s3qrh[i*ND*ND+j*ND+k]*ny*nz;
				xi[i*ND*ND+j*ND+k]=s1qih[i*ND*ND+j*ND+k]*ny*nx+s2qih[i*ND*ND+j*ND+k]*ny*ny+s3qih[i*ND*ND+j*ND+k]*ny*nz;
				//
				in[i*ND*ND+j*ND+k][0] = xr[i*ND*ND+j*ND+k]; //For IFFT (real)
				in[i*ND*ND+j*ND+k][1] = xi[i*ND*ND+j*ND+k]; //For IFFT (imag.)
			}
		}
	}
	fftw_execute(iplan); //For IFFT
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ s2k_dd[i][j][k]=xr[i][j][k]; } }
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k] = out[i*ND*ND+j*ND+k][0]/(fftsize*fftsize*fftsize); //For IFFT (real)
				//
				s2k_dd[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//***** Computation of the dipole-dipole interaction of s3 [convolution integral in Eq. (3.33)] ************************
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			for(k=0;k<=ndm;k++){
				if(k<=nd2-1){kk=k;}  if(k>=nd2){kk=k-nd;}
				//
				alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj+(double)kk*(double)kk);
				if(alnn==0.){alnn=1.;}
				nx=(double)ii/alnn;  ny=(double)jj/alnn;  nz=(double)kk/alnn;
				//xr[i][j][k]=s1qrh[i][j][k]*nz*nx+s2qrh[i][j][k]*nz*ny+s3qrh[i][j][k]*nz*nz;
				//xi[i][j][k]=s1qih[i][j][k]*nz*nx+s2qih[i][j][k]*nz*ny+s3qih[i][j][k]*nz*nz;
				xr[i*ND*ND+j*ND+k]=s1qrh[i*ND*ND+j*ND+k]*nz*nx+s2qrh[i*ND*ND+j*ND+k]*nz*ny+s3qrh[i*ND*ND+j*ND+k]*nz*nz;
				xi[i*ND*ND+j*ND+k]=s1qih[i*ND*ND+j*ND+k]*nz*nx+s2qih[i*ND*ND+j*ND+k]*nz*ny+s3qih[i*ND*ND+j*ND+k]*nz*nz;
				//
				in[i*ND*ND+j*ND+k][0] = xr[i*ND*ND+j*ND+k]; //For IFFT (real)
				in[i*ND*ND+j*ND+k][1] = xi[i*ND*ND+j*ND+k]; //For IFFT (imag.)
			}
		}
	}
	fftw_execute(iplan); //For IFFT
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ s3k_dd[i][j][k]=xr[i][j][k]; } }
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k] = out[i*ND*ND+j*ND+k][0]/(fftsize*fftsize*fftsize); //For IFFT (real)
				//
				s3k_dd[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}


//****** Numerical calculation of nonlinear evolution equations  ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//
				ip=i+1;  im=i-1;  jp=j+1;  jm=j-1;  kp=k+1;  km=k-1;
				if(i==ndm){ip=0;}  if(i==0){im=ndm;}//periodic boundary conditions
				if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
				if(k==ndm){kp=0;}  if(k==0){km=ndm;}
				//s1=s1h[i][j][k];s1ip=s1h[ip][j][k]; s1im=s1h[im][j][k]; s1jp=s1h[i][jp][k]; s1jm=s1h[i][jm][k]; s1jp=s1h[i][j][kp]; s1jm=s1h[i][j][km];
				//s2=s2h[i][j][k];s2ip=s2h[ip][j][k]; s2im=s2h[im][j][k]; s2jp=s2h[i][jp][k]; s2jm=s2h[i][jm][k]; s2jp=s2h[i][j][kp]; s3jm=s2h[i][j][km];
				//s3=s3h[i][j][k];s3ip=s3h[ip][j][k]; s3im=s3h[im][j][k]; s3jp=s3h[i][jp][k]; s3jm=s3h[i][jm][k]; s3jp=s3h[i][j][kp]; s4jm=s3h[i][j][km];
				s1=s1h[i*ND*ND+j*ND+k]; s1ip=s1h[ip*ND*ND+j*ND+k]; s1im=s1h[im*ND*ND+j*ND+k];
										s1jp=s1h[i*ND*ND+jp*ND+k]; s1jm=s1h[i*ND*ND+jm*ND+k];
										s1kp=s1h[i*ND*ND+j*ND+kp]; s1km=s1h[i*ND*ND+j*ND+km];
				s2=s2h[i*ND*ND+j*ND+k]; s2ip=s2h[ip*ND*ND+j*ND+k]; s2im=s2h[im*ND*ND+j*ND+k];
										s2jp=s2h[i*ND*ND+jp*ND+k]; s2jm=s2h[i*ND*ND+jm*ND+k];
										s2kp=s2h[i*ND*ND+j*ND+kp]; s2km=s2h[i*ND*ND+j*ND+km];
				s3=s3h[i*ND*ND+j*ND+k]; s3ip=s3h[ip*ND*ND+j*ND+k]; s3im=s3h[im*ND*ND+j*ND+k];
										s3jp=s3h[i*ND*ND+jp*ND+k]; s3jm=s3h[i*ND*ND+jm*ND+k];
										s3kp=s3h[i*ND*ND+j*ND+kp]; s3km=s3h[i*ND*ND+j*ND+km];

				//Computation of gradient potential
				s1k_surf=-kapaP*(s1ip+s1im + s1jp+s1jm + s1kp+s1km -6.0*s1);
				s2k_surf=-kapaP*(s2ip+s2im + s2jp+s2jm + s2kp+s2km -6.0*s2);
				s3k_surf=-kapaP*(s3ip+s3im + s3jp+s3jm + s3kp+s3km -6.0*s3);

				//Calculation of chemical potential [equation (4.55)]
					//Landau expansion form (BaTiO3) (structure phase transition) (cubic-tetragonal transition)
					//Gc = a1*(P1^2 + P2^2 + P3^2) + a11*(P1^4 + P2^4 + P3^4)
					//+ a12*( (P1^2 * P2^2) + (P2^2 * P3^2) + (P3^2 * P1^2) ) + a111*(P1^6 + P2^6 + P3^6)
					//+ a112*( (P1^2 * (P2^4 * P3^4)) + (P2^2 * (P1^4 * P3^4)) + (P3^2 * (P1^4 * P2^4))  )
					//+ a123*( P1^2 * P2^2 * P3^3 ) + a1111*(P1^8 + P2^8 + P3^8)
					//+ a1112*( (P1^6 * (P2^2 + P3^2)) + (P2^6 * (P1^2 + P3^2)) + (P3^6 * (P1^2 + P3^2)) )
					//+ a1122*( (P1^4 * P2^4) + (P2^4 * P3^4) + (P3^4 * P1^4) )
					//+ a1123*( (P1^4 * P2^2 * P3^2) + (P1^2 * P2^4 * P3^2) + (P1^2 * P2^2 * P3^4) )
				//dGc/dP1, P1=s1, P2=s2, P3=s3, A<-a
				s1k_chem=2.0*A1*(temp-Tc0)*s1 + 4.0*A11*(s1*s1*s1) + 2.0*A12*s1*(s2*s2 + s3*s3)
								+6.0*A111*pow(s1,5.0) + 2.0*A112*s1*( (pow(s2,4.0)+pow(s3,4.0)) + 2.0*(s1*s1)*(s2*s2 + 4.0*s3*s3) )
								+2.0*A123*s1*(s2*s2)*(s3*s3) + 8.0*A1111*pow(s1,7.0)
								+2.0*A1112*s1*( 3.0*pow(s1,4.0)*(s2*s2 + s3*s3) + (pow(s2,6.0)+pow(s3,6.0)) )
								+4.0*A1122*pow(s1,3.0)*( pow(s2,4.0) + pow(s3,4.0) )
								+2.0*A1123*s1*( 2.0*(s1*s1)*(s2*s2)*(s3*s3) + pow(s2,4.0)*(s3*s3) + (s2*s2)*pow(s3,4.0) );

				//dGc/dP2, P1=s1, P2=s2, P3=s3, A<-a
				s2k_chem=2.0*A1*(temp-Tc0)*s2 + 4.0*A11*(s2*s2*s2) + 2.0*A12*s2*(s1*s1 + s3*s3)
								+6.0*A111*pow(s2,5.0) + 2.0*A112*s2*( (pow(s1,4.0)+pow(s3,4.0)) + 2.0*(s2*s2)*(s1*s1 + 4.0*s3*s3) )
								+2.0*A123*s2*(s1*s1)*(s3*s3) +8.0*A1111*pow(s2,7.0)
								+2.0*A1112*s2*( 3.0*pow(s2,4.0)*(s1*s1 + s3*s3) + (pow(s1,6.0)+pow(s3,6.0)) )
								+4.0*A1122*pow(s2,3.0)*( pow(s1,4.0) + pow(s3,4.0) )
								+2.0*A1123*s2*( 2.0*(s1*s1)*(s2*s2)*(s3*s3) + pow(s1,4.0)*(s3*s3) + (s1*s1)*pow(s3,4.0) );

				//dGc/dP3, P1=s1, P2=s2, P3=s3, A<-a
				s3k_chem=2.0*A1*(temp-Tc0)*s3 + 4.0*A11*(s3*s3*s3) + 2.0*A12*s3*(s1*s1 + s2*s2)
								+6.0*A111*pow(s3,5.0) + 2.0*A112*s3*( (pow(s1,4.0)+pow(s2,4.0)) + 2.0*(s3*s3)*(s1*s1 + 4.0*s2*s2) )
								+2.0*A123*s3*(s2*s2)*(s3*s3) +8.0*A1111*pow(s3,7.0)
								+2.0*A1112*s3*( 3.0*pow(s3,4.0)*(s1*s1 + s2*s2) + (pow(s1,6.0)+pow(s2,6.0)) )
								+4.0*A1122*pow(s3,3.0)*( pow(s1,4.0) + pow(s2,4.0) )
								+2.0*A1123*s3*( 2.0*(s1*s1)*(s2*s2)*(s3*s3) + pow(s1,4.0)*(s2*s2) + (s1*s1)*pow(s2,4.0) );

				//Calculation of dipole-dipole potential [Formula (4.56) right side first term]
				//s1k_ddi=Add0*s1k_dd[i][j][k];
				//s2k_ddi=Add0*s2k_dd[i][j][k];
				//s3k_ddi=Add0*s3k_dd[i][j][k];
				s1k_ddi=Add0*s1k_dd[i*ND*ND+j*ND+k];
				s2k_ddi=Add0*s2k_dd[i*ND*ND+j*ND+k];
				s3k_ddi=Add0*s3k_dd[i*ND*ND+j*ND+k];

				//Computation of total potential
				s1k=s1k_chem+s1k_surf+s1k_ddi-E1_ex_0x;
				s2k=s2k_chem+s2k_surf+s2k_ddi-E1_ex_0y;
				s3k=s3k_chem+s3k_surf+s3k_ddi-E1_ex_0z;

				//Evolution equation [equation (4.57)] and time evolution of s-field (explicit method) calculation
				//s1ddtt=-smob1*s1k;  s1h2[i][j][k]=s1h[i][j][k]+s1ddtt*delt;
				//s2ddtt=-smob2*s2k;  s2h2[i][j][k]=s2h[i][j][k]+s2ddtt*delt;
				//s3ddtt=-smob3*s3k;  s3h2[i][j][k]=s3h[i][j][k]+s3ddtt*delt;
				s1ddtt=-smob1*s1k;  s1h2[i*ND*ND+j*ND+k]=s1h[i*ND*ND+j*ND+k]+s1ddtt*delt;
				s2ddtt=-smob2*s2k;  s2h2[i*ND*ND+j*ND+k]=s2h[i*ND*ND+j*ND+k]+s2ddtt*delt;
				s3ddtt=-smob3*s3k;  s3h2[i*ND*ND+j*ND+k]=s3h[i*ND*ND+j*ND+k]+s3ddtt*delt;
			}//k
		}//j
	}//i

	//Move the s-field to the main array and correct the numerical errors in the domain
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//s1h[i][j][k]=s1h2[i][j][k]; s2h[i][j][k]=s2h2[i][j][k]; s3h[i][j][k]=s3h2[i][j][k];
				s1h[i*ND*ND+j*ND+k]=s1h2[i*ND*ND+j*ND+k];
				s2h[i*ND*ND+j*ND+k]=s2h2[i*ND*ND+j*ND+k];
				s3h[i*ND*ND+j*ND+k]=s3h2[i*ND*ND+j*ND+k];
				//
				//if(s1h[i][j][k]>=Peq){s1h[i][j][k]=Peq;}  if(s1h[i][j][k]<=-Peq){s1h[i][j][k]=-Peq;}
				if(s1h[i*ND*ND+j*ND+k]>= Peq){s1h[i*ND*ND+j*ND+k]= Peq;}
				if(s1h[i*ND*ND+j*ND+k]<=-Peq){s1h[i*ND*ND+j*ND+k]=-Peq;}
				//if(s2h[i][j][k]>=Peq){s2h[i][j][k]=Peq;}  if(s2h[i][j][k]<=-Peq){s2h[i][j][k]=-Peq;}
				if(s2h[i*ND*ND+j*ND+k]>= Peq){s2h[i*ND*ND+j*ND+k]= Peq;}
				if(s2h[i*ND*ND+j*ND+k]<=-Peq){s2h[i*ND*ND+j*ND+k]=-Peq;}
				//if(s3h[i][j][k]>=Peq){s3h[i][j][k]=Peq;}  if(s3h[i][j][k]<=-Peq){s3h[i][j][k]=-Peq;}
				if(s3h[i*ND*ND+j*ND+k]>= Peq){s3h[i*ND*ND+j*ND+k]= Peq;}
				if(s3h[i*ND*ND+j*ND+k]<=-Peq){s3h[i*ND*ND+j*ND+k]=-Peq;}
				//
				//ss1ss2=sqrt(s1h[i][j][k]*s1h[i][j][k]+s2h[i][j][k]*s2h[i][j][k]);
				ss1ss2=sqrt(s1h[i*ND*ND+j*ND+k]*s1h[i*ND*ND+j*ND+k] + s2h[i*ND*ND+j*ND+k]*s2h[i*ND*ND+j*ND+k]);
				//if(ss1ss2>=Peq){s1h[i][j][k]=s1h[i][j][k]/ss1ss2*Peq; s2h[i][j][k]=s2h[i][j][k]/ss1ss2*Peq;}
				if(ss1ss2>=Peq){
					s1h[i*ND*ND+j*ND+k]=s1h[i*ND*ND+j*ND+k]/ss1ss2*Peq;
					s2h[i*ND*ND+j*ND+k]=s2h[i*ND*ND+j*ND+k]/ss1ss2*Peq;
				}
				//
				//ss1ss3=sqrt(s1h[i][j][k]*s1h[i][j][k]+s3h[i][j][k]*s3h[i][j][k]);
				ss1ss3=sqrt(s1h[i*ND*ND+j*ND+k]*s1h[i*ND*ND+j*ND+k] + s3h[i*ND*ND+j*ND+k]*s3h[i*ND*ND+j*ND+k]);
				//if(ss1ss3>=Peq){s1h[i][j][k]=s1h[i][j][k]/ss1ss3*Peq; s3h[i][j][k]=s3h[i][j][k]/ss1ss3*Peq;}
				if(ss1ss3>=Peq){
					s1h[i*ND*ND+j*ND+k]=s1h[i*ND*ND+j*ND+k]/ss1ss3*Peq;
					s3h[i*ND*ND+j*ND+k]=s3h[i*ND*ND+j*ND+k]/ss1ss3*Peq;
				}
				//
				//ss2ss3=sqrt(s2h[i][j][k]*s2h[i][j][k]+s3h[i][j][k]*s3h[i][j][k]);
				ss2ss3=sqrt(s2h[i*ND*ND+j*ND+k]*s2h[i*ND*ND+j*ND+k] + s3h[i*ND*ND+j*ND+k]*s3h[i*ND*ND+j*ND+k]);
				//if(ss2ss3>=Peq){s2h[i][j][k]=s2h[i][j][k]/ss2ss3*Peq; s3h[i][j][k]=s3h[i][j][k]/ss2ss3*Peq;}
				if(ss2ss3>=Peq){
					s2h[i*ND*ND+j*ND+k]=s2h[i*ND*ND+j*ND+k]/ss2ss3*Peq;
					s3h[i*ND*ND+j*ND+k]=s3h[i*ND*ND+j*ND+k]/ss2ss3*Peq;
				}
			}//k
		}//j
	}//i

	//Advance the time, and if it is not the end time, go to the start
	//if(keypress()){return 0;}//Waiting for key
	time1=time1+1.0;								//Add calculation count
	if(time1<time1max){goto start;}	//Determining if the maximum count has been reached
	printf("Finished \n");
	
end:;
	if(plan) fftw_destroy_plan(plan);		//For FFT
	if(iplan) fftw_destroy_plan(iplan);	//For IFFT
	fftw_free(in); fftw_free(out);		// For FFT and IFFT
	std::exit(0);
	//return 0;
}

//************ Initial field setup subroutine *************
void ini000(double *s1h, double *s2h, double *s3h, int ND)
{
	int i, j, k;	//integer
 	//srand(time(NULL)); // random number initialization
	int nd=ND, ndm=ND-1, nd2=ND/2;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//s1h[i][j][k]=0.01*(2.0*DRND(1)-1.0);//Set the field with random numbers up to }1%
				//s2h[i][j][k]=0.01*(2.0*DRND(1)-1.0);
				//s3h[i][j][k]=0.01*(2.0*DRND(1)-1.0);
				s1h[i*ND*ND+j*ND+k]=0.01*(2.0*DRND(1)-1.0);//Set the field with random numbers up to }1%
				s2h[i*ND*ND+j*ND+k]=0.01*(2.0*DRND(1)-1.0);
				s3h[i*ND*ND+j*ND+k]=0.01*(2.0*DRND(1)-1.0);
			}
		}
	}
}

//************ data save subroutine *******************************
void datsave(double *s1h, double *s2h, double *s3h, int ND)
{
	FILE		*stream;	//Stream pointer setting
	int 		i, j, k;			//integer
	int nd=ND, ndm=ND-1, nd2=ND/2;

	stream = fopen("test.dat", "a");	//Open the write destination file for appending
	fprintf(stream, "%e \n", time1);	//Saving calculation counts
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//fprintf(stream, "%e  %e  ", s1h[i][j][k], s2h[i][j][k]);	//Field data storage
				fprintf(stream, "%e  %e  %e  ", s1h[i*ND*ND+j*ND+k], s2h[i*ND*ND+j*ND+k], s3h[i*ND*ND+j*ND+k]);	//Field data storage
			}
		}
	}
	fprintf(stream, "\n");	//writing a newline
	fclose(stream);					//close file
}

void datsave_paraview(double *s1h, double *s2h, double *s3h, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j, k;
	int nd=ND, ndm=ND-1, nd2=ND/2;
	
	iout = iout + 1;
	printf("dip_result%06d.vtk \n",iout);
	sprintf(fName,"dip_result%06d.vtk",iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),(ndm+1));
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*(ndm+1)));
	fprintf(fp,"SCALARS Phase_field_1 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", s1h[i][j][k]); //Field data
				fprintf(fp,"%10.6f\n", s1h[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Phase_field_2 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", s2h[i][j][k]); //Field data
				fprintf(fp,"%10.6f\n", s2h[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Phase_field_3 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", s3h[i][j][k]); //Field data
				fprintf(fp,"%10.6f\n", s3h[i*ND*ND+j*ND+k]);
			}
		}
	}
	fclose(fp);
}
//************ data loading subroutine *******************************
void datin(double *s1h, double *s2h, double *s3h, int ND)
{
	FILE		*datin0;	//Stream pointer setting
	int 		i, j, k;			//integer
	double time_ext;	//calculation count number
	int nd=ND, ndm=ND-1, nd2=ND/2;

	datin0 = fopen("test.dat", "r");	//open source file
	fscanf(datin0, "%lf", &time_ext);		//Read calculation count
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//fscanf(datin0, "%lf  %lf  ", &s1h[i][j][k], &s2h[i][j][k]);	//Field data read
				fscanf(datin0, "%lf  %lf  %lf  ", &s1h[i*ND*ND+j*ND+k], &s2h[i*ND*ND+j*ND+k], &s3h[i*ND*ND+j*ND+k]);	//Field data read
			}
		}
	}
	fclose(datin0);//close file
}