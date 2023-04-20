#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//Random number function setting

//#define ND 128	//The number of divisions per side of the computational domain in difference computation
					// (a power of 2 due to the use of fast Fourier transform)
//#define IG 7		// 2^IG=ND

	//int nd=ND, ndm=ND-1; 			//Define the number of difference divisions
									// (number of difference blocks) on one side of the calculation area, ND-1
	//int nd2=ND/2;				 	/Define /ND/2: used in Fast Fourier Transform
	//int ig=IG;					//2^ig=ND
	double PI=3.14159;				//Pi
	double RR=8.3145;				//gas constant
	double time1;					//Number of calculation counts (proportional to time)
	double ep0=8.8541878e-12; 		//Vacuum permittivity (F/m)
	double Peq;						//equilibrium value of the moment of polarization
	int iout;

	//double s1h[ND][ND], s2h[ND][ND];//polarization moment in x direction, polarization moment in y direction

	double qs;						//Distinction between Fourier transform (qs:-1) and inverse Fourier transform (qs:1)
	//double xi[ND][ND], xr[ND][ND], xif[ND], xrf[ND];//array of real and imaginary parts of the Fourier transform
	//double s[ND],c[ND];			//sin and cos table
	//int ik[ND];					//bit reversal table

	void ini000_s12(double *s1h, double *s2h, int ND);	//initial profile of the polarization moment at time 0
	void table(double *s, double *c, int *ik, int ND, int ig);	//Creation subroutine for sin and cos table and bit reversal table
	void fft(double *xrf, double *xif, double *s, double *c, int ND, int ig);	//One-dimensional fast Fourier transform
	void rcfft(double *xrf, double *xif, double *xr, double *xi, int *ik, double *s, double *c, int ND, int ig);//2D Fast Fourier Transform
	void datsave(double *s1h, double *s2h, int ND);	//data save subroutine
	void datsave_paraview(double *s1h, double *s2h, int ND);	//data save subroutine
	void datin(double *s1h, double *s2h, int ND);	//Initial data loading

//******* main program ******************************************
int main(void)
{
	int ND;
	int nd, ndm, nd2, ig;
	
	int   i, j, k, l, ii, jj, Nstep;		//integer
	int   ip, im, jp, jm;					//integer

	//double s1qrh[ND][ND], s1qih[ND][ND];	//Fourier transform of s1 (real part, imaginary part)
	//double s2qrh[ND][ND], s2qih[ND][ND];	//Fourier transform of s2 (real part, imaginary part)
	//double s1h2[ND][ND], s2h2[ND][ND];	//Auxiliary arrays for s1 and s2

	//double ss1qrh[ND][ND], ss1qih[ND][ND];//Fourier transform of s1*s1 (real part, imaginary part)
	//double ss2qrh[ND][ND], ss2qih[ND][ND];//Fourier transform of s2*s2 (real part, imaginary part)
	//double s1s2qrh[ND][ND], s1s2qih[ND][ND];//Fourier transform of s1*s2 (real part, imaginary part)
	double ss1ss2;							//Work variables for correction of numerical errors in domain

	double a0_aa, a0_a, a0_c;				//Lattice constant of BaTiO3 (tetragonal)
	double al, temp, delt, vm0;				//Side length of calculation domain, temperature, time step, molar volume
	double time1max;						//Maximum calculation count (calculation end count)
	double b1;								//Length of one side of difference block

	double s1, s2;							//x and y components of the moment of polarization
	double s1ip, s1im, s1jp, s1jm;			//left, right, top, bottom values of s1
	double s2ip, s2im, s2jp, s2jm;			//left, right, top, bottom values of s2

	double s1k, s1k_chem, s1k_surf, s1k_str, s1k_ddi;	//potential for s1
	//double s1k_dd[ND][ND];				//dipole-dipole interaction potential
	double s2k, s2k_chem, s2k_surf, s2k_str, s2k_ddi;	//potential for s2
	//double s2k_dd[ND][ND];				//dipole-dipole interaction potential
	double smob1, smob2;					//Domain interface mobility
	double s1ddtt, s2ddtt;					//left-hand side of the evolution equation

	double A1, A11, A12;					//parameters in chemical free energy
	double A1e, A11e, A12e;
	double A111, A112, A123;
	double A111e, A112e, A123e;
	double A1111, A1112, A1122, A1123;
	double A1111e, A1112e, A1122e, A1123e;
	double Tc0;

	double nx, ny, alnn;					//Unit vector components of reciprocal space and reciprocal lattice vector length
	double kapaP, kapaPc;					//gradient energy factor

	double E1_ex, E1_ex_0, E1_ex_0c; 		//external electric field
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
	A1e     = data[9];
	A11e    = data[10];
	A12e    = data[11];
	A111e   = data[12];
	A112e   = data[13];
	A123e   = data[14];
	A1111e  = data[15];
	A1112e  = data[16];
	A1122e  = data[17];
	A1123e  = data[18];
	Tc0     = data[19];
	kapaPc  = data[20];
	E1_ex_0c= data[21];
	ep00    = data[22];
	Add0c   = data[23];
	readff  = int(data[24]);
	printf("---------------------------------\n");
	//
	ig = int(log2(ND));
	nd=ND;
	ndm=ND-1;
	nd2=ND/2;
	//
	double *s1h     = (double *)malloc(sizeof(double)*( ND*ND ));	//polarization moment in x direction
	double *s2h     = (double *)malloc(sizeof(double)*( ND*ND ));	//polarization moment in y direction
	//
	double *xi      = (double *)malloc(sizeof(double)*( ND*ND ));//array of real and imaginary parts of the Fourier transform
	double *xr      = (double *)malloc(sizeof(double)*( ND*ND ));//array of real and imaginary parts of the Fourier transform
	//
	double *xif     = (double *)malloc(sizeof(double)*( ND ));//array of real and imaginary parts of the Fourier transform
	double *xrf     = (double *)malloc(sizeof(double)*( ND ));//array of real and imaginary parts of the Fourier transform
	//
	double *s       = (double *)malloc(sizeof(double)*( ND ));//sin table
	double *c       = (double *)malloc(sizeof(double)*( ND ));//cos table
	//
	int *ik         = (int *)malloc(sizeof(int)*( ND ));//bit reversal table
	//
	double *s1qrh   = (double *)malloc(sizeof(double)*( ND*ND ));	//Fourier transform of s1 (real part)
	double *s1qih   = (double *)malloc(sizeof(double)*( ND*ND ));	//Fourier transform of s1 (imaginary part)
	double *s2qrh   = (double *)malloc(sizeof(double)*( ND*ND ));	//Fourier transform of s2 (real part)
	double *s2qih   = (double *)malloc(sizeof(double)*( ND*ND ));	//Fourier transform of s2 (imaginary part)
	double *s1h2    = (double *)malloc(sizeof(double)*( ND*ND ));	//auxiliary array of s1
	double *s2h2    = (double *)malloc(sizeof(double)*( ND*ND ));	//auxiliary array of s2
	double *ss1qrh  = (double *)malloc(sizeof(double)*( ND*ND ));	//Fourier transform of s1*s1 (real part)
	double *ss1qih  = (double *)malloc(sizeof(double)*( ND*ND ));	//Fourier transform of s1*s1 (imaginary part)
	double *ss2qrh  = (double *)malloc(sizeof(double)*( ND*ND ));	//Fourier transform of s2*s2 (real part)
	double *ss2qih  = (double *)malloc(sizeof(double)*( ND*ND ));	//Fourier transform of s2*s2 (imaginary part)
	double *s1s2qrh = (double *)malloc(sizeof(double)*( ND*ND ));	//Fourier transform of s1*s2 (real part)
	double *s1s2qih = (double *)malloc(sizeof(double)*( ND*ND ));	//Fourier transform of s1*s2 (imaginary part)
	//
	double *s1k_dd  = (double *)malloc(sizeof(double)*( ND*ND ));	//dipole-dipole interaction potential
	double *s2k_dd  = (double *)malloc(sizeof(double)*( ND*ND ));	//dipole-dipole interaction potential
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

	//smob1=1.; smob2=1.;					//Mobility at structural phase transitions (normalized and set to 1Åj

//--- Parameter value in chemical free energy [see Table 4.7]----------------------
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
	t1=3.*A111/4./A1111;
	t2=A11/2./A1111;
	t3=A1*(temp-Tc0)/4./A1111;
	tQ=(3.*t2-t1*t1)/9.;
	tR=(9.*t1*t2-27.*t3-2.*t1*t1*t1)/54.;
	tS=pow((tR+sqrt(tQ*tQ*tQ+tR*tR)),(1./3.));
	tT=pow((tR-sqrt(tQ*tQ*tQ+tR*tR)),(1./3.));
	Peq=sqrt(fabs(tS+tR-t1/3.));//equilibrium moment of polarization
	printf("Peq=  %f  \n", Peq);//Display of equilibrium moment of polarization

	//kapaP=1.0e-15/ep0*vm0/RR/temp/b1/b1;//gradient energy factor
	kapaP=kapaPc/ep0*vm0/RR/temp/b1/b1;//gradient energy factor
	
	//E1_ex_0=0.;//External electric field (x direction)
	//E1_ex_0=2.0e+6*vm0/RR/temp;//External electric field (x direction)
	E1_ex_0=E1_ex_0c*vm0/RR/temp;//External electric field (x direction)

	//ep00=200.;//dielectric constant

	//Add0=1.0/ep0/ep00*vm0/RR/temp;//Coefficient in dipole-dipole interaction calculation [see formula (4.7.8)]
	Add0=Add0c/ep0/ep00*vm0/RR/temp;//Coefficients in dipole-dipole interaction calculations

//*** Setting up sin and cos tables, bit reversal tables, and initial fields ***************

	table(s, c, ik, ND, ig);	//Setting up sin and cos tables and bit reversal tables
	if(readff == 0){
		ini000_s12(s1h, s2h, ND);	//initial profile of the polarization moment at time 0
	} else {
		datin(s1h, s2h, ND);		//Input initial tissue field
	}

//**** Simulation start ******************************
//Nstep = 10;
start: ;

	//if(time1<=100.){Nstep=10;} else{Nstep=100;}//Changing the time interval for saving data
	if((((int)(time1) % Nstep)==0)){datsave(s1h, s2h, ND);} //Save tissue data every fixed repeat count
	if((((int)(time1) % Nstep)==0)){datsave_paraview(s1h, s2h, ND);} //Save tissue data every fixed repeat count
	//if((int)(time1)==2000){datsave();} //Saving data at specific calculation counts

//**** Fourier transform of s1 [equation (3.37)] ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=s1h[i][j]; xi[i][j]=0.;
			xr[i*ND+j]=s1h[i*ND+j];
			xi[i*ND+j]=0.;
		}
	}
	qs=-1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s1qrh[i][j]=xr[i][j];  s1qih[i][j]=xi[i][j];
			s1qrh[i*ND+j]=xr[i*ND+j];
			s1qih[i*ND+j]=xi[i*ND+j];
		}
	}
	//s1qrh[0][0]=s1qih[0][0]=0.;
	s1qrh[0]=s1qih[0]=0.0;

//**** Fourier transform of s2 [equation (3.37)] ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=s2h[i][j]; xi[i][j]=0.;
			xr[i*ND+j]=s2h[i*ND+j];
			xi[i*ND+j]=0.0;
		}
	}
	qs=-1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s2qrh[i][j]=xr[i][j];  s2qih[i][j]=xi[i][j];
			s2qrh[i*ND+j]=xr[i*ND+j];
			s2qih[i*ND+j]=xi[i*ND+j];
		}
	}
	//s2qrh[0][0]=s2qih[0][0]=0.;
	s2qrh[0]=s2qih[0]=0.0;

//***** Computation of the dipole-dipole interaction of s1 [convolution integral in Eq. (3.33)] ***********************
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj);
			if(alnn==0.){alnn=1.;}
			nx=(double)ii/alnn;  ny=(double)jj/alnn;
			//xr[i][j]=s1qrh[i][j]*nx*nx+s2qrh[i][j]*nx*ny;
			//xi[i][j]=s1qih[i][j]*nx*nx+s2qih[i][j]*nx*ny;
			xr[i*ND+j]=s1qrh[i*ND+j]*nx*nx+s2qrh[i*ND+j]*nx*ny;
			xi[i*ND+j]=s1qih[i*ND+j]*nx*nx+s2qih[i*ND+j]*nx*ny;
		}
	}
	qs=1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ s1k_dd[i][j]=xr[i][j]; } }
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ s1k_dd[i*ND+j]=xr[i*ND+j]; } }

//***** Computation of the dipole-dipole interaction of s2 [convolution integral in Eq. (3.33)] ************************
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj);
			if(alnn==0.){alnn=1.;}
			nx=(double)ii/alnn;  ny=(double)jj/alnn;
			//xr[i][j]=s1qrh[i][j]*ny*nx+s2qrh[i][j]*ny*ny;
			//xi[i][j]=s1qih[i][j]*ny*nx+s2qih[i][j]*ny*ny;
			xr[i*ND+j]=s1qrh[i*ND+j]*ny*nx+s2qrh[i*ND+j]*ny*ny;
			xi[i*ND+j]=s1qih[i*ND+j]*ny*nx+s2qih[i*ND+j]*ny*ny;
		}
	}
	qs=1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ s2k_dd[i][j]=xr[i][j]; } }
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ s2k_dd[i*ND+j]=xr[i*ND+j]; } }

//****** Numerical calculation of nonlinear evolution equations ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1;  im=i-1;  jp=j+1;  jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//periodic boundary conditions
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
			//s1=s1h[i][j]; s1ip=s1h[ip][j]; s1im=s1h[im][j]; s1jp=s1h[i][jp]; s1jm=s1h[i][jm];
			//s2=s2h[i][j]; s2ip=s2h[ip][j]; s2im=s2h[im][j]; s2jp=s2h[i][jp]; s2jm=s2h[i][jm];
			s1=s1h[i*ND+j]; s1ip=s1h[ip*ND+j]; s1im=s1h[im*ND+j]; s1jp=s1h[i*ND+jp]; s1jm=s1h[i*ND+jm];
			s2=s2h[i*ND+j]; s2ip=s2h[ip*ND+j]; s2im=s2h[im*ND+j]; s2jp=s2h[i*ND+jp]; s2jm=s2h[i*ND+jm];

			//Computation of gradient potential
			s1k_surf=-kapaP*(s1ip+s1im+s1jp+s1jm-4.*s1);
			s2k_surf=-kapaP*(s2ip+s2im+s2jp+s2jm-4.*s2);

			//Calculation of chemical potential [equation (4.55)]
			s1k_chem=2.0*A1*(temp-Tc0)*s1
							+4.0*A11*s1*s1*s1+2.0*A12*s1*s2*s2
							+6.0*A111*pow(s1,5.0)+A112*(2.0*s1*pow(s2,4.0)+4.0*s2*s2*pow(s1,3.0))
							+8.0*A1111*pow(s1,7.0)
							+A1112*(6.0*pow(s1,5.0)*s2*s2+2.0*pow(s2,6.0)*s1)
							+4.0*A1122*pow(s1,3.0)*pow(s2,4.0);

			s2k_chem=2.0*A1*(temp-Tc0)*s2
							+4.0*A11*s2*s2*s2+2.0*A12*s2*s1*s1
							+6.0*A111*pow(s2,5.0)+A112*(2.0*s2*pow(s1,4.0)+4.0*s1*s1*pow(s2,3.0))
							+8.0*A1111*pow(s2,7.0)
							+A1112*(6.0*pow(s2,5.0)*s1*s1+2.0*pow(s1,6.0)*s2)
							+4.0*A1122*pow(s2,3.0)*pow(s1,4.0);

			//Calculation of dipole-dipole potential [Formula (4.56) right side first term]
			//s1k_ddi=Add0*s1k_dd[i][j];
			//s2k_ddi=Add0*s2k_dd[i][j];
			s1k_ddi=Add0*s1k_dd[i*ND+j];
			s2k_ddi=Add0*s2k_dd[i*ND+j];

			//Computation of total potential
			s1k=s1k_chem+s1k_surf+s1k_ddi-E1_ex_0;
			s2k=s2k_chem+s2k_surf+s2k_ddi;

			//Evolution equation [equation (4.57)] and time evolution of s-field (explicit method) calculation
			//s1ddtt=-smob1*s1k;  s1h2[i][j]=s1h[i][j]+s1ddtt*delt;
			//s2ddtt=-smob2*s2k;  s2h2[i][j]=s2h[i][j]+s2ddtt*delt;
			s1ddtt=-smob1*s1k;  s1h2[i*ND+j]=s1h[i*ND+j]+s1ddtt*delt;
			s2ddtt=-smob2*s2k;  s2h2[i*ND+j]=s2h[i*ND+j]+s2ddtt*delt;
		}
	}

	//Move the s-field to the main array and correct the numerical errors in the domain
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s1h[i][j]=s1h2[i][j]; s2h[i][j]=s2h2[i][j];
			s1h[i*ND+j]=s1h2[i*ND+j]; s2h[i*ND+j]=s2h2[i*ND+j];
			//if(s1h[i][j]>=Peq){s1h[i][j]=Peq;}  if(s1h[i][j]<=-Peq){s1h[i][j]=-Peq;}
			if(s1h[i*ND+j]>=Peq){s1h[i*ND+j]=Peq;}  if(s1h[i*ND+j]<=-Peq){s1h[i*ND+j]=-Peq;}
			//if(s2h[i][j]>=Peq){s2h[i][j]=Peq;}  if(s2h[i][j]<=-Peq){s2h[i][j]=-Peq;}
			if(s2h[i*ND+j]>=Peq){s2h[i*ND+j]=Peq;}  if(s2h[i*ND+j]<=-Peq){s2h[i*ND+j]=-Peq;}
			//ss1ss2=sqrt(s1h[i][j]*s1h[i][j]+s2h[i][j]*s2h[i][j]);
			ss1ss2=sqrt(s1h[i*ND+j]*s1h[i*ND+j]+s2h[i*ND+j]*s2h[i*ND+j]);
			//if(ss1ss2>=Peq){s1h[i][j]=s1h[i][j]/ss1ss2*Peq; s2h[i][j]=s2h[i][j]/ss1ss2*Peq;}
			if(ss1ss2>=Peq){s1h[i*ND+j]=s1h[i*ND+j]/ss1ss2*Peq; s2h[i*ND+j]=s2h[i*ND+j]/ss1ss2*Peq;}
		}
	}

	//Advance the time, and if it is not the end time, go to the start
	//if(keypress()){return 0;}//Waiting for key
	time1=time1+1.0;								//Add calculation count
	if(time1<time1max){goto start;}	//Determining if the maximum count has been reached
	printf("Finished \n");

end:;
	std::exit(0);
	//return 0;
}

//************ Initial field setup subroutine *************
void ini000_s12(double *s1h, double *s2h, int ND)
{
	int i, j;	//integer
	//srand(time(NULL)); // random number initialization
	int nd=ND, ndm=ND-1, nd2=ND/2;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s1h[i][j]=0.01*(2.0*DRND(1)-1.0);//Set the field with random numbers up to Å}1%
			//s2h[i][j]=0.01*(2.0*DRND(1)-1.0);
			s1h[i*ND+j]=0.01*(2.0*DRND(1)-1.0);//Set the field with random numbers up to Å}1%
			s2h[i*ND+j]=0.01*(2.0*DRND(1)-1.0);
		}
	}
}

//******* Sin, Cos table and bit reversal table settings ***************
void table(double *s, double *c, int *ik, int ND, int ig)
{
	int it, it1, it2, mc, mn;
	double q;
	int nd=ND, ndm=ND-1, nd2=ND/2;

	q=2.*PI/nd;
	for(it=0;it<=nd2-1;it++){
		c[it]=cos(q*it);  s[it]=sin(q*it);//Sin, Cos table
	}
	ik[0]=0;
	mn=nd2;
	mc=1;
	for(it1=1;it1<=ig;it1++){
		for(it2=0;it2<=mc-1;it2++){
			ik[it2+mc]=ik[it2]+mn;//bit reversal table
		}
		mn=mn/2;
		mc=2*mc;
	}
}

//********** One-dimensional fast Fourier transform **************************************
void fft(double *xrf, double *xif, double *s, double *c, int ND, int ig)
{
	int ix, ka, kb, l2, lf, mf, n2, nf;
	double tj, tr;
	int nd=ND, ndm=ND-1, nd2=ND/2;

	l2=1;
	for(lf=1;lf<=ig;lf++){
		n2=nd2/l2;
		for(mf=1;mf<=l2;mf++){
			for(nf=0;nf<=n2-1;nf++){
				ix=nf*l2;
				ka=nf+2*n2*(mf-1);
				kb=ka+n2;
				tr=xrf[ka]-xrf[kb];
				tj=xif[ka]-xif[kb];
				xrf[ka]=xrf[ka]+xrf[kb];
				xif[ka]=xif[ka]+xif[kb];
				xrf[kb]=tr*c[ix]-tj*qs*s[ix];
				xif[kb]=tj*c[ix]+tr*qs*s[ix];
			}
		}
		l2=l2*2;
	}

}

//************ 2D Fast Fourier Transform ***********************************
void rcfft(double *xrf, double *xif, double *xr, double *xi, int *ik, double *s, double *c, int ND, int ig)
{
	int i, ic, ir, j;
	int nd=ND, ndm=ND-1, nd2=ND/2;

	for(ir=0;ir<=ndm;ir++){
		for(ic=0;ic<=ndm;ic++){
			//xrf[ic]=xr[ir][ic];  xif[ic]=xi[ir][ic];
			xrf[ic]=xr[ir*ND+ic];
			xif[ic]=xi[ir*ND+ic];
		}
		fft(xrf, xif, s, c, ND, ig);
		for(ic=0;ic<=ndm;ic++){
			//xr[ir][ic]=xrf[ik[ic]];  xi[ir][ic]=xif[ik[ic]];
			xr[ir*ND+ic]=xrf[ik[ic]];
			xi[ir*ND+ic]=xif[ik[ic]];
		}
	}
	for(ic=0;ic<=ndm;ic++){
		for(ir=0;ir<=ndm;ir++){
			//xrf[ir]=xr[ir][ic];  xif[ir]=xi[ir][ic];
			xrf[ir]=xr[ir*ND+ic];
			xif[ir]=xi[ir*ND+ic];
		}
		fft(xrf, xif, s, c, ND, ig);
		for(ir=0;ir<=ndm;ir++){
			//xr[ir][ic]=xrf[ik[ir]];  xi[ir][ic]=xif[ik[ir]];
			xr[ir*ND+ic]=xrf[ik[ir]];
			xi[ir*ND+ic]=xif[ik[ir]];
		}
	}
	if(qs>0.){return;}
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=xr[i][j]/nd/nd;  xi[i][j]=xi[i][j]/nd/nd;
			xr[i*ND+j]=xr[i*ND+j]/nd/nd;
			xi[i*ND+j]=xi[i*ND+j]/nd/nd;
		}
	}

}

//************ data save subroutine *******************************
void datsave(double *s1h, double *s2h, int ND)
{
	FILE		*stream;	//Stream pointer setting
	int 		i, j;			//integer
	int nd=ND, ndm=ND-1, nd2=ND/2;

	stream = fopen("test_ext.dat", "a");	//Open the write destination file for appending
	fprintf(stream, "%e \n", time1);	//Saving calculation counts
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//fprintf(stream, "%e  %e  ", s1h[i][j], s2h[i][j]);	//Field data storage
			fprintf(stream, "%e  %e  ", s1h[i*ND+j], s2h[i*ND+j]);	//Field data storage
		}
	}
	fprintf(stream, "\n");	//writing a newline
	fclose(stream);					//close file
}


void datsave_paraview(double *s1h, double *s2h, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j;
	int nd=ND, ndm=ND-1, nd2=ND/2;
	
	iout = iout + 1;
	printf("dip_ext_result%06d.vtk \n",iout);
	sprintf(fName,"dip_ext_result%06d.vtk",iout);
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
			//fprintf(stream, "%e  %e  ", s1h[i][j], s2h[i][j]);	//Field data storage
			//fprintf(fp,"%10.6f\n", s1h[i][j]);
			fprintf(fp,"%10.6f\n", s1h[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Phase_field_2 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e  ", s1h[i][j], s2h[i][j]);	//Field data storage
			//fprintf(fp,"%10.6f\n", s2h[i][j]);
			fprintf(fp,"%10.6f\n", s2h[i*ND+j]);
		}
	}
	fclose(fp);
}

//************ data loading subroutine *******************************
void datin(double *s1h, double *s2h, int ND)
{
	FILE		*datin0;	//Stream pointer setting
	int 		i, j;			//integer
	double time_ext;	//calculation count number
	int nd=ND, ndm=ND-1, nd2=ND/2;

	datin0 = fopen("test.dat", "r");	//open source file
	fscanf(datin0, "%lf", &time_ext);	//Read calculation count
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//fscanf(datin0, "%lf  %lf  ", &s1h[i][j], &s2h[i][j]);	//Field data read
			fscanf(datin0, "%lf  %lf  ", &s1h[i*ND+j], &s2h[i*ND+j]);	//Field data read
		}
	}
	fclose(datin0);//close file
}

