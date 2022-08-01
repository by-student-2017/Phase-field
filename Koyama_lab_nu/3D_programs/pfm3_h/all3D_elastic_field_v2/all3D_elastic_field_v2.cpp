#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//#include "wingxa.h"

#define DRND(x) ((double)(x)/RAND_MAX*rand())
//#define ND 64
//#define IG 6	//2^IG=ND, IG=int(log2(ND))
//#define INX 400					//Pixel size (x direction)
//#define INY 400					//Pixel size (y direction)

	//int nd=ND; 					
	//int ndm=ND-1;				//define ND-1
	//int nd2=ND/2;				//define ND/2 for FFT
	//int ig=IG;					//2^ig=ND
	double PI=3.14159;			//pi
	double rr=8.3145;			//Gas constant
	double time1;				//count
	int iout=0;

	double c0;					//average of phase field
	//double ch[ND][ND][ND];		//phase field
	double eta_c[4][4];			//Lattice mismatch
	double cec[4][4][4][4];		//Elastic constant
	double sigma[4][4];			//Eigen stress (stress when elastically deformed by the amount of Eigen strain)
	//double ep11c[ND][ND][ND], ep22c[ND][ND][ND], ep33c[ND][ND][ND];	//Strain
	//double ep12c[ND][ND][ND], ep13c[ND][ND][ND], ep23c[ND][ND][ND];
	//double sig11[ND][ND][ND], sig22[ND][ND][ND], sig33[ND][ND][ND];	//Elastic stress
	//double sig12[ND][ND][ND], sig13[ND][ND][ND], sig23[ND][ND][ND];
	//double Estr[ND][ND][ND];	//Elastic strain energy density
	//double u1[ND][ND][ND], u2[ND][ND][ND], u3[ND][ND][ND];	//Displacement
	//double fld1[ND][ND][ND];

	//bag fix: move them from main program to here
	//double qrh1[ND][ND][ND], qih1[ND][ND][ND];	//Fourier transform of the field
	//double ec11[ND][ND][ND], ec22[ND][ND][ND], ec33[ND][ND][ND];	//Constrained strain array
	//double ec12[ND][ND][ND], ec13[ND][ND][ND], ec23[ND][ND][ND];

	int qs;						//FFT(qs:-1) or IFFT(qs:1)
	//double xi[ND][ND][ND], xr[ND][ND][ND], xif[ND], xrf[ND];//Real or imaginary array of Fourier transform
	//double s[ND],c[ND];			//sin and cos table
	//int ik[ND];				//Bit inversion table for FFT

	void table(double *s, double *c, int *ik, int ND, int ig);	//Creating sin and cos tables and bit inversion tables Subrutin
	void fft(double *xrf, double *xif, double *s, double *c, int ND, int ig);	//FFT1D
	void xyzfft(double *xrf, double *xif, double *xr, double *xi, int *ik, double *s, double *c, int ND, int ig);	//FFT3D
	double zcij(int i0, int ij0, int k0, int iii, int jjj, int ND);//Coefficient calculation of elastic function (Fourier space)
	double zuij(int i0, int ij0, int k0, int iii, int ND);//Displacement field coefficient calculation (Fourier space)

	void datin(double *ch, int ND);	//Subroutine for initial field reading
	void datsave(double *ch, double *Estr, 
			 double *ep11c, double *ep22c, double *ep33c, double *ep12c, double *ep13c, double *ep23c, 
			 double *sig11, double *sig22, double *sig33, double *sig12, double *sig13, double *sig23, 
			 double *u1, double *u2, double *u3, int ND);	//Data storage subroutine
	void datsave_paraview(double *ch, double *Estr, 
			 double *ep11c, double *ep22c, double *ep33c, double *ep12c, double *ep13c, double *ep23c, 
			 double *sig11, double *sig22, double *sig33, double *sig12, double *sig13, double *sig23, 
			 double *u1, double *u2, double *u3, int ND);	//Data storage subroutine

//******* main program ****************************************************
int main(void)
{
    int ND, IG;
	int nd, ndm, nd2, ig;
	
	//int i, j, k, l, ii, jj, kk, iii, jjj, ief;	//integer
	int i, j, k, l, iii, jjj;	//integer
	//int ip, im, jp, jm;			//integer

	//double c;					//phase field
	double al;		//Length of one side of the calculation area
	//double temp;				//Temperature [K]
	//double delt;				//Time increments
	double time1max;			//Maximum value of calculation count (used to stop calculation)
	double b1;					//The length of one side of the difference block
	double vm0;					//Molar volume
	double atom_n;				//Number of atoms in a unit cell
	double a1;					//Lattice constant of Fe
	
	double el_mag;				//Working variables related to elastic modulus
	double c11, c44, c12; 		//Modulus
	double ep000;				//Lattice mismatch
	//double epc11, epc22, epc33, epc12, epc13, epc23;	//Restraint strain
	//double ep011, ep022, ep033, ep012, ep013, ep023;	//Eigen strain
	//double dep011, dep022, dep033, dep012, dep023, dep013;
	//double epc011, epc022, epc033, epc012, epc023, epc013;
	//double ef1, ef2, ef3, ef11, ef22, ef33, ef12, ef13, ef23;	//Elastic function
	
	double Estr1, Estr2, Estr3, Estr4, Estr5;//Working variables when calculating elastic strain energy

//****** Setting of calculation conditions and material constants ****************************************
	printf("---------------------------------\n");
	printf("read parameters form parameters.txt\n");
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
	// pfm1 series
	ND      = int(data[0]);
	//delt    = data[1];
	//c0      = data[2];
	//Mx      = data[3];
	//My      = data[4];
	//Mz      = data[5];
	//al      = data[6];	// [nm]
	//amob_c  = data[7];
	//time1max= data[8];
	//temp    = data[9];	// [K]
	//L00     = data[10];	//L0=L00/RR/temp;
	//kapa_c0 = data[11];	//kapa_c=kapa_c0/b1/b1/RR/temp;
	//c_flu   = data[12];
	//flg     = int(data[13]);
	//Nstep   = int(data[14]);
	ep000   = data[15];
	el_mag  = data[16];
	c11     = data[17];
	c12     = data[18];
	c44     = data[19];
	printf("---------------------------------\n");
	//
	IG=int(log2(ND));
	nd=ND; 					
	ndm=ND-1;				//define ND-1
	nd2=ND/2;				//define ND/2 for FFT
	ig=IG;					//2^ig=ND
	
	//double ch[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	//a[z * ySize * xSize + y * xSize + x], a[i][j][k]=a[i*n*m + j*n + k]
	double *ch    = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//phase field
	//
	double *ep11c = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Strain
	double *ep22c = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep33c = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep12c = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep13c = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep23c = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	//
	double *sig11 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Elastic stress
	double *sig22 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sig33 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sig12 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sig13 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sig23 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	//
	double *Estr  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Elastic strain energy density
	//
	double *u1    = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Displacement
	double *u2    = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *u3    = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	//
	double *qrh1  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of the field
	double *qih1  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	//
	double *ec11  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Constrained strain array
	double *ec22  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ec33  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ec12  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ec13  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ec23  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	//
	double *xr    = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Real or imaginary array of Fourier transform
	double *xi    = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *xrf   = (double *)malloc(sizeof(double)*( ND ));
	double *xif   = (double *)malloc(sizeof(double)*( ND ));
	//
	double *s     = (double *)malloc(sizeof(double)*( ND ));	//sin and cos table
	double *c     = (double *)malloc(sizeof(double)*( ND ));
	//
	int *ik       = (int *)malloc(sizeof(int)*( ND ));	//Bit inversion table for FFT
	
	//printf("delt(0.005)=  "); scanf(" %lf",&delt);	//Time increments
	//printf("c0(0.4)=  "); scanf(" %lf",&c0);	//Average of phase field
	//printf("temp(900.0)=  "); scanf(" %lf",&temp);	//Temperature [K]
	//delt=0.002;					//time increments
	//c0 = 0.4;					//Average of phase field
	//temp=900.0;					//Temperature [K]
	//al=100.0*1.0E-09;			//calculation area [m]
	//al=al*1.0e-9;	// [m]
	//b1=al/nd;					//Difference block length

	//time1=0.0;					//Setting the initial calculation count number
	//time1max=1.0+1.0e+07;		//Setting the maximum calculation count

	//a1=2.8664E-10; 				//Lattice constant (experimental data) of Fe(bcc)
	//atom_n=2.0;
	//vm0=6.02E23*a1*a1*a1/atom_n;	//Calculation of molar volume (for bcc)

//*** Lattice mismatch settings ***
	//ep000=0.05;
	eta_c[1][1]=eta_c[2][2]=eta_c[3][3]=ep000; 	//Lattice mismatch
	eta_c[1][2]=eta_c[2][1]=eta_c[1][3]=eta_c[3][1]=eta_c[2][3]=eta_c[3][2]=0.0;

//***** Elastic constant of Fe (bcc) ****************************
	//el_mag=1.0E+11/1.0E+09;
	//c11=2.33*el_mag;
	//c12=1.35*el_mag;
	//c44=1.18*el_mag;
	c11=c11*el_mag;
	c12=c12*el_mag;
	c44=c44*el_mag;

	for(i=0;i<=3;i++){
		for(j=0;j<=3;j++){
			for(k=0;k<=3;k++){
				for(l=0;l<=3;l++){
					cec[i][j][k][l]=0.0;//Initialization of elastic constant array
				}
			}
		}
	}

	cec[1][1][1][1]=c11;
	cec[2][2][2][2]=c11;
	cec[3][3][3][3]=c11;
	cec[1][2][1][2]=cec[1][2][2][1]=cec[2][1][1][2]=cec[2][1][2][1]=c44;
	cec[2][3][2][3]=cec[2][3][3][2]=cec[3][2][2][3]=cec[3][2][3][2]=c44;
	cec[1][3][1][3]=cec[1][3][3][1]=cec[3][1][1][3]=cec[3][1][3][1]=c44;
	cec[1][1][2][2]=cec[2][2][1][1]=c12;
	cec[1][1][3][3]=cec[3][3][1][1]=c12;
	cec[2][2][3][3]=cec[3][3][2][2]=c12;

//--- Eigen stress (stress when elastically deformed by the amount of Eigen strain) --------------
	sigma[1][1]=cec[1][1][1][1]*eta_c[1][1]
			   +cec[1][1][2][2]*eta_c[2][2]
			   +cec[1][1][3][3]*eta_c[3][3]; //sigma1
	sigma[2][2]=cec[2][2][1][1]*eta_c[1][1]
			   +cec[2][2][2][2]*eta_c[2][2]
			   +cec[2][2][3][3]*eta_c[3][3]; //sigma2
	sigma[3][3]=cec[3][3][1][1]*eta_c[1][1]
			   +cec[3][3][2][2]*eta_c[2][2]
			   +cec[3][3][3][3]*eta_c[3][3]; //sigma3
	sigma[2][3]=cec[2][3][2][3]*eta_c[2][3]+cec[2][3][3][2]*eta_c[3][2]; //sigma4
	sigma[3][2]=cec[3][2][2][3]*eta_c[2][3]+cec[3][2][3][2]*eta_c[3][2]; //sigma4
	sigma[1][3]=cec[1][3][1][3]*eta_c[1][3]+cec[1][3][3][1]*eta_c[3][1]; //sigma5
	sigma[3][1]=cec[3][1][1][3]*eta_c[1][3]+cec[3][1][3][1]*eta_c[3][1]; //sigma5
	sigma[1][2]=cec[1][2][1][2]*eta_c[1][2]+cec[1][2][2][1]*eta_c[2][1]; //sigma6
	sigma[2][1]=cec[2][1][1][2]*eta_c[1][2]+cec[2][1][2][1]*eta_c[2][1]; //sigma6

//*** sin and cos table, bit inversion table, and initial field settings ***************

	table(s, c, ik, ND, ig);
	//ini000();
	datin(ch, ND); //Read phase field

 	//gwinsize(INX,INY); ginit(1); gsetorg(0,0);// drow figure

//**** Start of calculation of elastic field analysis ******************************
//start:;

//**** Fourier transform of the field (phase field) ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=ch[i][j][k]; xi[i][j][k]=0.;
				xr[i*ND*ND+j*ND+k]=ch[i*ND*ND+j*ND+k];
				xi[i*ND*ND+j*ND+k]=0.0;
			}
		}
	}
	qs=-1; xyzfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//qrh1[i][j][k]=xr[i][j][k]; qih1[i][j][k]=xi[i][j][k];
				qrh1[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
				qih1[i*ND*ND+j*ND+k]=xi[i*ND*ND+j*ND+k];
			}
		}
	}
	//qrh1[0][0][0]=qih1[0][0][0]=0.0;
	qrh1[0]=qih1[0]=0.0;

//***** Calculation of total strain fluctuation *************************************
//--- ec11 ---
	iii=1; jjj=1;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				//xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xr[i*ND*ND+j*ND+k]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)+qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
				xi[i*ND*ND+j*ND+k]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)+qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
			}
		}
	}
	qs=1; xyzfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//ec11[i][j][k]=xr[i][j][k];
				ec11[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- ec22 ---
	iii=2; jjj=2;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				//xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xr[i*ND*ND+j*ND+k]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)+qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
				xi[i*ND*ND+j*ND+k]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)+qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
			}
		}
	}
	qs=1; xyzfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//ec22[i][j][k]=xr[i][j][k];
				ec22[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- ec33 ---
	iii=3; jjj=3;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				//xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xr[i*ND*ND+j*ND+k]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)+qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
				xi[i*ND*ND+j*ND+k]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)+qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
			}
		}
	}
	qs=1; xyzfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//ec33[i][j][k]=xr[i][j][k];
				ec33[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- ec12 ---
	iii=1; jjj=2;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				//xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xr[i*ND*ND+j*ND+k]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)+qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
				xi[i*ND*ND+j*ND+k]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)+qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
			}
		}
	}
	qs=1; xyzfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//ec12[i][j][k]=xr[i][j][k];
				ec12[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- ec13 ---
	iii=1; jjj=3;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				//xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xr[i*ND*ND+j*ND+k]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)+qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
				xi[i*ND*ND+j*ND+k]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)+qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
			}
		}
	}
	qs=1; xyzfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//ec13[i][j][k]=xr[i][j][k];
				ec13[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- ec23 ---
	iii=2; jjj=3;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				//xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xr[i*ND*ND+j*ND+k]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)+qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
				xi[i*ND*ND+j*ND+k]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)+qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
			}
		}
	}
	qs=1; xyzfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//ec23[i][j][k]=xr[i][j][k];
				ec23[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//***** Calculation of total strain field, stress field, and elastic strain energy field *************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//Calculation of total strain field
				//ep11c[i][j][k]=ec11[i][j][k]+ep000*c0;
				//ep22c[i][j][k]=ec22[i][j][k]+ep000*c0;
				//ep33c[i][j][k]=ec33[i][j][k]+ep000*c0;
				//ep12c[i][j][k]=ec12[i][j][k];
				//ep13c[i][j][k]=ec13[i][j][k];
				//ep23c[i][j][k]=ec23[i][j][k];
				ep11c[i*ND*ND+j*ND+k]=ec11[i*ND*ND+j*ND+k]+ep000*c0;
				ep22c[i*ND*ND+j*ND+k]=ec22[i*ND*ND+j*ND+k]+ep000*c0;
				ep33c[i*ND*ND+j*ND+k]=ec33[i*ND*ND+j*ND+k]+ep000*c0;
				ep12c[i*ND*ND+j*ND+k]=ec12[i*ND*ND+j*ND+k];
				ep13c[i*ND*ND+j*ND+k]=ec13[i*ND*ND+j*ND+k];
				ep23c[i*ND*ND+j*ND+k]=ec23[i*ND*ND+j*ND+k];

				//Calculation of elastic stress field
				//sig11[i][j][k]=c11*ec11[i][j][k]+c12*ec22[i][j][k]+c12*ec33[i][j][k]
				//							   -(c11+2.*c12)*ep000*(ch[i][j][k]-c0);
				//sig22[i][j][k]=c12*ec11[i][j][k]+c11*ec22[i][j][k]+c12*ec33[i][j][k]
				//							   -(c11+2.*c12)*ep000*(ch[i][j][k]-c0);
				//sig33[i][j][k]=c12*ec11[i][j][k]+c12*ec22[i][j][k]+c11*ec33[i][j][k]
				//							   -(c11+2.*c12)*ep000*(ch[i][j][k]-c0);
				//sig12[i][j][k]=2.*c44*ec12[i][j][k];
				//sig13[i][j][k]=2.*c44*ec13[i][j][k];
				//sig23[i][j][k]=2.*c44*ec23[i][j][k];
				sig11[i*ND*ND+j*ND+k]=c11*ec11[i*ND*ND+j*ND+k]+c12*ec22[i*ND*ND+j*ND+k]+c12*ec33[i*ND*ND+j*ND+k]
											   -(c11+2.*c12)*ep000*(ch[i*ND*ND+j*ND+k]-c0);
				sig22[i*ND*ND+j*ND+k]=c12*ec11[i*ND*ND+j*ND+k]+c11*ec22[i*ND*ND+j*ND+k]+c12*ec33[i*ND*ND+j*ND+k]
											   -(c11+2.*c12)*ep000*(ch[i*ND*ND+j*ND+k]-c0);
				sig33[i*ND*ND+j*ND+k]=c12*ec11[i*ND*ND+j*ND+k]+c12*ec22[i*ND*ND+j*ND+k]+c11*ec33[i*ND*ND+j*ND+k]
											   -(c11+2.*c12)*ep000*(ch[i*ND*ND+j*ND+k]-c0);
				sig12[i*ND*ND+j*ND+k]=2.*c44*ec12[i*ND*ND+j*ND+k];
				sig13[i*ND*ND+j*ND+k]=2.*c44*ec13[i*ND*ND+j*ND+k];
				sig23[i*ND*ND+j*ND+k]=2.*c44*ec23[i*ND*ND+j*ND+k];

				//Calculation of elastic strain energy field
				//Estr1=1.5*(c11+2.0*c12)*ep000*ep000*ch[i][j][k]*ch[i][j][k];
				//Estr2=-(c11+2.0*c12)*(ep11c[i][j][k]
				//	                 +ep22c[i][j][k]
				//	                 +ep33c[i][j][k])*ep000*ch[i][j][k];
				//Estr3=0.5*c11*(ep11c[i][j][k]*ep11c[i][j][k]
				//			  +ep22c[i][j][k]*ep22c[i][j][k]
				//			  +ep33c[i][j][k]*ep33c[i][j][k]);
				//Estr4=c12*(ep11c[i][j][k]*ep22c[i][j][k]
				//		  +ep11c[i][j][k]*ep33c[i][j][k]
				//		  +ep22c[i][j][k]*ep33c[i][j][k]);
				//Estr5=2.*c44*(ep12c[i][j][k]*ep12c[i][j][k]
				//			 +ep13c[i][j][k]*ep13c[i][j][k]
				//			 +ep23c[i][j][k]*ep23c[i][j][k]);
				//Estr[i][j][k]=Estr1+Estr2+Estr3+Estr4+Estr5;
				Estr1=1.5*(c11+2.0*c12)*ep000*ep000*ch[i*ND*ND+j*ND+k]*ch[i*ND*ND+j*ND+k];
				Estr2=-(c11+2.0*c12)*(ep11c[i*ND*ND+j*ND+k]
					                 +ep22c[i*ND*ND+j*ND+k]
					                 +ep33c[i*ND*ND+j*ND+k])*ep000*ch[i*ND*ND+j*ND+k];
				Estr3=0.5*c11*(ep11c[i*ND*ND+j*ND+k]*ep11c[i*ND*ND+j*ND+k]
							  +ep22c[i*ND*ND+j*ND+k]*ep22c[i*ND*ND+j*ND+k]
							  +ep33c[i*ND*ND+j*ND+k]*ep33c[i*ND*ND+j*ND+k]);
				Estr4=c12*(ep11c[i*ND*ND+j*ND+k]*ep22c[i*ND*ND+j*ND+k]
						  +ep11c[i*ND*ND+j*ND+k]*ep33c[i*ND*ND+j*ND+k]
						  +ep22c[i*ND*ND+j*ND+k]*ep33c[i*ND*ND+j*ND+k]);
				Estr5=2.*c44*(ep12c[i*ND*ND+j*ND+k]*ep12c[i*ND*ND+j*ND+k]
							 +ep13c[i*ND*ND+j*ND+k]*ep13c[i*ND*ND+j*ND+k]
							 +ep23c[i*ND*ND+j*ND+k]*ep23c[i*ND*ND+j*ND+k]);
				Estr[i*ND*ND+j*ND+k]=Estr1+Estr2+Estr3+Estr4+Estr5;
			}
		}
	}


//***** Displacement field calculation *************************************
//--- u1 ---
	iii=1;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=-qrh1[i][j][k]*zuij(i, j, k, iii)-qih1[i][j][k]*zuij(i, j, k, iii);
				//xi[i][j][k]= qrh1[i][j][k]*zuij(i, j, k, iii)+qih1[i][j][k]*zuij(i, j, k, iii);
				xr[i*ND*ND+j*ND+k]=-qrh1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND)-qih1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND);
				xi[i*ND*ND+j*ND+k]= qrh1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND)+qih1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND);
			}
	 	}
	}
	qs=1; xyzfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//u1[i][j][k]=xr[i][j][k];
				u1[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- u2 ---
	iii=2;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=-qrh1[i][j][k]*zuij(i, j, k, iii)-qih1[i][j][k]*zuij(i, j, k, iii);
				//xi[i][j][k]= qrh1[i][j][k]*zuij(i, j, k, iii)+qih1[i][j][k]*zuij(i, j, k, iii);
				xr[i*ND*ND+j*ND+k]=-qrh1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND)-qih1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND);
				xi[i*ND*ND+j*ND+k]= qrh1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND)+qih1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND);
			}
	 	}
	}
	qs=1; xyzfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//u2[i][j][k]=xr[i][j][k];
				u2[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- u3 ---
	iii=3;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=-qrh1[i][j][k]*zuij(i, j, k, iii)-qih1[i][j][k]*zuij(i, j, k, iii);
				//xi[i][j][k]= qrh1[i][j][k]*zuij(i, j, k, iii)+qih1[i][j][k]*zuij(i, j, k, iii);
				xr[i*ND*ND+j*ND+k]=-qrh1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND)-qih1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND);
				xi[i*ND*ND+j*ND+k]= qrh1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND)+qih1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND);
			}
	 	}
	}
	qs=1; xyzfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//u3[i][j][k]=xr[i][j][k];
				u3[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//***** 全数値データ保存 ********************************************
	//datsave(ch, Estr, ep11c, ep22c, ep33c, ep12c, ep13c, ep23c, 
	//		 sig11, sig22, sig33, sig12, sig13, sig23, u1, u2, u3, ND);			//save data
	datsave_paraview(ch, Estr, ep11c, ep22c, ep33c, ep12c, ep13c, ep23c, 
			 sig11, sig22, sig33, sig12, sig13, sig23, u1, u2, u3, ND);	//save

	//if(keypress()){return 0;}//wait key input
	return 0;

}


//******* Sin, Cos table and bit inversion table settings ***************
void table(double *s, double *c, int *ik, int ND, int ig)
{
	int i, j, mc, mn;
	double q;
	int nd=ND, ndm=ND-1, nd2=ND/2;

	q=2.0*PI/nd;
	for(i=0;i<=nd2-1;i++){ c[i]=cos(q*i); s[i]=sin(q*i); }//Sin, Cos

	ik[0]=0; mn=nd2; mc=1;
	for(i=1;i<=ig;i++){
		for(j=0;j<=mc-1;j++){ ik[j+mc]=ik[j]+mn; }	//Bit inversion table
		mc*=2; mn/=2; 
	}
}

//********** FFT 1D **************************************
void fft(double *xrf, double *xif, double *s, double *c, int ND, int ig)
{
	int lg, lf, mf, nf, n2, ka, kb, ix;
	double tj, tr;
	int nd=ND, ndm=ND-1, nd2=ND/2;

	lg=1;
	for(lf=1;lf<=ig;lf++){
		n2=nd2/lg;
		for(mf=1;mf<=lg;mf++){
			for(nf=0;nf<=n2-1;nf++){
				ix=nf*lg;
				ka=nf+2*n2*(mf-1);
				kb=ka+n2;
				tr=xrf[ka]-xrf[kb];
				tj=xif[ka]-xif[kb];
				xrf[ka]+=xrf[kb];
				xif[ka]+=xif[kb];
				xrf[kb]=tr*c[ix]-tj*float(qs)*s[ix];
				xif[kb]=tj*c[ix]+tr*float(qs)*s[ix];
			}
		}
		lg*=2;
	}//lf
}


//************ FFT 3D ***********************************
void xyzfft(double *xrf, double *xif, double *xr, double *xi, int *ik, double *s, double *c, int ND, int ig)
{
	int i, j, k, iz, iy, ix;
	int nd=ND, ndm=ND-1, nd2=ND/2;

//*** xy-z ***
	for(ix=0;ix<=ndm;ix++){
		for(iy=0;iy<=ndm;iy++){
			for(iz=0;iz<=ndm;iz++){
				//xrf[iz]=xr[ix][iy][iz];
				//xif[iz]=xi[ix][iy][iz];
				xrf[iz]=xr[ix*ND*ND+iy*ND+iz];
				xif[iz]=xi[ix*ND*ND+iy*ND+iz];
			}
			fft(xrf, xif, s, c, ND, ig);
			for(iz=0;iz<=ndm;iz++){
				//xr[ix][iy][iz]=xrf[ik[iz]];
				//xi[ix][iy][iz]=xif[ik[iz]];
				xr[ix*ND*ND+iy*ND+iz]=xrf[ik[iz]];
				xi[ix*ND*ND+iy*ND+iz]=xif[ik[iz]];
			}
		}
	}

//*** xz-y ***
	for(ix=0;ix<=ndm;ix++){
		for(iz=0;iz<=ndm;iz++){
			for(iy=0;iy<=ndm;iy++){
				//xrf[iy]=xr[ix][iy][iz];
				//xif[iy]=xi[ix][iy][iz];
				xrf[iy]=xr[ix*ND*ND+iy*ND+iz];
				xif[iy]=xi[ix*ND*ND+iy*ND+iz];
			}
			fft(xrf, xif, s, c, ND, ig);
			for(iy=0;iy<=ndm;iy++){
				//xr[ix][iy][iz]=xrf[ik[iy]];
				//xi[ix][iy][iz]=xif[ik[iy]];
				xr[ix*ND*ND+iy*ND+iz]=xrf[ik[iy]];
				xi[ix*ND*ND+iy*ND+iz]=xif[ik[iy]];
			}
		}
	}

//*** yz-x ***
	for(iz=0;iz<=ndm;iz++){
		for(iy=0;iy<=ndm;iy++){
			for(ix=0;ix<=ndm;ix++){
				//xrf[ix]=xr[ix][iy][iz];
				//xif[ix]=xi[ix][iy][iz];
				xrf[ix]=xr[ix*ND*ND+iy*ND+iz];
				xif[ix]=xi[ix*ND*ND+iy*ND+iz];
			}
			fft(xrf, xif, s, c, ND, ig);
			for(ix=0;ix<=ndm;ix++){
				//xr[ix][iy][iz]=xrf[ik[ix]];
				//xi[ix][iy][iz]=xif[ik[ix]];
				xr[ix*ND*ND+iy*ND+iz]=xrf[ik[ix]];
				xi[ix*ND*ND+iy*ND+iz]=xif[ik[ix]];
			}
		}
	}

	if(qs > 0){return;}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=xr[i][j][k]/nd/nd/nd;
				//xi[i][j][k]=xi[i][j][k]/nd/nd/nd;
				xr[ix*ND*ND+iy*ND+iz]=xr[ix*ND*ND+iy*ND+iz]/nd/nd/nd;
				xi[ix*ND*ND+iy*ND+iz]=xi[ix*ND*ND+iy*ND+iz]/nd/nd/nd;
			}
		}
	}

}

//*** Zcij [eq.(5.26)] ****************************************
double zcij(int i0, int ij0, int k0, int iii, int jjj, int ND)
{
	//int i, j, k, m, n, p, q;
	int m, n;
	int ii, jj, kk;
	double nx, ny, nz, alnn;
	double nec[4];
	double a11, a22, a33, a12, a13, a23, a21, a31, a32;
	double b11, b22, b33, b12, b13, b23, b21, b31, b32;
	double det1, zij;
	double om[4][4];
	int nd=ND, ndm=ND-1, nd2=ND/2;

	if(i0<=nd2-1){ii=i0;}  if(i0>=nd2){ii=i0-nd;}
	if(ij0<=nd2-1){jj=ij0;}  if(ij0>=nd2){jj=ij0-nd;}
	if(k0<=nd2-1){kk=k0;}  if(k0>=nd2){kk=k0-nd;}
	alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj+(double)kk*(double)kk);
	if(alnn==0.){alnn=1.;}
	nec[1]=nx=(double)ii/alnn;
	nec[2]=ny=(double)jj/alnn;
	nec[3]=nz=(double)kk/alnn;

	//for(i=1;i<=3;i++){ for(j=1;j<=3;j++){ om[i][j]=0.; } }

	a11=cec[1][1][1][1]*nx*nx+cec[1][2][1][2]*ny*ny+cec[1][3][1][3]*nz*nz;
	a22=cec[1][2][1][2]*nx*nx+cec[2][2][2][2]*ny*ny+cec[2][3][2][3]*nz*nz;
	a33=cec[3][1][3][1]*nx*nx+cec[2][3][2][3]*ny*ny+cec[3][3][3][3]*nz*nz;
	a12=(cec[1][1][2][2]+cec[1][2][1][2])*nx*ny;
	a23=(cec[2][2][3][3]+cec[2][3][2][3])*ny*nz;
	a31=(cec[3][3][1][1]+cec[3][1][3][1])*nx*nz;
	a21=a12;
	a32=a23;
	a13=a31;

	b11=a22*a33-a23*a32;
	b22=a11*a33-a13*a31;
	b33=a11*a22-a12*a21;
	b12=-(a21*a33-a23*a31);
	b23=-(a11*a32-a12*a31);
	b31=-(a22*a31-a21*a32);
	b21=b12;
	b32=b23;
	b13=b31;

	det1=a11*a22*a33+a12*a23*a31+a13*a32*a21-a13*a31*a22-a11*a23*a32-a33*a12*a21;
	if(det1==0.){det1=1.;}

	om[1][1]=b11/det1;
	om[2][2]=b22/det1;
	om[3][3]=b33/det1;
	om[1][2]=om[2][1]=b12/det1;
	om[2][3]=om[3][2]=b23/det1;
	om[3][1]=om[1][3]=b31/det1;

	zij=0.0;
	for(m=1;m<=3;m++) {
		for(n=1;n<=3;n++) {
    		zij=zij+0.5*(sigma[m][n]*nec[jjj]*nec[n]*om[m][iii]
				  		+sigma[m][n]*nec[iii]*nec[n]*om[m][jjj]);
		}
	}
	return(zij);
}

//*** Zuij  [eq.(5.30)] ****************************************
double zuij(int i0, int ij0, int k0, int iii, int ND)
{
	//int i, j, k, m, n, p, q;
	int m, n;
	int ii, jj, kk;
	double nx, ny, nz, alnn;
	double nec[4];
	double a11, a22, a33, a12, a13, a23, a21, a31, a32;
	double b11, b22, b33, b12, b13, b23, b21, b31, b32;
	double det1, zij;
	double om[4][4];
	int nd=ND, ndm=ND-1, nd2=ND/2;

	if(i0<=nd2-1){ii=i0;}  if(i0>=nd2){ii=i0-nd;}
	if(ij0<=nd2-1){jj=ij0;}  if(ij0>=nd2){jj=ij0-nd;}
	if(k0<=nd2-1){kk=k0;}  if(k0>=nd2){kk=k0-nd;}
	alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj);
	if(alnn==0.){alnn=1.;}
	nec[1]=nx=(double)ii/alnn;
	nec[2]=ny=(double)jj/alnn;
	nec[3]=nz=(double)kk/alnn;

	//for(i=1;i<=3;i++){ for(j=1;j<=3;j++){ om[i][j]=0.; } }

	a11=cec[1][1][1][1]*nx*nx+cec[1][2][1][2]*ny*ny+cec[1][3][1][3]*nz*nz;
	a22=cec[1][2][1][2]*nx*nx+cec[2][2][2][2]*ny*ny+cec[2][3][2][3]*nz*nz;
	a33=cec[3][1][3][1]*nx*nx+cec[2][3][2][3]*ny*ny+cec[3][3][3][3]*nz*nz;
	a12=(cec[1][1][2][2]+cec[1][2][1][2])*nx*ny;
	a23=(cec[2][2][3][3]+cec[2][3][2][3])*ny*nz;
	a31=(cec[3][3][1][1]+cec[3][1][3][1])*nx*nz;
	a21=a12;
	a32=a23;
	a13=a31;

	b11=a22*a33-a23*a32;
	b22=a11*a33-a13*a31;
	b33=a11*a22-a12*a21;
	b12=-(a21*a33-a23*a31);
	b23=-(a11*a32-a12*a31);
	b31=-(a22*a31-a21*a32);
	b21=b12;
	b32=b23;
	b13=b31;

	det1=a11*a22*a33+a12*a23*a31+a13*a32*a21-a13*a31*a22-a11*a23*a32-a33*a12*a21;
	if(det1==0.){det1=1.;}

	om[1][1]=b11/det1;
	om[2][2]=b22/det1;
	om[3][3]=b33/det1;
	om[1][2]=om[2][1]=b12/det1;
	om[2][3]=om[3][2]=b23/det1;
	om[3][1]=om[1][3]=b31/det1;

	zij=0.;
	for(m=1;m<=3;m++) {
		for(n=1;n<=3;n++) {
    		zij=zij+sigma[m][n]*nec[n]*om[m][iii];
		}
	}
	zij=-zij/alnn;
	return(zij);
}

//************ Saving data of various places *******************************
void datsave(double *ch, double *Estr, 
			 double *ep11c, double *ep22c, double *ep33c, double *ep12c, double *ep13c, double *ep23c, 
			 double *sig11, double *sig22, double *sig33, double *sig12, double *sig13, double *sig23, 
			 double *u1, double *u2, double *u3, int ND)
{
	FILE		*stream;
	int 		i, j, k;
	int nd=ND, ndm=ND-1, nd2=ND/2;

	stream = fopen("el_field.dat", "w");

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				fprintf(stream, "%e  ", ch[i*ND*ND+j*ND+k]);
				fprintf(stream, "%e  ", Estr[i*ND*ND+j*ND+k]);
				fprintf(stream, "%e  ", ep11c[i*ND*ND+j*ND+k]);
				fprintf(stream, "%e  ", ep22c[i*ND*ND+j*ND+k]);
				fprintf(stream, "%e  ", ep33c[i*ND*ND+j*ND+k]);
				fprintf(stream, "%e  ", ep12c[i*ND*ND+j*ND+k]);
				fprintf(stream, "%e  ", ep13c[i*ND*ND+j*ND+k]);
				fprintf(stream, "%e  ", ep23c[i*ND*ND+j*ND+k]);
				fprintf(stream, "%e  ", sig11[i*ND*ND+j*ND+k]);
				fprintf(stream, "%e  ", sig22[i*ND*ND+j*ND+k]);
				fprintf(stream, "%e  ", sig33[i*ND*ND+j*ND+k]);
				fprintf(stream, "%e  ", sig12[i*ND*ND+j*ND+k]);
				fprintf(stream, "%e  ", sig13[i*ND*ND+j*ND+k]);
				fprintf(stream, "%e  ", sig23[i*ND*ND+j*ND+k]);
				fprintf(stream, "%e  ", u1[i*ND*ND+j*ND+k]);
				fprintf(stream, "%e  ", u2[i*ND*ND+j*ND+k]);
				fprintf(stream, "%e  ", u3[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(stream, "\n");
	fclose(stream);
}

void datsave_paraview(double *ch, double *Estr, 
			 double *ep11c, double *ep22c, double *ep33c, double *ep12c, double *ep13c, double *ep23c, 
			 double *sig11, double *sig22, double *sig33, double *sig12, double *sig13, double *sig23, 
			 double *u1, double *u2, double *u3, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j, k;
	int nd=ND, ndm=ND-1, nd2=ND/2;
	
	//iout = iout + 1;
	sprintf(fName,"el_result%06d.vtk",iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),(ndm+1));
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*(ndm+1)));
	fprintf(fp,"SCALARS concentration float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", ch[i][j][k]);
				fprintf(fp,"%10.6f\n", ch[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_strain_energy_density float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", Estr[i][j][k]);
				fprintf(fp,"%10.6f\n", Estr[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Strain_c11 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", ep11c[i][j][k]);
				fprintf(fp,"%10.6f\n", ep11c[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Strain_c22 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", ep22c[i][j][k]);
				fprintf(fp,"%10.6f\n", ep22c[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Strain_c33 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", ep33c[i][j][k]);
				fprintf(fp,"%10.6f\n", ep33c[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Strain_c12 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", ep12c[i][j][k]);
				fprintf(fp,"%10.6f\n", ep12c[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Strain_c13 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", ep13c[i][j][k]);
				fprintf(fp,"%10.6f\n", ep13c[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Strain_c23 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", ep23c[i][j][k]);
				fprintf(fp,"%10.6f\n", ep23c[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_11 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", sig11[i][j][k]);
				fprintf(fp,"%10.6f\n", sig11[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_22 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", sig22[i][j][k]);
				fprintf(fp,"%10.6f\n", sig22[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_33 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", sig33[i][j][k]);
				fprintf(fp,"%10.6f\n", sig33[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_12 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", sig12[i][j][k]);
				fprintf(fp,"%10.6f\n", sig12[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_13 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", sig13[i][j][k]);
				fprintf(fp,"%10.6f\n", sig13[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_23 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", sig23[i][j][k]);
				fprintf(fp,"%10.6f\n", sig23[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Displacement_u1 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", u1[i][j][k]);
				fprintf(fp,"%10.6f\n", u1[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Displacement_u2 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", u2[i][j][k]);
				fprintf(fp,"%10.6f\n", u2[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Displacement_u3 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", u3[i][j][k]);
				fprintf(fp,"%10.6f\n", u3[i*ND*ND+j*ND+k]);
			}
		}
	}
	fclose(fp);
}
//************ Reading field data *****************************************
void datin(double *ch, int ND)
{
	FILE *datin0;//Stream pointer setting
	int  i, j, k;//integer
	double c00;//Average value of the field
	int nd=ND, ndm=ND-1, nd2=ND/2;

	datin0 = fopen("data.dat", "r");//Open the source file

	fscanf(datin0, "%lf", &time1);

	c00=0.0;//Initial value of field mean
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				fscanf(datin0, "%lf", &ch[i*ND*ND+j*ND+k]);//Read field data
				c00+=ch[i*ND*ND+j*ND+k];
			}
		}
	}
	c0=c00/nd/nd/nd;//Average value of the field
	printf("c0=  %f  \n", c0);

	fclose(datin0);
}
