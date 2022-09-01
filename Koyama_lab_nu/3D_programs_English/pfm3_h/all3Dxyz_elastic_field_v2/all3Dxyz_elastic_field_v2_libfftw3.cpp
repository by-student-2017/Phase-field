#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include <fftw3.h>
//g++ martensite_v2_libfftw3.cpp -lfftw3
//#include <mpi.h> //mpi version
//#include <fftw3-mpi.h> //mpi version

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
	//double ch[NDX][NDY][NDZ];		//phase field
	double eta_c[4][4];			//Lattice mismatch
	double cec[4][4][4][4];		//Elastic constant
	double sigma[4][4];			//Eigen stress (stress when elastically deformed by the amount of Eigen strain)
	//double ep11c[NDX][NDY][NDZ], ep22c[NDX][NDY][NDZ], ep33c[NDX][NDY][NDZ];	//Strain
	//double ep12c[NDX][NDY][NDZ], ep13c[NDX][NDY][NDZ], ep23c[NDX][NDY][NDZ];
	//double sig11[NDX][NDY][NDZ], sig22[NDX][NDY][NDZ], sig33[NDX][NDY][NDZ];	//Elastic stress
	//double sig12[NDX][NDY][NDZ], sig13[NDX][NDY][NDZ], sig23[NDX][NDY][NDZ];
	//double Estr[NDX][NDY][NDZ];	//Elastic strain energy density
	//double u1[NDX][NDY][NDZ], u2[NDX][NDY][NDZ], u3[NDX][NDY][NDZ];	//Displacement
	//double fld1[NDX][NDY][NDZ];

	//bag fix: move them from main program to here
	//double qrh1[NDX][NDY][NDZ], qih1[NDX][NDY][NDZ];	//Fourier transform of the field
	//double ec11[NDX][NDY][NDZ], ec22[NDX][NDY][NDZ], ec33[NDX][NDY][NDZ];	//Constrained strain array
	//double ec12[NDX][NDY][NDZ], ec13[NDX][NDY][NDZ], ec23[NDX][NDY][NDZ];

	double zcij(int i0, int j0, int k0, int iii, int jjj, int NDX, int NDY, int NDZ);//Coefficient calculation of elastic function (Fourier space)
	double zuij(int i0, int j0, int k0, int iii, int NDX, int NDY, int NDZ);//Displacement field coefficient calculation (Fourier space)

	void datin(double *ch, int NDX, int NDY, int NDZ);	//Subroutine for initial field reading
	void datsave(double *ch, double *Estr, 
			 double *ep11c, double *ep22c, double *ep33c, double *ep12c, double *ep13c, double *ep23c, 
			 double *sig11, double *sig22, double *sig33, double *sig12, double *sig13, double *sig23, 
			 double *u1, double *u2, double *u3, int NDX, int NDY, int NDZ);	//Data storage subroutine
	void datsave_paraview(double *ch, double *Estr, 
			 double *ep11c, double *ep22c, double *ep33c, double *ep12c, double *ep13c, double *ep23c, 
			 double *sig11, double *sig22, double *sig33, double *sig12, double *sig13, double *sig23, 
			 double *u1, double *u2, double *u3, int NDX, int NDY, int NDZ);	//Data storage subroutine

//******* main program ****************************************************
int main(void)
{
    int NDX, NDY, NDZ;
	int ndx, ndxm, ndx2;
	int ndy, ndym, ndy2;
	int ndz, ndzm, ndz2;
	
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
	// pfm1 series
	NDX      = int(data[0]);
	NDY      = int(data[1]);
	NDZ      = int(data[2]);
	//N        = int(data[3]);
	//GNP      = int(data[4]);
	//delt    = data[5];
	//L       = data[6];	// [nm]
	//time1max= data[7];
	//temp    = data[8];	// [K]
	//vm0     = data[9];	//L0=L00/RR/temp;
	//gamma0  = data[10];	//kapa_c=kapa_c0/b1/b1/RR/temp;
	//delta   = data[11];
	//K0      = data[12];
	//W0      = data[13];
	//amobi   = data[14];
	//M0      = data[15];
	//E0      = data[16];
	//Nstep   = int(data[17]);
	ep000   = data[18];
	el_mag  = data[19];
	c11     = data[20];
	c12     = data[21];
	c44     = data[22];
	printf("---------------------------------\n");
	
	//double ch[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	//a[z * ySize * xSize + y * xSize + x], a[i][j][k]=a[i*n*m + j*n + k]
	double *ch    = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));	//phase field
	//
	double *ep11c = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));	//Strain
	double *ep22c = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));
	double *ep33c = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));
	double *ep12c = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));
	double *ep13c = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));
	double *ep23c = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));
	//
	double *sig11 = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));	//Elastic stress
	double *sig22 = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));
	double *sig33 = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));
	double *sig12 = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));
	double *sig13 = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));
	double *sig23 = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));
	//
	double *Estr  = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));	//Elastic strain energy density
	//
	double *u1    = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));	//Displacement
	double *u2    = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));
	double *u3    = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));
	//
	double *qrh1  = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));	//Fourier transform of the field
	double *qih1  = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));
	//
	double *ec11  = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));	//Constrained strain array
	double *ec22  = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));
	double *ec33  = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));
	double *ec12  = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));
	double *ec13  = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));
	double *ec23  = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));
	//
	double *xr    = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));	//Real or imaginary array of Fourier transform
	double *xi    = (double *)malloc(sizeof(double)*( NDX*NDY*NDZ + NDY*NDZ + NDZ ));
	//
	const int fftsizex = NDX, fftsizey = NDY, fftsizez = NDZ;
	fftw_complex *in, *out; // in[i][0] for real, in[i][1] for imag.
	fftw_plan plan, iplan;
	in = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey*fftsizez);
	out = (fftw_complex *)fftw_malloc(sizeof(fftw_complex) * fftsizex*fftsizey*fftsizez);
	//plan = fftw_plan_dft_1d(fftsize, in, out, FFTW_FORWARD, FFTW_ESTIMATE); //FFT
	//fftw_execute(plan); //FFT
	//iplan = fftw_plan_dft_1d(fftsize, in, out, FFTW_BACKWARD, FFTW_ESTIMATE); //IFFT
	//fftw_execute(iplan); //IFFT
	//
	//After calc.: fftw_destroy_plan(plan); fftw_free(in); fftw_free(out);
	//
	
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
	//
	// Memo 1 (cec[i][j][k][l]=cec[i][j][l][k]=cec[j][i][k][l]=cec[k][l][i][j])
	// cec[1][1][1][1]=c11;  cec[1][1][2][2]=c12;  cec[1][1][3][3]=c13;  c1123=c1131=c1112=0, c2311=c3111=c1211=0
	// cec[2][2][1][1]=c21;  cec[2][2][2][2]=c22;  cec[2][2][3][3]=c23;  c2223=c2231=c2212=0, c2322=c3122=c1222=0
	// cec[3][3][1][1]=c31;  cec[3][3][2][2]=c32;  cec[3][3][3][3]=c33;  c3323=c2231=c2212=0, c2333=c3133=c1232=0
	// cec[2][3][2][3]=c44;  c2331=c2312=0
	// cec[3][1][3][1]=c55;  c3123=c3112=0
	// cec[1][2][1][2]=c66;  c1223=c1231=0
	//
	// Memo 2
	// isotropic : c11=c22=c33 = lambda + 2*mu
	//             c12=c21=c13=c31=c23=c32 = lambda
	//             c44=c55=c66 = mu
	//             other=0
	//             (relation: c11-c12=2*c44)
	// cubic     : c11=c22=c33, c12=c21=c13=c31=c23=c32, c44=c55=c66, other=0
	//             (2*c44-c11+c12)>0 -> soft <100> and min Y<100>
	//             (2*c44-c11+c12)<0 -> soft <111> and min Y<111>
	//             Y<100>=(c11+2*c12)*(c11-c12)/c11
	//             Y<111>=6*c44*(c11+2*c12)/(c11+2*c12+4*c44)
	//             elastic anisotropy parameter, A = 2*c44/(c11-c12)
	// Tetragonal: c11=c22=c33, c12=c21, c13=c31=c23=c32, c44=c55, c66, other=0
	// hexagonal : c11=c22=c33, c12=c21, c13=c31=c23=c32, c44=c55, c66=(c11-c12)/2, other=0
	//
	// Memo 3 (sigma ij = C ijkl * epsilon k l)
	// sigma11 = c1111*epsilon11 + c1122*epsilon22 + c1133*epsilon33
	// sigma22 = c1122*epsilon11 + c2222*epsilon22 + c2233*epsilon33
	// sigma33 = c1133*epsilon11 + c2233*epsilon33 + c3333*epsilon33
	// sigma23 = c2323*2.0*epsilon23
	// sigma31 = c3131*2.0*epsilon31
	// sigma12 = c1212*2.0*epsilon12

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

//*** initial field settings ***************
	//ini000();
	datin(ch, NDX, NDY, NDZ); //Read phase field

 	//gwinsize(INX,INY); ginit(1); gsetorg(0,0);// drow figure

//**** Start of calculation of elastic field analysis ******************************
 plan = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);	//For FFT
iplan = fftw_plan_dft_3d(fftsizex, fftsizey, fftsizez, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);	//For IFFT
//start:;

//**** Fourier transform of the field (phase field) ********************************
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				//xr[i][j][k]=ch[i][j][k]; xi[i][j][k]=0.;
				xr[i*NDY*NDZ+j*NDZ+k]=ch[i*NDY*NDZ+j*NDZ+k];
				xi[i*NDY*NDZ+j*NDZ+k]=0.0;
				//
				in[i*NDY*NDZ+j*NDZ+k][0] = xr[i*NDY*NDZ+j*NDZ+k]; //For FFT (real)
				in[i*NDY*NDZ+j*NDZ+k][1] = xi[i*NDY*NDZ+j*NDZ+k]; //For FFT (imag.)
			}
		}
	}
	fftw_execute(plan); //For FFT
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				xr[i*NDY*NDZ+j*NDZ+k] = out[i*NDY*NDZ+j*NDZ+k][0]; //For FFT (real)
				xi[i*NDY*NDZ+j*NDZ+k] = out[i*NDY*NDZ+j*NDZ+k][1]; //For FFT (imag.)
				//
				//qrh1[i][j][k]=xr[i][j][k]; qih1[i][j][k]=xi[i][j][k];
				qrh1[i*NDY*NDZ+j*NDZ+k]=xr[i*NDY*NDZ+j*NDZ+k];
				qih1[i*NDY*NDZ+j*NDZ+k]=xi[i*NDY*NDZ+j*NDZ+k];
			}
		}
	}
	//qrh1[0][0][0]=qih1[0][0][0]=0.0;
	qrh1[0]=qih1[0]=0.0;

//***** Calculation of total strain fluctuation *************************************
//--- ec11 ---
	iii=1; jjj=1; // ec11 = ec iii jjj
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				//xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				//xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xr[i*NDY*NDZ+j*NDZ+k]=qrh1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ)+qih1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ);
				xi[i*NDY*NDZ+j*NDZ+k]=qrh1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ)+qih1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ);
				//
				in[i*NDY*NDZ+j*NDZ+k][0] = xr[i*NDY*NDZ+j*NDZ+k]; //For IFFT (real)
				in[i*NDY*NDZ+j*NDZ+k][1] = xi[i*NDY*NDZ+j*NDZ+k]; //For IFFT (imag.)
			}
		}
	}
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				xr[i*NDY*NDZ+j*NDZ+k] = out[i*NDY*NDZ+j*NDZ+k][0]/(fftsizex*fftsizey*fftsizez); //For IFFT (real)
				//
				//ec11[i][j][k]=xr[i][j][k];
				ec11[i*NDY*NDZ+j*NDZ+k]=xr[i*NDY*NDZ+j*NDZ+k];
			}
		}
	}

//--- ec22 ---
	iii=2; jjj=2; // ec22 = ec iii jjj
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				//xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				//xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xr[i*NDY*NDZ+j*NDZ+k]=qrh1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ)+qih1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ);
				xi[i*NDY*NDZ+j*NDZ+k]=qrh1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ)+qih1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ);
				//
				in[i*NDY*NDZ+j*NDZ+k][0] = xr[i*NDY*NDZ+j*NDZ+k]; //For IFFT (real)
				in[i*NDY*NDZ+j*NDZ+k][1] = xi[i*NDY*NDZ+j*NDZ+k]; //For IFFT (imag.)
			}
		}
	}
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				xr[i*NDY*NDZ+j*NDZ+k] = out[i*NDY*NDZ+j*NDZ+k][0]/(fftsizex*fftsizey*fftsizez); //For IFFT (real)
				//
				//ec22[i][j][k]=xr[i][j][k];
				ec22[i*NDY*NDZ+j*NDZ+k]=xr[i*NDY*NDZ+j*NDZ+k];
			}
		}
	}

//--- ec33 ---
	iii=3; jjj=3; // ec33 = ec iii jjj
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				//xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				//xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xr[i*NDY*NDZ+j*NDZ+k]=qrh1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ)+qih1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ);
				xi[i*NDY*NDZ+j*NDZ+k]=qrh1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ)+qih1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ);
				//
				in[i*NDY*NDZ+j*NDZ+k][0] = xr[i*NDY*NDZ+j*NDZ+k]; //For IFFT (real)
				in[i*NDY*NDZ+j*NDZ+k][1] = xi[i*NDY*NDZ+j*NDZ+k]; //For IFFT (imag.)
			}
		}
	}
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				xr[i*NDY*NDZ+j*NDZ+k] = out[i*NDY*NDZ+j*NDZ+k][0]/(fftsizex*fftsizey*fftsizez); //For IFFT (real)
				//
				//ec33[i][j][k]=xr[i][j][k];
				ec33[i*NDY*NDZ+j*NDZ+k]=xr[i*NDY*NDZ+j*NDZ+k];
			}
		}
	}

//--- ec12 ---
	iii=1; jjj=2; // ec12 = ec iii jjj
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				//xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				//xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xr[i*NDY*NDZ+j*NDZ+k]=qrh1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ)+qih1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ);
				xi[i*NDY*NDZ+j*NDZ+k]=qrh1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ)+qih1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ);
				//
				in[i*NDY*NDZ+j*NDZ+k][0] = xr[i*NDY*NDZ+j*NDZ+k]; //For IFFT (real)
				in[i*NDY*NDZ+j*NDZ+k][1] = xi[i*NDY*NDZ+j*NDZ+k]; //For IFFT (imag.)
			}
		}
	}
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				xr[i*NDY*NDZ+j*NDZ+k] = out[i*NDY*NDZ+j*NDZ+k][0]/(fftsizex*fftsizey*fftsizez); //For IFFT (real)
				//
				//ec12[i][j][k]=xr[i][j][k];
				ec12[i*NDY*NDZ+j*NDZ+k]=xr[i*NDY*NDZ+j*NDZ+k];
			}
		}
	}

//--- ec13 ---
	iii=1; jjj=3; // ec13 = ec iii jjj
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				//xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				//xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xr[i*NDY*NDZ+j*NDZ+k]=qrh1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ)+qih1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ);
				xi[i*NDY*NDZ+j*NDZ+k]=qrh1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ)+qih1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ);
				//
				in[i*NDY*NDZ+j*NDZ+k][0] = xr[i*NDY*NDZ+j*NDZ+k]; //For IFFT (real)
				in[i*NDY*NDZ+j*NDZ+k][1] = xi[i*NDY*NDZ+j*NDZ+k]; //For IFFT (imag.)
			}
		}
	}
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				xr[i*NDY*NDZ+j*NDZ+k] = out[i*NDY*NDZ+j*NDZ+k][0]/(fftsizex*fftsizey*fftsizez); //For IFFT (real)
				//
				//ec13[i][j][k]=xr[i][j][k];
				ec13[i*NDY*NDZ+j*NDZ+k]=xr[i*NDY*NDZ+j*NDZ+k];
			}
		}
	}

//--- ec23 ---
	iii=2; jjj=3; // ec23 = ec iii jjj
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				//xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				//xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xr[i*NDY*NDZ+j*NDZ+k]=qrh1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ)+qih1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ);
				xi[i*NDY*NDZ+j*NDZ+k]=qrh1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ)+qih1[i*NDY*NDZ+j*NDZ+k]*zcij(i, j, k, iii, jjj, NDX, NDY, NDZ);
				//
				in[i*NDY*NDZ+j*NDZ+k][0] = xr[i*NDY*NDZ+j*NDZ+k]; //For IFFT (real)
				in[i*NDY*NDZ+j*NDZ+k][1] = xi[i*NDY*NDZ+j*NDZ+k]; //For IFFT (imag.)
			}
		}
	}
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				xr[i*NDY*NDZ+j*NDZ+k] = out[i*NDY*NDZ+j*NDZ+k][0]/(fftsizex*fftsizey*fftsizez); //For IFFT (real)
				//
				//ec23[i][j][k]=xr[i][j][k];
				ec23[i*NDY*NDZ+j*NDZ+k]=xr[i*NDY*NDZ+j*NDZ+k];
			}
		}
	}

//***** Calculation of total strain field, stress field, and elastic strain energy field *************************************
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				//Calculation of total strain field
				//ep11c[i][j][k]=ec11[i][j][k]+ep000*c0;
				//ep22c[i][j][k]=ec22[i][j][k]+ep000*c0;
				//ep33c[i][j][k]=ec33[i][j][k]+ep000*c0;
				//ep12c[i][j][k]=ec12[i][j][k];
				//ep13c[i][j][k]=ec13[i][j][k];
				//ep23c[i][j][k]=ec23[i][j][k];
				ep11c[i*NDY*NDZ+j*NDZ+k]=ec11[i*NDY*NDZ+j*NDZ+k]+ep000*c0;
				ep22c[i*NDY*NDZ+j*NDZ+k]=ec22[i*NDY*NDZ+j*NDZ+k]+ep000*c0;
				ep33c[i*NDY*NDZ+j*NDZ+k]=ec33[i*NDY*NDZ+j*NDZ+k]+ep000*c0;
				ep12c[i*NDY*NDZ+j*NDZ+k]=ec12[i*NDY*NDZ+j*NDZ+k];
				ep13c[i*NDY*NDZ+j*NDZ+k]=ec13[i*NDY*NDZ+j*NDZ+k];
				ep23c[i*NDY*NDZ+j*NDZ+k]=ec23[i*NDY*NDZ+j*NDZ+k];

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
				sig11[i*NDY*NDZ+j*NDZ+k]=c11*ec11[i*NDY*NDZ+j*NDZ+k]+c12*ec22[i*NDY*NDZ+j*NDZ+k]+c12*ec33[i*NDY*NDZ+j*NDZ+k]
											   -(c11+2.*c12)*ep000*(ch[i*NDY*NDZ+j*NDZ+k]-c0);
				sig22[i*NDY*NDZ+j*NDZ+k]=c12*ec11[i*NDY*NDZ+j*NDZ+k]+c11*ec22[i*NDY*NDZ+j*NDZ+k]+c12*ec33[i*NDY*NDZ+j*NDZ+k]
											   -(c11+2.*c12)*ep000*(ch[i*NDY*NDZ+j*NDZ+k]-c0);
				sig33[i*NDY*NDZ+j*NDZ+k]=c12*ec11[i*NDY*NDZ+j*NDZ+k]+c12*ec22[i*NDY*NDZ+j*NDZ+k]+c11*ec33[i*NDY*NDZ+j*NDZ+k]
											   -(c11+2.*c12)*ep000*(ch[i*NDY*NDZ+j*NDZ+k]-c0);
				sig12[i*NDY*NDZ+j*NDZ+k]=2.*c44*ec12[i*NDY*NDZ+j*NDZ+k];
				sig13[i*NDY*NDZ+j*NDZ+k]=2.*c44*ec13[i*NDY*NDZ+j*NDZ+k];
				sig23[i*NDY*NDZ+j*NDZ+k]=2.*c44*ec23[i*NDY*NDZ+j*NDZ+k];

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
				Estr1=1.5*(c11+2.0*c12)*ep000*ep000*ch[i*NDY*NDZ+j*NDZ+k]*ch[i*NDY*NDZ+j*NDZ+k];
				Estr2=-(c11+2.0*c12)*(ep11c[i*NDY*NDZ+j*NDZ+k]
					                 +ep22c[i*NDY*NDZ+j*NDZ+k]
					                 +ep33c[i*NDY*NDZ+j*NDZ+k])*ep000*ch[i*NDY*NDZ+j*NDZ+k];
				Estr3=0.5*c11*(ep11c[i*NDY*NDZ+j*NDZ+k]*ep11c[i*NDY*NDZ+j*NDZ+k]
							  +ep22c[i*NDY*NDZ+j*NDZ+k]*ep22c[i*NDY*NDZ+j*NDZ+k]
							  +ep33c[i*NDY*NDZ+j*NDZ+k]*ep33c[i*NDY*NDZ+j*NDZ+k]);
				Estr4=c12*(ep11c[i*NDY*NDZ+j*NDZ+k]*ep22c[i*NDY*NDZ+j*NDZ+k]
						  +ep11c[i*NDY*NDZ+j*NDZ+k]*ep33c[i*NDY*NDZ+j*NDZ+k]
						  +ep22c[i*NDY*NDZ+j*NDZ+k]*ep33c[i*NDY*NDZ+j*NDZ+k]);
				Estr5=2.*c44*(ep12c[i*NDY*NDZ+j*NDZ+k]*ep12c[i*NDY*NDZ+j*NDZ+k]
							 +ep13c[i*NDY*NDZ+j*NDZ+k]*ep13c[i*NDY*NDZ+j*NDZ+k]
							 +ep23c[i*NDY*NDZ+j*NDZ+k]*ep23c[i*NDY*NDZ+j*NDZ+k]);
				Estr[i*NDY*NDZ+j*NDZ+k]=Estr1+Estr2+Estr3+Estr4+Estr5;
			}
		}
	}


//***** Displacement field calculation *************************************
//--- u1 ---
	iii=1; // u1 = u iii
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				//xr[i][j][k]=-qrh1[i][j][k]*zuij(i, j, k, iii)-qih1[i][j][k]*zuij(i, j, k, iii);
				//xi[i][j][k]= qrh1[i][j][k]*zuij(i, j, k, iii)+qih1[i][j][k]*zuij(i, j, k, iii);
				xr[i*NDY*NDZ+j*NDZ+k]=-qrh1[i*NDY*NDZ+j*NDZ+k]*zuij(i, j, k, iii, NDX, NDY, NDZ)-qih1[i*NDY*NDZ+j*NDZ+k]*zuij(i, j, k, iii, NDX, NDY, NDZ);
				xi[i*NDY*NDZ+j*NDZ+k]= qrh1[i*NDY*NDZ+j*NDZ+k]*zuij(i, j, k, iii, NDX, NDY, NDZ)+qih1[i*NDY*NDZ+j*NDZ+k]*zuij(i, j, k, iii, NDX, NDY, NDZ);
				//
				in[i*NDY*NDZ+j*NDZ+k][0] = xr[i*NDY*NDZ+j*NDZ+k]; //For IFFT (real)
				in[i*NDY*NDZ+j*NDZ+k][1] = xi[i*NDY*NDZ+j*NDZ+k]; //For IFFT (imag.)
			}
	 	}
	}
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				xr[i*NDY*NDZ+j*NDZ+k] = out[i*NDY*NDZ+j*NDZ+k][0]/(fftsizex*fftsizey*fftsizez); //For IFFT (real)
				//
				//u1[i][j][k]=xr[i][j][k];
				u1[i*NDY*NDZ+j*NDZ+k]=xr[i*NDY*NDZ+j*NDZ+k];
			}
		}
	}

//--- u2 ---
	iii=2; // u2 = u iii
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				//xr[i][j][k]=-qrh1[i][j][k]*zuij(i, j, k, iii)-qih1[i][j][k]*zuij(i, j, k, iii);
				//xi[i][j][k]= qrh1[i][j][k]*zuij(i, j, k, iii)+qih1[i][j][k]*zuij(i, j, k, iii);
				xr[i*NDY*NDZ+j*NDZ+k]=-qrh1[i*NDY*NDZ+j*NDZ+k]*zuij(i, j, k, iii, NDX, NDY, NDZ)-qih1[i*NDY*NDZ+j*NDZ+k]*zuij(i, j, k, iii, NDX, NDY, NDZ);
				xi[i*NDY*NDZ+j*NDZ+k]= qrh1[i*NDY*NDZ+j*NDZ+k]*zuij(i, j, k, iii, NDX, NDY, NDZ)+qih1[i*NDY*NDZ+j*NDZ+k]*zuij(i, j, k, iii, NDX, NDY, NDZ);
				//
				in[i*NDY*NDZ+j*NDZ+k][0] = xr[i*NDY*NDZ+j*NDZ+k]; //For IFFT (real)
				in[i*NDY*NDZ+j*NDZ+k][1] = xi[i*NDY*NDZ+j*NDZ+k]; //For IFFT (imag.)
			}
	 	}
	}
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				xr[i*NDY*NDZ+j*NDZ+k] = out[i*NDY*NDZ+j*NDZ+k][0]/(fftsizex*fftsizey*fftsizez); //For IFFT (real)
				//
				//u2[i][j][k]=xr[i][j][k];
				u2[i*NDY*NDZ+j*NDZ+k]=xr[i*NDY*NDZ+j*NDZ+k];
			}
		}
	}

//--- u3 ---
	iii=3; // u3 = u iii
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				//xr[i][j][k]=-qrh1[i][j][k]*zuij(i, j, k, iii)-qih1[i][j][k]*zuij(i, j, k, iii);
				//xi[i][j][k]= qrh1[i][j][k]*zuij(i, j, k, iii)+qih1[i][j][k]*zuij(i, j, k, iii);
				xr[i*NDY*NDZ+j*NDZ+k]=-qrh1[i*NDY*NDZ+j*NDZ+k]*zuij(i, j, k, iii, NDX, NDY, NDZ)-qih1[i*NDY*NDZ+j*NDZ+k]*zuij(i, j, k, iii, NDX, NDY, NDZ);
				xi[i*NDY*NDZ+j*NDZ+k]= qrh1[i*NDY*NDZ+j*NDZ+k]*zuij(i, j, k, iii, NDX, NDY, NDZ)+qih1[i*NDY*NDZ+j*NDZ+k]*zuij(i, j, k, iii, NDX, NDY, NDZ);
				//
				in[i*NDY*NDZ+j*NDZ+k][0] = xr[i*NDY*NDZ+j*NDZ+k]; //For IFFT (real)
				in[i*NDY*NDZ+j*NDZ+k][1] = xi[i*NDY*NDZ+j*NDZ+k]; //For IFFT (imag.)
			}
	 	}
	}
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				xr[i*NDY*NDZ+j*NDZ+k] = out[i*NDY*NDZ+j*NDZ+k][0]/(fftsizex*fftsizey*fftsizez); //For IFFT (real)
				//
				//u3[i][j][k]=xr[i][j][k];
				u3[i*NDY*NDZ+j*NDZ+k]=xr[i*NDY*NDZ+j*NDZ+k];
			}
		}
	}

//***** Save all numerical data ********************************************
	//datsave(ch, Estr, ep11c, ep22c, ep33c, ep12c, ep13c, ep23c, 
	//		 sig11, sig22, sig33, sig12, sig13, sig23, u1, u2, u3, NDX, NDY, NDZ);			//save data
	datsave_paraview(ch, Estr, ep11c, ep22c, ep33c, ep12c, ep13c, ep23c, 
			 sig11, sig22, sig33, sig12, sig13, sig23, u1, u2, u3, NDX, NDY, NDZ);	//save

	//if(keypress()){return 0;}//wait key input
	if(plan) fftw_destroy_plan(plan);		//For FFT
	if(iplan) fftw_destroy_plan(iplan);	//For IFFT
	fftw_free(in); fftw_free(out);		// For FFT and IFFT
	return 0;

}

//*** Zcij [eq.(5.26) or eq.(II 3.5)] ****************************************
double zcij(int i0, int j0, int k0, int iii, int jjj, int NDX, int NDY, int NDZ)
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
	int ndx=NDX, ndxm=NDX-1, ndx2=NDX/2;
	int ndy=NDY, ndym=NDY-1, ndy2=NDY/2;
	int ndz=NDZ, ndzm=NDZ-1, ndz2=NDZ/2;

	if(i0<=ndx2-1){ii=i0;}  if(i0>=ndx2){ii=i0-ndx;}
	if(j0<=ndy2-1){jj=j0;}  if(j0>=ndy2){jj=j0-ndy;}
	if(k0<=ndz2-1){kk=k0;}  if(k0>=ndz2){kk=k0-ndz;}
	alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj+(double)kk*(double)kk);
	if(alnn==0.){alnn=1.;}
	nec[1]=nx=(double)ii/alnn;	// n is the unit vector in the k direction, n=k/|k|
	nec[2]=ny=(double)jj/alnn;	// n is the unit vector in the k direction, n=k/|k|
	nec[3]=nz=(double)kk/alnn;	// n is the unit vector in the k direction, n=k/|k|

	//for(i=1;i<=3;i++){ for(j=1;j<=3;j++){ om[i][j]=0.; } }

	// C[i][k][j][l]*n[j]*n[l], C: elastic modulus, n: unit vector
	// C = cec
	a11=cec[1][1][1][1]*nx*nx+cec[1][2][1][2]*ny*ny+cec[1][3][1][3]*nz*nz;
	a22=cec[1][2][1][2]*nx*nx+cec[2][2][2][2]*ny*ny+cec[2][3][2][3]*nz*nz;
	a33=cec[3][1][3][1]*nx*nx+cec[2][3][2][3]*ny*ny+cec[3][3][3][3]*nz*nz;
	a12=(cec[1][1][2][2]+cec[1][2][1][2])*nx*ny;
	a23=(cec[2][2][3][3]+cec[2][3][2][3])*ny*nz;
	a31=(cec[3][3][1][1]+cec[3][1][3][1])*nx*nz;
	a21=a12;
	a32=a23;
	a13=a31;

	// cofactor
	b11=a22*a33-a23*a32;
	b22=a11*a33-a13*a31;
	b33=a11*a22-a12*a21;
	b12=-(a21*a33-a23*a31);
	b23=-(a11*a32-a12*a31);
	b31=-(a22*a31-a21*a32);
	b21=b12;
	b32=b23;
	b13=b31;

	// det (C[i][k][j][l]*n[j]*n[l])
	det1=a11*a22*a33+a12*a23*a31+a13*a32*a21-a13*a31*a22-a11*a23*a32-a33*a12*a21;
	if(det1==0.){det1=1.;}

	// inverse matrix
	om[1][1]=b11/det1;
	om[2][2]=b22/det1;
	om[3][3]=b33/det1;
	om[1][2]=om[2][1]=b12/det1;
	om[2][3]=om[3][2]=b23/det1;
	om[3][1]=om[1][3]=b31/det1;

	// sigma[i][j] = cec[i][j][k][l]*eta_c[k][l]*(Kronecker delta[k][l])
	// sigma: Eigen stress
	zij=0.0;
	for(m=1;m<=3;m++) {
		for(n=1;n<=3;n++) {
			//zij=zij+0.5*( om[m][iii]*nec[n]*nec[jjj] + om[m][jjj]*nec[n]*nec[iii] )*sigma[m][n]; // eq.(5.26) or eq.(II 3.5)
    		zij=zij+0.5*(sigma[m][n]*nec[jjj]*nec[n]*om[m][iii]
				  		+sigma[m][n]*nec[iii]*nec[n]*om[m][jjj]);
		}
	}
	return(zij);
}

//*** Zuij [eq.(5.30) or eq.(II 3.9)] ****************************************
double zuij(int i0, int j0, int k0, int iii, int NDX, int NDY, int NDZ)
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
	int ndx=NDX, ndxm=NDX-1, ndx2=NDX/2;
	int ndy=NDY, ndym=NDY-1, ndy2=NDY/2;
	int ndz=NDZ, ndzm=NDZ-1, ndz2=NDZ/2;

	if(i0<=ndx2-1){ii=i0;}  if(i0>=ndx2){ii=i0-ndx;}
	if(j0<=ndy2-1){jj=j0;}  if(j0>=ndy2){jj=j0-ndy;}
	if(k0<=ndz2-1){kk=k0;}  if(k0>=ndz2){kk=k0-ndz;}
	//alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj); // miss ?
	alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj+(double)kk*(double)kk);
	if(alnn==0.){alnn=1.;}
	nec[1]=nx=(double)ii/alnn;	// n is the unit vector in the k direction, n=k/|k|
	nec[2]=ny=(double)jj/alnn;	// n is the unit vector in the k direction, n=k/|k|
	nec[3]=nz=(double)kk/alnn;	// n is the unit vector in the k direction, n=k/|k|

	//for(i=1;i<=3;i++){ for(j=1;j<=3;j++){ om[i][j]=0.; } }

	// C[i][k][j][l]*n[j]*n[l], C: elastic modulus, n: unit vector
	// C = cec
	a11=cec[1][1][1][1]*nx*nx+cec[1][2][1][2]*ny*ny+cec[1][3][1][3]*nz*nz;
	a22=cec[1][2][1][2]*nx*nx+cec[2][2][2][2]*ny*ny+cec[2][3][2][3]*nz*nz;
	a33=cec[3][1][3][1]*nx*nx+cec[2][3][2][3]*ny*ny+cec[3][3][3][3]*nz*nz;
	a12=(cec[1][1][2][2]+cec[1][2][1][2])*nx*ny;
	a23=(cec[2][2][3][3]+cec[2][3][2][3])*ny*nz;
	a31=(cec[3][3][1][1]+cec[3][1][3][1])*nx*nz;
	a21=a12;
	a32=a23;
	a13=a31;

	// cofactor
	b11=a22*a33-a23*a32;
	b22=a11*a33-a13*a31;
	b33=a11*a22-a12*a21;
	b12=-(a21*a33-a23*a31);
	b23=-(a11*a32-a12*a31);
	b31=-(a22*a31-a21*a32);
	b21=b12;
	b32=b23;
	b13=b31;

	// det (C[i][k][j][l]*n[j]*n[l])
	det1=a11*a22*a33+a12*a23*a31+a13*a32*a21-a13*a31*a22-a11*a23*a32-a33*a12*a21;
	if(det1==0.){det1=1.;}

	// inverse matrix
	om[1][1]=b11/det1;
	om[2][2]=b22/det1;
	om[3][3]=b33/det1;
	om[1][2]=om[2][1]=b12/det1;
	om[2][3]=om[3][2]=b23/det1;
	om[3][1]=om[1][3]=b31/det1;

	// sigma[i][j] = cec[i][j][k][l]*eta_c[k][l]*(Kronecker delta[k][l])
	// sigma: Eigen stress
	zij=0.;
	for(m=1;m<=3;m++) {
		for(n=1;n<=3;n++) {
			//zij=zij-(om[m][iii]/alnn)*sigma[m][n]*nec[n]; // eq.(5.30) or eq.(II 3.9)
    		zij=zij+sigma[m][n]*nec[n]*om[m][iii];
		}
	}
	zij=-zij/alnn; // alnn=|k}
	return(zij);
}

//************ Saving data of various places *******************************
void datsave(double *ch, double *Estr, 
			 double *ep11c, double *ep22c, double *ep33c, double *ep12c, double *ep13c, double *ep23c, 
			 double *sig11, double *sig22, double *sig33, double *sig12, double *sig13, double *sig23, 
			 double *u1, double *u2, double *u3, int NDX, int NDY, int NDZ)
{
	FILE		*stream;
	int 		i, j, k;
	int ndx=NDX, ndxm=NDX-1, ndx2=NDX/2;
	int ndy=NDY, ndym=NDY-1, ndy2=NDY/2;
	int ndz=NDZ, ndzm=NDZ-1, ndz2=NDZ/2;

	stream = fopen("el_field.dat", "w");

	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				fprintf(stream, "%e  ", ch[i*NDY*NDZ+j*NDZ+k]);
				fprintf(stream, "%e  ", Estr[i*NDY*NDZ+j*NDZ+k]);
				fprintf(stream, "%e  ", ep11c[i*NDY*NDZ+j*NDZ+k]);
				fprintf(stream, "%e  ", ep22c[i*NDY*NDZ+j*NDZ+k]);
				fprintf(stream, "%e  ", ep33c[i*NDY*NDZ+j*NDZ+k]);
				fprintf(stream, "%e  ", ep12c[i*NDY*NDZ+j*NDZ+k]);
				fprintf(stream, "%e  ", ep13c[i*NDY*NDZ+j*NDZ+k]);
				fprintf(stream, "%e  ", ep23c[i*NDY*NDZ+j*NDZ+k]);
				fprintf(stream, "%e  ", sig11[i*NDY*NDZ+j*NDZ+k]);
				fprintf(stream, "%e  ", sig22[i*NDY*NDZ+j*NDZ+k]);
				fprintf(stream, "%e  ", sig33[i*NDY*NDZ+j*NDZ+k]);
				fprintf(stream, "%e  ", sig12[i*NDY*NDZ+j*NDZ+k]);
				fprintf(stream, "%e  ", sig13[i*NDY*NDZ+j*NDZ+k]);
				fprintf(stream, "%e  ", sig23[i*NDY*NDZ+j*NDZ+k]);
				fprintf(stream, "%e  ", u1[i*NDY*NDZ+j*NDZ+k]);
				fprintf(stream, "%e  ", u2[i*NDY*NDZ+j*NDZ+k]);
				fprintf(stream, "%e  ", u3[i*NDY*NDZ+j*NDZ+k]);
			}
		}
	}
	fprintf(stream, "\n");
	fclose(stream);
}

void datsave_paraview(double *ch, double *Estr, 
			 double *ep11c, double *ep22c, double *ep33c, double *ep12c, double *ep13c, double *ep23c, 
			 double *sig11, double *sig22, double *sig33, double *sig12, double *sig13, double *sig23, 
			 double *u1, double *u2, double *u3, int NDX, int NDY, int NDZ)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j, k;
	int ndx=NDX, ndxm=NDX-1, ndx2=NDX/2;
	int ndy=NDY, ndym=NDY-1, ndy2=NDY/2;
	int ndz=NDZ, ndzm=NDZ-1, ndz2=NDZ/2;
	
	//iout = iout + 1;
	sprintf(fName,"el_result%06d.vtk",iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndxm+1),(ndym+1),(ndzm+1));
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO %f %f %f \n",float(ndxm/ndxm),float(ndym/ndxm),float(ndzm/ndxm));
	fprintf(fp,"POINT_DATA %16d \n",((ndxm+1)*(ndym+1)*(ndzm+1)));
	fprintf(fp,"SCALARS concentration float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndzm;k++){
		for(j=0;j<=ndym;j++){
			for(i=0;i<=ndxm;i++){
				//fprintf(fp,"%10.6f\n", ch[i][j][k]);
				fprintf(fp,"%10.6f\n", ch[i*NDY*NDZ+j*NDZ+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_strain_energy_density float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndzm;k++){
		for(j=0;j<=ndym;j++){
			for(i=0;i<=ndxm;i++){
				//fprintf(fp,"%10.6f\n", Estr[i][j][k]);
				fprintf(fp,"%10.6f\n", Estr[i*NDY*NDZ+j*NDZ+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Strain_c11 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndzm;k++){
		for(j=0;j<=ndym;j++){
			for(i=0;i<=ndxm;i++){
				//fprintf(fp,"%10.6f\n", ep11c[i][j][k]);
				fprintf(fp,"%10.6f\n", ep11c[i*NDY*NDZ+j*NDZ+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Strain_c22 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndzm;k++){
		for(j=0;j<=ndym;j++){
			for(i=0;i<=ndxm;i++){
				//fprintf(fp,"%10.6f\n", ep22c[i][j][k]);
				fprintf(fp,"%10.6f\n", ep22c[i*NDY*NDZ+j*NDZ+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Strain_c33 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndzm;k++){
		for(j=0;j<=ndym;j++){
			for(i=0;i<=ndxm;i++){
				//fprintf(fp,"%10.6f\n", ep33c[i][j][k]);
				fprintf(fp,"%10.6f\n", ep33c[i*NDY*NDZ+j*NDZ+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Strain_c12 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndzm;k++){
		for(j=0;j<=ndym;j++){
			for(i=0;i<=ndxm;i++){
				//fprintf(fp,"%10.6f\n", ep12c[i][j][k]);
				fprintf(fp,"%10.6f\n", ep12c[i*NDY*NDZ+j*NDZ+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Strain_c13 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndzm;k++){
		for(j=0;j<=ndym;j++){
			for(i=0;i<=ndxm;i++){
				//fprintf(fp,"%10.6f\n", ep13c[i][j][k]);
				fprintf(fp,"%10.6f\n", ep13c[i*NDY*NDZ+j*NDZ+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Strain_c23 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndzm;k++){
		for(j=0;j<=ndym;j++){
			for(i=0;i<=ndxm;i++){
				//fprintf(fp,"%10.6f\n", ep23c[i][j][k]);
				fprintf(fp,"%10.6f\n", ep23c[i*NDY*NDZ+j*NDZ+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_11 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndzm;k++){
		for(j=0;j<=ndym;j++){
			for(i=0;i<=ndxm;i++){
				//fprintf(fp,"%10.6f\n", sig11[i][j][k]);
				fprintf(fp,"%10.6f\n", sig11[i*NDY*NDZ+j*NDZ+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_22 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndzm;k++){
		for(j=0;j<=ndym;j++){
			for(i=0;i<=ndxm;i++){
				//fprintf(fp,"%10.6f\n", sig22[i][j][k]);
				fprintf(fp,"%10.6f\n", sig22[i*NDY*NDZ+j*NDZ+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_33 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndzm;k++){
		for(j=0;j<=ndym;j++){
			for(i=0;i<=ndxm;i++){
				//fprintf(fp,"%10.6f\n", sig33[i][j][k]);
				fprintf(fp,"%10.6f\n", sig33[i*NDY*NDZ+j*NDZ+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_12 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndzm;k++){
		for(j=0;j<=ndym;j++){
			for(i=0;i<=ndxm;i++){
				//fprintf(fp,"%10.6f\n", sig12[i][j][k]);
				fprintf(fp,"%10.6f\n", sig12[i*NDY*NDZ+j*NDZ+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_13 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndzm;k++){
		for(j=0;j<=ndym;j++){
			for(i=0;i<=ndxm;i++){
				//fprintf(fp,"%10.6f\n", sig13[i][j][k]);
				fprintf(fp,"%10.6f\n", sig13[i*NDY*NDZ+j*NDZ+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_23 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndzm;k++){
		for(j=0;j<=ndym;j++){
			for(i=0;i<=ndxm;i++){
				//fprintf(fp,"%10.6f\n", sig23[i][j][k]);
				fprintf(fp,"%10.6f\n", sig23[i*NDY*NDZ+j*NDZ+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Displacement_u1 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndzm;k++){
		for(j=0;j<=ndym;j++){
			for(i=0;i<=ndxm;i++){
				//fprintf(fp,"%10.6f\n", u1[i][j][k]);
				fprintf(fp,"%10.6f\n", u1[i*NDY*NDZ+j*NDZ+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Displacement_u2 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndzm;k++){
		for(j=0;j<=ndym;j++){
			for(i=0;i<=ndxm;i++){
				//fprintf(fp,"%10.6f\n", u2[i][j][k]);
				fprintf(fp,"%10.6f\n", u2[i*NDY*NDZ+j*NDZ+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Displacement_u3 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndzm;k++){
		for(j=0;j<=ndym;j++){
			for(i=0;i<=ndxm;i++){
				//fprintf(fp,"%10.6f\n", u3[i][j][k]);
				fprintf(fp,"%10.6f\n", u3[i*NDY*NDZ+j*NDZ+k]);
			}
		}
	}
	fclose(fp);
}
//************ Reading field data *****************************************
void datin(double *ch, int NDX, int NDY, int NDZ)
{
	FILE *datin0;//Stream pointer setting
	int  i, j, k;//integer
	double c00;//Average value of the field
	int ndx=NDX, ndxm=NDX-1, ndx2=NDX/2;
	int ndy=NDY, ndym=NDY-1, ndy2=NDY/2;
	int ndz=NDZ, ndzm=NDZ-1, ndz2=NDZ/2;
	int rndxm, rndym, rndzm;

	datin0 = fopen("data.dat", "r");//Open the source file

	fscanf(datin0, "%d %d %d", &rndxm, &rndym, &rndzm);
	if (ndxm != rndxm || ndym != rndym || ndzm != rndzm){
		printf("data size is mismatch \n");
		printf("Please, change NDX, NDY or NDZ in parameters.txt \n");
	}

	fscanf(datin0, "%lf", &time1);

	c00=0.0;//Initial value of field mean
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				fscanf(datin0, "%lf", &ch[i*NDY*NDZ+j*NDZ+k]);//Read field data
				c00+=ch[i*NDY*NDZ+j*NDZ+k];
			}
		}
	}
	c0=c00/ndx/ndy/ndz;//Average value of the field
	printf("c0=  %f  \n", c0);

	fclose(datin0);
}
