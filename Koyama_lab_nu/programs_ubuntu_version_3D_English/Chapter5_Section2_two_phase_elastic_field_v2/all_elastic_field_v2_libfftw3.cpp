//Program to read phase field data (ph.dat) and calculate elastic field.
//Each elastic field is collectively saved in el_field.dat.
//At the same time, a two-dimensional image of each elastic field is also saved.

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#include <fftw3.h>
//g++ martensite_v2_libfftw3.cpp -lfftw3
//#include <mpi.h> //mpi version
//#include <fftw3-mpi.h> //mpi version

#define DRND(x) ((double)(x)/RAND_MAX*rand())//Random number function setting
//#define ND 128	//The number of divisions per side of the computational domain in
					// difference computation (a power of 2 due to the use of fast Fourier transform)
//#define IG 7		//2^IG=ND

	//int nd=ND, ndm=ND-1; 		//Define the number of difference divisions
								// (number of difference blocks) on one side of the calculation area, ND-1
	//int nd2=ND/2;				//Define ND/2: used in Fast Fourier Transform
	//int ig=IG;				//2^ig=ND
	double PI=3.14159;			//Pi
	double rr=8.3145;			//gas constant
	double time1;				//Number of calculation counts (proportional to time)
	int iout=-1;

	double c0;					//Phase field mean
	//double ch[ND][ND];		//field (phase field)
	double eta_c[4][4];			//lattice mismatch
	double cec[4][4][4][4];		//elastic constant
	double sigma[4][4];			//Eigen stress (stress when elastically deformed by Eigen strain)
	//double ep11c[ND][ND][ND], ep22c[ND][ND][ND], ep12c[ND][ND][ND];					//total strain
	//double sig11[ND][ND][ND], sig22[ND][ND][ND], sig33[ND][ND][ND], sig12[ND][ND][ND];	//elastic stress
	//double Estr[ND][ND][ND];		//elastic strain energy
	//double u1[ND][ND][ND], u2[ND][ND][ND], u3[ND][ND][ND];	//displacement
	//double fld1[ND][ND][ND];		//Work matrix for delivery of fields

	double zcij(int i0, int j0, int k0, int iii, int jjj, int ND);//Coefficient calculation of elastic energy function (Fourier space)
	double zuij(int i0, int j0, int k0, int iii, int ND);//Coefficient calculation of displacement field (Fourier space)

	void datin(double *ch, int ND);	//Initial field loading subroutine
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
	int ndm, ig;
	
	//int i, j, k, l, ii, jj, kk, iii, jjj, ief;		//integer
	int i, j, k, l, iii, jjj;							//integer
	//int ip, im, jp, jm;								//integer

	//double c;	//場（フェーズフィールド）
	//double qrh1[ND][ND][ND], qih1[ND][ND][ND];				//Fourier transform of the field
	//double ec11[ND][ND][ND], ec22[ND][ND][ND], ec12[ND][ND][ND];	//restraint strain array

	//double al, temp, delt;							//Side length of calculation domain, temperature, time step
	//double time1max;									//Maximum calculation count (used to stop calculation)
	//double b1, vm0, atom_n;							//Side length of difference block, molar volume, number of atoms in unit cell
	//double a1;										//Lattice parameter of Fe

	double el_mag;										//Working variables for elastic modulus
	double c11, c44, c12; 								//elastic modulus
	double ep000;										//lattice mismatch
	//double epc11, epc22, epc33, epc12, epc13, epc23;	//constraint strain
	//double ep011, ep022, ep033, ep012, ep013, ep023;	//eigen strain
	//double dep011, dep022, dep033, dep012, dep023, dep013;
	//double epc011, epc022, epc033, epc012, epc023, epc013;
	//double ef1, ef2, ef11, ef22, ef12;				//elastic function

	//double Estr1, Estr2, Estr3, Estr4, Estr5;			//Working variables for elastic strain energy calculation
	double Estr1, Estr2, Estr3, Estr4, Estr5;			//Working variables for elastic strain energy calculation

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
	// pfm1 series
	ND      = int(data[0]);
	//c2a     = data[1];
	//delt    = data[2];	// timestep
	//temp    = data[3];	// [K]
	//al      = data[4];	// [nm]
	//cmob22  = data[5];
	//om_12e  = data[6];
	//kapa_c2c= data[7];
	//time1max= data[8];
	//Nstep   = int(data[9]);
	ep000   = data[10];
	el_mag  = data[11];
	c11     = data[12];
	c12     = data[13];
	c44     = data[14];
	printf("---------------------------------\n");
	//
	IG=int(log2(ND));
	ndm=ND-1;				//define ND-1
	ig=IG;					//2^ig=ND
	//
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
	double *xr    = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//array of real and imaginary parts of the Fourier transform
	double *xi    = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//array of real and imaginary parts of the Fourier transform
	//
	double *ec11  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Constrained strain array
	double *ec22  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ec33  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ec12  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ec13  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ec23  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
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
	//temp=900.;			//temperature (K)
	//al=100.*1.0E-09;		//Computational area (m)
	//b1=al/nd;				//Diff block length

	//time1=0.0;			//Initial calculation count setting
	//time1max=1.0+1.0e+07;	//Setting the maximum calculation count

	//a1=2.8664E-10; 		//Lattice parameter of Fe(bcc) (experimental data)
	//atom_n=2.;	vm0=6.02E23*a1*a1*a1/atom_n;	//Calculation of molar volume (for bcc)

//*** Setting the lattice mismatch ***
	//ep000=0.05;
	eta_c[1][1]=eta_c[2][2]=eta_c[3][3]=ep000; 			//lattice mismatch
	eta_c[1][2]=eta_c[2][1]=eta_c[1][3]=eta_c[3][1]=eta_c[2][3]=eta_c[3][2]=0.0;

//***** Elastic constant of Fe(bcc) ****************************
  //el_mag=1.0E+11/1.0E+09;
  c11=2.33*el_mag;
  c12=1.35*el_mag;
  c44=1.18*el_mag;

	for(i=0;i<=3;i++){
		for(j=0;j<=3;j++){
			for(k=0;k<=3;k++){
				for(l=0;l<=3;l++){
					cec[i][j][k][l]=0.0;//Initialization of Elastic Constant Array
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

//--- Eigen stress (stress when elastically deformed by Eigen strain)）--------------
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

//***Setting up sin and cos tables, bit reversal tables, and initial fields ***************

	//ini000();
	datin(ch, ND); //Read phase field

//**** Calculation start of elastic field analysis ******************************
 plan = fftw_plan_dft_3d(fftsize, fftsize, fftsize, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);	//For FFT
iplan = fftw_plan_dft_3d(fftsize, fftsize, fftsize, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);	//For IFFT
//start: ;

//**** Fourier transform of field (phase field) ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=ch[i][j][k]; xi[i][j][k]=0.;
				xr[i*ND*ND+j*ND+k]=ch[i*ND*ND+j*ND+k];
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
				//qrh1[i][j][K]=xr[i][j][K]; qih1[i][j][K]=xi[i][j][K];
				qrh1[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
				qih1[i*ND*ND+j*ND+k]=xi[i*ND*ND+j*ND+k];
			}
		}
	}
	//qrh1[0][0]=qih1[0][0]=0.;
	qrh1[0]=qih1[0]=0.0;

//***** Calculation of total strain variation *************************************
//--- ec11 calculator ---
	iii=1; jjj=1; // ec11 = ec iii jjj
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][K]=qrh1[i][j][K]*zcij(i, j, k, iii, jjj)+qih1[i][j][K]*zcij(i, j, k, iii, jjj);
				//xi[i][j][K]=qrh1[i][j][K]*zcij(i, j, k, iii, jjj)+qih1[i][j][K]*zcij(i, j, k, iii, jjj);
				xr[i*ND*ND+j*ND+k]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)
								  +qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
				xi[i*ND*ND+j*ND+k]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)
								  +qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
				//
				in[i*ND*ND+j*ND+k][0] = xr[i*ND*ND+j*ND+k]; //For IFFT (real)
				in[i*ND*ND+j*ND+k][1] = xi[i*ND*ND+j*ND+k]; //For IFFT (imag.)
			}
		}
	}
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k] = out[i*ND*ND+j*ND+k][0]/(fftsize*fftsize*fftsize); //For IFFT (real)
				//
				//ec11[i][j][k]=xr[i][j][k];
				ec11[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- ec22 calculator ---
	iii=2; jjj=2; // ec22 = ec iii jjj
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				//xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xr[i*ND*ND+j*ND+k]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)
								  +qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
				xi[i*ND*ND+j*ND+k]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)
								  +qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
				//
				in[i*ND*ND+j*ND+k][0] = xr[i*ND*ND+j*ND+k]; //For IFFT (real)
				in[i*ND*ND+j*ND+k][1] = xi[i*ND*ND+j*ND+k]; //For IFFT (imag.)
			}
		}
	}
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k] = out[i*ND*ND+j*ND+k][0]/(fftsize*fftsize*fftsize); //For IFFT (real)
				//
				//ec22[i][j][k]=xr[i][j][k];
				ec22[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- ec33 calculator ---
	iii=3; jjj=3; // ec33 = ec iii jjj
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				//xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xr[i*ND*ND+j*ND+k]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)
								  +qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
				xi[i*ND*ND+j*ND+k]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)
								  +qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
				//
				in[i*ND*ND+j*ND+k][0] = xr[i*ND*ND+j*ND+k]; //For IFFT (real)
				in[i*ND*ND+j*ND+k][1] = xi[i*ND*ND+j*ND+k]; //For IFFT (imag.)
			}
		}
	}
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k] = out[i*ND*ND+j*ND+k][0]/(fftsize*fftsize*fftsize); //For IFFT (real)
				//
				//ec33[i][j][k]=xr[i][j][k];
				ec33[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- ec12 calculator ---
	iii=1; jjj=2; // ec12 = ec iii jjj
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				//xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xr[i*ND+j]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)
						  +qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
				xi[i*ND+j]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)
						  +qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
				//
				in[i*ND*ND+j*ND+k][0] = xr[i*ND*ND+j*ND+k]; //For IFFT (real)
				in[i*ND*ND+j*ND+k][1] = xi[i*ND*ND+j*ND+k]; //For IFFT (imag.)
			}
		}
	}
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k] = out[i*ND*ND+j*ND+k][0]/(fftsize*fftsize*fftsize); //For IFFT (real)
				//
				//ec12[i][j][k]=xr[i][j][k];
				ec12[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- ec13 calculator ---
	iii=1; jjj=3; // ec13 = ec iii jjj
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				//xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xr[i*ND+j]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)
						  +qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
				xi[i*ND+j]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)
						  +qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
				//
				in[i*ND*ND+j*ND+k][0] = xr[i*ND*ND+j*ND+k]; //For IFFT (real)
				in[i*ND*ND+j*ND+k][1] = xi[i*ND*ND+j*ND+k]; //For IFFT (imag.)
			}
		}
	}
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k] = out[i*ND*ND+j*ND+k][0]/(fftsize*fftsize*fftsize); //For IFFT (real)
				//
				//ec13[i][j][k]=xr[i][j][k];
				ec13[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- ec23 calculator ---
	iii=2; jjj=3; // ec23 = ec iii jjj
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				//xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xr[i*ND+j]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)
						  +qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
				xi[i*ND+j]=qrh1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND)
						  +qih1[i*ND*ND+j*ND+k]*zcij(i, j, k, iii, jjj, ND);
				//
				in[i*ND*ND+j*ND+k][0] = xr[i*ND*ND+j*ND+k]; //For IFFT (real)
				in[i*ND*ND+j*ND+k][1] = xi[i*ND*ND+j*ND+k]; //For IFFT (imag.)
			}
		}
	}
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k] = out[i*ND*ND+j*ND+k][0]/(fftsize*fftsize*fftsize); //For IFFT (real)
				//
				//ec23[i][j][k]=xr[i][j][k];
				ec23[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//***** Total strain field, stress field, and elastic strain energy field calculation *************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//Calculation of total strain field  [equation (5.1-5.25)]
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

				//Calculation of elastic stress field [equation (5.27)]
				// Hooke's law: sigma ij = C ijkl * (epsilon_c kl - epsilon_0 kl)
				// epsilon_0 kl = eta ij * phase-field
				// Vegard's law: epsilon = eta*(c-c0)
				//sig11[i][j][k]=c11*ec11[i][j][k]+c12*ec22[i][j][k]+c12*ec33[i][j][k]
				//							   -(c11+2.*c12)*ep000*(ch[i][j][k]-c0);
				//sig22[i][j][k]=c12*ec11[i][j][k]+c11*ec22[i][j][k]+c12*ec33[i][j][k]
				//							   -(c11+2.*c12)*ep000*(ch[i][j][k]-c0);
				//sig33[i][j][k]=c12*ec11[i][j][k]+c12*ec22[i][j][k]+c11*ec33[i][j][k]
				//							   -(c11+2.*c12)*ep000*(ch[i][j][k]-c0);
				//sig12[i][j][k]=2.*c44*ec12[i][j][k];
				//sig13[i][j][k]=2.*c44*ec13[i][j][k];
				//sig23[i][j][k]=2.*c44*ec23[i][j][k];
				//
				sig11[i*ND*ND+j*ND+k]=c11*ec11[i*ND*ND+j*ND+k]+c12*ec22[i*ND*ND+j*ND+k]+c12*ec33[i*ND*ND+j*ND+k]
											   -(c11+2.*c12)*ep000*(ch[i*ND*ND+j*ND+k]-c0);
				sig22[i*ND*ND+j*ND+k]=c12*ec11[i*ND*ND+j*ND+k]+c11*ec22[i*ND*ND+j*ND+k]+c12*ec33[i*ND*ND+j*ND+k]
											   -(c11+2.*c12)*ep000*(ch[i*ND*ND+j*ND+k]-c0);
				sig33[i*ND*ND+j*ND+k]=c12*ec11[i*ND*ND+j*ND+k]+c12*ec22[i*ND*ND+j*ND+k]+c11*ec33[i*ND*ND+j*ND+k]
											   -(c11+2.*c12)*ep000*(ch[i*ND*ND+j*ND+k]-c0);
				sig12[i*ND*ND+j*ND+k]=2.*c44*ec12[i*ND*ND+j*ND+k];
				sig13[i*ND*ND+j*ND+k]=2.*c44*ec13[i*ND*ND+j*ND+k];
				sig23[i*ND*ND+j*ND+k]=2.*c44*ec23[i*ND*ND+j*ND+k];

				//Calculation of elastic strain energy field [equation (5.28)]
				// Estr(r) = (1/2) * C ijkl * (epsilon_c kl - epsilon_0 kl) * (epsilon_c kl - epsilon_0 kl)
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


//***** Calculation of the displacement field [equation (5.31)] *************************************
//--- Calculation of u1 ---
	iii=1;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=-qrh1[i][j][k]*zuij(i, j, iii)-qih1[i][j][k]*zuij(i, j, iii);
				//xi[i][j][k]= qrh1[i][j][k]*zuij(i, j, iii)+qih1[i][j][k]*zuij(i, j, iii);
				xr[i*ND*ND+j*ND+k]=-qrh1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND)
								   -qih1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND);
				xi[i*ND*ND+j*ND+k]= qrh1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND)
								   +qih1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND);
				//
				in[i*ND*ND+j*ND+k][0] = xr[i*ND*ND+j*ND+k]; //For IFFT (real)
				in[i*ND*ND+j*ND+k][1] = xi[i*ND*ND+j*ND+k]; //For IFFT (imag.)
			}
		}
	}
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND+j] = out[i*ND+j][0]/(fftsize*fftsize); //For IFFT (real)
				//
				//u1[i][j][k]=xr[i][j][k];
				u1[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- Calculation of u2 ---
	iii=2;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=-qrh1[i][j][k]*zuij(i, j, k, iii)-qih1[i][j][k]*zuij(i, j, k, iii);
				//xi[i][j][k]= qrh1[i][j][k]*zuij(i, j, k, iii)+qih1[i][j][k]*zuij(i, j, k, iii);
				xr[i*ND*ND+j*ND+k]=-qrh1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND)
								   -qih1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND);
				xi[i*ND*ND+j*ND+k]= qrh1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND)
								   +qih1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND);
				//
				in[i*ND*ND+j*ND+k][0] = xr[i*ND*ND+j*ND+k]; //For IFFT (real)
				in[i*ND*ND+j*ND+k][1] = xi[i*ND*ND+j*ND+k]; //For IFFT (imag.)
			}
		}
	}
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND+j] = out[i*ND+j][0]/(fftsize*fftsize); //For IFFT (real)
				//
				//u2[i][j][k]=xr[i][j][k];
				u2[i*ND+j]=xr[i*ND+j];
			}
		}
	}

//--- Calculation of u3 ---
	iii=3;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//xr[i][j][k]=-qrh1[i][j][k]*zuij(i, j, k, iii)-qih1[i][j][k]*zuij(i, j, k, iii);
				//xi[i][j][k]= qrh1[i][j][k]*zuij(i, j, k, iii)+qih1[i][j][k]*zuij(i, j, k, iii);
				xr[i*ND*ND+j*ND+k]=-qrh1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND)
								   -qih1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND);
				xi[i*ND*ND+j*ND+k]= qrh1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND)
								   +qih1[i*ND*ND+j*ND+k]*zuij(i, j, k, iii, ND);
				//
				in[i*ND*ND+j*ND+k][0] = xr[i*ND*ND+j*ND+k]; //For IFFT (real)
				in[i*ND*ND+j*ND+k][1] = xi[i*ND*ND+j*ND+k]; //For IFFT (imag.)
			}
		}
	}
	fftw_execute(iplan); //For IFFT
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND+j] = out[i*ND+j][0]/(fftsize*fftsize); //For IFFT (real)
				//
				//u3[i][j][k]=xr[i][j][k];
				u3[i*ND+j]=xr[i*ND+j];;
			}
		}
	}


//***** Save all numerical data ********************************************
	//datsave(ch, Estr, ep11c, ep22c, ep33c, ep12c, ep13c, ep23c, 
	//		 sig11, sig22, sig33, sig12, sig13, sig23, u1, u2, u3, ND);			//save data
	datsave_paraview(ch, Estr, ep11c, ep22c, ep33c, ep12c, ep13c, ep23c, 
			 sig11, sig22, sig33, sig12, sig13, sig23, u1, u2, u3, ND);	//save
	//if(keypress()){return 0;}//キー待ち状態
	printf("Finished \n");
	if(plan) fftw_destroy_plan(plan);		//For FFT
	if(iplan) fftw_destroy_plan(iplan);	//For IFFT
	fftw_free(in); fftw_free(out);		// For FFT and IFFT
	std::exit(0);
	//return 0;
}

//*** Zcij [computation of formula (5.26)] ****************************************
double zcij(int i0, int j0, int k0, int iii, int jjj, int ND)
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
	int nd=ND, nd2=ND/2;

	if(i0<=nd2-1) {ii=i0;}  if(i0>=nd2) {ii=i0-nd;}
	if(j0<=nd2-1) {jj=j0;}  if(j0>=nd2) {jj=j0-nd;}
	if(k0<=nd2-1) {kk=k0;}  if(k0>=nd2) {kk=k0-nd;}
	//alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj); // 2D
	alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj+(double)kk*(double)kk);
	if(alnn==0.0){alnn=1.0;}
	nec[1]=nx=(double)ii/alnn;	// n is the unit vector in the k direction, n=k/|k|
	nec[2]=ny=(double)jj/alnn;	// n is the unit vector in the k direction, n=k/|k|
	nec[3]=nz=(double)kk/alnn;	// n is the unit vector in the k direction, n=k/|k|
	//nec[3]=nz=0.0; // 2D

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
	if(det1==0.0){det1=1.0;}

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

//*** Zuij [computation of formula (5.30)] ****************************************
double zuij(int i0, int j0, int k0, int iii, int ND)
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
	int nd=ND, nd2=ND/2;

	if(i0<=nd2-1) {ii=i0;}  if(i0>=nd2) {ii=i0-nd;}
	if(j0<=nd2-1) {jj=j0;}  if(j0>=nd2) {jj=j0-nd;}
	if(k0<=nd2-1) {kk=k0;}  if(k0>=nd2) {kk=k0-nd;}
	//alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj); // 3D
	alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj+(double)kk*(double)kk);
	if(alnn==0.0){alnn=1.0;}
	nec[1]=nx=(double)ii/alnn;	// n is the unit vector in the k direction, n=k/|k|
	nec[2]=ny=(double)jj/alnn;	// n is the unit vector in the k direction, n=k/|k|
	nec[3]=nz=(double)kk/alnn;	// n is the unit vector in the k direction, n=k/|k|
	//nec[3]=nz=0.0; // 2D

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
	if(det1==0.0){det1=1.0;}

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
			//zij=zij-(om[m][iii]/alnn)*sigma[m][n]*nec[n]; // eq.(5.30) or eq.(II 3.9), alnn=|k}
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
	int rndxm, rndym, rndzm;

	datin0 = fopen("data.dat", "r");//Open the source file
	
	fscanf(datin0, "%d %d %d", &rndxm, &rndym, &rndzm);
	if (ndm != rndxm || ndm != rndym || ndm != rndzm){
		printf("data size is mismatch \n");
		printf("Please, change ND in parameters.txt \n");
	}

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
