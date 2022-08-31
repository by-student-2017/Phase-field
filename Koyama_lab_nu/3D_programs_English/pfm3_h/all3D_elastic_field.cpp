#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

//#include "wingxa.h"

#define DRND(x) ((double)(x)/RAND_MAX*rand())
#define ND 64
#define IG 6	//2^IG=ND
//#define INX 400					//Pixel size (x direction)
//#define INY 400					//Pixel size (y direction)

	int nd=ND; 					
	int ndm=ND-1;				//define ND-1
	int nd2=ND/2;				//define ND/2 for FFT
	int ig=IG;					//2^ig=ND
	double PI=3.14159;			//pi
	double rr=8.3145;			//Gas constant
	double time1;				//count
	int iout=0;

	double c0;					//average of phase field
	double ch[ND][ND][ND];		//phase field
	double eta_c[4][4];			//Lattice mismatch
	double cec[4][4][4][4];		//Elastic constant
	double sigma[4][4];			//Eigen stress (stress when elastically deformed by the amount of Eigen strain)
	double ep11c[ND][ND][ND], ep22c[ND][ND][ND], ep33c[ND][ND][ND];	//Strain
	double ep12c[ND][ND][ND], ep13c[ND][ND][ND], ep23c[ND][ND][ND];
	double sig11[ND][ND][ND], sig22[ND][ND][ND], sig33[ND][ND][ND];	//Elastic stress
	double sig12[ND][ND][ND], sig13[ND][ND][ND], sig23[ND][ND][ND];
	double Estr[ND][ND][ND];	//Elastic strain energy density
	double u1[ND][ND][ND], u2[ND][ND][ND], u3[ND][ND][ND];	//Displacement
	//double fld1[ND][ND][ND];

	//bag fix: move them from main program to here
	double qrh1[ND][ND][ND], qih1[ND][ND][ND];	//Fourier transform of the field
	double ec11[ND][ND][ND], ec22[ND][ND][ND], ec33[ND][ND][ND];	//Constrained strain array
	double ec12[ND][ND][ND], ec13[ND][ND][ND], ec23[ND][ND][ND];

	int qs;						//FFT(qs:-1) or IFFT(qs:1)
	double xi[ND][ND][ND], xr[ND][ND][ND], xif[ND], xrf[ND];//Real or imaginary array of Fourier transform
	double s[ND],c[ND];			//sin and cos table
	int ik[ND];				//Bit inversion table for FFT

	void table();				//Creating sin and cos tables and bit inversion tables Subrutin
	void fft();					//FFT1D
	void xyzfft();				//FFT3D
	double zcij(int i0, int j0, int k0, int iii, int jjj);//Coefficient calculation of elastic function (Fourier space)
	double zuij(int i0, int j0, int k0, int iii);//Displacement field coefficient calculation (Fourier space)

	//void ini000();				//Initial concentration wave setting subroutine
	void datin();				//Subroutine for initial field reading
	void datsave();				//Data storage subroutine
	void datsave_paraview();	//Data storage subroutine
	//void graph();				//Scalar field drawing subroutine
	//void graph_u();				//Vector field drawing subroutine

//******* main program ****************************************************
int main(void)
{
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
	//printf("delt(0.005)=  "); scanf(" %lf",&delt);	//Time increments
	//printf("c0(0.4)=  "); scanf(" %lf",&c0);	//Average of phase field
	//printf("temp(900.0)=  "); scanf(" %lf",&temp);	//Temperature [K]
	//delt=0.002;					//time increments
	//c0 = 0.4;					//Average of phase field
	//temp=900.0;					//Temperature [K]
	al=100.0*1.0E-09;			//calculation area [m]
	b1=al/nd;					//Difference block length

	time1=0.0;					//Setting the initial calculation count number
	time1max=1.0+1.0e+07;		//Setting the maximum calculation count

	a1=2.8664E-10; 				//Lattice constant (experimental data) of Fe(bcc)
	atom_n=2.0;	vm0=6.02E23*a1*a1*a1/atom_n;	//Calculation of molar volume (for bcc)

//*** Lattice mismatch settings ***
	ep000=0.05;
	eta_c[1][1]=eta_c[2][2]=eta_c[3][3]=ep000; 	//Lattice mismatch
	eta_c[1][2]=eta_c[2][1]=eta_c[1][3]=eta_c[3][1]=eta_c[2][3]=eta_c[3][2]=0.0;

//***** Elastic constant of Fe (bcc) ****************************
	el_mag=1.0E+11/1.0E+09;
	c11=2.33*el_mag;
	c12=1.35*el_mag;
	c44=1.18*el_mag;

	for(i=0;i<=3;i++){
		for(j=0;j<=3;j++){
			for(k=0;k<=3;k++){
				for(l=0;l<=3;l++){
					cec[i][j][k][l]=0.;//Initialization of elastic constant array
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

//*** sin and cos table, bit inversion table, and initial field settings ***************

	table();
	//ini000();
	datin(); //Read phase field

 	//gwinsize(INX,INY); ginit(1); gsetorg(0,0);// drow figure

//**** Start of calculation of elastic field analysis ******************************
//start:;

//**** Fourier transform of the field (phase field) ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i][j][k]=ch[i][j][k]; xi[i][j][k]=0.;
			}
		}
	}
	qs=-1; xyzfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				qrh1[i][j][k]=xr[i][j][k]; qih1[i][j][k]=xi[i][j][k];
			}
		}
	}
	qrh1[0][0][0]=qih1[0][0][0]=0.0;

//***** Calculation of total strain fluctuation *************************************
//--- ec11 ---
	iii=1; jjj=1;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
			}
		}
	}
	qs=1; xyzfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ec11[i][j][k]=xr[i][j][k];
			}
		}
	}

//--- ec22 ---
	iii=2; jjj=2;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
			}
		}
	}
	qs=1; xyzfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ec22[i][j][k]=xr[i][j][k];
			}
		}
	}

//--- ec33 ---
	iii=3; jjj=3;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
			}
		}
	}
	qs=1; xyzfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ec33[i][j][k]=xr[i][j][k];
			}
		}
	}

//--- ec12 ---
	iii=1; jjj=2;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
			}
		}
	}
	qs=1; xyzfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ec12[i][j][k]=xr[i][j][k];
			}
		}
	}

//--- ec13 ---
	iii=1; jjj=3;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
			}
		}
	}
	qs=1; xyzfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ec13[i][j][k]=xr[i][j][k];
			}
		}
	}

//--- ec23 ---
	iii=2; jjj=3;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
				xi[i][j][k]=qrh1[i][j][k]*zcij(i, j, k, iii, jjj)+qih1[i][j][k]*zcij(i, j, k, iii, jjj);
			}
		}
	}
	qs=1; xyzfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ec23[i][j][k]=xr[i][j][k];
			}
		}
	}

//***** Calculation of total strain field, stress field, and elastic strain energy field *************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//Calculation of total strain field
				ep11c[i][j][k]=ec11[i][j][k]+ep000*c0;
				ep22c[i][j][k]=ec22[i][j][k]+ep000*c0;
				ep33c[i][j][k]=ec33[i][j][k]+ep000*c0;
				ep12c[i][j][k]=ec12[i][j][k];
				ep13c[i][j][k]=ec13[i][j][k];
				ep23c[i][j][k]=ec23[i][j][k];

				//Calculation of elastic stress field
				sig11[i][j][k]=c11*ec11[i][j][k]+c12*ec22[i][j][k]+c12*ec33[i][j][k]
											   -(c11+2.*c12)*ep000*(ch[i][j][k]-c0);
				sig22[i][j][k]=c12*ec11[i][j][k]+c11*ec22[i][j][k]+c12*ec33[i][j][k]
											   -(c11+2.*c12)*ep000*(ch[i][j][k]-c0);
				sig33[i][j][k]=c12*ec11[i][j][k]+c12*ec22[i][j][k]+c11*ec33[i][j][k]
											   -(c11+2.*c12)*ep000*(ch[i][j][k]-c0);
				sig12[i][j][k]=2.*c44*ec12[i][j][k];
				sig13[i][j][k]=2.*c44*ec13[i][j][k];
				sig23[i][j][k]=2.*c44*ec23[i][j][k];

				//Calculation of elastic strain energy field
				Estr1=1.5*(c11+2.0*c12)*ep000*ep000*ch[i][j][k]*ch[i][j][k];
				Estr2=-(c11+2.0*c12)*(ep11c[i][j][k]
					                 +ep22c[i][j][k]
					                 +ep33c[i][j][k])*ep000*ch[i][j][k];
				Estr3=0.5*c11*(ep11c[i][j][k]*ep11c[i][j][k]
							  +ep22c[i][j][k]*ep22c[i][j][k]
							  +ep33c[i][j][k]*ep33c[i][j][k]);
				Estr4=c12*(ep11c[i][j][k]*ep22c[i][j][k]
						  +ep11c[i][j][k]*ep33c[i][j][k]
						  +ep22c[i][j][k]*ep33c[i][j][k]);
				Estr5=2.*c44*(ep12c[i][j][k]*ep12c[i][j][k]
							 +ep13c[i][j][k]*ep13c[i][j][k]
							 +ep23c[i][j][k]*ep23c[i][j][k]);
				Estr[i][j][k]=Estr1+Estr2+Estr3+Estr4+Estr5;
			}
		}
	}


//***** Displacement field calculation *************************************
//--- u1 ---
	iii=1;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i][j][k]=-qrh1[i][j][k]*zuij(i, j, k, iii)-qih1[i][j][k]*zuij(i, j, k, iii);
				xi[i][j][k]= qrh1[i][j][k]*zuij(i, j, k, iii)+qih1[i][j][k]*zuij(i, j, k, iii);
			}
	 	}
	}
	qs=1; xyzfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				u1[i][j][k]=xr[i][j][k];
			}
		}
	}

//--- u2 ---
	iii=2;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i][j][k]=-qrh1[i][j][k]*zuij(i, j, k, iii)-qih1[i][j][k]*zuij(i, j, k, iii);
				xi[i][j][k]= qrh1[i][j][k]*zuij(i, j, k, iii)+qih1[i][j][k]*zuij(i, j, k, iii);
			}
	 	}
	}
	qs=1; xyzfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				u2[i][j][k]=xr[i][j][k];
			}
		}
	}

//--- u3 ---
	iii=3;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i][j][k]=-qrh1[i][j][k]*zuij(i, j, k, iii)-qih1[i][j][k]*zuij(i, j, k, iii);
				xi[i][j][k]= qrh1[i][j][k]*zuij(i, j, k, iii)+qih1[i][j][k]*zuij(i, j, k, iii);
			}
	 	}
	}
	qs=1; xyzfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				u3[i][j][k]=xr[i][j][k];
			}
		}
	}


//*****Display of the place, save of images, save of data  ********************************************
//--- phase field ----------
//	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ 
//		fld1[i][j][k]=ch[i][j][k]; 
//	} } }
//  graph();  save_screen("ch_field.bmp"); //show graph and save

//--- Elastic strain energy field (standardized) ----------
//	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ 
//		fld1[i][j][k]=(Estr[i][j][k]-0.0)/(0.4-0.0); 
//	} } }
//  graph();  save_screen("Estr_field.bmp");//show graph and save

//--- Total strain field (standardized) ----------
//	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ 
//		fld1[i][j][k]=(ep11c[i][j][k]-(-0.1))/(0.1-(-0.1)); 
//	} } }
//  graph();  save_screen("ep11c_field.bmp");//show graph and save

//	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ 
//		fld1[i][j][k]=(ep22c[i][j][k]-(-0.1))/(0.1-(-0.1)); 
//	} } }
//  graph();  save_screen("ep22c_field.bmp");//show graph and save

//	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ 
//		fld1[i][j][k]=(ep33c[i][j][k]-(-0.1))/(0.1-(-0.1)); 
//	} } }
//  graph();  save_screen("ep33c_field.bmp");//show graph and save

//	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ 
//		fld1[i][j][k]=(ep12c[i][j][k]-(-0.1))/(0.1-(-0.1)); 
//	} } }
//  graph();  save_screen("ep12c_field.bmp");//show graph and save

//	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ 
//		fld1[i][j][k]=(ep13c[i][j][k]-(-0.1))/(0.1-(-0.1)); 
//	} } }
//  graph();  save_screen("ep13c_field.bmp");//show graph and save

//	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ 
//		fld1[i][j][k]=(ep23c[i][j][k]-(-0.1))/(0.1-(-0.1)); 
//	} } }
//  graph();  save_screen("ep23c_field.bmp");//show graph and save

//--- Stress field (standardized) ----------
//	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ 
//		fld1[i][j][k]=(sig11[i][j][k]-(-10.0))/(10.0-(-10.0)); 
//	} } }
//  graph();  save_screen("sig11_field.bmp");//show graph and save

//	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ 
//		fld1[i][j][k]=(sig22[i][j][k]-(-10.0))/(10.0-(-10.0)); 
//	} } }
//  graph();  save_screen("sig22_field.bmp");//show graph and save

//	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ 
//		fld1[i][j][k]=(sig33[i][j][k]-(-10.0))/(10.0-(-10.0)); 
//	} } }
//  graph();  save_screen("sig33_field.bmp");//show graph and save

//	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ 
//		fld1[i][j][k]=(sig12[i][j][k]-(-10.0))/(10.0-(-10.0)); 
//	} } }
//  graph();  save_screen("sig12_field.bmp");//show graph and save

//	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ 
//		fld1[i][j][k]=(sig13[i][j][k]-(-10.0))/(10.0-(-10.0)); 
//	} } }
//  graph();  save_screen("sig13_field.bmp");//show graph and save

//	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ 
//		fld1[i][j][k]=(sig23[i][j][k]-(-10.0))/(10.0-(-10.0)); 
//	} } }
//  graph();  save_screen("sig23_field.bmp");//show graph and save

//--- Displacement field (display of vector field) ----------
	//graph_u();//Since the displacement field is a vector field, another subroutine graph_u is used for drawing.
	//save_screen("u_field.bmp");//save graph


//***** Save all numerical data ********************************************
	//datsave();			//save data
	datsave_paraview();	//save

	//if(keypress()){return 0;}//wait key input
	return 0;

}


//******* Sin, Cos table and bit inversion table settings ***************
void table()
{
	int i, j, mc, mn;
	double q;

	q=2.0*PI/nd;
	for(i=0;i<=nd2-1;i++){ c[i]=cos(q*i); s[i]=sin(q*i); }//Sin, Cos

	ik[0]=0; mn=nd2; mc=1;
	for(i=1;i<=ig;i++){
		for(j=0;j<=mc-1;j++){ ik[j+mc]=ik[j]+mn; }	//Bit inversion table
		mc*=2; mn/=2; 
	}
}

//********** FFT 1D **************************************
void fft()
{
	int lg, lf, mf, nf, n2, ka, kb, ix;
	double tj, tr;

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
void xyzfft()
{
	int i, j, k, iz, iy, ix;

//*** xy-z ***
	for(ix=0;ix<=ndm;ix++){
		for(iy=0;iy<=ndm;iy++){
			for(iz=0;iz<=ndm;iz++){
				xrf[iz]=xr[ix][iy][iz];
				xif[iz]=xi[ix][iy][iz];
			}
			fft();
			for(iz=0;iz<=ndm;iz++){
				xr[ix][iy][iz]=xrf[ik[iz]];
				xi[ix][iy][iz]=xif[ik[iz]];
			}
		}
	}

//*** xz-y ***
	for(ix=0;ix<=ndm;ix++){
		for(iz=0;iz<=ndm;iz++){
			for(iy=0;iy<=ndm;iy++){
				xrf[iy]=xr[ix][iy][iz];
				xif[iy]=xi[ix][iy][iz];
			}
			fft();
			for(iy=0;iy<=ndm;iy++){
				xr[ix][iy][iz]=xrf[ik[iy]];
				xi[ix][iy][iz]=xif[ik[iy]];
			}
		}
	}

//*** yz-x ***
	for(iz=0;iz<=ndm;iz++){
		for(iy=0;iy<=ndm;iy++){
			for(ix=0;ix<=ndm;ix++){
				xrf[ix]=xr[ix][iy][iz];
				xif[ix]=xi[ix][iy][iz];
			}
			fft();
			for(ix=0;ix<=ndm;ix++){
				xr[ix][iy][iz]=xrf[ik[ix]];
				xi[ix][iy][iz]=xif[ik[ix]];
			}
		}
	}

	if(qs > 0){return;}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
			xr[i][j][k]=xr[i][j][k]/nd/nd/nd;
			xi[i][j][k]=xi[i][j][k]/nd/nd/nd;
			}
		}
	}

}

//*** Zcij [eq.(5.26) or eq.(II 3.5)] ****************************************
double zcij(int i0, int j0, int k0, int iii, int jjj)
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

	if(i0<=nd2-1){ii=i0;}  if(i0>=nd2){ii=i0-nd;}
	if(j0<=nd2-1){jj=j0;}  if(j0>=nd2){jj=j0-nd;}
	if(k0<=nd2-1){kk=k0;}  if(k0>=nd2){kk=k0-nd;}
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
double zuij(int i0, int j0, int k0, int iii)
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

	if(i0<=nd2-1){ii=i0;}  if(i0>=nd2){ii=i0-nd;}
	if(j0<=nd2-1){jj=j0;}  if(j0>=nd2){jj=j0-nd;}
	if(k0<=nd2-1){kk=k0;}  if(k0>=nd2){kk=k0-nd;}
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
			//zij=zij-(om[m][iii]/alnn)*sigma[m][n]*nec[n]; // eq.(5.30) or eq.(II 3.9), alnn=|k}
    		zij=zij+sigma[m][n]*nec[n]*om[m][iii]; // eq.(5.30) or eq.(II 3.9)
		}
	}
	zij=-zij/alnn; //alnn=|k}
	return(zij);
}

//******* Scalar field display ***************************************
//void graph()
//{
//	int i, j, k, ii, jj, kk;					//integer
//	double col, col_R, col_G, col_B, col_RG;	//color
//	int ixmin=0, iymin=0, igx, igy, irad0;		//Screen coordinate system settings
//	double x, xmax, xmin, y, ymax, ymin, rad0, dia0;//Standardized coordinate system settings
//	int ixmax=INX, iymax=INY;					//show window

//  gcls();//clear
//	xmin=0.; xmax=1.; ymin=0.; ymax=1.;			//Drawing area (standardized)
//	printf("time %f\n", time1);					//show count
//	dia0=1./nd; 
//	rad0=dia0/2.;  	            irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;
	//Half the length of the diff block	//Convert to screen coordinate system (+1 is truncation correction when converting to integer)

//	kk=nd2;
//	for(i=0;i<=nd;i++){
//		for(j=0;j<=nd;j++){
//			x=rad0+dia0*i;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
//			y=rad0+dia0*j;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
//			//Coordinate calculation			//Convert to screen coordinate system
//			ii=i; jj=j;  if(i==nd){ii=0;}  if(j==nd){jj=0;}//Periodic boundary conditions
//			col_G=col_B=col_R=fld1[ii][jj][kk];//Set the field color in RGB
//			if(col_R>=1.){col_R=1.;}  if(col_R<=0.){col_R=0.;}//RGB domain correction
//			if(col_G>=1.){col_G=1.;}  if(col_G<=0.){col_G=0.;}
//			if(col_B>=1.){col_B=1.;}  if(col_B<=0.){col_B=0.;}
//			gcolor((int)(255*col_R),(int)(255*col_G),(int)(255*col_B));//color settings
//			grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);//Drawing a middle-painted rectangle
//		}
//	}
//	//swapbuffers();//Screen swap
//} 

//******* Displacement field display ***************************************
//void graph_u()
//{
//	int i, j, k, ii, jj, kk;//integer
//	double col;//color
//	double x, xmax, xmin, y, ymax, ymin, rad0, dia0;//Standardized coordinate system setting
//	int ixmin=0, iymin=0, igx, igy, irad0;//Screen coordinate system setting
//	double x1, y1, x2, y2, x3, y3, x4, y4, fact1, fact2, th0;//Arrow-related coordinate settings
//	int igx1, igy1, igx2, igy2, igx3, igy3, igx4, igy4;//Arrow-related screen coordinate settings
//	int ixmax=INX, iymax=INY;//Drawing Window range

//  gcls(); //Clear screen
//	xmin=0.; xmax=1.; ymin=0.; ymax=1.;//Drawing area (standardized)

//	printf("time %f\n",time1);//Display of calculated counts
//	dia0=1./nd;  
//	rad0=dia0/2.;               irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;
	//Half the length of the diff block	//Convert to screen coordinate system (+1 is truncation correction when converting to integer)

//**** Draw with a thin phase field ***************************
//	kk=nd2;
//	for(i=0;i<=nd;i++){
//		for(j=0;j<=nd;j++){
//			x=rad0+dia0*i;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
//			y=rad0+dia0*j;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
			//Coordinate calculation			//Convert to screen coordinate system
//			ii=i; jj=j;  if(i==nd){ii=0;}  if(j==nd){jj=0;}//Periodic boundary conditions
//			col=1.0-(1.0-ch[ii][jj][kk])*0.5;//Set the field color in RGB (multiply 0.5 to make it lighter)
//			if(col>=1.){col=1.;}  if(col<=0.){col=0.;}
//			gcolor((int)(255*col),(int)(255*col),(int)(255*col));//color settings
//			grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);//Drawing a middle-painted rectangle
//		}
//	}

//**** Draw the displacement field with an arrow ***************************
//	fact1=5.;//Arrow length
//	fact2=0.6;//Arrow tip position setting
//	th0=0.5*3.14159/2.;//Angle of the arrow part

//	kk=nd2;
//	for(i=0;i<=nd;i+=4){
//		for(j=0;j<=nd;j+=4){
//			x=rad0+dia0*i;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
//			y=rad0+dia0*j;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
//			//Coordinate calculation			//Convert to screen coordinate system
//			ii=i; jj=j; if(i==nd){ii=0;}  if(j==nd){jj=0;}//Periodic boundary conditions

//--- Coordinate calculation of each part of the arrow ---
//			x1=x-0.5*u1[ii][jj][kk]*fact1; y1=y-0.5*u2[ii][jj][kk]*fact1;
//			igx1=(ixmax-ixmin)/(xmax-xmin)*(x1-xmin)+ixmin;
//			igy1=(iymax-iymin)/(ymax-ymin)*(y1-ymin)+iymin;

//			x2=x+0.5*u1[ii][jj][kk]*fact1; y2=y+0.5*u2[ii][jj][kk]*fact1;
//			igx2=(ixmax-ixmin)/(xmax-xmin)*(x2-xmin)+ixmin;
//			igy2=(iymax-iymin)/(ymax-ymin)*(y2-ymin)+iymin;

//			x3=x+0.5*u1[ii][jj][kk]*fact1*fact2*cos(th0)-0.5*u2[ii][jj][kk]*fact1*fact2*sin(th0); 
//			y3=y+0.5*u1[ii][jj][kk]*fact1*fact2*sin(th0)+0.5*u2[ii][jj][kk]*fact1*fact2*cos(th0);
//			igx3=(ixmax-ixmin)/(xmax-xmin)*(x3-xmin)+ixmin;
//			igy3=(iymax-iymin)/(ymax-ymin)*(y3-ymin)+iymin;

//			x4=x+0.5*u1[ii][jj][kk]*fact1*fact2*cos(-th0)-0.5*u2[ii][jj][kk]*fact1*fact2*sin(-th0); 
//			y4=y+0.5*u1[ii][jj][kk]*fact1*fact2*sin(-th0)+0.5*u2[ii][jj][kk]*fact1*fact2*cos(-th0);
//			igx4=(ixmax-ixmin)/(xmax-xmin)*(x4-xmin)+ixmin;
//			igy4=(iymax-iymin)/(ymax-ymin)*(y4-ymin)+iymin;
//----------------------

			//circle(igx2, igy2, 1, 0);  paint(igx2, igy2,0,0);
//			gcolor(0,0,0);//color setting
//			gline(igx1, igy1, igx2, igy2);//Draw a line (arrow part)
//			gline(igx2, igy2, igx3, igy3);//Draw a line (arrow part)
//			gline(igx2, igy2, igx4, igy4);//Draw a line (arrow part)
//		}
//	}
//	//swapbuffers();//Screen swap
//}

//************ Saving data of various places *******************************
void datsave()
{
	FILE		*stream;
	int 		i, j, k;

	stream = fopen("el_field.dat", "w");

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				fprintf(stream, "%e  ", ch[i][j][k]);
				fprintf(stream, "%e  ", Estr[i][j][k]);
				fprintf(stream, "%e  ", ep11c[i][j][k]);
				fprintf(stream, "%e  ", ep22c[i][j][k]);
				fprintf(stream, "%e  ", ep33c[i][j][k]);
				fprintf(stream, "%e  ", ep12c[i][j][k]);
				fprintf(stream, "%e  ", ep13c[i][j][k]);
				fprintf(stream, "%e  ", ep23c[i][j][k]);
				fprintf(stream, "%e  ", sig11[i][j][k]);
				fprintf(stream, "%e  ", sig22[i][j][k]);
				fprintf(stream, "%e  ", sig33[i][j][k]);
				fprintf(stream, "%e  ", sig12[i][j][k]);
				fprintf(stream, "%e  ", sig13[i][j][k]);
				fprintf(stream, "%e  ", sig23[i][j][k]);
				fprintf(stream, "%e  ", u1[i][j][k]);
				fprintf(stream, "%e  ", u2[i][j][k]);
				fprintf(stream, "%e  ", u3[i][j][k]);
			}
		}
	}
	fprintf(stream, "\n");
	fclose(stream);
}

void datsave_paraview()
{
	FILE	*fp;
	char	fName[256];
	int 	i, j, k;
	
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
				fprintf(fp,"%10.6f\n", ch[i][j][k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_strain_energy_density float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", Estr[i][j][k]);
			}
		}
	}
	fprintf(fp,"SCALARS Strain_c11 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", ep11c[i][j][k]);
			}
		}
	}
	fprintf(fp,"SCALARS Strain_c22 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", ep22c[i][j][k]);
			}
		}
	}
	fprintf(fp,"SCALARS Strain_c33 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", ep33c[i][j][k]);
			}
		}
	}
	fprintf(fp,"SCALARS Strain_c12 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", ep12c[i][j][k]);
			}
		}
	}
	fprintf(fp,"SCALARS Strain_c13 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", ep13c[i][j][k]);
			}
		}
	}
	fprintf(fp,"SCALARS Strain_c23 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", ep23c[i][j][k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_11 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", sig11[i][j][k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_22 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", sig22[i][j][k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_33 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", sig33[i][j][k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_12 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", sig12[i][j][k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_13 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", sig13[i][j][k]);
			}
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_23 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", sig23[i][j][k]);
			}
		}
	}
	fprintf(fp,"SCALARS Displacement_u1 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", u1[i][j][k]);
			}
		}
	}
	fprintf(fp,"SCALARS Displacement_u2 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", u2[i][j][k]);
			}
		}
	}
	fprintf(fp,"SCALARS Displacement_u3 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", u3[i][j][k]);
			}
		}
	}
	fclose(fp);
}
//************ Reading field data *****************************************
void datin()
{
	FILE *datin0;//Stream pointer setting
	int  i, j, k;//integer
	double c00;//Average value of the field

	datin0 = fopen("data.dat", "r");//Open the source file

	fscanf(datin0, "%lf", &time1);

	c00=0.0;//Initial value of field mean
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				fscanf(datin0, "%lf", &ch[i][j][k]);//Read field data
				c00+=ch[i][j][k];
			}
		}
	}
	c0=c00/nd/nd/nd;//Average value of the field
	printf("c0=  %f  \n", c0);

	fclose(datin0);
}
//************[Initial concentration wave]*****************************************
//void ini000()
//{
//	int i, j, k;

//	for(i=0;i<=ndm;i++){
//		for(j=0;j<=ndm;j++){
//			for(k=0;k<=ndm;k++){
//				//ch[i][j][k]=0.0;
//				ch[i][j][k]=c0+0.1*(2.*DRND(1)-1.);//Spinodal decomposition
//			}
//		}
//	}
//
//}