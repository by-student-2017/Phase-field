#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#include <fftw3.h>
//g++ martensite_v2_libfftw3.cpp -lfftw3
//#include <mpi.h> //mpi version
//#include <fftw3-mpi.h> //mpi version

#define DRND(x) ((double)(x)/RAND_MAX*rand())//Random number function setting

	double PI=3.14159;				//Pi
	double RR=8.3145;				//gas constant
	double time1;					//Number of calculation counts (proportional to time)
	double ep0c=8.8541878e-12; 		//Vacuum permittivity (F/m)
	double Peq;						//equilibrium value of the moment of polarization
	int iout;
	
	double cec[4][4][4][4];		//Elastic constant
	double sigma[4][4];			//Eigen stress (stress when elastically deformed by the amount of Eigen strain)
	
	double zcij(int i0, int j0, int k0, int iii, int jjj, int ND);//Coefficient calculation of elastic function (Fourier space)
	
	void ini000(double *s1h, double *s2h, double *s3h, int ND);	//initial profile of the polarization moment at time 0	void datsave(double *s1h, double *s2h, int ND);	//data save subroutine
	void datsave(double *s1h, double *s2h, double *s3h, int ND);	//data save subroutine
	void datsave_paraview(double *s1h, double *s2h, double *s3h, int ND);	//data save subroutine
	void datin(double *s1h, double *s2h, double *s3h, int ND);	//Initial data loading

//******* main program ******************************************
int main(void)
{
	int ND;
	int nd, ndm, nd2, ig;
	
	int   i, j, k, l;
	int   ii, jj, kk;
	int   iii, jjj;
	int   Nstep;	//integer
	int   ip, im, jp, jm, kp, km;			//integer

	double ss1ss2, ss1ss3, ss2ss3;			//Work variables for correction of numerical errors in domain

	double a0_aa, a0_a, a0_c;				//Lattice constant of BaTiO3 (tetragonal)
	double al, temp, delt, vm0;				//Side length of calculation domain, temperature, time step, molar volume
	double time1max;						//Maximum calculation count (calculation end count)
	double b1;								//Length of one side of difference block

	double s1, s2, s3;						//x, y and z components of the moment of polarization
	double s1ip, s1im, s1jp, s1jm, s1kp, s1km;	//left, right, top, bottom, up, down, values of s1
	double s2ip, s2im, s2jp, s2jm, s2kp, s2km;	//left, right, top, bottom, up, down, values of s2
	double s3ip, s3im, s3jp, s3jm, s3kp, s3km;	//left, right, top, bottom, up, down, values of s3

	double s1k, s1k_chem, s1k_surf, s1k_str, s1k_ddi;	//potential for s1
	double s2k, s2k_chem, s2k_surf, s2k_str, s2k_ddi;	//potential for s2
	double s3k, s3k_chem, s3k_surf, s3k_str, s3k_ddi;	//potential for s3
	double smob1, smob2, smob3;				//Domain interface mobility
	double s1ddtt, s2ddtt, s3ddtt;			//left-hand side of the evolution equation

	double A1, A11, A12;					//parameters in chemical free energy
	double A1e, A11e, A12e;
	double A111, A112, A123;
	double A111e, A112e, A123e;
	double A1111, A1112, A1122, A1123;
	double A1111e, A1112e, A1122e, A1123e;
	double Tc0;

	double nx, ny, nz, alnn;				//Unit vector components of reciprocal space and reciprocal lattice vector length
	double kapaP, kapaPc;					//gradient energy factor

	double E1_ex; 							//external electric field
	double E1_ex_0x, E1_ex_x; 				//external electric field (x)
	double E1_ex_0y, E1_ex_y; 				//external electric field (y)
	double E1_ex_0z, E1_ex_z; 				//external electric field (z)
	double ep00;							//dielectric constant
	double Add0, Add0c;						//Coefficients in dipole-dipole interaction calculations

	double t1, t2, t3, tQ, tR, tS, tT;		//Working variables for calculation of equilibrium moment of polarization
	double dtemp;
	int readff;
	
	double C11, C22, C33;					//elastic constant
	double C12, C21, C13, C31, C23, C32;
	double C44, C55, C66;
	double Q11, Q22, Q33;					//electrostrictive coefficients [m^4/C^2]
	double Q12, Q21, Q13, Q31, Q23, Q32;
	double Q44, Q55, Q66;
	double ep_a[4][4];							//elastic strain related external force
	//
	double sum11, sum22, sum33;					//spatial integral of s1 and s2
	double sum12, sum13, sum23;					//spatial integral of s1 and s2
	double ep0[4][4];							//Mean value of transformation strain in the structure
	
	double epT[4][4];
	double el_fac;								//Normalization constant of elastic constant

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
	dtemp   = data[5];
	al      = data[6];
	vm0     = data[7];
	smob1   = data[8];
	smob2   = data[9];
	smob3   = data[10];
	A1e     = data[11]; // Landau expansion form
	A11e    = data[12]; // Landau expansion form
	A12e    = data[13]; // Landau expansion form
	A111e   = data[14]; // Landau expansion form
	A112e   = data[15]; // Landau expansion form
	A123e   = data[16]; // Landau expansion form
	A1111e  = data[17]; // Landau expansion form
	A1112e  = data[18]; // Landau expansion form
	A1122e  = data[19]; // Landau expansion form
	A1123e  = data[20]; // Landau expansion form
	Tc0     = data[21];
	kapaPc  = data[22];
	E1_ex_x = data[23];
	E1_ex_y = data[24];
	E1_ex_z = data[25];
	ep00    = data[26];
	Add0c   = data[27];
	readff  = int(data[28]);
	el_fac=1.0E+11/1.0E+09;
	C11     = data[29]*el_fac; C22=C33=C11;
	C12     = data[30]*el_fac; C21=C12; C31=C13=C12;
	C44     = data[31]*el_fac; C55=C66=C44;
	ep_a[1][1]=data[32];
	ep_a[2][2]=data[33];
	ep_a[3][3]=data[34];
	Q11     = data[35]; Q22=Q33=Q11;
	Q12     = data[36]; Q21=Q12; Q31=Q13=Q12;
	Q44     = data[37]; Q55=Q66=Q44;
	printf("---------------------------------\n");
	//
	ig = int(log2(ND));
	nd=ND;
	ndm=ND-1;
	nd2=ND/2;
	//
	double *s1h     = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//polarization moment in x direction
	double *s2h     = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//polarization moment in y direction
	double *s3h     = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//polarization moment in z direction	
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
	//After calc.: fftw_destroy_plan(plan); fftw_free(in); fftw_free(out);
	//
	double *s1qrh   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s1 (real part) (x)
	double *s2qrh   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s2 (real part) (y)
	double *s3qrh   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s2 (real part) (z)
	//
	double *s1qih   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s1 (imaginary part) (x)
	double *s2qih   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s2 (imaginary part) (y)
	double *s3qih   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Fourier transform of s2 (imaginary part) (z)
	//
	double *s1h2    = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//auxiliary array of s1
	double *s2h2    = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//auxiliary array of s2
	double *s3h2    = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//auxiliary array of s3
	//
	double *s1k_dd  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//dipole-dipole interaction potential (x)
	double *s2k_dd  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//dipole-dipole interaction potential (y)
	double *s3k_dd  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//dipole-dipole interaction potential (z)
	//
	double *ep11qrh0= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep22qrh0= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep33qrh0= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep12qrh0= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep21qrh0= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep13qrh0= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep31qrh0= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep23qrh0= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep32qrh0= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	//
	double *ep11qih0= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep22qih0= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep33qih0= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep12qih0= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep21qih0= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep13qih0= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep31qih0= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep23qih0= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep32qih0= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	//
	double *ep11h0  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep22h0  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep33h0  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep12h0  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep21h0  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep13h0  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep31h0  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep23h0  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ep32h0  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	//
	double *sigma11_r = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma22_r = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma33_r = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma12_r = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma21_r = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma13_r = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma31_r = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma23_r = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma32_r = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	//
	double *sigma11_i = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma22_i = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma33_i = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma12_i = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma21_i = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma13_i = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma31_i = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma23_i = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma32_i = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	//
	double *ec11     = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ec22     = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ec33     = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ec12     = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ec13     = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ec23     = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	
	time1=0.0;								//Initial calculation count setting
	
	b1=al*1.0E-06/nd;						//Length of one side of difference block (m)

//--- Parameter value in chemical free energy [see Table 4.7]----------------------
	A1=A1e*vm0/RR/temp;
	A11=A11e*vm0/RR/temp;
	A12=A12e*vm0/RR/temp;
	A111=A111e*vm0/RR/temp;
	A112=A112e*vm0/RR/temp;
	A123=A123e*vm0/RR/temp;
	A1111=A1111e*vm0/RR/temp;
	A1112=A1112e*vm0/RR/temp;
	A1122=A1122e*vm0/RR/temp;
	A1123=A1123e*vm0/RR/temp;

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

	kapaP=kapaPc/ep0c*vm0/RR/temp/b1/b1;//gradient energy factor
	
	E1_ex_0x=E1_ex_x*vm0/RR/temp;//External electric field (x direction)
	E1_ex_0y=E1_ex_y*vm0/RR/temp;//External electric field (y direction)
	E1_ex_0z=E1_ex_z*vm0/RR/temp;//External electric field (z direction)

	Add0=Add0c/ep0c/ep00*vm0/RR/temp;//Coefficients in dipole-dipole interaction calculations [see formula (4.52)]

//--- elastic constant ----------------------
	for(i=0;i<=3;i++){
		for(j=0;j<=3;j++){
			for(k=0;k<=3;k++){
				for(l=0;l<=3;l++){
					cec[i][j][k][l]=0.0;//Initialization of elastic constant array
				}
			}
		}
	}
	
	cec[1][1][1][1]=C11;
	cec[2][2][2][2]=C22;
	cec[3][3][3][3]=C33;
	cec[2][3][2][3]=cec[2][3][3][2]=cec[3][2][2][3]=cec[3][2][3][2]=C44;
	cec[1][3][1][3]=cec[1][3][3][1]=cec[3][1][1][3]=cec[3][1][3][1]=C55;
	cec[1][2][1][2]=cec[1][2][2][1]=cec[2][1][1][2]=cec[2][1][2][1]=C66;
	cec[1][1][2][2]=cec[2][2][1][1]=C12;
	cec[1][1][3][3]=cec[3][3][1][1]=C13;
	cec[2][2][3][3]=cec[3][3][2][2]=C23;

//*** Setting up sin and cos tables, bit reversal tables, and initial fields ***************
	
	if(readff == 0){
		ini000(s1h, s2h, s3h, ND);	//initial profile of the polarization moment at time 0
	} else {
		datin(s1h, s2h, s3h, ND);	//Input initial tissue field
	}

//**** Simulation start ******************************
//Nstep = 10;
iout = -1;
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
				s1qrh[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
				s1qih[i*ND*ND+j*ND+k]=xi[i*ND*ND+j*ND+k];
			}
		}
	}
	s1qrh[0]=s1qih[0]=0.0;

//**** Fourier transform of s2 [equation (3.27)] ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
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
				s3qrh[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
				s3qih[i*ND*ND+j*ND+k]=xi[i*ND*ND+j*ND+k];
			}
		}
	}
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
				//
				xr[i*ND*ND+j*ND+k]=s1qrh[i*ND*ND+j*ND+k]*nx*nx+s2qrh[i*ND*ND+j*ND+k]*nx*ny+s3qrh[i*ND*ND+j*ND+k]*nx*nz;
				xi[i*ND*ND+j*ND+k]=s1qih[i*ND*ND+j*ND+k]*nx*nx+s2qih[i*ND*ND+j*ND+k]*nx*ny+s3qih[i*ND*ND+j*ND+k]*nx*nz;
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
				//
				xr[i*ND*ND+j*ND+k]=s1qrh[i*ND*ND+j*ND+k]*ny*nx+s2qrh[i*ND*ND+j*ND+k]*ny*ny+s3qrh[i*ND*ND+j*ND+k]*ny*nz;
				xi[i*ND*ND+j*ND+k]=s1qih[i*ND*ND+j*ND+k]*ny*nx+s2qih[i*ND*ND+j*ND+k]*ny*ny+s3qih[i*ND*ND+j*ND+k]*ny*nz;
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
				//
				xr[i*ND*ND+j*ND+k]=s1qrh[i*ND*ND+j*ND+k]*nz*nx+s2qrh[i*ND*ND+j*ND+k]*nz*ny+s3qrh[i*ND*ND+j*ND+k]*nz*nz;
				xi[i*ND*ND+j*ND+k]=s1qih[i*ND*ND+j*ND+k]*nz*nx+s2qih[i*ND*ND+j*ND+k]*nz*ny+s3qih[i*ND*ND+j*ND+k]*nz*nz;
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
				s3k_dd[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//**** Fourier transform of the Eigen strain field [equation (5)] ep11 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k]=ep11h0[i*ND*ND+j*ND+k]
								  =Q11*s1h[i*ND*ND+j*ND+k]*s1h[i*ND*ND+j*ND+k]
								  +Q12*s2h[i*ND*ND+j*ND+k]*s2h[i*ND*ND+j*ND+k]
								  +Q12*s3h[i*ND*ND+j*ND+k]*s3h[i*ND*ND+j*ND+k];
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
				ep11qrh0[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
				ep11qih0[i*ND*ND+j*ND+k]=xi[i*ND*ND+j*ND+k];
			}
		}
	}
	ep11qrh0[0]=ep11qih0[0]=0.0;

//**** Fourier transform of the Eigen strain field [equation (5)] ep22 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k]=ep22h0[i*ND*ND+j*ND+k]
								  =Q12*s1h[i*ND*ND+j*ND+k]*s1h[i*ND*ND+j*ND+k]
								  +Q11*s2h[i*ND*ND+j*ND+k]*s2h[i*ND*ND+j*ND+k]
								  +Q12*s3h[i*ND*ND+j*ND+k]*s3h[i*ND*ND+j*ND+k];
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
				ep22qrh0[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
				ep22qih0[i*ND*ND+j*ND+k]=xi[i*ND*ND+j*ND+k];
			}
		}
	}
	ep22qrh0[0]=ep22qih0[0]=0.0;
	
//**** Fourier transform of the Eigen strain field ep33 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k]=ep33h0[i*ND*ND+j*ND+k]
								  =Q12*s1h[i*ND*ND+j*ND+k]*s1h[i*ND*ND+j*ND+k]
								  +Q12*s2h[i*ND*ND+j*ND+k]*s2h[i*ND*ND+j*ND+k]
								  +Q11*s3h[i*ND*ND+j*ND+k]*s3h[i*ND*ND+j*ND+k];
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
				ep33qrh0[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
				ep33qih0[i*ND*ND+j*ND+k]=xi[i*ND*ND+j*ND+k];
			}
		}
	}
	ep33qrh0[0]=ep33qih0[0]=0.0;

//**** Fourier transform of the Eigen strain field ep12 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k]=ep12h0[i*ND*ND+j*ND+k]
								  =Q66*s1h[i*ND*ND+j*ND+k]*s2h[i*ND*ND+j*ND+k];
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
				ep12qrh0[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
				ep12qih0[i*ND*ND+j*ND+k]=xi[i*ND*ND+j*ND+k];
			}
		}
	}
	ep12qrh0[0]=ep12qih0[0]=0.0;

//**** Fourier transform of the Eigen strain field ep13 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k]=ep13h0[i*ND*ND+j*ND+k]
								  =Q55*s1h[i*ND*ND+j*ND+k]*s3h[i*ND*ND+j*ND+k];
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
				ep13qrh0[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
				ep13qih0[i*ND*ND+j*ND+k]=xi[i*ND*ND+j*ND+k];
			}
		}
	}
	ep13qrh0[0]=ep13qih0[0]=0.0;

//**** Fourier transform of the Eigen strain field ep23 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k]=ep23h0[i*ND*ND+j*ND+k]
								  =Q44*s2h[i*ND*ND+j*ND+k]*s3h[i*ND*ND+j*ND+k];
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
				ep23qrh0[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
				ep23qih0[i*ND*ND+j*ND+k]=xi[i*ND*ND+j*ND+k];
			}
		}
	}
	ep23qrh0[0]=ep23qih0[0]=0.0;

//*** Calculation of the average value of the Eigen strain field ***
	sum11=sum22=sum33=0.0;
	sum12=sum13=sum23=0.0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//
				sum11+=ep11h0[i*ND*ND+j*ND+k];
				sum22+=ep22h0[i*ND*ND+j*ND+k];
				sum33+=ep33h0[i*ND*ND+j*ND+k];
				//
				sum12+=ep12h0[i*ND*ND+j*ND+k];
				sum13+=ep13h0[i*ND*ND+j*ND+k];
				sum23+=ep23h0[i*ND*ND+j*ND+k];
			}
		}
	}
	ep0[1][1]=sum11/nd/nd/nd; ep0[2][2]=sum22/nd/nd/nd; ep0[3][3]=sum33/nd/nd/nd;
	ep0[1][2]=ep0[2][1]=sum12/nd/nd/nd;
	ep0[1][3]=ep0[3][1]=sum13/nd/nd/nd;
	ep0[2][3]=ep0[3][2]=sum23/nd/nd/nd;
	
//***** Calculation of sigma *************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//
				sigma11_r[i*ND*ND+j*ND+k]=cec[1][1][1][1]*ep11qrh0[i*ND*ND+j*ND+k]
										 +cec[1][1][2][2]*ep22qrh0[i*ND*ND+j*ND+k]
										 +cec[1][1][3][3]*ep33qrh0[i*ND*ND+j*ND+k];
				sigma22_r[i*ND*ND+j*ND+k]=cec[2][2][1][1]*ep11qrh0[i*ND*ND+j*ND+k]
										 +cec[2][2][2][2]*ep22qrh0[i*ND*ND+j*ND+k]
										 +cec[2][2][3][3]*ep33qrh0[i*ND*ND+j*ND+k];
				sigma33_r[i*ND*ND+j*ND+k]=cec[3][3][1][1]*ep11qrh0[i*ND*ND+j*ND+k]
										 +cec[3][3][2][2]*ep22qrh0[i*ND*ND+j*ND+k]
										 +cec[3][3][3][3]*ep33qrh0[i*ND*ND+j*ND+k];
				sigma12_r[i*ND*ND+j*ND+k]=cec[1][2][1][2]*ep12qrh0[i*ND*ND+j*ND+k]
										 +cec[1][2][2][1]*ep21qrh0[i*ND*ND+j*ND+k];
				sigma21_r[i*ND*ND+j*ND+k]=cec[2][1][1][2]*ep12qrh0[i*ND*ND+j*ND+k]
										 +cec[2][1][2][1]*ep21qrh0[i*ND*ND+j*ND+k];
				sigma13_r[i*ND*ND+j*ND+k]=cec[1][3][1][3]*ep13qrh0[i*ND*ND+j*ND+k]
										 +cec[1][3][3][1]*ep31qrh0[i*ND*ND+j*ND+k];
				sigma31_r[i*ND*ND+j*ND+k]=cec[3][1][1][3]*ep13qrh0[i*ND*ND+j*ND+k]
										 +cec[3][1][3][1]*ep31qrh0[i*ND*ND+j*ND+k];
				sigma23_r[i*ND*ND+j*ND+k]=cec[2][3][2][3]*ep23qrh0[i*ND*ND+j*ND+k]
										 +cec[2][3][3][2]*ep32qrh0[i*ND*ND+j*ND+k];
				sigma32_r[i*ND*ND+j*ND+k]=cec[3][2][2][3]*ep23qrh0[i*ND*ND+j*ND+k]
										 +cec[3][2][3][2]*ep32qrh0[i*ND*ND+j*ND+k];
				//
				sigma11_i[i*ND*ND+j*ND+k]=cec[1][1][1][1]*ep11qih0[i*ND*ND+j*ND+k]
										 +cec[1][1][2][2]*ep22qih0[i*ND*ND+j*ND+k]
										 +cec[1][1][3][3]*ep33qih0[i*ND*ND+j*ND+k];
				sigma22_i[i*ND*ND+j*ND+k]=cec[2][2][1][1]*ep11qih0[i*ND*ND+j*ND+k]
										 +cec[2][2][2][2]*ep22qih0[i*ND*ND+j*ND+k]
										 +cec[2][2][3][3]*ep33qih0[i*ND*ND+j*ND+k];
				sigma33_i[i*ND*ND+j*ND+k]=cec[3][3][1][1]*ep11qih0[i*ND*ND+j*ND+k]
										 +cec[3][3][2][2]*ep22qih0[i*ND*ND+j*ND+k]
										 +cec[3][3][3][3]*ep33qih0[i*ND*ND+j*ND+k];
				sigma12_i[i*ND*ND+j*ND+k]=cec[1][2][1][2]*ep12qih0[i*ND*ND+j*ND+k]
										 +cec[1][2][2][1]*ep21qih0[i*ND*ND+j*ND+k];
				sigma21_i[i*ND*ND+j*ND+k]=cec[2][1][1][2]*ep12qih0[i*ND*ND+j*ND+k]
										 +cec[2][1][2][1]*ep21qih0[i*ND*ND+j*ND+k];
				sigma13_i[i*ND*ND+j*ND+k]=cec[1][3][1][3]*ep13qih0[i*ND*ND+j*ND+k]
										 +cec[1][3][3][1]*ep31qih0[i*ND*ND+j*ND+k];
				sigma31_i[i*ND*ND+j*ND+k]=cec[3][1][1][3]*ep13qih0[i*ND*ND+j*ND+k]
										 +cec[3][1][3][1]*ep31qih0[i*ND*ND+j*ND+k];
				sigma23_i[i*ND*ND+j*ND+k]=cec[2][3][2][3]*ep23qih0[i*ND*ND+j*ND+k]
										 +cec[2][3][3][2]*ep32qih0[i*ND*ND+j*ND+k];
				sigma32_i[i*ND*ND+j*ND+k]=cec[3][2][2][3]*ep23qih0[i*ND*ND+j*ND+k]
										 +cec[3][2][3][2]*ep32qih0[i*ND*ND+j*ND+k];
			}
		}
	}

//***** Calculation of constraint strain variation ec11 [formula (4.9)] *************************************
	iii=1; jjj=1; // ec11 = ec iii jjj
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			for(k=0;k<=ndm;k++){
				if(k<=nd2-1){kk=k;}  if(k>=nd2){kk=k-nd;}
				//
				sigma[1][1]=sigma11_r[i*ND*ND+j*ND+k]; sigma[2][2]=sigma22_r[i*ND*ND+j*ND+k]; sigma[3][3]=sigma33_r[i*ND*ND+j*ND+k];
				sigma[1][2]=sigma[2][1]=sigma12_r[i*ND*ND+j*ND+k];
				sigma[1][3]=sigma[3][1]=sigma13_r[i*ND*ND+j*ND+k];
				sigma[2][3]=sigma[3][2]=sigma23_r[i*ND*ND+j*ND+k];
				xr[i*ND*ND+j*ND+k]=zcij(i, j, k, iii, jjj, ND); // The zcij function includes sigma[][] and cec[][][][]
				//
				sigma[1][1]=sigma11_i[i*ND*ND+j*ND+k]; sigma[2][2]=sigma22_i[i*ND*ND+j*ND+k]; sigma[3][3]=sigma33_i[i*ND*ND+j*ND+k];
				sigma[1][2]=sigma[2][1]=sigma12_i[i*ND*ND+j*ND+k];
				sigma[1][3]=sigma[3][1]=sigma13_i[i*ND*ND+j*ND+k];
				sigma[2][3]=sigma[3][2]=sigma23_i[i*ND*ND+j*ND+k];
				xi[i*ND*ND+j*ND+k]=zcij(i, j, k, iii, jjj, ND); // The zcij function includes sigma[][] and cec[][][][]
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
				ec11[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//***** Calculation of constraint strain variation ec22 [Equation (4.9)] *****************************
	iii=2; jjj=2; // ec11 = ec iii jjj
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			for(k=0;k<=ndm;k++){
				if(k<=nd2-1){kk=k;}  if(k>=nd2){kk=k-nd;}
				//
				sigma[1][1]=sigma11_r[i*ND*ND+j*ND+k]; sigma[2][2]=sigma22_r[i*ND*ND+j*ND+k]; sigma[3][3]=sigma33_r[i*ND*ND+j*ND+k];
				sigma[1][2]=sigma[2][1]=sigma12_r[i*ND*ND+j*ND+k];
				sigma[1][3]=sigma[3][1]=sigma13_r[i*ND*ND+j*ND+k];
				sigma[2][3]=sigma[3][2]=sigma23_r[i*ND*ND+j*ND+k];
				xr[i*ND*ND+j*ND+k]=zcij(i, j, k, iii, jjj, ND); // The zcij function includes sigma[][] and cec[][][][]
				//
				sigma[1][1]=sigma11_i[i*ND*ND+j*ND+k]; sigma[2][2]=sigma22_i[i*ND*ND+j*ND+k]; sigma[3][3]=sigma33_i[i*ND*ND+j*ND+k];
				sigma[1][2]=sigma[2][1]=sigma12_i[i*ND*ND+j*ND+k];
				sigma[1][3]=sigma[3][1]=sigma13_i[i*ND*ND+j*ND+k];
				sigma[2][3]=sigma[3][2]=sigma23_i[i*ND*ND+j*ND+k];
				xi[i*ND*ND+j*ND+k]=zcij(i, j, k, iii, jjj, ND); // The zcij function includes sigma[][] and cec[][][][]
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
				ec22[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//***** Calculation of constraint strain variation ec33 *****************************
	iii=3; jjj=3; // ec11 = ec iii jjj
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			for(k=0;k<=ndm;k++){
				if(k<=nd2-1){kk=k;}  if(k>=nd2){kk=k-nd;}
				//
				sigma[1][1]=sigma11_r[i*ND*ND+j*ND+k]; sigma[2][2]=sigma22_r[i*ND*ND+j*ND+k]; sigma[3][3]=sigma33_r[i*ND*ND+j*ND+k];
				sigma[1][2]=sigma[2][1]=sigma12_r[i*ND*ND+j*ND+k];
				sigma[1][3]=sigma[3][1]=sigma13_r[i*ND*ND+j*ND+k];
				sigma[2][3]=sigma[3][2]=sigma23_r[i*ND*ND+j*ND+k];
				xr[i*ND*ND+j*ND+k]=zcij(i, j, k, iii, jjj, ND); // The zcij function includes sigma[][] and cec[][][][]
				//
				sigma[1][1]=sigma11_i[i*ND*ND+j*ND+k]; sigma[2][2]=sigma22_i[i*ND*ND+j*ND+k]; sigma[3][3]=sigma33_i[i*ND*ND+j*ND+k];
				sigma[1][2]=sigma[2][1]=sigma12_i[i*ND*ND+j*ND+k];
				sigma[1][3]=sigma[3][1]=sigma13_i[i*ND*ND+j*ND+k];
				sigma[2][3]=sigma[3][2]=sigma23_i[i*ND*ND+j*ND+k];
				xi[i*ND*ND+j*ND+k]=zcij(i, j, k, iii, jjj, ND); // The zcij function includes sigma[][] and cec[][][][]
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
				ec33[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//***** Calculation of constraint strain variation ec12 *****************************
	iii=1; jjj=2; // ec11 = ec iii jjj
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			for(k=0;k<=ndm;k++){
				if(k<=nd2-1){kk=k;}  if(k>=nd2){kk=k-nd;}
				//
				sigma[1][1]=sigma11_r[i*ND*ND+j*ND+k]; sigma[2][2]=sigma22_r[i*ND*ND+j*ND+k]; sigma[3][3]=sigma33_r[i*ND*ND+j*ND+k];
				sigma[1][2]=sigma[2][1]=sigma12_r[i*ND*ND+j*ND+k];
				sigma[1][3]=sigma[3][1]=sigma13_r[i*ND*ND+j*ND+k];
				sigma[2][3]=sigma[3][2]=sigma23_r[i*ND*ND+j*ND+k];
				xr[i*ND*ND+j*ND+k]=zcij(i, j, k, iii, jjj, ND); // The zcij function includes sigma[][] and cec[][][][]
				//
				sigma[1][1]=sigma11_i[i*ND*ND+j*ND+k]; sigma[2][2]=sigma22_i[i*ND*ND+j*ND+k]; sigma[3][3]=sigma33_i[i*ND*ND+j*ND+k];
				sigma[1][2]=sigma[2][1]=sigma12_i[i*ND*ND+j*ND+k];
				sigma[1][3]=sigma[3][1]=sigma13_i[i*ND*ND+j*ND+k];
				sigma[2][3]=sigma[3][2]=sigma23_i[i*ND*ND+j*ND+k];
				xi[i*ND*ND+j*ND+k]=zcij(i, j, k, iii, jjj, ND); // The zcij function includes sigma[][] and cec[][][][]
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
				ec12[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//***** Calculation of constraint strain variation ec13 *****************************
	iii=1; jjj=3; // ec11 = ec iii jjj
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			for(k=0;k<=ndm;k++){
				if(k<=nd2-1){kk=k;}  if(k>=nd2){kk=k-nd;}
				//
				sigma[1][1]=sigma11_r[i*ND*ND+j*ND+k]; sigma[2][2]=sigma22_r[i*ND*ND+j*ND+k]; sigma[3][3]=sigma33_r[i*ND*ND+j*ND+k];
				sigma[1][2]=sigma[2][1]=sigma12_r[i*ND*ND+j*ND+k];
				sigma[1][3]=sigma[3][1]=sigma13_r[i*ND*ND+j*ND+k];
				sigma[2][3]=sigma[3][2]=sigma23_r[i*ND*ND+j*ND+k];
				xr[i*ND*ND+j*ND+k]=zcij(i, j, k, iii, jjj, ND); // The zcij function includes sigma[][] and cec[][][][]
				//
				sigma[1][1]=sigma11_i[i*ND*ND+j*ND+k]; sigma[2][2]=sigma22_i[i*ND*ND+j*ND+k]; sigma[3][3]=sigma33_i[i*ND*ND+j*ND+k];
				sigma[1][2]=sigma[2][1]=sigma12_i[i*ND*ND+j*ND+k];
				sigma[1][3]=sigma[3][1]=sigma13_i[i*ND*ND+j*ND+k];
				sigma[2][3]=sigma[3][2]=sigma23_i[i*ND*ND+j*ND+k];
				xi[i*ND*ND+j*ND+k]=zcij(i, j, k, iii, jjj, ND); // The zcij function includes sigma[][] and cec[][][][]
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
				ec13[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//***** Calculation of constraint strain variation ec23 *****************************
	iii=2; jjj=3; // ec11 = ec iii jjj
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			for(k=0;k<=ndm;k++){
				if(k<=nd2-1){kk=k;}  if(k>=nd2){kk=k-nd;}
				//
				sigma[1][1]=sigma11_r[i*ND*ND+j*ND+k]; sigma[2][2]=sigma22_r[i*ND*ND+j*ND+k]; sigma[3][3]=sigma33_r[i*ND*ND+j*ND+k];
				sigma[1][2]=sigma[2][1]=sigma12_r[i*ND*ND+j*ND+k];
				sigma[1][3]=sigma[3][1]=sigma13_r[i*ND*ND+j*ND+k];
				sigma[2][3]=sigma[3][2]=sigma23_r[i*ND*ND+j*ND+k];
				xr[i*ND*ND+j*ND+k]=zcij(i, j, k, iii, jjj, ND); // The zcij function includes sigma[][] and cec[][][][]
				//
				sigma[1][1]=sigma11_i[i*ND*ND+j*ND+k]; sigma[2][2]=sigma22_i[i*ND*ND+j*ND+k]; sigma[3][3]=sigma33_i[i*ND*ND+j*ND+k];
				sigma[1][2]=sigma[2][1]=sigma12_i[i*ND*ND+j*ND+k];
				sigma[1][3]=sigma[3][1]=sigma13_i[i*ND*ND+j*ND+k];
				sigma[2][3]=sigma[3][2]=sigma23_i[i*ND*ND+j*ND+k];
				xi[i*ND*ND+j*ND+k]=zcij(i, j, k, iii, jjj, ND); // The zcij function includes sigma[][] and cec[][][][]
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
				ec23[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//****** Numerical calculation of nonlinear evolution equations ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//
				ip=i+1;  im=i-1;  jp=j+1;  jm=j-1;  kp=k+1;  km=k-1;
				if(i==ndm){ip=0;}  if(i==0){im=ndm;}//periodic boundary conditions
				if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
				if(k==ndm){kp=0;}  if(k==0){km=ndm;}
				//
				s1=s1h[i*ND*ND+j*ND+k]; s1ip=s1h[ip*ND*ND+j*ND+k]; s1im=s1h[im*ND*ND+j*ND+k];
										s1jp=s1h[i*ND*ND+jp*ND+k]; s1jm=s1h[i*ND*ND+jm*ND+k];
										s1kp=s1h[i*ND*ND+j*ND+kp]; s1km=s1h[i*ND*ND+j*ND+km];
				s2=s2h[i*ND*ND+j*ND+k]; s2ip=s2h[ip*ND*ND+j*ND+k]; s2im=s2h[im*ND*ND+j*ND+k];
										s2jp=s2h[i*ND*ND+jp*ND+k]; s2jm=s2h[i*ND*ND+jm*ND+k];
										s2kp=s2h[i*ND*ND+j*ND+kp]; s2km=s2h[i*ND*ND+j*ND+km];
				s3=s3h[i*ND*ND+j*ND+k]; s3ip=s3h[ip*ND*ND+j*ND+k]; s3im=s3h[im*ND*ND+j*ND+k];
										s3jp=s3h[i*ND*ND+jp*ND+k]; s3jm=s3h[i*ND*ND+jm*ND+k];
										s3kp=s3h[i*ND*ND+j*ND+kp]; s3km=s3h[i*ND*ND+j*ND+km];
				//
				//Computation of gradient potential
				s1k_surf = -kapaP*(s1ip+s1im + s1jp+s1jm + s1kp+s1km -6.0*s1);
				s2k_surf = -kapaP*(s2ip+s2im + s2jp+s2jm + s2kp+s2km -6.0*s2);
				s3k_surf = -kapaP*(s3ip+s3im + s3jp+s3jm + s3kp+s3km -6.0*s3);
				//
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
				s1k_ddi = Add0*s1k_dd[i*ND*ND+j*ND+k];
				s2k_ddi = Add0*s2k_dd[i*ND*ND+j*ND+k];
				s3k_ddi = Add0*s3k_dd[i*ND*ND+j*ND+k];
				
				// epsilonT ij = epsilon0 ij - average epsilon0 ij - d(epsilonC ij) - epsilonA ij
				epT[1][1] = ep11h0[i*ND*ND+j*ND+k] - ep0[1][1] - ec11[i*ND*ND+j*ND+k] - ep_a[1][1];
				epT[2][2] = ep22h0[i*ND*ND+j*ND+k] - ep0[2][2] - ec22[i*ND*ND+j*ND+k] - ep_a[2][2];
				epT[3][3] = ep33h0[i*ND*ND+j*ND+k] - ep0[3][3] - ec33[i*ND*ND+j*ND+k] - ep_a[3][3];
	epT[1][2] = epT[2][1] = ep12h0[i*ND*ND+j*ND+k] - ep0[1][2] - ec12[i*ND*ND+j*ND+k] - ep_a[1][2];
	epT[1][3] = epT[3][1] = ep13h0[i*ND*ND+j*ND+k] - ep0[1][3] - ec13[i*ND*ND+j*ND+k] - ep_a[1][3];
	epT[2][3] = epT[3][2] = ep23h0[i*ND*ND+j*ND+k] - ep0[2][3] - ec23[i*ND*ND+j*ND+k] - ep_a[3][3];
				
				//stress-free strain, epsilon 0 ij
				//ep0p[1][1] = Q11*s1*s1 + Q12*s2*s2 + Q13*s3*s3;
				//ep0p[2][2] = Q21*s1*s1 + Q22*s2*s2 + Q23*s3*s3;
				//ep0p[3][3] = Q31*s1*s1 + Q32*s2*s2 + Q33*s3*s3;
				//ep0p[1][2] = Q66*s1*s2;
				//ep0p[1][3] = Q55*s1*s3;
				//ep0p[2][3] = Q44*s2*s3;
				
				//d(Estr)/d(s1) = C(i,j,k,l) * [el_strain(k,l)] * eigen(i,j)_s1
				//d(Estr)/d(s1) = C(i,j,k,l) * [epT(k,l)] * d(epsilon0(polar) ij)/d(s1)
				//d(Estr)/d(s1) = C(i,j,k,l) * d(epsilon0(polar) ij)/d(s1) * [epT(k,l)]
				s1k_str= cec[1][1][1][1]*(2.0*Q11*s1)*epT[1][1]
						+cec[2][2][2][2]*(2.0*Q21*s1)*epT[2][2]
						+cec[3][3][3][3]*(2.0*Q31*s1)*epT[3][3]
						+cec[1][1][2][2]*(2.0*Q11*s1)*epT[2][2]*2.0
						+cec[1][1][3][3]*(2.0*Q11*s1)*epT[3][3]*2.0
						+cec[2][2][3][3]*(2.0*Q21*s2)*epT[3][3]*2.0
						+cec[2][3][2][3]*(0.0       )*epT[2][3]*4.0
						+cec[1][3][1][3]*(    Q55*s3)*epT[1][3]*4.0
						+cec[1][2][1][2]*(    Q66*s2)*epT[1][2]*4.0;
				//
				s2k_str= cec[1][1][1][1]*(2.0*Q12*s2)*epT[1][1]
						+cec[2][2][2][2]*(2.0*Q22*s2)*epT[2][2]
						+cec[3][3][3][3]*(2.0*Q32*s2)*epT[3][3]
						+cec[1][1][2][2]*(2.0*Q12*s2)*epT[2][2]*2.0
						+cec[1][1][3][3]*(2.0*Q12*s2)*epT[3][3]*2.0
						+cec[2][2][3][3]*(2.0*Q22*s2)*epT[3][3]*2.0
						+cec[2][3][2][3]*(    Q44*s3)*epT[2][3]*4.0
						+cec[1][3][1][3]*(0.0       )*epT[1][3]*4.0
						+cec[1][2][1][2]*(    Q66*s1)*epT[1][2]*4.0;
				//
				s3k_str= cec[1][1][1][1]*(2.0*Q13*s3)*epT[1][1]
						+cec[2][2][2][2]*(2.0*Q23*s3)*epT[2][2]
						+cec[3][3][3][3]*(2.0*Q33*s3)*epT[3][3]
						+cec[1][1][2][2]*(2.0*Q13*s3)*epT[2][2]*2.0
						+cec[1][1][3][3]*(2.0*Q13*s3)*epT[3][3]*2.0
						+cec[2][2][3][3]*(2.0*Q23*s3)*epT[3][3]*2.0
						+cec[2][3][2][3]*(    Q44*s2)*epT[2][3]*4.0
						+cec[1][3][1][3]*(    Q55*s1)*epT[1][3]*4.0
						+cec[1][2][1][2]*(0.0       )*epT[1][2]*4.0;

				//Computation of total potential
				// dGc(Pi,T) + d(kappa_p*sum(dPi)) + dEdipole + dEappel + dEstr
				s1k = s1k_chem + s1k_surf + s1k_ddi - E1_ex_0x + s1k_str;
				s2k = s2k_chem + s2k_surf + s2k_ddi - E1_ex_0y + s2k_str;
				s3k = s3k_chem + s3k_surf + s3k_ddi - E1_ex_0z + s3k_str;

				//Evolution equation [equation (4.57)] and time evolution of s-field (explicit method) calculation
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
				//
				s1h[i*ND*ND+j*ND+k]=s1h2[i*ND*ND+j*ND+k];
				s2h[i*ND*ND+j*ND+k]=s2h2[i*ND*ND+j*ND+k];
				s3h[i*ND*ND+j*ND+k]=s3h2[i*ND*ND+j*ND+k];
				//
				if(s1h[i*ND*ND+j*ND+k]>= Peq){s1h[i*ND*ND+j*ND+k]= Peq;}
				if(s1h[i*ND*ND+j*ND+k]<=-Peq){s1h[i*ND*ND+j*ND+k]=-Peq;}
				//
				if(s2h[i*ND*ND+j*ND+k]>= Peq){s2h[i*ND*ND+j*ND+k]= Peq;}
				if(s2h[i*ND*ND+j*ND+k]<=-Peq){s2h[i*ND*ND+j*ND+k]=-Peq;}
				//
				if(s3h[i*ND*ND+j*ND+k]>= Peq){s3h[i*ND*ND+j*ND+k]= Peq;}
				if(s3h[i*ND*ND+j*ND+k]<=-Peq){s3h[i*ND*ND+j*ND+k]=-Peq;}
				//
				ss1ss2=sqrt(s1h[i*ND*ND+j*ND+k]*s1h[i*ND*ND+j*ND+k] + s2h[i*ND*ND+j*ND+k]*s2h[i*ND*ND+j*ND+k]);
				if(ss1ss2>=Peq){
					s1h[i*ND*ND+j*ND+k]=s1h[i*ND*ND+j*ND+k]/ss1ss2*Peq;
					s2h[i*ND*ND+j*ND+k]=s2h[i*ND*ND+j*ND+k]/ss1ss2*Peq;
				}
				//
				ss1ss3=sqrt(s1h[i*ND*ND+j*ND+k]*s1h[i*ND*ND+j*ND+k] + s3h[i*ND*ND+j*ND+k]*s3h[i*ND*ND+j*ND+k]);
				if(ss1ss3>=Peq){
					s1h[i*ND*ND+j*ND+k]=s1h[i*ND*ND+j*ND+k]/ss1ss3*Peq;
					s3h[i*ND*ND+j*ND+k]=s3h[i*ND*ND+j*ND+k]/ss1ss3*Peq;
				}
				//
				ss2ss3=sqrt(s2h[i*ND*ND+j*ND+k]*s2h[i*ND*ND+j*ND+k] + s3h[i*ND*ND+j*ND+k]*s3h[i*ND*ND+j*ND+k]);
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

	temp += dtemp * delt;
end:;
	if(plan) fftw_destroy_plan(plan);	//For FFT
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


//*** Zcij [eq.(5.26) or eq.(II 3.5)] ****************************************
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
	int nd=ND, ndm=ND-1, nd2=ND/2;

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

//************ data save subroutine *******************************
void datsave(double *s1h, double *s2h, double *s3h, int ND)
{
	FILE		*stream;	//Stream pointer setting
	int 		i, j, k;			//integer
	int nd=ND, ndm=ND-1, nd2=ND/2;

	stream = fopen("test_ext.dat", "a");	//Open the write destination file for appending
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
	printf("dip_ext_result%06d.vtk \n",iout);
	sprintf(fName,"dip_ext_result%06d.vtk",iout);
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
	fscanf(datin0, "%lf", &time_ext);	//Read calculation count
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

