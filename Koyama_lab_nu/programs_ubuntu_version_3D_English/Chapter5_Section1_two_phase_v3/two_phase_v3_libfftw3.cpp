#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#include <fftw3.h>
//g++ martensite_v2_libfftw3.cpp -lfftw3
//#include <mpi.h> //mpi version
//#include <fftw3-mpi.h> //mpi version

#define DRND(x) ((double)(x)/RAND_MAX*rand())	//Random number setting

	double PI=3.141592654;		//Pi
	double rr=8.3145;			//gas constant
	double temp;				//absolute temperature
	double time1;				//number of calculation counts (proportional to time)
	double c2a;					//Average composition (B component)
	int Nstep, iout;

	double c0;					//Phase field mean
	double eta_c[4][4];			//lattice mismatch
	double cec[4][4][4][4];		//elastic constant
	double sigma[4][4];			//Eigen stress (stress when elastically deformed by Eigen strain)

	double zcij(int i0, int j0, int k0, int iii, int jjj, int ND);//Coefficient calculation of elastic energy function (Fourier space)
	double zuij(int i0, int j0, int k0, int iii, int ND);//Coefficient calculation of displacement field (Fourier space)

	void ini000(double *c2h, int ND);				//Initial concentration profile setting subroutine
	void datin(double *c2h, int ND);	//Initial field loading subroutine
	void datsave2(double *c2h, int ND);
	void datsave(double *ph, double *c2h, double *Estr, 
			 double *ep11c, double *ep22c, double *ep33c, double *ep12c, double *ep13c, double *ep23c, 
			 double *sig11, double *sig22, double *sig33, double *sig12, double *sig13, double *sig23, 
			 double *u1, double *u2, double *u3, int ND);	//Data storage subroutine
	void datsave_paraview(double *ph, double *c2h, double *Estr, 
			 double *ep11c, double *ep22c, double *ep33c, double *ep12c, double *ep13c, double *ep23c, 
			 double *sig11, double *sig22, double *sig33, double *sig12, double *sig13, double *sig23, 
			 double *u1, double *u2, double *u3, int ND);	//Data storage subroutine

//******* main program ******************************************
int main(void)
{
	int ND, nd, ndm;
	
	double c1, c2;				//local concentration
	double c2max, c2min;		//Maximum and minimum density
	double c2k_chem, c2k_su;	//local potential
	double dc2a, sumc2;			//Variables used for concentration field balance calculation
	double dakd2;				//Second derivative of the diffusion potential
	double c2ddtt;				//Time variation of concentration field

	int   i, j, k, l, iii, jjj;;	//integer
	int   ip, im, jp, jm, kp, km;//(i+1),(i-1),(j+1),(j-1),(k+1),(k-1)
	double al, b1, rtemp, delt;	//Length of one side of computational domain, length of one side of difference block, RT, time step
	double time1max;			//Maximum calculation count (calculation end count)
	double cmob22;				//mobility

	double om_12, om_12e;		//interaction parameter
	double kapa_c2, kapa_c2c;	//concentration gradient energy coefficient

	double dtemp;
	int    readff;
	
	double el_mag;									//Working variables for elastic modulus
	double c11, c22, c33;							//elastic constant
	double c12, c21, c13, c31, c23, c32;
	double c44, c55, c66;
	//double ep000;									//lattice mismatch
	
	double Estr1, Estr2, Estr3, Estr4, Estr5;		//Working variables for elastic strain energy calculation
	
	double ep_a[4][4];								//elastic strain related external force
	//
	double s2k_str;
	double sum11, sum22, sum33;						//spatial integral of s1 and s2
	double sum12, sum13, sum23;						//spatial integral of s1 and s2
	double ep0[4][4];								//Mean value of transformation strain in the structure
	//
	double epT[4][4];
	double s1ddtt, s2ddtt;							//Time variation of s1 and s2 (left side of evolution equation)
	
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
	c2a     = data[1];
	delt    = data[2];	// timestep
	temp    = data[3];	// [K]
	dtemp   = data[4];
	al      = data[5];	// [nm]
	cmob22  = data[6];
	om_12e  = data[7];
	kapa_c2c= data[8];
	time1max= int(data[9]);
	Nstep   = int(data[10]);
	readff  = int(data[11]);
	eta_c[1][1] = data[12];
	eta_c[2][2] = data[13];
	eta_c[3][3] = data[14];
	el_mag  = data[15];
	c11     = data[16]*el_mag;
	c22     = data[17]*el_mag;
	c33     = data[18]*el_mag;
	c12     = data[19]*el_mag; c21=c12;
	c13     = data[20]*el_mag; c31=c13;
	c23     = data[21]*el_mag; c32=c23;
	c44     = data[22]*el_mag;
	c55     = data[23]*el_mag;
	c66     = data[24]*el_mag;
	ep_a[1][1]=data[25];
	ep_a[2][2]=data[26];
	ep_a[3][3]=data[27];
	ep_a[1][2]=ep_a[1][3]=ep_a[2][3]=0.0;
	printf("---------------------------------\n");
	//
	nd=ND;
	ndm=ND-1;
	//
	//double ch[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	//a[z * ySize * xSize + y * xSize + x], a[i][j][k]=a[i*n*m + j*n + k]
	double *ph   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//local volume fraction
	double *c2h  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//topical composition
	double *c2h2 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//auxiliary matrix for the local concentration field
	double *c2k  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//local potential
	//
	double *ep11h0   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Transformational strain in tissue
	double *ep22h0   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Transformational strain in tissue
	double *ep33h0   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Transformational strain in tissue
	double *ep12h0   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Transformational strain in tissue
	double *ep13h0   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Transformational strain in tissue
	double *ep23h0   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Transformational strain in tissue
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
	//After calc.: fftw_destroy_plan(plan); fftw_free(in); fftw_free(out);

	rtemp=rr*temp;				//RT

	al=al*1.0E-09;				//Side length of calculation area (nm)
	b1=al/(double)ND;			//Length of one side of difference block

	om_12=om_12e/rtemp; 		//Interaction parameters (in J/mol, non-dimensionalized at RT)

	kapa_c2=kapa_c2c/b1/b1/rtemp;//gradient energy coefficient (in Jm^2/mol, non-dimensionalized by RT and b1^2)

	time1=0.0;

//*** Setting the lattice mismatch ****************************
	//ep000=0.05;
	//eta_c[1][1]=eta_c[2][2]=eta_c[3][3]=ep000; 			//lattice mismatch
	eta_c[1][2]=eta_c[2][1]=0.0;
	eta_c[1][3]=eta_c[3][1]=0.0;
	eta_c[2][3]=eta_c[3][2]=0.0;

//***** Elastic constant of Fe(bcc) ****************************
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
	cec[2][2][2][2]=c22;
	cec[3][3][3][3]=c33;
	cec[2][3][2][3]=cec[2][3][3][2]=cec[3][2][2][3]=cec[3][2][3][2]=c44;
	cec[1][3][1][3]=cec[1][3][3][1]=cec[3][1][1][3]=cec[3][1][3][1]=c55;
	cec[1][2][1][2]=cec[1][2][2][1]=cec[2][1][1][2]=cec[2][1][2][1]=c66;
	cec[1][1][2][2]=cec[2][2][1][1]=c12;
	cec[1][1][3][3]=cec[3][3][1][1]=c13;
	cec[2][2][3][3]=cec[3][3][2][2]=c23;

//--- Eigen stress (stress when elastically deformed by Eigen strain)j--------------
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

//*** Initial concentration field setting and drawing window display *****************************************

	if(readff==0){
		printf("make initial concentration (random) \n");
		ini000(c2h, ND);	//Setting the initial concentration field
	} else {
		printf("read data.dat file \n");
		datin(c2h, ND);
	}
	
//**** Simulation start ******************************
 plan = fftw_plan_dft_3d(fftsize, fftsize, fftsize, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);	//For FFT
iplan = fftw_plan_dft_3d(fftsize, fftsize, fftsize, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);	//For IFFT
iout = -1;
start: ;

	if((((int)(time1) % Nstep)==0)) {iout = iout + 1;}
	if((((int)(time1) % Nstep)==0)) {datsave2(c2h, ND);}
	if((((int)(time1) % Nstep)==0)) {datsave(ph, c2h, Estr, 
			 ep11c, ep22c, ep33c, ep12c, ep13c, ep23c, 
			 sig11, sig22, sig33, sig12, sig13, sig23, 
			 u1, u2, u3, ND);}	//Save the concentration field every fixed repetition count
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(ph, c2h, Estr, 
			 ep11c, ep22c, ep33c, ep12c, ep13c, ep23c, 
			 sig11, sig22, sig33, sig12, sig13, sig23, 
			 u1, u2, u3, ND);}	//Save the concentration field every fixed repetition count

//**** Fourier transform of field (phase field) ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k]=c2h[i*ND*ND+j*ND+k];
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
				qrh1[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
				qih1[i*ND*ND+j*ND+k]=xi[i*ND*ND+j*ND+k];
			}
		}
	}
	qrh1[0]=qih1[0]=0.0;

//***** Calculation of total strain variation *************************************
//--- ec11 calculator ---
	iii=1; jjj=1; // ec11 = ec iii jjj
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ep11h0[i*ND*ND+j*ND+k]=eta_c[iii][jjj]*c2h[i*ND*ND+j*ND+k];
				//
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
				ec11[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- ec22 calculator ---
	iii=2; jjj=2; // ec22 = ec iii jjj
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ep22h0[i*ND*ND+j*ND+k]=eta_c[iii][jjj]*c2h[i*ND*ND+j*ND+k];
				//
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
				ec22[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- ec33 calculator ---
	iii=3; jjj=3; // ec33 = ec iii jjj
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ep33h0[i*ND*ND+j*ND+k]=eta_c[iii][jjj]*c2h[i*ND*ND+j*ND+k];
				//
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
				ec33[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- ec12 calculator ---
	iii=1; jjj=2; // ec12 = ec iii jjj
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ep12h0[i*ND*ND+j*ND+k]=eta_c[iii][jjj]*c2h[i*ND*ND+j*ND+k];
				//
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
				ec12[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- ec13 calculator ---
	iii=1; jjj=3; // ec13 = ec iii jjj
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ep13h0[i*ND*ND+j*ND+k]=eta_c[iii][jjj]*c2h[i*ND*ND+j*ND+k];
				//
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
				ec13[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- ec23 calculator ---
	iii=2; jjj=3; // ec23 = ec iii jjj
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ep23h0[i*ND*ND+j*ND+k]=eta_c[iii][jjj]*c2h[i*ND*ND+j*ND+k];
				//
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
				ec23[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

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
	ep0[1][2]=sum12/nd/nd/nd; ep0[1][3]=sum13/nd/nd/nd; ep0[2][3]=sum23/nd/nd/nd;

//***** Total strain field, stress field, and elastic strain energy field calculation *************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//
				ep11c[i*ND*ND+j*ND+k]=ec11[i*ND*ND+j*ND+k]+eta_c[1][1]*c0;
				ep22c[i*ND*ND+j*ND+k]=ec22[i*ND*ND+j*ND+k]+eta_c[2][2]*c0;
				ep33c[i*ND*ND+j*ND+k]=ec33[i*ND*ND+j*ND+k]+eta_c[3][3]*c0;
				ep12c[i*ND*ND+j*ND+k]=ec12[i*ND*ND+j*ND+k];
				ep13c[i*ND*ND+j*ND+k]=ec13[i*ND*ND+j*ND+k];
				ep23c[i*ND*ND+j*ND+k]=ec23[i*ND*ND+j*ND+k];
				//
				sig11[i*ND*ND+j*ND+k]=c11*ec11[i*ND*ND+j*ND+k]+c12*ec22[i*ND*ND+j*ND+k]+c13*ec33[i*ND*ND+j*ND+k]
											   -(c11+c12+c13)*eta_c[1][1]*(c2h[i*ND*ND+j*ND+k]-c0);
				sig22[i*ND*ND+j*ND+k]=c12*ec11[i*ND*ND+j*ND+k]+c22*ec22[i*ND*ND+j*ND+k]+c23*ec33[i*ND*ND+j*ND+k]
											   -(c11+c12+c13)*eta_c[2][2]*(c2h[i*ND*ND+j*ND+k]-c0);
				sig33[i*ND*ND+j*ND+k]=c13*ec11[i*ND*ND+j*ND+k]+c23*ec22[i*ND*ND+j*ND+k]+c33*ec33[i*ND*ND+j*ND+k]
											   -(c11+c12+c13)*eta_c[3][3]*(c2h[i*ND*ND+j*ND+k]-c0);
				sig12[i*ND*ND+j*ND+k]=2.0*c66*ec12[i*ND*ND+j*ND+k];
				sig13[i*ND*ND+j*ND+k]=2.0*c55*ec13[i*ND*ND+j*ND+k];
				sig23[i*ND*ND+j*ND+k]=2.0*c44*ec23[i*ND*ND+j*ND+k];
				//
				Estr1=1.5*(c11*eta_c[1][1]*eta_c[1][1]
						  +c12*eta_c[2][2]*eta_c[2][2]
						  +c13*eta_c[3][3]*eta_c[3][3])*c2h[i*ND*ND+j*ND+k]*c2h[i*ND*ND+j*ND+k];
				Estr2=-(c11*ep11c[i*ND*ND+j*ND+k]*eta_c[1][1]
					   +c12*ep22c[i*ND*ND+j*ND+k]*eta_c[2][2]
					   +c13*ep33c[i*ND*ND+j*ND+k]*eta_c[3][3])*c2h[i*ND*ND+j*ND+k];
				Estr3=0.5*c11*ep11c[i*ND*ND+j*ND+k]*ep11c[i*ND*ND+j*ND+k]
					 +0.5*c22*ep22c[i*ND*ND+j*ND+k]*ep22c[i*ND*ND+j*ND+k]
					 +0.5*c33*ep33c[i*ND*ND+j*ND+k]*ep33c[i*ND*ND+j*ND+k];
				Estr4=c12*ep11c[i*ND*ND+j*ND+k]*ep22c[i*ND*ND+j*ND+k]
					 +c13*ep11c[i*ND*ND+j*ND+k]*ep33c[i*ND*ND+j*ND+k]
					 +c23*ep22c[i*ND*ND+j*ND+k]*ep33c[i*ND*ND+j*ND+k];
				Estr5=2.0*c66*ep12c[i*ND*ND+j*ND+k]*ep12c[i*ND*ND+j*ND+k]
					 +2.0*c55*ep13c[i*ND*ND+j*ND+k]*ep13c[i*ND*ND+j*ND+k]
					 +2.0*c44*ep23c[i*ND*ND+j*ND+k]*ep23c[i*ND*ND+j*ND+k];
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
				xr[i*ND*ND+j*ND+k] = out[i*ND*ND+j*ND+k][0]/(fftsize*fftsize*fftsize); //For IFFT (real)
				//
				u1[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- Calculation of u2 ---
	iii=2;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
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
				xr[i*ND*ND+j*ND+k] = out[i*ND*ND+j*ND+k][0]/(fftsize*fftsize*fftsize); //For IFFT (real)
				//
				u2[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//--- Calculation of u3 ---
	iii=3;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
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
				xr[i*ND*ND+j*ND+k] = out[i*ND*ND+j*ND+k][0]/(fftsize*fftsize*fftsize); //For IFFT (real)
				//
				u3[i*ND*ND+j*ND+k]=xr[i*ND*ND+j*ND+k];
			}
		}
	}

//***** Computation of the potential field (see Section 5.2) ***********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				
				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==ndm){ip=0;} if(i==0){im=ndm;}	//periodic boundary conditions
				if(j==ndm){jp=0;} if(j==0){jm=ndm;}	//periodic boundary conditions
				if(k==ndm){kp=0;} if(k==0){km=ndm;}	//periodic boundary conditions

				c2=c2h[i*ND*ND+j*ND+k];  c1=1.0-c2;//local concentration field

				c2k_chem=om_12*(c1-c2)+(log(c2)-log(c1));//chemical diffusion potential
				//
				c2k_su=-2.*kapa_c2*(c2h[ip*ND*ND+j*ND+k]+c2h[im*ND*ND+j*ND+k]
								   +c2h[i*ND*ND+jp*ND+k]+c2h[i*ND*ND+jm*ND+k]
      							   +c2h[i*ND*ND+j*ND+kp]+c2h[i*ND*ND+j*ND+km]
								   -6.0*c2);//gradient potential

				c2k[i*ND*ND+j*ND+k]=c2k_chem+c2k_su;//diffusion potential
			}
		}
	}

//***** Evolution equation calculation (see section 5.2) **********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//
				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==ndm){ip=0;} if(i==0){im=ndm;}	//periodic boundary conditions
				if(j==ndm){jp=0;} if(j==0){jm=ndm;}	//periodic boundary conditions
				if(k==ndm){kp=0;} if(k==0){km=ndm;}	//periodic boundary conditions
				//
				dakd2=c2k[ip*ND*ND+j*ND+k]+c2k[im*ND*ND+j*ND+k]
					 +c2k[i*ND*ND+jp*ND+k]+c2k[i*ND*ND+jm*ND+k]
      				 +c2k[i*ND*ND+j*ND+kp]+c2k[i*ND*ND+j*ND+km]
					 -6.0*c2k[i*ND*ND+j*ND+k];	//Second derivative of the diffusion potential
				//
				//****** Calculation of elastic potential [equation (4.8)] ********************************
				// epsilonT ij = epsilon0 ij - average epsilon0 ij - d(epsilonC ij) - epsilonA ij
				epT[1][1] = ep11h0[i*ND*ND+j*ND+k] - ep0[1][1] - ec11[i*ND*ND+j*ND+k] - ep_a[1][1];
				epT[2][2] = ep22h0[i*ND*ND+j*ND+k] - ep0[2][2] - ec22[i*ND*ND+j*ND+k] - ep_a[2][2];
				epT[3][3] = ep33h0[i*ND*ND+j*ND+k] - ep0[3][3] - ec33[i*ND*ND+j*ND+k] - ep_a[3][3];
	epT[1][2] = epT[2][1] = ep12h0[i*ND*ND+j*ND+k] - ep0[1][2] - ec12[i*ND*ND+j*ND+k] - ep_a[1][2];
	epT[1][3] = epT[3][1] = ep13h0[i*ND*ND+j*ND+k] - ep0[1][3] - ec13[i*ND*ND+j*ND+k] - ep_a[1][3];
	epT[2][3] = epT[3][2] = ep23h0[i*ND*ND+j*ND+k] - ep0[2][3] - ec23[i*ND*ND+j*ND+k] - ep_a[3][3];
				//
				s2k_str= cec[1][1][1][1]*eta_c[1][1]*epT[1][1]
						+cec[2][2][2][2]*eta_c[2][2]*epT[2][2]
						+cec[3][3][3][3]*eta_c[3][3]*epT[3][3]
						+cec[1][1][2][2]*eta_c[1][1]*epT[2][2]*2.0
						+cec[1][1][3][3]*eta_c[1][1]*epT[3][3]*2.0
						+cec[2][2][3][3]*eta_c[2][2]*epT[3][3]*2.0
						+cec[2][3][2][3]*eta_c[2][3]*epT[2][3]*4.0
						+cec[1][3][1][3]*eta_c[1][3]*epT[1][3]*4.0
						+cec[1][2][1][2]*eta_c[1][2]*epT[1][2]*4.0;
				//
				c2ddtt = cmob22*(dakd2 - s2k_str);//diffusion equation
				//
				c2h2[i*ND*ND+j*ND+k]=c2h[i*ND*ND+j*ND+k]+c2ddtt*delt;//Time evolution of the concentration field
				//
				if(c2h[i*ND*ND+j*ND+k]>=1.0){c2h[i*ND*ND+j*ND+k]=1.0-1.0e-06;}//Concentration field correction
				if(c2h[i*ND*ND+j*ND+k]<=0.0){c2h[i*ND*ND+j*ND+k]=1.0e-06;}
			}
		}
	}

//*** Concentration field balance correction ***********************************************
	sumc2=0.;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				sumc2+=c2h2[i*ND*ND+j*ND+k];
			}
		}
	}
	dc2a=sumc2/ND/ND/ND-c2a;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				c2h[i*ND*ND+j*ND+k]=c2h2[i*ND*ND+j*ND+k]-dc2a;
			}
		}
	}

//*** Setting the phase field (see section 5.3.1) ***************************************
	c2max=0.0; c2min=1.0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				if(c2max<=c2h[i*ND*ND+j*ND+k]){c2max=c2h[i*ND*ND+j*ND+k];}
				if(c2min>=c2h[i*ND*ND+j*ND+k]){c2min=c2h[i*ND*ND+j*ND+k];}
			}
		} 
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ph[i*ND*ND+j*ND+k]=(c2h[i*ND*ND+j*ND+k]-c2min)/(c2max-c2min);
			}
		} 
	}

//*********************************************************************
	//if(keypress()){return 0;}	//Waiting for key
	time1=time1+1.0;
	temp += dtemp * delt;
	if(time1<time1max){goto start;}//Determining if the maximum count has been reached
	
end:;
	printf("Finished \n");
	if(plan) fftw_destroy_plan(plan);	//For FFT
	if(iplan) fftw_destroy_plan(iplan);	//For IFFT
	fftw_free(in); fftw_free(out);		// For FFT and IFFT
	std::exit(0);
	//return 0;
}

//************ Initial concentration field setting subroutine *************
void ini000(double *c2h, int ND)
{
	int i, j, k, id;
 	//srand(time(NULL));//random number seeding
	int ndm=ND-1;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//Set the concentration field with random numbers up to +/- 1%
				c2h[i*ND*ND+j*ND+k]=c2a+0.01*(2.0*DRND(1)-1.0);
			}
		}
	}

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
void datsave(double *ph, double *c2h, double *Estr, 
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
				fprintf(stream, "%e  ", ph[i*ND*ND+j*ND+k]);
				fprintf(stream, "%e  ", c2h[i*ND*ND+j*ND+k]);
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

void datsave_paraview(double *ph, double *c2h, double *Estr, 
			 double *ep11c, double *ep22c, double *ep33c, double *ep12c, double *ep13c, double *ep23c, 
			 double *sig11, double *sig22, double *sig33, double *sig12, double *sig13, double *sig23, 
			 double *u1, double *u2, double *u3, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j, k;
	int nd=ND, ndm=ND-1, nd2=ND/2;
	
	//iout = iout + 1;
	printf("paraview output no.%06d \n",iout);
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
	fprintf(fp,"SCALARS Phase_field float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", ph[i][j][k]);
				fprintf(fp,"%10.6f\n", ph[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS concentration float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", c2h[i][j][k]);
				fprintf(fp,"%10.6f\n", c2h[i*ND*ND+j*ND+k]);
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

//************ data save subroutine *******************************
void datsave2(double *c2h, int ND)
{
	FILE	*stream;		//Stream pointer setting
	char	fName[256];
	int 	i, j, k;		//integer
	int nd=ND, ndm=ND-1;
	int ndxm=ndm, ndym=ndm, ndzm=ndm;

	sprintf(fName,"data_%06d.dat",iout);
	//stream = fopen("test.dat", "a");	//Open the file to write to in append mode
	stream = fopen(fName, "w");
	fprintf(stream, "%d %d %d \n", ndm, ndm, ndm);
	fprintf(stream, "%e  \n", time1);	//Saving calculation counts
	for(i=0;i<=ndxm;i++){
		for(j=0;j<=ndym;j++){
			for(k=0;k<=ndzm;k++){
				fprintf(stream, "%e   ", c2h[i*ND*ND+j*ND+k]);//Save phase field
			}
		}
	}
	fprintf(stream, "\n");	//writing a newline
	fclose(stream);					//close file
}

//*********** data entry subroutine **************************
void datin(double *c2h, int ND)
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
				fscanf(datin0, "%lf", &c2h[i*ND*ND+j*ND+k]);//Read field data
				c00+=c2h[i*ND*ND+j*ND+k];
			}
		}
	}
	c0=c00/nd/nd/nd;//Average value of the field
	printf("c0=  %f  \n", c0);

	fclose(datin0);
}
