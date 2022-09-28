#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#include <fftw3.h>
//g++ martensite_v2_libfftw3.cpp -lfftw3
//#include <mpi.h> //mpi version
//#include <fftw3-mpi.h> //mpi version

#define DRND(x) ((double)(x)/RAND_MAX*rand())

	double RR=8.3145;			//Gas constant
	double PI=3.141592654;		//pi
	
	double c0;					//average of phase field
	double temp;
	double time1;				//count

	int iout;
	int flg;
	
	double zcij(int i0, int j0, int k0, int iii, int jjj, int ND);//Coefficient calculation of elastic function (Fourier space)
	double zuij(int i0, int j0, int k0, int iii, int ND);//Displacement field coefficient calculation (Fourier space)
	
	double eta_c[4][4];			//Lattice mismatch
	double cec[4][4][4][4];		//Elastic constant
	double sigma[4][4];			//Eigen stress (stress when elastically deformed by the amount of Eigen strain)

	void ini_field(double *ch, int ND);	//Initial concentration wave setting subroutine
	void datin(double *ch, int ND);	//Subroutine for initial field reading
	//
	void datsave_old(double *ch, int ND);	//Concentration data save subroutine
	//
	void datsave(double *ch, double *Estr, 
			 double *ep11c, double *ep22c, double *ep33c, double *ep12c, double *ep13c, double *ep23c, 
			 double *sig11, double *sig22, double *sig33, double *sig12, double *sig13, double *sig23, 
			 double *u1, double *u2, double *u3, int ND);	//Data storage subroutine
	void datsave_paraview(double *ch, double *Estr, 
			 double *ep11c, double *ep22c, double *ep33c, double *ep12c, double *ep13c, double *ep23c, 
			 double *sig11, double *sig22, double *sig33, double *sig12, double *sig13, double *sig23, 
			 double *u1, double *u2, double *u3, int ND);	//Data storage subroutine

//******* Main program ******************************************************
int main(void)
{
	int ND;
	int nd, ndm, nd2; 				//Organization split, organization split-1

	double mu_chem, mu_str, mu_surf;
	int    i, j, k, l, iii, jjj;	//integer
	int    Nstep;
	double time1max; 				//maximum time
	double amob_c;					//atomic mobility constant
	double c;						//concentration
	double cddtt;					//concentration increment
	double kapa_c, kapa_c0;			//concentration gradient energy constant
	double al;						//Computational domain
	double b1;						//Differential block size
	double L0, L00;
	double delt;					//ticking
	double Mx, My, Mz;				//mobility

	double cip, cim, cjp, cjm, ckp, ckm;
	double ck1dev, ck2dev;
	double ck0, ckip, ckim, ckjp, ckjm, ckkp, ckkm;
	int   ip, im, jp, jm, kp, km;
	double sumc, dc0;
	double c_flu, c_flu0;			//Magnitude of concentration field fluctuations

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
	double s2k_str_x, s2k_str_y, s2k_str_z;
	double sum11, sum22, sum33;						//spatial integral of s1 and s2
	double sum12, sum13, sum23;						//spatial integral of s1 and s2
	double ep0[4][4];								//Mean value of transformation strain in the structure
	//
	double epT[4][4];
	double s1ddtt, s2ddtt;							//Time variation of s1 and s2 (left side of evolution equation)
	
//********************************************************************************
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
	//
	ND      = int(data[0]);
	delt    = data[1];
	c0      = data[2];
	Mx      = data[3];
	My      = data[4];
	Mz      = data[5];
	al      = data[6];	// [nm]
	amob_c  = data[7];
	time1max= int(data[8]);
	temp    = data[9];	// [K]
	dtemp   = data[10];
	L00     = data[11];	//L0=L00/RR/temp;
	kapa_c0 = data[12];	//kapa_c=kapa_c0/b1/b1/RR/temp;
	c_flu   = data[13];
	flg     = int(data[14]);
	Nstep   = int(data[15]);
	readff  = int(data[16]);
	eta_c[1][1] = data[17];
	eta_c[2][2] = data[18];
	eta_c[3][3] = data[19];
	el_mag  = data[20];
	c11     = data[21]*el_mag;
	c22     = data[22]*el_mag;
	c33     = data[23]*el_mag;
	c12     = data[24]*el_mag; c21=c12;
	c13     = data[25]*el_mag; c31=c13;
	c23     = data[26]*el_mag; c32=c23;
	c44     = data[27]*el_mag;
	c55     = data[28]*el_mag;
	c66     = data[29]*el_mag;
	ep_a[1][1]=data[30];
	ep_a[2][2]=data[31];
	ep_a[3][3]=data[32];
	ep_a[1][2]=ep_a[1][3]=ep_a[2][3]=0.0;
	printf("---------------------------------\n");
	//
	nd=ND;
	ndm=ND-1;
	nd2=ND/2;
	
	//double ch[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	//a[z * ySize * xSize + y * xSize + x], a[i][j][k]=a[i*n*m + j*n + k]
	double *ck  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//diffusion potential
	double *ch  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//Concentration data array in tissue
	double *ch2 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
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
	
	al=al*1.0e-9;	//m

	b1=al/nd;

	time1=0.0;
	
	L0=L00/RR/temp;
	kapa_c=kapa_c0/b1/b1/RR/temp;
	
	c_flu0 = c_flu;

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

//*************************************************************************
	if(readff==0){
		printf("make initial concentration (random) \n");
		ini_field(ch, ND);
	} else {
		printf("read data.dat file \n");
		datin(ch, ND);
	}

//**** start **************************************************************
 plan = fftw_plan_dft_3d(fftsize, fftsize, fftsize, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);	//For FFT
iplan = fftw_plan_dft_3d(fftsize, fftsize, fftsize, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);	//For IFFT
iout = -1;
start: ;
	printf("time: %f \n", time1);
	if((((int)(time1) % Nstep)==0)) {iout = iout + 1;}
	if((((int)(time1) % Nstep)==0)) {datsave_old(ch,ND);}
	if((((int)(time1) % Nstep)==0)) {datsave(ch, Estr, 
			 ep11c, ep22c, ep33c, ep12c, ep13c, ep23c, 
			 sig11, sig22, sig33, sig12, sig13, sig23, 
			 u1, u2, u3, ND);}	//Save the concentration field every fixed repetition count
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(ch, Estr, 
			 ep11c, ep22c, ep33c, ep12c, ep13c, ep23c, 
			 sig11, sig22, sig33, sig12, sig13, sig23, 
			 u1, u2, u3, ND);}	//Save the concentration field every fixed repetition count

//**** Fourier transform of field (phase field) ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
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
				ep11h0[i*ND*ND+j*ND+k]=eta_c[iii][jjj]*ch[i*ND*ND+j*ND+k];
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
				ep22h0[i*ND*ND+j*ND+k]=eta_c[iii][jjj]*ch[i*ND*ND+j*ND+k];
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
				ep33h0[i*ND*ND+j*ND+k]=eta_c[iii][jjj]*ch[i*ND*ND+j*ND+k];
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
				ep12h0[i*ND*ND+j*ND+k]=eta_c[iii][jjj]*ch[i*ND*ND+j*ND+k];
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
				ep13h0[i*ND*ND+j*ND+k]=eta_c[iii][jjj]*ch[i*ND*ND+j*ND+k];
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
				ep23h0[i*ND*ND+j*ND+k]=eta_c[iii][jjj]*ch[i*ND*ND+j*ND+k];
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
											   -(c11+c12+c13)*eta_c[1][1]*(ch[i*ND*ND+j*ND+k]-c0);
				sig22[i*ND*ND+j*ND+k]=c12*ec11[i*ND*ND+j*ND+k]+c22*ec22[i*ND*ND+j*ND+k]+c23*ec33[i*ND*ND+j*ND+k]
											   -(c11+c12+c13)*eta_c[2][2]*(ch[i*ND*ND+j*ND+k]-c0);
				sig33[i*ND*ND+j*ND+k]=c13*ec11[i*ND*ND+j*ND+k]+c23*ec22[i*ND*ND+j*ND+k]+c33*ec33[i*ND*ND+j*ND+k]
											   -(c11+c12+c13)*eta_c[3][3]*(ch[i*ND*ND+j*ND+k]-c0);
				sig12[i*ND*ND+j*ND+k]=2.0*c66*ec12[i*ND*ND+j*ND+k];
				sig13[i*ND*ND+j*ND+k]=2.0*c55*ec13[i*ND*ND+j*ND+k];
				sig23[i*ND*ND+j*ND+k]=2.0*c44*ec23[i*ND*ND+j*ND+k];
				//
				Estr1=1.5*(c11*eta_c[1][1]*eta_c[1][1]
						  +c12*eta_c[2][2]*eta_c[2][2]
						  +c13*eta_c[3][3]*eta_c[3][3])*ch[i*ND*ND+j*ND+k]*ch[i*ND*ND+j*ND+k];
				Estr2=-(c11*ep11c[i*ND*ND+j*ND+k]*eta_c[1][1]
					   +c12*ep22c[i*ND*ND+j*ND+k]*eta_c[2][2]
					   +c13*ep33c[i*ND*ND+j*ND+k]*eta_c[3][3])*ch[i*ND*ND+j*ND+k];
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

//******[Calculation of diffusion potential]********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//
				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==ndm) {ip=0;}  if(i==0) {im=ndm;}
				if(j==ndm) {jp=0;}  if(j==0) {jm=ndm;}
				if(k==ndm) {kp=0;}  if(k==0) {km=ndm;}

				c=ch[i*ND*ND+j*ND+k];

			 	mu_chem=32.0*L0*c*(1.0-c)*(1.0-2.0*c); 		//chemical potential
			 	//mu_chem=L0*(1.0-2.0*c)+log(c)-log(1.0-c); //chemical potential

				mu_surf=-2.*kapa_c*(ch[ip*ND*ND+j*ND+k]+ch[im*ND*ND+j*ND+k]
								   +ch[i*ND*ND+jp*ND+k]+ch[i*ND*ND+jm*ND+k]
      							   +ch[i*ND*ND+j*ND+kp]+ch[i*ND*ND+j*ND+km]
								   -6.0*c);

				ck[i*ND*ND+j*ND+k]= mu_chem+mu_surf;

			}
		}
	}

//******[Time variation of concentration field]**************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//
				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==ndm) {ip=0;}  if(i==0) {im=ndm;}
				if(j==ndm) {jp=0;}  if(j==0) {jm=ndm;}
				if(k==ndm) {kp=0;}  if(k==0) {km=ndm;}
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
				s2k_str_x = cec[1][1][1][1]*eta_c[1][1]*epT[1][1]
						   +cec[1][1][2][2]*eta_c[1][1]*epT[2][2]*1.0
						   +cec[1][1][3][3]*eta_c[1][1]*epT[3][3]*1.0
						   +cec[1][3][1][3]*eta_c[1][3]*epT[1][3]*2.0
						   +cec[1][2][1][2]*eta_c[1][2]*epT[1][2]*2.0;
				
				s2k_str_y = cec[2][2][2][2]*eta_c[2][2]*epT[2][2]
						   +cec[1][1][2][2]*eta_c[1][1]*epT[2][2]*1.0
						   +cec[2][2][3][3]*eta_c[2][2]*epT[3][3]*1.0
						   +cec[2][3][2][3]*eta_c[2][3]*epT[2][3]*2.0
						   +cec[1][2][1][2]*eta_c[1][2]*epT[1][2]*2.0;
				
				s2k_str_z= cec[3][3][3][3]*eta_c[3][3]*epT[3][3]
						  +cec[1][1][3][3]*eta_c[1][1]*epT[3][3]*1.0
						  +cec[2][2][3][3]*eta_c[2][2]*epT[3][3]*1.0
						  +cec[2][3][2][3]*eta_c[2][3]*epT[2][3]*2.0
						  +cec[1][3][1][3]*eta_c[1][3]*epT[1][3]*2.0;
				//
				cddtt=Mx*(ck[ip*ND*ND+j*ND+k]+ck[im*ND*ND+j*ND+k]-2.0*ck[i*ND*ND+j*ND+k] - s2k_str_x)
					 +My*(ck[i*ND*ND+jp*ND+k]+ck[i*ND*ND+jm*ND+k]-2.0*ck[i*ND*ND+j*ND+k] - s2k_str_y)
					 +Mz*(ck[i*ND*ND+j*ND+kp]+ck[i*ND*ND+j*ND+km]-2.0*ck[i*ND*ND+j*ND+k] - s2k_str_z);

				ch2[i*ND*ND+j*ND+k] = ch[i*ND*ND+j*ND+k] + (cddtt+c_flu*(2.0*DRND(1.0)-1.0))*delt;
			}
		}
	}

//******[Concentration field balance correction]**********************************
	sumc=0.;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				sumc+=ch2[i*ND*ND+j*ND+k];
			}
	   }
	}
    dc0=sumc/nd/nd/nd-c0;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
	    	for(k=0;k<=ndm;k++){
				ch[i*ND*ND+j*ND+k]=ch2[i*ND*ND+j*ND+k]-dc0;
				if(ch[i*ND*ND+j*ND+k]>=1.0){ch[i*ND*ND+j*ND+k]=1.0-1.0e-06;}
				if(ch[i*ND*ND+j*ND+k]<=0.0){ch[i*ND*ND+j*ND+k]=1.0e-06;}
	    	}
	    }
	}

//******[time increase]*************************************************
	//if(keypress()){return 0;}
	time1=time1+1.;
	temp += dtemp * delt;
	if (time1<time1max) {goto start;}//Determining if the maximum count has been reached

end:;
	printf("Finished \n");
	if(plan) fftw_destroy_plan(plan);	//For FFT
	if(iplan) fftw_destroy_plan(iplan);	//For IFFT
	fftw_free(in); fftw_free(out);		// For FFT and IFFT
	std::exit(0);
	//return 0;
}

//************[initial concentration wave]*****************************************
void ini_field(double *ch, int ND)
{
	int i, j, k, ii, jj, kk, PN, ipn;
	int x1, y1, z1, r0;
	//int flg;
	int nd=ND, ndm=ND-1, nd2=ND/2; 	//Organization split, organization split-1
	double rnd0, sumc, r; 
 	srand(time(NULL)); // random number initialization
	int nd01, nd02, dd1;

//Initialization
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//ch[i][j][k]=0.0;
				ch[i*ND*ND+j*ND+k]=0.0;
			}
		}
	}


// 1: Spinodal decomposition
// 2: Nucleation-growth decomposition
// 3: 1 sphere (volume fraction is calculated here)
// 4: 1 torus (volume fraction is calculated here)

	//flg=1;//Select one of the above.

	switch(flg){

	case 1: //Spinodal decomposition
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//ch[i][j][k]=c0+0.1*(2.*DRND(1)-1.);
				ch[i*ND*ND+j*ND+k]=c0+0.1*(2.*DRND(1)-1.);
			}
		}
	}
	break;


	case 2: //nucleation
	r0=nd/20;  //put a nucleus of diameter r0
	//PN=50; //place PN nuclei
	PN=c0/( 4./3.*PI*(1./20.)*(1./20.)*(1./20.) );

	for(ipn=1;ipn<=PN;ipn++){
		x1=nd*DRND(1); y1=nd*DRND(1); z1=nd*DRND(1);
		//while(ch[x1][y1][z1]!=0.0){x1=nd*DRND(1); y1=nd*DRND(1); z1=nd*DRND(1);}
		while(ch[x1*ND*ND+y1*ND+z1]!=0.0){x1=nd*DRND(1); y1=nd*DRND(1); z1=nd*DRND(1);}
		for(i=-r0;i<=(ndm+r0);i++){
			ii=i; if(i<0){ii=nd+i;}  if(i>ndm){ii=i-nd;}
			for(j=-r0;j<=(ndm+r0);j++){
				jj=j; if(j<0){jj=nd+j;}  if(j>ndm){jj=j-nd;}
				for(k=-r0;k<=(ndm+r0);k++){
					kk=k; if(k<0){kk=nd+k;}  if(k>ndm){kk=k-nd;}
					r=sqrt( (double(i-x1))*(double(i-x1))
						   +(double(j-y1))*(double(j-y1))
						   +(double(k-z1))*(double(k-z1)) );
					//if(r<=r0){ ch[ii][jj][kk]=1.0; } 
					if(r<=r0){ ch[ii*ND*ND+jj*ND+kk]=1.0; } 
				}
			}
		}
	}

	sumc=0.;
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i][j][k]; } } }
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i*ND*ND+j*ND+k]; } } }
		//c0=sumc/nd/nd/nd;
	break;


	case 3: //1 ball
	nd01=nd/3;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				dd1=(i-nd2)*(i-nd2)+(j-nd2)*(j-nd2)+(k-nd2)*(k-nd2);
				//if(dd1<(nd01*nd01)){ch[i][j][k]=0.9;} else{ch[i][j][k]=0.2;}
				if(dd1<(nd01*nd01)){ch[i*ND*ND+j*ND+k]=0.9;} else{ch[i*ND*ND+j*ND+k]=0.2;}
			}
		}
	}
	sumc=0.;
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i][j][k]; } } }
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i*ND*ND+j*ND+k]; } } }
		c0=sumc/nd/nd/nd;
	break;


	case 4: //1 torus
	nd01=nd/3; nd02=nd/6;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=nd2-4;k<=nd2+4;k++){
				dd1=(i-nd2)*(i-nd2)+(j-nd2)*(j-nd2);
				//if(dd1<(nd01*nd01)){ch[i][j][k]=0.9;} else{ch[i][j][k]=0.2;}
				//if(dd1<(nd02*nd02)){ch[i][j][k]=0.2;}
				if(dd1<(nd01*nd01)){ch[i*ND*ND+j*ND+k]=0.9;} else{ch[i*ND*ND+j*ND+k]=0.2;}
				if(dd1<(nd02*nd02)){ch[i*ND*ND+j*ND+k]=0.2;}
			}
		}
	}
	sumc=0.;
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i][j][k]; } } }
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i*ND*ND+j*ND+k]; } } }
		c0=sumc/nd/nd/nd;
	break;
	}//switch
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

//************[Saving concentration wave data]************************************
void datsave_old(double *ch, int ND)
{
	FILE		*stream;
	char	fName[256];
	int 		i, j, k;
	int ndm=ND-1;

	sprintf(fName,"data_%06d.dat",iout);
	//stream = fopen("test3D_rod.dat", "a");
	stream = fopen(fName, "w");
	fprintf(stream, "%d %d %d \n", ndm, ndm, ndm);
	fprintf(stream, "%e\n", time1);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//fprintf(stream, "%e  ", ch[i][j][k]);
				fprintf(stream, "%e  ", ch[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(stream, "\n");
	fclose(stream);
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
	if (ndm != rndxm){
		printf("data size is mismatch \n");
		printf("Please, change ND= %d in parameters.txt \n", rndxm);
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