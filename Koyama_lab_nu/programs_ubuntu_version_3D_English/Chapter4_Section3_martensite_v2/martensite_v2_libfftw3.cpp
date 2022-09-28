#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#include <fftw3.h>
//g++ martensite_v2_libfftw3.cpp -lfftw3
//#include <mpi.h> //mpi version
//#include <fftw3-mpi.h> //mpi version

#define DRND(x) ((double)(x)/RAND_MAX*rand())//Random number function setting
	
	double PI=3.14159;			//Pi
	double rr=8.3145;			//gas constant
	double time1;				//Number of calculation counts (proportional to time)
	int Nstep, iout;
	
	double cec[4][4][4][4];		//Elastic constant
	double sigma[4][4];			//Eigen stress (stress when elastically deformed by the amount of Eigen strain)
	
	double zcij(int i0, int j0, int k0, int iii, int jjj, int ND);//Coefficient calculation of elastic function (Fourier space)
	
	void ini000(double *s1h, double *s2h, int ND);	//Initial field setup subroutine
	void datsave(double *s1h, double *s2h, int ND);			//data save subroutine
	void datsave_paraview(double *s1h, double *s2h, int ND);//data save subroutine

//******* main program ******************************************
int main(void)
{
    int ND, IG;
	int nd, ndm, nd2, ig;
	
	double s1, s2;									//Phase field of martensite
	//
	double s1k_chem, s1k_str;						//potential
	double s2k_chem, s2k_str;						//potential
	//
	double c11, c22, c33, c12, c13, c23, c44, c55, c66; //elastic constant
	double c11c, c22c, c33c, c12c, c13c, c23c, c44c, c55c, c66c; //elastic constant
	//
	double eta_s1[4][4], eta_s2[4][4];				//Eigen distortion component
	double ep_a[4][4];								//elastic strain related external force
	//
	double sum11, sum22, sum33;						//spatial integral of s1 and s2
	double sum12, sum13, sum23;						//spatial integral of s1 and s2
	double ep0[4][4];								//Mean value of transformation strain in the structure
	//
	double epT[4][4];
	double s1ddtt, s2ddtt;							//Time variation of s1 and s2 (left side of evolution equation)
	double el_fac;									//Normalization constant of elastic constant

	int   i, j, k, l, ii, jj, kk, iii, jjj;			//integer
	int   p, q, m, n;								//integer
	int   ip, im, jp, jm, kp, km;					//integer
	int   Nstep;									//integer
	//
	double al, temp, delt;							//Computational domain, temperature, time step
	double time1max;								//Maximum calculation count (calculation end count)
	double b1, vm0, atom_n;							//Side length of difference block, molar volume, number of atoms in unit cell
	double smob;									//Mobility (relaxation coefficient of crystal transformation)
	double nxx, nyy, nxy, alnn;						//product of basis vectors in Fourier space, norm

	double AA0, AA1, AA2, AA3;						//factor in Gibbs energy
	double AA0e;									//factor in Gibbs energy
	//
	double a1_c, b1_c, c1_c;						//lattice constant
	double a1_t, b1_t, c1_t;						//lattice constant
	//
	double kappa_s1, kappa_s2;						//gradient energy factor
	double kappa_s1c, kappa_s2c;					//gradient energy factor
	//
	double ds_fac;									//Fluctuation coefficient of crystal transformation

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
	c22c    = data[16];
	c33c    = data[17];
	c12c    = data[18];
	c13c    = data[19];
	c23c    = data[20];
	c44c    = data[21];
	c55c    = data[22];
	c66c    = data[23];
	ep_a[1][1]=data[24];
	ep_a[2][2]=data[25];
	ep_a[3][3]=data[26];
	eta_s1[1][1]=data[27];
	eta_s1[2][2]=data[28];
	eta_s1[3][3]=data[29];
	eta_s2[1][1]=data[30];
	eta_s2[2][2]=data[31];
	eta_s2[3][3]=data[32];
	printf("---------------------------------\n");
	//
	IG=int(log2(ND));
	nd=ND; 					
	ndm=ND-1;				//define ND-1
	nd2=ND/2;				//define ND/2 for FFT
	ig=IG;					//2^ig=ND
	//
	double *s1h      = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Phase field of martensite
	double *s2h      = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Phase field of martensite
	//
	double *ep11h0   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Transformational strain in tissue
	double *ep22h0   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Transformational strain in tissue
	double *ep33h0   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Transformational strain in tissue
	double *ep12h0   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Transformational strain in tissue
	double *ep13h0   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Transformational strain in tissue
	double *ep23h0   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Transformational strain in tissue
	//
	double *ep11qrh0 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Fourier transform of constraint strain variation
	double *ep22qrh0 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Fourier transform of constraint strain variation
	double *ep33qrh0 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Fourier transform of constraint strain variation
	double *ep12qrh0 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Fourier transform of constraint strain variation
	double *ep21qrh0 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Fourier transform of constraint strain variation
	double *ep13qrh0 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Fourier transform of constraint strain variation
	double *ep31qrh0 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Fourier transform of constraint strain variation
	double *ep23qrh0 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Fourier transform of constraint strain variation
	double *ep32qrh0 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Fourier transform of constraint strain variation
	//
	double *ep11qih0 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Fourier transform of constraint strain variation
	double *ep22qih0 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Fourier transform of constraint strain variation
	double *ep33qih0 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Fourier transform of constraint strain variation
	double *ep12qih0 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Fourier transform of constraint strain variation
	double *ep21qih0 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Fourier transform of constraint strain variation
	double *ep13qih0 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Fourier transform of constraint strain variation
	double *ep31qih0 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Fourier transform of constraint strain variation
	double *ep23qih0 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Fourier transform of constraint strain variation
	double *ep32qih0 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Fourier transform of constraint strain variation
	//
	double *s1k_su   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//potential
	double *s2k_su   = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//potential
	//
	double *ec11     = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Constraint strain variation (real space)
	double *ec22     = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Constraint strain variation (real space)
	double *ec33     = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Constraint strain variation (real space)
	double *ec12     = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Constraint strain variation (real space)
	double *ec13     = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Constraint strain variation (real space)
	double *ec23     = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//Constraint strain variation (real space)
	//
	double *xi       = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//array of real and imaginary parts of the Fourier transform
	double *xr       = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));//array of real and imaginary parts of the Fourier transform
	//
	double *sigma11_r= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma22_r= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma33_r= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma12_r= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma21_r= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma13_r= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma31_r= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma23_r= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma32_r= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	//
	double *sigma11_i= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma22_i= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma33_i= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma12_i= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma21_i= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma13_i= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma31_i= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma23_i= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *sigma32_i= (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
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
	al=al*1.0E-09;				//Computational area (m)
	b1=al/nd;					//Diff block length

	time1=0.0;					//Initial calculation count setting

	AA0=AA0e/rr/temp;			//Chemical driving force for martensite transformation

	kappa_s1=kappa_s1c/rr/temp/b1/b1;	//gradient energy factor
	kappa_s2=kappa_s2c/rr/temp/b1/b1;	//gradient energy factor

//*** Setting the Eigen distortion of the s1 field ***************
	//eta_s1[1][1]=0.05; eta_s1[2][2]=-0.05;
	//eta_s1[1][1]=eta11_s1;
	//eta_s1[2][2]=eta22_s1;
	//eta_s1[3][3]=eta33_s1;
	eta_s1[1][2]=eta_s1[2][1]=0.0;
	eta_s1[1][3]=eta_s1[3][1]=0.0;
	eta_s1[2][3]=eta_s1[3][2]=0.0;

//*** Setting the Eigen distortion of the s2 field ***************
	//eta_s2[1][1]=eta11_s2;
	//eta_s2[2][2]=eta22_s2;
	//eta_s2[3][3]=eta33_s2;
	eta_s2[1][2]=eta_s2[2][1]=0.0;
	eta_s2[1][3]=eta_s2[3][1]=0.0;
	eta_s2[2][3]=eta_s2[3][2]=0.0;

//***** elastic constant ****************************
	for(i=0;i<=3;i++){
		for(j=0;j<=3;j++){
			for(k=0;k<=3;k++){
				for(l=0;l<=3;l++){
					cec[i][j][k][l]=0.0;//Initialization of elastic constant array
				}
			}
		}
	}
	
	el_fac=1.0E+11*vm0/rr/temp;
	c11=c11c*el_fac; c22=c22c*el_fac; c33=c33c*el_fac;
	c12=c12c*el_fac; c13=c13c*el_fac; c23=c23c*el_fac;
	c44=c44c*el_fac; c55=c55c*el_fac; c66=c66c*el_fac;
	
	cec[1][1][1][1]=c11;
	cec[2][2][2][2]=c22;
	cec[3][3][3][3]=c33;
	cec[2][3][2][3]=cec[2][3][3][2]=cec[3][2][2][3]=cec[3][2][3][2]=c44;
	cec[1][3][1][3]=cec[1][3][3][1]=cec[3][1][1][3]=cec[3][1][3][1]=c55;
	cec[1][2][1][2]=cec[1][2][2][1]=cec[2][1][1][2]=cec[2][1][2][1]=c66;
	cec[1][1][2][2]=cec[2][2][1][1]=c12;
	cec[1][1][3][3]=cec[3][3][1][1]=c13;
	cec[2][2][3][3]=cec[3][3][2][2]=c23;
	
	//c12=c11-2.0*c44;
	//if(lam0==0.0){
	//	lam0=c12;//Set Lahme constants from cij data
	//}
	//if(mu0==0.0){
	//	mu0=c44;//Set Lahme constants from cij data
	//}
	//nu0=lam0/2.0/(lam0+mu0);//Poisson's ratio
	//printf("nu0= %f  \n", nu0);//Display the value of Poisson's ratio (to check the validity of the Lame constant)

//*** External force setting *******************************
 	//sig22_a=0.0;//In this calculation, set 0 because external force is not considered.
	//ep11_a=-0.5*lam0/mu0/(3.0*lam0+2.0*mu0)*sig22_a;
	//ep22_a=(lam0+mu0)/mu0/(3.0*lam0+2.0*mu0)*sig22_a;
	//ep12_a=ep21_a=0.0;

//*** Setting up sin and cos tables, bit reversal tables, and initial fields ***************

	ini000(s1h, s2h, ND);		//Initial field setting

//**** Simulation start ******************************
//Nstep = 10;
iout = -1;
 plan = fftw_plan_dft_3d(fftsize, fftsize, fftsize, in, out, FFTW_FORWARD,  FFTW_ESTIMATE);	//For FFT
iplan = fftw_plan_dft_3d(fftsize, fftsize, fftsize, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);	//For IFFT
start:;

	//if(time1<=100.){Nstep=10;} else{Nstep=200;}		//Changing the time interval for saving data
	//if((((int)(time1) % Nstep)==0)) {datsave(s1h, s2h, ND);} 	//Save tissue data every fixed repeat count
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(s1h, s2h, ND);} 	//Save tissue data every fixed repeat count

//***** gradient potential ***********************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==ndm){ip=0;}  if(i==0){im=ndm;}	//periodic boundary conditions
				if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
				if(k==ndm){kp=0;}  if(k==0){km=ndm;}
				//
				s1k_su[i*ND*ND+j*ND+k]=-kappa_s1*(
					 s1h[ip*ND*ND+j*ND+k]+s1h[im*ND*ND+j*ND+k]
					+s1h[i*ND*ND+jp*ND+k]+s1h[i*ND*ND+jm*ND+k]
					+s1h[i*ND*ND+j*ND+kp]+s1h[i*ND*ND+j*ND+km]
					-6.0*s1h[i*ND*ND+j*ND+k]);//eq.(4.2.4)
				//
				s2k_su[i*ND*ND+j*ND+k]=-kappa_s2*(
					 s2h[ip*ND*ND+j*ND+k]+s2h[im*ND*ND+j*ND+k]
					+s2h[i*ND*ND+jp*ND+k]+s2h[i*ND*ND+jm*ND+k]
					+s2h[i*ND*ND+j*ND+kp]+s2h[i*ND*ND+j*ND+km]
					-6.0*s2h[i*ND*ND+j*ND+k]);
			}
		}
	}

//**** Fourier transform of the Eigen strain field [equation (4.7)] ep11 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k]=ep11h0[i*ND*ND+j*ND+k]=eta_s1[1][1]*s1h[i*ND*ND+j*ND+k]+eta_s2[1][1]*s2h[i*ND*ND+j*ND+k];
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

//**** Fourier transform of the Eigen strain field [equation (4.7)] ep22 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				xr[i*ND*ND+j*ND+k]=ep22h0[i*ND*ND+j*ND+k]=eta_s1[2][2]*s1h[i*ND*ND+j*ND+k]+eta_s2[2][2]*s2h[i*ND*ND+j*ND+k];
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
				xr[i*ND*ND+j*ND+k]=ep33h0[i*ND*ND+j*ND+k]=eta_s1[3][3]*s1h[i*ND*ND+j*ND+k]+eta_s2[3][3]*s2h[i*ND*ND+j*ND+k];
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
				xr[i*ND*ND+j*ND+k]=ep12h0[i*ND*ND+j*ND+k]=eta_s1[1][2]*s1h[i*ND*ND+j*ND+k]+eta_s2[1][2]*s2h[i*ND*ND+j*ND+k];
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
				xr[i*ND*ND+j*ND+k]=ep13h0[i*ND*ND+j*ND+k]=eta_s1[1][3]*s1h[i*ND*ND+j*ND+k]+eta_s2[1][3]*s2h[i*ND*ND+j*ND+k];
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
				xr[i*ND*ND+j*ND+k]=ep23h0[i*ND*ND+j*ND+k]=eta_s1[2][3]*s1h[i*ND*ND+j*ND+k]+eta_s2[2][3]*s2h[i*ND*ND+j*ND+k];
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
	ep0[1][2]=sum12/nd/nd/nd; ep0[1][3]=sum13/nd/nd/nd; ep0[2][3]=sum23/nd/nd/nd;
	
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

//****** Potential calculation ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//
				s1=s1h[i*ND*ND+j*ND+k];
				s2=s2h[i*ND*ND+j*ND+k];

//****** Calculation of chemical potential [equation (4.4)] ********************************
				// Landau polynomial expansion
				// fchem = df*{(a/2)*sum(sp^2)-(b/3)*sum(sp^3)+(c/4)*(sum(sp^2))^2}
				// chem_pot = d(fchem)/d(sp)
				s1k_chem=AA0*s1*(AA1-AA2*s1+AA3*(s1*s1+s2*s2));
				s2k_chem=AA0*s2*(AA1-AA2*s2+AA3*(s1*s1+s2*s2));

//****** Calculation of elastic potential [equation (4.8)] ********************************
				// epsilonT ij = epsilon0 ij - average epsilon0 ij - d(epsilonC ij) - epsilonA ij
				epT[1][1] = ep11h0[i*ND*ND+j*ND+k] - ep0[1][1] - ec11[i*ND*ND+j*ND+k] - ep_a[1][1];
				epT[2][2] = ep22h0[i*ND*ND+j*ND+k] - ep0[2][2] - ec22[i*ND*ND+j*ND+k] - ep_a[2][2];
				epT[3][3] = ep33h0[i*ND*ND+j*ND+k] - ep0[3][3] - ec33[i*ND*ND+j*ND+k] - ep_a[3][3];
	epT[1][2] = epT[2][1] = ep12h0[i*ND*ND+j*ND+k] - ep0[1][2] - ec12[i*ND*ND+j*ND+k] - ep_a[1][2];
	epT[1][3] = epT[3][1] = ep13h0[i*ND*ND+j*ND+k] - ep0[1][3] - ec13[i*ND*ND+j*ND+k] - ep_a[1][3];
	epT[2][3] = epT[3][2] = ep23h0[i*ND*ND+j*ND+k] - ep0[2][3] - ec23[i*ND*ND+j*ND+k] - ep_a[3][3];
				//
				//d(Estr)/d(s1) = C(i,j,k,l) * [el_strain(k,l)] * eigen(i,j)_s1
				//d(Estr)/d(s1) = C(i,j,k,l) * [epT(k,l)] * eta_s1(i,j)
				s1k_str= cec[1][1][1][1]*eta_s1[1][1]*epT[1][1]
						+cec[2][2][2][2]*eta_s1[2][2]*epT[2][2]
						+cec[3][3][3][3]*eta_s1[3][3]*epT[3][3]
						+cec[1][1][2][2]*eta_s1[1][1]*epT[2][2]*2.0
						+cec[1][1][3][3]*eta_s1[1][1]*epT[3][3]*2.0
						+cec[2][2][3][3]*eta_s1[2][2]*epT[3][3]*2.0
						+cec[2][3][2][3]*eta_s1[2][3]*epT[2][3]*4.0
						+cec[1][3][1][3]*eta_s1[1][3]*epT[1][3]*4.0
						+cec[1][2][1][2]*eta_s1[1][2]*epT[1][2]*4.0;
				//
				s2k_str= cec[1][1][1][1]*eta_s2[1][1]*epT[1][1]
						+cec[2][2][2][2]*eta_s2[2][2]*epT[2][2]
						+cec[3][3][3][3]*eta_s2[3][3]*epT[3][3]
						+cec[1][1][2][2]*eta_s2[1][1]*epT[2][2]*2.0
						+cec[1][1][3][3]*eta_s2[1][1]*epT[3][3]*2.0
						+cec[2][2][3][3]*eta_s2[2][2]*epT[3][3]*2.0
						+cec[2][3][2][3]*eta_s2[2][3]*epT[2][3]*4.0
						+cec[1][3][1][3]*eta_s2[1][3]*epT[1][3]*4.0
						+cec[1][2][1][2]*eta_s2[1][2]*epT[1][2]*4.0;

//****** Calculate the time evolution of the phase field [equation (4.10)] ********************************
				s1ddtt = -smob*( s1k_chem + s1k_su[i*ND*ND+j*ND+k] + s1k_str );
				s2ddtt = -smob*( s2k_chem + s2k_su[i*ND*ND+j*ND+k] + s2k_str );
				//
				s1h[i*ND*ND+j*ND+k] = s1h[i*ND*ND+j*ND+k] + ( s1ddtt+ds_fac*(2.0*DRND(1.0)-1.0) )*delt;	//Explicit method
				s2h[i*ND*ND+j*ND+k] = s2h[i*ND*ND+j*ND+k] + ( s2ddtt+ds_fac*(2.0*DRND(1.0)-1.0) )*delt;

//*** Correction for s domain (0<=s<=1) ***
				if(s1h[i*ND*ND+j*ND+k]>=1.0){s1h[i*ND*ND+j*ND+k]=1.0;}
				if(s1h[i*ND*ND+j*ND+k]<=0.0){s1h[i*ND*ND+j*ND+k]=0.0;}
				if(s2h[i*ND*ND+j*ND+k]>=1.0){s2h[i*ND*ND+j*ND+k]=1.0;}
				if(s2h[i*ND*ND+j*ND+k]<=0.0){s2h[i*ND*ND+j*ND+k]=0.0;}
			}
		}
	}

	//if(keypress()){return 0;}//Waiting for key

	time1=time1+1.0;								//Add calculation count
	if(time1<time1max){goto start;}	//Determining if the maximum count has been reached
	
end:;
	printf("Finished \n");
	if(plan) fftw_destroy_plan(plan);	//For FFT
	if(iplan) fftw_destroy_plan(iplan);	//For IFFT
	fftw_free(in); fftw_free(out);		// For FFT and IFFT
	std::exit(0);
	//return 0;
}

//************ Initial field setup subroutine *************
void ini000(double *s1h, double *s2h, int ND)
{
	int i, j, k;
	//srand(time(NULL)); // random number initialization
	int nd=ND, ndm=ND-1, nd2=ND/2;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				s1h[i*ND*ND+j*ND+k]=0.2*DRND(1.0);//Set the field with random numbers up to 20%
				s2h[i*ND*ND+j*ND+k]=0.2*DRND(1.0);//Set the field with random numbers up to 20%
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
void datsave(double *s1h, double *s2h, int ND)
{
	FILE *stream;		//Stream pointer setting
	int i, j, k;			//integer
	int nd=ND, ndm=ND-1, nd2=ND/2;

	stream = fopen("test.dat", "a");	//Open the write destination file for appending
	fprintf(stream, "%e\n", time1);		//Saving calculation counts
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				fprintf(stream, "%e  %e ", s1h[i*ND*ND+j*ND+k], s2h[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(stream, "\n");	//writing a newline
	fclose(stream);			//close file
}

void datsave_paraview(double *s1h, double *s2h, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j, k;
	int nd=ND, ndm=ND-1, nd2=ND/2;
	
	iout = iout + 1;
	printf("pf_result%06d.vtk \n",iout);
	sprintf(fName,"pf_result%06d.vtk",iout);
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
				fprintf(fp,"%10.6f\n", s1h[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Phase_field_2 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				fprintf(fp,"%10.6f\n", s2h[i*ND*ND+j*ND+k]);
			}
		}
	}
	fclose(fp);
}
