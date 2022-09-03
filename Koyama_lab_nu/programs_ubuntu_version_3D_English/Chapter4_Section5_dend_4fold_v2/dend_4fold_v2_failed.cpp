#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//Random number function setting

//#define NDP 351							//Number of divisions per side of calculation area in difference calculation + 1

	//int nd=NDP-1, ndm=NDP-2, nd2=(NDP-1)/2;
	//Define the number of difference blocks on one side of the calculation area, define nd-1, define nd/2
	double delt;						//ticking
	double PI=3.14159;					//Pi
	double RR=8.3145;					//gas constant
	double Tini, Tm, time1;				//Initial temperature, melting point, number of calculation counts
	//double s1h[NDP][NDP][NDP], Th[NDP][NDP][NDP];	//phase field, temperature field
	int Nstep, iout;

	void ini000(double *s1h, double *Th, int NDP);	//Initial field setup subroutine
	void datsave(double *s1h, double *Th, int NDP);	//data save subroutine
	void datsave_paraview(double *s1h, double *Th, int NDP);//data save subroutine
	

//******* main program ******************************************
int main(void)
{
	int NDP;			//Number of divisions per side of calculation area in difference calculation + 1
	int nd, ndm, nd2;	//Define the number of difference blocks on one side of the calculation area, define nd-1, define nd/2
	
	//double s1h2[NDP][NDP], Th2[NDP][NDP];	//Auxiliary array of fields
	double s1kai, s1kais;			//potential
	int    i, j, k, l;				//integer
	int    ip, im, jp, jm, kp, km;	//integer
	double time1max;				//Maximum calculation count (calculation end count)

	double s1;						//phase field
	double s1ip, s1im;				//s1 left, right, top, bottom phase field
	double s1jp, s1jm;
	double s1kp, s1km;
	double s1ipjp, s1ipjm;			//Diagonal side phase field of s1
	double s1imjp, s1imjm;
	double s1ipkp, s1ipkm;
	double s1imkp, s1imkm;
	double s1jpkp, s1jpkm;
	double s1jmkp, s1jmkm;
	double s1ddtt;					//Time variation of s1 (left side of evolution equation)

	double TT;						//temperature field
	double Tip, Tim, Tjp, Tjm, Tkp, Tkm;	//Left, right, top, bottom temperature of TT
	double Tddtt;					//Temporal change in temperature (left side of thermal diffusion equation)

	double ep, ep1p, ep2p;			//Gradient energy coefficient related variables
	double dx_s1, dy_s1, dz_s1;		//Variables related to spatial differentiation of the phase field
	double l_dxdy;
	double dxx_s1, dyy_s1, dzz_s1;	//Variables related to spatial differentiation of the phase field
	double dxy_s1, dxz_s1, dyz_s1;

	double al;				//length of one side of computational domain
	double dx, dy, dz;		//Difference grid size (x direction)
	double gamma;			//interfacial energy density
	double delta;			//interface width
	double ram;				//lambda
	double bbb;				//Coefficient related to interface width
	double j_fold;			//Anisotropic mode
	double astre;			//anisotropic strength
	double cndct;			//Thermal conductivity
	double speht;			//specific heat
	double rlate;			//latent heat
	double skine;			//interface kinetic coefficient
	double th0;				//Angle of preferred growth direction
	double th;				//the angle normal to the interface
	double aaa, aaac;		//gradient energy factor 0
	double www, wwwc;		//Energy Barrier for Penalty Term
	double pmobi, pmobic;	//Phasefield Mobility
	double dtp;				//ticking
	double dtt;				//ticking
	double anois;			//Amplitude of noise
	double dF;				//driving force
	double dami1, dami2;	//Variables used for setting conditions to skip field calculations

	double zeta;			//zeta = astre;
	double div_s1;			//div(s1) = div(ph)
	double dx_ep, dy_ep, dz_ep;
	double ddphdx_ep, ddphdy_ep, ddphdz_ep;	// ddphdx_ep = d(ep)/d(d(s1)/dx)
	double dx_ddphdx_ep, dy_ddphdy_ep, dz_ddphdz_ep;	// ddphdx_ep = d(ep)/d(d(s1)/dx)
	double term1, term2x, term2y, term2z, term3;
	double dx_abs_div_ph, dy_abs_div_ph, dz_abs_div_ph;	// d|div(ph)|/dx, d|div(ph)|/dy, d|div(ph)|/dz

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
	NDP     = int(data[0]);
	delt    = data[1];
	dx      = data[2];
	dy      = data[3];
	dz      = data[4];
	gamma   = data[5];
	ram     = data[6];
	j_fold  = data[7];
	astre   = data[8];
	cndct   = data[9];
	speht   = data[10];
	rlate   = data[11];
	Tm      = data[12];
	Tini    = data[13];
	skine   = data[14];
	th0     = data[15];
	aaac    = data[16];
	wwwc    = data[17];
	pmobic  = data[18];
	anois   = data[19];
	time1max= int(data[20]);
	Nstep   = int(data[21]);
	printf("---------------------------------\n");
	//
	nd=NDP-1;
	ndm=NDP-2;
	nd2=(NDP-1)/2;
	//
	double *s1h  = (double *)malloc(sizeof(double)*( NDP*NDP*NDP + NDP*NDP + NDP ));	//phase field
	double *Th   = (double *)malloc(sizeof(double)*( NDP*NDP*NDP + NDP*NDP + NDP ));	//temperature field
	double *s1h2 = (double *)malloc(sizeof(double)*( NDP*NDP*NDP + NDP*NDP + NDP ));
	double *Th2  = (double *)malloc(sizeof(double)*( NDP*NDP*NDP + NDP*NDP + NDP ));
	//
	//printf("dtemp(1.0)=  ");	scanf(" %lf",&dtemp);	//dtemp=1.0;
	//printf("DELT(1.5)=  "); scanf(" %lf",&delt);	//delt=1.5;

	//dx=dy=30.0e-9;		//Difference grid size (x direction) [m]
	delta=3.0*dx;      		//interface width [m]
	//dx=dy=20.0e-9;		//Difference grid size (x direction) [m]
	//delta=4.0*dx;    		//interface width [m]
	al=dx*(double)nd;		//Computational domain [m]
	//gamma=0.37;       	//interfacial energy [J/m^2]
	//ram=0.1;        		//lambda
	bbb=2.0*log((1.0+(1.0-2.0*ram))/(1.0-(1.0-2.0*ram)))/2.0;//Coefficient related to interface width
	//j_fold=4.0;         	//Anisotropic mode

	//astre=0.005;       	//anisotropic strength
	//astre=0.01;       	//anisotropic strength
	//astre=0.03;       	//Anisotropic strength (Fig. 4.8)
	//astre=0.05;       	//anisotropic strength
	zeta = astre;

//    << pure Ni >>
	//cndct=84.01;      	//Thermal conductivity [W/mK]
	//speht=5.42e+06;   	//Specific heat [J/Km^3]
	//rlate=2.350e+09;  	//latent heat [J/m^3]
	//Tm=1728.0;      		//Melting point [K]
	//Tini=1511.2;     		//initial temperature [K]
	//Tini=1424.5;     		//initial temperature [K]
	//skine=2.0;        	//interface kinetic coefficient [m/Ks]
	//th0=0.0;         		//Angle of preferred growth direction

//---- See Section 4.5.2 ----------------------------------------------------
	//aaa=sqrt(3.0*delta*gamma/bbb);         	//gradient energy factor
	aaa=sqrt(aaac*delta*gamma/bbb);         	//gradient energy factor
	//www=6.0*gamma*bbb/delta;               	//Energy Barrier for Penalty Term
	www=wwwc*gamma*bbb/delta;               	//Energy Barrier for Penalty Term
	//pmobi=bbb*Tm*skine/(3.0*delta*rlate); 	//Phasefield Mobility
	pmobi=bbb*Tm*skine/(pmobic*delta*rlate); 	//Phasefield Mobility

	dtp=dx*dx/(5.0*pmobi*aaa*aaa);	//time increments
	dtt=dx*dx/(5.0*cndct/speht);	//time increments
	if(dtp>dtt){delt=dtt;} else{delt=dtp;}
	printf("delt= %e \n", delt);
//-----------------------------------------------------------------

	//anois=0.1;	//Amplitude of noise

	time1=0.0;		//Initial calculation time
	//time1max=1.0+1.0e+08;	//Maximum calculation time

//*** Initial concentration field setting and drawing window display *****************************************

	ini000(s1h, Th, NDP);//Initial field setting

//**** Simulation start ******************************
//Nstep = 100;
iout = -1;
start: ;

	if((((int)(time1) % Nstep)==0)){datsave(s1h, Th, NDP);}		//Save the field every fixed repetition count
	if((((int)(time1) % Nstep)==0)){datsave_paraview(s1h, Th, NDP);}//Save the field every fixed repetition count

//****** Time evolution of phase field and temperature field  **************
	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			for(k=0;k<=nd;k++){
				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==nd){ip=ndm;}  if(i==0){im=1;}
				if(j==nd){jp=ndm;}  if(j==0){jm=1;}
				if(k==nd){kp=ndm;}  if(k==0){km=1;}

				//s1=s1h[i][j];//phase field
				//s1ip=s1h[ip][j];  s1im=s1h[im][j];  s1jp=s1h[i][jp];  s1jm=s1h[i][jm];
				//s1ipjp=s1h[ip][jp]; s1ipjm=s1h[ip][jm];  s1imjp=s1h[im][jp];  s1imjm=s1h[im][jm];
				s1=s1h[i*NDP*NDP+j*NDP+k];//phase field
				s1ip=s1h[ip*NDP*NDP+j*NDP+k];  s1im=s1h[im*NDP*NDP+j*NDP+k];
				s1jp=s1h[i*NDP*NDP+jp*NDP+k];  s1jm=s1h[i*NDP*NDP+jm*NDP+k];
				s1kp=s1h[i*NDP*NDP+j*NDP+kp];  s1km=s1h[i*NDP*NDP+j*NDP+km];
				// i and j, x and y
				s1ipjp=s1h[ip*NDP*NDP+ip*NDP+k]; s1ipjm=s1h[ip*NDP*NDP+im*NDP+k];
				s1imjp=s1h[im*NDP*NDP+jp*NDP+k]; s1imjm=s1h[im*NDP*NDP+jm*NDP+k];
				// i and k, x and z
				s1ipkp=s1h[ip*NDP*NDP+i*NDP+kp]; s1ipkm=s1h[ip*NDP*NDP+i*NDP+km];
				s1imkp=s1h[im*NDP*NDP+j*NDP+kp]; s1imkm=s1h[im*NDP*NDP+j*NDP+km];
				// j and k, y and z
				s1jpkp=s1h[i*NDP*NDP+ip*NDP+kp]; s1jpkm=s1h[i*NDP*NDP+ip*NDP+km];
				s1jmkp=s1h[i*NDP*NDP+jm*NDP+kp]; s1jmkm=s1h[i*NDP*NDP+jm*NDP+km];

				//TT=Th[i][j][k];//temperature field
				//Tip=Th[ip][j][k];  Tim=Th[im][j][k];
				//Tjp=Th[i][jp][k];  Tjm=Th[i][jm][k];
				//Tkp=Th[i][j][kp];  Tkm=Th[i][j][km];
				TT=Th[i*NDP*NDP+j*NDP+k];//temperature field
				Tip=Th[ip*NDP*NDP+j*NDP+k];  Tim=Th[im*NDP*NDP+j*NDP+k];
				Tjp=Th[i*NDP*NDP+jp*NDP+k];  Tjm=Th[i*NDP*NDP+jm*NDP+k];
				Tkp=Th[i*NDP*NDP+j*NDP+kp];  Tkm=Th[i*NDP*NDP+j*NDP+km];

//----- Deciding When to Skip Calculations ----------------------------------------------
				dami1=fabs(s1+s1ip+s1im+s1jp+s1jm+s1kp+s1km);
				dami2=fabs(TT+Tip+Tim+Tjp+Tjm+Tkp+Tkm-7.0*Tini);
				if( (dami1<=1.0e-20)&&(dami2<=1.0e-20) ){
					//s1h2[i][j][K]=s1h[i][j][K];  Th2[i][j][K]=Th[i][j][K];
					s1h2[i*NDP*NDP+j*NDP+k]=s1h[i*NDP*NDP+j*NDP+k];
					 Th2[i*NDP*NDP+j*NDP+k]= Th[i*NDP*NDP+j*NDP+k];
					goto dami;
				}
//---------------------------------------------------------------------------------

				dx_s1=(s1ip-s1im)/(2.0*dx);  			//Spatial first derivative of the phase field
				dy_s1=(s1jp-s1jm)/(2.0*dy);
				dz_s1=(s1kp-s1km)/(2.0*dz);
				dxx_s1=(s1ip+s1im-2.0*s1)/(dx*dx);		//Spatial second derivative of the phase field
				dyy_s1=(s1jp+s1jm-2.0*s1)/(dy*dy);
				dzz_s1=(s1kp+s1km-2.0*s1)/(dz*dz);
				//d(dx)/dy=( (s1ip-s1im)jp - (s1ip-s1im)jm )/(4.0*dx*dy)
				//		=(s1ipjp-s1imjp-s1ipjm+s1imjm)/(4.0*dx*dy)
				//		=(s1ipjp+s1imjm-s1imjp-s1ipjm)/(4.0*dx*dy)
				dxy_s1=(s1ipjp+s1imjm-s1imjp-s1ipjm)/(4.0*dx*dy);
				//d(dx)/dz=( (s1ip-s1im)kp - (s1ip-s1im)km )/(4.0*dx*dy)
				//		=(s1ipkp-s1imkp-s1ipkm+s1imkm)/(4.0*dx*dy)
				//		=(s1ipkp+s1imkm-s1imkp-s1ipkm)/(4.0*dx*dy)
				dxz_s1=(s1ipkp+s1imkm-s1imkp-s1ipkm)/(4.0*dx*dz);
				//d(dy)/dz=( (s1jp-s1jm)kp - (s1jp-s1jm)km )/(4.0*dx*dy)
				//		=(s1jpkp-s1jmkp-s1jpkm+s1jmkm)/(4.0*dx*dy)
				//		=(s1jpkp+s1jmkm-s1jmkp-s1jpkm)/(4.0*dx*dy)
				dyz_s1=(s1jpkp+s1jmkm-s1jmkp-s1jpkm)/(4.0*dy*dz);
				//
				//Eq.(4.24): theta = arctan((dph/dy)/(dph/dx))
				//th=atan(dy_s1/(dx_s1+1.0e-20));			//angle normal to interface [equation (4.24)]
				
				//Eq.(4.23): epsion(theta)=epsilon0*(1-zeta*cos(k*(theta-theta0)))
				//ep=aaa*(1.0+astre*cos(j_fold*(th-th0)));	//square root of gradient energy coefficient [equation (4.23)]
				//ep1p=-aaa*astre*j_fold*sin(j_fold*(th-th0));	//first derivative of ep with respect to angle
				//ep2p=-aaa*astre*j_fold*j_fold*cos(j_fold*(th-th0));	//Second derivative of ep with respect to angle
				
				//other test: Eq.(4.18): a = ep
				//zeta = astre;
				div_s1 = sqrt(dx_s1*dx_s1 + dy_s1*dy_s1 + dz_s1*dz_s1);
				// Processing to skip to avoid occurrence of "nan"
				if( div_s1<=1.0e-20 ){
					s1h2[i*NDP*NDP+j*NDP+k]=s1h[i*NDP*NDP+j*NDP+k];
					 Th2[i*NDP*NDP+j*NDP+k]= Th[i*NDP*NDP+j*NDP+k];
					goto dami;
				}
				ep = aaa*(1.0-3.0*zeta)*(
					1.0 + (4.0*zeta)/(1.0-3.0*zeta)*(pow(dx_s1,4.0)+pow(dy_s1,4.0)+pow(dz_s1,4.0))/(pow(div_s1,4.0))
					//1.0 + (4.0*zeta)/(1.0-3.0*zeta)*(pow(dx_s1,4.0)+pow(dy_s1,4.0)+pow(dz_s1,4.0))*(pow(div_s1,-4.0))
				);
				//
				//dx_abs_div_ph = (dx_s1*dx_s1 + dy_s1*dy_s1 + dz_s1*dz_s1)*(dx_s1*dxx_s1 + dy_s1*dxy_s1 + dz_s1*dxz_s1);
				//dy_abs_div_ph = (dx_s1*dx_s1 + dy_s1*dy_s1 + dz_s1*dz_s1)*(dx_s1*dxy_s1 + dy_s1*dyy_s1 + dz_s1*dyz_s1);
				//dz_abs_div_ph = (dx_s1*dx_s1 + dy_s1*dy_s1 + dz_s1*dz_s1)*(dx_s1*dxz_s1 + dy_s1*dyz_s1 + dz_s1*dzz_s1);
				dx_abs_div_ph = div_s1*div_s1*(dx_s1*dxx_s1 + dy_s1*dxy_s1 + dz_s1*dxz_s1);
				dy_abs_div_ph = div_s1*div_s1*(dx_s1*dxy_s1 + dy_s1*dyy_s1 + dz_s1*dyz_s1);
				dz_abs_div_ph = div_s1*div_s1*(dx_s1*dxz_s1 + dy_s1*dyz_s1 + dz_s1*dzz_s1);
				//
				dx_ep = aaa*4.0*zeta*(
					 (4.0*pow(dx_s1,3.0)*dxx_s1+4.0*pow(dy_s1,3.0)*dxy_s1+4.0*pow(dz_s1,3.0)*dxz_s1)*pow(div_s1,-4.0)
					//+(pow(dx_s1,4.0)+pow(dy_s1,4.0)+pow(dz_s1,4.0))*-4.0*pow(div_s1,-5.0)*dx_abs_div_ph
					+(pow(dx_s1,4.0)+pow(dy_s1,4.0)+pow(dz_s1,4.0))*
						-4.0*pow(div_s1,-3.0)*(dx_s1*dxx_s1 + dy_s1*dxy_s1 + dz_s1*dxz_s1)
				);
				dy_ep = aaa*4.0*zeta*(
					 (4.0*pow(dx_s1,3.0)*dxy_s1+4.0*pow(dy_s1,3.0)*dyy_s1+4.0*pow(dz_s1,3.0)*dyz_s1)*pow(div_s1,-4.0)
					//+(pow(dx_s1,4.0)+pow(dy_s1,4.0)+pow(dz_s1,4.0))*-4.0*pow(div_s1,-5.0)*dy_abs_div_ph
					+(pow(dx_s1,4.0)+pow(dy_s1,4.0)+pow(dz_s1,4.0))*
						-4.0*pow(div_s1,-3.0)*(dx_s1*dxy_s1 + dy_s1*dyy_s1 + dz_s1*dyz_s1)
				);
				dz_ep = aaa*4.0*zeta*(
					 (4.0*pow(dx_s1,3.0)*dxz_s1+4.0*pow(dy_s1,3.0)*dyz_s1+4.0*pow(dz_s1,3.0)*dzz_s1)*pow(div_s1,-4.0)
					//+(pow(dx_s1,4.0)+pow(dy_s1,4.0)+pow(dz_s1,4.0))*-4.0*pow(div_s1,-5.0)*dz_abs_div_ph
					+(pow(dx_s1,4.0)+pow(dy_s1,4.0)+pow(dz_s1,4.0))*
						-4.0*pow(div_s1,-3.0)*(dx_s1*dxz_s1 + dy_s1*dyz_s1 + dz_s1*dzz_s1)
				);
				term1 = 2.0*ep*(dx_ep*dx_s1 + dy_ep*dy_s1 + dz_ep*dz_s1) 
					  + ep*ep*(dxx_s1 + dyy_s1 + dzz_s1);
				
				// ddphdx_ep = d(ep)/d(d(s1)/dx)
				ddphdx_ep = aaa*4.0*zeta*(
							4.0*pow(dx_s1,3.0)*pow(div_s1,-4.0)
							+(pow(dx_s1,4.0)+pow(dy_s1,4.0)+pow(dz_s1,4.0))*
								//-4.0*pow(div_s1,-5.0)*0.5*(div_s1*div_s1)*2.0*dx_s1
								-4.0*pow(div_s1,-3.0)*dx_s1
				);
				ddphdy_ep = aaa*4.0*zeta*(
							4.0*pow(dy_s1,3.0)*pow(div_s1,-4.0)
							+(pow(dx_s1,4.0)+pow(dy_s1,4.0)+pow(dz_s1,4.0))*
								//-4.0*pow(div_s1,-5.0)*0.5*(div_s1*div_s1)*2.0*dy_s1
								-4.0*pow(div_s1,-3.0)*dy_s1
				);
				ddphdz_ep = aaa*4.0*zeta*(
							4.0*pow(dz_s1,3.0)*pow(div_s1,-4.0)
							+(pow(dx_s1,4.0)+pow(dy_s1,4.0)+pow(dz_s1,4.0))*
								//-4.0*pow(div_s1,-5.0)*0.5*(div_s1*div_s1)*2.0*dz_s1
								-4.0*pow(div_s1,-3.0)*dz_s1
				);
				//
				dx_ddphdx_ep = aaa*4.0*zeta*(
					12.0*pow(dx_s1,2.0)*dxx_s1*pow(div_s1,-4.0)
					//+4.0*pow(dx_s1,3.0)*-4.0*pow(div_s1,-5.0)*dx_abs_div_ph
					+4.0*pow(dx_s1,3.0)*-4.0*pow(div_s1,-3.0)*(dx_s1*dxx_s1 + dy_s1*dxy_s1 + dz_s1*dxz_s1)
					//
					+(4.0*pow(dx_s1,3.0)*dxx_s1+4.0*pow(dy_s1,3.0)*dxy_s1+4.0*pow(dz_s1,3.0)*dxz_s1)*
						-4.0*pow(div_s1,-3.0)*dx_s1
					+(pow(dx_s1,4.0)+pow(dy_s1,4.0)+pow(dz_s1,4.0))*(
						//12.0*pow(div_s1,-4.0)*dx_abs_div_ph*dx_s1
						12.0*pow(div_s1,-2.0)*(dx_s1*dxx_s1 + dy_s1*dxy_s1 + dz_s1*dxz_s1)*dx_s1
						-4.0*pow(div_s1,-3.0)*dxx_s1
						)
				);
				dy_ddphdy_ep = aaa*4.0*zeta*(
					12.0*pow(dy_s1,2.0)*dyy_s1*pow(div_s1,-4.0)
					//+4.0*pow(dy_s1,3.0)*-4.0*pow(div_s1,-5.0)*dy_abs_div_ph
					+4.0*pow(dy_s1,3.0)*-4.0*pow(div_s1,-3.0)*(dx_s1*dxy_s1 + dy_s1*dyy_s1 + dz_s1*dyz_s1)
					//
					+(4.0*pow(dx_s1,3.0)*dxy_s1+4.0*pow(dy_s1,3.0)*dyy_s1+4.0*pow(dz_s1,3.0)*dyz_s1)*
						-4.0*pow(div_s1,-3.0)*dy_s1
					+(pow(dx_s1,4.0)+pow(dy_s1,4.0)+pow(dz_s1,4.0))*(
						//12.0*pow(div_s1,-4.0)*dy_abs_div_ph*dy_s1
						12.0*pow(div_s1,-2.0)*(dx_s1*dxy_s1 + dy_s1*dyy_s1 + dz_s1*dyz_s1)*dy_s1
						-4.0*pow(div_s1,-3.0)*dyy_s1
						)
				);
				dz_ddphdz_ep = aaa*4.0*zeta*(
					12.0*pow(dz_s1,2.0)*dzz_s1*pow(div_s1,-4.0)
					//+4.0*pow(dz_s1,3.0)*-4.0*pow(div_s1,-5.0)*dz_abs_div_ph
					+4.0*pow(dz_s1,3.0)*-4.0*pow(div_s1,-3.0)*(dx_s1*dxz_s1 + dy_s1*dyz_s1 + dz_s1*dzz_s1)
					//
					+(4.0*pow(dx_s1,3.0)*dxz_s1+4.0*pow(dy_s1,3.0)*dyz_s1+4.0*pow(dz_s1,3.0)*dzz_s1)*
						-4.0*pow(div_s1,-3.0)*dz_s1
					+(pow(dx_s1,4.0)+pow(dy_s1,4.0)+pow(dz_s1,4.0))*(
						//12.0*pow(div_s1,-4.0)*dz_abs_div_ph*dz_s1
						12.0*pow(div_s1,-2.0)*(dx_s1*dxz_s1 + dy_s1*dyz_s1 + dz_s1*dzz_s1)*dz_s1
						-4.0*pow(div_s1,-3.0)*dzz_s1
						)
				);
				//
				term2x = dx_ep*ddphdx_ep*(div_s1*div_s1)
					   + ep*dx_ddphdx_ep*(div_s1*div_s1)
					   + ep*ddphdx_ep*2.0*div_s1*dx_abs_div_ph;
				term2y = dy_ep*ddphdy_ep*(div_s1*div_s1)
					   + ep*dy_ddphdy_ep*(div_s1*div_s1)
					   + ep*ddphdy_ep*2.0*div_s1*dy_abs_div_ph;
				term2z = dz_ep*ddphdz_ep*(div_s1*div_s1)
					   + ep*dz_ddphdz_ep*(div_s1*div_s1)
					   + ep*ddphdz_ep*2.0*div_s1*dz_abs_div_ph;

				//gradient potential [part of equation (4.29)]
				//s1kais=-ep*ep*(dxx_s1+dyy_s1)
				//	   -ep*ep1p*((dyy_s1-dxx_s1)*sin(2.0*th)+2.0*dxy_s1*cos(2.0*th))
				//	   +0.5*(ep1p*ep1p+ep*ep2p)*(2.0*dxy_s1*sin(2.0*th)
				//	         -dxx_s1-dyy_s1-(dyy_s1-dxx_s1)*cos(2.0*th));

				//chemical driving force [part of formula (4.29)]
				// 15*L*(T-Tm)/(2*Wsl*Tm)*ph*(1-ph)
				dF=15.0/(2.0*www)*rlate*(TT-Tm)/Tm*s1*(1.0-s1);
				
				//chemical potential [part of equation (4.29)] (other test: Eq.(4.21))
				//s1kai=4.0*www*s1*(1.0-s1)*(0.5-s1+dF+anois*(DRND(1)-0.5));
				//other test: Eq.(4.21)
				term3=4.0*www*s1*(1.0-s1)*(s1-0.5-dF+anois*(DRND(1)-0.5));

				//Phase-field evolution equation [equation (4.25)]
				//s1ddtt=-pmobi*(s1kai+s1kais);
				//other test: Eq.(4.21)
				s1ddtt=pmobi*(term1 + term2x + term2y + term2z + term3);

				//s1h2[i][j][k]=s1+s1ddtt*delt;	//Phase field time evolution (explicit method)
				//s1h2[i][j][k]=s1+s1ddtt*delt+anois*(DRND(1)-0.5)*s1*(1.0-s1);
				//if(s1h2[i][j][k]>=1.0){s1h2[i][j][k]=1.0;}//Phase field domain correction
				//if(s1h2[i][j][k]<=0.0){s1h2[i][j][k]=0.0;}
				s1h2[i*NDP*NDP+j*NDP+k]=s1+s1ddtt*delt;	//Phase field time evolution (explicit method)
				if(s1h2[i*NDP*NDP+j*NDP+k]>=1.0){s1h2[i*NDP*NDP+j*NDP+k]=1.0;}//Phase field domain correction
				if(s1h2[i*NDP*NDP+j*NDP+k]<=0.0){s1h2[i*NDP*NDP+j*NDP+k]=0.0;}

				//Thermal diffusion equation [equation (4.30)]
				Tddtt=( cndct*( 
						(Tip+Tim-2.0*TT)/dx/dx
					   +(Tjp+Tjm-2.0*TT)/dy/dy
					   +(Tkp+Tkm-2.0*TT)/dz/dz
					  ) +30.0*s1*(1.0-s1)*s1*(1.0-s1)*rlate*s1ddtt )/speht;

				//Th2[i][j][k]=Th[i][j][k]+Tddtt*delt;		//Time evolution of temperature field (explicit method)
				Th2[i*NDP*NDP+j*NDP+k]=Th[i*NDP*NDP+j*NDP+k]+Tddtt*delt;	//Time evolution of temperature field (explicit method)

				dami:;
			}
		}
	}

//**** move field auxiliary array to main array ************************************
	for(i=0;i<=nd;i++){
	    //for(j=0;j<=nd;j++){  s1h[i][j][k]=s1h2[i][j][k]; Th[i][j][k]=Th2[i][j][k]; }
		for(j=0;j<=nd;j++){
			for(k=0;k<=nd;k++){
				s1h[i*NDP*NDP+j*NDP+k]=s1h2[i*NDP*NDP+j*NDP+k];
				Th[i*NDP*NDP+j*NDP+k]=Th2[i*NDP*NDP+j*NDP+k];
			}
		}
	}
//*********************************************************************

	//if(keypress()){return 0;}//Waiting for key
	time1=time1+1.;
	if(time1<time1max){goto start;}//Time Increment and Determining Whether Maximum Time Has Been Reached
	printf("Finished \n");

	end:;
	std::exit(0);
	//return 0;

}

//************ Initial field setup subroutine *************
void ini000(double *s1h, double *Th, int NDP)
{
	int i, j, k;
 	srand(time(NULL)); // random number initialization
	int nd=NDP-1, ndm=NDP-2, nd2=(NDP-1)/2;

	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			for(k=0;k<=nd;k++){
				//s1h[i][j]=0.0;
				//if((i<10.)&&(j<10.)){s1h[i][j][k]=0.9;}
				//if((i*i+j*j)<20.){s1h[i][j][k]=0.9;}
				s1h[i*NDP*NDP+j*NDP+k]=0.0;
				if((i*i+j*j+k*k)<20.0){
					s1h[i*NDP*NDP+j*NDP+k]=0.9;
				}
			}
		}
	}

	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			for(k=0;k<=nd;k++){
				//Th[i][j][k]=Tini+s1h[i][j][k]*(Tm-Tini);
				Th[i*NDP*NDP+j*NDP+k]=Tini+s1h[i*NDP*NDP+j*NDP+k]*(Tm-Tini);
			}
		}
	}

}

//************ data save subroutine *******************************
void datsave(double *s1h, double *Th, int NDP)
{
	FILE		*stream;	//Stream pointer setting
	int 		i, j, k;			//integer
	int nd=NDP-1, ndm=NDP-2, nd2=(NDP-1)/2;

	stream = fopen("test.dat", "a");	//Open the file to write to in append mode
	fprintf(stream, "%f\n", time1);		//Save repeat count
	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			for(k=0;k<=nd;k++){
				//fprintf(stream, "%e  %e  ", s1h[i][j][k], Th[i][j][k]);//Conservation of local fields
				fprintf(stream, "%e  %e  ", s1h[i*NDP*NDP+j*NDP+k], Th[i*NDP*NDP+j*NDP+k]);//Conservation of local fields
			}
		}
	}
	fprintf(stream, "\n");	//writing a newline
	fclose(stream);					//close file
}

void datsave_paraview(double *s1h, double *Th, int NDP)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j, k;
	int nd=NDP-1, ndm=NDP-2, nd2=(NDP-1)/2;
	
	iout = iout + 1;
	printf("pf_result%06d.vtk \n",iout);
	sprintf(fName,"pf_result%06d.vtk",iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(nd+1),(nd+1),(nd+1));
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %16d \n",((nd+1)*(nd+1)*(nd+1)));
	fprintf(fp,"SCALARS Phase_field float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=nd;k++){
		for(j=0;j<=nd;j++){
			for(i=0;i<=nd;i++){
				fprintf(fp,"%10.6f\n", s1h[i*NDP*NDP+j*NDP+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Temperature float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=nd;k++){
		for(j=0;j<=nd;j++){
			for(i=0;i<=nd;i++){
				fprintf(fp,"%10.6f\n", Th[i*NDP*NDP+j*NDP+k]);
			}
		}
	}
	fclose(fp);
}
