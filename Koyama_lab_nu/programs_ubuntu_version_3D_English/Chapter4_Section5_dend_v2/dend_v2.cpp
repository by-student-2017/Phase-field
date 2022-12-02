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

	void ini000(double *s1h, double *Th, int NDPX, int NDPY, int NDPZ);	//Initial field setup subroutine
	void datsave(double *s1h, double *Th, int NDPX, int NDPY, int NDPZ);	//data save subroutine
	void datsave_paraview(double *s1h, double *Th, int NDPX, int NDPY, int NDPZ);//data save subroutine
	

//******* main program ******************************************
int main(void)
{
	int NDPX, NDPY, NDPZ;	//Number of divisions per side of calculation area in difference calculation + 1
	int ndx, ndy, ndz;
	int ndxm, ndym, ndzm;
	int ndx2, ndy2, ndz2;	//Define the number of difference blocks on one side of the calculation area, define nd-1, define nd/2
	
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
	//
	double abs_dph, abs_dph_p2;
	double dx_s1_p4, dy_s1_p4, dz_s1_p4;
	double alphap, alpha2;
	double partx, party, partz;
	double term1x, term1y, term1z, term2;
	//
	double interd;

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
	NDPX    = int(data[0]);
	NDPY    = int(data[1]);
	NDPZ    = int(data[2]);
	delt    = data[3];
	dx      = data[4];
	dy      = data[5];
	dz      = data[6];
	gamma   = data[7];
	ram     = data[8];
	j_fold  = data[9];
	astre   = data[10];
	cndct   = data[11];
	speht   = data[12];
	rlate   = data[13];
	Tm      = data[14];
	Tini    = data[15];
	skine   = data[16];
	th0     = data[17];
	aaac    = data[18];
	wwwc    = data[19];
	pmobic  = data[20];
	anois   = data[21];
	interd  = data[22];
	time1max= int(data[23]);
	Nstep   = int(data[24]);
	printf("---------------------------------\n");
	//
	ndx=NDPX-1;
	ndxm=NDPX-2;
	ndx2=(NDPX-1)/2;
	//
	ndy=NDPY-1;
	ndym=NDPY-2;
	ndy2=(NDPY-1)/2;
	//
	ndz=NDPZ-1;
	ndzm=NDPZ-2;
	ndz2=(NDPZ-1)/2;
	//
	double *s1h  = (double *)malloc(sizeof(double)*( NDPX*NDPY*NDPZ + NDPY*NDPZ + NDPZ ));	//phase field
	double *Th   = (double *)malloc(sizeof(double)*( NDPX*NDPY*NDPZ + NDPY*NDPZ + NDPZ ));	//temperature field
	double *s1h2 = (double *)malloc(sizeof(double)*( NDPX*NDPY*NDPZ + NDPY*NDPZ + NDPZ ));
	double *Th2  = (double *)malloc(sizeof(double)*( NDPX*NDPY*NDPZ + NDPY*NDPZ + NDPZ ));
	//
	double *s1x  = (double *)malloc(sizeof(double)*( NDPX*NDPY*NDPZ + NDPY*NDPZ + NDPZ ));
	double *s1y  = (double *)malloc(sizeof(double)*( NDPX*NDPY*NDPZ + NDPY*NDPZ + NDPZ ));
	double *s1z  = (double *)malloc(sizeof(double)*( NDPX*NDPY*NDPZ + NDPY*NDPZ + NDPZ ));
	//
	//printf("dtemp(1.0)=  ");	scanf(" %lf",&dtemp);	//dtemp=1.0;
	//printf("DELT(1.5)=  "); scanf(" %lf",&delt);	//delt=1.5;

	//dx=dy=30.0e-9;		//Difference grid size (x direction) [m]
	delta=interd*dx;      		//interface width [m]
	//dx=dy=20.0e-9;		//Difference grid size (x direction) [m]
	//delta=4.0*dx;    		//interface width [m]
	//al=dx*(double)nd;		//Computational domain [m]
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
	//pmobi=bbb*Tm*skine/(pmobic*delta*rlate); 	//Phasefield Mobility
	pmobi=sqrt(2.0*www)/(6.0*aaa)*(Tm/rlate*skine);

	//dtp=dx*dx/(n*pmobi*aaa*aaa);	//time increments, 1d:n=1, 2d:n=4, 3d:n=6
	//dtt=dx*dx/(n*cndct/speht);	//time increments, 1d:n=1, 2d:n=4, 3d:n=6
	dtp=dx*dx/(7.0*pmobi*aaa*aaa);	//time increments, aniso case: larger than n
	dtt=dx*dx/(7.0*cndct/speht);	//time increments
	if(dtp>dtt){delt=dtt;} else{delt=dtp;}
	printf("delt= %e \n", delt);
//-----------------------------------------------------------------

	//anois=0.1;	//Amplitude of noise

	time1=0.0;		//Initial calculation time
	//time1max=1.0+1.0e+08;	//Maximum calculation time

//*** Initial concentration field setting and drawing window display *****************************************

	ini000(s1h, Th, NDPX, NDPY, NDPZ);//Initial field setting

//**** Simulation start ******************************
//Nstep = 100;
iout = -1;
start: ;

	if((((int)(time1) % Nstep)==0)){datsave(s1h, Th, NDPX, NDPY, NDPZ);}		//Save the field every fixed repetition count
	if((((int)(time1) % Nstep)==0)){datsave_paraview(s1h, Th, NDPX, NDPY, NDPZ);}//Save the field every fixed repetition count

//****** Time evolution of phase field and temperature field  **************
	for(i=0;i<=ndx;i++){
		for(j=0;j<=ndy;j++){
			for(k=0;k<=ndz;k++){
				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==ndx){ip=ndxm;}  if(i==0){im=1;}
				if(j==ndy){jp=ndym;}  if(j==0){jm=1;}
				if(k==ndz){kp=ndzm;}  if(k==0){km=1;}
				//
				s1=s1h[i*NDPY*NDPZ+j*NDPZ+k];//phase field
				s1ip=s1h[ip*NDPY*NDPZ+j*NDPZ+k];  s1im=s1h[im*NDPY*NDPZ+j*NDPZ+k];
				s1jp=s1h[i*NDPY*NDPZ+jp*NDPZ+k];  s1jm=s1h[i*NDPY*NDPZ+jm*NDPZ+k];
				s1kp=s1h[i*NDPY*NDPZ+j*NDPZ+kp];  s1km=s1h[i*NDPY*NDPZ+j*NDPZ+km];
				//
				dx_s1=(s1ip-s1im)/(2.0*dx);	//Spatial first derivative of the phase field
				dy_s1=(s1jp-s1jm)/(2.0*dy);
				dz_s1=(s1kp-s1km)/(2.0*dz);

//----- Deciding When to Skip Calculations ----------------------------------------------
				dami1=fabs(s1+s1ip+s1im+s1jp+s1jm+s1kp+s1km);
				if( dami1<=1.0e-20 ){
					//s1h2[i][j][K]=s1h[i][j][K];  Th2[i][j][K]=Th[i][j][K];
					s1h2[i*NDPY*NDPZ+j*NDPZ+k]=s1h[i*NDPY*NDPZ+j*NDPZ+k];
					goto damif;
				}
				
				abs_dph=( dx_s1*dx_s1 + dy_s1*dy_s1 + dz_s1*dz_s1 );
				if( abs_dph==0.0 ){
					s1x[i*NDPY*NDPZ+j*NDPZ+k] = 0.0;
					s1y[i*NDPY*NDPZ+j*NDPZ+k] = 0.0;
					s1z[i*NDPY*NDPZ+j*NDPZ+k] = 0.0;
					goto damif;
				}
				//
				abs_dph_p2=abs_dph*abs_dph;
				//
				dx_s1_p4=dx_s1*dx_s1*dx_s1*dx_s1; //pow(dx_s1,4.0)
				dy_s1_p4=dy_s1*dy_s1*dy_s1*dy_s1; //pow(dy_s1,4.0)
				dz_s1_p4=dz_s1*dz_s1*dz_s1*dz_s1; //pow(dz_s1,4.0)
				//
				alphap=aaa*(1.0-3.0*zeta)*(1.0 + (4.0*zeta)/(1.0-3.0*zeta)*(dx_s1_p4 + dy_s1_p4 + dz_s1_p4)/abs_dph_p2); //Eq.(4.18)
				alpha2=alphap*alphap;
				//
				partx=alphap*aaa*(4.0*zeta)*(
					//+(4.0*dx_s1*dx_s1*dx_s1)/abs_dph_p2
					//+(dx_s1_p4 + dy_s1_p4 + dz_s1_p4)*-2.0/(abs_dph*abs_dph*abs_dph)*2.0*dx_s1
					//)*abs_dph;
					//
					4.0*dx_s1*( (dx_s1*dx_s1)*(dy_s1*dy_s1 + dz_s1*dz_s1) -(dy_s1_p4 + dz_s1_p4) )/abs_dph_p2);
					//
				party=alphap*aaa*(4.0*zeta)*(
					//+(4.0*dy_s1*dy_s1*dy_s1)/abs_dph_p2
					//+(dx_s1_p4 + dy_s1_p4 + dz_s1_p4)*-2.0/(abs_dph*abs_dph*abs_dph)*2.0*dy_s1
					//)*abs_dph;
					//
					4.0*dy_s1*( (dy_s1*dy_s1)*(dx_s1*dx_s1 + dz_s1*dz_s1) -(dx_s1_p4 + dz_s1_p4) )/abs_dph_p2);
					//
				partz=alphap*aaa*(4.0*zeta)*(
					//+(4.0*dz_s1*dz_s1*dz_s1)/abs_dph_p2
					//+(dx_s1_p4 + dy_s1_p4 + dz_s1_p4)*-2.0/(abs_dph*abs_dph*abs_dph)*2.0*dz_s1
					//)*abs_dph;
					//
					4.0*dz_s1*( (dz_s1*dz_s1)*(dx_s1*dx_s1 + dy_s1*dy_s1) -(dx_s1_p4 + dy_s1_p4) )/abs_dph_p2);
					//
				//
				//parts of Eq.(4.21) before differentiation
				s1x[i*NDPY*NDPZ+j*NDPZ+k] = alpha2*dx_s1 + partx;
				s1y[i*NDPY*NDPZ+j*NDPZ+k] = alpha2*dy_s1 + party;
				s1z[i*NDPY*NDPZ+j*NDPZ+k] = alpha2*dz_s1 + partz;
				//
				damif:;
			}
		}
	}
	//
	for(i=0;i<=ndx;i++){
		for(j=0;j<=ndy;j++){
			for(k=0;k<=ndz;k++){
				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==ndx){ip=ndxm;}  if(i==0){im=1;}
				if(j==ndy){jp=ndym;}  if(j==0){jm=1;}
				if(k==ndz){kp=ndzm;}  if(k==0){km=1;}
				
				s1=s1h[i*NDPY*NDPZ+j*NDPZ+k];//phase field
				s1ip=s1h[ip*NDPY*NDPZ+j*NDPZ+k];  s1im=s1h[im*NDPY*NDPZ+j*NDPZ+k];
				s1jp=s1h[i*NDPY*NDPZ+jp*NDPZ+k];  s1jm=s1h[i*NDPY*NDPZ+jm*NDPZ+k];
				s1kp=s1h[i*NDPY*NDPZ+j*NDPZ+kp];  s1km=s1h[i*NDPY*NDPZ+j*NDPZ+km];

				//TT=Th[i][j][k];//temperature field
				//Tip=Th[ip][j][k];  Tim=Th[im][j][k];
				//Tjp=Th[i][jp][k];  Tjm=Th[i][jm][k];
				//Tkp=Th[i][j][kp];  Tkm=Th[i][j][km];
				TT=Th[i*NDPY*NDPZ+j*NDPZ+k];//temperature field
				Tip=Th[ip*NDPY*NDPZ+j*NDPZ+k];  Tim=Th[im*NDPY*NDPZ+j*NDPZ+k];
				Tjp=Th[i*NDPY*NDPZ+jp*NDPZ+k];  Tjm=Th[i*NDPY*NDPZ+jm*NDPZ+k];
				Tkp=Th[i*NDPY*NDPZ+j*NDPZ+kp];  Tkm=Th[i*NDPY*NDPZ+j*NDPZ+km];

//----- Deciding When to Skip Calculations ----------------------------------------------
				dami1=fabs(s1+s1ip+s1im+s1jp+s1jm+s1kp+s1km);
				dami2=fabs(TT+Tip+Tim+Tjp+Tjm+Tkp+Tkm-7.0*Tini);
				if( (dami1<=1.0e-20) && (dami2<=1.0e-20) ){
					//s1h2[i][j][K]=s1h[i][j][K];  Th2[i][j][K]=Th[i][j][K];
					s1h2[i*NDPY*NDPZ+j*NDPZ+k]=s1h[i*NDPY*NDPZ+j*NDPZ+k];
					 Th2[i*NDPY*NDPZ+j*NDPZ+k]= Th[i*NDPY*NDPZ+j*NDPZ+k];
					goto dami;
				}
//---------------------------------------------------------------------------------
				//
				//Eq.(4.24): theta = arctan((dph/dy)/(dph/dx))
				//th=atan(dy_s1/(dx_s1+1.0e-20));			//angle normal to interface [equation (4.24)]
				
				//Eq.(4.23): epsion(theta)=epsilon0*(1-zeta*cos(k*(theta-theta0)))
				//ep=aaa*(1.0+astre*cos(j_fold*(th-th0)));	//square root of gradient energy coefficient [equation (4.23)]
				//ep1p=-aaa*astre*j_fold*sin(j_fold*(th-th0));	//first derivative of ep with respect to angle
				//ep2p=-aaa*astre*j_fold*j_fold*cos(j_fold*(th-th0));	//Second derivative of ep with respect to angle
				
				//other test: Eq.(4.18): a = ep
				//zeta = astre;
				//div_s1 = sqrt(dx_s1*dx_s1 + dy_s1*dy_s1 + dz_s1*dz_s1);
				// Processing to skip to avoid occurrence of "nan"
				//if( div_s1<=1.0e-20 ){
				//	s1h2[i*NDP*NDP+j*NDP+k]=s1h[i*NDP*NDP+j*NDP+k];
				//	 Th2[i*NDP*NDP+j*NDP+k]= Th[i*NDP*NDP+j*NDP+k];
				//	goto dami;
				//}
				
				//ep = alphap
				//ep = aaa*(1.0-3.0*zeta)*(
				//	1.0 + (4.0*zeta)/(1.0-3.0*zeta)*(pow(dx_s1,4.0)+pow(dy_s1,4.0)+pow(dz_s1,4.0))/(pow(div_s1,4.0))
				//	//1.0 + (4.0*zeta)/(1.0-3.0*zeta)*(pow(dx_s1,4.0)+pow(dy_s1,4.0)+pow(dz_s1,4.0))*(pow(div_s1,-4.0))
				//);
				
				term1x = (s1x[ip*NDPY*NDPZ+j*NDPZ+k] - s1x[im*NDPY*NDPZ+j*NDPZ+k])/(2.0*dx);
				term1y = (s1y[i*NDPY*NDPZ+jp*NDPZ+k] - s1y[i*NDPY*NDPZ+jm*NDPZ+k])/(2.0*dy);
				term1z = (s1z[i*NDPY*NDPZ+j*NDPZ+kp] - s1z[i*NDPY*NDPZ+j*NDPZ+km])/(2.0*dz);

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
				term2=4.0*www*s1*(1.0-s1)*(s1-0.5-dF+anois*(DRND(1)-0.5));

				//Phase-field evolution equation [equation (4.25)]
				//s1ddtt=-pmobi*(s1kai+s1kais);
				//other test: Eq.(4.21)
				s1ddtt=pmobi*(term1x + term1y + term1z + term2);

				//s1h2[i][j][k]=s1+s1ddtt*delt;	//Phase field time evolution (explicit method)
				//s1h2[i][j][k]=s1+s1ddtt*delt+anois*(DRND(1)-0.5)*s1*(1.0-s1);
				//if(s1h2[i][j][k]>=1.0){s1h2[i][j][k]=1.0;}//Phase field domain correction
				//if(s1h2[i][j][k]<=0.0){s1h2[i][j][k]=0.0;}
				s1h2[i*NDPY*NDPZ+j*NDPZ+k]=s1+s1ddtt*delt;	//Phase field time evolution (explicit method)
				if(s1h2[i*NDPY*NDPZ+j*NDPZ+k]>=1.0){s1h2[i*NDPY*NDPZ+j*NDPZ+k]=1.0;}//Phase field domain correction
				if(s1h2[i*NDPY*NDPZ+j*NDPZ+k]<=0.0){s1h2[i*NDPY*NDPZ+j*NDPZ+k]=0.0;}

				//Thermal diffusion equation [equation (4.30)]
				Tddtt=( cndct*( 
						(Tip+Tim-2.0*TT)/dx/dx
					   +(Tjp+Tjm-2.0*TT)/dy/dy
					   +(Tkp+Tkm-2.0*TT)/dz/dz
					  ) +30.0*s1*(1.0-s1)*s1*(1.0-s1)*rlate*s1ddtt )/speht;

				//Th2[i][j][k]=Th[i][j][k]+Tddtt*delt;		//Time evolution of temperature field (explicit method)
				Th2[i*NDPY*NDPZ+j*NDPZ+k]=Th[i*NDPY*NDPZ+j*NDPZ+k]+Tddtt*delt;	//Time evolution of temperature field (explicit method)

				dami:;
			}
		}
	}

//**** move field auxiliary array to main array ************************************
	for(i=0;i<=ndx;i++){
	    //for(j=0;j<=nd;j++){  s1h[i][j][k]=s1h2[i][j][k]; Th[i][j][k]=Th2[i][j][k]; }
		for(j=0;j<=ndy;j++){
			for(k=0;k<=ndz;k++){
				s1h[i*NDPY*NDPZ+j*NDPZ+k]=s1h2[i*NDPY*NDPZ+j*NDPZ+k];
				Th[i*NDPY*NDPZ+j*NDPZ+k]=Th2[i*NDPY*NDPZ+j*NDPZ+k];
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
void ini000(double *s1h, double *Th, int NDPX, int NDPY, int NDPZ)
{
	int i, j, k;
 	srand(time(NULL)); // random number initialization
	int ndx=NDPX-1, ndxm=NDPX-2, ndx2=(NDPX-1)/2;
	int ndy=NDPY-1, ndym=NDPY-2, ndy2=(NDPY-1)/2;
	int ndz=NDPZ-1, ndzm=NDPZ-2, ndz2=(NDPZ-1)/2;

	for(i=0;i<=ndx;i++){
		for(j=0;j<=ndy;j++){
			for(k=0;k<=ndz;k++){
				//s1h[i][j]=0.0;
				//if((i<10.)&&(j<10.)){s1h[i][j][k]=0.9;}
				//if((i*i+j*j)<20.){s1h[i][j][k]=0.9;}
				s1h[i*NDPY*NDPZ+j*NDPZ+k]=0.0;
				if((i*i+j*j+k*k)<20.0){
					s1h[i*NDPY*NDPZ+j*NDPZ+k]=0.9;
				}
			}
		}
	}

	for(i=0;i<=ndx;i++){
		for(j=0;j<=ndy;j++){
			for(k=0;k<=ndz;k++){
				//Th[i][j][k]=Tini+s1h[i][j][k]*(Tm-Tini);
				Th[i*NDPY*NDPZ+j*NDPZ+k]=Tini+s1h[i*NDPY*NDPZ+j*NDPZ+k]*(Tm-Tini);
			}
		}
	}

}

//************ data save subroutine *******************************
void datsave(double *s1h, double *Th, int NDPX, int NDPY, int NDPZ)
{
	FILE		*stream;	//Stream pointer setting
	int 		i, j, k;			//integer
	int ndx=NDPX-1, ndxm=NDPX-2, ndx2=(NDPX-1)/2;
	int ndy=NDPY-1, ndym=NDPY-2, ndy2=(NDPY-1)/2;
	int ndz=NDPZ-1, ndzm=NDPZ-2, ndz2=(NDPZ-1)/2;

	stream = fopen("test.dat", "w");	//Open the file to write to in append mode
	fprintf(stream, "%f\n", time1);		//Save repeat count
	for(i=0;i<=ndx;i++){
		for(j=0;j<=ndy;j++){
			for(k=0;k<=ndz;k++){
				//fprintf(stream, "%e  %e  ", s1h[i][j][k], Th[i][j][k]);//Conservation of local fields
				fprintf(stream, "%e  %e  ", s1h[i*NDPY*NDPZ+j*NDPZ+k], Th[i*NDPY*NDPZ+j*NDPZ+k]);//Conservation of local fields
			}
		}
	}
	fprintf(stream, "\n");	//writing a newline
	fclose(stream);					//close file
}

void datsave_paraview(double *s1h, double *Th, int NDPX, int NDPY, int NDPZ)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j, k;
	int ndx=NDPX-1, ndxm=NDPX-2, ndx2=(NDPX-1)/2;
	int ndy=NDPY-1, ndym=NDPY-2, ndy2=(NDPY-1)/2;
	int ndz=NDPZ-1, ndzm=NDPZ-2, ndz2=(NDPZ-1)/2;
	
	iout = iout + 1;
	printf("pf_result%06d.vtk \n",iout);
	sprintf(fName,"pf_result%06d.vtk",iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndx+1),(ndy+1),(ndz+1));
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %16d \n",((ndx+1)*(ndy+1)*(ndz+1)));
	fprintf(fp,"SCALARS Phase_field float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndz;k++){
		for(j=0;j<=ndy;j++){
			for(i=0;i<=ndx;i++){
				fprintf(fp,"%10.6f\n", s1h[i*NDPY*NDPZ+j*NDPZ+k]);
			}
		}
	}
	fprintf(fp,"SCALARS Temperature float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndz;k++){
		for(j=0;j<=ndy;j++){
			for(i=0;i<=ndx;i++){
				fprintf(fp,"%10.6f\n", Th[i*NDPY*NDPZ+j*NDPZ+k]);
			}
		}
	}
	fclose(fp);
}
