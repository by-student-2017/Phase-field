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
	//double s1h[NDP][NDP], Th[NDP][NDP];	//phase field, temperature field
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
	int   i, j, k, l;				//integer
	int   ip, im, jp, jm;			//integer
	double time1max;				//Maximum calculation count (calculation end count)

	double s1;						//phase field
	double s1ip, s1im, s1jp, s1jm;	//s1 left, right, top, bottom phase field
	double s1ipjp, s1ipjm, s1imjp, s1imjm;//Diagonal side phase field of s1
	double s1ddtt;					//Time variation of s1 (left side of evolution equation)

	double TT;						//temperature field
	double Tip, Tim, Tjp, Tjm;		//Left, right, top, bottom temperature of TT
	double Tddtt;					//Temporal change in temperature (left side of thermal diffusion equation)

	double ep, ep1p, ep2p;			//Gradient energy coefficient related variables
	double dx_s1, dy_s1, l_dxdy;	//Variables related to spatial differentiation of the phase field
	double dxx_s1, dyy_s1, dxy_s1;	//Variables related to spatial differentiation of the phase field

	double al;				//length of one side of computational domain
	double dx, dy;			//Difference grid size (x direction)
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
	gamma   = data[4];
	ram     = data[5];
	j_fold  = data[6];
	astre   = data[7];
	cndct   = data[8];
	speht   = data[9];
	rlate   = data[10];
	Tm      = data[11];
	Tini    = data[12];
	skine   = data[13];
	th0     = data[14];
	aaac    = data[15];
	wwwc    = data[16];
	pmobic  = data[17];
	anois   = data[18];
	time1max= int(data[19]);
	Nstep   = int(data[20]);
	printf("---------------------------------\n");
	//
	nd=NDP-1;
	ndm=NDP-2;
	nd2=(NDP-1)/2;
	//
	double *s1h  = (double *)malloc(sizeof(double)*( NDP*NDP ));	//phase field
	double *Th   = (double *)malloc(sizeof(double)*( NDP*NDP ));	//temperature field
	double *s1h2 = (double *)malloc(sizeof(double)*( NDP*NDP ));
	double *Th2  = (double *)malloc(sizeof(double)*( NDP*NDP ));
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
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==nd){ip=ndm;}  if(i==0){im=1;}
			if(j==nd){jp=ndm;}  if(j==0){jm=1;}

			//s1=s1h[i][j];//phase field
			//s1ip=s1h[ip][j];  s1im=s1h[im][j];  s1jp=s1h[i][jp];  s1jm=s1h[i][jm];
			//s1ipjp=s1h[ip][jp]; s1ipjm=s1h[ip][jm];  s1imjp=s1h[im][jp];  s1imjm=s1h[im][jm];
			s1=s1h[i*NDP+j];//phase field
			s1ip=s1h[ip*NDP+j];  s1im=s1h[im*NDP+j];  s1jp=s1h[i*NDP+jp];  s1jm=s1h[i*NDP+jm];
			s1ipjp=s1h[ip*NDP+jp]; s1ipjm=s1h[ip*NDP+jm];  s1imjp=s1h[im*NDP+jp];  s1imjm=s1h[im*NDP+jm];

			//TT=Th[i][j];//temperature field
			//Tip=Th[ip][j];  Tim=Th[im][j];  Tjp=Th[i][jp];  Tjm=Th[i][jm];
			TT=Th[i*NDP+j];//temperature field
			Tip=Th[ip*NDP+j];  Tim=Th[im*NDP+j];  Tjp=Th[i*NDP+jp];  Tjm=Th[i*NDP+jm];

//----- Deciding When to Skip Calculations ----------------------------------------------
			dami1=fabs(s1+s1ip+s1im+s1jp+s1jm);  dami2=fabs(TT+Tip+Tim+Tjp+Tjm-5.0*Tini);
			if( (dami1<=1.0e-20)&&(dami2<=1.0e-20) ){
				//s1h2[i][j]=s1h[i][j];  Th2[i][j]=Th[i][j];
				s1h2[i*NDP+j]=s1h[i*NDP+j];  Th2[i*NDP+j]=Th[i*NDP+j];
				goto dami;
			}
//---------------------------------------------------------------------------------

			dx_s1=(s1ip-s1im)/2.0/dx;  				//Spatial first derivative of the phase field
			dy_s1=(s1jp-s1jm)/2.0/dy;
			dxx_s1=(s1ip+s1im-2.0*s1)/dx/dx;		//Spatial second derivative of the phase field
			dyy_s1=(s1jp+s1jm-2.0*s1)/dy/dy;
			dxy_s1=(s1ipjp+s1imjm-s1imjp-s1ipjm)/4.0/dx/dy;
			th=atan(dy_s1/(dx_s1+1.0e-20));			//angle normal to interface [equation (4.24)]

			ep=aaa*(1.0+astre*cos(j_fold*(th-th0)));	//square root of gradient energy coefficient [equation (4.23)]
			ep1p=-aaa*astre*j_fold*sin(j_fold*(th-th0));	//first derivative of ep with respect to angle
			ep2p=-aaa*astre*j_fold*j_fold*cos(j_fold*(th-th0));	//Second derivative of ep with respect to angle

			s1kais=-ep*ep*(dxx_s1+dyy_s1)
					-ep*ep1p*((dyy_s1-dxx_s1)*sin(2.0*th)+2.0*dxy_s1*cos(2.0*th))
					 +0.5*(ep1p*ep1p+ep*ep2p)*(2.0*dxy_s1*sin(2.0*th)
					-dxx_s1-dyy_s1-(dyy_s1-dxx_s1)*cos(2.0*th));
			//gradient potential [part of equation (4.29)]

			dF=15.0/(2.0*www)*rlate*(TT-Tm)/Tm*s1*(1.0-s1);//chemical driving force [part of formula (4.29)]
			s1kai=4.0*www*s1*(1.0-s1)*(0.5-s1+dF+anois*(DRND(1)-0.5));
			//chemical potential [part of equation (4.29)]

			s1ddtt=-pmobi*(s1kai+s1kais);	//Phase-field evolution equation [equation (4.25)]

			//s1h2[i][j]=s1+s1ddtt*delt;	//Phase field time evolution (explicit method)
			//s1h2[i][j]=s1+s1ddtt*delt+anois*(DRND(1)-0.5)*s1*(1.0-s1);
			//if(s1h2[i][j]>=1.0){s1h2[i][j]=1.0;}//Phase field domain correction
			//if(s1h2[i][j]<=0.0){s1h2[i][j]=0.0;}
			s1h2[i*NDP+j]=s1+s1ddtt*delt;	//Phase field time evolution (explicit method)
			if(s1h2[i*NDP+j]>=1.0){s1h2[i*NDP+j]=1.0;}//Phase field domain correction
			if(s1h2[i*NDP+j]<=0.0){s1h2[i*NDP+j]=0.0;}

			Tddtt=( cndct*( (Tip+Tim-2.0*TT)/dx/dx+(Tjp+Tjm-2.0*TT)/dy/dy )
      	                    +30.0*s1*(1.0-s1)*s1*(1.0-s1)*rlate*s1ddtt )/speht;
			//Thermal diffusion equation [equation (4.30)]
			//Th2[i][j]=Th[i][j]+Tddtt*delt;		//Time evolution of temperature field (explicit method)
			Th2[i*NDP+j]=Th[i*NDP+j]+Tddtt*delt;	//Time evolution of temperature field (explicit method)

			dami:;
		}
	}

//**** move field auxiliary array to main array ************************************
	for(i=0;i<=nd;i++){
	    //for(j=0;j<=nd;j++){  s1h[i][j]=s1h2[i][j]; Th[i][j]=Th2[i][j]; }
		for(j=0;j<=nd;j++){  s1h[i*NDP+j]=s1h2[i*NDP+j]; Th[i*NDP+j]=Th2[i*NDP+j]; }
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
	int i, j;
 	srand(time(NULL)); // random number initialization
	int nd=NDP-1, ndm=NDP-2, nd2=(NDP-1)/2;

	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			//s1h[i][j]=0.0;
			//if((i<10.)&&(j<10.)){s1h[i][j]=0.9;}
			//if((i*i+j*j)<20.){s1h[i][j]=0.9;}
			s1h[i*NDP+j]=0.0;
			if((i*i+j*j)<20.0){s1h[i*NDP+j]=0.9;}
		}
	}

	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			//Th[i][j]=Tini+s1h[i][j]*(Tm-Tini);
			Th[i*NDP+j]=Tini+s1h[i*NDP+j]*(Tm-Tini);
		}
	}

}

//************ data save subroutine *******************************
void datsave(double *s1h, double *Th, int NDP)
{
	FILE		*stream;	//Stream pointer setting
	int 		i, j;			//integer
	int nd=NDP-1, ndm=NDP-2, nd2=(NDP-1)/2;

	stream = fopen("test.dat", "a");	//Open the file to write to in append mode
	fprintf(stream, "%f\n", time1);		//Save repeat count
	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			//fprintf(stream, "%e  %e  ", s1h[i][j], Th[i][j]);//Conservation of local fields
			fprintf(stream, "%e  %e  ", s1h[i*NDP+j], Th[i*NDP+j]);//Conservation of local fields
		}
	}
	fprintf(stream, "\n");	//writing a newline
	fclose(stream);					//close file
}

void datsave_paraview(double *s1h, double *Th, int NDP)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j;
	int nd=NDP-1, ndm=NDP-2, nd2=(NDP-1)/2;
	
	iout = iout + 1;
	printf("pf_result%06d.vtk \n",iout);
	sprintf(fName,"pf_result%06d.vtk",iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),(1));
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*1));
	fprintf(fp,"SCALARS Phase_field float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e  ", s1h[i][j], Th[i][j]);//Conservation of local fields
			//fprintf(fp,"%10.6f\n", s1h[i][j]);
			fprintf(fp,"%10.6f\n", s1h[i*NDP+j]);
		}
	}
	fprintf(fp,"SCALARS Temperature float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e  ", s1h[i][j], Th[i][j]);//Conservation of local fields
			//fprintf(fp,"%10.6f\n", Th[i][j]);
			fprintf(fp,"%10.6f\n", Th[i*NDP+j]);
		}
	}
	fclose(fp);
}
