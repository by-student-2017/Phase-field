#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

//#include "wingxa.h"

#define DRND(x) ((double)(x)/RAND_MAX*rand())

//#define ND 30
//#define INXY 400				//•`‰æwindow‚P•Ó‚ÌƒsƒNƒZƒ‹ƒTƒCƒY

	//int nd=ND, ndm=ND-1; 		//‘gD‚Ì•ªŠ„A‘gD‚Ì•ªŠ„-1
	double rr=8.3145;
	double ca, temp, time1;
	//double ch[ND][ND][ND];		//‘gD“à‚Ì”Z“xƒf|ƒ^”z—ñ
	
	int iout;

	void shokiha_c(double *ch, int ND);	//‰Šú”Z“x”gİ’èƒTƒuƒ‹|ƒ`ƒ“
	//void graph_c();			//ƒOƒ‰ƒt•\¦ƒTƒuƒ‹|ƒ`ƒ“
	void datsave(double *ch, int ND);	//”Z“xƒf|ƒ^•Û‘¶ƒTƒuƒ‹|ƒ`ƒ“
	void datsave_paraview(double *ch, int ND);	//”Z“xƒf|ƒ^•Û‘¶ƒTƒuƒ‹|ƒ`ƒ“
	//void datin();				//‰Šú”Z“x”g“Ç‚İ‚İƒTƒuƒ‹|ƒ`ƒ“

//******* Main program ******************************************************
int main(void)
{
	int ND;
	int nd, ndm; 
	
	//double ck[ND][ND][ND];		//ŠgUƒ|ƒeƒ“ƒVƒƒƒ‹
	//double ch2[ND][ND][ND];
	double mu_chem, mu_str, mu_surf;
	int   i, j, k, l, Nstep;
	double time1max; 			//Å‘åŠÔ
	double amob_c;				//Œ´q‚ÌˆÕ“®“x’è”
	double c;					//”Z“x
	double cddtt;				//”Z“x‚Ì‘•ª
	double kapa_c;				//”Z“xŒù”zƒGƒlƒ‹ƒM|’è”
	double kapa_c0;
	double om0, om1, om2;		//Œ´qŠÔ‘ŠŒİì—pƒpƒ‰ƒ|ƒ^’è”
	double al;					//ŒvZ—Ìˆæ
	double b1;					//·•ªƒuƒƒbƒNƒTƒCƒY
	double delt;				//ŠÔ‚İ
	double a0;					//Fe‚ÌŠiq’è”
	double natom;				
	double vm0;					//ƒ‚ƒ‹‘ÌÏ
	double c11a, c12a, c44a;	//Fe‚Ì’e«—¦
	double c11a0, c12a0, c44a0;
	double c11b, c12b, c44b;	//Cr‚Ì’e«—¦
	double c11b0, c12b0, c44b0;
	double c11, c12, c44;		//‡‹à‚Ì•½‹Ï’e«—¦
	double y100;				//’e«ŠÖ”
	double eta;					//Šiqƒ~ƒXƒ}ƒbƒ`
	double Da, Db, Dab;			//©ŒÈŠgUŒW”
	double Dac, Dae, Dbc, Dbe;
	double Mc, dMc;				//ƒ‚ƒrƒŠƒeƒB|ŠÖ”

	double c1, c2, c2_cn=0.858;
	double Tc, Tn;
	double Tc1, Tc2, Tc3, Tc4;
	double Tn1, Tn2;
	double beta_c, beta_n;
	double beta_c1, beta_c2, beta_c3;
	double beta_n1, beta_n2;
	double c_alp;
	double c_alp1_c, c_alp2_c;
	double c_alp1_n, c_alp2_n;
	double d_Tc, d_Tn;
	double d_Tc1, d_Tc2, d_Tc3, d_Tc4;
	double d_Tn1;
	double d_beta_c, d_beta_n;
	double d_beta_c1, d_beta_c2;
	double d_beta_n1;
	double d_c_alp;
	double d_c_alp_c;
	double d_c_alp_n;
	double L0;
	double L01, L02;
	double fs=0.285;
	double kkf, kkp, sm;
	double d_kf, d_kp;
	double d_dg_cpm_eqm, d_dg_nm, d_gc;

	double cip, cim, cjp, cjm, ckp, ckm;
	double ck1dev, ck2dev;
	double ck0, ckip, ckim, ckjp, ckjm, ckkp, ckkm;
	int   ip, im, jp, jm, kp, km;
	double sumc, dca;

//********************************************************************************
	printf("---------------------------------\n");
	printf("read parameters form parameters.txt\n");
	FILE *fp;
	char name[40], comment[72];
	float param;
	float data[100];
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
	ca      = data[2];
	temp    = data[3];	// [K]
	al      = data[4];	// [nm]
	amob_c  = data[5];
	time1max= data[6];
	a0      = data[7];
	natom   = data[8];
	L01     = data[9];
	L02     = data[10];
	kapa_c0 = data[11];	//kapa_c=kapa_c0/b1/b1/RR/temp;
	Nstep   = int(data[12]);
	eta     = data[13];
	c11a0   = data[14];
	c12a0   = data[15];
	c44a0   = data[16];
	c11b0   = data[17];
	c12b0   = data[18];
	c44b0   = data[19];
	Dac     = data[20];
	Dae     = data[21];
	Dbc     = data[22];
	Dbe     = data[23];
	Tc1     = data[24];
	Tc2     = data[25];
	Tc3     = data[26];
	Tc4     = data[27];
	beta_c1 = data[28];
	beta_c2 = data[29];
	beta_c3 = data[30];
	c_alp1_c= data[31];
	c_alp2_c= data[32];
	d_Tc1   = data[33];
	d_Tc2   = data[34];
	d_Tc3   = data[35];
	d_Tc4   = data[36];
	d_beta_c1= data[37];
	d_beta_c2= data[38];
	d_c_alp_c= data[39];
	Tn1     = data[40];
	Tn2     = data[41];
	beta_n1 = data[42];
	beta_n2 = data[43];
	c_alp1_n= data[44];
	c_alp2_n= data[45];
	d_Tn1   = data[46];
	d_beta_n1=data[47];
	d_c_alp_n=data[48];
	printf("---------------------------------\n");
	//
	nd=ND;
	ndm=ND-1;
	
	//double ch[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	//a[z * ySize * xSize + y * xSize + x], a[i][j][k]=a[i*n*m + j*n + k]
	double *ck  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//ŠgUƒ|ƒeƒ“ƒVƒƒƒ‹
	double *ch  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//‘gD“à‚Ì”Z“xƒf|ƒ^”z—ñ
	double *ch2 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	
	//printf("DELT(0.05)=  "); scanf(" %lf",&delt);	//	delt=0.05;
	//printf("ca (0.4)= "); scanf(" %lf",&ca);	//	ca=0.4;
	//printf("Temp(K) (773.)= "); scanf(" %lf",&temp);	//	temp=773.;

	//al=30.;		// [nm]
	al=al*1.0e-9;	// [m]
	b1=al/nd;
	//amob_c=1.;
	time1=0.;
	//time1max=100001.;

	//a0=2.8664E-10;
	a0=a0*1.0e-9;	// [m]
	//vm0=6.02E23*a0*a0*a0/2.;
	vm0=6.02E23*a0*a0*a0/natom;

	//L0=(21020.8-9.31889*temp)/rr/temp;
	L0=(L01 + L02*temp)/rr/temp;

	//kapa_c=6.0e-15/b1/b1/rr/temp;
	kapa_c=kapa_c0/b1/b1/rr/temp;

	//eta=0.00614;

	//c11a=2.331e11*vm0/rr/temp;
	//c12a=1.3544e11*vm0/rr/temp;
	//c44a=1.1783e11*vm0/rr/temp;
	c11a=c11a0*vm0/rr/temp;
	c12a=c12a0*vm0/rr/temp;
	c44a=c44a0*vm0/rr/temp;

	//c11b=3.5e11*vm0/rr/temp;
	//c12b=0.678e11*vm0/rr/temp;
	//c44b=1.008e11*vm0/rr/temp;
	c11b=c11b0*vm0/rr/temp;
	c12b=c12b0*vm0/rr/temp;
	c44b=c44b0*vm0/rr/temp;

	c11=(1.0-ca)*c11a+ca*c11b;
	c12=(1.0-ca)*c12a+ca*c12b;
	c44=(1.0-ca)*c44a+ca*c44b;

	y100=c11+c12-2.*(c12*c12/c11);

	//Da=1.0e-4*exp(-294000./rr/temp);	//Fe‚Ì©ŒÈŠgUŒW”
	//Db=2.0e-5*exp(-308000./rr/temp);	//Cr‚Ì©ŒÈŠgUŒW”
	Da=Dac*exp(Dae/rr/temp);	//Fe‚Ì©ŒÈŠgUŒW”
	Db=Dbc*exp(Dbe/rr/temp);	//Cr‚Ì©ŒÈŠgUŒW”
	Dab=Db/Da;					//Fe‚Ì©ŒÈŠgUŒW”‚Å‹KŠi‰»

	//printf(" %f, %f %f\n", gc, dg_nm, dg_cpm_eqm);

//*************************************************************************
	shokiha_c(ch, ND);
	//datin();
 	//gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);//•`‰æWindow•\¦

//**** start **************************************************************
//Nstep = 200;
iout = -1;
start: ;
	printf("time: %f \n", time1);
	//if((((int)(time1) % Nstep)==0)) {datsave(ch, ND);}
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(ch, ND);}
	//if(time1==300.){datsave();}
	//if((((int)(time1) % 10)==0)) {graph_c();} 

//******[ŠgUƒ|ƒeƒ“ƒVƒƒƒ‹‚ÌŒvZ]********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){

				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==ndm) {ip=0;}  if(i==0) {im=ndm;}
				if(j==ndm) {jp=0;}  if(j==0) {jm=ndm;}
				if(k==ndm) {kp=0;}  if(k==0) {km=ndm;}

				//c=ch[i][j][k];
				//cip=ch[ip][j][k];  cim=ch[im][j][k];
				//cjp=ch[i][jp][k];  cjm=ch[i][jm][k];
				//ckp=ch[i][j][kp];  ckm=ch[i][j][km];
				c=ch[i*ND*ND+j*ND+k];
				cip=ch[ip*ND*ND+j*ND+k];  cim=ch[im*ND*ND+j*ND+k];
				cjp=ch[i*ND*ND+jp*ND+k];  cjm=ch[i*ND*ND+jm*ND+k];
				ckp=ch[i*ND*ND+j*ND+kp];  ckm=ch[i*ND*ND+j*ND+km];

				mu_surf=-2.*kapa_c*(cip+cim+cjp+cjm+ckp+ckm-6.*c);
				
				c1=1.-c ; c2=c;
				
			if(isnan(c)){printf("%f, %f, %f, %f \n",Tc, Tn, beta_c, beta_n); exit(0);} // check nan
			if(c<=c2_cn){
				//Tc=1043.*c1-310.*c2+(1207.3-321.2*(c2-c1))*c1*c2;
				Tc = Tc1*c1 + Tc2*c2 + (Tc3 + Tc4*(c2-c1))*c1*c2;
				if(Tc <= 0.0001){Tc=0.0001;} //for nan
				//beta_c=2.216*c1-0.4*c2+0.2525*c1*c2;
				beta_c = beta_c1*c1 + beta_c2*c2 + beta_c3*c1*c2;
				if(beta_c <= -0.9999){beta_c=-0.9999;} //for nan
				//c_alp=0.96*c1+1.0*c2;
				c_alp = c_alp1_c*c1 + c_alp2_c*c2;
				sm=c_alp*rr*log(beta_c+1.);
				kkf=4.*(1.-fs)*sm/(1.-exp(-4.));
				kkp=8.*fs*sm;
				//d_Tc=-1353.+(1207.3-321.2*(c2-c1))*(c1-c2)-642.4*c1*c2;
				d_Tc = d_Tc1 + (d_Tc2 + d_Tc3*(c2-c1))*(c1-c2) + d_Tc4*c1*c2;
				//d_beta_c=-2.616+0.2525*(c1-c2);
				d_beta_c = d_beta_c1 + d_beta_c2*(c1-c2);
				//d_c_alp=0.04;
				d_c_alp = d_c_alp_c;
				d_kf=4.*(1.-fs)/(1.-exp(-4.))*d_c_alp*rr*log(beta_c+1.)
					+4.*(1.-fs)/(1.-exp(-4.))*d_beta_c*c_alp*rr/(beta_c+1.);
				d_kp=8.*fs*d_c_alp*rr*log(beta_c+1.)
					+8.*fs*d_beta_c*c_alp*rr/(beta_c+1.);
				if(temp>=Tc){d_dg_cpm_eqm=-d_kp*Tc/64.*exp(8.*(1.-temp/Tc))
										  -d_Tc*kkp/64.*exp(8.*(1.-temp/Tc))
										  -d_Tc*kkp*temp/8./Tc*exp(8.*(1.-temp/Tc));}
					else    {d_dg_cpm_eqm=1./8.*d_kp*(temp-9.*Tc/8.)
										  -9.*kkp/64.*d_Tc
								  		  -0.25*d_kf*(-temp+3./4.*Tc
										   +Tc/4.*exp(-4.*(1.-temp/Tc)))
							 	 	  	  -kkf/4.*d_Tc*(0.75+0.25*exp(-4.*(1.-temp/Tc))
										 -temp/Tc*exp(-4.*(1.-temp/Tc)));}
			}
			else{
				//Tn=-1873.1*c1+310.*c2;
				Tn = Tn1*c1 + Tn2*c2;
				if(Tn <= 0.0001){Tn=0.0001;} //for nan
				//beta_n=-2.417*c1+0.4*c2;
				beta_n = beta_n1*c1 + beta_n2*c2;
				if(beta_n <= -0.9999){beta_n=-0.9999;} //for nan
				//c_alp=0.96*c1+1.0*c2;
				c_alp = c_alp1_n*c1 + c_alp2_n*c2;
				sm=c_alp*rr*log(beta_n+1.);
				kkf=4.*(1.-fs)*sm/(1.-exp(-4.));
				kkp=8.*fs*sm;
				//d_Tn=2183.1;
				d_Tn = d_Tn1;
				//d_beta_n=2.817;
				d_beta_n = d_beta_n1;
				//d_c_alp=0.04;
				d_c_alp = d_c_alp_n;
				d_kf=4.*(1.-fs)/(1.-exp(-4.))*d_c_alp*rr*log(beta_n+1.)
					+4.*(1.-fs)/(1.-exp(-4.))*d_beta_n*c_alp*rr/(beta_n+1.);
				d_kp=8.*fs*d_c_alp*rr*log(beta_n+1.)
					+8.*fs*d_beta_n*c_alp*rr/(beta_n+1.);
				if(temp>=Tn){d_dg_cpm_eqm=-d_kp*Tn/64.*exp(8.*(1.-temp/Tn))
										  -d_Tn*kkp/64.*exp(8.*(1.-temp/Tn))
										  -d_Tn*kkp*temp/8./Tn*exp(8.*(1.-temp/Tn));}
					else    {d_dg_cpm_eqm=1./8.*d_kp*(temp-9.*Tn/8.)
										  -9.*kkp/64.*d_Tn
								  		  -0.25*d_kf*(-temp+3./4.*Tn
										   +Tn/4.*exp(-4.*(1.-temp/Tn)))
								  	  	  -kkf/4.*d_Tn*(0.75+0.25*exp(-4.*(1.-temp/Tn))
										 -temp/Tn*exp(-4.*(1.-temp/Tn)));}
				}
				d_dg_nm=L0*(c1-c2);
				mu_chem=log(c2)-log(c1)+d_dg_nm+d_dg_cpm_eqm/rr/temp;
				mu_str=2.*eta*eta*y100*(c-ca);
				//ck[i][j][k]=mu_chem+mu_str+mu_surf;
				ck[i*ND*ND+j*ND+k]=mu_chem+mu_str+mu_surf;
			}
		}
	}

//******[”Z“xê‚ÌŠÔ•Ï‰»]**************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==ndm) {ip=0;}  if(i==0) {im=ndm;}
				if(j==ndm) {jp=0;}  if(j==0) {jm=ndm;}
				if(k==ndm) {kp=0;}  if(k==0) {km=ndm;}
				
				//cip=ch[ip][j][k];  cim=ch[im][j][k];
				//cjp=ch[i][jp][k];  cjm=ch[i][jm][k];
				//ckp=ch[i][j][kp];  ckm=ch[i][j][km];
				cip=ch[ip*ND*ND+j*ND+k];  cim=ch[im*ND*ND+j*ND+k];
				cjp=ch[i*ND*ND+jp*ND+k];  cjm=ch[i*ND*ND+jm*ND+k];
				ckp=ch[i*ND*ND+j*ND+kp];  ckm=ch[i*ND*ND+j*ND+km];
				
				//ck0=ck[i][j][k];
				//ckip=ck[ip][j][k];  ckim=ck[im][j][k];
				//ckjp=ck[i][jp][k];  ckjm=ck[i][jm][k];
				//ckkp=ck[i][j][kp];  ckkm=ck[i][j][km];
				ck0=ck[i*ND*ND+j*ND+k];
				ckip=ck[ip*ND*ND+j*ND+k];  ckim=ck[im*ND*ND+j*ND+k];
				ckjp=ck[i*ND*ND+jp*ND+k];  ckjm=ck[i*ND*ND+jm*ND+k];
				ckkp=ck[i*ND*ND+j*ND+kp];  ckkm=ck[i*ND*ND+j*ND+km];
				
				ck2dev=(ckip+ckim+ckjp+ckjm+ckkp+ckkm)-6.*ck0;
				Mc=(ca+Dab*(1.-ca))*ca*(1.-ca);
				cddtt=amob_c*Mc*ck2dev;
				//ch2[i][j][k]=ch[i][j][k]+cddtt*delt;
				ch2[i*ND*ND+j*ND+k]=ch[i*ND*ND+j*ND+k]+cddtt*delt;
			}
		}
	}

//******[”Z“xê‚Ìûx‚Ì•â³]**********************************
	sumc=0.;
	for(i=1;i<=ndm;i++){
	    for(j=1;j<=ndm;j++){
	    	for(k=1;k<=ndm;k++){
				//sumc+=ch2[i][j][k];
	    		sumc+=ch2[i*ND*ND+j*ND+k];
	    	}
	    }
	}
    dca=sumc/nd/nd/nd-ca;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
	    	for(k=0;k<=ndm;k++){
				//ch[i][j][k]=ch2[i][j][k]-dca;
				//if(ch[i][j][k]>=1.){ch[i][j][k]=0.9999;}
				//if(ch[i][j][k]<=0.){ch[i][j][k]=0.0001;}
				ch[i*ND*ND+j*ND+k]=ch2[i*ND*ND+j*ND+k]-dca;
				if(ch[i*ND*ND+j*ND+k]>=1.0){ch[i*ND*ND+j*ND+k]=0.9999;}
				if(ch[i*ND*ND+j*ND+k]<=0.0){ch[i*ND*ND+j*ND+k]=0.0001;}
	    	}
	    }
	}

//******[ŠÔ‘‰Á]*************************************************
	//if(keypress()){return 0;}
	time1=time1+1.;
	if (time1<time1max) {goto start;}

end:;
  return 0;
}

//************[‰Šú”Z“x”g]*****************************************
void shokiha_c(double *ch, int ND)
{
	int i, j, k;
	double rnd0; 
  	srand(time(NULL)); // —”‰Šú‰»
	int ndm=ND-1;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				rnd0=2.*DRND(1)-1.;
				//ch[i][j][k]=ca+0.01*rnd0;
				ch[i*ND*ND+j*ND+k]=ca+0.01*rnd0;
			}
		}
	}
}

//************[”Z“x”gƒf|ƒ^‚Ì•Û‘¶]************************************
void datsave(double *ch, int ND)
{
	FILE		*stream;
	int 		i, j, k;
	int ndm=ND-1;

	stream = fopen("test.dat", "a");
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

void datsave_paraview(double *ch, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j, k;
	int ndm=ND-1;
	
	iout = iout + 1;
	sprintf(fName,"sp_result%06d.vtk",iout);
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
	fclose(fp);
}
//***********[ƒtƒ@ƒCƒ‹‚©‚ç‚Ì‰Šú”Z“x”g“ü—Í]****************************
//void datin()
//{
//	FILE		*datin0;
//	int 		i, j, k;

//	datin0 = fopen("test.dat", "r");
//	fscanf(datin0, "%lf", &time1);
//	for(i=0;i<=ndm;i++){
//		for(j=0;j<=ndm;j++){
//			for(j=0;j<=ndm;j++){
//				fscanf(datin0, "%lf", &ch[i][j][k]);
//			}
//		}
//	}
//	fclose(datin0);
//}
