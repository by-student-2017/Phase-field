#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

//#include "wingxa.h"

#define DRND(x) ((double)(x)/RAND_MAX*rand())

#define ND 64
#define INXY 400					//�`��window�P�ӂ̃s�N�Z���T�C�Y

	int nd=ND, ndm=ND-1, nd2=ND/2; 	//�g�D�̕����A�g�D�̕���-1
	double RR=8.3145;
	double PI=3.141592654;
	double c0, temp, time1;
	double ch[ND][ND][ND];			//�g�D���̔Z�x�f�|�^�z��
	int iout;

	void ini_field();				//�����Z�x�g�ݒ�T�u���|�`��
	//void graph_c();				//�O���t�\���T�u���|�`��
	void datsave();				//�Z�x�f�|�^�ۑ��T�u���|�`��
	void datsave_paraview();		//�Z�x�f�|�^�ۑ��T�u���|�`��
	//void datin();					//�����Z�x�g�ǂݍ��݃T�u���|�`��

//******* Main program ******************************************************
int main(void)
{

	double ck[ND][ND][ND];			//�g�U�|�e���V����
	double ch2[ND][ND][ND];
	double mu_chem, mu_str, mu_surf;
	int   i, j, k, l, Nstep;
	double time1max; 				//�ő厞��
	double amob_c;					//���q�̈Փ��x�萔
	double c;						//�Z�x
	double cddtt;					//�Z�x�̑���
	double kapa_c;					//�Z�x���z�G�l���M�|�萔
	double al;						//�v�Z�̈�
	double b1;						//�����u���b�N�T�C�Y
	double L0;
	double delt;					//���ԍ���
	double Mx, My, Mz;				//�Փ��x

	double cip, cim, cjp, cjm, ckp, ckm;
	double ck1dev, ck2dev;
	double ck0, ckip, ckim, ckjp, ckjm, ckkp, ckkm;
	int   ip, im, jp, jm, kp, km;
	double sumc, dc0;
	double c_flu, c_flu0;			//�Z�x��̗h�炬�̑傫��

//********************************************************************************

	printf("DELT(0.005)=  "); scanf(" %lf",&delt);
//	delt=0.005;

	//printf("ca (0.4)= "); scanf(" %lf",&c0);
	c0=0.4;

	//printf("Mx(1.0) = ");	scanf(" %lf",&Mx); 	//�ړ��x
	Mx=1.; 						//�W���l

	//printf("My(1.0) = ");	scanf(" %lf",&My); 	//�ړ��x
	My=1.; 						//�W���l

	//printf("Mz(1.0) = ");	scanf(" %lf",&Mz); 	//�ړ��x
	Mz=0.01; 						//�W���l

	al=64.0;		//nm
	//al=50.0;			//nm
	al=al*1.0e-9;	//m

	b1=al/nd;
	amob_c=1.;

	time1=0.; time1max=30001.;

	temp=1000.0;//K

	//L0=2.5e+04/RR/temp;
	L0=1.0e+03/RR/temp;
	kapa_c=5.0e-15/b1/b1/RR/temp;
	//Mx=1.0;  								//�Փ��x
	//My=Mx;  								//�Փ��x
	//Mz=Mx;  								//�Փ��x
	//My=0.0;  								//�Փ��x
	//Mz=0.0;  								//�Փ��x

	c_flu=c_flu0=0.0;						//�Z�x��̗h�炬�W��
	//c_flu=c_flu0=0.1;						//�Z�x��̗h�炬�W��

//*************************************************************************
	ini_field();
 	//gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);//�`��Window�\��

//**** start **************************************************************
Nstep = 200;
iout = -1;
start: ;
	printf("time: %f \n", time1);
	if((((int)(time1) % Nstep)==0)) {iout = iout + 1;}
	if((((int)(time1) % Nstep)==0)) {datsave();}
	if((((int)(time1) % Nstep)==0)) {datsave_paraview();}
	//if(time1==2000.){datsave();}
	//if((((int)(time1) % 50)==0)) {graph_c();} 

	//if(time1<1.0e4){c_flu=c_flu0;} else{c_flu=0.0;}

//******[�g�U�|�e���V�����̌v�Z]********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){

				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==ndm) {ip=0;}  if(i==0) {im=ndm;}
				if(j==ndm) {jp=0;}  if(j==0) {jm=ndm;}
				if(k==ndm) {kp=0;}  if(k==0) {km=ndm;}

				c=ch[i][j][k];

			 	mu_chem=32.0*L0*c*(1.0-c)*(1.0-2.0*c); 		//���w�|�e���V����
			 	//mu_chem=L0*(1.0-2.0*c)+log(c)-log(1.0-c); //���w�|�e���V����

				mu_surf=-2.*kapa_c*(ch[ip][j][k]+ch[im][j][k]
      							   +ch[i][jp][k]+ch[i][jm][k]
								   +ch[i][j][kp]+ch[i][j][km]-6.0*c);
				ck[i][j][k]=mu_chem+mu_surf;

			}
		}
	}

//******[�Z�x��̎��ԕω�]**************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==ndm) {ip=0;}  if(i==0) {im=ndm;}
				if(j==ndm) {jp=0;}  if(j==0) {jm=ndm;}
				if(k==ndm) {kp=0;}  if(k==0) {km=ndm;}

				cddtt=Mx*(ck[ip][j][k]+ck[im][j][k]-2.0*ck[i][j][k])
					 +My*(ck[i][jp][k]+ck[i][jm][k]-2.0*ck[i][j][k])
					 +Mz*(ck[i][j][kp]+ck[i][j][km]-2.0*ck[i][j][k]);

				ch2[i][j][k]=ch[i][j][k]+(cddtt+c_flu*(2.0*DRND(1.0)-1.0))*delt;
			}
		}
	}

//******[�Z�x��̎��x�̕␳]**********************************
	sumc=0.;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				sumc+=ch2[i][j][k];
			}
	   }
	}
    dc0=sumc/nd/nd/nd-c0;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
	    	for(k=0;k<=ndm;k++){
				ch[i][j][k]=ch2[i][j][k]-dc0;
				if(ch[i][j][k]>=1.0){ch[i][j][k]=1.0-1.0e-06;}
				if(ch[i][j][k]<=0.0){ch[i][j][k]=1.0e-06;}
	    	}
	    }
	}

//******[���ԑ���]*************************************************
	//if(keypress()){return 0;}
	time1=time1+1.;
	if (time1<time1max) {goto start;}

end:;
  return 0;
}

//************[�����Z�x�g]*****************************************
void ini_field()
{
	int i, j, k, ii, jj, kk, PN, ipn;
	int x1, y1, z1, r0;
	int flg;
	double rnd0, sumc, r; 
 	srand(time(NULL)); // ����������
	int nd01, nd02, dd1;

//������
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
				for(k=0;k<=ndm;k++){
					ch[i][j][k]=0.0;
			}
		}
	}


// 1: �X�s�[�m�[�_������
// 2: �j�`��-�����^����
// 3: �P�̋�(�̐ϕ����͂�����ŋN�Z���Ă���)
// 4: �P�̃g�[���X(�̐ϕ����͂�����ŋN�Z���Ă���)

	flg=1;//�ȏ�L�̂����ꂩ��I������B

	switch(flg){

	case 1: //�X�s�[�m�[�_������
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
				for(k=0;k<=ndm;k++){
					ch[i][j][k]=c0+0.1*(2.*DRND(1)-1.);
			}
		}
	}
	break;


	case 2: //�j�`��
	r0=nd/20;  //�ar0�̊j��u��
	//PN=50; //�j��PN�u��
	PN=c0/( 4./3.*PI*(1./20.)*(1./20.)*(1./20.) );

	for(ipn=1;ipn<=PN;ipn++){
		x1=nd*DRND(1); y1=nd*DRND(1); z1=nd*DRND(1);
		while(ch[x1][y1][z1]!=0.0){x1=nd*DRND(1); y1=nd*DRND(1); z1=nd*DRND(1);}
		for(i=-r0;i<=(ndm+r0);i++){
			ii=i; if(i<0){ii=nd+i;}  if(i>ndm){ii=i-nd;}
			for(j=-r0;j<=(ndm+r0);j++){
				jj=j; if(j<0){jj=nd+j;}  if(j>ndm){jj=j-nd;}
				for(k=-r0;k<=(ndm+r0);k++){
					kk=k; if(k<0){kk=nd+k;}  if(k>ndm){kk=k-nd;}
					r=sqrt( (double(i-x1))*(double(i-x1))
								 +(double(j-y1))*(double(j-y1))
								 +(double(k-z1))*(double(k-z1)) );
					if(r<=r0){ ch[ii][jj][kk]=1.0; } 
				}
			}
		}
	}

	sumc=0.;
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i][j][k]; } } }
  //c0=sumc/nd/nd/nd;
	break;


	case 3: //�P�̋�
	nd01=nd/3;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				dd1=(i-nd2)*(i-nd2)+(j-nd2)*(j-nd2)+(k-nd2)*(k-nd2);
				if(dd1<(nd01*nd01)){ch[i][j][k]=0.9;} else{ch[i][j][k]=0.2;}
			}
		}
	}
	sumc=0.;
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i][j][k]; } } }
  c0=sumc/nd/nd/nd;
	break;


	case 4: //�P�̃g�[���X
	nd01=nd/3; nd02=nd/6;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=nd2-4;k<=nd2+4;k++){
				dd1=(i-nd2)*(i-nd2)+(j-nd2)*(j-nd2);
				if(dd1<(nd01*nd01)){ch[i][j][k]=0.9;} else{ch[i][j][k]=0.2;}
				if(dd1<(nd02*nd02)){ch[i][j][k]=0.2;}
			}
		}
	}
	sumc=0.;
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i][j][k]; } } }
  c0=sumc/nd/nd/nd;
	break;

	}//switch

}

//*******[�O���t]**************************************************
//void graph_c()
//{
//		int i, j, k, ii, jj, kk;
//		double col;
//		int ixmin=0, iymin=0, igx, igy, irad0;
//		int ixmax=INXY;
		int iymax=INXY;
//		double c, x, xmax, xmin, y, ymax, ymin, rad0;

   	//gcls(); //��ʃN���A
//		xmin=0.; xmax=1.;
//		ymin=0.; ymax=1.;

//		printf("time %f\n",time1);
//		rad0=1./nd/2.;
//		irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;

//	for(k=0;k<=nd;k++){
//		for(i=0;i<=nd;i++){
//			for(j=0;j<=nd;j++){
//				x=1./nd*i+rad0;
//				igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
//				y=1./nd*j+rad0;
//				igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
//				ii=i; jj=j, kk=k;
//				if(i==nd){ii=0;}  if(j==nd){jj=0;}  if(k==nd){kk=0;}
//				col=1.-ch[ii][jj][kk];
//				if(col>=1.){col=1.;}  if(col<=0.){col=0.;}
//				gcolor((int)(255*col),(int)(255*col),(int)(255*col));
//				grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);
//			}
//		}
//		swapbuffers();
//	}
//} 

//************[�Z�x�g�f�|�^�̕ۑ�]************************************
void datsave()
{
	FILE		*stream;
	char	fName[256];
	int 		i, j, k;

	sprintf(fName,"data_%06d.dat",iout);
	//stream = fopen("test3D_rod.dat", "a");
	stream = fopen(fName, "w");
	fprintf(stream, "%d %d %d \n", ndm, ndm, ndm);
	fprintf(stream, "%e\n", time1);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				fprintf(stream, "%e  ", ch[i][j][k]);
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
	
	sprintf(fName,"rodms_3D_result%06d.vtk",iout);
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
	fclose(fp);
}
//***********[�t�@�C������̏����Z�x�g����]****************************
//void datin()
//{
//	FILE		*datin0;
//	int 		i, j;

//	datin0 = fopen("test.dat", "r");
//	fscanf(datin0, "%lf", &time1);
//	for(i=0;i<=ndm;i++){
//		for(j=0;j<=ndm;j++){
//			fscanf(datin0, "%lf", &ch[i][j]);
//		}
//	}
//	fclose(datin0);
//}
