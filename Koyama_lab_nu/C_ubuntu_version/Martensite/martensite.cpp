#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

//#include "wingxa.h"

#define DRND(x) ((double)(x)/RAND_MAX*rand())
#define ND 128
#define IG 7
#define INXY 400			//�`��window�P�ӂ̃s�N�Z���T�C�Y

	int nd=ND, ndm=ND-1; 	//�g�D�̕�����
	int nd2=ND/2;			//(�g�D�̕�����)�^2�F�t�|���G�ϊ����Ŏg�p
	int ig=IG;				//2^ig=ND
	double PI=3.14159;		//�~����
	double rr=8.3145;		//�K�X�萔
	double time1;			//����

	double s1h[ND][ND], s2h[ND][ND];	//�}���e���T�C�g��phase field

	double qs;				//�t�|���G�ϊ�(qs:-1)�Ƌt�t�|���G�ϊ�(qs:1)�̋��
	double xi[ND][ND], xr[ND][ND], xif[ND], xrf[ND];//�t�|���G�ϊ��̎����E�����z��
	double s[ND],c[ND];	//sin��cos�̃e�[�u��
	int ik[ND];				//�r�b�g���]�z��
	int iout;

	void shokiha();			//�����g�ݒ�T�u���|�`��
	void graph_s1();		//�����x��`��T�u���|�`��
	void table();			//sin��cos�̃e�[�u���쐬�T�u���|�`��
	void fft();				//�P�����e�e�s
	void rcfft();			//�Q�����e�e�s
	void datsave_gnuplot();	//�f�|�^�ۑ��T�u���|�`��
	void datsave_paraview();//�f�|�^�ۑ��T�u���|�`��
	void datin();			//�����g�ǂݍ��ݗp�T�u���|�`��

//******* Main program ******************************************
int main(void)
{
	double s1, s2;			//�}���e���T�C�g��phase field
	double ep11h0[ND][ND], ep22h0[ND][ND];		//�g�D���̕ϑԘc
	double ep11qrh0[ND][ND],	ep11qih0[ND][ND];		//�S���c�ϓ��ʂ̃t�[���G�ϊ�
	double ep22qrh0[ND][ND],	ep22qih0[ND][ND];		//�S���c�ϓ��ʂ̃t�[���G�ϊ�
	double s1k_chem, s1k_str, s1k_su[ND][ND];			//�|�e���V����
	double s2k_chem, s2k_str, s2k_su[ND][ND];			//�|�e���V����
	double c11, c12, c44, lam0, mu0, nu0; //�e���萔
	double eta_s1[4][4], eta_s2[4][4];	//eigen�c����
	double ec11[ND][ND], ec22[ND][ND];	//�S���c�ϓ��ʁi����ԁj
	double ep11T, ep22T;
	double ep11_0, ep22_0;				//�g�D���̕ϑԘc�̕��ϒl
	double ep11_a, ep22_a, ep12_a, ep21_a;//�O�͂ɋN������c
	double sig11_a, sig22_a;			//�O��
	double Z11ep, Z12ep, Z21ep, Z22ep;	//�t�[���G�ϊ����̌W��
	double sum11, sum22;				//s1��s2�̋�Ԑϕ�
	double s1ddtt, s2ddtt;				//s1��s2�̎��ԕω���
	double el_fac;						//�e���萔�K�i���ϐ�

	int   i, j, k, l, ii, jj, kk, iii, jjj;		//����
	int   p, q, m, n;					//����
	int   ip, im, jp, jm, Nstep;		//����
	double al, temp, delt;				//�v�Z�̈�A���ԁA���Ԃ�����
	double time1max;					//�ő厞�ԁi�v�Z���~�߂�ۂɎg�p�j
	double b1, vm0, atom_n;				//�K�i�������A�����̐ρA�P�ʖE���̌��q��
	double smob;						//�����ϑԂ̊ɘa�W��
	double nxx, nyy, nxy, alnn;			//�t��Ԃ̊�{�x�N�g���̐ρA�m����

	double AA0, AA1, AA2, AA3;			//���w�I�쓮�͒萔
	double a1_c, b1_c, c1_c;			//�i�q�萔
	double a1_t, b1_t, c1_t;			//�i�q�萔
	double kappa_s1, kappa_s2;			//���z�G�l���M�|�萔
	double ds_fac;						//�����ϑԂ̗h�炬�W��

//****** reg data ****************************************

	printf("DELT(0.1)=  ");	scanf(" %lf",&delt);	//���Ԃ����ݓ���
	//delt=0.1;

	temp=500.0;			//���x
	al=250.0*1.0E-09;	//�v�Z�̈�[m]
	b1=al/nd;			//�����u���b�N�̒���

	time1=-10.0;		//�����ݒ莞��
	time1max=1.0+1.0e+07;	//�v�Z���Ԃ̍ő�l

	smob=1.0;			//�����ϑԂ̊ɘa�W��
	ds_fac=0.01;		//�����ϑԂ̗h�炬�W��

	AA0=1000.0/rr/temp;	//�}���e���T�C�g�ϑԂ̉��w�I�쓮��
	AA1=10.0;  AA2=3.0*AA1+12.0;  AA3=2.0*AA1+12.0;//���w�I�쓮�͒萔

	kappa_s1=kappa_s2=5.0e-15/rr/temp/b1/b1;//���z�G�l���M�|�萔

	a1_c=b1_c=c1_c=3.563E-10;//�i�q�萔
	//a1_t=b1_t=3.541E-10; c1_t=0.5*7.216E-10;
	atom_n=4.0;  vm0=6.02E23*a1_c*b1_c*c1_c/atom_n;//�����̐ς̌v�Z�ifcc������j

//*** s1��̕ϑԘc�̐ݒ� ***
	eta_s1[1][1]=0.08; eta_s1[2][2]=-0.04;
	eta_s1[3][3]=0.;
	eta_s1[1][2]=eta_s1[2][1]=eta_s1[1][3]=eta_s1[3][1]=eta_s1[2][3]=eta_s1[3][2]=0.;

//*** s2��̕ϑԘc�̐ݒ� ***
	eta_s2[1][1]=eta_s1[2][2];
	eta_s2[2][2]=eta_s1[1][1];
	eta_s2[3][3]=0.;
	eta_s2[1][2]=eta_s2[2][1]=eta_s2[1][3]=eta_s2[3][1]=eta_s2[2][3]=eta_s2[3][2]=0.;

//***** Ni�̒e���萔 ****************************
	el_fac=1.0E+11*vm0/rr/temp;
  c11=2.508*el_fac;
  c44=1.235*el_fac;
  c12=1.500*el_fac;
  //c12=c11-2.0*c44;
	lam0=c12;		mu0=c44;//���[���̒萔
	nu0=lam0/2.0/(lam0+mu0);//�|�A�\����
	printf("nu0= %f  \n", nu0);

//*** �O�͂̐ݒ� ***
 	sig22_a=0.;//0�ɐݒ�
	ep11_a=-0.5*lam0/mu0/(3.0*lam0+2.0*mu0)*sig22_a;
	ep22_a=(lam0+mu0)/mu0/(3.0*lam0+2.0*mu0)*sig22_a;
	ep12_a=ep21_a=0.0;

//*** sin�����cos�e�|�u������я����g�̐ݒ� ***************

 	//gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);//�`��Window�\��

	table();
	shokiha();

//**** ���C���v���O���� ******************************
iout = -1;
Nstep = 10;
start: ;

	//if(time1<=100.){Nstep=5;} else{Nstep=200;}
	printf("time: %f \n", time1);
	//if((((int)(time1) % Nstep)==0)) {datsave_gnuplot();} //phase field�̕ۑ�
	if((((int)(time1) % Nstep)==0)) {datsave_paraview();} //phase field�̕ۑ�
	//if((((int)(time1) % 10)==0)) {graph_s1();} //phase field�̕`��
	//if((((int)(time1) % 100)==0)) {datsave();} //phase field�̕ۑ�

//***** ���z�|�e���V���� ***********************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
			s1k_su[i][j]=-kappa_s1*(s1h[ip][j]+s1h[im][j]+s1h[i][jp]+s1h[i][jm]-4.0*s1h[i][j]);
			s2k_su[i][j]=-kappa_s2*(s2h[ip][j]+s2h[im][j]+s2h[i][jp]+s2h[i][jm]-4.0*s2h[i][j]);
		}
	}

//**** �ϑԘc��̃t�|���G�ϊ� ep11 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			xr[i][j]=ep11h0[i][j]=eta_s1[1][1]*s1h[i][j]+eta_s2[1][1]*s2h[i][j];
			xi[i][j]=0.;
		}
	}
	qs=-1.; rcfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){ 
			ep11qrh0[i][j]=xr[i][j];
			ep11qih0[i][j]=xi[i][j];
		}
	}
	ep11qrh0[0][0]=ep11qih0[0][0]=0.;

//**** �ϑԘc��̃t�|���G�ϊ� ep22 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			xr[i][j]=ep22h0[i][j]=eta_s1[2][2]*s1h[i][j]+eta_s2[2][2]*s2h[i][j];
			xi[i][j]=0.;
		}
	}
	qs=-1.; rcfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){ 
			ep22qrh0[i][j]=xr[i][j];
			ep22qih0[i][j]=xi[i][j];
		}
	}
	ep22qrh0[0][0]=ep22qih0[0][0]=0.;

//*** �ϑԘc��̕��ϒl�̎Z�o ***
	sum11=sum22=0.;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			sum11+=ep11h0[i][j];  sum22+=ep22h0[i][j];
	  }
	}
  ep11_0=sum11/nd/nd;  ep22_0=sum22/nd/nd;

//***** �S���c�ϓ���ec11�̌v�Z *************************************
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj);
			if(alnn==0.){alnn=1.;}
			nxx=(double)ii/alnn*(double)ii/alnn;
			nyy=(double)jj/alnn*(double)jj/alnn;
			Z11ep=nxx*(2.0*(1.0-nu0)-nxx-nu0/(1.0-nu0)*nyy)/(1.0-2.0*nu0);
			Z12ep=nxx*(2.0*nu0      -nyy-nu0/(1.0-nu0)*nxx)/(1.0-2.0*nu0);
			xr[i][j]=Z11ep*ep11qrh0[i][j]+Z12ep*ep22qrh0[i][j];
			xi[i][j]=Z11ep*ep11qih0[i][j]+Z12ep*ep22qih0[i][j];
	 }
	}
	qs=1.; rcfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ec11[i][j]=xr[i][j];
		}
	}

//***** �S���c�ϓ���ec22�̌v�Z *****************************
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj);
			if(alnn==0.){alnn=1.;}
			nxx=(double)ii/alnn*(double)ii/alnn;
			nyy=(double)jj/alnn*(double)jj/alnn;
			Z21ep=nyy*(2.0*nu0      -nxx-nu0/(1.0-nu0)*nyy)/(1.0-2.0*nu0);
			Z22ep=nyy*(2.0*(1.0-nu0)-nyy-nu0/(1.0-nu0)*nxx)/(1.0-2.0*nu0);
			xr[i][j]=Z21ep*ep11qrh0[i][j]+Z22ep*ep22qrh0[i][j];
			xi[i][j]=Z21ep*ep11qih0[i][j]+Z22ep*ep22qih0[i][j];
	 }
	}
	qs=1.; rcfft();
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ec22[i][j]=xr[i][j];
		}
	}

//******  �|�e���V�����̌v�Z ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){

			s1=s1h[i][j];  	s2=s2h[i][j];

//******  ���w�|�e���V�����̌v�Z ********************************
			s1k_chem=AA0*s1*(AA1-AA2*s1+AA3*(s1*s1+s2*s2));
			s2k_chem=AA0*s2*(AA1-AA2*s2+AA3*(s1*s1+s2*s2));

//******  �e���|�e���V�����̌v�Z ********************************

			ep11T=ep11h0[i][j]-ep11_0-ec11[i][j]-ep11_a;
			ep22T=ep22h0[i][j]-ep22_0-ec22[i][j]-ep22_a;

			s1k_str=ep11T*((lam0+2.0*mu0)*eta_s1[1][1]+lam0*eta_s1[2][2])
						 +ep22T*((lam0+2.0*mu0)*eta_s1[2][2]+lam0*eta_s1[1][1]);
			s2k_str=ep11T*((lam0+2.0*mu0)*eta_s2[1][1]+lam0*eta_s2[2][2])
						 +ep22T*((lam0+2.0*mu0)*eta_s2[2][2]+lam0*eta_s2[1][1]);

//****** phase field�̎��Ԕ��W�̌v�Z ********************************
			s1ddtt=-smob*(s1k_chem+s1k_su[i][j]+s1k_str);
			s2ddtt=-smob*(s2k_chem+s2k_su[i][j]+s2k_str);
			s1h[i][j]=s1h[i][j]+( s1ddtt+ds_fac*(2.0*DRND(1.)-1.0) )*delt;
			s2h[i][j]=s2h[i][j]+( s2ddtt+ds_fac*(2.0*DRND(1.)-1.0) )*delt;

//--- s�̕ψ�(0<=s<=1)�̕␳ ---
			if(s1h[i][j]>=1.0){s1h[i][j]=1.0;}  if(s1h[i][j]<=0.0){s1h[i][j]=0.0;}
			if(s2h[i][j]>=1.0){s2h[i][j]=1.0;}  if(s2h[i][j]<=0.0){s2h[i][j]=0.0;}
		}
	}

	//if(keypress()){return 0;}

	time1=time1+1.0;//���Ԃ̑���
	if(time1<time1max){goto start;}

end:;
  return 0;
}

//************ �����g�ݒ�T�u���|�`�� *************
void shokiha()
{
	int i, j;
  srand(time(NULL)); // ����������

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s1h[i][j]=DRND(1.0)/2.0; s2h[i][j]=DRND(1.0)/2.0;
			if(abs(j-nd2)<(nd/40)){s1h[i][j]=DRND(1.0); s2h[i][j]=DRND(1.0);}
		}
	}
}

//******* phase field�̕`�� ***************************************
//void graph_s1()
//{
//	int i, j, ii, jj;
//	double col, col_R, col_G, col_B, col_RG;
//	int ixmin=0, iymin=0, igx, igy, irad0;
//	double c, x, xmax, xmin, y, ymax, ymin, rad0, dia0;
//	int ixmax=INXY;
//	int iymax=INXY;

//	gcls(); //��ʃN���A
//	xmin=0.; xmax=1.; ymin=0.; ymax=1.;

//	printf("time %f\n",time1);
//	dia0=1.0/nd;
//	rad0=dia0/2.0;
//	irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;

//	for(i=0;i<=nd;i++){
//		for(j=0;j<=nd;j++){
//			x=rad0+dia0*i;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
//			y=rad0+dia0*j;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
//			ii=i; jj=j; if(i==nd){ii=0;} if(j==nd){jj=0;}

//			col_R=s1h[ii][jj];
//			col_G=s2h[ii][jj];
//			col_RG=col_R+col_G;
//			if(col_RG>1.){col_RG=1.;}
//			col_B=1.-col_RG;

//			if(col_R>=0.999){col_R=1.;} if(col_R<=0.001){col_R=0.;}
//			if(col_G>=0.999){col_G=1.;} if(col_G<=0.001){col_G=0.;}
//			if(col_B>=0.999){col_B=1.;} if(col_B<=0.001){col_B=0.;}

//			gcolor((int)(255*col_R),(int)(255*col_G),(int)(255*col_B));
//			grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);
//		}
//	}
//	swapbuffers();
//}

//******* Sin, Cos �̃e�[�u���ݒ� *************************************
void table()
{
	int it, it1, it2, mc, mn;
	double q;

	q=2.0*PI/nd;
	for(it=0;it<=nd2-1;it++){ c[it]=cos(q*it); s[it]=sin(q*it); }
	ik[0]=0; mn=nd2; mc=1;
	for(it1=1;it1<=ig;it1++){
		for(it2=0;it2<=mc-1;it2++){
			ik[it2+mc]=ik[it2]+mn;
		}
		mn=mn/2; mc=2*mc;
	}
}

//********** �P���������t�[���G�ϊ� **************************************
void fft()
{
	int ix, ka, kb, l2, lf, mf, n2, nf;
	double tj, tr;

	l2=1;
	for(lf=1;lf<=ig;lf++){
		n2=nd2/l2;
		for(mf=1;mf<=l2;mf++){
			for(nf=0;nf<=n2-1;nf++){
				ix=nf*l2;
				ka=nf+2*n2*(mf-1);
				kb=ka+n2;
				tr=xrf[ka]-xrf[kb];  					tj=xif[ka]-xif[kb];
				xrf[ka]=xrf[ka]+xrf[kb]; 			xif[ka]=xif[ka]+xif[kb];
				xrf[kb]=tr*c[ix]+tj*qs*s[ix];	xif[kb]=tj*c[ix]-tr*qs*s[ix];
			}
		}
		l2=l2*2;
	}

}

//************ �Q���������t�[���G�ϊ� ***********************************
void rcfft()
{
	int i, ic, ir, j;

	for(ir=0;ir<=ndm;ir++){
		for(ic=0;ic<=ndm;ic++){
			xrf[ic]=xr[ir][ic];	xif[ic]=xi[ir][ic];
		}
	fft();
		for(ic=0;ic<=ndm;ic++){
			xr[ir][ic]=xrf[ik[ic]];	xi[ir][ic]=xif[ik[ic]];
		}
	}
	for(ic=0;ic<=ndm;ic++){
		for(ir=0;ir<=ndm;ir++){
			xrf[ir]=xr[ir][ic];	xif[ir]=xi[ir][ic];
		}
	fft();
		for(ir=0;ir<=ndm;ir++){
			xr[ir][ic]=xrf[ik[ir]];	xi[ir][ic]=xif[ik[ir]];
		}
	}
	if(qs>0.){return;}
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			xr[i][j]=xr[i][j]/nd/nd;	xi[i][j]=xi[i][j]/nd/nd;
		}
	}

}

//************ �f�[�^�̕ۑ� *******************************
void datsave_gnuplot()
{
	FILE	*stream;
	char	fName[256];
	int 	i, j;

	iout = iout + 1;
	sprintf(fName,"mt_%06d.dat",iout);
	stream = fopen(fName, "w");

	fprintf(stream, "%e\n", time1);
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(stream, "%5d  %5d  %e  %e \n", i, j, s1h[i][j], s2h[i][j]);
		}
		fprintf(stream, "\n");
	}
	fclose(stream);
}

void datsave_paraview()
{
	FILE	*fp;
	char	fName[256];
	int 	i, j;
	
	iout = iout + 1;
	sprintf(fName,"mt_result%05d.vtk",iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d 1 \n",(ndm+1),(ndm+1));
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %11d \n",((ndm+1)*(ndm+1)*1));
	fprintf(fp,"SCALARS phase_field_1 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n",s1h[i][j]);
		}
	}
	fprintf(fp,"SCALARS phase_field_2 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n",s2h[i][j]);
		}
	}
	fclose(fp);
}
