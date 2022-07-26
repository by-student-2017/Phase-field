#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#include "gint.fx" 			//��ʏ������֐�
#include "gcls.fx" 			//��ʏ����֐�
#include "rgb.fx"				//RGB�֐�
#include "getscreenx.fx"		//���X�T�C�Y�֐�
#include "getscreeny.fx"		//���Y�T�C�Y�֐�
#include "rectfill.fx"			//�����`�h�׊֐�
#include "line.fx"
#include "pset.fx"
#include "circle.fx"
#include "rectangle.fx"
#include "getkey.fx"			//�L�|���͑҂��֐�

#define DRND(x) ((double)(x)/RAND_MAX*rand()) 	//�����֐��̐ݒ�

#define ND 512								//�P�ӂ̍���������

	int nd=ND, ndm=ND-1; 			//�g�D�̍����������A�g�D�̍���������-1
	double rr=8.3145;					//�K�X�萔
	double c0, temp, time1;			//�����g���A���x�A�v�Z����
	double ch[ND];			//�g�D���̔Z�x�f�|�^�z��

	void ini_comp();					//�����Z�x��ݒ�T�u���|�`��
	void graph_c();					//�Z�x��\���T�u���|�`��
	void datsave();					//�Z�x�f�|�^�ۑ��T�u���|�`��

//******* ���C���v���O���� ******************************************************
main()
{
 	gint(1000,200,"graph1",rgb(1,1,1));		//�`��pWindow�ݒ�

	double ck[ND];					//�g�U�|�e���V����
	double ch2[ND]; 					//�g�D���̔Z�x�f�|�^�\���z��
	double mu_chem, mu_surf;		//�e�|�e���V����
	int  i, ii; 								//����
	double time1max; 							//�ő厞��
	double c;										//�Z�x
	double cddtt;								//�Z�x�̑���
	double kappa_c;								//�Z�x���z�G�l���M�|�萔
	double al;									//�v�Z�̈�
	double b1;									//�����u���b�N�T�C�Y
	double delt;									//���ԍ���
	double Mc;							//�Փ��x�֐��Ƃ��̔���
	double L0;									//���q�ԑ��ݍ�p�p�����[�^
	double c_flu;								//�Z�x��̗h�炬�̑傫��

	double cip, cim; 	//�����u���b�N�ɂ�����c�𒆐S�ɁA���̏㉺���E�̔Z�x
	int   ip, im; 		//�����ii+1, i-1, j+1, j-1�j
	double sumc, dc0; 			// �Z�x��̑��a�A���ϑg������̂���

//****************************************************************************

	printf("DELT(0.01)=  "); scanf(" %lf",&delt);		//���Ԃ����ݓ���
	//delt=0.01; 					//�W���l

	printf("c0(0.3) = "); scanf(" %lf",&c0); 		//�����g������
	//c0=0.3; 					//�W���l

	//printf("Temp(K) = "); scanf(" %lf",&temp); 	//�������x����
	temp=1000.;				//�W���l

	al=500.;					//�v�Z�̈�̈�ӂ̒���(nm)
	al=al*1.0e-9;				//(m)�ɕϊ�
	b1=al/nd;					//�����P�u���b�N�̃T�C�Y
	time1=0.0;					//�X�^�|�g����(���|�v�̉�)
	time1max=1.0+5.0e+05;		//�v�Z�Ő؎���(���|�v�̉�)

	L0=2.5e+04/rr/temp; //���q�ԑ��ݍ�p�p�����|�^�i���łɖ��������ς݁j
	kappa_c=5.0e-15/b1/b1/rr/temp;	//�Z�x���z�G�l���M�|�W���i �Ŗ��������j
	Mc=c0*(1.0-c0);  			//�Փ��x
	c_flu=0.1;								//�Z�x��̗h�炬�W��

//**********************************************************************************
	ini_comp();								//�����Z�x��̐ݒ�

//*** �������̎��Ԕ��W�ߒ��̌v�Z�X�^�|�g *******************************************
start: ;

	if((((int)(time1) % 1000)==0)) {datsave();}	//�Z�x���ۑ�
	if((((int)(time1) % 500)==0)){ graph_c(); }	//�Z�x���\��

//******[�g�U�|�e���V�����̌v�Z]****************************************************
	for(i=0;i<=ndm;i++){
		ip=i+1; im=i-1;  if(i==ndm){ip=0;}  if(i==0){im=ndm;}//�����I���E����
		c=ch[i]; cip=ch[ip]; cim=ch[im];
	 	mu_chem=L0*(1.-2.*c)+log(c)-log(1.-c); 		//���w�|�e���V����
		mu_surf=-2.*kappa_c*(cip+cim-2.*c);	//���z�|�e���V����
		ck[i]=mu_chem+mu_surf; 				//�g�U�|�e���V����
	}

//******[�Z�x��̎��ԕω�]**************************************
	for(i=0;i<=ndm;i++){
		ip=i+1; im=i-1;  if(i==ndm){ip=0;} 	if(i==0){im=ndm;}//�����I���E����
		cddtt=Mc*(ck[ip]+ck[im]-2.0*ck[i]); 							//����`�g�U������
		ch2[i]=ch[i]+(cddtt+c_flu*(2.0*DRND(1.0)-1.0))*delt;
	}

//*** [�Z�x��̎��x�̕␳] *******************************************************
//*** ���l�v�Z�ł���̂ŁA���� �Z�x��̎��x�̕␳���s���K�v������B] *************
  sumc=0.; for(i=0;i<=ndm;i++){ sumc=sumc+ch2[i]; }
	dc0=sumc/nd-c0;
	for(i=0;i<=ndm;i++){
		ch[i]=ch2[i]-dc0; if(ch[i]>=1.0){ch[i]=0.99999;}  if(ch[i]<=0.0){ch[i]=0.00001;}
	}

//******[���ԑ���]*************************************************
	time1=time1+1.0;  if(time1<time1max){goto start;}

end:;
  return EXIT_SUCCESS;
}

//************[�����Z�x�g]*****************************************
void ini_comp()
{
	int i, j;
  srand(time(NULL)); // ����������
	for(i=0;i<=ndm;i++){ ch[i]=c0+0.001*(2.0*DRND(1)-1.0); }
}

//*******[������g�D�̕`��]**************************************************
void graph_c()
{
	int i, ii, i1, ii1, i2, ii2;
	double col;
	double c, x, xmax, xmin, y, ymax, ymin, rad0, dx, dy;
	double gx1, gy1, gx2, gy2;
	int ixmin=0, iymin=0, igx1, igy1, igx2, igy2, irad0;
	int ixmax=getscreenx()-1; // �ő�x���W  Maximum-X
	int iymax=getscreeny()-1; // �ő�y���W  Maximum-Y
	int idx, idy;

	xmin=0.0; xmax=1.0; dx=0.1;
	ymin=0.0; ymax=1.0; dy=0.1;
	idx=ixmax*(dx/(xmax-xmin))+0.5;
	idy=iymax*(dy/(ymax-ymin))+0.5;

	gcls(); //��ʃN���A

	printf("time %f\n",time1);
	rad0=1.0/nd/2.0;
	irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;

	for(i=0;i<=nd-1;i++){
		i1=i; i2=i+1;
		gx1=1.0/nd*i1+rad0;  gx2=1.0/nd*i2+rad0;
		ii1=i1; ii2=i2; if(ii2==nd){ii2=0;}
		gy1=ch[ii1];        gy2=ch[ii2];
		igx1=(ixmax-ixmin)/(xmax-xmin)*(gx1-xmin)+ixmin;
		igy1=iymin+iymax-((iymax-iymin)/(ymax-ymin)*(gy1-ymin)+iymin);
		igx2=(ixmax-ixmin)/(xmax-xmin)*(gx2-xmin)+ixmin;
		igy2=iymin+iymax-((iymax-iymin)/(ymax-ymin)*(gy2-ymin)+iymin);
		line(igx1, igy1, igx2, igy2, 2);			//circle(igx, igy, 1, 2);
	}

	rectangle(ixmin, iymin, ixmax, iymax, 0);
	for(i=0;i<ixmax;i+=idx){line(i, iymin, i, iymax, 0);}
	for(i=0;i<iymax;i+=idy){line(ixmin, i, ixmax, i, 0);}

}

//************[�Z�x�g�f�|�^�̕ۑ�]************************************
void datsave()
{
	FILE		*stream;
	int 		i;

	stream = fopen("test.dat", "a");		//�ۑ��t�@�C������test.dat�Ƃ��Ă���B
	fprintf(stream, "%e \n", time1);		//���Ԃ̏�������
	for(i=0;i<=ndm;i++){ fprintf(stream, "%e  ", ch[i]); }		//�Z�x��̏�������
	fprintf(stream, "\n");
	fclose(stream);
}

