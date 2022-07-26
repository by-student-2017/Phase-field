#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#include "wingxa.h"

#define DRND(x) ((double)(x)/RAND_MAX*rand()) 	//�����֐��̐ݒ�

#define ND 80								//�P�ӂ̍���������
#define INXY 400						//�`��window�P�ӂ̃s�N�Z���T�C�Y

	int nd=ND, ndm=ND-1; 			//�g�D�̍����������A�g�D�̍���������-1
	double rr=8.3145;					//�K�X�萔
	double ca, temp, time1;			//�����g���A���x�A�v�Z����
	double ch[ND][ND];			//�g�D���̔Z�x�f�|�^�z��

	void shokiha_c();					//�����Z�x��ݒ�T�u���|�`��
	void graph_c();					//�Z�x��\���T�u���|�`��
	void datsave();					//�Z�x�f�|�^�ۑ��T�u���|�`��

//******* ���C���v���O���� ******************************************************
int main(void){

	double ck[ND][ND];					//�g�U�|�e���V����
	double ch2[ND][ND]; 					//�g�D���̔Z�x�f�|�^�\���z��
	double mu_chem, mu_str, mu_surf;		//�e�|�e���V����
	int  i, j, k, l; 								//����
	double time1max; 							//�ő厞��
	double amob_c;								//���q�̈Փ��x�萔
	double c;										//�Z�x
	double cddtt;								//�Z�x�̑���
	double kapa_c;								//�Z�x���z�G�l���M�|�萔
	double al;									//�v�Z�̈�
	double b1;									//�����u���b�N�T�C�Y
	double delt;									//���ԍ���
	double a0;									//Fe�̊i�q�萔
	double vm0;									//�����̐�
	double c11a, c12a, c44a;					//Fe�̒e����
	double c11b, c12b, c44b;					//Cr�̒e����
	double c11, c12, c44;						//�����̕��ϒe����
	double y100;								//�e���֐�
	double eta;									//�i�q�~�X�}�b�`
	double Da, Db, Dab;						//���Ȋg�U�W���Ƃ��̔�
	double Mc, dMc;							//�Փ��x�֐��Ƃ��̔���
	double L0;									//���q�ԑ��ݍ�p�p�����[�^

	double cip, cim, cjp, cjm; 	//�����u���b�N�ɂ�����c�𒆐S�ɁA���̏㉺���E�̔Z�x
	double ck1dev, ck2dev;  	//�g�U�|�e���V�����̂P�K�����ƂQ�K����
	double ck0, ckip, ckim, ckjp, ckjm; //�����u���b�N�ɂ�ck0�𒆐S�ɂ��̏㉺���E�̊g�U�|�e���V����
	int   ip, im, jp, jm; 		//�����ii+1, i-1, j+1, j-1�j
	double sumc, dca; 			// �Z�x��̑��a�A���ϑg������̂���

//****************************************************************************

	printf("DELT(0.2)=  "); scanf(" %lf",&delt);		//���Ԃ����ݓ���
//	delt=0.2; 					//�W���l

	printf("ca(0.4) = "); scanf(" %lf",&ca); 		//�����g������
//	ca=0.4; 					//�W���l

	printf("Temp(K)(673.) = "); scanf(" %lf",&temp); 	//�������x����
//	temp=673.;				//�W���l

	al=150.;					//�v�Z�̈�̈�ӂ̒���(nm)
	al=al*1.0e-9;				//(m)�ɕϊ�
	b1=al/nd;					//�����P�u���b�N�̃T�C�Y
	amob_c=1.;				//�Փ��x�i�P�ɋK�i���j
	time1=0.;					//�X�^�|�g����(���|�v�̉�)
	time1max=100001.;		//�v�Z�Ő؎���(���|�v�̉�)

	a0=2.8664E-10;  				//Fe�̊i�q�萔
	vm0=6.02E23*a0*a0*a0/2.; 	//�����̐�
	L0=(21020.8-9.31889*temp)/rr/temp; //���q�ԑ��ݍ�p�p�����|�^�i���łɖ��������ς݁j[��(2-2)]
	kapa_c=6.0e-15/b1/b1/rr/temp;	//�Z�x���z�G�l���M�|�W���i �Ŗ��������j
	eta=0.00614;						//�i�q�~�X�}�b�`

	c11a=2.331e11*vm0/rr/temp; 		//Fe�̒e�����A �Ŗ�������
	c12a=1.3544e11*vm0/rr/temp;
	c44a=1.1783e11*vm0/rr/temp;
	c11b=3.5e11*vm0/rr/temp; 		//Cr�̒e����
	c12b=0.678e11*vm0/rr/temp;
	c44b=1.008e11*vm0/rr/temp;

	c11=(1.-ca)*c11a+ca*c11b;			//�����̒e����
	c12=(1.-ca)*c12a+ca*c12b;
	c44=(1.-ca)*c44a+ca*c44b;

	y100=c11+c12-2.*(c12*c12/c11);		//�e���֐�Y<100>

	Da=1.0e-4*exp(-294000./rr/temp);		//Fe�̎��Ȋg�U�W��
	Db=2.0e-5*exp(-308000./rr/temp);		//Cr�̎��Ȋg�U�W��
	Dab=Db/Da;								//Fe�̎��Ȋg�U�W���ŋK�i��

//**********************************************************************************
	shokiha_c();								//�����Z�x��̐ݒ�
 	gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);//�`��Window�\��

//*** �������̎��Ԕ��W�ߒ��̌v�Z�X�^�|�g *******************************************
start: ;

//	if((((int)(time1) % 200)==0)) {datsave();} 	//200���|�v���ɔZ�x��f�|�^���f�B�X�N�ɕۑ�
//  ���݂��̍s�̑O��//�����蕶�S�̂��R�����g���ɂ��Ă���B//���O���΃f�|�^�̕ۑ�
//  ���s���邪�A���l�f�|�^�͔��ɑ傫���Ȃ�ꍇ������̂ŁA//���O���ꍇ�ɂ́A
//  �n�|�h�f�B�X�N�̋󂫗e�ʂ�f�|�^�ۑ��ʓ���ǂ��m�F������A//���O���ꂽ���B

	if((((int)(time1) % 100)==0)) {graph_c();}	//20���|�v���ɔZ�x����r�b�g�}�b�v�摜�Ƃ��ĉ�ʂɕ\��
//******[�g�U�|�e���V�����̌v�Z]****************************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm) {ip=0;}  //�����I���E����
			if(i==0) {im=ndm;}
			if(j==ndm) {jp=0;}
			if(j==0) {jm=ndm;}

			c=ch[i][j];
			cip=ch[ip][j]; cim=ch[im][j];
			cjp=ch[i][jp]; cjm=ch[i][jm];

		 	mu_chem=L0*(1.-2.*c)+log(c)-log(1.-c); 		//���w�g�U�|�e���V����[��(2-3)]
			mu_surf=-2.*kapa_c*(cip+cim+cjp+cjm-4.*c);	//���z�g�U�|�e���V����[��(2-5)]
			mu_str=2.*eta*eta*y100*(c-ca); 					//�e���g�U�|�e���V����[��(2-8)]

			ck[i][j]=mu_chem+mu_str+mu_surf; 				//�g�U�|�e���V����

		}
	}

//******[�Z�x��̎��ԕω�]**************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm) {ip=0;} 	//�����I���E����
			if(i==0) {im=ndm;}
			if(j==ndm) {jp=0;}
			if(j==0) {jm=ndm;}

			c=ch[i][j];
			cip=ch[ip][j]; cim=ch[im][j];
			cjp=ch[i][jp]; cjm=ch[i][jm];

			ck0=ck[i][j];
			ckip=ck[ip][j]; ckim=ck[im][j];
			ckjp=ck[i][jp]; ckjm=ck[i][jm];

			ck2dev=(ckip+ckim+ckjp+ckjm)-4.*ck0; 							//�g�U�|�e���V�����̂Q�K����
			Mc=(ca+Dab*(1.-ca))*ca*(1.-ca);  										//�Փ��x�֐�
			cddtt=amob_c*Mc*ck2dev; 							//����`�g�U������
			ch2[i][j]=ch[i][j]+cddtt*delt; 											//�Z�x��̎��Ԕ��W

		}
	}

//*** [�Z�x��̎��x�̕␳] *******************************************************
//*** ���l�v�Z�ł���̂ŁA���� �Z�x��̎��x�̕␳���s���K�v������B] *************
  sumc=0.;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			sumc+=ch2[i][j];
		}
	}
	dca=sumc/nd/nd-ca;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ch[i][j]=ch2[i][j]-dca;
			if(ch[i][j]>=1.){ch[i][j]=0.9999;}
			if(ch[i][j]<=0.){ch[i][j]=0.0001;}
		}
	}

//******[���ԑ���]*************************************************
	if(keypress()){return 0;}
	time1=time1+1.;
	if (time1<time1max) {goto start;}

end:;
  return 0;
}

//************[�����Z�x�g]*****************************************
void shokiha_c()
{
	int i, j;
	double rnd0; 
  srand(time(NULL)); // ����������

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			rnd0=2.*DRND(1)-1.;
			ch[i][j]=ca+rnd0/50.;
		}
	}
}

//*******[������g�D�̕`��]**************************************************
void graph_c()
{
	int i, j, ii, jj;
	double col;
	double c, x, xmax, xmin, y, ymax, ymin, rad0;
	int ixmin=0, iymin=0, igx, igy, irad0;
	int ixmax=INXY;
	int iymax=INXY;

	gcls(); //��ʃN���A
	xmin=0.; xmax=1.;
	ymin=0.; ymax=1.;

	printf("time %f\n",time1);
	rad0=1./nd/2.;
	irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;

	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			x=1./nd*i+rad0;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
			y=1./nd*j+rad0;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
			ii=i; jj=j;  if(i==nd){ii=0;}  if(j==nd){jj=0;}
			col=1.-ch[ii][jj];
			if(col>=1.0){col=1.0;}  if(col<=0.0){col=0.0;}
			gcolor((int)(255*col),(int)(255*col),(int)(255*col));
			grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);
		}
	}
	swapbuffers();

} 

//************[�Z�x�g�f�|�^�̕ۑ�]************************************
void datsave()
{
	FILE		*stream;
	int 		i, j;

	stream = fopen("test.dat", "a");		//�ۑ��t�@�C������test.dat�Ƃ��Ă���B
	fprintf(stream, "%e\n", time1);		//���Ԃ̏�������
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fprintf(stream, "%e  ", ch[i][j]);		//�Z�x��̏�������
		}
	}
	fprintf(stream, "\n");
	fclose(stream);
}
//*********************************************************************

