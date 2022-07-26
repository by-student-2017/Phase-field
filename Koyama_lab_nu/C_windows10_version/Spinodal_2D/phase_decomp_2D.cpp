#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#include "wingxa.h"

#define DRND(x) ((double)(x)/RAND_MAX*rand()) 	//�����֐��̐ݒ�
#define ND 64								//�P�ӂ̍���������
#define INXY 400						//�`��window�P�ӂ̃s�N�Z���T�C�Y

	int nd=ND, ndm=ND-1; 			//�g�D�̍����������A�g�D�̍���������-1
	double rr=8.3145;					//�K�X�萔
	double c0, temp, time1;		//�����g���A���x�A�v�Z����
	double ch[ND][ND];				//�g�D���̔Z�x�f�|�^�z��

	void ini_comp();					//�����Z�x��ݒ�T�u���|�`��
	void graph_c();						//�Z�x��\���T�u���|�`��
	void datsave();						//�Z�x�f�|�^�ۑ��T�u���|�`��

int main(void)
{
	int  i, j, k, l; 						//����
	int ip, im, jp, jm; 				//�����ii+1, i-1, j+1, j-1�j
	double ck[ND][ND];					//�g�U�|�e���V����
	double ch2[ND][ND]; 				//�g�D���̔Z�x�f�|�^�\���z��
	double mu_chem, mu_surf;		//�e�|�e���V����
	double time1max; 						//�ő厞��
	double amob_c;							//���q�̈Փ��x�萔
	double c;										//�Z�x
	double cddtt;								//�Z�x�̑���
	double kappa_c;							//�Z�x���z�G�l���M�|�萔
	double al;									//�v�Z�̈�
	double b1;									//�����u���b�N�T�C�Y
	double delt;								//���ԍ���
	double a0;									//�i�q�萔
	double Dab;									//���Ȋg�U�W���Ƃ��̔�
	double Mc;									//�Փ��x�֐��Ƃ��̔���
	double L0;									//���q�ԑ��ݍ�p�p�����[�^
	double c_flu;								//�Z�x��̗h�炬�̑傫��
	double cip, cim, cjp, cjm; 	//�����u���b�N�ɂ�����c�𒆐S�ɁA���̏㉺���E�̔Z�x
	double sumc, dc0; 					//�Z�x��̑��a�A���ϑg������̂���

//****************************************************************************

	printf("delt(0.01)=  "); scanf(" %lf",&delt);		//���Ԃ����ݓ���
	//delt=0.01; 					//�W���l

	//printf("Temp(K) = "); scanf(" %lf",&temp); 		//�������x����
	temp=1000.;						//�W���l

	printf("c0 = ");	scanf(" %lf",&c0); 						//�����g������
	//c0=0.3; 						//�W���l
	//1000K�ɂ�����X�s�m�[�_���g���͂��悻0.21	
	//1000K�ɂ�����o�C�m�[�_���g���͂��悻0.07


	al=60.0;						//�v�Z�̈�̈�ӂ̒���(nm)
	al=al*1.0e-9;				//(m)�ɕϊ�
	b1=al/nd;						//�����P�u���b�N�̃T�C�Y
	time1=0.0;					//�X�^�|�g����(���|�v�̉�)
	time1max=1.0e+06;		//�v�Z�Ő؎���(���|�v�̉�)

	L0=25000./rr/temp; //���q�ԑ��ݍ�p�p�����|�^�i���łɖ��������ς݁j[��(2-2)]
	kappa_c=5.0e-15/b1/b1/rr/temp;	//�Z�x���z�G�l���M�|�W���i �Ŗ��������j
	Mc=c0*(1.0-c0);  								//�Փ��x
	c_flu=0.1;											//�Z�x��̗h�炬�W��

//**********************************************************************************

 	gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);//�`��Window�\��
	ini_comp();										//�����Z�x��̐ݒ�

//*** �������̎��Ԕ��W�ߒ��̌v�Z�X�^�|�g *******************************************
start: ;
	//if((((int)(time1) % 200)==0)) {datsave();}	//200���|�v���ɔZ�x���ۑ�
	if((((int)(time1) % 200)==0)) {graph_c();}		//200���|�v���ɔZ�x���`��

//******[�g�U�|�e���V�����̌v�Z]****************************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//�����I���E����
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
			c=ch[i][j]; cip=ch[ip][j]; cim=ch[im][j]; cjp=ch[i][jp]; cjm=ch[i][jm];

		 	mu_chem=L0*(1.0-2.0*c)+log(c)-log(1.-c); 			//���w�|�e���V����
			mu_surf=-2.0*kappa_c*(cip+cim+cjp+cjm-4.0*c);	//���z�|�e���V����
			ck[i][j]=mu_chem+mu_surf; 										//�g�U�|�e���V����
		}
	}

//******[�Z�x��̎��ԕω�]**************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}	//�����I���E����
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
			cddtt=Mc*(ck[ip][j]+ck[im][j]+ck[i][jp]+ck[i][jm]-4.0*ck[i][j]); 	//����`�g�U������
			ch2[i][j]=ch[i][j]+(cddtt+c_flu*(2.0*DRND(1.0)-1.0))*delt; 				//�Z�x��̎��Ԕ��W
		}
	}

//*** [�Z�x��̎��x�̕␳] *******************************************************
//*** ���l�v�Z�ł���̂ŁA���� �Z�x��̎��x�̕␳���s���K�v������B] *************
  sumc=0.;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){ sumc+=ch2[i][j]; }
	}
	dc0=sumc/nd/nd-c0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ch[i][j]=ch2[i][j]-dc0;
			if(ch[i][j]>=1.0){ch[i][j]=0.9999;}  if(ch[i][j]<=0.0){ch[i][j]=0.0001;}
		}
	}

//******[���ԑ���]*************************************************

	if(keypress()){return 0;}

	time1=time1+1.0;  if(time1<time1max){goto start;}

end:;
  return 0;
}

//************[�����Z�x�g]*****************************************
void ini_comp()
{
	int i, j;
  srand(time(NULL)); // ����������

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ch[i][j]=c0+0.01*(2.0*DRND(1)-1.0);
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

	stream = fopen("test.dat", "a");	//�ۑ��t�@�C������test.dat�Ƃ��Ă���B
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

