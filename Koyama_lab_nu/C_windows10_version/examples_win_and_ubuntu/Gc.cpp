#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//#include "wingxa.h"

#define NDP 101
#define PI 3.14159
#define RR 8.3145
#define INX 600						//�`��window�P��(x)�̃s�N�Z���T�C�Y
#define INY 400						//�`��window�P�ӂ�(y)�s�N�Z���T�C�Y

	int nd=NDP-1;//�ϐ�s�̕ψ�(0=<s=<1)������
	double sh[NDP], Gch[NDP];//s��Gc�̔z��
	double s2h[NDP], Gc2h[NDP];//s��Gc�̕⏕�z��

	double Gc(double a1, double s);//Gc�̊֐�
	void frame0(double xmin,double xmax,double dx,double ymin,double ymax,double dy);//�O���t�g�`��
	void datsave();//�f�[�^�ۑ��T�u���[�`��
	void datload();//�f�[�^�ǂݏo���T�u���[�`��

//*** ���R�G�l���M�|�`�� ***********************************************************
int main(void){

	int i, i1, i2;
	double a1, b1, c1;
	double s, s_max, Gc_max;

	double xmax, xmin, dx, ymax, ymin, dy;
	int ixmin=0, ixmax=INX; // �`��window�̕�
	int iymin=0, iymax=INY; // �`��window�̍���
	double gx, gy, gx1, gy1, gx2, gy2;
	int igx, igy, igx1, igy1, igx2, igy2;

	xmin=0.0; xmax=1.0 ; dx=0.1;//s�̕`��͈�
	ymin=-1.0; ymax=1.0; dy=0.1;//�f���̕`��͈�

//--- �f�|�^�̌v�Z -------------------------------------------------------
	printf("a1(10.0) =  "); scanf(" %lf",&a1);//a�̓���

	s_max=a1/(2.0*a1+12.0);//Gc�̋ɑ�l��^����s�̒l
	Gc_max=a1*a1*a1*(a1+8.0)/32.0/pow((a1+6.0),3.0);//Gc�̋ɑ�l
	printf("s_max, Gc_max= %f, %e \n", s_max, Gc_max);//�W�����o�͂ɐ��l��\��

	for(i=0;i<=nd;i++){
		s=(double)i/(double)nd; sh[i]=s;  Gch[i]=Gc(a1, s);//�f�[�^(s,Gc)�̌v�Z
	}

//--- �f�|�^�̕`�� -------------------------------------------------
	//frame0(xmin,xmax,dx,ymin,ymax,dy);//window�̕`��

	//for(i=0;i<nd;i++){
	//	i1=i; i2=i+1;			//�X�ׂ̗荇���Q�_���A�����ŘA���I�Ɍ���

	//	//���ۂ̒l
	//	gx1=sh[i1];  gy1=Gch[i1];  gx2=sh[i2];  gy2=Gch[i2];

	//	//�X�N���[�����W�ɕϊ�
	//	igx1=(ixmax-ixmin)/(xmax-xmin)*(gx1-xmin)+ixmin;
	//	igy1=(iymax-iymin)/(ymax-ymin)*(gy1-ymin)+iymin;
	//	igx2=(ixmax-ixmin)/(xmax-xmin)*(gx2-xmin)+ixmin;
	//	igy2=(iymax-iymin)/(ymax-ymin)*(gy2-ymin)+iymin;
	//	glinewidth(3);//���̕�
	//	gcolor(0,0,255); gline(igx1, igy1, igx2, igy2);//�ׂ荇���Q�_���A�����Ō���
	//}

	datsave();//�f�[�^(s,Gc)���n�[�h�f�B�X�N�ɕۑ�

	datload();//�f�[�^(s,Gc)���n�[�h�f�B�X�N����ǂݏo��

	sleep(2000);//2�b�x��

	//for(i=0;i<=nd;i++){
	//	gx=s2h[i]; gy=Gc2h[i];//�n�[�h�f�B�X�N����ǂݏo�����f�[�^(s,Gc)
	//	//printf("s, Gc= %f,  %e \n", gx, gy);
	//	igx=(ixmax-ixmin)/(xmax-xmin)*(gx-xmin)+ixmin;
	//	igy=(iymax-iymin)/(ymax-ymin)*(gy-ymin)+iymin;
	//	glinewidth(2);//���̕�
	//	gcolor(255,0,0); gcircle(igx, igy, 4);//�f�[�^�_�𒆐S�ɏ����ȉ~��`��
	//}

	//if(keypress()){return 0;}
	return 1;
}

//*** ���R�G�l���M�[�֐� *****************************************************************
double Gc(double a1, double s)
{
	double b1, c1, gc;

	b1=3.0*a1+12.0;  c1=2.0*a1+12.0;
	gc=a1*s*s/2.0-b1*s*s*s/3.0+c1*pow(s,4.0)/4.0;
	return(gc);
}

//**** window�̐ݒ� ************************************************************
//void frame0(double xmin,double xmax,double dx,double ymin,double ymax,double dy){
//	int idx, idy, i;
//	int ixmin=0, ixmax=INX; // �`��window�̕�
//	int iymin=0, iymax=INY; // �`��window�̍���

//	gwinsize(ixmax,iymax); ginit(1); gsetorg(0,0); //�`��window�̕\��
//	gcolor(255,255,255); grect(0,0, ixmax,iymax);  //window���𔒓h��

//	idx=ixmax*(dx/(xmax-xmin))+0.5;  idy=iymax*(dy/(ymax-ymin))+0.5;
//	gcolor(0,0,0); grectangle(ixmin, iymin, ixmax, iymax);//window�̊O�g
//	for(i=0;i<ixmax;i+=idx){gcolor(0,0,0); gline(i, iymin, i, iymax);}//window���̌r��
//	for(i=0;i<iymax;i+=idy){gcolor(0,0,0); gline(ixmin, i, ixmax, i);}//window���̌r��
//}


//**** �f�[�^�̕ۑ� ************************************************
void datsave()
{
	FILE		*stream;
	int 		i;

	stream = fopen("test.dat", "w");//�㏑���̏ꍇ
	//stream = fopen("test.dat", "a");//�ǋL�̏ꍇ
	for(i=0;i<=nd;i++){
		fprintf(stream, "%e  %e  \n", sh[i], Gch[i]);//�f�[�^�̕ۑ�
	}
	fprintf(stream, "\n");
	fclose(stream);
}
//**** �f�[�^�̓ǂݏo�� ***********************************************
void datload()
{
	FILE	*datload0;
	int i;

	datload0 = fopen("test.dat", "r");
	for(i=0;i<=nd;i++){
		fscanf(datload0, "%lf  %lf  \n", &s2h[i], &Gc2h[i]);
	}
	fclose(datload0);
}

