#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())  //�����̐ݒ�

//#define ND 100					//�����v�Z�ɂ�����v�Z�̈��ӂ̕�����

	//int nd=ND;					//�v�Z�̈�̈�ӂ̍���������(�����u���b�N��)
	//int ndm=ND-1;				//ND-1���`
	double PI=3.141592654;		//��
	double rr=8.3145;			//�K�X�萔
	double temp;				//��Ή��x
	double time1;				//�v�Z�J�E���g��(���Ԃɔ��)
	double c2a, c3a;			//���ϑg��(1:A����, 2:B����, 3:C����)
	//double c2h[ND][ND], c3h[ND][ND];//�Ǐ��g��
	int Nstep, iout;

	void ini000(double *c2h, double *c3h, int ND);	//�����Z�x�v���t�@�C���̐ݒ�T�u���[�`��
	void datsave(double *c2h, double *c3h, int ND);	//�f�[�^�ۑ��T�u���[�`��
	void datsave_paraview(double *c2h, double *c3h, int ND);//�f�[�^�ۑ��T�u���[�`��

//******* ���C���v���O���� ******************************************
int main(void)
{
	int ND, nd, ndm;
	
	double c1, c2, c3;							//�Ǐ��Z�x
	//double c2h2[ND][ND], c3h2[ND][ND];			//�Ǐ��Z�x��̕⏕�s��
	double c2k_chem, c2k_su;					//�Ǐ��|�e���V����
	//double c2k[ND][ND];						//�Ǐ��|�e���V����
	double c3k_chem, c3k_su;					//�Ǐ��|�e���V����
	//double c3k[ND][ND];						//�Ǐ��|�e���V����
	double dc2a, sumc2, dc3a, sumc3, sumc23;	//�Z�x��̎��x�v�Z�Ɏg�p���Ă���ϐ�
	double dakd2, dakd3;						//�g�U�|�e���V�����̓�K����
	double c2ddtt, c3ddtt;						//�Z�x��̎��ԕϓ���

	int   i, j;									//����
	int   ip, im, jp, jm;						//(i+1),(i-1),(j+1),(j-1)
	double al, b1, rtemp, delt;					//�v�Z�̈��ӂ̒����A�����v���b�N�P�ӂ̒����ART�A���Ԃ�����
	double time1max;							//�v�Z�J�E���g���̍ő�l�i�v�Z�I���J�E���g�j
	double cmob22, cmob33, cmob23, cmob32;		//�Փ��x

	double om_12, om_23, om_13;					//���ݍ�p�p�����[�^
	double om_12e, om_23e, om_13e;				//���ݍ�p�p�����[�^
	double kapa_c2, kapa_c3;					//�Z�x���z�G�l���M�[�W��
	double kapa_c2c, kapa_c3c;					//�Z�x���z�G�l���M�[�W��

//****** �v�Z��������ѕ����萔�̐ݒ� ****************************************
	printf("---------------------------------\n");
	printf("read parameters from parameters.txt\n");
	FILE *fp;
	char name[40], comment[72];
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
	ND      = int(data[0]);
	c2a     = data[1];
	c3a     = data[2];
	delt    = data[3];
	temp    = data[4];
	al      = data[5];	// [nm]
	cmob22  = data[6];
	cmob33  = data[7];
	cmob23  = data[8];
	om_12e  = data[9];
	om_13e  = data[10];
	om_23e  = data[11];
	kapa_c2c= data[12];
	kapa_c3c= data[13];
	time1max= int(data[14]);
	Nstep   = int(data[15]);
	printf("---------------------------------\n");
	//
	nd=ND;					//�v�Z�̈�̈�ӂ̍���������(�����u���b�N��)
	ndm=ND-1;				//ND-1���`
	//
	double *c2h  = (double *)malloc(sizeof(double)*( ND*ND ));	//�Ǐ��g��
	double *c3h  = (double *)malloc(sizeof(double)*( ND*ND ));	//�Ǐ��g��
	double *c2h2 = (double *)malloc(sizeof(double)*( ND*ND ));	//�Ǐ��Z�x��̕⏕�s��
	double *c3h2 = (double *)malloc(sizeof(double)*( ND*ND ));	//�Ǐ��Z�x��̕⏕�s��
	double *c2k  = (double *)malloc(sizeof(double)*( ND*ND ));	//�Ǐ��|�e���V����
	double *c3k  = (double *)malloc(sizeof(double)*( ND*ND ));	//�Ǐ��|�e���V����
	//
	//printf("C2B(0.333) =  ");	scanf(" %lf",&c2a);//�W�����o�͂��畽�ϑg��(B����)�����	//	c2a=1./3.;
	//printf("C3C(0.333) =  ");	scanf(" %lf",&c3a);//�W�����o�͂��畽�ϑg��(C����)�����	//	c3a=1./3.;
	//printf("delt(0.005)=  ");	scanf(" %lf",&delt);//���ԍ���	//	delt=0.005;

	//temp=900.0;					//���x�iK�j
	rtemp=rr*temp;				//RT
	//al=100.0*1.0E-09;			//�v�Z�̈�̈�ӂ̒���(m)
	al=al*1.0E-09;			//�v�Z�̈�̈�ӂ̒���(m)
	b1=al/(double)ND;			//�����v���b�N�P�ӂ̒���

	//cmob22=1.0;					//�Փ��x
	//cmob33=1.0;					//�Փ��x
	//cmob23=cmob32=-0.5;			//�Փ��x
	cmob32=cmob23;

	//om_12=25000./rtemp; 		//���ݍ�p�p�����[�^(J/mol�ŁART�Ŗ�������)
	//om_13=25000./rtemp;
	//om_23=25000./rtemp;
	om_12=om_12e/rtemp; 		//���ݍ�p�p�����[�^(J/mol�ŁART�Ŗ�������)
	om_13=om_13e/rtemp;
	om_23=om_23e/rtemp;

	//kapa_c2=5.0e-15/b1/b1/rtemp;//�Z�x���z�G�l���M�[�W��(Jm^2/mol�ŁART��b1^2�Ŗ�������)
	//kapa_c3=5.0e-15/b1/b1/rtemp;//�Z�x���z�G�l���M�[�W��(Jm^2/mol�ŁART��b1^2�Ŗ�������)
	kapa_c2=kapa_c2c/b1/b1/rtemp;//�Z�x���z�G�l���M�[�W��(Jm^2/mol�ŁART��b1^2�Ŗ�������)
	kapa_c3=kapa_c3c/b1/b1/rtemp;//�Z�x���z�G�l���M�[�W��(Jm^2/mol�ŁART��b1^2�Ŗ�������)

	time1=0.0;					//�v�Z�J�E���g���̏����l
	//time1max=1.0e05+1.0;		//�v�Z�J�E���g���̍ő�l

//*** �����Z�x��̐ݒ�ƕ`��Window�\�� *****************************************

	ini000(c2h, c3h, ND);//�����Z�x��̐ݒ�

//**** �V�~�����[�V�����X�^�[�g ******************************
//Nstep = 2000;
iout = -1;
start: ;

	//if((((int)(time1) % Nstep)==0)) {datsave(c2h, c3h, ND);} //���J�Ԃ��J�E���g���ɔZ�x���ۑ�
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(c2h, c3h, ND);} //���J�Ԃ��J�E���g���ɔZ�x���ۑ�

//***** �|�e���V������̌v�Z ***********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//�����I���E����
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}

			//c2=c2h[i][j]; c3=c3h[i][j]; c1=1.0-c2-c3;//�Ǐ��Z�x��
			c2=c2h[i*ND+j]; c3=c3h[i*ND+j]; c1=1.0-c2-c3;//�Ǐ��Z�x��

			c2k_chem=om_12*(c1-c2)-om_13*c3+om_23*c3+(log(c2)-log(c1));//���w�g�U�|�e���V����
			//c2k_su=-2.*kapa_c2*(c2h[ip][j]+c2h[im][j]+c2h[i][jp]+c2h[i][jm]-4.0*c2)
			//				  -kapa_c3*(c3h[ip][j]+c3h[im][j]+c3h[i][jp]+c3h[i][jm]-4.0*c3);//���z�|�e���V����
			c2k_su=-2.*kapa_c2*(c2h[ip*ND+j]+c2h[im*ND+j]+c2h[i*ND+jp]+c2h[i*ND+jm]-4.0*c2)
					  -kapa_c3*(c3h[ip*ND+j]+c3h[im*ND+j]+c3h[i*ND+jp]+c3h[i*ND+jm]-4.0*c3);//���z�|�e���V����

			c3k_chem=om_13*(c1-c3)-om_12*c2+om_23*c2+(log(c3)-log(c1));//���w�g�U�|�e���V����
			//c3k_su=-2.*kapa_c3*(c3h[ip][j]+c3h[im][j]+c3h[i][jp]+c3h[i][jm]-4.0*c3)
			//				  -kapa_c2*(c2h[ip][j]+c2h[im][j]+c2h[i][jp]+c2h[i][jm]-4.0*c2);//���z�|�e���V����
			c3k_su=-2.*kapa_c3*(c3h[ip*ND+j]+c3h[im*ND+j]+c3h[i*ND+jp]+c3h[i*ND+jm]-4.0*c3)
					  -kapa_c2*(c2h[ip*ND+j]+c2h[im*ND+j]+c2h[i*ND+jp]+c2h[i*ND+jm]-4.0*c2);//���z�|�e���V����

			//c2k[i][j]=c2k_chem+c2k_su;//�g�U�|�e���V����(��(4.1))
			//c3k[i][j]=c3k_chem+c3k_su;
			c2k[i*ND+j]=c2k_chem+c2k_su;//�g�U�|�e���V����(��(4.1))
			c3k[i*ND+j]=c3k_chem+c3k_su;
		}
	}

//***** ���W�������̌v�Z **********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;} if(i==0){im=ndm;}//�����I���E����
			if(j==ndm){jp=0;} if(j==0){jm=ndm;}

			//dakd2=c2k[ip][j]+c2k[im][j]+c2k[i][jp]+c2k[i][jm]-4.0*c2k[i][j];//�g�U�|�e���V�����̓�K����
			//dakd3=c3k[ip][j]+c3k[im][j]+c3k[i][jp]+c3k[i][jm]-4.0*c3k[i][j];
			dakd2=c2k[ip*ND+j]+c2k[im*ND+j]+c2k[i*ND+jp]+c2k[i*ND+jm]-4.0*c2k[i*ND+j];//�g�U�|�e���V�����̓�K����
			dakd3=c3k[ip*ND+j]+c3k[im*ND+j]+c3k[i*ND+jp]+c3k[i*ND+jm]-4.0*c3k[i*ND+j];

			c2ddtt=cmob22*dakd2+cmob23*dakd3;//�g�U������(��(4.2))
			c3ddtt=cmob32*dakd2+cmob33*dakd3;

			//c2h2[i][j]=c2h[i][j]+c2ddtt*delt;//�Z�x��̎��Ԕ��W
			//c3h2[i][j]=c3h[i][j]+c3ddtt*delt;
			c2h2[i*ND+j]=c2h[i*ND+j]+c2ddtt*delt;//�Z�x��̎��Ԕ��W
			c3h2[i*ND+j]=c3h[i*ND+j]+c3ddtt*delt;

			//if(c2h[i][j]>=1.0){c2h[i][j]=1.0-1.0e-06;}//�Z�x��̕ψ�␳
			//if(c2h[i][j]<=0.0){c2h[i][j]=1.0e-06;}
			//if(c3h[i][j]>=1.0){c3h[i][j]=1.0-1.0e-06;}
			//if(c3h[i][j]<=0.0){c3h[i][j]=1.0e-06;}
			if(c2h[i*ND+j]>=1.0){c2h[i*ND+j]=1.0-1.0e-06;}//�Z�x��̕ψ�␳
			if(c2h[i*ND+j]<=0.0){c2h[i*ND+j]=1.0e-06;}
			if(c3h[i*ND+j]>=1.0){c3h[i*ND+j]=1.0-1.0e-06;}
			if(c3h[i*ND+j]<=0.0){c3h[i*ND+j]=1.0e-06;}
		}
	}

//*** �Z�x��̎��x�␳ ***********************************************
	//sumc2=0.; for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ sumc2+=c2h2[i][j]; } }
	//dc2a=sumc2/ND/ND-c2a;
	//for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ c2h[i][j]=c2h2[i][j]-dc2a; } }
	sumc2=0.; for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ sumc2+=c2h2[i*ND+j]; } }
	dc2a=sumc2/ND/ND-c2a;
	for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ c2h[i*ND+j]=c2h2[i*ND+j]-dc2a; } }

	//sumc3=0.; for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ sumc3+=c3h2[i][j]; } }
	//dc3a=sumc3/ND/ND-c3a;
	//for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ c3h[i][j]=c3h2[i][j]-dc3a; } }
	sumc3=0.; for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ sumc3+=c3h2[i*ND+j]; } }
	dc3a=sumc3/ND/ND-c3a;
	for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ c3h[i*ND+j]=c3h2[i*ND+j]-dc3a; } }

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//sumc23=c2h[i][j]+c3h[i][j];
			//if(sumc23>=1.0){
			//	c2h[i][j]=c2h[i][j]/sumc23-1.0e-06;
			//	c3h[i][j]=c3h[i][j]/sumc23-1.0e-06;
			//}
			sumc23=c2h[i*ND+j]+c3h[i*ND+j];
			if(sumc23>=1.0){
				c2h[i*ND+j]=c2h[i*ND+j]/sumc23-1.0e-06;
				c3h[i*ND+j]=c3h[i*ND+j]/sumc23-1.0e-06;
			}
		}
	}
//*********************************************************************

	//if(keypress()){return 0;}	//�L�[�҂����
	time1=time1+1.0;  if(time1<time1max){goto start;}//�ő�J�E���g���ɓ��B�������ǂ����̔��f

	end:;
  return 0;
}


//************ �����Z�x��̐ݒ�T�u���[�`�� *************
void ini000(double *c2h, double *c3h, int ND)
{
	int i, j, id;
 	//srand(time(NULL));//�����̎�ݒ�
	int ndm=ND-1;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//c2h[i][j]=c2a+0.01*(2.0*DRND(1)-1.0);//�Z�x����ő�}1%�̗����ɂĐݒ�
			//c3h[i][j]=c3a+0.01*(2.0*DRND(1)-1.0);
			c2h[i*ND+j]=c2a+0.01*(2.0*DRND(1)-1.0);//�Z�x����ő�}1%�̗����ɂĐݒ�
			c3h[i*ND+j]=c3a+0.01*(2.0*DRND(1)-1.0);
		}
	}

}

//************ �f�[�^�ۑ��T�u���[�`�� *******************************
void datsave(double *c2h, double *c3h, int ND)
{
	FILE *stream;	//�X�g���[���̃|�C���^�ݒ�
	int i, j;			//����
	int ndm=ND-1;

	stream = fopen("test.dat", "a");	//�������ސ�̃t�@�C����ǋL�����ŃI�[�v��
	fprintf(stream, "%f\n", time1);		//�v�Z�J�E���g���̕ۑ�
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//fprintf(stream, "%e  %e  ", c2h[i][j], c3h[i][j]);//�Ǐ��Z�x��̕ۑ�
			fprintf(stream, "%e  %e  ", c2h[i*ND+j], c3h[i*ND+j]);//�Ǐ��Z�x��̕ۑ�
		}
	}
	fprintf(stream, "\n");	//���s�̏�������
	fclose(stream);					//�t�@�C�����N���[�Y
}

void datsave_paraview(double *c2h, double *c3h, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j;
	int ndm=ND-1;
	
	iout = iout + 1;
	printf("sp_result%06d.vtk \n",iout);
	sprintf(fName,"sp_result%06d.vtk",iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),(1));
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*1));
	fprintf(fp,"SCALARS concentration_A float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e  ", c2h[i][j], c3h[i][j]);//�Ǐ��Z�x��̕ۑ�
			//fprintf(fp,"%10.6f\n", (1.0-c2h[i][j]-c3h[i][j]));
			fprintf(fp,"%10.6f\n", (1.0-c2h[i*ND+j]-c3h[i*ND+j]));
		}
	}
	fprintf(fp,"SCALARS concentration_B float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e  ", c2h[i][j], c3h[i][j]);//�Ǐ��Z�x��̕ۑ�
			//fprintf(fp,"%10.6f\n", c2h[i][j]);
			fprintf(fp,"%10.6f\n", c2h[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS concentration_C float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e  ", c2h[i][j], c3h[i][j]);//�Ǐ��Z�x��̕ۑ�
			//fprintf(fp,"%10.6f\n", c3h[i][j]);
			fprintf(fp,"%10.6f\n", c3h[i*ND+j]);
		}
	}
	fclose(fp);
}