#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())	//�����̐ݒ�

//#define ND 128					//�����v�Z�ɂ�����v�Z�̈��ӂ̕�����

	//int nd=ND;					//�v�Z�̈�̈�ӂ̍���������(�����u���b�N��)
	//int ndm=ND-1;				//ND-1���`
	double PI=3.141592654, rr=8.3145, temp, time1;//��,�K�X�萔,��Ή��x,�v�Z�J�E���g��(���Ԃɔ��)
	double c2a;					//���ϑg��(B����)
	//double c2h[ND][ND];			//�Ǐ��g��
	//double ph[ND][ND];			//�Ǐ��̐ϕ���
	int Nstep, iout;

	void ini000(double *c2h, int ND);				//�����Z�x�v���t�@�C���̐ݒ�T�u���[�`��
	void datsave(double *ph, double *c2h, int ND);				//�f�[�^�ۑ��T�u���[�`��
	void datsave_paraview(double *ph, double *c2h, int ND);	//�f�[�^�ۑ��T�u���[�`��

//******* ���C���v���O���� ******************************************
int main(void)
{
	int ND, nd, ndm;
	
	double c1, c2;				//�Ǐ��Z�x
	double c2max, c2min;		//�ő�Z�x�ƍŏ��Z�x
	//double c2h2[ND][ND];		//�Ǐ��Z�x��̕⏕�s��
	//double c2k_chem, c2k_su, c2k[ND][ND];	//�Ǐ��|�e���V����
	double c2k_chem, c2k_su;	//�Ǐ��|�e���V����
	double dc2a, sumc2;			//�Z�x��̎��x�v�Z�Ɏg�p���Ă���ϐ�
	double dakd2;				//�g�U�|�e���V�����̓�K����
	double c2ddtt;				//�Z�x��̎��ԕϓ���

	int   i, j;					//����
	int   ip, im, jp, jm;		//(i+1),(i-1),(j+1),(j-1)
	double al, b1, rtemp, delt;	//�v�Z�̈��ӂ̒����A�����v���b�N�P�ӂ̒����ART�A���Ԃ�����
	double time1max;			//�v�Z�J�E���g���̍ő�l�i�v�Z�I���J�E���g�j
	double cmob22;				//�Փ��x

	double om_12, om_12e;		//���ݍ�p�p�����[�^
	double kapa_c2, kapa_c2c;	//�Z�x���z�G�l���M�[�W��

//****** �v�Z��������ѕ����萔�̐ݒ� ****************************************
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
	ND      = int(data[0]);
	c2a     = data[1];
	delt    = data[2];	// timestep
	temp    = data[3];	// [K]
	al      = data[4];	// [nm]
	cmob22  = data[5];
	om_12e  = data[6];
	kapa_c2c= data[7];
	time1max= int(data[8]);
	Nstep   = int(data[9]);
	printf("---------------------------------\n");
	//
	nd=ND;
	ndm=ND-1;
	//
	//double ch[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	//a[z * ySize * xSize + y * xSize + x], a[i][j][k]=a[i*n*m + j*n + k]
	double *c2h  = (double *)malloc(sizeof(double)*( ND*ND ));	//�Ǐ��g��
	double *ph   = (double *)malloc(sizeof(double)*( ND*ND ));	//�Ǐ��̐ϕ���
	double *c2h2 = (double *)malloc(sizeof(double)*( ND*ND ));	//�Ǐ��Z�x��̕⏕�s��
	double *c2k  = (double *)malloc(sizeof(double)*( ND*ND ));	//�Ǐ��|�e���V����

	//printf("C2B(0.45) =  ");	scanf(" %lf",&c2a);//�W�����o�͂��畽�ϑg��(B����)�����	//	c2a=0.45;
	//printf("delt(0.005)=  ");	scanf(" %lf",&delt);//���ԍ���	//	delt=0.005;

	//temp=900.0;					//���x�iK�j
	rtemp=rr*temp;				//RT
	//al=100.0*1.0E-09;			//�v�Z�̈�̈�ӂ̒���(nm)
	al=al*1.0E-09;			//�v�Z�̈�̈�ӂ̒���(nm)
	b1=al/(double)ND;			//�����v���b�N�P�ӂ̒���

	//cmob22=1.0;					//�Փ��x

	//om_12=25000./rtemp; 		//���ݍ�p�p�����[�^(J/mol�ŁART�Ŗ�������)
	om_12=om_12e/rtemp; 		//���ݍ�p�p�����[�^(J/mol�ŁART�Ŗ�������)

	//kapa_c2=5.0e-15/b1/b1/rtemp;//�Z�x���z�G�l���M�[�W��(Jm^2/mol�ŁART��b1^2�Ŗ�������)
	kapa_c2=kapa_c2c/b1/b1/rtemp;//�Z�x���z�G�l���M�[�W��(Jm^2/mol�ŁART��b1^2�Ŗ�������)

	time1=0.0;
	//time1max=1.0e05+1.0;//�v�Z�J�E���g���̏����l�ƍő�l

//*** �����Z�x��̐ݒ�ƕ`��Window�\�� *****************************************

	ini000(c2h, ND);					//�����Z�x��̐ݒ�

//**** �V�~�����[�V�����X�^�[�g ******************************
//Nstep = 1000;
iout = -1;
start: ;

	if((((int)(time1) % Nstep)==0)) {datsave(ph, c2h, ND);} //���J�Ԃ��J�E���g���ɔZ�x���ۑ�
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(ph, c2h, ND);} //���J�Ԃ��J�E���g���ɔZ�x���ۑ�

//***** �|�e���V������̌v�Z(5.2�ߎQ��) ***********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//�����I���E����
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}

			//c2=c2h[i][j];  c1=1.0-c2;//�Ǐ��Z�x��
			c2=c2h[i*ND+j];  c1=1.0-c2;//�Ǐ��Z�x��

			c2k_chem=om_12*(c1-c2)+(log(c2)-log(c1));//���w�g�U�|�e���V����
			//c2k_su=-2.*kapa_c2*(c2h[ip][j]+c2h[im][j]+c2h[i][jp]+c2h[i][jm]-4.0*c2);//���z�|�e���V����
			c2k_su=-2.*kapa_c2*(c2h[ip*ND+j]+c2h[im*ND+j]+c2h[i*ND+jp]+c2h[i*ND+jm]-4.0*c2);//���z�|�e���V����

			//c2k[i][j]=c2k_chem+c2k_su;//�g�U�|�e���V����
			c2k[i*ND+j]=c2k_chem+c2k_su;//�g�U�|�e���V����
		}
	}

//***** ���W�������̌v�Z(5.2�ߎQ��) **********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;} if(i==0){im=ndm;}//�����I���E����
			if(j==ndm){jp=0;} if(j==0){jm=ndm;}

			//dakd2=c2k[ip][j]+c2k[im][j]+c2k[i][jp]+c2k[i][jm]-4.0*c2k[i][j];//�g�U�|�e���V�����̓�K����
			dakd2=c2k[ip*ND+j]+c2k[im*ND+j]+c2k[i*ND+jp]+c2k[i*ND+jm]-4.0*c2k[i*ND+j];//�g�U�|�e���V�����̓�K����

			c2ddtt=cmob22*dakd2;//�g�U������

			//c2h2[i][j]=c2h[i][j]+c2ddtt*delt;//�Z�x��̎��Ԕ��W
			c2h2[i*ND+j]=c2h[i*ND+j]+c2ddtt*delt;//�Z�x��̎��Ԕ��W

			//if(c2h[i][j]>=1.0){c2h[i][j]=1.0-1.0e-06;}//�Z�x��̕␳
			//if(c2h[i][j]<=0.0){c2h[i][j]=1.0e-06;}
			if(c2h[i*ND+j]>=1.0){c2h[i*ND+j]=1.0-1.0e-06;}//�Z�x��̕␳
			if(c2h[i*ND+j]<=0.0){c2h[i*ND+j]=1.0e-06;}
		}
	}

//*** �Z�x��̎��x�␳ ***********************************************
	//sumc2=0.; for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ sumc2+=c2h2[i][j]; } }
	sumc2=0.; for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ sumc2+=c2h2[i*ND+j]; } }
	dc2a=sumc2/ND/ND-c2a;
	//for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ c2h[i][j]=c2h2[i][j]-dc2a; } }
	for(i=0;i<=ndm;i++){  for(j=0;j<=ndm;j++){ c2h[i*ND+j]=c2h2[i*ND+j]-dc2a; } }

//*** �t�F�[�Y�t�B�[���h�̐ݒ�(5.3.1���Q��) ***************************************
	c2max=0.0; c2min=1.0;
	for(i=0;i<=ndm;i++){  
		for(j=0;j<=ndm;j++){ 
			//if(c2max<=c2h[i][j]){c2max=c2h[i][j];}
			//if(c2min>=c2h[i][j]){c2min=c2h[i][j];}
			if(c2max<=c2h[i*ND+j]){c2max=c2h[i*ND+j];}
			if(c2min>=c2h[i*ND+j]){c2min=c2h[i*ND+j];}
		} 
	}

	for(i=0;i<=ndm;i++){  
		for(j=0;j<=ndm;j++){ 
			//ph[i][j]=(c2h[i][j]-c2min)/(c2max-c2min);
			ph[i*ND+j]=(c2h[i*ND+j]-c2min)/(c2max-c2min);
		} 
	}

//*********************************************************************
	//if(keypress()){return 0;}	//�L�[�҂����
	time1=time1+1.0;  if(time1<time1max){goto start;}//�ő�J�E���g���ɓ��B�������ǂ����̔��f

	end:;
  return 0;
}

//************ �����Z�x��̐ݒ�T�u���[�`�� *************
void ini000(double *c2h, int ND)
{
	int i, j, id;
 	//srand(time(NULL));//�����̎�ݒ�
	int ndm=ND-1;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//c2h[i][j]=c2a+0.01*(2.0*DRND(1)-1.0);//�Z�x����ő�}1%�̗����ɂĐݒ�
			c2h[i*ND+j]=c2a+0.01*(2.0*DRND(1)-1.0);
		}
	}

}

//************ �f�[�^�ۑ��T�u���[�`�� *******************************
void datsave(double *ph, double *c2h, int ND)
{
	FILE *stream;	//�X�g���[���̃|�C���^�ݒ�
	int i, j;			//����
	int ndm=ND-1;

	stream = fopen("ph.dat", "a");	//�������ސ�̃t�@�C����ǋL�����ŃI�[�v��
	fprintf(stream, "%f\n", time1);		//�J�Ԃ��J�E���g���̕ۑ�
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//fprintf(stream, "%e  ", ph[i][j]);//�t�F�[�Y�t�B�[���h�̕ۑ�
			fprintf(stream, "%e  ", ph[i*ND+j]);
			//fprintf(stream, "%e  ", c2h[i][j]);//�Ǐ��Z�x��̕ۑ�
		}
	}
	fprintf(stream, "\n");	//���s�̏�������
	fclose(stream);					//�t�@�C�����N���[�Y
}

void datsave_paraview(double *ph, double *c2h, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j, k;
	int ndm=ND-1;
	
	iout = iout + 1;
	printf("pf_result%06d.vtk \n",iout);
	sprintf(fName,"pf_result%06d.vtk",iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),1);
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %11d \n",((ndm+1)*(ndm+1)*1));
	fprintf(fp,"SCALARS Phase_field float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(fp, "%e  ", ph[i][j]);//�t�F�[�Y�t�B�[���h�̕ۑ�
			//fprintf(fp,"%10.6f\n", ph[i][j]);
			fprintf(fp,"%10.6f\n", ph[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS concentration float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(fp, "%e  ", c2h[i][j]);//�Ǐ��Z�x��̕ۑ�
			//fprintf(fp,"%10.6f\n", c2h[i][j]);
			fprintf(fp,"%10.6f\n", c2h[i*ND+j]);
		}
	}
	fclose(fp);
}