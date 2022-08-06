#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//�����̊֐��ݒ�

//#define ND 400						//�����v�Z�ɂ�����v�Z�̈��ӂ̕�����

	//int nd=ND, ndm=ND-1; 			//�v�Z�̈�̈�ӂ̍���������(�����u���b�N��)�AND-1���`
	//int nd2=ND/2;				 	//ND/2���`�g�p
	double PI=3.14159;				//�~����
	double rr=8.3145;				//�K�X�萔
	double time1;					//����
	int iout;

	//double s1h[ND][ND];				//���͈͋K���x��

	void ini000(double *s1h, int ND);	//������̐ݒ�T�u���|�`��
	void datsave(double *s1h, int ND);	//�f�|�^�ۑ��T�u���|�`��
	void datsave_paraview(double *s1h, int ND);//�f�|�^�ۑ��T�u���|�`��

//******* ���C���v���O���� ******************************************
int main(void)
{
	int ND, nd, ndm, nd2;
	
	double s1;						//���͈͋K���x��
	//double s1h2[ND][ND];			//��̕⏕�z��
	double s1k_chem, s1k_su;		//�|�e���V����
	double sum11;					//s1�̋�Ԑϕ�
	double s1ddtt;					//s1�̎��ԕω��ʁi���W�������̍��Ӂj

	int   i, j, k, l, ii, jj, kk, iii, jjj;	//����
	int   p, q, m, n;						//����
	int   ip, im, jp, jm, Nstep;			//����
	double al, temp, delt;					//�v�Z�̈�A���ԁA���Ԃ�����
	double time1max;						//�ő厞�ԁi�v�Z���~�߂�ۂɎg�p�j
	double b1, vm0, atom_n;					//�K�i�������A�����̐ρA�P�ʖE���̌��q��
	double smob;							//�����ϑԂ̊ɘa�W��

	double AA0, AA1, AA2, AA3;				//�M�Y�u�G�l���M�[���̌W��
	double AA0e;
	double a1_c, b1_c, c1_c;				//�i�q�萔
	double kappa_s1;						//���z�G�l���M�|�W��
	double kappa_s1c;
	//double ds_fac;							//�����ϑԂ̗h�炬�W��

//****** �v�Z��������ѕ����萔�̐ݒ� ****************************************
	printf("---------------------------------\n");
	printf("read parameters form parameters.txt\n");
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
	delt    = data[1];
	temp    = data[2];
	al      = data[3];	// [nm]
	smob    = data[4];
	AA0e    = data[5];
	kappa_s1c=data[6];
	vm0     = data[7];
	time1max= int(data[8]);
	Nstep   = int(data[9]);
	printf("---------------------------------\n");
	
	//
	nd=ND, ndm=ND-1; 			//�v�Z�̈�̈�ӂ̍���������(�����u���b�N��)�AND-1���`
	nd2=ND/2;				 	//ND/2���`�g�p
	//
	double *s1h  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//���͈͋K���x��
	double *s1h2 = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//��̕⏕�z��
	//
	//printf("DELT(1.0)=  ");	scanf(" %lf",&delt);//���Ԃ����ݓ���	//delt=1.0;

	//temp=500.0;					//���x(K)
	//al=1000.0*1.0E-09;			//�v�Z�̈�[m]
	al=al*1.0E-09;				//�v�Z�̈�[m]
	b1=al/nd;					//�����u���b�N�̒���

	time1=0.0;					//�����v�Z�J�E���g���̐ݒ�
	//time1max=1.0+1.0e+07;		//�ő�v�Z�J�E���g���̐ݒ�

	//smob=1.0;					//���r���e�B�[�i�ϑԂ̊ɘa�W���j

	//AA0=200.0/rr/temp;			//�K��-�s�K���ϑԂ̉��w�I�쓮��
	AA0=AA0e/rr/temp;			//�K��-�s�K���ϑԂ̉��w�I�쓮��

	//kappa_s1=5.0e-15/rr/temp/b1/b1;//���z�G�l���M�|�W��
	kappa_s1=kappa_s1c/rr/temp/b1/b1;//���z�G�l���M�|�W��

	//a1_c=b1_c=c1_c=3.563E-10;	//�i�q�萔
	//atom_n=4.0;  vm0=6.02E23*a1_c*b1_c*c1_c/atom_n;	//�����̐ς̌v�Z�ifcc������j

//*** ������̐ݒ� ***************

	ini000(s1h, ND);

//**** �V�~�����[�V�����X�^�[�g ******************************
//Nstep = 10;
iout = -1;
start: ;

	//if(time1<=200.){Nstep=10;} else{Nstep=100;}		//�f�[�^�ۑ����鎞�ԊԊu�̕ύX
	if((((int)(time1) % Nstep)==0)) {datsave(s1h, ND);} 	//���J�Ԃ��J�E���g���ɑg�D�f�[�^��ۑ�
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(s1h, ND);} 	//���J�Ԃ��J�E���g���ɑg�D�f�[�^��ۑ�

//******  �|�e���V�����̌v�Z ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}

			//s1=s1h[i][j];
			//s1k_su=-kappa_s1*(s1h[ip][j]+s1h[im][j]+s1h[i][jp]+s1h[i][jm]-4.0*s1);
			s1=s1h[i*ND+j];
			s1k_su=-kappa_s1*(s1h[ip*ND+j]+s1h[im*ND+j]+s1h[i*ND+jp]+s1h[i*ND+jm]-4.0*s1);
															//���z�|�e���V����[��(4.16)]
			s1k_chem=-4.0*AA0*s1*(1.0-s1*s1);				//���w�|�e���V�����̌v�Z[��(4.14)]

			s1ddtt=-smob*(s1k_chem+s1k_su);					//��̎��Ԕ��W�̌v�Z[��(4.17)]
			//s1h2[i][j]=s1h[i][j]+s1ddtt*delt;				//�z��@
			s1h2[i*ND+j]=s1h[i*ND+j]+s1ddtt*delt;				//�z��@

			//if(s1h2[i][j]>=1.0){s1h2[i][j]=1.0;}  
			//if(s1h2[i][j]<=-1.0){s1h2[i][j]=-1.0;}			//s�̕ψ�(-1<=s<=1)�̕␳
			if(s1h2[i*ND+j]>= 1.0){s1h2[i*ND+j]= 1.0;}  
			if(s1h2[i*ND+j]<=-1.0){s1h2[i*ND+j]=-1.0;}			//s�̕ψ�(-1<=s<=1)�̕␳
		}
	}

	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){  s1h[i][j]=s1h2[i][j]; } }
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){  s1h[i*ND+j]=s1h2[i*ND+j]; } }
	//�⏕�z�񂩂��z��Ƀf�[�^�R�s�[

	//if(keypress()){return 0;}				//�L�[�҂����

	time1=time1+1.0;								//�v�Z�J�E���g���̉��Z
	if(time1<time1max){goto start;}	//�ő�J�E���g���ɓ��B�������ǂ����̔��f

end:;
  return 0;
}

//************ �����g�ݒ�T�u���|�`�� *******************************
void ini000(double *s1h, int ND)
{
	int i, j;
	srand(time(NULL)); // ����������
	int ndm=ND-1;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s1h[i][j]=0.01*(2.0*DRND(1.0)-1.0);//����ő�}1%�̗����ɂĐݒ�
			s1h[i*ND+j]=0.01*(2.0*DRND(1.0)-1.0);//����ő�}1%�̗����ɂĐݒ�
		}
	}
}

//************ �f�[�^�̕ۑ��T�u���[�`�� *******************************
void datsave(double *s1h, int ND)
{
	FILE		*stream;	//�X�g���[���̃|�C���^�ݒ�
	int 		i, j;			//����
	int ndm=ND-1;

	stream = fopen("test.dat", "a");	//�������ݐ�̃t�@�C����ǋL�����ŃI�[�v��

	fprintf(stream, "%e\n", time1);		//�J�Ԃ��J�E���g���̕ۑ�
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//fprintf(stream, "%e  ", s1h[i][j]);//��̃f�[�^�ۑ�
			fprintf(stream, "%e  ", s1h[i*ND+j]);//��̃f�[�^�ۑ�
		}
	}
	fprintf(stream, "\n");	//���s�̏�������
	fclose(stream);					//�t�@�C�����N���[�Y
}

void datsave_paraview(double *s1h, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j;
	int ndm=ND-1;
	
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
			//fprintf(stream, "%e  ", s1h[i][j]);//��̃f�[�^�ۑ�
			//fprintf(fp,"%10.6f\n", s1h[i][j]);
			fprintf(fp,"%10.6f\n", s1h[i*ND+j]);
		}
	}
	fclose(fp);
}
