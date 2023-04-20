#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//�����̊֐��ݒ�

//#define ND 128						//�����v�Z�ɂ�����v�Z�̈��ӂ̕�����(�����t�[���G�ϊ���p���邽�߂Q�ׂ̂���)
//#define IG 7						// 2^IG=ND

	//int nd=ND, ndm=ND-1; 			//�v�Z�̈�̈�ӂ̍���������(�����u���b�N��)�AND-1���`
	//int nd2=ND/2;				 	//ND/2���`�F�����t�|���G�ϊ��Ŏg�p
	//int ig=IG;						//2^ig=ND
	double PI=3.14159;				//�~����
	double RR=8.3145;				//�K�X�萔
	double time1;					//�v�Z�J�E���g��(���Ԃɔ��)
	double ep0=8.8541878e-12; 		//�^��̗U�d��(F/m)
	double Peq;						//���Ƀ��[�����g�̕��t�l

	//double s1h[ND][ND], s2h[ND][ND];//x�����̕��Ƀ��[�����g,y�����̕��Ƀ��[�����g

	double qs;						//�t�|���G�ϊ�(qs:-1)�ƃt�|���G�t�ϊ�(qs:1)�̋��
	//double xi[ND][ND], xr[ND][ND], xif[ND], xrf[ND];//�t�|���G�ϊ��̎����E�����z��
	//double s[ND],c[ND];				//sin��cos�̃e�[�u��
	//int ik[ND];						//�r�b�g���]�e�[�u��
	int iout;

	void ini000_s12(double *s1h, double *s2h, int ND);			//����0�ɂ����镪�Ƀ��[�����g�̏����v���t�@�C��
	void table(double *s, double *c, int *ik, int ND, int ig);	//sin��cos�̃e�[�u���ƃr�b�g���]�e�[�u���̍쐬�T�u���|�`��
	void fft(double *xrf, double *xif, double *s, double *c, int ND, int ig);	//�P���������t�[���G�ϊ�
	void rcfft(double *xrf, double *xif, double *xr, double *xi, int *ik, double *s, double *c, int ND, int ig);//�Q���������t�[���G�ϊ�
	void datsave(double *s1h, double *s2h, int ND);	//�f�|�^�ۑ��T�u���|�`��
	void datsave_paraview(double *s1h, double *s2h, int ND);//�f�|�^�ۑ��T�u���|�`��
	void datin(double *s1h, double *s2h, int ND);	//�����f�[�^�ǂݍ���

//******* ���C���v���O���� ******************************************
int main(void)
{
    int ND;
	int nd, ndm, nd2, ig;
	
	int   i, j, k, l, ii, jj, Nstep;		//����
	int   ip, im, jp, jm;					//����

	//double s1qrh[ND][ND], s1qih[ND][ND];	//s1�̃t�[���G�ϊ��i�����A�����j
	//double s2qrh[ND][ND], s2qih[ND][ND];	//s2�̃t�[���G�ϊ��i�����A�����j
	//double s1h2[ND][ND], s2h2[ND][ND];		//s1��s2�̕⏕�z��

	//double ss1qrh[ND][ND], ss1qih[ND][ND];	//s1*s1�̃t�[���G�ϊ��i�����A�����j
	//double ss2qrh[ND][ND], ss2qih[ND][ND];	//s2*s2�̃t�[���G�ϊ��i�����A�����j
	//double s1s2qrh[ND][ND], s1s2qih[ND][ND];//s1*s2�̃t�[���G�ϊ��i�����A�����j
	double ss1ss2;							//�ψ�Ɋւ��鐔�l�v�Z�덷�̕␳�p�̍�ƕϐ�

	double a0_aa, a0_a, a0_c;				//BaTiO3�̊i�q�萔(������)
	double al, temp, delt, vm0;				//�v�Z�̈�P�҂̒����A���x�A���Ԃ����݁A�����̐�
	double time1max;						//�v�Z�J�E���g���̍ő�l�i�v�Z�I���J�E���g�j
	double b1;								//�����u���b�N�P�ӂ̒���

	double s1, s2;							//���Ƀ��[�����g��x�����y����
	double s1ip, s1im, s1jp, s1jm;			//s1�̍��E�㉺�̒l
	double s2ip, s2im, s2jp, s2jm;			//s2�̍��E�㉺�̒l

	double s1k, s1k_chem, s1k_surf, s1k_str, s1k_ddi;	//s1�Ɋւ���|�e���V����
	//double s1k_dd[ND][ND];					//�o�Ɏq-�o�Ɏq���ݍ�p�|�e���V����
	double s2k, s2k_chem, s2k_surf, s2k_str, s2k_ddi;	//s2�Ɋւ���|�e���V����
	//double s2k_dd[ND][ND];					//�o�Ɏq-�o�Ɏq���ݍ�p�|�e���V����
	double smob1, smob2;					//�h���C���E�ʂ̈Փ��x
	double s1ddtt, s2ddtt;					//���W�������̍���

	double A1, A11, A12;					//���w�I���R�G�l���M�[���̃p�����[�^
	double A1e, A11e, A12e;
	double A111, A112, A123;
	double A111e, A112e, A123e;
	double A1111, A1112, A1122, A1123;
	double A1111e, A1112e, A1122e, A1123e;
	double Tc0;

	double nx, ny, alnn;					//�t��Ԃ̒P�ʃx�N�g�������Ƌt�i�q�x�N�g����
	double kapaP, kapaPc;					//���z�G�l���M�[�W��

	double E1_ex, E1_ex_0; 					//�O���d��
	double ep00;							//��U�d��
	double Add0, Add0c;						//�o�Ɏq-�o�Ɏq���ݍ�p�v�Z�ɂ�����W��

	double t1, t2, t3, tQ, tR, tS, tT;		//���t���Ƀ��[�����g�Z�o���̍�ƕϐ�

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
	delt    = data[1];
	time1max= int(data[2]);
	Nstep   = int(data[3]);
	temp    = data[4];
	al      = data[5];
	vm0     = data[6];
	smob1   = data[7];
	smob2   = data[8];
	A1e     = data[9];
	A11e    = data[10];
	A12e    = data[11];
	A111e   = data[12];
	A112e   = data[13];
	A123e   = data[14];
	A1111e  = data[15];
	A1112e  = data[16];
	A1122e  = data[17];
	A1123e  = data[18];
	Tc0     = data[19];
	kapaPc  = data[20];
	E1_ex_0 = data[21];
	ep00    = data[22];
	Add0c   = data[23];
	printf("---------------------------------\n");
	//
	ig = int(log2(ND));
	nd=ND;
	ndm=ND-1;
	nd2=ND/2;
	//
	double *s1h     = (double *)malloc(sizeof(double)*( ND*ND ));	//x�����̕��Ƀ��[�����g
	double *s2h     = (double *)malloc(sizeof(double)*( ND*ND ));	//y�����̕��Ƀ��[�����g
	//
	double *xi      = (double *)malloc(sizeof(double)*( ND*ND ));//�t�|���G�ϊ��̎����E�����z��
	double *xr      = (double *)malloc(sizeof(double)*( ND*ND ));//�t�|���G�ϊ��̎����E�����z��
	//
	double *xif     = (double *)malloc(sizeof(double)*( ND ));//�t�|���G�ϊ��̎����E�����z��
	double *xrf     = (double *)malloc(sizeof(double)*( ND ));//�t�|���G�ϊ��̎����E�����z��
	//
	double *s       = (double *)malloc(sizeof(double)*( ND ));//sin�̃e�[�u��
	double *c       = (double *)malloc(sizeof(double)*( ND ));//cos�̃e�[�u��
	//
	int *ik         = (int *)malloc(sizeof(int)*( ND ));//�r�b�g���]�e�[�u��
	//
	double *s1qrh   = (double *)malloc(sizeof(double)*( ND*ND ));	//s1�̃t�[���G�ϊ��i�����j
	double *s1qih   = (double *)malloc(sizeof(double)*( ND*ND ));	//s1�̃t�[���G�ϊ��i�����j
	double *s2qrh   = (double *)malloc(sizeof(double)*( ND*ND ));	//s2�̃t�[���G�ϊ��i�����j
	double *s2qih   = (double *)malloc(sizeof(double)*( ND*ND ));	//s2�̃t�[���G�ϊ��i�����j
	double *s1h2    = (double *)malloc(sizeof(double)*( ND*ND ));	//s1�̕⏕�z��
	double *s2h2    = (double *)malloc(sizeof(double)*( ND*ND ));	//s2�̕⏕�z��
	double *ss1qrh  = (double *)malloc(sizeof(double)*( ND*ND ));	//s1*s1�̃t�[���G�ϊ��i�����j
	double *ss1qih  = (double *)malloc(sizeof(double)*( ND*ND ));	//s1*s1�̃t�[���G�ϊ��i�����j
	double *ss2qrh  = (double *)malloc(sizeof(double)*( ND*ND ));	//s2*s2�̃t�[���G�ϊ��i�����j
	double *ss2qih  = (double *)malloc(sizeof(double)*( ND*ND ));	//s2*s2�̃t�[���G�ϊ��i�����j
	double *s1s2qrh = (double *)malloc(sizeof(double)*( ND*ND ));	//s1*s2�̃t�[���G�ϊ��i�����j
	double *s1s2qih = (double *)malloc(sizeof(double)*( ND*ND ));	//s1*s2�̃t�[���G�ϊ��i�����j
	//
	double *s1k_dd  = (double *)malloc(sizeof(double)*( ND*ND ));	//�o�Ɏq-�o�Ɏq���ݍ�p�|�e���V����
	double *s2k_dd  = (double *)malloc(sizeof(double)*( ND*ND ));	//�o�Ɏq-�o�Ɏq���ݍ�p�|�e���V����
	//
	//printf("DELT(0.1)=  "); scanf(" %lf",&delt);//���ԍ��ݐݒ�	//delt=0.1;

	time1=0.0;								//�����v�Z�J�E���g���̐ݒ�
	//time1max=1.0+1.0e+07;					//�ő�v�Z�J�E���g���̐ݒ�
	//���Ԃ͑S�Ė�����������Ă���̂Œ��ӁB

	//temp=298.0;								//���x(K)

	//al=250.;								//�v�Z�̈�̂P�ӂ̒����i��m�j
	b1=al*1.0E-06/nd;						//�����u���b�N�P�ӂ̒����im�j

	//a0_aa=0.4;  a0_a=0.3992;  a0_c=0.4036;	//�i�q�萔(nm)
	//vm0=a0_a*a0_a*a0_c*1.0e-27*6.02e+23;	//�����̐ρi���q�P�����j

	//smob1=1.; smob2=1.;						//�\�����]�ڂɂ�����Փ��x�i�K�i�����ĂP�ɐݒ�j

//--- ���w�I���R�G�l���M�[���̃p�����[�^�l [�\4.7�Q��]----------------------
	//A1=4.124e+05*vm0/RR/temp;
	A1=A1e*vm0/RR/temp;
	//A11=-2.097e+08*vm0/RR/temp;
	A11=A11e*vm0/RR/temp;
	//A12=7.974e+08*vm0/RR/temp;
	A12=A12e*vm0/RR/temp;
	//A111=1.294e+09*vm0/RR/temp;
	A111=A111e*vm0/RR/temp;
	//A112=-1.950e+09*vm0/RR/temp;
	A112=A112e*vm0/RR/temp;
	//A123=-2.500e+09*vm0/RR/temp;
	A123=A123e*vm0/RR/temp;
	//A1111=3.863e+10*vm0/RR/temp;
	A1111=A1111e*vm0/RR/temp;
	//A1112=2.529e+10*vm0/RR/temp;
	A1112=A1112e*vm0/RR/temp;
	//A1122=1.637e+10*vm0/RR/temp;
	A1122=A1122e*vm0/RR/temp;
	//A1123=1.367e+10*vm0/RR/temp;
	A1123=A1123e*vm0/RR/temp;

	//Tc0=115.0+273.0;  //K

//--- ���t���Ƀ��[�����g�l�̌v�Z ----------------------
	t1=3.0*A111/4.0/A1111;
	t2=A11/2.0/A1111;
	t3=A1*(temp-Tc0)/4.0/A1111;
	tQ=(3.0*t2-t1*t1)/9.0;
	tR=(9.0*t1*t2-27.0*t3-2.0*t1*t1*t1)/54.0;
	tS=pow((tR+sqrt(tQ*tQ*tQ+tR*tR)),(1.0/3.0));
	tT=pow((tR-sqrt(tQ*tQ*tQ+tR*tR)),(1.0/3.0));
	Peq=sqrt(fabs(tS+tR-t1/3.));//���t���Ƀ��[�����g�l
	printf("Peq=  %f  \n", Peq);//���t���Ƀ��[�����g�l�̕\��

	//kapaP=1.0e-15/ep0*vm0/RR/temp/b1/b1;//���z�G�l���M�[�W��
	kapaP=kapaPc/ep0*vm0/RR/temp/b1/b1;//���z�G�l���M�[�W��

	//E1_ex_0=0.0;//�O���d��(x����)
	//E1_ex_0=2.0e+6*vm0/RR/temp;//�O���d��(x����)

	//ep00=200.0;//��U�d��

	//Add0=1.0/ep0/ep00*vm0/RR/temp;//�o�Ɏq-�o�Ɏq���ݍ�p�v�Z�ɂ�����W��[��(4.52)�Q��]
	Add0=Add0c/ep0/ep00*vm0/RR/temp;//�o�Ɏq-�o�Ɏq���ݍ�p�v�Z�ɂ�����W��

//*** sin�����cos�e�|�u���A�r�b�g���]�e�[�u���A����я�����̐ݒ� ***************

	table(s, c, ik, ND, ig);		//sin�����cos�e�|�u���ƃr�b�g���]�e�[�u���̐ݒ�
	ini000_s12(s1h, s2h, ND);		//����0s�ɂ����镪�Ƀ��[�����g�̏����v���t�@�C��
	//datin(s1h, s2h, ND);			//�����g�D��̓���

//**** �V�~�����[�V�����X�^�[�g ******************************
//Nstep = 10;
start: ;

	//if(time1<=100.){Nstep=10;} else{Nstep=100;}//�f�[�^�ۑ����鎞�ԊԊu�̕ύX
	if((((int)(time1) % Nstep)==0)){datsave(s1h, s2h, ND);} //���J�Ԃ��J�E���g���ɑg�D�f�[�^��ۑ�
	if((((int)(time1) % Nstep)==0)){datsave_paraview(s1h, s2h, ND);} //���J�Ԃ��J�E���g���ɑg�D�f�[�^��ۑ�
	//if((int)(time1)==2000){datsave(s1h, s2h, ND);} //����v�Z�J�E���g�ɂ�����f�|�^�̕ۑ�

//**** s1�̃t�[���G�ϊ�[��(3.37)] ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=s1h[i][j]; xi[i][j]=0.;
			xr[i*ND+j]=s1h[i*ND+j];
			xi[i*ND+j]=0.;
		}
	}
	qs=-1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s1qrh[i][j]=xr[i][j];  s1qih[i][j]=xi[i][j];
			s1qrh[i*ND+j]=xr[i*ND+j];
			s1qih[i*ND+j]=xi[i*ND+j];
		}
	}
	//s1qrh[0][0]=s1qih[0][0]=0.;
	s1qrh[0]=s1qih[0]=0.0;

//**** s2�̃t�[���G�ϊ�[��(3.37)] ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=s2h[i][j]; xi[i][j]=0.;
			xr[i*ND+j]=s2h[i*ND+j];
			xi[i*ND+j]=0.0;
		}
	}
	qs=-1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s2qrh[i][j]=xr[i][j];  s2qih[i][j]=xi[i][j];
			s2qrh[i*ND+j]=xr[i*ND+j];
			s2qih[i*ND+j]=xi[i*ND+j];
		}
	}
	//s2qrh[0][0]=s2qih[0][0]=0.;
	s2qrh[0]=s2qih[0]=0.0;

//***** s1�̑o�Ɏq-�o�Ɏq���ݍ�p�̌v�Z[��(3.33)���̏􍞂ݐϕ�] ***********************
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj);
			if(alnn==0.){alnn=1.;}
			nx=(double)ii/alnn;  ny=(double)jj/alnn;
			//xr[i][j]=s1qrh[i][j]*nx*nx+s2qrh[i][j]*nx*ny;
			//xi[i][j]=s1qih[i][j]*nx*nx+s2qih[i][j]*nx*ny;
			xr[i*ND+j]=s1qrh[i*ND+j]*nx*nx+s2qrh[i*ND+j]*nx*ny;
			xi[i*ND+j]=s1qih[i*ND+j]*nx*nx+s2qih[i*ND+j]*nx*ny;
		}
	}
	qs=1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ s1k_dd[i][j]=xr[i][j]; } }
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ s1k_dd[i*ND+j]=xr[i*ND+j]; } }

//***** s2�̑o�Ɏq-�o�Ɏq���ݍ�p�̌v�Z[��(3.33)���̏􍞂ݐϕ�] ************************
	for(i=0;i<=ndm;i++){
		if(i<=nd2-1){ii=i;}  if(i>=nd2){ii=i-nd;}
		for(j=0;j<=ndm;j++){
			if(j<=nd2-1){jj=j;}  if(j>=nd2){jj=j-nd;}
			alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj);
			if(alnn==0.){alnn=1.;}
			nx=(double)ii/alnn;  ny=(double)jj/alnn;
			//xr[i][j]=s1qrh[i][j]*ny*nx+s2qrh[i][j]*ny*ny;
			//xi[i][j]=s1qih[i][j]*ny*nx+s2qih[i][j]*ny*ny;
			xr[i*ND+j]=s1qrh[i*ND+j]*ny*nx+s2qrh[i*ND+j]*ny*ny;
			xi[i*ND+j]=s1qih[i*ND+j]*ny*nx+s2qih[i*ND+j]*ny*ny;
		}
	}
	qs=1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ s2k_dd[i][j]=xr[i][j]; } }
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ s2k_dd[i*ND+j]=xr[i*ND+j]; } }

//******  ����`���W�������̐��l�v�Z  ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1;  im=i-1;  jp=j+1;  jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//�����I���E����
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
			//s1=s1h[i][j]; s1ip=s1h[ip][j]; s1im=s1h[im][j]; s1jp=s1h[i][jp]; s1jm=s1h[i][jm];
			//s2=s2h[i][j]; s2ip=s2h[ip][j]; s2im=s2h[im][j]; s2jp=s2h[i][jp]; s2jm=s2h[i][jm];
			s1=s1h[i*ND+j]; s1ip=s1h[ip*ND+j]; s1im=s1h[im*ND+j]; s1jp=s1h[i*ND+jp]; s1jm=s1h[i*ND+jm];
			s2=s2h[i*ND+j]; s2ip=s2h[ip*ND+j]; s2im=s2h[im*ND+j]; s2jp=s2h[i*ND+jp]; s2jm=s2h[i*ND+jm];

			//���z�|�e���V�����̌v�Z
			s1k_surf=-kapaP*(s1ip+s1im+s1jp+s1jm-4.*s1);
			s2k_surf=-kapaP*(s2ip+s2im+s2jp+s2jm-4.*s2);

			//���w�|�e���V�����̌v�Z[��(4.55)]
			s1k_chem=2.0*A1*(temp-Tc0)*s1
							+4.0*A11*s1*s1*s1+2.0*A12*s1*s2*s2
							+6.0*A111*pow(s1,5.0)+A112*(2.0*s1*pow(s2,4.0)+4.0*s2*s2*pow(s1,3.0))
							+8.0*A1111*pow(s1,7.0)
							+A1112*(6.0*pow(s1,5.0)*s2*s2+2.0*pow(s2,6.0)*s1)
							+4.0*A1122*pow(s1,3.0)*pow(s2,4.0);

			s2k_chem=2.0*A1*(temp-Tc0)*s2
							+4.0*A11*s2*s2*s2+2.0*A12*s2*s1*s1
							+6.0*A111*pow(s2,5.0)+A112*(2.0*s2*pow(s1,4.0)+4.0*s1*s1*pow(s2,3.0))
							+8.0*A1111*pow(s2,7.0)
							+A1112*(6.0*pow(s2,5.0)*s1*s1+2.0*pow(s1,6.0)*s2)
							+4.0*A1122*pow(s2,3.0)*pow(s1,4.0);

			//�o�Ɏq-�o�Ɏq�|�e���V�����̌v�Z[��(4.56)�E�ӑ�P��] 
			//s1k_ddi=Add0*s1k_dd[i][j];
			//s2k_ddi=Add0*s2k_dd[i][j];
			s1k_ddi=Add0*s1k_dd[i*ND+j];
			s2k_ddi=Add0*s2k_dd[i*ND+j];

			//�S�|�e���V�����̌v�Z
			s1k=s1k_chem+s1k_surf+s1k_ddi-E1_ex_0;
			s2k=s2k_chem+s2k_surf+s2k_ddi;

			//���W������[��(4.57)]��s��̎��Ԕ��W�i�z��@�j�̌v�Z
			//s1ddtt=-smob1*s1k;  s1h2[i][j]=s1h[i][j]+s1ddtt*delt;
			//s2ddtt=-smob2*s2k;  s2h2[i][j]=s2h[i][j]+s2ddtt*delt;
			s1ddtt=-smob1*s1k;  s1h2[i*ND+j]=s1h[i*ND+j]+s1ddtt*delt;
			s2ddtt=-smob2*s2k;  s2h2[i*ND+j]=s2h[i*ND+j]+s2ddtt*delt;
		}
	}

	//s�����z��Ɉړ�����ѕψ�Ɋւ��鐔�l�v�Z�덷�̕␳
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s1h[i][j]=s1h2[i][j]; s2h[i][j]=s2h2[i][j];
			s1h[i*ND+j]=s1h2[i*ND+j]; s2h[i*ND+j]=s2h2[i*ND+j];
			//if(s1h[i][j]>=Peq){s1h[i][j]=Peq;}  if(s1h[i][j]<=-Peq){s1h[i][j]=-Peq;}
			if(s1h[i*ND+j]>=Peq){s1h[i*ND+j]=Peq;}  if(s1h[i*ND+j]<=-Peq){s1h[i*ND+j]=-Peq;}
			//if(s2h[i][j]>=Peq){s2h[i][j]=Peq;}  if(s2h[i][j]<=-Peq){s2h[i][j]=-Peq;}
			if(s2h[i*ND+j]>=Peq){s2h[i*ND+j]=Peq;}  if(s2h[i*ND+j]<=-Peq){s2h[i*ND+j]=-Peq;}
			//ss1ss2=sqrt(s1h[i][j]*s1h[i][j]+s2h[i][j]*s2h[i][j]);
			ss1ss2=sqrt(s1h[i*ND+j]*s1h[i*ND+j]+s2h[i*ND+j]*s2h[i*ND+j]);
			//if(ss1ss2>=Peq){s1h[i][j]=s1h[i][j]/ss1ss2*Peq; s2h[i][j]=s2h[i][j]/ss1ss2*Peq;}
			if(ss1ss2>=Peq){s1h[i*ND+j]=s1h[i*ND+j]/ss1ss2*Peq; s2h[i*ND+j]=s2h[i*ND+j]/ss1ss2*Peq;}
		}
	}

	//���Ԃ�i�߁A�I�����ԂłȂ����ɂ̓X�^�[�g��
	//if(keypress()){return 0;}//�L�[�҂����
	time1=time1+1.0;								//�v�Z�J�E���g���̉��Z
	if(time1<time1max){goto start;}	//�ő�J�E���g���ɓ��B�������ǂ����̔��f

	end:;
  return 0;
}

//************ ������̐ݒ�T�u���|�`�� *************
void ini000_s12(double *s1h, double *s2h, int ND)
{
	int i, j;	//����
 	//srand(time(NULL)); // ����������
	int nd=ND, ndm=ND-1, nd2=ND/2;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s1h[i][j]=0.01*(2.0*DRND(1)-1.0);//����ő�}1%�̗����ɂĐݒ�
			//s2h[i][j]=0.01*(2.0*DRND(1)-1.0);
			s1h[i*ND+j]=0.01*(2.0*DRND(1)-1.0);//����ő�}1%�̗����ɂĐݒ�
			s2h[i*ND+j]=0.01*(2.0*DRND(1)-1.0);
		}
	}
}

//******* Sin, Cos �̃e�[�u������уr�b�g���]�e�[�u���̐ݒ� ***************
void table(double *s, double *c, int *ik, int ND, int ig)
{
	int it, it1, it2, mc, mn;
	double q;
	int nd=ND, ndm=ND-1, nd2=ND/2;

	q=2.*PI/nd;
	for(it=0;it<=nd2-1;it++){
		c[it]=cos(q*it);  s[it]=sin(q*it);//Sin, Cos �̃e�[�u��
	}
	ik[0]=0;
	mn=nd2;
	mc=1;
	for(it1=1;it1<=ig;it1++){
		for(it2=0;it2<=mc-1;it2++){
			ik[it2+mc]=ik[it2]+mn;//�r�b�g���]�e�[�u��
		}
		mn=mn/2;
		mc=2*mc;
	}
}

//********** �P���������t�[���G�ϊ� **************************************
void fft(double *xrf, double *xif, double *s, double *c, int ND, int ig)
{
	int ix, ka, kb, l2, lf, mf, n2, nf;
	double tj, tr;
	int nd=ND, ndm=ND-1, nd2=ND/2;

	l2=1;
	for(lf=1;lf<=ig;lf++){
		n2=nd2/l2;
		for(mf=1;mf<=l2;mf++){
			for(nf=0;nf<=n2-1;nf++){
				ix=nf*l2;
				ka=nf+2*n2*(mf-1);
				kb=ka+n2;
				tr=xrf[ka]-xrf[kb];
				tj=xif[ka]-xif[kb];
				xrf[ka]=xrf[ka]+xrf[kb];
				xif[ka]=xif[ka]+xif[kb];
				xrf[kb]=tr*c[ix]-tj*qs*s[ix];
				xif[kb]=tj*c[ix]+tr*qs*s[ix];
			}
		}
		l2=l2*2;
	}

}

//************ �Q���������t�[���G�ϊ� ***********************************
void rcfft(double *xrf, double *xif, double *xr, double *xi, int *ik, double *s, double *c, int ND, int ig)
{
	int i, ic, ir, j;
	int nd=ND, ndm=ND-1, nd2=ND/2;

	for(ir=0;ir<=ndm;ir++){
		for(ic=0;ic<=ndm;ic++){
			//xrf[ic]=xr[ir][ic];  xif[ic]=xi[ir][ic];
			xrf[ic]=xr[ir*ND+ic];
			xif[ic]=xi[ir*ND+ic];
		}
		fft(xrf, xif, s, c, ND, ig);
		for(ic=0;ic<=ndm;ic++){
			//xr[ir][ic]=xrf[ik[ic]];  xi[ir][ic]=xif[ik[ic]];
			xr[ir*ND+ic]=xrf[ik[ic]];
			xi[ir*ND+ic]=xif[ik[ic]];
		}
	}
	for(ic=0;ic<=ndm;ic++){
		for(ir=0;ir<=ndm;ir++){
			//xrf[ir]=xr[ir][ic];  xif[ir]=xi[ir][ic];
			xrf[ir]=xr[ir*ND+ic];
			xif[ir]=xi[ir*ND+ic];
		}
		fft(xrf, xif, s, c, ND, ig);
		for(ir=0;ir<=ndm;ir++){
			//xr[ir][ic]=xrf[ik[ir]];  xi[ir][ic]=xif[ik[ir]];
			xr[ir*ND+ic]=xrf[ik[ir]];
			xi[ir*ND+ic]=xif[ik[ir]];
		}
	}
	if(qs>0.){return;}
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=xr[i][j]/nd/nd;  xi[i][j]=xi[i][j]/nd/nd;
			xr[i*ND+j]=xr[i*ND+j]/nd/nd;
			xi[i*ND+j]=xi[i*ND+j]/nd/nd;
		}
	}

}

//************ �f�[�^�ۑ��T�u���[�`�� *******************************
void datsave(double *s1h, double *s2h, int ND)
{
	FILE		*stream;	//�X�g���[���̃|�C���^�ݒ�
	int 		i, j;			//����
	int nd=ND, ndm=ND-1, nd2=ND/2;

	stream = fopen("test.dat", "a");	//�������ݐ�̃t�@�C����ǋL�����ŃI�[�v��
	fprintf(stream, "%e \n", time1);	//�v�Z�J�E���g���̕ۑ�
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//fprintf(stream, "%e  %e  ", s1h[i][j], s2h[i][j]);	//��̃f�[�^�ۑ�
			fprintf(stream, "%e  %e  ", s1h[i*ND+j], s2h[i*ND+j]);	//��̃f�[�^�ۑ�
		}
	}
	fprintf(stream, "\n");	//���s�̏�������
	fclose(stream);					//�t�@�C�����N���[�Y
}

void datsave_paraview(double *s1h, double *s2h, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j;
	int nd=ND, ndm=ND-1, nd2=ND/2;
	
	iout = iout + 1;
	printf("dip_result%06d.vtk \n",iout);
	sprintf(fName,"dip_result%06d.vtk",iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),(1));
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*1));
	fprintf(fp,"SCALARS Phase_field_1 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e  ", s1h[i][j], s2h[i][j]);	//��̃f�[�^�ۑ�
			//fprintf(fp,"%10.6f\n", s1h[i][j]);
			fprintf(fp,"%10.6f\n", s1h[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Phase_field_2 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e  ", s1h[i][j], s2h[i][j]);	//��̃f�[�^�ۑ�
			//fprintf(fp,"%10.6f\n", s2h[i][j]);
			fprintf(fp,"%10.6f\n", s2h[i*ND+j]);
		}
	}
	fclose(fp);
}
//************ �f�[�^�ǂݍ��݃T�u���[�`�� *******************************
void datin(double *s1h, double *s2h, int ND)
{
	FILE		*datin0;	//�X�g���[���̃|�C���^�ݒ�
	int 		i, j;			//����
	double time_ext;	//�v�Z�J�E���g��
	int nd=ND, ndm=ND-1, nd2=ND/2;

	datin0 = fopen("test.dat", "r");	//�ǂݍ��݌��̃t�@�C�����I�[�v��
	fscanf(datin0, "%lf", &time_ext);		//�v�Z�J�E���g���̓Ǎ�
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//fscanf(datin0, "%lf  %lf  ", &s1h[i][j], &s2h[i][j]);	//��̃f�[�^�Ǎ�
			fscanf(datin0, "%lf  %lf  ", &s1h[i*ND+j], &s2h[i*ND+j]);	//��̃f�[�^�Ǎ�
		}
	}
	fclose(datin0);//�t�@�C�����N���[�Y
}