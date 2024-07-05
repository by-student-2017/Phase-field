//�t�F�[�Y�t�B�[���h�̃f�[�^(ph.dat)��ǂ�ŁA
//�e������v�Z����v���O�����B
//�e�e����́A�ꊇ����el_field.dat�ɕۑ������B
//�����ɁA�e�e����̂Q�����摜���ۑ������B

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//�����̊֐��ݒ�
//#define ND 128					//�����v�Z�ɂ�����v�Z�̈��ӂ̕�����(�����t�[���G�ϊ���p���邽�߂Q�ׂ̂���)
//#define IG 7					//2^IG=ND

	//int nd=ND, ndm=ND-1; 		//�v�Z�̈�̈�ӂ̍���������(�����u���b�N��)�AND-1���`
	//int nd2=ND/2;				//ND/2���`�F�����t�|���G�ϊ��Ŏg�p
	//int ig=IG;					//2^ig=ND
	double PI=3.14159;			//�~����
	double rr=8.3145;			//�K�X�萔
	//double time1;				//�v�Z�J�E���g��(���Ԃɔ��)
	int iout=-1;

	double c0;					//�t�F�[�Y�t�B�[���h�̕��ϒl
	//double ch[ND][ND];			//��i�t�F�[�Y�t�B�[���h�j
	double eta_c[4][4];			//�i�q�~�X�}�b�`
	double cec[4][4][4][4];		//�e���萔
	double sigma[4][4];			//�A�C�Q�����́i�A�C�Q���Ђ��ݕ������e���ό`�������̉��́j
	//double ep11c[ND][ND], ep22c[ND][ND], ep12c[ND][ND];									//�S�Ђ���
	//double sig11[ND][ND], sig22[ND][ND], sig33[ND][ND], sig12[ND][ND];	//�e������
	//double Estr[ND][ND];		//�e���Ђ��݃G�l���M�[
	//double u1[ND][ND], u2[ND][ND], u3[ND][ND];	//�ψ�
	//double fld1[ND][ND];		//��̎󂯓n���p�̍�ƍs��

	int qs;					//�t�|���G�ϊ�(qs:-1)�ƃt�|���G�t�ϊ�(qs:1)�̋��
	//double xi[ND][ND], xr[ND][ND], xif[ND], xrf[ND];//�t�|���G�ϊ��̎����E�����z��
	//double s[ND],c[ND];	//sin��cos�̃e�[�u��
	//int ik[ND];					//�r�b�g���]�e�[�u��

	void table(double *s, double *c, int *ik, int ND, int ig);	//sin��cos�̃e�[�u���ƃr�b�g���]�e�[�u���̍쐬�T�u���|�`��
	void fft(double *xrf, double *xif, double *s, double *c, int ND, int ig);	//�P���������t�[���G�ϊ�
	void rcfft(double *xrf, double *xif, double *xr, double *xi, int *ik, double *s, double *c, int ND, int ig);	//�Q���������t�[���G�ϊ�
	double zcij(int i0, int j0, int iii, int jjj, int ND);//�e���G�l���M�[�֐��̌W���v�Z�i�t�[���G��ԁj
	double zuij(int i0, int j0, int iii, int ND);//�ψʏ�̌W���v�Z�i�t�[���G��ԁj

	void datin(double *ch, int ND);	//������ǂݍ��ݗp�T�u���|�`��
	void datsave(double *ch, double *Estr, 
			 double *ep11c, double *ep22c, double *ep12c,
			 double *sig11, double *sig22, double *sig33, double *sig12, 
			 double *u1, double *u2, int ND);	//�f�|�^�ۑ��T�u���|�`��
	void datsave_paraview(double *ch, double *Estr, 
			 double *ep11c, double *ep22c, double *ep12c,
			 double *sig11, double *sig22, double *sig33, double *sig12, 
			 double *u1, double *u2, int ND);	//�f�|�^�ۑ��T�u���|�`��

//******* ���C���v���O���� ****************************************************
int main(void)
{
    int ND, IG;
	int ndm, ig;
	
	//int i, j, k, l, ii, jj, kk, iii, jjj, ief;			//����
	int i, j, k, l, iii, jjj;	//integer
	//int ip, im, jp, jm;									//����

	//double c;	//��i�t�F�[�Y�t�B�[���h�j
	//double qrh1[ND][ND], qih1[ND][ND];					//��̃t�[���G�ϊ�
	//double ec11[ND][ND], ec22[ND][ND], ec12[ND][ND];	//�S���c�z��

	//double al, temp, delt;								//�v�Z�̈�̂P�ӂ̒����A���x�A���Ԃ�����
	//double time1max;									//�v�Z�J�E���g���̍ő�l�i�v�Z���~�߂�ۂɎg�p�j
	//double b1, vm0, atom_n;								//�����u���b�N�̂P�ӂ̒����A�����̐ρA�P�ʖE���̌��q��
	//double a1;											//Fe�̊i�q�萔

	double el_mag;										//�e�����Ɋւ����ƕϐ�
	double c11, c44, c12; 								//�e����
	double ep000;										//�i�q�~�X�}�b�`
	//double epc11, epc22, epc33, epc12, epc13, epc23;	//�S���Ђ���
	//double ep011, ep022, ep033, ep012, ep013, ep023;	//�A�C�Q���Ђ���
	//double dep011, dep022, dep033, dep012, dep023, dep013;
	//double epc011, epc022, epc033, epc012, epc023, epc013;
	//double ef1, ef2, ef11, ef22, ef12;					//�e���֐�

	//double Estr1, Estr2, Estr3, Estr4, Estr5;			//�e���Ђ��݃G�l���M�[�v�Z�̍ۂ̍�ƕϐ�
	double Estr1, Estr2, Estr3;							//�e���Ђ��݃G�l���M�[�v�Z�̍ۂ̍�ƕϐ�

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
	// pfm1 series
	ND      = int(data[0]);
	//c2a     = data[1];
	//delt    = data[2];	// timestep
	//temp    = data[3];	// [K]
	//al      = data[4];	// [nm]
	//cmob22  = data[5];
	//om_12e  = data[6];
	//kapa_c2c= data[7];
	//time1max= data[8];
	//Nstep   = int(data[9]);
	ep000   = data[10];
	el_mag  = data[11];
	c11     = data[12];
	c12     = data[13];
	c44     = data[14];
	printf("---------------------------------\n");
	//
	IG=int(log2(ND));
	ndm=ND-1;				//define ND-1
	ig=IG;					//2^ig=ND
	//
	//double ch[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	//a[z * ySize * xSize + y * xSize + x], a[i][j][k]=a[i*n*m + j*n + k]
	double *ch    = (double *)malloc(sizeof(double)*( ND*ND ));	//phase field
	//
	double *ep11c = (double *)malloc(sizeof(double)*( ND*ND ));	//Strain
	double *ep22c = (double *)malloc(sizeof(double)*( ND*ND ));
	//double *ep33c = (double *)malloc(sizeof(double)*( ND*ND ));
	double *ep12c = (double *)malloc(sizeof(double)*( ND*ND ));
	//double *ep13c = (double *)malloc(sizeof(double)*( ND*ND ));
	//double *ep23c = (double *)malloc(sizeof(double)*( ND*ND ));
	//
	double *sig11 = (double *)malloc(sizeof(double)*( ND*ND ));	//Elastic stress
	double *sig22 = (double *)malloc(sizeof(double)*( ND*ND ));
	double *sig33 = (double *)malloc(sizeof(double)*( ND*ND ));
	double *sig12 = (double *)malloc(sizeof(double)*( ND*ND ));
	//double *sig13 = (double *)malloc(sizeof(double)*( ND*ND ));
	//double *sig23 = (double *)malloc(sizeof(double)*( ND*ND ));
	//
	double *Estr  = (double *)malloc(sizeof(double)*( ND*ND ));	//Elastic strain energy density
	//
	double *u1    = (double *)malloc(sizeof(double)*( ND*ND ));	//Displacement
	double *u2    = (double *)malloc(sizeof(double)*( ND*ND ));
	//double *u3    = (double *)malloc(sizeof(double)*( ND*ND ));
	//
	double *qrh1  = (double *)malloc(sizeof(double)*( ND*ND ));	//Fourier transform of the field
	double *qih1  = (double *)malloc(sizeof(double)*( ND*ND ));
	//
	double *ec11  = (double *)malloc(sizeof(double)*( ND*ND ));	//Constrained strain array
	double *ec22  = (double *)malloc(sizeof(double)*( ND*ND ));
	//double *ec33  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	double *ec12  = (double *)malloc(sizeof(double)*( ND*ND ));
	//double *ec13  = (double *)malloc(sizeof(double)*( ND*ND ));
	//double *ec23  = (double *)malloc(sizeof(double)*( ND*ND ));
	//
	double *xr    = (double *)malloc(sizeof(double)*( ND*ND ));	//Real or imaginary array of Fourier transform
	double *xi    = (double *)malloc(sizeof(double)*( ND*ND ));
	double *xrf   = (double *)malloc(sizeof(double)*( ND ));
	double *xif   = (double *)malloc(sizeof(double)*( ND ));
	//
	double *s     = (double *)malloc(sizeof(double)*( ND ));	//sin and cos table
	double *c     = (double *)malloc(sizeof(double)*( ND ));
	//
	int *ik       = (int *)malloc(sizeof(int)*( ND ));	//Bit inversion table for FFT
	
	//
	//temp=900.;				//���x(K)
	//al=100.*1.0E-09;	//�v�Z�̈�(m)
	//b1=al/nd;					//�����u���b�N�̒���

	//time1=0.0;						//�����v�Z�J�E���g���̐ݒ�
	//time1max=1.0+1.0e+07;	//�ő�v�Z�J�E���g���̐ݒ�

	//a1=2.8664E-10; 	//Fe(bcc)�̊i�q�萔�i�����f�|�^�j
	//atom_n=2.;	vm0=6.02E23*a1*a1*a1/atom_n;	//�����̐ς̌v�Z�ibcc���̏ꍇ�j

//*** �i�q�~�X�}�b�`�̐ݒ� ***
	//ep000=0.05;
	eta_c[1][1]=eta_c[2][2]=eta_c[3][3]=ep000; 			//�i�q�~�X�}�b�`
	eta_c[1][2]=eta_c[2][1]=eta_c[1][3]=eta_c[3][1]=eta_c[2][3]=eta_c[3][2]=0.0;

//***** Fe(bcc)�̒e���萔 ****************************
  //el_mag=1.0E+11/1.0E+09;
  c11=2.33*el_mag;
  c12=1.35*el_mag;
  c44=1.18*el_mag;

	for(i=0;i<=3;i++){
		for(j=0;j<=3;j++){
			for(k=0;k<=3;k++){
				for(l=0;l<=3;l++){
					cec[i][j][k][l]=0.0;//�e���萔�z��̏�����
				}
			}
		}
	}

	cec[1][1][1][1]=c11;
	cec[2][2][2][2]=c11;
	cec[3][3][3][3]=c11;
	cec[1][2][1][2]=cec[1][2][2][1]=cec[2][1][1][2]=cec[2][1][2][1]=c44;
	cec[2][3][2][3]=cec[2][3][3][2]=cec[3][2][2][3]=cec[3][2][3][2]=c44;
	cec[1][3][1][3]=cec[1][3][3][1]=cec[3][1][1][3]=cec[3][1][3][1]=c44;
	cec[1][1][2][2]=cec[2][2][1][1]=c12;
	cec[1][1][3][3]=cec[3][3][1][1]=c12;
	cec[2][2][3][3]=cec[3][3][2][2]=c12;

//--- �A�C�Q�����́i�A�C�Q���Ђ��ݕ������e���ό`�������̉��́j--------------
	sigma[1][1]=cec[1][1][1][1]*eta_c[1][1]
			   +cec[1][1][2][2]*eta_c[2][2]
			   +cec[1][1][3][3]*eta_c[3][3]; //sigma1
	sigma[2][2]=cec[2][2][1][1]*eta_c[1][1]
			   +cec[2][2][2][2]*eta_c[2][2]
			   +cec[2][2][3][3]*eta_c[3][3]; //sigma2
	sigma[3][3]=cec[3][3][1][1]*eta_c[1][1]
			   +cec[3][3][2][2]*eta_c[2][2]
			   +cec[3][3][3][3]*eta_c[3][3]; //sigma3
	sigma[2][3]=cec[2][3][2][3]*eta_c[2][3]+cec[2][3][3][2]*eta_c[3][2]; //sigma4
	sigma[3][2]=cec[3][2][2][3]*eta_c[2][3]+cec[3][2][3][2]*eta_c[3][2]; //sigma4
	sigma[1][3]=cec[1][3][1][3]*eta_c[1][3]+cec[1][3][3][1]*eta_c[3][1]; //sigma5
	sigma[3][1]=cec[3][1][1][3]*eta_c[1][3]+cec[3][1][3][1]*eta_c[3][1]; //sigma5
	sigma[1][2]=cec[1][2][1][2]*eta_c[1][2]+cec[1][2][2][1]*eta_c[2][1]; //sigma6
	sigma[2][1]=cec[2][1][1][2]*eta_c[1][2]+cec[2][1][2][1]*eta_c[2][1]; //sigma6

//*** sin�����cos�e�|�u���A�r�b�g���]�e�[�u���A����я�����̐ݒ� ***************

	table(s, c, ik, ND, ig);
	//ini000();
	datin(ch, ND); //�t�F�[�Y�t�B�[���h�̓Ǎ�

//**** �e�����͂̌v�Z�X�^�[�g ******************************
//start: ;

//**** ��i�t�F�[�Y�t�B�[���h�j�̃t�|���G�ϊ� ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=ch[i][j]; xi[i][j]=0.;
			xr[i*ND+j]=ch[i*ND+j];
			xi[i*ND+j]=0.0;
		}
	}
	qs=-1; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//qrh1[i][j]=xr[i][j]; qih1[i][j]=xi[i][j];
			qrh1[i*ND+j]=xr[i*ND+j];
			qih1[i*ND+j]=xi[i*ND+j];
		}
	}
	//qrh1[0][0]=qih1[0][0]=0.;
	qrh1[0]=qih1[0]=0.0;

//***** �S�Ђ��ݕϓ��ʂ̌v�Z *************************************
//--- ec11�̌v�Z ---
	iii=1; jjj=1;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=qrh1[i][j]*zcij(i, j, iii, jjj)+qih1[i][j]*zcij(i, j, iii, jjj);
			//xi[i][j]=qrh1[i][j]*zcij(i, j, iii, jjj)+qih1[i][j]*zcij(i, j, iii, jjj);
			xr[i*ND+j]=qrh1[i*ND+j]*zcij(i, j, iii, jjj, ND)+qih1[i*ND+j]*zcij(i, j, iii, jjj, ND);
			xi[i*ND+j]=qrh1[i*ND+j]*zcij(i, j, iii, jjj, ND)+qih1[i*ND+j]*zcij(i, j, iii, jjj, ND);
		}
	}
	qs=1; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//ec11[i][j]=xr[i][j];
			ec11[i*ND+j]=xr[i*ND+j];
		}
	}

//--- ec22�̌v�Z ---
	iii=2; jjj=2;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=qrh1[i][j]*zcij(i, j, iii, jjj)+qih1[i][j]*zcij(i, j, iii, jjj);
			//xi[i][j]=qrh1[i][j]*zcij(i, j, iii, jjj)+qih1[i][j]*zcij(i, j, iii, jjj);
			xr[i*ND+j]=qrh1[i*ND+j]*zcij(i, j, iii, jjj, ND)+qih1[i*ND+j]*zcij(i, j, iii, jjj, ND);
			xi[i*ND+j]=qrh1[i*ND+j]*zcij(i, j, iii, jjj, ND)+qih1[i*ND+j]*zcij(i, j, iii, jjj, ND);
		}
	}
	qs=1; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//ec22[i][j]=xr[i][j];
			ec22[i*ND+j]=xr[i*ND+j];
		}
	}

//--- ec12�̌v�Z ---
	iii=1; jjj=2;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=qrh1[i][j]*zcij(i, j, iii, jjj)+qih1[i][j]*zcij(i, j, iii, jjj);
			//xi[i][j]=qrh1[i][j]*zcij(i, j, iii, jjj)+qih1[i][j]*zcij(i, j, iii, jjj);
			xr[i*ND+j]=qrh1[i*ND+j]*zcij(i, j, iii, jjj, ND)+qih1[i*ND+j]*zcij(i, j, iii, jjj, ND);
			xi[i*ND+j]=qrh1[i*ND+j]*zcij(i, j, iii, jjj, ND)+qih1[i*ND+j]*zcij(i, j, iii, jjj, ND);
		}
	}
	qs=1; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//ec12[i][j]=xr[i][j];
			ec12[i*ND+j]=xr[i*ND+j];
		}
	}


//***** �S�Ђ��ݏ�A���͏�A����ђe���c�G�l���M�[��̌v�Z *************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//�S�Ђ��ݏ�̌v�Z[��(5.25)]
			//ep11c[i][j]=ec11[i][j]+ep000*c0;
			//ep22c[i][j]=ec22[i][j]+ep000*c0;
			//ep12c[i][j]=ec12[i][j];
			ep11c[i*ND+j]=ec11[i*ND+j]+ep000*c0;
			ep22c[i*ND+j]=ec22[i*ND+j]+ep000*c0;
			ep12c[i*ND+j]=ec12[i*ND+j];

			//�e�����͏�̌v�Z[��(5.27)]
			//sig11[i][j]=c11*ec11[i][j]+c12*ec22[i][j]-(c11+2.*c12)*ep000*(ch[i][j]-c0);
			//sig22[i][j]=c12*ec11[i][j]+c11*ec22[i][j]-(c11+2.*c12)*ep000*(ch[i][j]-c0);
			//sig33[i][j]=c12*ec11[i][j]+c12*ec22[i][j]-(c11+2.*c12)*ep000*(ch[i][j]-c0);
			//sig12[i][j]=2.*c44*ec12[i][j];
			sig11[i*ND+j]=c11*ec11[i*ND+j]+c12*ec22[i*ND+j]-(c11+2.0*c12)*ep000*(ch[i*ND+j]-c0);
			sig22[i*ND+j]=c12*ec11[i*ND+j]+c11*ec22[i*ND+j]-(c11+2.0*c12)*ep000*(ch[i*ND+j]-c0);
			sig33[i*ND+j]=c12*ec11[i*ND+j]+c12*ec22[i*ND+j]-(c11+2.0*c12)*ep000*(ch[i*ND+j]-c0);
			sig12[i*ND+j]=2.0*c44*ec12[i*ND+j];

			//�e���Ђ��݃G�l���M�[��̌v�Z[��(5.28)]
			//Estr1=1.5*(c11+2.0*c12)*ep000*ep000*ch[i][j]*ch[i][j];
			//Estr2=-(c11+2.0*c12)*(ep11c[i][j]+ep22c[i][j])*ep000*ch[i][j];
			//Estr3=0.5*c11*(ep11c[i][j]*ep11c[i][j]+ep22c[i][j]*ep22c[i][j])
            //+c12*ep11c[i][j]*ep22c[i][j]+2.*c44*ep12c[i][j]*ep12c[i][j];
			//Estr[i][j]=Estr1+Estr2+Estr3;
			Estr1=1.5*(c11+2.0*c12)*ep000*ep000*ch[i*ND+j]*ch[i*ND+j];
			Estr2=-(c11+2.0*c12)*(ep11c[i*ND+j]+ep22c[i*ND+j])*ep000*ch[i*ND+j];
			Estr3=0.5*c11*(ep11c[i*ND+j]*ep11c[i*ND+j]+ep22c[i*ND+j]*ep22c[i*ND+j])
                      +c12*ep11c[i*ND+j]*ep22c[i*ND+j]+2.0*c44*ep12c[i*ND+j]*ep12c[i*ND+j];
			Estr[i*ND+j]=Estr1+Estr2+Estr3;
		}
	}


//***** �ψʏ�̌v�Z[��(5.31)] *************************************
//--- u1�̌v�Z ---
	iii=1;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=-qrh1[i][j]*zuij(i, j, iii)-qih1[i][j]*zuij(i, j, iii);
			//xi[i][j]=qrh1[i][j]*zuij(i, j, iii)+qih1[i][j]*zuij(i, j, iii);
			xr[i*ND+j]=-qrh1[i*ND+j]*zuij(i, j, iii, ND)-qih1[i*ND+j]*zuij(i, j, iii, ND);
			xi[i*ND+j]= qrh1[i*ND+j]*zuij(i, j, iii, ND)+qih1[i*ND+j]*zuij(i, j, iii, ND);
		}
	}
	qs=1; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//u1[i][j]=xr[i][j];
			u1[i*ND+j]=xr[i*ND+j];
		}
	}

//--- u2�̌v�Z ---
	iii=2;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=-qrh1[i][j]*zuij(i, j, iii)-qih1[i][j]*zuij(i, j, iii);
			//xi[i][j]=qrh1[i][j]*zuij(i, j, iii)+qih1[i][j]*zuij(i, j, iii);
			xr[i*ND+j]=-qrh1[i*ND+j]*zuij(i, j, iii, ND)-qih1[i*ND+j]*zuij(i, j, iii, ND);
			xi[i*ND+j]= qrh1[i*ND+j]*zuij(i, j, iii, ND)+qih1[i*ND+j]*zuij(i, j, iii, ND);
		}
	}
	qs=1; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//u2[i][j]=xr[i][j];
			u2[i*ND+j]=xr[i*ND+j];
		}
	}

//--- u3�̌v�Z ---
//	iii=3;
//	for(i=0;i<=ndm;i++){
//		for(j=0;j<=ndm;j++){
//			//xr[i][j]=-qrh1[i][j]*zuij(i, j, iii)-qih1[i][j]*zuij(i, j, iii);
//			//xi[i][j]=qrh1[i][j]*zuij(i, j, iii)+qih1[i][j]*zuij(i, j, iii);
//			xr[i*ND+j]=-qrh1[i*ND+j]*zuij(i, j, iii, ND)-qih1[i*ND+j]*zuij(i, j, iii, ND);
//			xi[i*ND+j]= qrh1[i*ND+j]*zuij(i, j, iii, ND)+qih1[i*ND+j]*zuij(i, j, iii, ND);
//		}
//	}
//	qs=1; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
//	for(i=0;i<=ndm;i++){
//		for(j=0;j<=ndm;j++){
//			//u3[i][j]=xr[i][j];
//			u3[i*ND+j]=xr[i*ND+j];
//		}
//	}


//***** �S���l�f�[�^�ۑ� ********************************************
	//datsave(ch, Estr, ep11c, ep22c, ep12c,
	//				sig11, sig22, sig33, sig12, u1, u2, ND);	//save
	datsave_paraview(ch, Estr, ep11c, ep22c, ep12c,
					sig11, sig22, sig33, sig12, u1, u2, ND);	//save
	//if(keypress()){return 0;}//�L�[�҂����
	return 0;
}


//******* Sin, Cos �̃e�[�u������уr�b�g���]�e�[�u���̐ݒ� ***************
void table(double *s, double *c, int *ik, int ND, int ig)
{
	int i, j, mc, mn;
	double q;
	int nd=ND, nd2=ND/2;

	q=2.0*PI/nd;
	for(i=0;i<=nd2-1;i++){ c[i]=cos(q*i); s[i]=sin(q*i); }//Sin, Cos �̃e�[�u��

	ik[0]=0; mn=nd2; mc=1;
	for(i=1;i<=ig;i++){
		for(j=0;j<=mc-1;j++){ ik[j+mc]=ik[j]+mn; }	//�r�b�g���]�e�[�u��
		mc*=2; mn/=2; 
	}
}

//********** �P���������t�[���G�ϊ� **************************************
void fft(double *xrf, double *xif, double *s, double *c, int ND, int ig)
{
	int lg, lf, mf, nf, n2, ka, kb, ix;
	double tj, tr;
	int nd2=ND/2;

	lg=1;
	for(lf=1;lf<=ig;lf++){
		n2=nd2/lg;
		for(mf=1;mf<=lg;mf++){
			for(nf=0;nf<=n2-1;nf++){
				ix=nf*lg;
				ka=nf+2*n2*(mf-1);
				kb=ka+n2;
				tr=xrf[ka]-xrf[kb];
				tj=xif[ka]-xif[kb];
				xrf[ka]+=xrf[kb];
				xif[ka]+=xif[kb];
				xrf[kb]=tr*c[ix]-tj*float(qs)*s[ix];
				xif[kb]=tj*c[ix]+tr*float(qs)*s[ix];
			}
		}
		lg*=2;
	}//lf
}

//************ �Q���������t�[���G�ϊ� ***********************************
void rcfft(double *xrf, double *xif, double *xr, double *xi, int *ik, double *s, double *c, int ND, int ig)
{
	int i, ic, ir, j;
	int nd=ND, ndm=ND-1;

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

	if(qs > 0){return;}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=xr[i][j]/nd/nd;  xi[i][j]=xi[i][j]/nd/nd;
			xr[i*ND+j]=xr[i*ND+j]/nd/nd;
			xi[i*ND+j]=xi[i*ND+j]/nd/nd;
		}
	}

}

//*** Zcij [��(5.26)�̌v�Z] ****************************************
double zcij(int i0, int j0, int iii, int jjj, int ND)
{
	//int i, j, k, m, n, p, q;
	int m, n;
	int ii, jj;
	double nx, ny, nz, alnn;
	double nec[4];
	double a11, a22, a33, a12, a13, a23, a21, a31, a32;
	double b11, b22, b33, b12, b13, b23, b21, b31, b32;
	double det1, zij;
	double om[4][4];
	int nd=ND, nd2=ND/2;

	if(i0<=nd2-1) {ii=i0;}
	if(i0>=nd2) {ii=i0-nd;}
	if(j0<=nd2-1) {jj=j0;}
	if(j0>=nd2) {jj=j0-nd;}
	alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj);
	if(alnn==0.0){alnn=1.0;}
	nec[1]=nx=(double)ii/alnn;
	nec[2]=ny=(double)jj/alnn;
	nec[3]=nz=0.0;

  //for(i=1;i<=3;i++){ for(j=1;j<=3;j++){ om[i][j]=0.; } }

	a11=cec[1][1][1][1]*nx*nx+cec[1][2][1][2]*ny*ny+cec[1][3][1][3]*nz*nz;
	a22=cec[1][2][1][2]*nx*nx+cec[2][2][2][2]*ny*ny+cec[2][3][2][3]*nz*nz;
	a33=cec[3][1][3][1]*nx*nx+cec[2][3][2][3]*ny*ny+cec[3][3][3][3]*nz*nz;
	a12=(cec[1][1][2][2]+cec[1][2][1][2])*nx*ny;
	a23=(cec[2][2][3][3]+cec[2][3][2][3])*ny*nz;
	a31=(cec[3][3][1][1]+cec[3][1][3][1])*nx*nz;
	a21=a12;
	a32=a23;
	a13=a31;

	b11=a22*a33-a23*a32;
	b22=a11*a33-a13*a31;
	b33=a11*a22-a12*a21;
	b12=-(a21*a33-a23*a31);
	b23=-(a11*a32-a12*a31);
	b31=-(a22*a31-a21*a32);
	b21=b12;
	b32=b23;
	b13=b31;

	det1=a11*a22*a33+a12*a23*a31+a13*a32*a21-a13*a31*a22-a11*a23*a32-a33*a12*a21;
	if(det1==0.0){det1=1.0;}

	om[1][1]=b11/det1;
	om[2][2]=b22/det1;
	om[3][3]=b33/det1;
	om[1][2]=om[2][1]=b12/det1;
	om[2][3]=om[3][2]=b23/det1;
	om[3][1]=om[1][3]=b31/det1;

	zij=0.0;
	for(m=1;m<=3;m++) {
		for(n=1;n<=3;n++) {
			zij=zij+0.5*(sigma[m][n]*nec[jjj]*nec[n]*om[m][iii]
						+sigma[m][n]*nec[iii]*nec[n]*om[m][jjj]);
		}
   }
	return(zij);
}

//*** Zuij  [��(5.30)�̌v�Z] ****************************************
double zuij(int i0, int j0, int iii, int ND)
{
	//int i, j, k, m, n, p, q;
	int m, n;
	int ii, jj;
	double nx, ny, nz, alnn;
	double nec[4];
	double a11, a22, a33, a12, a13, a23, a21, a31, a32;
	double b11, b22, b33, b12, b13, b23, b21, b31, b32;
	double det1, zij;
	double om[4][4];
	int nd=ND, nd2=ND/2;

	if(i0<=nd2-1) {ii=i0;}
	if(i0>=nd2) {ii=i0-nd;}
	if(j0<=nd2-1) {jj=j0;}
	if(j0>=nd2) {jj=j0-nd;}
	alnn=sqrt((double)ii*(double)ii+(double)jj*(double)jj);
	if(alnn==0.0){alnn=1.0;}
	nec[1]=nx=(double)ii/alnn;
	nec[2]=ny=(double)jj/alnn;
	nec[3]=nz=0.0;

  //for(i=1;i<=3;i++){ for(j=1;j<=3;j++){ om[i][j]=0.; } }

	a11=cec[1][1][1][1]*nx*nx+cec[1][2][1][2]*ny*ny+cec[1][3][1][3]*nz*nz;
	a22=cec[1][2][1][2]*nx*nx+cec[2][2][2][2]*ny*ny+cec[2][3][2][3]*nz*nz;
	a33=cec[3][1][3][1]*nx*nx+cec[2][3][2][3]*ny*ny+cec[3][3][3][3]*nz*nz;
	a12=(cec[1][1][2][2]+cec[1][2][1][2])*nx*ny;
	a23=(cec[2][2][3][3]+cec[2][3][2][3])*ny*nz;
	a31=(cec[3][3][1][1]+cec[3][1][3][1])*nx*nz;
	a21=a12;
	a32=a23;
	a13=a31;

	b11=a22*a33-a23*a32;
	b22=a11*a33-a13*a31;
	b33=a11*a22-a12*a21;
	b12=-(a21*a33-a23*a31);
	b23=-(a11*a32-a12*a31);
	b31=-(a22*a31-a21*a32);
	b21=b12;
	b32=b23;
	b13=b31;

	det1=a11*a22*a33+a12*a23*a31+a13*a32*a21-a13*a31*a22-a11*a23*a32-a33*a12*a21;
	if(det1==0.0){det1=1.0;}

	om[1][1]=b11/det1;
	om[2][2]=b22/det1;
	om[3][3]=b33/det1;
	om[1][2]=om[2][1]=b12/det1;
	om[2][3]=om[3][2]=b23/det1;
	om[3][1]=om[1][3]=b31/det1;

	zij=0.0;
	for(m=1;m<=3;m++) {
		for(n=1;n<=3;n++) {
			zij=zij+sigma[m][n]*nec[n]*om[m][iii];
		}
	}
	zij=-zij/alnn;
	return(zij);
}

//************ �e��̏�̃f�[�^�̕ۑ� *******************************
void datsave(double *ch, double *Estr, 
			 double *ep11c, double *ep22c, double *ep12c,
			 double *sig11, double *sig22, double *sig33, double *sig12, 
			 double *u1, double *u2, int ND)
{
	FILE		*stream;//�X�g���[���̃|�C���^�ݒ�
	int 		i, j;//����
	int ndm=ND-1;

	stream = fopen("el_field.dat", "w");//�������ݐ�̃t�@�C�����㏑�������ŃI�[�v��

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fprintf(stream, "%e  ", ch[i*ND+j]);//�e���̃f�[�^�ۑ�
			fprintf(stream, "%e  ", Estr[i*ND+j]);
			fprintf(stream, "%e  ", ep11c[i*ND+j]);
			fprintf(stream, "%e  ", ep22c[i*ND+j]);
			fprintf(stream, "%e  ", ep12c[i*ND+j]);
			fprintf(stream, "%e  ", sig11[i*ND+j]);
			fprintf(stream, "%e  ", sig22[i*ND+j]);
			fprintf(stream, "%e  ", sig33[i*ND+j]);
			fprintf(stream, "%e  ", sig12[i*ND+j]);
			fprintf(stream, "%e  ", u1[i*ND+j]);
			fprintf(stream, "%e  ", u2[i*ND+j]);
		}
	}
	fprintf(stream, "\n");//���s�̏�������
	fclose(stream);//�t�@�C�����N���[�Y
}

void datsave_paraview(double *ch, double *Estr, 
			 double *ep11c, double *ep22c, double *ep12c,
			 double *sig11, double *sig22, double *sig33, double *sig12, 
			 double *u1, double *u2, int ND)
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
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),1);
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %11d \n",((ndm+1)*(ndm+1)*1));
	fprintf(fp,"SCALARS concentration float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", ch[i*ND+j]);//�e���̃f�[�^�ۑ�
		}
	}
	fprintf(fp,"SCALARS Elastic_strain_energy_density float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", Estr[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Strain_c11 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", ep11c[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Strain_c22 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", ep22c[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Strain_c12 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", ep12c[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_11 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", sig11[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_22 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", sig22[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_33 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", sig33[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Elastic_stress_12 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", sig12[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Displacement_u1 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", u1[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Displacement_u2 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", u2[i*ND+j]);
		}
	}
	fprintf(fp,"VECTORS vectors float \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f %10.6f %10.6f\n", u1[i*ND+j], u2[i*ND+j], 0.0); // �ψʏ����ɂĕ`��
		}
	}
	fclose(fp);
}
//************ ��̃f�[�^�̓ǂݍ��� *****************************************
void datin(double *ch, int ND)
{
	FILE		*datin0;//�X�g���[���̃|�C���^�ݒ�
	int 		i, j;//����
	double c00;//��̕��ϒl
	int nd=ND, ndm=ND-1;

	datin0 = fopen("ph.dat", "r");//�Ǎ��݌��̃t�@�C�����I�[�v��
	c00=0.0;//��̕��ϒl�̏����l
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fscanf(datin0, "%lf", &ch[i*ND+j]);//��̃f�[�^�Ǎ�
			c00+=ch[i*ND+j];
		}
	}
	c0=c00/nd/nd;//��̕��ϒl
	printf("c0=  %f  \n", c0);

	fclose(datin0);//�t�@�C�����N���[�Y
}
