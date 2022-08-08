#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//�����̊֐��ݒ�
//#define ND 256					//�����v�Z�ɂ�����v�Z�̈��ӂ̕�����(�����t�[���G�ϊ���p���邽�߂Q�ׂ̂���)
//#define IG 8					//2^IG=ND

	//int nd=ND, ndm=ND-1; 		//�v�Z�̈�̈�ӂ̍���������(�����u���b�N��)�AND-1���`
	//int nd2=ND/2;				 //ND/2���`�F�����t�|���G�ϊ��Ŏg�p
	//int ig=IG;					//2^ig=ND
	double PI=3.14159;			//�~����
	double rr=8.3145;			//�K�X�萔
	double time1;				//�v�Z�J�E���g��(���Ԃɔ��)
	int Nstep, iout;

	//double s1h[ND][ND], s2h[ND][ND];	//�}���e���T�C�g�̃t�F�[�Y�t�B�[���h

	double qs;					//�t�|���G�ϊ�(qs:-1)�ƃt�|���G�t�ϊ�(qs:1)�̋��
	//double xi[ND][ND], xr[ND][ND], xif[ND], xrf[ND];//�t�|���G�ϊ��̎����E�����z��
	//double s[ND],c[ND];			//sin��cos�̃e�[�u��
	//int ik[ND];					//�r�b�g���]�e�[�u��

	void ini000(double *s1h, double *s2h, int ND);	//������̐ݒ�T�u���|�`��
	void table(double *s, double *c, int *ik, int ND, int ig);//sin��cos�̃e�[�u���ƃr�b�g���]�e�[�u���̍쐬�T�u���|�`��
	void fft(double *xrf, double *xif, double *s, double *c, int ND, int ig);//�P���������t�[���G�ϊ�
	void rcfft(double *xrf, double *xif, double *xr, double *xi, int *ik, double *s, double *c, int ND, int ig);//�Q���������t�[���G�ϊ�
	void datsave(double *s1h, double *s2h, int ND);			//�f�|�^�ۑ��T�u���|�`��
	void datsave_paraview(double *s1h, double *s2h, int ND);//�f�|�^�ۑ��T�u���|�`��

//******* ���C���v���O���� ******************************************
int main(void)
{
    int ND, IG;
	int nd, ndm, nd2, ig;
	
	double s1, s2;									//�}���e���T�C�g�̃t�F�[�Y�t�B�[���h
	//double ep11h0[ND][ND], ep22h0[ND][ND];			//�g�D���̕ϑԘc
	//double ep11qrh0[ND][ND],	ep11qih0[ND][ND];	//�S���c�ϓ��ʂ̃t�[���G�ϊ�
	//double ep22qrh0[ND][ND],	ep22qih0[ND][ND];	//�S���c�ϓ��ʂ̃t�[���G�ϊ�
	double s1k_chem, s1k_str;						//�|�e���V����
	//double s1k_su[ND][ND];							//�|�e���V����
	double s2k_chem, s2k_str;						//�|�e���V����
	//double s2k_su[ND][ND];							//�|�e���V����
	double c11, c12, c44, lam0, mu0, nu0; 			//�e���萔
	double c11c, c12c, c44c;			 			//�e���萔
	double eta_s1[4][4], eta_s2[4][4];				//�A�C�Q���c����
	//double ec11[ND][ND], ec22[ND][ND];				//�S���c�ϓ��ʁi����ԁj
	double ep11T, ep22T;
	double ep11_0, ep22_0;							//�g�D���̕ϑԘc�̕��ϒl
	double ep11_a, ep22_a, ep12_a, ep21_a;			//�O�͂ɋN������c
	double sig11_a, sig22_a;						//�O��
	double Z11ep, Z12ep, Z21ep, Z22ep;				//�t�[���G�t�ϊ����̌W��
	double sum11, sum22;							//s1��s2�̋�Ԑϕ�
	double s1ddtt, s2ddtt;							//s1��s2�̎��ԕω��ʁi���W�������̍��Ӂj
	double el_fac;									//�e���萔�̋K�i���萔

	int   i, j, k, l, ii, jj, kk, iii, jjj;			//����
	int   p, q, m, n;								//����
	int   ip, im, jp, jm, Nstep;					//����
	double al, temp, delt;							//�v�Z�̈�A���x�A���Ԃ�����
	double time1max;								//�v�Z�J�E���g���̍ő�l�i�v�Z�I���J�E���g�j
	double b1, vm0, atom_n;							//�����v���b�N�P�ӂ̒����A�����̐ρA�P�ʖE���̌��q��
	double smob;									//���r���e�B�[�i�����ϑԂ̊ɘa�W���j
	double nxx, nyy, nxy, alnn;						//�t�[���G��Ԃ̊�{�x�N�g���̐ρA�m����

	double AA0, AA1, AA2, AA3;						//�M�u�X�G�l���M�[���̌W��
	double AA0e;									//�M�u�X�G�l���M�[���̌W��
	double a1_c, b1_c, c1_c;						//�i�q�萔
	double a1_t, b1_t, c1_t;						//�i�q�萔
	double kappa_s1, kappa_s2;						//���z�G�l���M�|�W��
	double kappa_s1c, kappa_s2c;					//���z�G�l���M�|�W��
	double ds_fac;									//�����ϑԂ̗h�炬�W��
	double eta1, eta2;

//****** �v�Z��������ѕ����萔�̐ݒ� ****************************************
	printf("---------------------------------\n");
	printf("read parameters form parameters.txt\n");
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
	delt    = data[1];
	temp    = data[2];
	al      = data[3];	// [nm]
	smob    = data[4];
	ds_fac  = data[5];
	AA0e    = data[6];
	AA1     = data[7];
	AA2     = data[8];
	AA3     = data[9];
	kappa_s1c=data[10];
	kappa_s2c=data[11];
	vm0     = data[12];
	time1max= int(data[13]);
	Nstep   = int(data[14]);
	c11c    = data[15];
	c12c    = data[16];
	c44c    = data[17];
	lam0    = data[18];
	mu0     = data[19];
	eta1    = data[20];
	eta2    = data[21];
	sig22_a = data[22];
	printf("---------------------------------\n");
	//
	IG=int(log2(ND));
	nd=ND; 					
	ndm=ND-1;				//define ND-1
	nd2=ND/2;				//define ND/2 for FFT
	ig=IG;					//2^ig=ND
	//
	double *s1h      = (double *)malloc(sizeof(double)*( ND*ND + ND ));//�}���e���T�C�g�̃t�F�[�Y�t�B�[���h
	double *s2h      = (double *)malloc(sizeof(double)*( ND*ND + ND ));//�}���e���T�C�g�̃t�F�[�Y�t�B�[���h
	//
	double *ep11h0   = (double *)malloc(sizeof(double)*( ND*ND + ND ));//�g�D���̕ϑԘc
	double *ep22h0   = (double *)malloc(sizeof(double)*( ND*ND + ND ));//�g�D���̕ϑԘc
	double *ep11qrh0 = (double *)malloc(sizeof(double)*( ND*ND + ND ));//�S���c�ϓ��ʂ̃t�[���G�ϊ�
	double *ep11qih0 = (double *)malloc(sizeof(double)*( ND*ND + ND ));//�S���c�ϓ��ʂ̃t�[���G�ϊ�
	double *ep22qrh0 = (double *)malloc(sizeof(double)*( ND*ND + ND ));//�S���c�ϓ��ʂ̃t�[���G�ϊ�
	double *ep22qih0 = (double *)malloc(sizeof(double)*( ND*ND + ND ));//�S���c�ϓ��ʂ̃t�[���G�ϊ�
	//
	double *s1k_su   = (double *)malloc(sizeof(double)*( ND*ND + ND ));//�|�e���V����
	double *s2k_su   = (double *)malloc(sizeof(double)*( ND*ND + ND ));//�|�e���V����
	//
	double *ec11     = (double *)malloc(sizeof(double)*( ND*ND + ND ));//�S���c�ϓ��ʁi����ԁj
	double *ec22     = (double *)malloc(sizeof(double)*( ND*ND + ND ));//�S���c�ϓ��ʁi����ԁj
	//
	double *xi       = (double *)malloc(sizeof(double)*( ND*ND + ND ));//�t�|���G�ϊ��̎����E�����z��
	double *xr       = (double *)malloc(sizeof(double)*( ND*ND + ND ));//�t�|���G�ϊ��̎����E�����z��
	//
	double *xif      = (double *)malloc(sizeof(double)*( ND ));//�t�|���G�ϊ��̎����E�����z��
	double *xrf      = (double *)malloc(sizeof(double)*( ND ));//�t�|���G�ϊ��̎����E�����z��
	//
	double *s        = (double *)malloc(sizeof(double)*( ND ));//sin�̃e�[�u��
	double *c        = (double *)malloc(sizeof(double)*( ND ));//cos�̃e�[�u��
	//
	int *ik       = (int *)malloc(sizeof(int)*( ND ));//�r�b�g���]�e�[�u��
	//
	//printf("DELT(0.2)=  ");	scanf(" %lf",&delt);	//���Ԃ����ݓ���	//delt=0.2;

	//temp=500.0;					//���x(K)
	//al=500.0*1.0E-09;			//�v�Z�̈�(m)
	al=al*1.0E-09;				//�v�Z�̈�(m)
	b1=al/nd;					//�����u���b�N�̒���

	time1=0.0;					//�����v�Z�J�E���g���̐ݒ�
	//time1max=1.0+1.0e+07;		//�ő�v�Z�J�E���g���̐ݒ�

	//smob=1.0;					//���r���e�B�[�i�����ϑԂ̊ɘa�W���j
	//ds_fac=0.01;				//�����ϑԂ̗h�炬�W��

	//AA0=1000.0/rr/temp;			//�}���e���T�C�g�ϑԂ̉��w�I�쓮��
	AA0=AA0e/rr/temp;			//�}���e���T�C�g�ϑԂ̉��w�I�쓮��
	//AA1=1.0;  AA2=3.0*AA1+12.0;  AA3=2.0*AA1+12.0;	//�M�u�X�G�l���M�[���̌W��

	//kappa_s1=kappa_s2=5.0e-15/rr/temp/b1/b1;	//���z�G�l���M�|�W��
	kappa_s1=kappa_s1c/rr/temp/b1/b1;	//���z�G�l���M�|�W��
	kappa_s2=kappa_s2c/rr/temp/b1/b1;	//���z�G�l���M�|�W��

	//a1_c=b1_c=c1_c=3.5E-10;						//�i�q�萔(m)
	//atom_n=4.0;  vm0=6.02E23*a1_c*b1_c*c1_c/atom_n;	//�����̐ς̌v�Z�ifcc������j

//*** s1��̃A�C�Q���c�̐ݒ� ***************
	//eta_s1[1][1]=0.05; eta_s1[2][2]=-0.05;
	eta_s1[1][1]=eta1; eta_s1[2][2]=eta2;
	eta_s1[3][3]=0.;
	eta_s1[1][2]=eta_s1[2][1]=eta_s1[1][3]=eta_s1[3][1]=eta_s1[2][3]=eta_s1[3][2]=0.;

//*** s2��̃A�C�Q���c�̐ݒ� ***************
	eta_s2[1][1]=eta_s1[2][2];
	eta_s2[2][2]=eta_s1[1][1];
	eta_s2[3][3]=0.;
	eta_s2[1][2]=eta_s2[2][1]=eta_s2[1][3]=eta_s2[3][1]=eta_s2[2][3]=eta_s2[3][2]=0.;

//***** Ni�̒e���萔 ****************************
	el_fac=1.0E+11*vm0/rr/temp;
	//c11=2.508*el_fac;
	//c44=1.235*el_fac;
	//c12=1.500*el_fac;
	c11=c11c*el_fac;
	c12=c12c*el_fac;
	c44=c44c*el_fac;
	//c12=c11-2.0*c44;
	if(lam0==0.0){
		lam0=c12;//cij�̃f�[�^���烉�[���萔��ݒ�
	}
	if(mu0==0.0){
		mu0=c44;//cij�̃f�[�^���烉�[���萔��ݒ�
	}
	nu0=lam0/2.0/(lam0+mu0);//�|�A�\����
	printf("nu0= %f  \n", nu0);//�|�A�\����̒l��\���i���[���萔�̑Ó����̊m�F�̂��߁j

//*** �O�͂̐ݒ� *******************************
 	//sig22_a=0.0;//�{�v�Z�ł́A�O�͂��l�����Ă��Ȃ��̂łO��ݒ�
	ep11_a=-0.5*lam0/mu0/(3.0*lam0+2.0*mu0)*sig22_a;
	ep22_a=(lam0+mu0)/mu0/(3.0*lam0+2.0*mu0)*sig22_a;
	ep12_a=ep21_a=0.0;

//*** sin�����cos�e�|�u���A�r�b�g���]�e�[�u���A����я�����̐ݒ� ***************

	table(s, c, ik, ND, ig);	//sin�����cos�e�|�u���ƃr�b�g���]�e�[�u���̐ݒ�

	ini000(s1h, s2h, ND);		//������̐ݒ�

//**** �V�~�����[�V�����X�^�[�g ******************************
//Nstep = 10;
iout = -1;
start: ;

	//if(time1<=100.){Nstep=10;} else{Nstep=200;}		//�f�[�^�ۑ����鎞�ԊԊu�̕ύX
	//if((((int)(time1) % Nstep)==0)) {datsave(s1h, s2h, ND);} 	//���J�Ԃ��J�E���g���ɑg�D�f�[�^��ۑ�
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(s1h, s2h, ND);} 	//���J�Ԃ��J�E���g���ɑg�D�f�[�^��ۑ�

//***** ���z�|�e���V���� ***********************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}	//�����I���E����
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}
			//s1k_su[i][j]=-kappa_s1*(s1h[ip][j]+s1h[im][j]+s1h[i][jp]+s1h[i][jm]-4.0*s1h[i][j]);//��(4.2.4)
			//s2k_su[i][j]=-kappa_s2*(s2h[ip][j]+s2h[im][j]+s2h[i][jp]+s2h[i][jm]-4.0*s2h[i][j]);
			s1k_su[i*ND+j]=-kappa_s1*(s1h[ip*ND+j]+s1h[im*ND+j]+s1h[i*ND+jp]+s1h[i*ND+jm]-4.0*s1h[i*ND+j]);//��(4.2.4)
			s2k_su[i*ND+j]=-kappa_s2*(s2h[ip*ND+j]+s2h[im*ND+j]+s2h[i*ND+jp]+s2h[i*ND+jm]-4.0*s2h[i*ND+j]);
		}
	}

//**** �A�C�Q���c��[��(4.7)]�̃t�|���G�ϊ� ep11 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=ep11h0[i][j]=eta_s1[1][1]*s1h[i][j]+eta_s2[1][1]*s2h[i][j];
			//xi[i][j]=0.0;
			xr[i*ND+j]=ep11h0[i*ND+j]=eta_s1[1][1]*s1h[i*ND+j]+eta_s2[1][1]*s2h[i*ND+j];
			xi[i*ND+j]=0.0;
		}
	}
	qs=-1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){ 
			//ep11qrh0[i][j]=xr[i][j];
			//ep11qih0[i][j]=xi[i][j];
			ep11qrh0[i*ND+j]=xr[i*ND+j];
			ep11qih0[i*ND+j]=xi[i*ND+j];
		}
	}
	//ep11qrh0[0][0]=ep11qih0[0][0]=0.;
	ep11qrh0[0]=ep11qih0[0]=0.0;

//**** �A�C�Q���c��[��(4.7)]�̃t�|���G�ϊ� ep22 ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=ep22h0[i][j]=eta_s1[2][2]*s1h[i][j]+eta_s2[2][2]*s2h[i][j];
			//xi[i][j]=0.0;
			xr[i*ND+j]=ep22h0[i*ND+j]=eta_s1[2][2]*s1h[i*ND+j]+eta_s2[2][2]*s2h[i*ND+j];
			xi[i*ND+j]=0.0;
		}
	}
	qs=-1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){ 
			//ep22qrh0[i][j]=xr[i][j];
			//ep22qih0[i][j]=xi[i][j];
			ep22qrh0[i*ND+j]=xr[i*ND+j];
			ep22qih0[i*ND+j]=xi[i*ND+j];
		}
	}
	//ep22qrh0[0][0]=ep22qih0[0][0]=0.;
	ep22qrh0[0]=ep22qih0[0]=0.0;

//*** �A�C�Q���c��̕��ϒl�̎Z�o ***
	sum11=sum22=0.;
	for(i=0;i<=ndm;i++){
		//for(j=0;j<=ndm;j++){ sum11+=ep11h0[i][j];  sum22+=ep22h0[i][j]; }
		for(j=0;j<=ndm;j++){ sum11+=ep11h0[i*ND+j];  sum22+=ep22h0[i*ND+j]; }
	}
	ep11_0=sum11/nd/nd;  ep22_0=sum22/nd/nd;

//***** �S���c�ϓ���ec11�̌v�Z[��(4.9)] *************************************
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
			//xr[i][j]=Z11ep*ep11qrh0[i][j]+Z12ep*ep22qrh0[i][j];
			//xi[i][j]=Z11ep*ep11qih0[i][j]+Z12ep*ep22qih0[i][j];
			xr[i*ND+j]=Z11ep*ep11qrh0[i*ND+j]+Z12ep*ep22qrh0[i*ND+j];
			xi[i*ND+j]=Z11ep*ep11qih0[i*ND+j]+Z12ep*ep22qih0[i*ND+j];
		}
	}
	qs=1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//ec11[i][j]=xr[i][j];
			ec11[i*ND+j]=xr[i*ND+j];
		}
	}

//***** �S���c�ϓ���ec22�̌v�Z[��(4.9)] *****************************
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
			//xr[i][j]=Z21ep*ep11qrh0[i][j]+Z22ep*ep22qrh0[i][j];
			//xi[i][j]=Z21ep*ep11qih0[i][j]+Z22ep*ep22qih0[i][j];
			xr[i*ND+j]=Z21ep*ep11qrh0[i*ND+j]+Z22ep*ep22qrh0[i*ND+j];
			xi[i*ND+j]=Z21ep*ep11qih0[i*ND+j]+Z22ep*ep22qih0[i*ND+j];
		}
	}
	qs=1.; rcfft(xrf, xif, xr, xi, ik, s, c, ND, ig);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//ec22[i][j]=xr[i][j];
			ec22[i*ND+j]=xr[i*ND+j];
		}
	}

//******  �|�e���V�����̌v�Z ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){

			//s1=s1h[i][j];  	s2=s2h[i][j];
			s1=s1h[i*ND+j];  	s2=s2h[i*ND+j];

//******  ���w�|�e���V�����̌v�Z[��(4.4)] ********************************
			s1k_chem=AA0*s1*(AA1-AA2*s1+AA3*(s1*s1+s2*s2));
			s2k_chem=AA0*s2*(AA1-AA2*s2+AA3*(s1*s1+s2*s2));

//******  �e���|�e���V�����̌v�Z[��(4.8)] ********************************

			//ep11T=ep11h0[i][j]-ep11_0-ec11[i][j]-ep11_a;
			//ep22T=ep22h0[i][j]-ep22_0-ec22[i][j]-ep22_a;
			ep11T=ep11h0[i*ND+j]-ep11_0-ec11[i*ND+j]-ep11_a;
			ep22T=ep22h0[i*ND+j]-ep22_0-ec22[i*ND+j]-ep22_a;

			s1k_str=ep11T*((lam0+2.0*mu0)*eta_s1[1][1]+lam0*eta_s1[2][2])
				   +ep22T*((lam0+2.0*mu0)*eta_s1[2][2]+lam0*eta_s1[1][1]);
			s2k_str=ep11T*((lam0+2.0*mu0)*eta_s2[1][1]+lam0*eta_s2[2][2])
				   +ep22T*((lam0+2.0*mu0)*eta_s2[2][2]+lam0*eta_s2[1][1]);

//****** �t�F�[�Y�t�B�[���h�̎��Ԕ��W�̌v�Z[��(4.10)] ********************************
			//s1ddtt=-smob*(s1k_chem+s1k_su[i][j]+s1k_str);
			//s2ddtt=-smob*(s2k_chem+s2k_su[i][j]+s2k_str);
			//s1h[i][j]=s1h[i][j]+( s1ddtt+ds_fac*(2.0*DRND(1.)-1.0) )*delt;//�z��@
			//s2h[i][j]=s2h[i][j]+( s2ddtt+ds_fac*(2.0*DRND(1.)-1.0) )*delt;
			s1ddtt=-smob*(s1k_chem+s1k_su[i*ND+j]+s1k_str);
			s2ddtt=-smob*(s2k_chem+s2k_su[i*ND+j]+s2k_str);
			s1h[i*ND+j]=s1h[i*ND+j]+( s1ddtt+ds_fac*(2.0*DRND(1.)-1.0) )*delt;//�z��@
			s2h[i*ND+j]=s2h[i*ND+j]+( s2ddtt+ds_fac*(2.0*DRND(1.)-1.0) )*delt;

//*** s�̕ψ�(0<=s<=1)�̕␳ ***
			//if(s1h[i][j]>=1.0){s1h[i][j]=1.0;}  if(s1h[i][j]<=0.0){s1h[i][j]=0.0;}
			//if(s2h[i][j]>=1.0){s2h[i][j]=1.0;}  if(s2h[i][j]<=0.0){s2h[i][j]=0.0;}
			if(s1h[i*ND+j]>=1.0){s1h[i*ND+j]=1.0;}
			if(s1h[i*ND+j]<=0.0){s1h[i*ND+j]=0.0;}
			if(s2h[i*ND+j]>=1.0){s2h[i*ND+j]=1.0;}
			if(s2h[i*ND+j]<=0.0){s2h[i*ND+j]=0.0;}
		}
	}

	//if(keypress()){return 0;}//�L�[�҂����

	time1=time1+1.0;								//�v�Z�J�E���g���̉��Z
	if(time1<time1max){goto start;}	//�ő�J�E���g���ɓ��B�������ǂ����̔��f

end:;
  return 0;
}

//************ ������̐ݒ�T�u���|�`�� *************
void ini000(double *s1h, double *s2h, int ND)
{
	int i, j;
	//srand(time(NULL)); // ����������
	int nd=ND, ndm=ND-1, nd2=ND/2;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//s1h[i][j]=0.2*DRND(1.0); s2h[i][j]=0.2*DRND(1.0);//����ő�20%�̗����ɂĐݒ�
			//if(abs(j-nd2)<(nd/40)){s1h[i][j]=DRND(1.0); s2h[i][j]=DRND(1.0);}
			s1h[i*ND+j]=0.2*DRND(1.0);//����ő�20%�̗����ɂĐݒ�
			s2h[i*ND+j]=0.2*DRND(1.0);//����ő�20%�̗����ɂĐݒ�
		}
	}
}

//******* Sin, Cos �̃e�[�u������уr�b�g���]�e�[�u���̐ݒ� ***************
void table(double *s, double *c, int *ik, int ND, int ig)
{
	int it, it1, it2, mc, mn;
	double q;
	int nd=ND, ndm=ND-1, nd2=ND/2;

	q=2.0*PI/nd;
	for(it=0;it<=nd2-1;it++){ c[it]=cos(q*it); s[it]=sin(q*it); }//Sin, Cos �̃e�[�u��

	ik[0]=0; mn=nd2; mc=1;
	for(it1=1;it1<=ig;it1++){
		for(it2=0;it2<=mc-1;it2++){
			ik[it2+mc]=ik[it2]+mn;				//�r�b�g���]�e�[�u��
		}
		mn=mn/2; mc=2*mc;
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
				tr=xrf[ka]-xrf[kb];  					tj=xif[ka]-xif[kb];
				xrf[ka]=xrf[ka]+xrf[kb]; 			xif[ka]=xif[ka]+xif[kb];
				xrf[kb]=tr*c[ix]-tj*qs*s[ix];	xif[kb]=tj*c[ix]+tr*qs*s[ix];
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
			//xrf[ic]=xr[ir][ic];	xif[ic]=xi[ir][ic];
			xrf[ic]=xr[ir*ND+ic];
			xif[ic]=xi[ir*ND+ic];
		}
		fft(xrf, xif, s, c, ND, ig);
		for(ic=0;ic<=ndm;ic++){
			//xr[ir][ic]=xrf[ik[ic]];	xi[ir][ic]=xif[ik[ic]];
			xr[ir*ND+ic]=xrf[ik[ic]];
			xi[ir*ND+ic]=xif[ik[ic]];
		}
	}
	for(ic=0;ic<=ndm;ic++){
		for(ir=0;ir<=ndm;ir++){
			//xrf[ir]=xr[ir][ic];	xif[ir]=xi[ir][ic];
			xrf[ir]=xr[ir*ND+ic];
			xif[ir]=xi[ir*ND+ic];
		}
		fft(xrf, xif, s, c, ND, ig);
		for(ir=0;ir<=ndm;ir++){
			//xr[ir][ic]=xrf[ik[ir]];	xi[ir][ic]=xif[ik[ir]];
			xr[ir*ND+ic]=xrf[ik[ir]];
			xi[ir*ND+ic]=xif[ik[ir]];
		}
	}
	if(qs>0.0){return;}
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//xr[i][j]=xr[i][j]/nd/nd;	xi[i][j]=xi[i][j]/nd/nd;
			xr[i*ND+j]=xr[i*ND+j]/nd/nd;
			xi[i*ND+j]=xi[i*ND+j]/nd/nd;
		}
	}

}

//************ �f�[�^�ۑ��T�u���[�`�� *******************************
void datsave(double *s1h, double *s2h, int ND)
{
	FILE *stream;		//�X�g���[���̃|�C���^�ݒ�
	int i, j;			//����
	int nd=ND, ndm=ND-1, nd2=ND/2;

	stream = fopen("test.dat", "a");	//�������ݐ�̃t�@�C����ǋL�����ŃI�[�v��
	fprintf(stream, "%e\n", time1);		//�v�Z�J�E���g���̕ۑ�
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//fprintf(stream, "%e  %e ", s1h[i][j], s2h[i][j]);//��̃f�[�^�ۑ�
			fprintf(stream, "%e  %e ", s1h[i*ND+j], s2h[i*ND+j]);
		}
	}
	fprintf(stream, "\n");	//���s�̏�������
	fclose(stream);			//�t�@�C�����N���[�Y
}

void datsave_paraview(double *s1h, double *s2h, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j;
	int nd=ND, ndm=ND-1, nd2=ND/2;
	
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
	fprintf(fp,"SCALARS Phase_field_1 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e ", s1h[i][j], s2h[i][j]);//��̃f�[�^�ۑ�
			//fprintf(fp,"%10.6f\n", s1h[i][j]);
			fprintf(fp,"%10.6f\n", s1h[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS Phase_field_2 float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e ", s1h[i][j], s2h[i][j]);//��̃f�[�^�ۑ�
			//fprintf(fp,"%10.6f\n", s2h[i][j]);
			fprintf(fp,"%10.6f\n", s2h[i*ND+j]);
		}
	}
	fclose(fp);
}