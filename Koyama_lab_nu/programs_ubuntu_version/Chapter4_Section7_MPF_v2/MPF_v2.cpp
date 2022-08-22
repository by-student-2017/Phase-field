//�f�[�^��ǂݍ���No.1-5
//���ԍ���1-5
//�E�ʃG�l���M�[���x�ɕ��ʍ��ˑ������l���i�ł��邪���Ă��Ȃ��j
//�Z�x����l���i�Q�����j�A�����̌������͂P�� No.5

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//�����̐ݒ�

//#define ND 100				//�����v�Z�ɂ�����v�Z�̈��ӂ̕�����
//#define N 6					//�l�����錋�����ʂ̐��{�P

	//int nd=ND, ndm=ND-1;	//�v�Z�̈�̈�ӂ̍���������(�����u���b�N��), ND-1���`
	//int nm=N-1, nmm=N-2;	//�l�����錋�����ʂ̐��AN-2�i�l�����錋�����ʂ̐��|�P�j���`
	double PI=3.141592;		//pi
	double RR=8.3145;		//�K�X�萔
	double time1;			//�v�Z�J�E���g��
	//double ph[N][ND][ND], ph2[N][ND][ND];	//�t�F�[�Y�t�B�[���h�A�t�F�[�Y�t�B�[���h�⏕�z��
	double c0;				//���ϑg��
	//double ch[ND][ND], ch2[ND][ND];	//�Z�x��A�Z�x��̕⏕�z��
	//double aij[N][N];		//���z�G�l���M�[�W��
	//double wij[N][N];		//�y�i���e�B�[���̌W��
	//double tij[N][N];		//���E�̈Փ��x
	//double eij[N][N];		//���E�ړ��̋쓮��
	//int m00h[N][ND][ND];	//�ʒu(i,j)����т��̎���(i�}1,j�}1)�ɂ����āAp���O�ł͂Ȃ����ʂ̔ԍ�
	//int n00h[ND][ND];		//�ʒu(i,j)����т��̎���(i�}1,j�}1)�ɂ����āAp���O�ł͂Ȃ����ʂ̌�
	int Nstep, iout;

	void ini000(double *ph, double *ch, int ND, int N);	//������̐ݒ�T�u���[�`��
	void datsave(double *ph, double *ch, int ND, int N);//�f�[�^�ۑ��T�u���[�`��
	void datsave_paraview(double *ph, double *ch, int ND, int N);	//�f�[�^�ۑ��T�u���[�`��
	void datin(double *ph, int ND, int N);	//�f�[�^���̓T�u���[�`��

//******* ���C���v���O���� ******************************************
int main(void)
{
	int ND, N;
	int nd, ndm, nm, nmm;
	
	int i, j, k, l, ii, jj, kk, ll, it;	//����
	int ip, im, jp, jm;									//����
	int n1, n2, n3;											//����
	int nalph;
	int n00;//�ʒu(i,j)����т��̎���(i�}1,j�}1)�ɂ����āAp���O�ł͂Ȃ����ʂ̌�
	//int n000;//�ʒu(i,j)�ɂ����āAp���O�ł͂Ȃ����ʂ̌��in00>=n000�j
	//double cah[ND][ND], cbh[ND][ND];//�Ǐ����t�i���s�ڐ����j�̔Z�x��
	//double sah[ND][ND], sbh[ND][ND];//�����ƃ����̑��݊m���iS����S���j
	double time1max;//�v�Z�J�E���g���̍ő�l�i�v�Z�I���J�E���g�j
	double delt, L, b1;//���Ԃ����݁A�v�Z�̈�P�ӂ̒����A�����u���b�N��ӂ̒���
	double M1, M0;		//���E�̈Փ��x
	double W1, W0;		//�y�i���e�B�[���̌W��
	double K1, K0;		//���z�G�l���M�[�W��
	double E1, E0;		//���E�ړ��̋쓮��
	double temp;		//���x
	double sum1, sum2, sum3;//�e��̘a�̍�ƕϐ�
	double pddtt;		//�t�F�[�Y�t�B�[���h�̎��ԕω���

	double gamma, gamma0;//���E�G�l���M���x
	double delta;	//���E���i�����u���b�N���ɂĕ\���j
	double amobi;	//���E�̈Փ��x
	double vm0;		//�����̐�

	double A_alph, A_beta, c_alph0, c_beta0;//���w�I���R�G�l���M�[���̃p�����[�^
	double A_alph0, A_beta0;		//���w�I���R�G�l���M�[���̃p�����[�^
	double c, sa, sb, Da, Db, D;	//�Z�x��AS���AS���A�g�U�W���iD���AD���j�A�d�ݕt�g�U�W��D
	double gii, gjj, cii, cjj;		//���w�I���R�G�l���M�[�A�Ǐ����t�Z�x
	double dgdc;					//���w�I���R�G�l���M�[�̔Z�x����
	double cddtt;					//�Z�x�̕ω���
	double dc0;						//�Z�x��␳�p�̍�ƕϐ�
	double dev1_a, dev1_b, dev2_a, dev2_b;//�g�U���������ɂ���������v�Z�̍ۂ̍�ƕϐ�


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
	N       = int(data[1]);
	delt    = data[2];
	temp    = data[3];
	L       = data[4];
	vm0     = data[5];
	gamma0  = data[6];
	delta   = data[7];
	K0      = data[8];
	W0      = data[9];
	amobi   = data[10];
	M0      = data[11];
	E0      = data[12];
	time1max= int(data[13]);
	Nstep   = int(data[14]);
	nalph   = data[15];
	A_alph0 = data[16];
	c_alph0 = data[17];
	A_beta0 = data[18];
	c_beta0 = data[19];
	Da      = data[20];
	Db      = data[21];
	D       = data[22];
	printf("---------------------------------\n");
	//
	nd=ND;
	ndm=ND-1;
	nm=N-1;
	nmm=N-2;
	//
	double *ph   = (double *)malloc(sizeof(double)*( N*ND*ND + ND*ND + ND ));	//�t�F�[�Y�t�B�[���h
	double *ph2  = (double *)malloc(sizeof(double)*( N*ND*ND + ND*ND + ND ));	//�t�F�[�Y�t�B�[���h�⏕�z��
	double *ch   = (double *)malloc(sizeof(double)*( ND*ND + ND ));
	double *ch2  = (double *)malloc(sizeof(double)*( ND*ND + ND ));
	double *aij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//���z�G�l���M�[�W��
	double *wij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//�y�i���e�B�[���̌W��
	double *tij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//���E�̈Փ��x
	double *eij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//���E�ړ��̋쓮��
	double *m00h = (double *)malloc(sizeof(double)*( N*ND*ND + ND*ND + ND ));	//�ʒu(i,j)����т��̎���(i�}1,j�}1)�ɂ����āAp��0�ł͂Ȃ����ʂ̔ԍ�
	double *n00h = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//�ʒu(i,j)����т��̎���(i�}1,j�}1)�ɂ����āAp��0�ł͂Ȃ����ʂ̌�
	//
	double *cah  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//�Ǐ����t�i���s�ڐ����j�̔Z�x��
	double *cbh  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//�Ǐ����t�i���s�ڐ����j�̔Z�x��
	double *sah  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//�����ƃ����̑��݊m���iS����S���j
	double *sbh  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//�����ƃ����̑��݊m���iS����S���j
	//
	//printf("delt(0.2)=  "); scanf(" %lf",&delt);//���ԍ��݂̓���	//	delt=0.2;
	//printf("c0(0.4)=  "); scanf(" %lf",&c0);//���ϑg���̓���	//c0=0.2;

	//temp=1000.0; 							//���x(K)
	//L=2000.0;									//�v�Z�̈�̈�ӂ̒���(nm)
	b1=L/(double)ND*1.0e-9; 	//�����v���b�N�P�ӂ̒���(m)

	//vm0=7.0e-6;								//�����̐�
	//gamma=0.5*vm0/RR/temp/b1;	//���E�G�l���M���x�i0.5J/m^2�j�𖳎�����
	gamma=gamma0*vm0/RR/temp/b1;	//���E�G�l���M���x�i0.5J/m^2�j�𖳎�����
	//delta=5.0;								//���E���i�����u���b�N���ɂĕ\���j

	//K1=8.0*delta*gamma/PI/PI;	//���z�G�l���M�[�W��[��(3.23)]
	K1=K0*delta*gamma/PI/PI;	//���z�G�l���M�[�W��[��(3.23)]
	//W1=4.0*gamma/delta;				//�y�i���e�B�[���̌W��[��(3.23)]
	W1=W0*gamma/delta;				//�y�i���e�B�[���̌W��[��(3.23)]
	//amobi=1.0;
	//M1=amobi*PI*PI/(8.0*delta);	//���E�̈Փ��x[��(3.23)]
	M1=amobi*PI*PI/(M0*delta);	//���E�̈Փ��x[��(3.23)]
	//E1=50.0/RR/temp;					//���E�ړ��̋쓮��
	E1=E0/RR/temp;					//���E�ړ��̋쓮��

	time1=0.0;
	//time1max=10000001.;//�v�Z�J�E���g���̏����l�ƍő�l

	//nalph=4;	//�����̌������ʂ̐��i�l�����Ă��郿���̌������̐��j
	//A_alph=1.0e+02/RR/temp;//���w�I���R�G�l���M�[���̃p�����[�^
	A_alph=A_alph0/RR/temp;//���w�I���R�G�l���M�[���̃p�����[�^
	//c_alph0=0.1;
	//A_beta=1.0e+02/RR/temp;
	A_beta=A_beta0/RR/temp;
	//c_beta0=0.9;
	//Da=Db=D=1.0;	//�e���̊g�U�W���͑S�ē������Ɖ���

//*** ��(4.36)-��(4.39)�̔z��iK,W,M,E�j�̐ݒ� **************************
	for(i=1;i<=nm;i++){
		for(j=i+1;j<=nm;j++){
			//wij[i][j]=wij[j][i]=W1;
			//aij[i][j]=aij[j][i]=K1;
			//tij[i][j]=tij[j][i]=M1;
			//eij[i][j]=eij[j][i]=0.0;
			wij[i*ND+j]=wij[j*ND+i]=W1;
			aij[i*ND+j]=aij[j*ND+i]=K1;
			tij[i*ND+j]=tij[j*ND+i]=M1;
			eij[i*ND+j]=eij[j*ND+i]=0.0;
		}
		//wij[i][i]=0.0;  aij[i][i]=0.0;  tij[i][i]=0.0;  eij[i][i]=0.0;
		wij[i*ND+i]=0.0; aij[i*ND+i]=0.0; tij[i*ND+i]=0.0; eij[i*ND+i]=0.0;
	}

//*** ������̐ݒ�ƕ`��Window�\�� *****************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=1;k<=nm;k++){
				//ph[k][i][j]=0.0;  ph2[k][i][j]=0.0;//�t�F�[�Y�t�B�[���h�A����ѕ⏕�z��̏�����
				 ph[k*ND*ND+i*ND+j]=0.0;	//�t�F�[�Y�t�B�[���h�̏�����
				ph2[k*ND*ND+i*ND+j]=0.0;	//�⏕�z��̏�����
			}
			//ch[i][j]=0.0;  ch2[i][j]=0.0;//�Z�x��A����ѕ⏕�z��̏�����
			 ch[i*ND+j]=0.0;	//�Z�x��̏�����
			ch2[i*ND+j]=0.0;	//�⏕�z��̏�����
		}
	}

	ini000(ph, ch, ND, N);//������̐ݒ�
	//datin(ph, ND, N);//������̓���

//**** �V�~�����[�V�����X�^�[�g ******************************
//Nstep = 200;
iout = -1;
start: ;

	if((((int)(time1) % Nstep)==0)) {datsave(ph, ch, ND, N);}//���J�Ԃ��J�E���g���ɏ��ۑ�
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(ph, ch, ND, N);}//���J�Ԃ��J�E���g���ɏ��ۑ�
	//if(time1==200.) {datsave();}//����̎��Ԃ̏��ۑ�

//******  �Ǐ����t�g���ic����c���j�̌v�Z[��(4.47)]  ********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//c=ch[i][j];	//�Z�x��
			c=ch[i*ND+j];	//�Z�x��
			//sum1=0.0; for(ii=1;ii<=nalph;ii++){ sum1+=ph[ii][i][j]; }
			sum1=0.0; for(ii=1;ii<=nalph;ii++){ sum1+=ph[ii*ND*ND+i*ND+j]; }
			//sa=sah[i][j]=sum1;	//S���̌v�Z[��(4.46)]
			sa=sah[i*ND+j]=sum1;	//S���̌v�Z[��(4.46)]
			//sum1=0.0; for(ii=nalph+1;ii<=nm;ii++){ sum1+=ph[ii][i][j]; }
			sum1=0.0; for(ii=nalph+1;ii<=nm;ii++){ sum1+=ph[ii*ND*ND+i*ND+j]; }
			//sb=sbh[i][j]=sum1;	//S���̌v�Z[��(4.46)]
			sb=sbh[i*ND+j]=sum1;	//S���̌v�Z[��(4.46)]

			//�Ǐ����t�g���̌v�Z[��(4.47)]
			//cah[i][j]=(A_beta*c+(A_alph*c_alph0-A_beta*c_beta0)*sb)/(A_beta*sa+A_alph*sb);
			cah[i*ND+j]=(A_beta*c+(A_alph*c_alph0-A_beta*c_beta0)*sb)/(A_beta*sa+A_alph*sb);
			//cbh[i][j]=(A_alph*c+(A_beta*c_beta0-A_alph*c_alph0)*sa)/(A_beta*sa+A_alph*sb);
			cbh[i*ND+j]=(A_alph*c+(A_beta*c_beta0-A_alph*c_alph0)*sa)/(A_beta*sa+A_alph*sb);
			//if(cah[i][j]>=1.0){cah[i][j]=1.0;}  if(cah[i][j]<=0.0){cah[i][j]=0.0;}//�Z�x��̕ψ�␳
			if(cah[i*ND+j]>=1.0){cah[i*ND+j]=1.0;}  if(cah[i*ND+j]<=0.0){cah[i*ND+j]=0.0;}//�Z�x��̕ψ�␳
			//if(cbh[i][j]>=1.0){cbh[i][j]=1.0;}  if(cbh[i][j]<=0.0){cbh[i][j]=0.0;}
			if(cbh[i*ND+j]>=1.0){cbh[i*ND+j]=1.0;}  if(cbh[i*ND+j]<=0.0){cbh[i*ND+j]=0.0;}
		}
	}

//**** �e�����v���b�N�ɂ�����n00h[i][j]��m00h[n00][i][j]�𒲍� *********************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1;  jp=j+1;  jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//�����I���E����
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}

//--- �ʒu(i,j)����т��̎���(i�}1,j�}1)�ɂ����āAp���O�ł͂Ȃ����ʂ̌�---
			n00=0;
			for(ii=1;ii<=nm;ii++){
				//if( (ph[ii][i][j]>0.0)||
				//  ( (ph[ii][i][j]==0.0)&&(ph[ii][ip][j]>0.0)||
				//  	(ph[ii][im][j]>0.0)||
				//  	(ph[ii][i][jp]>0.0)||
				//  	(ph[ii][i][jm]>0.0)   ) ){
				//  		n00++; m00h[n00][i][j]=ii;
				//		//printf("%d  ", n00);
				//}
				if( (ph[ii*ND*ND+i*ND+j]>0.0)||
				  ( (ph[ii*ND*ND+i*ND+j]==0.0)&&(ph[ii*ND*ND+ip*ND+j]>0.0)||
				  	(ph[ii*ND*ND+im*ND+j]>0.0)||
				  	(ph[ii*ND*ND+i*ND+jp]>0.0)||
				  	(ph[ii*ND*ND+i*ND+jm]>0.0)   ) ){
				  		n00++; m00h[n00*ND*ND+i*ND+j]=ii;
						//printf("%d  ", n00);
				}
			}
			//n00h[i][j]=n00;
			n00h[i*ND+j]=n00;
//--------------------------------------------------------------------------
		}
	}


//***** �t�F�[�Y�t�B�[���h�̔��W�������̌v�Z **********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1; im=i-1;  jp=j+1;  jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//�����I���E����
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}

			//for(n1=1; n1<=n00h[i][j]; n1++){
			for(n1=1; n1<=n00h[i*ND+j]; n1++){
				//ii=m00h[n1][i][j];  pddtt=0.0;
				ii=m00h[n1*ND*ND+i*ND+j];  pddtt=0.0;
				//for(n2=1; n2<=n00h[i][j]; n2++){
				for(n2=1; n2<=n00h[i*ND+j]; n2++){
					//jj=m00h[n2][i][j];  sum1=0.0;
					jj=m00h[n2*ND*ND+i*ND+j];  sum1=0.0;
					//for(n3=1; n3<=n00h[i][j]; n3++){
					for(n3=1; n3<=n00h[i*ND+j]; n3++){
						//kk=m00h[n3][i][j];
						kk=m00h[n3*ND*ND+i*ND+j];
						//sum1+=0.5*(aij[ii][kk]-aij[jj][kk])*(ph[kk][ip][j]+ph[kk][im][j]
						//									+ph[kk][i][jp]+ph[kk][i][jm]-4.0*ph[kk][i][j])
						//		 +(wij[ii][kk]-wij[jj][kk])*ph[kk][i][j];//[��(4.41)�̈ꕔ]
						sum1+=0.5*(aij[ii*ND+kk]-aij[jj*ND+kk])*(ph[kk*ND*ND+ip*ND+j]+ph[kk*ND*ND+im*ND+j]
																+ph[kk*ND*ND+i*ND+jp]+ph[kk*ND*ND+i*ND+jm]-4.0*ph[kk*ND*ND+i*ND+j])
								 +(wij[ii*ND+kk]-wij[jj*ND+kk])*(ph[kk*ND*ND+i*ND+j]);//[��(4.41)�̈ꕔ]
					}

					//if(ii<=nalph){gii=A_alph*(cah[i][j]-c_alph0)*(cah[i][j]-c_alph0); cii=cah[i][j]; }
					if(ii<=nalph){gii=A_alph*(cah[i*ND+j]-c_alph0)*(cah[i*ND+j]-c_alph0); cii=cah[i*ND+j]; }
					//else{gii=A_beta*(cbh[i][j]-c_beta0)*(cbh[i][j]-c_beta0); cii=cbh[i][j]; }//��(4.43)
					else{gii=A_beta*(cbh[i*ND+j]-c_beta0)*(cbh[i*ND+j]-c_beta0); cii=cbh[i*ND+j]; }//��(4.43)
					//if(jj<=nalph){gjj=A_alph*(cah[i][j]-c_alph0)*(cah[i][j]-c_alph0); cjj=cah[i][j]; }
					if(jj<=nalph){gjj=A_alph*(cah[i*ND+j]-c_alph0)*(cah[i*ND+j]-c_alph0); cjj=cah[i*ND+j]; }
					//else{gjj=A_beta*(cbh[i][j]-c_beta0)*(cbh[i][j]-c_beta0); cjj=cbh[i][j]; }//��(4.43)
					else{gjj=A_beta*(cbh[i*ND+j]-c_beta0)*(cbh[i*ND+j]-c_beta0); cjj=cbh[i*ND+j]; }//��(4.43)
						//dgdc=2.0*A_alph*(cah[i][j]-c_alph0);//���w�I���R�G�l���M�[�̔Z�x����[��(4.48)]�̌v�Z
						dgdc=2.0*A_alph*(cah[i*ND+j]-c_alph0);//���w�I���R�G�l���M�[�̔Z�x����[��(4.48)]�̌v�Z
						//dgdc=2.0*A_beta*(cbh[i][j]-c_beta0);
						
						//pddtt+=-2.0*tij[ii][jj]/double(n00h[i][j])*(sum1+gii-gjj-(cii-cjj)*dgdc );
						pddtt+=-2.0*tij[ii*ND+jj]/double(n00h[i*ND+j])*(sum1+gii-gjj-(cii-cjj)*dgdc );
						//�t�F�[�Y�t�B�[���h�̔��W������[��(4.41)]
				}
				//ph2[ii][i][j]=ph[ii][i][j]+pddtt*delt;//�t�F�[�Y�t�B�[���h�̎��Ԕ��W�i�z��@�j
				//if(ph2[ii][i][j]>=1.0){ph2[ii][i][j]=1.0;}//�t�F�[�Y�t�B�[���h�̕ψ�␳
				//if(ph2[ii][i][j]<=0.0){ph2[ii][i][j]=0.0;}
				ph2[ii*ND*ND+i*ND+j]=ph[ii*ND*ND+i*ND+j]+pddtt*delt;//�t�F�[�Y�t�B�[���h�̎��Ԕ��W�i�z��@�j
				if(ph2[ii*ND*ND+i*ND+j]>=1.0){ph2[ii*ND*ND+i*ND+j]=1.0;}//�t�F�[�Y�t�B�[���h�̕ψ�␳
				if(ph2[ii*ND*ND+i*ND+j]<=0.0){ph2[ii*ND*ND+i*ND+j]=0.0;}
			}
		}//j
	}//i

//*****  �Z�x��̎��Ԕ��W�̌v�Z�i�g�U�������j**********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			ip=i+1;  im=i-1;  jp=j+1;  jm=j-1;
			if(i==ndm){ip=0;}  if(i==0){im=ndm;}//�����I���E����
			if(j==ndm){jp=0;}  if(j==0){jm=ndm;}

			//�g�U���������ɂ���������v�Z
			//dev1_a=0.25*( (sah[ip][j]-sah[im][j])*(cah[ip][j]-cah[im][j])
			//			 +(sah[i][jp]-sah[i][jm])*(cah[i][jp]-cah[i][jm]) );
			dev1_a=0.25*( (sah[ip*ND+j]-sah[im*ND+j])*(cah[ip*ND+j]-cah[im*ND+j])
						 +(sah[i*ND+jp]-sah[i*ND+jm])*(cah[i*ND+jp]-cah[i*ND+jm]) );
			//dev1_b=0.25*( (sbh[ip][j]-sbh[im][j])*(cbh[ip][j]-cbh[im][j])
			//			 +(sbh[i][jp]-sbh[i][jm])*(cbh[i][jp]-cbh[i][jm]) );
			dev1_b=0.25*( (sbh[ip*ND+j]-sbh[im*ND+j])*(cbh[ip*ND+j]-cbh[im*ND+j])
						 +(sbh[i*ND+jp]-sbh[i*ND+jm])*(cbh[i*ND+jp]-cbh[i*ND+jm]) );
			//dev2_a=sah[i][j]*(cah[ip][j]+cah[im][j]+cah[i][jp]+cah[i][jm]-4.0*cah[i][j]);
			dev2_a=sah[i*ND+j]*(cah[ip*ND+j]+cah[im*ND+j]+cah[i*ND+jp]+cah[i*ND+jm]-4.0*cah[i*ND+j]);
			//dev2_b=sbh[i][j]*(cbh[ip][j]+cbh[im][j]+cbh[i][jp]+cbh[i][jm]-4.0*cbh[i][j]);
			dev2_b=sbh[i*ND+j]*(cbh[ip*ND+j]+cbh[im*ND+j]+cbh[i*ND+jp]+cbh[i*ND+jm]-4.0*cbh[i*ND+j]);

			cddtt=Da*(dev1_a+dev2_a)+Db*(dev1_b+dev2_b);	//�g�U������[��(4.42)]
			//ch2[i][j]=ch[i][j]+cddtt*delt;	//�Z�x��̎��Ԕ��W(�z��@)
			//ch2[i][j]=ch[i][j]+cddtt*delt+(2.*DRND(1.)-1.)*0.001;	//�Z�x��̎��Ԕ��W(�z��@)
			ch2[i*ND+j]=ch[i*ND+j]+cddtt*delt+(2.*DRND(1.)-1.)*0.001;	//�Z�x��̎��Ԕ��W(�z��@)
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=1;k<=nm;k++){
				//ph[k][i][j]=ph2[k][i][j];//�⏕�z�����z��Ɉړ��i�t�F�[�Y�t�B�[���h�j
				ph[k*ND*ND+i*ND+j]=ph2[k*ND*ND+i*ND+j];//�⏕�z�����z��Ɉړ��i�t�F�[�Y�t�B�[���h�j
			}
			//ch[i][j]=ch2[i][j];//�⏕�z�����z��Ɉړ��i�Z�x��j
			ch[i*ND+j]=ch2[i*ND+j];//�⏕�z�����z��Ɉړ��i�Z�x��j
		}
	}

//*** �t�F�[�Y�t�B�[���h�̋K�i���␳ ***********************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k][i][j]; }
			//for(k=1;k<=nm;k++){ ph[k][i][j]=ph[k][i][j]/sum1; }
			sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k*ND*ND+i*ND+j]; }
			for(k=1;k<=nm;k++){ ph[k*ND*ND+i*ND+j]=ph[k*ND*ND+i*ND+j]/sum1; }
		}
	}

//*** �Z�x��̎��x�␳ *************************************************************
	//sum1=0.0; for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ sum1+=ch[i][j]; } }
	sum1=0.0; for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ sum1+=ch[i*ND+j]; } }
	dc0=sum1/nd/nd-c0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//ch[i][j]=ch[i][j]-dc0;
			//if(ch[i][j]>1.0){ch[i][j]=1.0;}   if(ch[i][j]<0.0){ch[i][j]=0.0;}
			ch[i*ND+j]=ch[i*ND+j]-dc0;
			if(ch[i*ND+j]>1.0){ch[i*ND+j]=1.0;}
			if(ch[i*ND+j]<0.0){ch[i*ND+j]=0.0;}
	  }
	}

//*********************************************************************
	//if(keypress()){return 0;}//�L�[�҂����
	time1=time1+1.;  if(time1<time1max){goto start;}//�ő�J�E���g���ɓ��B�������ǂ����̔��f

	end:;
  return 0;
}

//************ ������(�t�F�[�Y�t�B�[���h�ƔZ�x��)�̐ݒ�T�u���[�`�� *************
void ini000(double *ph, double *ch, int ND, int N)
{
	int i, j, k, l, it;		//����
	int ii, jj, kk;				//����
	int ip, im, jp, jm;		//����
	int igx, igy, ixmin, ixmax, iymin, iymax;	//�����u���b�N���W�n�̐ݒ�
	double x, y, xmin, xmax, ymin, ymax;			//�K�i�����W�n�̐ݒ�
	double sum1, sum2, t, r0, phi, r;					//��ƕϐ�
	//srand(time(NULL)); //����������
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//ch[i][j]=0.5;//�����̔Z�x��i�ߖO�a�ŗn�́j��ݒ�
			ch[i*ND+j]=0.5;//�����̔Z�x��i�ߖO�a�ŗn�́j��ݒ�
		}
	}

	datin(ph, ND, N);//�����̌����g�D�f�[�^��ǂݍ���

	xmin=-1.0; xmax=1.0; ymin=-1.0;  ymax=1.0;	//�K�i�����W�n
	ixmin=0; ixmax=ND; iymin=0;  iymax=ND;			//�����u���b�N���W�n

	r0=2.0;//�����̃T�C�Y

//�ȉ��A�����̗��E�̂R�d�_�i�U�ӏ��j�Ƀ����̊j��u��
	x=0.5;           igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
	y=0.5/sqrt(3.);  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			r=sqrt( (double(i-igx))*(double(i-igx))+(double(j-igy))*(double(j-igy)) );
			//if(r<=r0){ ph[nm][i][j]=1.0; 	for(k=1;j<=nm-1;k++){ph[k][i][j]=0.0;} ch[i][j]=0.9; }
			if(r<=r0){ ph[nm*ND*ND+i*ND+j]=1.0;
				for(k=1;j<=nm-1;k++){ph[k*ND*ND+i*ND+j]=0.0;} 
				ch[i*ND+j]=0.9;
			}
		}
	}

	x=-0.5;          igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
	y=0.5/sqrt(3.);  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			r=sqrt( (double(i-igx))*(double(i-igx))+(double(j-igy))*(double(j-igy)) );
			//if(r<=r0){ ph[nm][i][j]=1.0; 	for(k=1;j<=nm-1;k++){ph[k][i][j]=0.0;} ch[i][j]=0.9;  }
			if(r<=r0){
				ph[nm*ND*ND+i*ND+j]=1.0;
				for(k=1;j<=nm-1;k++){ph[k*ND*ND+i*ND+j]=0.0;}
				ch[i*ND+j]=0.9;
			}
		}
	}

	x=0.5;           igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
	y=-0.5/sqrt(3.); igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			r=sqrt( (double(i-igx))*(double(i-igx))+(double(j-igy))*(double(j-igy)) );
			//if(r<=r0){ ph[nm][i][j]=1.0; 	for(k=1;j<=nm-1;k++){ph[k][i][j]=0.0;} ch[i][j]=0.9;  }
			if(r<=r0){
				ph[nm*ND*ND+i*ND+j]=1.0;
				for(k=1;j<=nm-1;k++){ph[k*ND*ND+i*ND+j]=0.0;}
				ch[i*ND+j]=0.9;
			}
		}
	}

	x=-0.5;           igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
	y=-0.5/sqrt(3.);  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			r=sqrt( (double(i-igx))*(double(i-igx))+(double(j-igy))*(double(j-igy)) );
			//if(r<=r0){ ph[nm][i][j]=1.0; 	for(k=1;j<=nm-1;k++){ph[k][i][j]=0.0;} ch[i][j]=0.9;  }
			if(r<=r0){
				ph[nm*ND*ND+i*ND+j]=1.0;
				for(k=1;j<=nm-1;k++){ph[k*ND*ND+i*ND+j]=0.0;}
				ch[i*ND+j]=0.9;
			}
		}
	}

	x=0.0;          igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
	y=1./sqrt(3.);  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			r=sqrt( (double(i-igx))*(double(i-igx))+(double(j-igy))*(double(j-igy)) );
			//if(r<=r0){ ph[nm][i][j]=1.0; 	for(k=1;j<=nm-1;k++){ph[k][i][j]=0.0;} ch[i][j]=0.9;  }
			if(r<=r0){
				ph[nm*ND*ND+i*ND+j]=1.0;
				for(k=1;j<=nm-1;k++){ph[k*ND*ND+i*ND+j]=0.0;}
				ch[i*ND+j]=0.9;
			}
		}
	}

	x=0.0;           igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
	y=-1./sqrt(3.);  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			r=sqrt( (double(i-igx))*(double(i-igx))+(double(j-igy))*(double(j-igy)) );
			//if(r<=r0){ ph[nm][i][j]=1.0; 	for(k=1;j<=nm-1;k++){ph[k][i][j]=0.0;} ch[i][j]=0.9;  }
			if(r<=r0){
				ph[nm*ND*ND+i*ND+j]=1.0;
				for(k=1;j<=nm-1;k++){ph[k*ND*ND+i*ND+j]=0.0;}
				ch[i*ND+j]=0.9;
			}
		}
	}

//--- �t�F�[�Y�t�B�[���h�̋K�i���␳�ƕ��ϑg���̎Z�o ------------------------
	sum2=0.0;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k][i][j]; }
			//for(k=1;k<=nm;k++){ ph[k][i][j]=ph[k][i][j]/sum1; }
			//sum2+=ch[i][j];
			sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k*ND*ND+i*ND+j]; }
			for(k=1;k<=nm;k++){ ph[k*ND*ND+i*ND+j]=ph[k*ND*ND+i*ND+j]/sum1; }
			sum2+=ch[i*ND+j];
		}
	}
	c0=sum2/nd/nd;

}

//************ �f�[�^�ۑ��T�u���[�`�� *******************************
void datsave(double *ph, double *ch, int ND, int N)
{
	FILE		*stream;//�X�g���[���̃|�C���^�ݒ�
	int 		i, j, k;//����
	double 	col;
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;

	stream = fopen("test.dat", "a");//�������ސ�̃t�@�C����ǋL�����ŃI�[�v��
	fprintf(stream, "%e  \n", time1);//�v�Z�J�E���g���̕ۑ�
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=nm;k++){
				//fprintf(stream, "%e   ", ph[k][i][j]);//�t�F�[�Y�t�B�[���h�̕ۑ�
				fprintf(stream, "%e   ", ph[k*ND*ND+i*ND+j]);//�t�F�[�Y�t�B�[���h�̕ۑ�
			}
			//fprintf(stream, "%e   ", ch[i][j]);//�Z�x��̕ۑ�
			fprintf(stream, "%e   ", ch[i*ND+j]);//�Z�x��̕ۑ�
		}
	}
	fprintf(stream, "\n");//���s�̏�������
	fclose(stream);//�t�@�C�����N���[�Y
}

void datsave_paraview(double *ph, double *ch, int ND, int N)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j, k;
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;
	double *pht  = (double *)malloc(sizeof(double)*( ND*ND + ND ));
	
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			pht[i*ND+j]=0.0;
		}
	}
	
	iout = iout + 1;
	printf("paraview output no.%06d \n",iout);
	for(k=0;k<=nm;k++){
		sprintf(fName,"mpf_N%03d_result%06d.vtk",k,iout);
		fp = fopen(fName, "w");
		fprintf(fp,"# vtk DataFile Version 3.0 \n");
		fprintf(fp,"output.vtk \n");
		fprintf(fp,"ASCII \n");
		fprintf(fp,"DATASET STRUCTURED_POINTS \n");
		fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),1);
		fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
		fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
		fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*1));
		fprintf(fp,"SCALARS phase_field float \n");
		fprintf(fp,"LOOKUP_TABLE default \n");
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(stream, "%e   ", ph[k][i][j]);//�t�F�[�Y�t�B�[���h�̕ۑ�
				//fprintf(fp,"%10.6f\n", ph[k][i][j]);
				fprintf(fp,"%10.6f\n", ph[k*ND*ND+i*ND+j]);
				pht[i*ND+j]+=ph[k*ND*ND+i*ND+j]*float(k+1.0);
			}
		}
		fclose(fp);
	}
	
	sprintf(fName,"mpf_N%03d_result%06d.vtk",(nm+1),iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),1);
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*1));
	fprintf(fp,"SCALARS phase_field float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			fprintf(fp,"%10.6f\n", pht[i*ND+j]);
		}
	}
	fprintf(fp,"SCALARS concentration float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e   ", ch[i][j]);//�Z�x��̕ۑ�
			//fprintf(fp,"%10.6f\n", ch[i][j]);
			fprintf(fp,"%10.6f\n", ch[i*ND+j]);
		}
	}
	fclose(fp);
	free(pht);
}
//*********** �f�[�^���̓T�u���[�`�� **************************
void datin(double *ph, int ND, int N)
{
	FILE		*datin0;//�X�g���[���̃|�C���^�ݒ�
	int 		i, j, k;//����
	double dami1;//��ƕϐ�
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;

	datin0 = fopen("test.dat", "r");//�ǂݍ��݌��̃t�@�C�����I�[�v��
	fscanf(datin0, "%lf  ", &dami1);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=1;k<=nm;k++){
				//fscanf(datin0, "%lf  ", &ph[k][i][j]);//�t�F�[�Y�t�B�[���h�̓ǂݍ���
				fscanf(datin0, "%lf  ", &ph[k*ND*ND+i*ND+j]);//�t�F�[�Y�t�B�[���h�̓ǂݍ���
			}
		}
	}
	fclose(datin0);//�t�@�C�����N���[�Y

}

