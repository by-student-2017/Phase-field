//���������g�D�`���̃v���O����
//���ԍ���1-nm
//���ԍ�nm���t��

#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//�����̐ݒ�

//#define ND 100					//�����v�Z�ɂ�����v�Z�̈��ӂ̕�����
//#define N 5							//�l�����錋�����ʂ̐��{�P

	//int nd=ND, ndm=ND-1;	//�v�Z�̈�̈�ӂ̍���������(�����u���b�N��), ND-1���`
	//int nm=N-1, nmm=N-2;	//�l�����錋�����ʂ̐��AN-2�i�l�����錋�����ʂ̐��|�P�j���`
	double PI=3.141592;		//pi
	double RR=8.3145;		//�K�X�萔
	double time1;			//�v�Z�J�E���g��
	//double ph[N][ND][ND], ph2[N][ND][ND];	//�t�F�[�Y�t�B�[���h�A�t�F�[�Y�t�B�[���h�⏕�z��
	//double aij[N][N];//���z�G�l���M�[�W��
	//double wij[N][N];//�y�i���e�B�[���̌W��
	//double tij[N][N];//���E�̈Փ��x
	//double eij[N][N];//���E�ړ��̋쓮��
	//int m00h[N][ND][ND];//�ʒu(i,j)����т��̎���(i�}1,j�}1)�ɂ����āAp���O�ł͂Ȃ����ʂ̔ԍ�
	//int n00h[ND][ND];//�ʒu(i,j)����т��̎���(i�}1,j�}1)�ɂ����āAp���O�ł͂Ȃ����ʂ̌�
	int Nstep, iout;

	void ini000(double *ph, int ND, int N);	//������̐ݒ�T�u���[�`��
	void datsave(double *ph, int ND, int N);	//�f�[�^�ۑ��T�u���[�`��
	void datsave_paraview(double *ph, int ND, int N);	//�f�[�^�ۑ��T�u���[�`��

//******* ���C���v���O���� ******************************************
int main(void)
{
	int ND, N;
	int nd, ndm, nm, nmm;

	int i, j, k, l, ii, jj, kk, ll, it;	//����
	int ip, im, jp, jm;									//����
	int n1, n2, n3;											//����
	int n00;		//�ʒu(i,j)����т��̎���(i�}1,j�}1)�ɂ����āAp���O�ł͂Ȃ����ʂ̌�
	//int n000;		//�ʒu(i,j)�ɂ����āAp���O�ł͂Ȃ����ʂ̌��in00>=n000�j
	double time1max;		//�v�Z�J�E���g���̍ő�l�i�v�Z�I���J�E���g�j
	double delt, L, b1;		//���Ԃ����݁A�v�Z�̈�P�ӂ̒����A�����u���b�N��ӂ̒���
	double M1, M0;			//���E�̈Փ��x
	double W1, W0;			//�y�i���e�B�[���̌W��
	double K1, K0;			//���z�G�l���M�[�W��
	double E1, E0;			//���E�ړ��̋쓮��
	double temp;			//���x
	double sum1, sum2, sum3;//�e��̘a�̍�ƕϐ�
	double pddtt;			//�t�F�[�Y�t�B�[���h�̎��ԕω���

	double gamma, gamma0;	//���E�G�l���M���x
	double delta;			//���E���i�����u���b�N���ɂĕ\���j
	double amobi;			//���E�̈Փ��x
	double vm0;				//�����̐�

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
	printf("---------------------------------\n");
	//
	nd=ND;
	ndm=ND-1;
	nm=N-1;
	nmm=N-2;
	//
	double *ph   = (double *)malloc(sizeof(double)*( N*ND*ND + ND*ND + ND ));	//�t�F�[�Y�t�B�[���h
	double *ph2  = (double *)malloc(sizeof(double)*( N*ND*ND + ND*ND + ND ));	//�t�F�[�Y�t�B�[���h�⏕�z��
	double *aij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//���z�G�l���M�[�W��
	double *wij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//�y�i���e�B�[���̌W��
	double *tij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//���E�̈Փ��x
	double *eij  = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//���E�ړ��̋쓮��
	double *m00h = (double *)malloc(sizeof(double)*( N*ND*ND + ND*ND + ND ));	//�ʒu(i,j)����т��̎���(i�}1,j�}1)�ɂ����āAp��0�ł͂Ȃ����ʂ̔ԍ�
	double *n00h = (double *)malloc(sizeof(double)*( ND*ND + ND ));	//�ʒu(i,j)����т��̎���(i�}1,j�}1)�ɂ����āAp��0�ł͂Ȃ����ʂ̌�
	//
	//printf("delt(5.0)=  "); scanf(" %lf",&delt);//���ԍ��݂̓���	//	delt=5.0;

	//temp=1000.0; 					//���x(K)
	//L=2000.0;					//�v�Z�̈�̈�ӂ̒���(nm)
	b1=L/(double)ND*1.0e-9; 	//�����v���b�N�P�ӂ̒���(m)

	//vm0=7.0e-6;					/�����̐�
	//gamma=0.5*vm0/RR/temp/b1;	//���E�G�l���M���x�i0.5J/m^2�j�𖳎�����
	gamma=gamma0*vm0/RR/temp/b1;	//���E�G�l���M���x�i0.5J/m^2�j�𖳎�����
	//delta=7.0;					//���E���i�����u���b�N���ɂĕ\���j

	//K1=8.0*delta*gamma/PI/PI;	//���z�G�l���M�[�W��[��(4.40)]
	K1=K0*delta*gamma/PI/PI;	//���z�G�l���M�[�W��[��(4.40)]
	//W1=4.0*gamma/delta;				//�y�i���e�B�[���̌W��[��(4.40)]
	W1=W0*gamma/delta;				//�y�i���e�B�[���̌W��[��(4.40)]
	//amobi=1.0;
	//M1=amobi*PI*PI/(8.0*delta);	//���E�̈Փ��x[��(4.40)]
	M1=amobi*PI*PI/(M0*delta);	//���E�̈Փ��x[��(4.40)]
	//E1=50.0/RR/temp;				//���E�ړ��̋쓮��
	E1=E0/RR/temp;				//���E�ړ��̋쓮��

	time1=0.;					//�v�Z�J�E���g���̏����l
	time1max=10000001.0;		//�v�Z�J�E���g���̍ő�l

//*** ��(4.32) - ��(4.35)�̔z��iK,W,M,E�j�̐ݒ� *****************************
	for(i=1;i<=nm;i++){
		for(j=1;j<=nm;j++){
			//wij[i][j]=W1;  aij[i][j]=K1;  tij[i][j]=0.0;  eij[i][j]=0.0;
			//if( (i==nm)||(j==nm) ){eij[i][j]=E1; tij[i][j]=M1;}
			//if(i>j){eij[i][j]=-eij[i][j];}
			//if(i==j){	wij[i][j]=0.0;  aij[i][j]=0.0;  tij[i][j]=0.0; eij[i][j]=0.0;}
			wij[i*ND+j]=W1;  aij[i*ND+j]=K1;  tij[i*ND+j]=0.0;  eij[i*ND+j]=0.0;
			if( (i==nm)||(j==nm) ){eij[i*ND+j]=E1; tij[i*ND+j]=M1;}
			if(i>j){eij[i*ND+j]=-eij[i*ND+j];}
			if(i==j){	wij[i*ND+j]=0.0;  aij[i*ND+j]=0.0;  tij[i*ND+j]=0.0; eij[i*ND+j]=0.0;}
		}
	}

//*** ������̐ݒ�ƕ`��Window�\�� *****************************************
	for(k=1;k<=nm;k++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				//ph[k][i][j]=0.0;  ph2[k][i][j]=0.0;//�t�F�[�Y�t�B�[���h�A����ѕ⏕�z��̏�����
				 ph[k*ND*ND+i*ND+j]=0.0;	//�t�F�[�Y�t�B�[���h�̏�����
				ph2[k*ND*ND+i*ND+j]=0.0;	//�⏕�z��̏�����
			}
		}
	}
	ini000(ph, ND, N);	//������̐ݒ�

//**** �V�~�����[�V�����X�^�[�g ******************************
//Nstep = 20;
iout = -1;
start: ;

	if((((int)(time1) % Nstep)==0)) {datsave(ph, ND, N);}	//���J�Ԃ��J�E���g���ɏ��ۑ�
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(ph, ND, N);}	//���J�Ԃ��J�E���g���ɏ��ۑ�
	//if(time1==200.) {datsave();}					//����̎��Ԃ̏��ۑ�

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
				//}
				if( (ph[ii*ND*ND+i*ND+j]>0.0)||
					( (ph[ii*ND*ND+i*ND+j]==0.0)&&(ph[ii*ND*ND+ip*ND+j]>0.0)||
					  (ph[ii*ND*ND+im*ND+j]>0.0)||
					  (ph[ii*ND*ND+i*ND+jp]>0.0)||
					  (ph[ii*ND*ND+i*ND+jm]>0.0) ) ){
					  	n00++; m00h[n00*ND*ND+i*ND+j]=ii;
				}
			}
			//n00h[i][j]=n00;
			n00h[i*ND+j]=n00;
//--------------------------------------------------------------------------
		}
	}

//***** ���W�������̌v�Z **********************************
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
						//								  +(wij[ii][kk]-wij[jj][kk])*ph[kk][i][j];	//[��(4.31)�̈ꕔ]
						sum1+=0.5*(aij[ii*ND+kk]-aij[jj*ND+kk])*(ph[kk*ND*ND+ip*ND+j]+ph[kk*ND*ND+im*ND+j]
																+ph[kk*ND*ND+i*ND+jp]+ph[kk*ND*ND+i*ND+jm]-4.0*ph[kk*ND*ND+i*ND+j])
								 +(wij[ii*ND+kk]-wij[jj*ND+kk])*(ph[kk*ND*ND+i*ND+j]);	//[��(4.31)�̈ꕔ]
					}
					//pddtt+=-2.0*tij[ii][jj]/double(n00h[i][j])
					//	*(sum1-8.0/PI*eij[ii][jj]*sqrt(ph[ii][i][j]*ph[jj][i][j]));
					pddtt+=-2.0*tij[ii*ND+jj]/double(n00h[i*ND+j])
				  *(sum1-8.0/PI*eij[ii*ND+jj]*sqrt(ph[ii*ND*ND+i*ND+j]*ph[jj*ND*ND+i*ND+j]));
					//�t�F�[�Y�t�B�[���h�̔��W������[��(4.31)]
				}
				//ph2[ii][i][j]=ph[ii][i][j]+pddtt*delt;		//�t�F�[�Y�t�B�[���h�̎��Ԕ��W�i�z��@�j
				//if(ph2[ii][i][j]>=1.0){ph2[ii][i][j]=1.0;}//�t�F�[�Y�t�B�[���h�̕ψ�␳
				//if(ph2[ii][i][j]<=0.0){ph2[ii][i][j]=0.0;}
				ph2[ii*ND*ND+i*ND+j]=ph[ii*ND*ND+i*ND+j]+pddtt*delt;		//�t�F�[�Y�t�B�[���h�̎��Ԕ��W�i�z��@�j
				if(ph2[ii*ND*ND+i*ND+j]>=1.0){ph2[ii*ND*ND+i*ND+j]=1.0;}//�t�F�[�Y�t�B�[���h�̕ψ�␳
				if(ph2[ii*ND*ND+i*ND+j]<=0.0){ph2[ii*ND*ND+i*ND+j]=0.0;}
			}
		}//j
	}//i

	for(k=1;k<=nm;k++){
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				//ph[k][i][j]=ph2[k][i][j];	//�⏕�z�����z��Ɉړ�
				ph[k*ND*ND+i*ND+j]=ph2[k*ND*ND+i*ND+j];	//�⏕�z�����z��Ɉړ�
			}
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

//*********************************************************************
	//if(keypress()){return 0;}//�L�[�҂����
	time1=time1+1.;  if(time1<time1max){goto start;}//�ő�J�E���g���ɓ��B�������ǂ����̔��f

	end:;
  return 0;
}

//************ ������(�t�F�[�Y�t�B�[���h)�̐ݒ�T�u���[�`�� *************
void ini000(double *ph, int ND, int N)
{
	int i, j, k, l, it;		//����
	int ii, jj, kk;				//����
	int ip, im, jp, jm;		//����
	int x1, y1, x1h[10], y1h[10];//�����j�̍��W
	double sum1, t, r0, phi, r;
	//srand(time(NULL)); // ����������
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;

	x1h[1]=0.2*nd;   y1h[1]=0.2*nd;		//�����j�P�̍��W�ݒ�
	x1h[2]=0.75*nd;  y1h[2]=0.4*nd;		//�����j�Q�̍��W�ݒ�
	x1h[3]=0.5*nd;   y1h[3]=0.75*nd;	//�����j�R�̍��W�ݒ�

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//for(ii=1;ii<=nm-1;ii++){ph[ii][i][j]=0.0;}
			//ph[nm][i][j]=1.0;//nm�Ԗڂ̃t�F�[�Y�t�B�[���h���P�ɏ�����
			for(ii=1;ii<=nm-1;ii++){ph[ii*ND*ND+i*ND+j]=0.0;}
			ph[nm*ND*ND+i*ND+j]=1.0;//nm�Ԗڂ̃t�F�[�Y�t�B�[���h���P�ɏ�����
		}
	}

	r0=10.0;
	for(ii=1;ii<=nm-1;ii++){
		x1=x1h[ii]; y1=y1h[ii];
		//x1=nd*DRND(1); y1=nd*DRND(1);
		for(i=0;i<=ndm;i++){
			for(j=0;j<=ndm;j++){
				r=sqrt( (double(i-x1))*(double(i-x1))+(double(j-y1))*(double(j-y1)) );
				//if(r<=r0){ ph[ii][i][j]=1.0;  ph[nm][i][j]=0.0; } //�����j�ʒu�̃t�F�[�Y�t�B�[���h��ݒ�
				if(r<=r0){
					ph[ii*ND*ND+i*ND+j]=1.0;
					ph[nm*ND*ND+i*ND+j]=0.0;
				} //�����j�ʒu�̃t�F�[�Y�t�B�[���h��ݒ�
			}
		}
	}

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			//sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k][i][j]; }
			//for(k=1;k<=nm;k++){ ph[k][i][j]=ph[k][i][j]/sum1; }//�t�F�[�Y�t�B�[���h�̋K�i���␳
			sum1=0.0; for(k=1;k<=nm;k++){ sum1+=ph[k*ND*ND+i*ND+j]; }
			for(k=1;k<=nm;k++){ ph[k*ND*ND+i*ND+j]=ph[k*ND*ND+i*ND+j]/sum1; }//�t�F�[�Y�t�B�[���h�̋K�i���␳
		}
	}

}

//************ �f�[�^�ۑ��T�u���[�`�� *******************************
void datsave(double *ph, int ND, int N)
{
	FILE		*stream;		//�X�g���[���̃|�C���^�ݒ�
	int 		i, j, k;		//����
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;

	stream = fopen("test.dat", "a");	//�������ސ�̃t�@�C����ǋL�����ŃI�[�v��
	fprintf(stream, "%e  \n", time1);	//�v�Z�J�E���g���̕ۑ�
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=nm;k++){
				//fprintf(stream, "%e   ", ph[k][i][j]);//�t�F�[�Y�t�B�[���h�̕ۑ�
				fprintf(stream, "%e   ", ph[k*ND*ND+i*ND+j]);//�t�F�[�Y�t�B�[���h�̕ۑ�
			}
		}
	}
	fprintf(stream, "\n");	//���s�̏�������
	fclose(stream);					//�t�@�C�����N���[�Y
}

void datsave_paraview(double *ph, int ND, int N)
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
	for(k=1;k<=nm;k++){
		sprintf(fName,"mpf0_N%03d_result%06d.vtk",k,iout);
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
				//fprintf(stream, "%e   ", ph[k][i][j]);//Save phase field
				//fprintf(fp,"%10.6f\n", ph[k][i][j]);
				fprintf(fp,"%10.6f\n", ph[k*ND*ND+i*ND+j]);
				pht[i*ND+j]+=ph[k*ND*ND+i*ND+j]*float(k);
			}
		}
		fclose(fp);
	}
	
	for(k=0;k<=0;k++){
		sprintf(fName,"mpf0_N%03d_result%06d.vtk",k,iout);
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
		fclose(fp);
	}
	free(pht);
}
//*********** �f�[�^���̓T�u���[�`�� **************************
void datin(double *ph, int ND, int N)
{
	FILE		*datin0;//�X�g���[���̃|�C���^�ݒ�
	int 		i, j, k;//����
	int nd=ND, ndm=ND-1, nm=N-1, nmm=N-2;

	datin0 = fopen("test.dat", "r");//�ǂݍ��݌��̃t�@�C�����I�[�v��
	fscanf(datin0, "%lf", &time1);	//�v�Z�J�E���g���̓ǂݍ���
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=nm;k++){
				//fscanf(datin0, "%lf  ", &ph[k][i][j]);//�t�F�[�Y�t�B�[���h�̓ǂݍ���
				fscanf(datin0, "%lf  ", &ph[k*ND*ND+i*ND+j]);//�t�F�[�Y�t�B�[���h�̓ǂݍ���
			}
		}
	}
	fclose(datin0);//�t�@�C�����N���[�Y

}