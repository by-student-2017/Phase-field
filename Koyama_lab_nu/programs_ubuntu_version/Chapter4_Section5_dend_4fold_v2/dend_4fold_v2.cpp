#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#define DRND(x) ((double)(x)/RAND_MAX*rand())//�����̊֐��ݒ�

//#define NDP 351							//�����v�Z�ɂ�����v�Z�̈��ӂ̕������{�P

	//int nd=NDP-1, ndm=NDP-2, nd2=(NDP-1)/2;//�v�Z�̈�̈�ӂ̍����u���b�N���And-1���`�And/2���`
	double delt;						//���ԍ���
	double PI=3.14159;					//�~����
	double RR=8.3145;					//�K�X�萔
	double Tini, Tm, time1;				//�������x�A�Z�_�A�v�Z�J�E���g��
	//double s1h[NDP][NDP], Th[NDP][NDP];	//�t�F�[�Y�t�B�[���h�A���x��
	int Nstep, iout;

	void ini000(double *s1h, double *Th, int NDP);	//������̐ݒ�T�u���|�`��
	void datsave(double *s1h, double *Th, int NDP);	//�f�|�^�ۑ��T�u���|�`��
	void datsave_paraview(double *s1h, double *Th, int NDP);//�f�|�^�ۑ��T�u���|�`��
	

//******* ���C���v���O���� ******************************************
int main(void)
{
	int NDP;						//�����v�Z�ɂ�����v�Z�̈��ӂ̕������{�P
	int nd, ndm, nd2;				//�v�Z�̈�̈�ӂ̍����u���b�N���And-1���`�And/2���`
	
	//double s1h2[NDP][NDP], Th2[NDP][NDP];	//��̕⏕�z��
	double s1kai, s1kais;			//�|�e���V����
	int   i, j, k, l;				//����
	int   ip, im, jp, jm;			//����
	double time1max;				//�v�Z�J�E���g���̍ő�l�i�v�Z�I���J�E���g�j

	double s1;						//�t�F�[�Y�t�B�[���h
	double s1ip, s1im, s1jp, s1jm;	//s1�̍��E�㉺�̃t�F�[�Y�t�B�[���h
	double s1ipjp, s1ipjm, s1imjp, s1imjm;//s1�̎΂߉��̃t�F�[�Y�t�B�[���h
	double s1ddtt;					//s1�̎��ԕω��ʁi���W�������̍��Ӂj

	double TT;						//���x��
	double Tip, Tim, Tjp, Tjm;		//TT�̍��E�㉺�̉��x
	double Tddtt;					//���x�̎��ԕω��ʁi�M�g�U�������̍��Ӂj

	double ep, ep1p, ep2p;			//���z�G�l���M�[�W���֘A�̕ϐ�
	double dx_s1, dy_s1, l_dxdy;	//�t�F�[�Y�t�B�[���h�̋�Ԕ����֘A�̕ϐ�
	double dxx_s1, dyy_s1, dxy_s1;	//�t�F�[�Y�t�B�[���h�̋�Ԕ����֘A�̕ϐ�

	double al;				//�v�Z�̈�̂P�ӂ̒���
	double dx, dy;			//�����i�q�T�C�Y(x����)
	double gamma;			//�E�ʃG�l���M�[���x
	double delta;			//�E�ʕ�
	double ram;				//��
	double bbb;				//�E�ʕ��Ɋ֌W����W��
	double j_fold;			//�ٕ������[�h
	double astre;			//�ٕ������x
	double cndct;			//�M�`����
	double speht;			//��M
	double rlate;			//���M
	double skine;			//�E�ʃJ�C�l�e�B�b�N�W��
	double th0;				//�D�搬�������̊p�x
	double th;				//�E�ʂ̖@�������̊p�x
	double aaa, aaac;		//���z�G�l���M�[�W��0
	double www, wwwc;		//�y�i���e�B�[���̃G�l���M�[���
	double pmobi, pmobic;	//�t�F�[�Y�t�B�[���h�̃��r���e�B
	double dtp;				//���ԍ���
	double dtt;				//���ԍ���
	double anois;			//�m�C�Y�̐U��
	double dF;				//�쓮��
	double dami1, dami2;	//��̌v�Z���X�L�b�v��������ݒ�Ɏg�p����ϐ�

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
	NDP     = int(data[0]);
	delt    = data[1];
	dx      = data[2];
	dy      = data[3];
	gamma   = data[4];
	ram     = data[5];
	j_fold  = data[6];
	astre   = data[7];
	cndct   = data[8];
	speht   = data[9];
	rlate   = data[10];
	Tm      = data[11];
	Tini    = data[12];
	skine   = data[13];
	th0     = data[14];
	aaac    = data[15];
	wwwc    = data[16];
	pmobic  = data[17];
	anois   = data[18];
	time1max= int(data[19]);
	Nstep   = int(data[20]);
	printf("---------------------------------\n");
	//
	nd=NDP-1;
	ndm=NDP-2;
	nd2=(NDP-1)/2;
	//
	double *s1h  = (double *)malloc(sizeof(double)*( NDP*NDP + NDP ));	//�t�F�[�Y�t�B�[���h
	double *Th   = (double *)malloc(sizeof(double)*( NDP*NDP + NDP ));	//���x��
	double *s1h2 = (double *)malloc(sizeof(double)*( NDP*NDP + NDP ));
	double *Th2  = (double *)malloc(sizeof(double)*( NDP*NDP + NDP ));
	//
	//printf("dtemp(1.0)=  ");	scanf(" %lf",&dtemp);	//dtemp=1.0;
	//printf("DELT(1.5)=  "); scanf(" %lf",&delt);	//delt=1.5;

	//dx=dy=30.0e-9;				//�����i�q�T�C�Y(x����)[m]
	delta=3.0*dx;      		//�E�ʕ�[m]
	//dx=dy=20.0e-9;			//�����i�q�T�C�Y(x����)[m]
	//delta=4.0*dx;    		//�E�ʕ�[m]
	al=dx*(double)nd;			//�v�Z�̈�[m]
	//gamma=0.37;       		//�E�ʃG�l���M�[[J/m^2]
	//ram=0.1;        			//��
	bbb=2.0*log((1.0+(1.0-2.0*ram))/(1.0-(1.0-2.0*ram)))/2.0;//�E�ʕ��Ɋ֌W����W��
	//j_fold=4.0;         	//�ٕ������[�h

	//astre=0.005;       	//�ٕ������x
	//astre=0.01;       	//�ٕ������x
	//astre=0.03;       		//�ٕ������x(�}4.8)
	//astre=0.05;       	//�ٕ������x

//    << pure Ni >>
	//cndct=84.01;      		//�M�`����[W/mK]
	//speht=5.42e+06;   		//��M[J/Km^3]
	//rlate=2.350e+09;  		//���M[J/m^3]
	//Tm=1728.0;      		//�Z�_[K]
	//Tini=1511.2;     		//�������x[K]
	//Tini=1424.5;     		//�������x[K]
	//skine=2.0;        		//�E�ʃJ�C�l�e�B�b�N�W��[m/Ks]
	//th0=0.0;         		//�D�搬�������̊p�x

//---- 4.5.2���Q�� ----------------------------------------------------
	//aaa=sqrt(3.0*delta*gamma/bbb);         	//���z�G�l���M�[�W��0
	aaa=sqrt(aaac*delta*gamma/bbb);         	//���z�G�l���M�[�W��
	//www=6.0*gamma*bbb/delta;               	//�y�i���e�B�[���̃G�l���M�[���
	www=wwwc*gamma*bbb/delta;               	//�y�i���e�B�[���̃G�l���M�[���
	//pmobi=bbb*Tm*skine/(3.0*delta*rlate); 	//�t�F�[�Y�t�B�[���h�̃��r���e�B
	pmobi=bbb*Tm*skine/(pmobic*delta*rlate); 	//�t�F�[�Y�t�B�[���h�̃��r���e�B

	dtp=dx*dx/(5.0*pmobi*aaa*aaa);	//���Ԃ�����
	dtt=dx*dx/(5.0*cndct/speht);	//���Ԃ�����
	if(dtp>dtt){delt=dtt;} else{delt=dtp;}
	printf("delt= %e \n", delt);
//-----------------------------------------------------------------

	//anois=0.1;	//�m�C�Y�̐U��

	time1=0.0;		//�v�Z���Ԃ̏����l
	//time1max=1.0+1.0e+08;	//�v�Z���Ԃ̍ő�l

//*** �����Z�x��̐ݒ�ƕ`��Window�\�� *****************************************

	ini000(s1h, Th, NDP);//������̐ݒ�

//**** �V�~�����[�V�����X�^�[�g ******************************
//Nstep = 100;
iout = -1;
start: ;

	if((((int)(time1) % Nstep)==0)){datsave(s1h, Th, NDP);}		//���J�Ԃ��J�E���g���ɏ��ۑ�
	if((((int)(time1) % Nstep)==0)){datsave_paraview(s1h, Th, NDP);}//���J�Ԃ��J�E���g���ɏ��ۑ�

//****** �t�F�[�Y�t�B�[���h����щ��x��̎��Ԕ��W  **************
	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			ip=i+1; im=i-1; jp=j+1; jm=j-1;
			if(i==nd){ip=ndm;}  if(i==0){im=1;}
			if(j==nd){jp=ndm;}  if(j==0){jm=1;}

			//s1=s1h[i][j];//�t�F�[�Y�t�B�[���h
			//s1ip=s1h[ip][j];  s1im=s1h[im][j];  s1jp=s1h[i][jp];  s1jm=s1h[i][jm];
			//s1ipjp=s1h[ip][jp]; s1ipjm=s1h[ip][jm];  s1imjp=s1h[im][jp];  s1imjm=s1h[im][jm];
			s1=s1h[i*NDP+j];//�t�F�[�Y�t�B�[���h
			s1ip=s1h[ip*NDP+j];  s1im=s1h[im*NDP+j];  s1jp=s1h[i*NDP+jp];  s1jm=s1h[i*NDP+jm];
			s1ipjp=s1h[ip*NDP+jp]; s1ipjm=s1h[ip*NDP+jm];  s1imjp=s1h[im*NDP+jp];  s1imjm=s1h[im*NDP+jm];

			//TT=Th[i][j];//���x��
			//Tip=Th[ip][j];  Tim=Th[im][j];  Tjp=Th[i][jp];  Tjm=Th[i][jm];
			TT=Th[i*NDP+j];//���x��
			Tip=Th[ip*NDP+j];  Tim=Th[im*NDP+j];  Tjp=Th[i*NDP+jp];  Tjm=Th[i*NDP+jm];

//----- �v�Z���X�L�b�v����ꍇ�̔��f ----------------------------------------------
			dami1=fabs(s1+s1ip+s1im+s1jp+s1jm);  dami2=fabs(TT+Tip+Tim+Tjp+Tjm-5.0*Tini);
			if( (dami1<=1.0e-20)&&(dami2<=1.0e-20) ){
				//s1h2[i][j]=s1h[i][j];  Th2[i][j]=Th[i][j];
				s1h2[i*NDP+j]=s1h[i*NDP+j];  Th2[i*NDP+j]=Th[i*NDP+j];
				goto dami;
			}
//---------------------------------------------------------------------------------

			dx_s1=(s1ip-s1im)/2.0/dx;  				//�t�F�[�Y�t�B�[���h�̋�ԂP�K����
			dy_s1=(s1jp-s1jm)/2.0/dy;
			dxx_s1=(s1ip+s1im-2.0*s1)/dx/dx;		//�t�F�[�Y�t�B�[���h�̋�ԂQ�K����
			dyy_s1=(s1jp+s1jm-2.0*s1)/dy/dy;
			dxy_s1=(s1ipjp+s1imjm-s1imjp-s1ipjm)/4.0/dx/dy;
			th=atan(dy_s1/(dx_s1+1.0e-20));			//�E�ʂ̖@�������̊p�x[��(4.24)]

			ep=aaa*(1.0+astre*cos(j_fold*(th-th0)));	//���z�G�l���M�[�W���̕�����[��(4.23)]
			ep1p=-aaa*astre*j_fold*sin(j_fold*(th-th0));	//ep�̊p�x�ɂ��P�K����
			ep2p=-aaa*astre*j_fold*j_fold*cos(j_fold*(th-th0));	//ep�̊p�x�ɂ��Q�K����

			s1kais=-ep*ep*(dxx_s1+dyy_s1)
					-ep*ep1p*((dyy_s1-dxx_s1)*sin(2.0*th)+2.0*dxy_s1*cos(2.0*th))
					 +0.5*(ep1p*ep1p+ep*ep2p)*(2.0*dxy_s1*sin(2.0*th)
					-dxx_s1-dyy_s1-(dyy_s1-dxx_s1)*cos(2.0*th));
			//���z�|�e���V����[��(4.29)�̈ꕔ]

			dF=15.0/(2.0*www)*rlate*(TT-Tm)/Tm*s1*(1.0-s1);//���w�I�쓮��[��(4.29)�̈ꕔ]
			s1kai=4.0*www*s1*(1.0-s1)*(0.5-s1+dF+anois*(DRND(1)-0.5));
			//���w�|�e���V����[��(4.29)�̈ꕔ]

			s1ddtt=-pmobi*(s1kai+s1kais);	//�t�F�[�Y�t�B�[���h�̔��W������[��(4.25)]

			//s1h2[i][j]=s1+s1ddtt*delt;	//�t�F�[�Y�t�B�[���h�̎��Ԕ��W(�z��@)
			//s1h2[i][j]=s1+s1ddtt*delt+anois*(DRND(1)-0.5)*s1*(1.0-s1);
			//if(s1h2[i][j]>=1.0){s1h2[i][j]=1.0;}//�t�F�[�Y�t�B�[���h�̕ψ�␳
			//if(s1h2[i][j]<=0.0){s1h2[i][j]=0.0;}
			s1h2[i*NDP+j]=s1+s1ddtt*delt;	//�t�F�[�Y�t�B�[���h�̎��Ԕ��W(�z��@)
			if(s1h2[i*NDP+j]>=1.0){s1h2[i*NDP+j]=1.0;}//�t�F�[�Y�t�B�[���h�̕ψ�␳
			if(s1h2[i*NDP+j]<=0.0){s1h2[i*NDP+j]=0.0;}

			Tddtt=( cndct*( (Tip+Tim-2.0*TT)/dx/dx+(Tjp+Tjm-2.0*TT)/dy/dy )
      	                    +30.0*s1*(1.0-s1)*s1*(1.0-s1)*rlate*s1ddtt )/speht;
			//�M�g�U������[��(4.30)]
			//Th2[i][j]=Th[i][j]+Tddtt*delt;		//���x��̎��Ԕ��W(�z��@)
			Th2[i*NDP+j]=Th[i*NDP+j]+Tddtt*delt;		//���x��̎��Ԕ��W(�z��@)

			dami:;
		}
	}

//**** ��̕⏕�z�����z��Ɉړ� ************************************
	for(i=0;i<=nd;i++){
	    //for(j=0;j<=nd;j++){  s1h[i][j]=s1h2[i][j]; Th[i][j]=Th2[i][j]; }
		for(j=0;j<=nd;j++){  s1h[i*NDP+j]=s1h2[i*NDP+j]; Th[i*NDP+j]=Th2[i*NDP+j]; }
	}
//*********************************************************************

	//if(keypress()){return 0;}//�L�[�҂����
	time1=time1+1.;  if(time1<time1max){goto start;}//���ԑ�������эő厞�Ԃɓ��B�������ǂ����̔��f

	end:;
  return 0;

}

//************ ������̐ݒ�T�u���[�`�� *************
void ini000(double *s1h, double *Th, int NDP)
{
	int i, j;
 	srand(time(NULL)); // ����������
	int nd=NDP-1, ndm=NDP-2, nd2=(NDP-1)/2;

	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			//s1h[i][j]=0.0;
			//if((i<10.)&&(j<10.)){s1h[i][j]=0.9;}
			//if((i*i+j*j)<20.){s1h[i][j]=0.9;}
			s1h[i*NDP+j]=0.0;
			if((i*i+j*j)<20.0){s1h[i*NDP+j]=0.9;}
		}
	}

	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			//Th[i][j]=Tini+s1h[i][j]*(Tm-Tini);
			Th[i*NDP+j]=Tini+s1h[i*NDP+j]*(Tm-Tini);
		}
	}

}

//************ �f�[�^�ۑ��T�u���[�`�� *******************************
void datsave(double *s1h, double *Th, int NDP)
{
	FILE		*stream;	//�X�g���[���̃|�C���^�ݒ�
	int 		i, j;			//����
	int nd=NDP-1, ndm=NDP-2, nd2=(NDP-1)/2;

	stream = fopen("test.dat", "a");	//�������ސ�̃t�@�C����ǋL�����ŃI�[�v��
	fprintf(stream, "%f\n", time1);		//�J�Ԃ��J�E���g���̕ۑ�
	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			//fprintf(stream, "%e  %e  ", s1h[i][j], Th[i][j]);//�Ǐ���̕ۑ�
			fprintf(stream, "%e  %e  ", s1h[i*NDP+j], Th[i*NDP+j]);//�Ǐ���̕ۑ�
		}
	}
	fprintf(stream, "\n");	//���s�̏�������
	fclose(stream);					//�t�@�C�����N���[�Y
}

void datsave_paraview(double *s1h, double *Th, int NDP)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j;
	int nd=NDP-1, ndm=NDP-2, nd2=(NDP-1)/2;
	
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
			//fprintf(stream, "%e  %e  ", s1h[i][j], Th[i][j]);//�Ǐ���̕ۑ�
			//fprintf(fp,"%10.6f\n", s1h[i][j]);
			fprintf(fp,"%10.6f\n", s1h[i*NDP+j]);
		}
	}
	fprintf(fp,"SCALARS Temperature float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(j=0;j<=ndm;j++){
		for(i=0;i<=ndm;i++){
			//fprintf(stream, "%e  %e  ", s1h[i][j], Th[i][j]);//�Ǐ���̕ۑ�
			//fprintf(fp,"%10.6f\n", Th[i][j]);
			fprintf(fp,"%10.6f\n", Th[i*NDP+j]);
		}
	}
	fclose(fp);
}
