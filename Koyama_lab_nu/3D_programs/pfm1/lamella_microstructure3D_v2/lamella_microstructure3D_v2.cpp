#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

//#include "wingxa.h"

#define DRND(x) ((double)(x)/RAND_MAX*rand())

//#define ND 64
//#define INXY 400					//描画window１辺のピクセルサイズ

	//int nd=ND, ndm=ND-1, nd2=ND/2; 	//組織の分割、組織の分割-1
	double RR=8.3145;
	double PI=3.141592654;
	double c0, temp, time1;
	//double ch[ND][ND][ND];			//組織内の濃度デ－タ配列
	int iout;
	int flg;

	void ini_field(double *ch, int ND);	//初期濃度波設定サブル－チン
	//void graph_c();					//グラフ表示サブル－チン
	void datsave(double *ch, int ND);	//濃度デ－タ保存サブル－チン
	void datsave_paraview(double *ch, int ND);	//濃度デ－タ保存サブル－チン
	//void datin();						//初期濃度波読み込みサブル－チン

//******* Main program ******************************************************
int main(void)
{
	int ND;
	int nd, ndm, nd2; 				//組織の分割、組織の分割-1
	
	//double ck[ND][ND][ND];			//拡散ポテンシャル
	//double ch[ND][ND][ND];			//組織内の濃度デ－タ配列
	//double ch2[ND][ND][ND];
	double mu_chem, mu_str, mu_surf;
	int    i, j, k, l;
	int    Nstep;
	double time1max; 				//最大時間
	double amob_c;					//原子の易動度定数
	double c;						//濃度
	double cddtt;					//濃度の増分
	double kapa_c, kapa_c0;			//濃度勾配エネルギ－定数
	double al;						//計算領域
	double b1;						//差分ブロックサイズ
	double L0, L00;
	double delt;					//時間刻み
	double Mx, My, Mz;				//易動度

	double cip, cim, cjp, cjm, ckp, ckm;
	double ck1dev, ck2dev;
	double ck0, ckip, ckim, ckjp, ckjm, ckkp, ckkm;
	int    ip, im, jp, jm, kp, km;
	double sumc, dc0;
	double c_flu, c_flu0;			//濃度場の揺らぎの大きさ

//********************************************************************************
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
	c0      = data[2];
	Mx      = data[3];
	My      = data[4];
	Mz      = data[5];
	al      = data[6];	// [nm]
	amob_c  = data[7];
	time1max= data[8];
	temp    = data[9];	// [K]
	L00     = data[10];	//L0=L00/RR/temp;
	kapa_c0 = data[11];	//kapa_c=kapa_c0/b1/b1/RR/temp;
	c_flu   = data[12];
	flg     = int(data[13]);
	Nstep   = int(data[14]);
	printf("---------------------------------\n");
	//
	nd=ND;
	ndm=ND-1;
	nd2=ND/2;
	
	//double ch[(nstep+1)*nstep]; // [(nstep+1)*nstep]=[nstep*nstep+nstep], a[i][j]= a[i*n+j] for a[][n]
	//a[z * ySize * xSize + y * xSize + x], a[i][j][k]=a[i*n*m + j*n + k]
	double *ck  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//拡散ポテンシャル
	double *ch  = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));	//組織内の濃度デ－タ配列
	double *ch2 = (double *)malloc(sizeof(double)*( ND*ND*ND + ND*ND + ND ));
	
	//printf("DELT(0.005)=  "); scanf(" %lf",&delt);	//delt=0.005;
	//printf("ca (0.4)= "); scanf(" %lf",&c0);	//c0=0.4;
	//printf("Mx(1.0) = ");	scanf(" %lf",&Mx);	//移動度, Mx=0.01;	//標準値
	//printf("My(1.0) = ");	scanf(" %lf",&My);	//移動度, My=0.01;	/標準値
	//printf("Mz(1.0) = ");	scanf(" %lf",&Mz);	//移動度, Mz=1.0; 	//標準値

	//al=64.0;		// [nm]
	al=al*1.0e-9;	// [m]

	b1=al/nd;
	//amob_c=1.0;

	time1=0.0; 
	//time1max=30001.0;

	//temp=1000.0;	// [K]

	////L0=2.5e+04/RR/temp;
	//L0=1.0e+03/RR/temp;
	//kapa_c=5.0e-15/b1/b1/RR/temp;
	L0=L00/RR/temp;
	kapa_c=kapa_c0/b1/b1/RR/temp;
	
	//Mx=1.0;  						//易動度
	//My=Mx;  						//易動度
	//Mz=Mx;  						//易動度
	//My=0.0;  						//易動度
	//Mz=0.0;  						//易動度

	//c_flu=c_flu0=0.0;				//濃度場の揺らぎ係数, e.g., 0.0 or s0.1
	c_flu0 = c_flu;

//*************************************************************************
	ini_field(ch, ND);
	//gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);//描画Window表示

//**** start **************************************************************
//Nstep = 200;
iout = -1;
start: ;
	printf("time: %f \n", time1);
	if((((int)(time1) % Nstep)==0)) {datsave(ch,ND);}
	if((((int)(time1) % Nstep)==0)) {datsave_paraview(ch,ND);}

	//if(time1<1.0e4){c_flu=c_flu0;} else{c_flu=0.0;}

//******[拡散ポテンシャルの計算]********************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){

				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==ndm) {ip=0;}  if(i==0) {im=ndm;}
				if(j==ndm) {jp=0;}  if(j==0) {jm=ndm;}
				if(k==ndm) {kp=0;}  if(k==0) {km=ndm;}

				//c=ch[i][j][k];
				c=ch[i*ND*ND+j*ND+k];

			 	mu_chem=32.0*L0*c*(1.0-c)*(1.0-2.0*c); 		//化学ポテンシャル
			 	//mu_chem=L0*(1.0-2.0*c)+log(c)-log(1.0-c); //化学ポテンシャル

				//mu_surf=-2.*kapa_c*(ch[ip][j][k]+ch[im][j][k]
      			//				   +ch[i][jp][k]+ch[i][jm][k]
				//				   +ch[i][j][kp]+ch[i][j][km]-6.0*c);
				mu_surf=-2.*kapa_c*(ch[ip*ND*ND+j*ND+k]+ch[im*ND*ND+j*ND+k]
								   +ch[i*ND*ND+jp*ND+k]+ch[i*ND*ND+jm*ND+k]
      							   +ch[i*ND*ND+j*ND+kp]+ch[i*ND*ND+j*ND+km]
								   -6.0*c);
				
				//ck[i][j][k]=mu_chem+mu_surf;
				ck[i*ND*ND+j*ND+k]= mu_chem+mu_surf;
			}
		}
	}

//******[濃度場の時間変化]**************************************
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				ip=i+1; im=i-1; jp=j+1; jm=j-1; kp=k+1; km=k-1;
				if(i==ndm) {ip=0;}  if(i==0) {im=ndm;}
				if(j==ndm) {jp=0;}  if(j==0) {jm=ndm;}
				if(k==ndm) {kp=0;}  if(k==0) {km=ndm;}

				//cddtt=Mx*(ck[ip][j][k]+ck[im][j][k]-2.0*ck[i][j][k])
				//     +My*(ck[i][jp][k]+ck[i][jm][k]-2.0*ck[i][j][k])
				//     +Mz*(ck[i][j][kp]+ck[i][j][km]-2.0*ck[i][j][k]);
				cddtt=Mx*(ck[ip*ND*ND+j*ND+k]+ck[im*ND*ND+j*ND+k]-2.0*ck[i*ND*ND+j*ND+k])
					 +My*(ck[i*ND*ND+jp*ND+k]+ck[i*ND*ND+jm*ND+k]-2.0*ck[i*ND*ND+j*ND+k])
					 +Mz*(ck[i*ND*ND+j*ND+kp]+ck[i*ND*ND+j*ND+km]-2.0*ck[i*ND*ND+j*ND+k]);

				//ch2[i][j][k]=ch[i][j][k]+(cddtt+c_flu*(2.0*DRND(1.0)-1.0))*delt;
				ch2[i*ND*ND+j*ND+k] = ch[i*ND*ND+j*ND+k] + (cddtt+c_flu*(2.0*DRND(1.0)-1.0))*delt;
			}
		}
	}

//******[濃度場の収支の補正]**********************************
	sumc=0.;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//sumc+=ch2[i][j][k];
				sumc+=ch2[i*ND*ND+j*ND+k];
			}
	   }
	}
    dc0=sumc/nd/nd/nd-c0;

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
	    	for(k=0;k<=ndm;k++){
				//ch[i][j][k]=ch2[i][j][k]-dc0;
				//if(ch[i][j][k]>=1.0){ch[i][j][k]=1.0-1.0e-06;}
				//if(ch[i][j][k]<=0.0){ch[i][j][k]=1.0e-06;}
				ch[i*ND*ND+j*ND+k]=ch2[i*ND*ND+j*ND+k]-dc0;
				if(ch[i*ND*ND+j*ND+k]>=1.0){ch[i*ND*ND+j*ND+k]=1.0-1.0e-06;}
				if(ch[i*ND*ND+j*ND+k]<=0.0){ch[i*ND*ND+j*ND+k]=1.0e-06;}
	    	}
	    }
	}

//******[時間増加]*************************************************
	//if(keypress()){return 0;}
	time1=time1+1.;
	if (time1<time1max) {goto start;}

end:;
  return 0;
}

//************[初期濃度波]*****************************************
void ini_field(double *ch, int ND)
{
	int i, j, k, ii, jj, kk, PN, ipn;
	int x1, y1, z1, r0;
	//int flg;
	int nd=ND, ndm=ND-1, nd2=ND/2; 	//組織の分割、組織の分割-1
	double rnd0, sumc, r; 
 	srand(time(NULL)); // 乱数初期化
	int nd01, nd02, dd1;

//初期化
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
				for(k=0;k<=ndm;k++){
					//ch[i][j][k]=0.0;
					ch[i*ND*ND+j*ND+k]=0.0;
			}
		}
	}


// 1: スピーノーダル分解
// 2: 核形成-成長型分解
// 3: １個の球(体積分率はこちらで起算している)
// 4: １個のトーラス(体積分率はこちらで起算している)

	//flg=1;//以上記のいずれかを選択する。

	switch(flg){

	case 1:	//スピーノーダル分解
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
				for(k=0;k<=ndm;k++){
					//ch[i][j][k]=c0+0.1*(2.*DRND(1)-1.);
					ch[i*ND*ND+j*ND+k]=c0+0.1*(2.*DRND(1)-1.);
			}
		}
	}
	break;


	case 2:	//核形成
	r0=nd/20;	//径r0の核を置く
	//PN=50;	//核をPN個置く
	PN=c0/( 4./3.*PI*(1./20.)*(1./20.)*(1./20.) );

	for(ipn=1;ipn<=PN;ipn++){
		x1=nd*DRND(1); y1=nd*DRND(1); z1=nd*DRND(1);
		//while(ch[x1][y1][z1]!=0.0){x1=nd*DRND(1); y1=nd*DRND(1); z1=nd*DRND(1);}
		while(ch[x1*ND*ND+y1*ND+z1]!=0.0){x1=nd*DRND(1); y1=nd*DRND(1); z1=nd*DRND(1);}
		for(i=-r0;i<=(ndm+r0);i++){
			ii=i; if(i<0){ii=nd+i;}  if(i>ndm){ii=i-nd;}
			for(j=-r0;j<=(ndm+r0);j++){
				jj=j; if(j<0){jj=nd+j;}  if(j>ndm){jj=j-nd;}
				for(k=-r0;k<=(ndm+r0);k++){
					kk=k; if(k<0){kk=nd+k;}  if(k>ndm){kk=k-nd;}
					r=sqrt( (double(i-x1))*(double(i-x1))
						   +(double(j-y1))*(double(j-y1))
						   +(double(k-z1))*(double(k-z1)) );
					//if(r<=r0){ ch[ii][jj][kk]=1.0; } 
					if(r<=r0){ ch[ii*ND*ND+jj*ND+kk]=1.0; } 
				}
			}
		}
	}

	sumc=0.;
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i][j][k]; } } }
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i*ND*ND+j*ND+k]; } } }
		//c0=sumc/nd/nd/nd;
	break;


	case 3: //１個の球
	nd01=nd/3;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				dd1=(i-nd2)*(i-nd2)+(j-nd2)*(j-nd2)+(k-nd2)*(k-nd2);
				//if(dd1<(nd01*nd01)){ch[i][j][k]=0.9;} else{ch[i][j][k]=0.2;}
				if(dd1<(nd01*nd01)){ch[i*ND*ND+j*ND+k]=0.9;} else{ch[i*ND*ND+j*ND+k]=0.2;}
			}
		}
	}
	sumc=0.;
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i][j][k]; } } }
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i*ND*ND+j*ND+k]; } } }
		c0=sumc/nd/nd/nd;
	break;


	case 4: //１個のトーラス
	nd01=nd/3; nd02=nd/6;
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=nd2-4;k<=nd2+4;k++){
				dd1=(i-nd2)*(i-nd2)+(j-nd2)*(j-nd2);
				//if(dd1<(nd01*nd01)){ch[i][j][k]=0.9;} else{ch[i][j][k]=0.2;}
				//if(dd1<(nd02*nd02)){ch[i][j][k]=0.2;}
				if(dd1<(nd01*nd01)){ch[i*ND*ND+j*ND+k]=0.9;} else{ch[i*ND*ND+j*ND+k]=0.2;}
				if(dd1<(nd02*nd02)){ch[i*ND*ND+j*ND+k]=0.2;}
			}
		}
	}
	sumc=0.;
	//for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i][j][k]; } } }
	for(i=0;i<=ndm;i++){ for(j=0;j<=ndm;j++){ for(k=0;k<=ndm;k++){ sumc+=ch[i*ND*ND+j*ND+k]; } } }
		c0=sumc/nd/nd/nd;
	break;
	}//switch
}

//************[濃度波デ－タの保存]************************************
void datsave(double *ch, int ND)
{
	FILE		*stream;
	int 		i, j, k;
	int ndm=ND-1;

	//stream = fopen("test3D_lamella.dat", "w");
	stream = fopen("test3D_lamella.dat", "a");
	fprintf(stream, "%e\n", time1);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			for(k=0;k<=ndm;k++){
				//fprintf(stream, "%e  ", ch[i][j][k]);
				fprintf(stream, "%e  ", ch[i*ND*ND+j*ND+k]);
			}
		}
	}
	fprintf(stream, "\n");
	fclose(stream);
}

void datsave_paraview(double *ch, int ND)
{
	FILE	*fp;
	char	fName[256];
	int 	i, j, k;
	int ndm=ND-1;
	
	iout = iout + 1;
	sprintf(fName,"lamella_3D_result%06d.vtk",iout);
	fp = fopen(fName, "w");
	fprintf(fp,"# vtk DataFile Version 3.0 \n");
	fprintf(fp,"output.vtk \n");
	fprintf(fp,"ASCII \n");
	fprintf(fp,"DATASET STRUCTURED_POINTS \n");
	fprintf(fp,"DIMENSIONS %5d %5d %5d \n",(ndm+1),(ndm+1),(ndm+1));
	fprintf(fp,"ORIGIN 0.0 0.0 0.0 \n");
	fprintf(fp,"ASPECT_RATIO 1 1 1 \n");
	fprintf(fp,"POINT_DATA %16d \n",((ndm+1)*(ndm+1)*(ndm+1)));
	fprintf(fp,"SCALARS concentration float \n");
	fprintf(fp,"LOOKUP_TABLE default \n");
	for(k=0;k<=ndm;k++){
		for(j=0;j<=ndm;j++){
			for(i=0;i<=ndm;i++){
				//fprintf(fp,"%10.6f\n", ch[i][j][k]);
				fprintf(fp,"%10.6f\n", ch[i*ND*ND+j*ND+k]);
			}
		}
	}
	fclose(fp);
}
//***********[ファイルからの初期濃度波入力]****************************
//void datin()
//{
//	FILE		*datin0;
//	int 		i, j;

//	datin0 = fopen("test.dat", "r");
//	fscanf(datin0, "%lf", &time1);
//	for(i=0;i<=ndm;i++){
//		for(j=0;j<=ndm;j++){
//			fscanf(datin0, "%lf", &ch[i][j]);
//		}
//	}
//	fclose(datin0);
//}
