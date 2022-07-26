#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <string.h>

#include "wingxa.h"

#define ND 128
#define INXY 400

	int nd=ND, ndm=ND-1, nd2=ND/2;
	double PI=3.14159, time1;
	double s1h[ND][ND], s2h[ND][ND];

	void graph_s1();			//結晶度場描画サブル－チン

/******* Main program ****************************************************/
int main(void)
{
	int   	i, j;
	char	fname[100];
	FILE	*datin0;

	int   	ii=1;
	char	bmp_file0[100], bmp_file1[100], bmp_file2[100];
	int 	bmp_n[100];

//---data---------------------------------------

	strcpy(fname, "test2.dat");
	strcpy(bmp_file1, "mar2_");
	strcpy(bmp_file0, bmp_file1);

	bmp_n[1]=0;
	bmp_n[2]=100;
	bmp_n[3]=200;
	bmp_n[4]=300;
	bmp_n[5]=400;
	bmp_n[6]=500;
	bmp_n[7]=600;
	bmp_n[8]=700;
	bmp_n[9]=800;
	//bmp_n[10]=1600;
	//bmp_n[11]=1800;
	//bmp_n[12]=2000;

	//bmp_n[1]=0;
	//bmp_n[2]=50;
	//bmp_n[3]=100;
	//bmp_n[4]=200;
	//bmp_n[5]=500;
	//bmp_n[6]=1000;
	//bmp_n[7]=2000;
	//bmp_n[8]=5000;
	//bmp_n[9]=10000;
	//bmp_n[10]=20000;
	//bmp_n[11]=40000;
	//bmp_n[12]=50000;
	//bmp_n[13]=80000;
	//bmp_n[14]=100000;

	//printf("input file name = ? ");
	//scanf("%s", fname);

	datin0 = fopen(fname, "r");

 	gwinsize(INXY,INXY); ginit(1); gsetorg(0,0);

start: ;

	fscanf(datin0, "%lf", &time1);

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fscanf(datin0, "%lf %lf", &s1h[i][j], &s2h[i][j]);
		}
	}

	if((int)time1==bmp_n[ii]){
		graph_s1();
		itoa(bmp_n[ii], bmp_file2, 10);
		strcat(bmp_file1,bmp_file2);
		strcat(bmp_file1,".bmp");
		save_screen(bmp_file1);
		ii=ii+1;
		strcpy(bmp_file1, bmp_file0);

		if(keypress()){return 0;}

	}

	/*datsave();*/

	if (feof(datin0)==0) {goto start;}
   fclose(datin0);

  return 0;

end:;
}

//******* phase fieldの描画 ***************************************
void graph_s1()
{
	int i, j, ii, jj;
	double col, col_R, col_G, col_B, col_RG;
	int ixmin=0, iymin=0, igx, igy, irad0;
	double c, x, xmax, xmin, y, ymax, ymin, rad0, dia0;
	int ixmax=INXY;
	int iymax=INXY;

	gcls(); //画面クリア
	xmin=0.; xmax=1.; ymin=0.; ymax=1.;

	printf("time %f\n",time1);
	dia0=1.0/nd;
	rad0=dia0/2.0;
	irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;

	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			x=rad0+dia0*i;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
			y=rad0+dia0*j;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
			ii=i; jj=j; if(i==nd){ii=0;} if(j==nd){jj=0;}

			col_R=s1h[ii][jj];
			col_G=s2h[ii][jj];
			col_RG=col_R+col_G;
			if(col_RG>1.){col_RG=1.;}
			col_B=1.-col_RG;

			if(col_R>=0.999){col_R=1.;} if(col_R<=0.001){col_R=0.;}
			if(col_G>=0.999){col_G=1.;} if(col_G<=0.001){col_G=0.;}
			if(col_B>=0.999){col_B=1.;} if(col_B<=0.001){col_B=0.;}

			gcolor((int)(255*col_R),(int)(255*col_G),(int)(255*col_B));
			grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);
		}
	}
	//swapbuffers();
}
