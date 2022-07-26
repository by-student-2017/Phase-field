#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <string.h>

#include "wingxa.h"

#define ND 512
#define INX 1000						//描画window１辺(x)のピクセルサイズ
#define INY 200						//描画window１辺の(y)ピクセルサイズ

	int nd=ND, nd2=ND/2, ndm=ND-1;
	double PI=3.14159, time1;
	double ch[ND];

	void graph_c();

/******* Main program ****************************************************/
int main(void)
{
	int   	i;
	int   	ii=0;
	char	fname[100];
	FILE	*datin0, *stream;
	int 	idat[20];
	char	ida_file0[30], ida_file1[30], ida_file2[30];
	double al, x;

	strcpy(ida_file1, "./tmp/c1D40_");
	strcpy(ida_file0, ida_file1);

	al=500.0;//[nm]

	//idat[0]=300000;

	idat[0]=0;
	idat[1]=1000;
	idat[2]=2000;
	idat[3]=5000;
	idat[4]=10000;
	idat[5]=20000;
	idat[6]=30000;
	idat[7]=50000;
	idat[8]=100000;
	idat[9]=200000;
	idat[10]=300000;

	//idat[10]=6000;
	//idat[11]=8000;
	//idat[12]=10000;

	printf("input file name = ? ");	scanf("%s", fname);

 	gwinsize(INX,INY); ginit(1); gsetorg(0,0);//描画Window表示

	datin0 = fopen(fname, "r");

start: ;

	fscanf(datin0, "%lf", &time1);
	printf("time %e \n",time1);
	for(i=0;i<=ndm;i++){ fscanf(datin0, "%lf ", &ch[i]); }

	if(fabs(time1-(float)idat[ii])<=1.){
		graph_c(); sleep(500);
		itoa(idat[ii], ida_file2, 10);
		strcat(ida_file1,ida_file2);
		strcat(ida_file1,".dat");

		stream = fopen(ida_file1, "a");
			//fprintf(stream, "%e \n", time1);
		for(i=0;i<=ndm;i++){
			x=al/double(nd)*(double)i;
			fprintf(stream, "%e  %e  \n", x, ch[i]); 
		}
		fclose(stream);

		ii=ii+1;
		strcpy(ida_file1, ida_file0);
	}

	if(keypress()){return 0;}
	if (feof(datin0)==0) {goto start;}
   	fclose(datin0);

  return 0;

end:;
}

//*******[相分解組織の描画]**************************************************
void graph_c()
{
	int i, ii, i1, ii1, i2, ii2;
	double col;
	double c, x, xmax, xmin, y, ymax, ymin, rad0, dx, dy;
	double gx1, gy1, gx2, gy2;
	int ixmin=0, iymin=0, igx1, igy1, igx2, igy2, irad0;
	int ixmax=INX;
	int iymax=INY;
	int idx, idy;

	xmin=0.0; xmax=1.0; dx=0.1;
	ymin=0.0; ymax=1.0; dy=0.1;
	idx=ixmax*(dx/(xmax-xmin))+0.5;
	idy=iymax*(dy/(ymax-ymin))+0.5;

	//gcls(); //画面クリア
	gcolor(255,255,255); grect(0,0,ixmax,iymax);//全画面白塗り

	printf("time %f\n",time1);
	rad0=1.0/nd/2.0;
	irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;

	for(i=0;i<=nd-1;i++){
		i1=i; i2=i+1;
		gx1=1.0/nd*i1+rad0;  gx2=1.0/nd*i2+rad0;
		ii1=i1; ii2=i2; if(ii2==nd){ii2=0;}
		gy1=ch[ii1];        gy2=ch[ii2];
		igx1=(ixmax-ixmin)/(xmax-xmin)*(gx1-xmin)+ixmin;
		igy1=(iymax-iymin)/(ymax-ymin)*(gy1-ymin)+iymin;
		igx2=(ixmax-ixmin)/(xmax-xmin)*(gx2-xmin)+ixmin;
		igy2=(iymax-iymin)/(ymax-ymin)*(gy2-ymin)+iymin;
		gcolor(255,0,0); gline(igx1, igy1, igx2, igy2);		//gcircle(igx, igy, 1);
	}

	gcolor(0,0,0); grectangle(ixmin, iymin, ixmax, iymax);
	for(i=0;i<ixmax;i+=idx){gline(i, iymin, i, iymax);}
	for(i=0;i<iymax;i+=idy){gline(ixmin, i, ixmax, i);}
	//swapbuffers();

}

