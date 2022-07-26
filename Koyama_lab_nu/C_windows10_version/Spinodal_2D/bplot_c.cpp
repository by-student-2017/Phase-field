#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <string.h>

#include "wingxa.h"

#define ND 64
#define INXY 400

	int nd=ND, ndm=ND-1, nd2=ND/2;
	double PI=3.14159, time1;
	double ch[ND][ND];

	void graph_c();			//濃度場描画サブル－チン

/******* Main program ****************************************************/
main()
{
 	gwinsize(INXY,INXY); ginit(1); gsetorg(0,0);

	int   	i, j;
	char	fname[100];
	FILE	*datin0;

	int ii=1;
	char bmp_file0[100], bmp_file1[100], bmp_file2[100];
	int bmp_n[100];

//---data---------------------------------------
	strcpy(fname, "test.dat");
	strcpy(bmp_file1, "c_field_");
	strcpy(bmp_file0, bmp_file1);

	bmp_n[1]=0;
	bmp_n[2]=200;
	bmp_n[3]=400;
	bmp_n[4]=600;
	bmp_n[5]=800;
	bmp_n[6]=1000;
	bmp_n[7]=2000;
	bmp_n[8]=5000;
	bmp_n[9]=10000;
	bmp_n[10]=15000;
	bmp_n[11]=20000;
	bmp_n[12]=30000;
	bmp_n[13]=40000;
	bmp_n[14]=50000;
	//bmp_n[14]=60000;

	datin0 = fopen(fname, "r");

start: ;

	fscanf(datin0, "%lf", &time1);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fscanf(datin0, "%lf ", &ch[i][j]);
		}
	}

	if((int)time1==bmp_n[ii]){
		graph_c();
		//graph_s1();
		itoa(bmp_n[ii], bmp_file2, 10);
		strcat(bmp_file1,bmp_file2);
		strcat(bmp_file1,".bmp");
		save_screen(bmp_file1);
		ii=ii+1;
		strcpy(bmp_file1, bmp_file0);

		if(keypress()){return 0;}

	}

	if (feof(datin0)==0) {goto start;}
   fclose(datin0);

  return EXIT_SUCCESS;

end:;
}

//*******[相分解組織の描画]**************************************************
void graph_c()
{
	int i, j, ii, jj;
	double col;
	double c, x, xmax, xmin, y, ymax, ymin, rad0;
	int ixmin=0, iymin=0, igx, igy, irad0;
	int ixmax=INXY; // 最大x座標  Maximum-X
	int iymax=INXY; // 最大y座標  Maximum-Y

	gcls(); //画面クリア
	xmin=0.; xmax=1.;
	ymin=0.; ymax=1.;

	printf("time %f\n",time1);
	rad0=1./nd/2.;
	irad0=(ixmax-ixmin)/(xmax-xmin)*rad0+1;

	for(i=0;i<=nd;i++){
		for(j=0;j<=nd;j++){
			x=1./nd*i+rad0;  igx=(ixmax-ixmin)/(xmax-xmin)*(x-xmin)+ixmin;
			y=1./nd*j+rad0;  igy=(iymax-iymin)/(ymax-ymin)*(y-ymin)+iymin;
			ii=i; jj=j;  if(i==nd){ii=0;}  if(j==nd){jj=0;}
			col=1.-ch[ii][jj];
			if(col>=1.0){col=1.0;}  if(col<=0.0){col=0.0;}
			gcolor((int)(255*col),(int)(255*col),(int)(255*col));
			grect(igx-irad0,igy-irad0,igx+irad0,igy+irad0);
		}
	}
	//swapbuffers();

} 

