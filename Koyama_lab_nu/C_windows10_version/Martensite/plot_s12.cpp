#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#include "wingxa.h"

#define ND 128
#define INXY 400

	int nd=ND, ndm=ND-1, nd2=ND/2;
	double PI=3.14159, time1;
	double s1h[ND][ND], s2h[ND][ND];

	void graph_s1();			//場描画サブル－チン

/******* Main program ****************************************************/
int main(void)
{
	int   	i, j;
	char	fname[100];
	FILE	*datin0;

	printf("input file name = ? "); scanf("%s", fname);

	datin0 = fopen(fname, "r");

 	gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);

start: ;

	fscanf(datin0, "%lf", &time1);

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fscanf(datin0, "%lf  %lf", &s1h[i][j], &s2h[i][j]);
		}
	}

	graph_s1();

	if(keypress()){return 0;}

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
	swapbuffers();
}

