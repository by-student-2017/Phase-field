#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>

#include "wingxa.h"

#define ND 64
#define INXY 400

	int nd=ND, ndm=ND-1, nd2=ND/2;
	double PI=3.14159, time1;
	double ch[ND][ND];

	void graph_c();			//”Z“xê•`‰æƒTƒuƒ‹|ƒ`ƒ“

/******* Main program ****************************************************/
int main(void)
{
	int i, j;
	char	fname[100];
	FILE	*datin0;

	printf("input file name = ? ");	scanf("%s", fname);
	datin0 = fopen(fname, "r");

 	gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);

start: ;

	fscanf(datin0, "%lf", &time1);
	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
			fscanf(datin0, "%lf  ", &ch[i][j]);
		}
	}

	graph_c();

	if(keypress()){return 0;}

	if (feof(datin0)==0) {goto start;}
  fclose(datin0);

  return 0;
end:;
}

//*******[‘Š•ª‰ğ‘gD‚Ì•`‰æ]**************************************************
void graph_c()
{
	int i, j, ii, jj;
	double col;
	double c, x, xmax, xmin, y, ymax, ymin, rad0;
	int ixmin=0, iymin=0, igx, igy, irad0;
	int ixmax=INXY;
	int iymax=INXY;

	gcls(); //‰æ–ÊƒNƒŠƒA
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
	swapbuffers();

} 

