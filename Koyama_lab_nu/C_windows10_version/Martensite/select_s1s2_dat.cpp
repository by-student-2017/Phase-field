#include <stdio.h>
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <string.h>

#include "wingxa.h"

#define ND 128
#define IG 7
#define INXY 400

	int nd=ND, nd2=ND/2, ndm=ND-1;
	double PI=3.14159, time1;
	double s1h[ND][ND], s2h[ND][ND];

	void graph_s1();

/******* Main program ****************************************************/
int main(void)
{

	int   	i, j;
	int   	ii=0;
	char	fname[100];
	FILE	*datin0, *stream;
	int 	idat[30];
	char	ida_file0[30], ida_file1[30], ida_file2[30];

	strcpy(ida_file1, "./tmp/Mini_");
	strcpy(ida_file0, ida_file1);

	idat[0]=200;
	idat[1]=400;
	idat[2]=2000;
	//idat[3]=2800;
	//idat[4]=3000;
	//idat[5]=3000;
	//idat[6]=8000;
	//idat[7]=3000;
	//idat[8]=4000;
	//idat[9]=5000;
	//idat[10]=6000;
	//idat[11]=8000;
	//idat[12]=10000;

	printf("input file name = ? "); scanf("%s", fname);

	datin0 = fopen(fname, "r");

 	gwinsize(INXY,INXY); ginit(2); gsetorg(0,0);

start: ;

	fscanf(datin0, "%lf", &time1);
	printf("time %e \n",time1);

	for(i=0;i<=ndm;i++){
		for(j=0;j<=ndm;j++){
				fscanf(datin0, "%lf %lf", &s1h[i][j], &s2h[i][j]);
		}
	}

	if(fabs(time1-(float)idat[ii])<=1.){
		graph_s1();
		itoa(idat[ii], ida_file2, 10);
		strcat(ida_file1,ida_file2);
		strcat(ida_file1,".dat");

		stream = fopen(ida_file1, "a");
			fprintf(stream, "%e \n", time1);
			for(i=0;i<=ndm;i++){
				for(j=0;j<=ndm;j++){
					fprintf(stream, "%e  %e  ", s1h[i][j], s2h[i][j]);
				}
			}
		fclose(stream);

		ii=ii+1;
		strcpy(ida_file1, ida_file0);

		if(keypress()){return 0;}

	}

	if (feof(datin0)==0) {goto start;}
   	fclose(datin0);

  return 0;

end:;
}

//******* phase field‚Ì•`‰æ ***************************************
void graph_s1()
{
	int i, j, ii, jj;
	double col, col_R, col_G, col_B, col_RG;
	int ixmin=0, iymin=0, igx, igy, irad0;
	double c, x, xmax, xmin, y, ymax, ymin, rad0, dia0;
	int ixmax=INXY;
	int iymax=INXY;

	gcls(); //‰æ–ÊƒNƒŠƒA
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


