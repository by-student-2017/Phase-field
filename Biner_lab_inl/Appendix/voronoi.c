#include <stdio.h> //printf()
#include <math.h> //mod() and -lm

/* This program generates random polycrystalline microstructure by using
   Voronoi tessellation for the phase-field models */

void Discrete_Voronoi_diagram_2d();

void main(){
	
	//open output file
	//----- Assign unit names for the output files
	
	//Output file name for Voronoi tessellation
	FILE  *out=fopen("Voronoi_vertices.out","w");
	
	//Intermediate plot file name
	FILE *out1=fopen("plot_1.out","w");
	
	//Final graphical output file name
	FILE *out2=fopen("final_plot.p","w");
	
	/* File name that contains the coordinates of
	   the randomly generated points */
	FILE *out3=fopen("original_points.out","w");
	
	/* File name that conatins the Cartesian coordinates of
	   the corners of the simulation cell */
	FILE *out4=fopen("cell_1.out","w");
	
	/* Tabulated output file that will be used in
	   phase-field simulations */
	FILE *out5=fopen("grain_25.inp","w");
	//----- ----- ----- ---- -----
	
	int npoints=25;   //Number of grains in the simulation cell
	double xmax=32.0; //The length of the simulation cell in the x-direction
	double ymax=32.0; //The length of the simulation cell in the y-direction
	double x0=0.0;    //The origin of the coordinate system (x-direction)
	double y0=0.0;    //The origin of the coordinate system (y-direction)
	double extra=2.0; //A parameter will be used for trimming the Voronoi tessellation results
	
	
	//----- ----- ----- ---- ----- ----- ----- ----- ----
	//  generate random points for voronoi tessellation
	//----- ----- ----- ---- ----- ----- ----- ----- ----
	// Randomly generate the x- and y-coordinates of the points
	
	srand(12345);
	
	double x[npoints*9];
	double y[npoints*9];
	for(int ipoint=0;i<npoints;ipoint++){
		x[ipoint] = xmax*( (double)rand()/RAND_MAX );
		y[ipoint] = ymax*( (double)rand()/RAND_MAX );
	}
	
	//----- ----- ----- ---- ----- -----
	//  Duplicate the points for symmetry
	//----- ----- ----- ---- ----- -----
	/* Replicate the points for their nine periodic images,
	   for generation of periodic simulation cell in phase-field simulations */
	
	int jpoint;
	
	for(int ipoint=0;i<npoints;ipoint++){
		jpoint = npoints*1 + ipoint;
		x[jpoint]=x[ipoint];
		y[jpoint]=y[ipoint]-ymax;
	}
	//
	for(int ipoint=0;i<npoints;ipoint++){
		jpoint = npoints*2 + ipoint;
		x[jpoint]=x[ipoint]+xmax;
		y[jpoint]=y[ipoint]-ymax;
	}
	//
	for(int ipoint=0;i<npoints;ipoint++){
		jpoint = npoints*3 + ipoint;
		x[jpoint]=x[ipoint]+xmax;
		y[jpoint]=y[ipoint];
	}
	//
	for(int ipoint=0;i<npoints;ipoint++){
		jpoint = npoints*4 + ipoint;
		x[jpoint]=x[ipoint]+xmax;
		y[jpoint]=y[ipoint]+ymax;
	}
	//
	for(int ipoint=0;i<npoints;ipoint++){
		jpoint = npoints*5 + ipoint;
		x[jpoint]=x[ipoint];
		y[jpoint]=y[ipoint]+ymax;
	}
	//
	for(int ipoint=0;i<npoints;ipoint++){
		jpoint = npoints*6 + ipoint;
		x[jpoint]=x[ipoint]-xmax;
		y[jpoint]=y[ipoint]+ymax;
	}
	//
	for(int ipoint=0;i<npoints;ipoint++){
		jpoint = npoints*7 + ipoint;
		x[jpoint]=x[ipoint]-xmax;
		y[jpoint]=y[ipoint];
	}
	//
	for(int ipoint=0;i<npoints;ipoint++){
		jpoint = npoints*8 + ipoint;
		x[jpoint]=x[ipoint]-xmax;
		y[jpoint]=y[ipoint]-ymax;
	}
	
	//----- ----- ----- ---- ----- -----
	// Print-out original random points
	//----- ----- ----- ---- ----- -----
	//Output the coordinates of the points to file
	
	for(int i=0;i<npoints;ipoint++){
		fprintf(out3,"%4.6e %14.6e \n",x[i],y[i]);
	}
	
	//----- ----- ----- ---- ----- -----
	// Print the simulation cell coordinates
	//----- ----- ----- ---- ----- -----
	//Note: FILE *out4=fopen("cell_1.out","w");
	// Output the corner coordinates of the simulation cell to file
	
	fprintf(out4,"%4.6e %14.6e \n",x0,y0);
	fprintf(out4,"%4.6e %14.6e \n",xmax,y0);
	fprintf(out4,"%4.6e %14.6e \n",xmax,ymax);
	fprintf(out4,"%4.6e %14.6e \n",x0,ymax);
	fprintf(out4,"%4.6e %14.6e \n",x0,y0);
	
	
	//----- ----- ----- ---- ----- -----
	// generate voronoi diagram
	//----- ----- ----- ---- ----- -----
	/* Generate Voronoi tessellation by using
	   Matlab/Octave function voronoi */
	
	//voronoin(x,y,c,f);
	//
	double distance_x=1.0;
	double distance_y=1.0;
	int Nx=(int)(xmax-x0)/distance_x);
	int Ny=(int)(ymax-y0)/distance_y);
	//
	Discrete_Voronoi_diagram_2d(npoints, x, y,
		Nx, Ny, C, F,
		distance_x, distance_y);
	//
	for(int i=0;i<Nx;i++){
		for(int j=0;j<Ny;j++){
			//----- ----- ----- ---- ----- -----
			//Duplicate the c and the f
			//non:0, +:1, -:2
			//----- ----- ----- ---- ----- -----#1
			c[i+Nx*0][j+Ny*2]=c[i][j]; //non,-ymax
			f[i+Nx*0][j+Ny*2]=f[i][j]; //non,-ymax
			//----- ----- ----- ---- ----- -----#2
			c[i+Nx*1][j+Ny*2]=c[i][j]; //+xmax,-ymax
			f[i+Nx*1][j+Ny*2]=f[i][j]; //+xmax,-ymax
			//----- ----- ----- ---- ----- -----#3
			c[i+Nx*1][j+Ny*0]=c[i][j]; //+xmax,non
			f[i+Nx*1][j+Ny*0]=f[i][j]; //+xmax,non
			//----- ----- ----- ---- ----- -----#4
			c[i+Nx*1][j+Ny*1]=c[i][j]; //+xmax,+ymax
			f[i+Nx*1][j+Ny*1]=f[i][j]; //+xmax,+ymax
			//----- ----- ----- ---- ----- -----#5
			c[i+Nx*0][j+Ny*1]=c[i][j]; //non,+ymax
			f[i+Nx*0][j+Ny*1]=f[i][j]; //non,+ymax
			//----- ----- ----- ---- ----- -----#6
			c[i+Nx*2][j+Ny*1]=c[i][j]; //-xmax,+ymax
			f[i+Nx*2][j+Ny*1]=f[i][j]; //-xmax,+ymax
			//----- ----- ----- ---- ----- -----#7
			c[i+Nx*2][j+Ny*0]=c[i][j]; //-xmax,non
			f[i+Nx*2][j+Ny*0]=f[i][j]; //-xmax,non
			//----- ----- ----- ---- ----- -----#8
			c[i+Nx*2][j+Ny*2]=c[i][j]; //-xmax,-ymax
			f[i+Nx*2][j+Ny*2]=f[i][j]; //-xmax,-ymax
			//----- ----- ----- ---- ----- -----
		}
	}
	
	//----- ----- ----- ---- ----- -----
	// rearrange voronoin output
	// for their connectivity
	//----- ----- ----- ---- ----- -----
	
	nvelem = npionts;
	
	ncount=0;
	for(int i=0;i<nvelem;i++){
		flag=1;
		
		vnodes=f{i,;};
		nnode=size(vnodes,2);
		
		for(int j=0;j<nnode;j++){
			if(vnodes[j]==1){
				flag=0;
			}
		}
		
		if(flag==1){
			ncount = ncount + 1;
			for(int j=0;j<nnode;j++){
				lnods[ncount][j]=vnodes[j];
			}
		}
	}
	
	
	//----- ----- ----- ---- ----- -----
	// print voronoi results to file
	// to be viewed later in Voroni_vertices.out
	//----- ----- ----- ---- ----- -----
	
	for(int i=0;i<ncount;i++){
		fprintf(out, "# i %14.6e %d \n",i);
		
		// double lnods[ncount][nnode];
		// int ncount = sizeof(lnods)/sizeof(lnods[0]); //rows
		nnode  = sizeof(lnods[0)/sizeof(lnods[0][0]);  //colums
		
		for(int j=0;j<nnode;j++){
			kk = lnods[i][j];
			if(kk != 0){
				fprintf(out, "%14.6e %14.6e \n",c[kk][0],c[kk][1]);
			}
		}
		
		kk=lnods[i][0];
		fprintf(out, "%14.6e %14.6e \n",c[kk][0],c[kk][1]);
		
		fprintf(out, "\n");
	}
	
	
	//----- ----- ----- ---- ----- -----
	// Clip far outside voronoi elements
	// from the simulation cell
	//----- ----- ----- ---- ----- -----
	
	nelem=0;
	for(int i=0;i<ncount;i++){
		flag=0;
		for(int j=0;j<nnode;j++){
			kk=lnods[i][j];
			if(kk != 0){
				if( (c[kk][0] >= -extra) && (c[kk][0] <= xmax+extra) ){
					if( (c[kk][1] >= -extra) && (c[kk][1] <= ymax+extra) ){
						flag=1;
					}
				}
			}
		}//end for(j
		
		if(flag==1){
			nelem = nelem + 1;
			jnode = 0;
			for(int j=0;j<nnode;j++){
				kk=lnods[i][j];
				if(kk != 0){
					jnode = jnode + 1;
					lnods2[nelem][jnode] = lnods[i][j];
				}
			}//end for(j
			
			nnode2[nelem]=jnode;
		}//end if(flag
	}//end for(icount
	
	
	//----- ----- ----- ---- ----- -----
	// Assign grain numbers to
	// voronoi elements
	//----- ----- ----- ---- ----- -----
	
	twopi=8.0*atan(1.0);
	epsilon=1.0e-4;
	
	for(int isector=0;isector<9;isector++){
		for(int ipoint=0;ipoint<npoint;ipoint++){
			jpoint = isector*npoint + ipoint;
			
			for(int ielem=0;ielem<nelem;ielem++){
				theta=0.0;
				nnode=nnode2[ielem];
				
				for(int inode=0;inode<nnode;inode++){
					kk=lnods2[ielem][inode];
					
					xv1=c[kk][0];
					yv1=c[kk][1];
					
					jnode=inode+1;
					if(inode == nnode){
						jnode=1;
					}
					
					jj=lnods2[ielem][jnode];
					xv2=c[jj][0];
					yv2=c[jj][1];
					
					p2x=(xv1-x[jpoint]);
					p2y=(yv1-y[jpoint]);
					
					p1x=(xv2-x[jpoint]);
					p2y=(yv2-y[jpoint]);
					
					x1=sqrt(p1x*p1x + p1y*p1y);
					x2=sqrt(p2x*p2x + p2y+p2y);
					
					if(x1*x2 <= epsilon){
						theta=twopi;
					}else{
						tx1 = (p1x*p2x + p1y*p2y)/(x1*x2);
						
						if(fabs(tx1) >= 1.0){
							tx1=0.9999999999;
							
						}
						
						theta=theta+acos(tx1);
					}
				}//end for(inode
				
				if(fabs(theta-twopi) <= epsilon){
					igrain[ielem]=ipoint;
				}
			}//end for(ielem
		}//end for(ipoint
		
		//----- ----- ----- ---- ----- -----
		// print out as input file
		//----- ----- ----- ---- ----- -----
		
		nn1 = sizeof(c)/sizeof(c[0]); //rows
		nn2  = sizeof(c[0)/sizeof(c[0][0]);  //colums
		nnode=sizeof(lnods[0)/sizeof(lnods[0][0]);
		fprintf(out5, "%5d %5d %5d \n",(nn1-1),nnode,nelem);
		
		for(int i=1;i<nn1;i++){
			fprintf(out5, "%5d %14.6e %14.6e \n",i,c[i][0],c[i][1]);
		}
		
		for(int i=0;i<nelem;i++){
			fprintf(out5, "%5d ",i);
			for(int j=0;j<nnode;j++){
				fprintf(out5, "%5d ",lnods2[i][j]);
			}
			fprintf(out5, "%5d ",igrain[i]);
			fprintf(out5, "%5d ",lnods2[i][j]);
			fprintf(out5, "\n");
		}
		
		//----- ----- ----- ---- ----- -----
		//graphic output
		//----- ----- ----- ---- ----- -----
		
		for(int i=0;i<nelem;i++){
			fprintf(out1, "# i %d %d \n",i,nnode2[i]);
		}
		
		nnode  = sizeof(lnods2[0)/sizeof(lnods2[0][0]);  //colums
		
		ncount=0;
		xcod=0.0;
		ycod=0.0;
		
		for(int j=0;j<nnode;j++){
			kk=lnods2[i][j];
			fprintf(out1, "# %5d %5d \n",j,kk);
			
			if(kk != 0){
				fprintf(out1, "%14.6e %14.6e \n",c[kk][0],c[kk][1]);
				ncount=ncount+1;
				xcod=xcod+c[kk][0];
				ycod=ycod+c[kk][1];
			}
		}
		
		kk=lnods2[i][0];
		fprintf(out1, "%14.6e %14.6e \n",c[kk][0],c[kk][1]);
		fprintf(out1, "\n");
		
		xcod=xcod/ncount;
		ycod=ycod/ncount;
		
		fprintf(out2, "set label");
		fprintf(out2, '"');
		fprintf(out2, "%d",i);
		fprintf(out2, '"at');
		fprintf(out2, "%14.6e %14.6e \n",xcod,ycod);
		fprintf(out2, "\n";
	}
	
	fprintf(out2, 'plot "plot_1.out" wl, "cell_1.out wl, "original_points.out" \n');
	
	return;
}