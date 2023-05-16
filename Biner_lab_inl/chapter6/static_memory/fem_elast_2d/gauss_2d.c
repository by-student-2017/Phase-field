/* This function evaluates the position of the sampling points and
   the associated weights for chosen order of numerical integration. */

/* Variable and array list
  ngaus: The order of numerical integration
  nnode: Number of nodes per element
  posgp[ngaus]: Position of sampling points
  weigp[ngaus]: Weighting factors at the sampling points
*/

void gauss_2d(int ngaus, int nnode,
	double *posgp, double *weigp){
	
	/* For three-node isoparametric triangular elements,
	   the position of sampling points and weights Radau rule. */
	if(nnode == 3){
		//
		if(ngaus == 1){
			posgp[0]=1.0/3.0;
			posgp[1]=1.0/3.0;
			weigp[0]=0.5;
		}
		
		if(ngaus == 3){
			posgp[0]=0.5;
			posgp[1]=0.5;
			posgp[2]=0.0;
			//
			posgp[3]=0.0;
			posgp[4]=0.5;
			posgp[5]=0.5;
			//
			weigp[0]=1.0/6.0;
			weigp[1]=1.0/6.0;
			weigp[2]=1.0/6.0;
		}
		
		if(ngaus == 7){
			posgp[0]=0.0;
			posgp[1]=0.5;
			posgp[2]=1.0;
			posgp[3]=0.5;
			posgp[4]=0.0;
			posgp[5]=0.0;
			posgp[6]=1.0/3.0;
			//
			posgp[7]=0.0;
			posgp[8]=0.0;
			posgp[9]=0.0;
			posgp[10]=0.5;
			posgp[11]=1.0;
			posgp[12]=0.5;
			posgp[13]=1.0/3.0;
			//
			weigp[0]=1.0/40.0;
			weigp[1]=1.0/15.0;
			weigp[2]=1.0/40.0;
			weigp[3]=1.0/15.0;
			weigp[4]=1.0/40.0;
			weigp[5]=1.0/15.0;
			weigp[6]=9.0/10.0;
		}
		//
	}//end if(nnode
	
	/* For four- and eight-node isoparametric elements,
	   the position of sampling points and weights by using
	   Gauss-Legendre rule. */
	if(nnode != 3){
		//
		if(ngaus == 2){
			posgp[0]=-0.57735026918963;
			weigp[0]=1.0;
		}
		if(ngaus != 2){
			posgp[0]=-0.7745966241483;
			posgp[1]=0.0;
			weigp[0]=0.55555555555556;
			weigp[1]=0.88888888888889;
		}
		
		kgaus = ngaus / 2;
		for(int igash=0;igash<kgaus;igash++){
			jgash = ngaus + 1 - igash;
			posgp[jgash] = -posgp[igash];
			weigp[jgash] =  weigp[igash];
		}
		//
	}
	
	return;
}