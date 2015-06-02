//=========================================
// Main driver for the proper time integral
//=========================================

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <builtin_types.h> 
#include "intT.h" 



double* transformT(double t, void* p, int* order)
//Transforms the integral [0,infinity) to [0,1]
//x=(1-t)/t
{
	struct Wparams *params = (struct Wparams *) p;
	double* tfunc;
	double x = (1.0-t)/t;
	int i;
	
	tfunc = Tfunc(x,p, order);
	//loop over groups of worldlines
	for(i = 0; i < params->ng; i++)
	{
		tfunc[i] = (tfunc[i]/t)/t;

	}
	return tfunc;

}

void testval(double *val, int n, double t)
//used for outputting debugging information
{
  int i;
  double avg = 0.0;
  for(i = 0; i < n; i++)
  {
	if(isnan(val[i])) printf("for t=%f val[%d]=nan\n",t, i);
	avg += val[i];
  }
  avg /= n;


  printf("mean for Tt=%f, %f\n",t,avg);
}

void Tintout(double *val, int n, double rho, int fermion)
//Produces output about the T-integral
{
	int i;
	double avg = 0.0, SEM = 0.0;
  	for(i = 0; i < n; i++)
  	{
		if(isnan(val[i])) printf("for rho=%f val[%d]=nan\n",rho, i);
		avg += val[i];
  	}
  	avg /= n;
	for(i = 0; i < n; i++)
 	{
		SEM += (val[i]-avg)*(val[i]-avg);
  	}
  	SEM = sqrt(SEM/((double)n*((double)n-1.0)));
  	//rho Leff SEM-Leff
	
	if(fermion == 1)
  		printf("Leffvsrho: %e %e %e\n", rho, avg/(4.0*pi), SEM/(4.0*pi));
	else
  		printf("Leffvsrho: %e %e %e\n", rho, avg/(-2.0*pi), SEM/(2.0*pi));
}


int* shuffleWL(int Nl)
//Returns a vector representing a shuffled order of worldlines
{
  gsl_rng * r;
  int i;
  
  // select random number generator 
  r = gsl_rng_alloc (gsl_rng_mt19937);

  int* order=(int*)malloc(Nl*sizeof(*order));
          
  for (i = 0; i < Nl; i++)
  {
  	order[i] = i;
  }
          
  //gsl_ran_shuffle (r, order, Nl, sizeof (int));
  free(r);

  return order;
}



double* integrateT(int n, void * p)
//Performs the proper time integral using Simpson's Method
//n = number of points to use for Simpson's method
//p = Wparams parameters for the integrand
{
        double* val;    //value of the integral for each group
	double* trt;    //Integrands evaluated at three points
	double* trtphalf;
	double* trtp1;
	int* order;     //Stores the order of worldlines
	struct Wparams *params = (struct Wparams *) p;
	int inx, ig, i; 
	val = (double *)calloc(params->ng,sizeof(*val));
	double h = 1.0/((double) n);   //
	double t;
	if(!val) printf("integrate(): malloc failed\n");
	//initialize val to zero
	for(ig = 0; ig < params->ng; ig++)
	{
		val[ig] = 0.0;
	}
	//main integration loop
	for(inx = 0; inx < n; inx++)
	{
		printf("T integrate inx=%d\n",inx);
		//shuffle worldlines
		order=shuffleWL(params->Nl);
		t=h*((double)inx+1.0);
		trt=transformT(t,p,order);
		trtphalf=transformT(h*((double)inx + 0.5), p, order);
		trtp1=transformT(h*((double) inx+1.0), p, order);
		
		//A kludge to ignore possible infinities or nans
		for(ig = 0; ig < params->ng; ig++)
		{
			if(isinf(trt[ig])) 
			{
				printf("Warning: trt[%d]=%f.  Setting to zero.\n",ig,trt[ig]);
				trt[ig]=0.0;
			}
			if(isinf(trtphalf[ig])) 
			{
				printf("Warning: trtphalf[%d]=%f.  Setting to zero.\n",ig,trtphalf[ig]);
				trtphalf[ig]=0.0;
			}
			if(isinf(trtp1[ig])) 
			{
				printf("Warning: trtp1[%d]=%f.  Setting to zero.\n",ig,trtp1[ig]);
				trtp1[ig]=0.0;
			}

        		val[ig]+=h/6.0*(
				trt[ig] 
				+ 4.0*trtphalf[ig]	
				+trtp1[ig]); 
			if(isinf(val[ig])){
				printf("trt[%d]=%e\n",ig, trt[ig]);
				printf("trtphalf[%d]=%e\n",ig,trtphalf[ig]);
				printf("trtp1[%d]=%e\n",ig,trtp1[ig]);
				printf("val[%d]=%e\n",ig,val[ig]);
				
			}
			//printf("val[%d]=%e\n",ig,val[ig]);
		}
		testval(val, params->ng, t);
		free(order);
		
	}
	Tintout(val, params->ng, params->xcm.x, params->fermion);
	free(trt);
	free(trtphalf);
	free(trtp1);
	
        return val;
            
} 
