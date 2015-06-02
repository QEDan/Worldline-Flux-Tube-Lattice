#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <builtin_types.h> 
#include "intT.h" 

double e, m2;

double MagField( double rho2, void* p );

float4* getspline( float* rho2, double lambda2, double a, int Npoints );
double testspline( float rho2, double lambda2, double a, 
	int Npoints, float *rho2sp, float4 *flcoefs) ;

double* rhogrand(double rho, void* p)
//The Lagrangian Density Function
{
  int i, Npoints = 40;
  double *intout;

  struct Wparams *params = (struct Wparams *) p;
  //by symmetry, we may integrate rho along the x-axis
  params->xcm.x = rho;
  params->xcm.y = 0.0;
  params->xcm.z = 0.0;
  //perform the proper time integration
  if(verbosity >= 2)
  	printf("Beginning T integration\n");
  //cudaSetDevice();
  intout=integrateT(Npoints,p);
  if(verbosity >= 2)
  	printf("Returning from T integration\n");
  //for each group of worldlines
  for(i=0;i<params->ng;i++)
  {
	//multiply the integrand by rho
	intout[i] = rho*intout[i];
  }
  return intout;
}

void testvalrho(double *val, int n, double t)
//Outputs diagnostic information about the rho integrand
{
  int i;
  double avg=0.0;
  double SEM=0.0; //standard error in the mean
  //compute the average
  for(i=0;i<n;i++)
  {
	if(isnan(val[i])) printf("for rho=%f val[%d]=nan\n",(1.0-t)/t, i);
	avg+=val[i];
  }
  avg/=n;
  //compute the standard error
  for(i=0;i<n;i++)
  {
	SEM+=(val[i]-avg)*(val[i]-avg);
  }
  SEM=sqrt(SEM/((double)n*((double)n-1)));
  //rho Leff SEM-Leff
  printf("testvalrho: %f %f %f\n", (1.0-t)/t, avg, SEM);
}

double* transformRho(double rhot, void* p)
//Maps the interval [0,inf) to [0,1] using rho=(1-rhot)/rhot
{
	struct Wparams *params = (struct Wparams *) p;
	double* rhofunc;
	double rho=(1.0-rhot)/rhot;
	int i;
	
	rhofunc=rhogrand(rho,p);
	for(i=0;i<params->ng;i++)
	{
		rhofunc[i]=(rhofunc[i]/rhot)/rhot;

	}
	return rhofunc;

}

double* integrateRhotoa(int n, void * p)
//integrates rho from 0 to a/2 for periodic flux tubes
//This algorithm is the composite Simpson's rule as
//described in Burden and Faires "Numerical Analysis", 7th ed.
{
        double* val;
	double* XI0;
	double* XI1;
	double* XI2;
	double* rg;
	double* rg2;
	struct Wparams *params = (struct Wparams *) p;
	double rho;
	int inx, ig, i;
	val = (double *)calloc(params->ng, sizeof(*val));
	XI0 = (double *)calloc(params->ng, sizeof(*XI0));
	XI1 = (double *)calloc(params->ng, sizeof(*XI1));
	XI2 = (double *)calloc(params->ng, sizeof(*XI2));
	double h = tubedist/((double) 2*n);
	//temp replacement -> double h = 4.0*lmin/((double) 2*n);
	double rhot;
	if(!val) printf("integrate(): malloc failed\n");
	rg = rhogrand(0.0, p);
	rg2 = rhogrand(tubedist/2.0, p);
	//temp replacement -> rg2 = rhogrand(4.0*lmin/2.0, p);
	for(ig = 0; ig < params->ng; ig++)
	{
		val[ig] = 0.0;
		XI1[ig] = 0.0;
		XI2[ig] = 0.0;
		XI0[ig] = rg[ig] + rg2[ig];
	}
	for(inx = 0; inx < n; inx++)
	{
		printf("rho integrate inx=%d\n", inx);
		rho = h*((double)inx + 1.0);
		rg = rhogrand(rho, p);
		for(ig = 0; ig < params->ng; ig++)
		{
			if(inx%2 == 0) XI2[ig] += rg[ig];
			else XI1[ig] += rg[ig];
		}
		free(rg);
		
	}
	for(ig = 0; ig < params->ng; ig++)
	{
		val[ig] = h/3.0*(XI0[ig] + 2.0*XI2[ig] + 4.0*XI1[ig]);
	}	
	free(XI0);free(XI1);free(XI2);free(rg2);
        return val;
            
} 

double* integrateRho(int n, void * p)
//integrates rho from 0 to infinity for isolated flux tubes
//This algorithm is the composite Simpson's rule as
//described in Burden and Faires "Numerical Analysis", 7th ed.
{
        double* val;
	double* rhort;
	double* rhortphalf;
	double* rhortp1;
	struct Wparams *params = (struct Wparams *) p;
	int inx, ig, i;
	val=(double *)calloc(params->ng, sizeof(*val));
	double h = 1.0/((double) n);
	double rhot;
	if(!val) printf("integrate(): malloc failed\n");
	for(ig = 0; ig<params->ng; ig++)
	{
		val[ig] = 0.0;
	}
	for(inx = 0; inx < n; inx++)
	{
		printf("rho integrate inx=%d\n", inx);
		rhot = h*((double)inx+1.0);
		rhort = transformRho(rhot, p);
		rhortphalf = transformRho(h*((double)inx + 0.5), p);
		rhortp1 = transformRho(h*((double) inx+1.0), p);
		
		
		for(ig = 0; ig < params->ng; ig++)
		{
			//check for infinities and nans
			if(isinf(rhort[ig])) 
			{
				printf("Warning: rhort[%d]=%f.  Setting to zero.\n",ig,rhort[ig]);
				rhort[ig]=0.0;
			}
			if(isinf(rhortphalf[ig])) 
			{
				printf("Warning: rhortphalf[%d]=%f.  Setting to zero.\n",ig,rhortphalf[ig]);
				rhortphalf[ig]=0.0;
			}
			if(isinf(rhortp1[ig])) 
			{
				printf("Warning: rhortp1[%d]=%f.  Setting to zero.\n",ig,rhortp1[ig]);
				rhortp1[ig]=0.0;
			}
			//This equation extrapolates the three different Nppl spacings to Nppl=infinity
        		val[ig] += h/6.0*(
				rhort[ig] 
				+ 4.0*rhortphalf[ig]	
				+rhortp1[ig]); 
			if(isinf(val[ig])){
				printf("rhort[%d]=%f\n", ig, rhort[ig]);
				printf("rhortphalf[%d]=%f\n", ig, rhortphalf[ig]);
				printf("rhortp1[%d]=%f\n", ig, rhortp1[ig]);
				printf("val[%d]=%f\n", ig, val[ig]);
				
			}
			//printf("val[%d]=%f\n",ig,val[ig]);
		}
		//testvalrho(val, params->ng, rhot);
		
	}
	free(rhort);
	free(rhortphalf);
	free(rhortp1);
	
        return val;
            
} 

double PeriodicCA(double l2)
//Computes the spline Classical action
{
  double a = tubedist;
  double l = sqrt(l2);
  double expal = exp(a/l);
  double expalp1 = 1.0 + expal;
  double PCA;
  PCA = 2.0*a + 4.0*a/(expalp1*expalp1*expalp1) - l + 4.0*l/expalp1
	-2.0*(3.0*a+2.0*l)/(expalp1*expalp1) 
	+ l*log(16.0)-4.0*l*log(expalp1);
  return PCA;
}

double FixFluxCA(double l2)
//Computes the fixflux (periodic) classical action
{
  double q3 = 0.0187671;
  double a = tubedist;
  double l = sqrt(l2);
  double A = q3*(a*a/(l2*q2*q2));
  double lmlmin = (l-lmin)/(a-lmin);
  double FFCA;
  FFCA = 4.0*A + lmlmin*((9.0/4.0*A-9.0/2.0)*lmlmin - 6.0*A + 12.0);
  return FFCA;
}

double ClassAction(double F, double l2)
//Returns the Classical action for the given flux tube profile
{
  double CA;
  switch(profile){
	case step:
	  CA=-2.0*pi*F*F/(e*e*l2);
	  break;
	case smooth:
	  CA=-2.0*pi*F*F/(3.0*l2*e*e);
	  break;
	case quadratic:
	  CA=-8.0*pi*F*F/(3.0*l2*e*e);
	  break;
	case gaussian:
	  CA=-pi*F*F/(l2*e*e);
	  break;
	case periodic:
	  CA=-pi*sqrt(l2)/(24.0*l2*l2*log(2.0)*log(2.0))*PeriodicCA(l2);
	  break;
	case spline:
	  CA=1.0;
	  break;
	case fixflux:
	  CA = -pi*F*F/(e*e*tubedist*tubedist)*FixFluxCA(l2);
  }	
   //return -2.0/(3.0*pi*l2*e*e)*F*F;
   return CA;
}

double *EAction(int Npoints, void *p)
//Calls the integral over rho and returns 
//the effective action
//Npoints is the number of points used in rho integral
{
  int i;
  double *intout;
  //The average and standard error for the effective action
  double result, resultSE;
  struct Wparams *params = (struct Wparams *) p;

  //Determine which class of rho integral is appropriate
  if(profile == periodic || profile == spline || profile == fixflux)
	intout=integrateRhotoa(2*Npoints, p);
  else
  	intout=integrateRho(Npoints, p);
  
  result = 0.0;
  for(i = 0; i < params->ng; i++)
  {
	if(params->fermion == 1)
        	intout[i] /= 4.0*pi;
	else
		intout[i] /= -2.0*pi;
	result += intout[i];
	printf("rho: intout[%d]=%e\n", i, intout[i]);	
  }
  result /= params->ng;
  resultSE = 0.0;
  for(i = 0; i < params->ng; i++)
  {
	resultSE += (intout[i]-result)*(intout[i]-result);
  }
  resultSE = sqrt(resultSE/((double)params->ng*((double)params->ng-1.0)));
  if (verbosity >= 1)
  {
    printf("EAvsLambda: %e %e %e %e \n",sqrt(params->l2), ClassAction(params->F,params->l2), result, resultSE);
    printf("ActRatio: %e %e %e\n",sqrt(params->l2), result/ClassAction(params->F,params->l2), 
	  resultSE/ClassAction(params->F,params->l2));
  }
  return intout;
}

int main (int argc, char *argv[])
{
  //Number of worldlines, points per line, and dimensions in the worldline file
  const int Nl = 5120, Nppl = 1024, Dim = 2;
  const int fermion = 0; //1 => fermion, 0 => scalar
  int i, nBlocks, nThreads = 512;  //nThreads set to 256 for periodic fields.
  char *filename = "worldlines.5120.1024.2.dat"; //worldline text file
  //Arrays to hold values of the Wilson loops on the device and host
  float4 *Wsscal_d, *Wsferm_d;
  float4 *Wsscal_h, *Wsferm_h;
  float4 xcm;
  double exactX;
  //F=0.5*e*B_0*l2
  double F = 1.0, l2=lmin*lmin;//l2 = 0.02*tubedist*tubedist;
  double texe, T,X;
  float4 *worldlines_h; //worldlines on host
  float4 *worldlines_d; //worldlines on CUDA device
  cudaError_t errorcode;
  double result, resultSE; //mean and standard error of effective action
  double* intout;
  float4 *flcoefs_d; //spline coefficients on device
  float4 *flcoefs_h; //spline coefficients on host
  float *rho2sp_d;  //splined rho values on device
  float *rho2sp_h;  //splined rho values on host

  if(fermion==0)  printf("Warning: performing scQED calculation\n");
  switch(profile){
	case step: printf("F_lambda(rho2) = rho^2/lambda^2*theta(lambda^2-rho^2) + theta(rho^2-lambda^2)\n");
		break;
	case smooth: printf("F_lambda(rho2) = rho^2/(lambda^2 + rho^2)\n");
		break;
	case quadratic: printf("F_lambda(rho2) = 2*rho^2/lambda^2*(1-0.5*rho^2/lambda^2)*theta(lambda^2-rho^2)\n");
	        printf("+theta(rho^2-lambda^2)\n");
		break;
	case gaussian: printf("f_lambda(rho2) = 1-exp(-rho2/lambda2)\n");
		break;
	case periodic: printf("f_lambda(rho2) = Periodic B field\n");
		printf("a=%f\n",tubedist);
		nThreads = 256; //Periodic kernel requires more registers, so cannot support large nThreads.
		break;
	case spline: printf("f_lambda(rho2) = Spline data for periodic field\n");
	        printf("a=%f\n",tubedist);
		nThreads = 256; //Spline kernel requires more registers, so cannot support large nThreads.
	case fixflux: printf("f_lambda(rho2) = A*bump(2rho/lambda) + B\n ");
		printf("a=%f\n", tubedist);
		nThreads = 128; //Fixed Flux kernel requires more registers, so cannot support large nThreads.
		break;
	default: printf("Warning: Invalid field profile set\n");
		break;	
  }

  if(Nl%nThreads == 0)
	nBlocks = Nl/nThreads; 
  else 
	nBlocks = Nl/nThreads + 1;

  //parameters for the integration function
  struct Wparams params;

  //allocate memory for worldlines on host and device
  worldlines_h = (float4*)malloc(Nppl*nBlocks*nThreads*sizeof(*worldlines_h));
  cudaMalloc((void **)&worldlines_d,Nppl*nBlocks*nThreads*sizeof(*worldlines_d));
  printf("allocating memory \n");
  if(worldlines_h == NULL | worldlines_d == NULL)
  {
	fprintf(stderr, "out of memory 1\n");
	return(1);
  }
  //read worldlines in from file
  getwl(worldlines_h, filename, Nl, Nppl, Dim);
  printf("Copying worldlines to GPU device \n");
  //Copy worldlines to device
  errorcode=cudaMemcpy(worldlines_d, worldlines_h,
		nThreads*nBlocks*Nppl*sizeof(*worldlines_h), cudaMemcpyHostToDevice);
  if(errorcode > 0) printf("cudaMemcpy WLs: %s\n",cudaGetErrorString(errorcode));
  //allocate memory on the host
  Wsscal_h=(float4 *)malloc(nThreads*nBlocks*sizeof(*Wsscal_h));
  Wsferm_h=(float4 *)malloc(nThreads*nBlocks*sizeof(*Wsferm_h));
  //allocate memory on the device
  errorcode=cudaMalloc((float4**)&Wsscal_d, nThreads*nBlocks*sizeof(*Wsscal_d));
  if(errorcode > 0) printf("cudaMalloc Ws: %s\n",cudaGetErrorString(errorcode));
  errorcode=cudaMalloc((float4**)&Wsferm_d, nThreads*nBlocks*sizeof(*Wsferm_d));
  if(errorcode > 0) printf("cudaMalloc Ws: %s\n",cudaGetErrorString(errorcode));
  if(Wsscal_h==NULL | Wsscal_d==NULL |Wsferm_h==NULL | Wsferm_d==NULL)
  {
	fprintf(stderr,"out of memory 1\n");
	return(1);
  }


  //define some integration parameters
  params.Nl = Nl;
  params.Nppl = Nppl;
  params.F = F; 
  params.l2 = l2;
  params.xcm = xcm;
  params.Wsscal_h = Wsscal_h;
  params.Wsscal_d = Wsscal_d;
  params.Wsferm_h = Wsferm_h;
  params.Wsferm_d = Wsferm_d;
  params.worldlines = worldlines_d; 
  params.SE = 0.0;
  params.nBlocks = nBlocks;
  params.nThreads = nThreads;
  params.ng = 40;
  params.fermion = fermion;
  //intout=(float *)malloc(params.ng*sizeof(intout[0]));
  //if(!intout) printf("main(): malloc failed\n");
  printf("parameter dump:\n");
  printf("F=%f B_0=%f Bk\n l2=%f \n fermion=%d \n a=%f \n", F, 
	e*MagField(0.0,&params), l2, fermion, tubedist); 

  if(profile == spline)
  {
	rho2sp_h = (float *)malloc(Nspline*sizeof(*rho2sp_h));	
        flcoefs_h = getspline(rho2sp_h,l2,tubedist,Nspline);
	for(i = 0; i < Nspline; i++)
	{
		if(isnan(rho2sp_h[i]) || isnan(flcoefs_h[i].x) || 
		  isnan(flcoefs_h[i].y) ||isnan(flcoefs_h[i].y) ||isnan(flcoefs_h[i].y) )
			printf("is nan detected\n");
	        else if(isinf(rho2sp_h[i]) || isinf(flcoefs_h[i].x) || 
		  isinf(flcoefs_h[i].y) ||isinf(flcoefs_h[i].y) ||isinf(flcoefs_h[i].y) )
			printf("isinf detected\n");

	}
	errorcode = cudaMalloc((float**)&rho2sp_d, Nspline*sizeof(*rho2sp_d));
	if(errorcode > 0) printf("cudaMalloc rho2sp: %s\n", cudaGetErrorString(errorcode));
	errorcode = cudaMemcpy(rho2sp_d, rho2sp_h,
		Nspline*sizeof(*rho2sp_d),cudaMemcpyHostToDevice);
  	if(errorcode > 0) printf("cudaMemcpy rho2sp: %s\n", cudaGetErrorString(errorcode));

	errorcode = cudaMalloc((float4**)&flcoefs_d, Nspline*sizeof(*flcoefs_d));
	if(errorcode > 0) printf("cudaMalloc flcoefs: %s\n", cudaGetErrorString(errorcode));
	errorcode = cudaMemcpy(flcoefs_d, flcoefs_h,
		Nspline*sizeof(*flcoefs_d), cudaMemcpyHostToDevice);
  	if(errorcode > 0) printf("cudaMemcpy flcoefs: %s\n", cudaGetErrorString(errorcode));
	printf("TESTSPLINE: %f\n", testspline(1.0f, l2, tubedist, Nspline, rho2sp_h, flcoefs_h));

  }
  else
  //Spline data will not be used
  {
	flcoefs_d = NULL;
	rho2sp_d = NULL;
  }
  params.flcoefs=flcoefs_d;
  params.rho2sp=rho2sp_d;

  //compute the effective action in groups of worldlines
  //first argument is the number of points for the rho integral
  intout=EAction(20, &params);
  
  //compute the average
  result=0.0;
  for(i=0;i<params.ng;i++)
  {
	result+=intout[i];
	printf("rho: intout[%d]=%e\n",i,intout[i]);	
  }
  result/=params.ng;
  //Compute the standard error
  resultSE = 0.0;
  for(i = 0; i < params.ng; i++)
  {
	resultSE += (intout[i]-result)*(intout[i]-result);
  }
  resultSE = sqrt( resultSE/((double)params.ng*((double)params.ng-1.0)) );
  printf("EAction: %e +/- %e\n", result, resultSE);

  //Free the device memory
  cudaFree(Wsscal_d);
  cudaFree(Wsferm_d);
  cudaFree(worldlines_d);
  //Free the host memory
  free(Wsscal_h);
  free(Wsferm_h);
  free(worldlines_h);
  free(intout);

  if(profile == spline)
  {
   	free(flcoefs_h);
	free(rho2sp_h);
	cudaFree(flcoefs_d);
	cudaFree(rho2sp_d);
  }
  
     
  return 0;
}
