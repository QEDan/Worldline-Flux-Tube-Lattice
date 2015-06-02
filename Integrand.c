#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <builtin_types.h> 
#include "intT.h" 
#include <fenv.h>
#include <signal.h>

  //electron mass squared
  double m2 = 1.0;
  //electron charge
  double e = 0.30282212;
  //double m2=1.0;
  //double e=1.0;

double EV(double T, void * p, int* WLlist);

double Exact(double T, double B, int fermion)
//Return the exact integrand for constant fields
{
  double TB = e*T*B;
  double Exact;
  if( B < 1.0e-40)
	Exact = 0.0;
  else
    if(fermion == 1)
	Exact = exp(-m2*T)/(T*T*T)*(TB/tanh(TB)-1.0-1.0/3.0*TB*TB);
    else
	if (TB < 50.0)
		Exact = exp(-m2*T)/(T*T*T)*(TB/sinh(TB)-1.0+1.0/6.0*TB*TB);
	else
		Exact = exp(-m2*T)/(T*T*T)*(1.0/6.0*TB*TB-1.0);
  //printf("%f %e\n",T, exp(-T));
  return Exact;
}

double bumpd(double x)
//The bump function
{
  if (abs(x) < 1.0)
	return exp(-1.0/(1.0-x*x));
  else
	return 0.0;	  
}

double fixfluxB(double rho2, double lambda2)
//Magnetic field strength for the fixflux profile
{
  float l = sqrt(lambda2);
  float lmlmin = (l - lmin)/(tubedist - lmin);
  if(sqrt(rho2) > tubedist/2.0) printf("Warning: inappropriate magnetic field being used. rho=%f\n", sqrt(rho2));
  return  2.0/(lambda2*q2)*(1.0-0.75*lmlmin)*bumpd(2.0*sqrt(rho2/lambda2))
	+ 3.0/(tubedist*tubedist)*lmlmin;
}

double periodicB(double rho2, double lambda2)
//magnetic field for the periodic flux tube profile
{
  double expr, rho, lambda, flraf, f, ceilraf;
  double expral, exprapl;
  rho = sqrt(rho2);
  lambda = sqrt(lambda2); 
  expr = exp(-rho/lambda);
  flraf = floor(rho/tubedist);
  ceilraf = ceil(rho/tubedist);
  if(rho < tubedist)
  {
	expral = 0.0;
	if((tubedist-rho)/lambda < 10.0)
		expral = exp(-(rho-tubedist)/lambda);
	f = expr/(tubedist*(1.0 + expr)*(1.0 + expr)); 
	if(rho > 0.0)
		f += expral/(rho*(1.0+expral)*(1.0+expral));
  }
  else
  {
        if(abs(fmod(rho,tubedist)) < 1.0e-6)
 	{
		f = flraf/(4.0*rho);
	}
	else
	{
	  expral = exp(-(rho-flraf*tubedist)/lambda);
	  exprapl = exp(-(rho-ceilraf*tubedist)/lambda);
          f = flraf*expral/(rho*(1.0+expral)*(1.0+expral)) +
		ceilraf*exprapl/(rho*(1.0+exprapl)*(1.0+exprapl));
        }
  }
  //f=1.0f/lambda2*exp(-1.0f*rho2/lambda2);
  //f=lambda2/8.0f;
  f=tubedist/(2.0*log(2.0)*lambda2)*f;
  //f=1.0;
  return f;
}

double MagField(double rho2, void* p)
//Compute the magnetic field. 
{
  double B;
  double sinrho;
  struct Wparams *params = (struct Wparams *) p;
  switch(profile){
	case step:
	  if(rho2<params->l2)
		B=1.0/params->l2;
	  else
		B=0.0;
	  break;
	case smooth:
	  B=params->l2/((params->l2+rho2)*(params->l2+rho2));
	  break;
	case quadratic:
	  if(rho2<params->l2)
		B=2.0/params->l2*(1.0-rho2/params->l2);
	  else
		B=0.0;
	  break;
	case gaussian:
	  if(rho2<100.0*params->l2)
	  	B=1.0/(params->l2*exp(rho2/params->l2));
	  else
		B=0.0;
	  break;
	case periodic:
	  B=periodicB(rho2,params->l2);
	  break;
	case spline:
	  sinrho=sin(pi*sqrt(rho2)/tubedist);
	  B=1.0/(params->l2*exp(sinrho*sinrho/params->l2));
	  break;
	case fixflux:
	  B = fixfluxB(rho2,params->l2);
	  break;
  }
  B=2.0*params->F/e*B;
  
  return B;
}

double smTscal(double B, double T)
//Effective action for T<<1.0 for scalar QED
{
	double B2=e*e*B*B;
	double B4=B2*B2;
	double B6=B2*B4;
	double T2=T*T;
	double T3=T2*T;
	double T4=T3*T;
	double result;
	result = 7.0/360.0*B4*T*(1.0-m2*T)+(147.0*B4*m2*m2-31.0*B6)/15120.0*T3
		+(31.0*B6*m2-49.0*B4*m2*m2*m2)/15120.0*T4;
	return result;
}

double smTferm(double B, double T)
//Effective action for T<<1.0 for fermionic QED
{
	double B2=e*e*B*B;
	double B4=B2*B2;
	double B6=B2*B4;
	double T2=T*T;
	double T3=T2*T;
	double T4=T3*T;
	double result;
	result = 1.0/45.0*B4*T*(m2*T-1.0) + (2.0*B6/945.0 - B4*m2*m2/90.0)*T3
		+(7.0*B4*m2*m2*m2-4.0*B6*m2)/1890.0*T4;
	//result = 1.0/45.0*B4*T*(T-1.0) + (2.0*B6/945.0 - B4/90.0)*T3
	//	+(7.0*B4-4.0*B6)/1890.0*T4;
	return result;
}

double meanigrand(double* igrand, int ngroups, double* igranderr)
//compute the mean and std. err. for igrand
{
  double meanig=0.0;
  int i;
  *igranderr=0.0;
  for(i=0;i<ngroups;i++)
  {
	meanig+=igrand[i];
  }
  meanig/=ngroups;
  for(i=0; i < ngroups; i++)
  {
	*igranderr+=(igrand[i]-meanig)*(igrand[i]-meanig);
  }
  *igranderr = sqrt(*igranderr);
  *igranderr /= ngroups;
  return meanig;
}


double* Tfunc(double T, void* p, int* order)
//The integrand of the proper time, T, integral
{
  struct Wparams *params = (struct Wparams *) p;
  double rho2 = params->xcm.x*params->xcm.x + params->xcm.y*params->xcm.y;
  double smallT;
  double B;
  double TB;
  double* igrand = malloc(params->ng*sizeof(*igrand));
  int groupsize = 128; //loops per group
  int* WLlist = malloc(groupsize*sizeof(*WLlist));
  int i, j;
  double rho = params->xcm.x;
  double l2 = params->l2;
  double sin2rho = sin(2.0*pi*rho/tubedist);
  double denominator;
  double igmean; double igerr;
  if(!igrand || !WLlist) printf("error Tfunc():  malloc failed\n");

  //printf("func: T=%f\n",T);
  B=MagField(rho2,p);
  TB=e*T*B;
  //printf("T=%f f'/f''=%f\n",T,(params->l2+params->xcm.x*params->xcm.x));
  //determine a suitable definition for small T
  switch(profile)
  {
	case step:
	  smallT=abs(rho*rho-l2);
          break;
        case smooth:
	  smallT= 0.5*(l2+rho*rho);
          break;
        case quadratic:
	  smallT= 0.5*(l2+rho*rho);
          break;
        case gaussian:
	  smallT=l2;
	  break;
        case periodic:
	  smallT=0.5*(l2+rho*rho);
          break;
	case spline:
	  if(rho/tubedist>1.0e-8)
	  	smallT=tubedist*tubedist*l2/(pi*pi*sin2rho*sin2rho);
	  else
		smallT=1.0e12;
	  break;
	case fixflux:
	  //small if T << l2, or if the loop is far from any bump functions.
	  if(rho*rho<l2/4.0)
		smallT = 0.25*l2;
	  else
		smallT = 0.25*l2+(rho-sqrt(l2)/2.0)*(rho-sqrt(l2)/2.0);
  }
  for(i = 0; i < params->ng; i++)
  {      
	//printf("func igrand[%d]=%f\n",i,igrand[i]);
	if(T < 0.1 && T < 0.01*smallT)
	//if(T<0.7)
  	{
		if(verbosity >= 2)
		  printf("using small T expression: %f %f %f\n", rho2, T, smallT);
		if(params->fermion == 1)
			igrand[i] = smTferm(B, T);
		else
			igrand[i] = smTscal(B, T);
  	}
        else if(T < 0.01*smallT)
	//use the exact expression
	{
		if(verbosity >= 2)
		  printf("using const. field expression %f %f %f\n", rho2, T, smallT);
		if(TB > 1.0e-6)
			igrand[i] = Exact(T, B, params->fermion);
		else
			igrand[i] = 0.0;	
	}
	else if(T > 49.0)
	//for large T, the contribution is exponentially supressed.
	{
		if(verbosity >= 2)
                  printf("using large T expression %f %f %f\n", rho2, T);
		if(params->fermion==1)
		{
			igrand[i]=-1.0*exp(-m2*T)/(3.0*T)*e*e*B*B;
			//if(B>1e-6)
			//	igrand[i]=exp(-m2*T)/(T*T*T)*(TB/tanh(TB)-1.0-1.0/3.0*TB*TB);
			//else 
			//	igrand[i]=0.0;
			//printf("igrand(%f)=%f\n",T,igrand[i]);
		}
		else
			igrand[i]=exp(-m2*T)/(6.0*T)*e*e*B*B;
	}
  	else
  	{
		for(j=0;j<groupsize;j++)
        	{
			//printf("func j=%d\n",j);
			WLlist[j]=order[j + i*groupsize];
        	}
		if(verbosity >= 3)
		  printf("EV= %e Exact = %e\n",EV(T,p, WLlist),T/sinh(T));
		if(params->fermion==1)
			igrand[i]= exp(-m2*T)/(T*T*T)*(EV(T, p, WLlist)-1.0-1.0/3.0*TB*TB);
		else
			igrand[i]= exp(-m2*T)/(T*T*T)*(EV(T, p, WLlist)-1.0+1.0/6.0*TB*TB);
		//igrand[i]=Exact(T,1.0);
		//return Exact(T,1.0);

  	}
 
  }
  if(T > 0.0f && verbosity >= 2){
	igmean = meanigrand(igrand, params->ng, &igerr);
	printf("WLvconstExpression: ");
  		printf("%e %e %e %e %e %e\n",
  		  rho2, T, Exact(T, B, params->fermion), smTscal(B, T), igmean, igerr);
  }
  free(WLlist);
  return igrand;
}

void signal_handler(int sig)
//Floating point signal handler
{
  printf( "Error: Floating point signal detected.\n");
  exit(1);
}

void EVMean(double *EV, float4 *Wsscal_h, float4 *Wsferm_h, int n, int *WL, double T, int fermion)
//Compute the expectation value and error of the Wilson loops
//WL is a list representing a group of wordline indices 
//to compute the expectation value of
{
  int i, WLi;
  double EVWlNb2, EVWlN, EVWlNb4;
  double SEWlNb2, SEWlN, SEWlNb4;
  double *EVpart1;
  double *EVpart2;
  double *EVpart4;
  double Nby4part;

  //Turn on floating point exceptions:                                                                                                                      
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);

  // Set the signal handler:                                                                                                                                 
  signal(SIGFPE, signal_handler);

  *EV = 0.0;
  EVWlNb2 = 0.0; EVWlN = 0.0; EVWlNb4 = 0.0;
  //allocate arrays to store the wilson loop of each worldline
  EVpart1 = (double *) malloc(n*sizeof(*EVpart1));
  EVpart2 = (double *) malloc(n*sizeof(*EVpart2));
  EVpart4 = (double *) malloc(n*sizeof(*EVpart4));
  if(!EVpart1 || !EVpart2 || !EVpart4) printf("EVonly(): malloc failed\n");
  //determine the sum of the contributions
  for(i = 0; i < n; i++){
	WLi=WL[i];
	//Compute the scalar part
	Nby4part = Wsscal_h[WLi].y + Wsscal_h[WLi].z;
	EVpart1[i] = cos(0.5*Wsscal_h[WLi].x+0.25*(Nby4part));
	EVpart2[i] = cos(Wsscal_h[WLi].x);
	EVpart4[i] = cos(0.5*(Nby4part));
	if(fermion == 1)
	{
		if(isinf(Wsferm_h[WLi].x) !=0 ||isinf(Wsferm_h[WLi].y) !=0 ||isinf(Wsferm_h[WLi].z) !=0 
		  || isnan(Wsferm_h[WLi].x) !=0 || isnan(Wsferm_h[WLi].y) !=0 ||isnan(Wsferm_h[WLi].z) !=0 )
		{
			printf("Warning: WSferm is infinite. Worldline=%d\n",WLi);
		}
		//Compute the fermion part
		//printf("Wsferm_h[%d].x=%f\n",WLi,Wsferm_h[WLi].x);
		//printf("Wsferm_h[%d].y=%f\n",WLi,Wsferm_h[WLi].y);
		//printf("Wsferm_h[%d].z=%f\n",WLi,Wsferm_h[WLi].z);
		Nby4part = Wsferm_h[WLi].y+Wsferm_h[WLi].z;
		if(abs(Nby4part) > 600.0 | abs(Wsferm_h[WLi].x) > 600.0)
		{
			printf("Warning: Large cosh argument.  T=%f, WLi=%d\n",T,WLi);
			EVpart1[i] *= 1.0e10;
			EVpart2[i] *= 1.0e10;
			EVpart4[i] *= 1.0e10;
		}
                else
		{
			EVpart1[i] *= cosh(0.5*(double)Wsferm_h[WLi].x+0.25*((double)Nby4part));
			EVpart2[i] *= cosh((double)Wsferm_h[WLi].x);
			EVpart4[i] *= cosh(0.5*((double)Nby4part));
			//EVpart1[i]*=cosh(T);
			//EVpart2[i]*=cosh(T);
			//EVpart4[i]*=cosh(T);
		}
		if(isnan(EVpart1[i]) != 0 || isnan(EVpart2[i]) != 0 || isnan(EVpart4[i]) != 0 
			|| isinf(EVpart1[i]) != 0 || isinf(EVpart2[i]) != 0 || isinf(EVpart4[i]) != 0) 
			printf("Warning: Infinity detected: WL# %d\n",WLi);
	}
	EVWlNb4 += EVpart4[i];
	EVWlNb2 += EVpart2[i];
	EVWlN   += EVpart1[i];
        //printf("EVpart1[%d]=%f\n",i,EVpart1[i]);
  }
  EVWlNb4 /=n;
  EVWlNb2 /=n;
  EVWlN  /=n;
  //Extrapolate to N=infinity
  *EV=8.0/3.0*EVWlN - 2.0*EVWlNb2 + 1.0/3.0*EVWlNb4;
  

  free(EVpart1);
  free(EVpart2);
  free(EVpart4);
  
}


