//This file contains mostly untested code for computing splined
//periodic flux tube configurations

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <builtin_types.h> 
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include "intT.h" 

//extern "C"
//__device__ __host__ float interp(float rho2, float *rho2sp, float4 *coefs);


int func (double t, const double y[], double f[],
           void *params)
{
    double *p = (double *)params;
    double lambda2 = p[0];
    double a = p[1];
    f[0] = 2.0/lambda2*t*exp(-pow(sin(pi*t/a),2)/lambda2);
     return GSL_SUCCESS;
}

int jac(double t, const double y[], double *dfdy, double dfdt[], void *params)
{
    double *p = (double *)params;
    double lambda2 = p[0];
    double a = p[1];
    
    gsl_matrix_view dfdy_mat 
        = gsl_matrix_view_array(dfdy, 1, 1);
    gsl_matrix * m = &dfdy_mat.matrix;

    double der=0.0;

    gsl_matrix_set (m, 0, 0, der);
    dfdt[0] = 2.0/lambda2*exp(-pow(sin(pi*t/a),2)/lambda2)*(1.0-
	2.0*pi/(a*lambda2)*t*sin(pi*t/a)*cos(pi/a*t));
    return GSL_SUCCESS;
}

float4* flspline(double *fli, float *rho2, int N)
//Computes the cubic spline coefficients for f_lambda(rho^2)
//uses Burden and Faires algorithm 3.5 (Clamped Cubic Spline)
//Step numbers refer to this textbook
{
  int i;
  double* h;
  double* l;
  double* mu;
  double* z;
  double* alpha;
  float4* flscoefs;
  h=(double *)malloc(N*sizeof(*h));
  l=(double *)malloc(N*sizeof(*l));
  mu=(double *)malloc(N*sizeof(*mu));
  z=(double *)malloc(N*sizeof(*z));
  alpha = (double *)malloc(N*sizeof(*alpha));
  flscoefs=(float4 *)malloc(N*sizeof(*flscoefs));

  
//Step 1
  for(i=0;i<N;i++)
  {
	flscoefs[i].x = fli[i];
  }
  for(i=0;i<N-1;i++) h[i]=rho2[i+1]-rho2[i];
//Step 2
  alpha[0]=3.0*(fli[1]-fli[0])/h[0];
  alpha[N-1]=-3.0*(fli[N-1]-fli[N-2])/h[N-2];
//Step 3
  for(i=1;i<N-1;i++)
  {
	//printf("%d \n",i);
	alpha[i]=3.0/h[i]*(fli[i+1]-fli[i])-3.0/h[i-1]*(fli[i]-fli[i-1]);
  }
//Step 4 - 
//Steps 4-7 solve a tridiagonal system using Crout factorization (Burden and Faires alg 6.7)
  l[0]=2.0*h[0];
  mu[0]=0.5;
  z[0]=alpha[0]/l[0];
//Step 5
  for(i=1;i<N-1;i++)
  {
	l[i]=2.0*(rho2[i+1]-rho2[i-1])-h[i-1]*mu[i-1];
	mu[i]=h[i]/l[i];
	z[i]=(alpha[i]-h[i-1]*z[i-1])/l[i];
  }
//Step 6
  l[N-1]=h[N-2]*(2.0-mu[N-2]);
  z[N-1]=(alpha[N-1]-h[N-2]*z[N-2])/l[N-1];
  flscoefs[N-1].z = z[N-1];
//Step 7
  for(i=N-2;i>-1;i--)
  {
	flscoefs[i].z = z[i]-mu[i]*flscoefs[i+1].z;
	flscoefs[i].y = (fli[i+1]-fli[i])/h[i] - h[i]*(flscoefs[i+1].z+2.0*flscoefs[i].z)/3.0;
	flscoefs[i].w = (flscoefs[i+1].z-flscoefs[i].z)/(3.0*h[i]);
  }
//Free allocated memory
  free(h);free(l);free(mu);free(z);free(alpha);
//Step 8
  return flscoefs;

  

}

float4* getspline(float* rho2, double lambda2, double a, int Npoints)
{

    double params[2];
    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;
    float4 *coefs;

    gsl_odeiv2_step * s  = gsl_odeiv2_step_alloc (T, 1);
    gsl_odeiv2_control * c  = gsl_odeiv2_control_y_new (1e-9,1e-7);
    gsl_odeiv2_evolve * e  = gsl_odeiv2_evolve_alloc (1);

    params[0] = lambda2;
    params[1] = a;


    gsl_odeiv2_system sys = {func, jac, 1, params};

    double rho = 0.0;
    double rhof = 20.0;
    double h = rhof/((double)Npoints);
    double y[1] = { 0.00 };
    double *flambda;
    int i=0;


    flambda=(double *)malloc(Npoints*sizeof(*flambda));
    while (i<Npoints)
    {
        int status = gsl_odeiv2_evolve_apply_fixed_step (e, c, s,
                                        &sys, 
                                        &rho, h, y);
 	if(status == GSL_FAILURE)
	{
		printf("ERROR:GSL_FAILURE\n");
		break;
	}
        else if (status != GSL_SUCCESS)
            break;
        rho2[i]=rho*rho;
	flambda[i]=y[0];
	i+=1;
	//rho+=h;
	
    }

    gsl_odeiv2_evolve_free (e);
    gsl_odeiv2_control_free (c);
    gsl_odeiv2_step_free (s);

    coefs=flspline(flambda,rho2,Npoints);
    free(flambda);
    return coefs;
}


double getinterp(float rho2, float *rho2sp, float4 *coefs)
{
  int j;
  int upperi=Nspline-1, loweri=0;
  float rho2diff=0.0;
  double flambda;

    //Discover which interval to look in using a binary search
  if(rho2<rho2sp[Nspline-1] && rho2>rho2sp[0])
//  if(0)
  {
   	while(upperi-loweri > 1)
	{
		if(rho2 >= rho2sp[(upperi+loweri)/2]) loweri=(upperi+loweri)/2;
		else upperi=(upperi+loweri)/2;
  	}
  	//interpolate using the jth interval
  	j=loweri;
	rho2diff=rho2-rho2sp[j];
	//rho2diff=0.0;
  	flambda= coefs[j].x+rho2diff*(coefs[j].y+rho2diff*(coefs[j].z+rho2diff*coefs[j].w));
  }
  else
  {
	j=0;
	rho2diff=0.0f;
	flambda=0.0f;
  }
  
  //*flambda= coefs[j].x+rho2diff*(coefs[j].y+rho2diff*(coefs[j].z+rho2diff*coefs[j].w));
  //*fprime = coefs[j].y+rho2diff*(2.0f*coefs[j].z + rho2diff*(3.0f*coefs[j].w));
  //*flambda=1.0f-exp(-1.0f*rho2);
  //*fprime=1.0f*exp(-1.0f*rho2);
  return flambda;
}


double testspline(float rho2, double lambda2, double a, int Npoints, float *rho2sp, float4 *flcoefs)
{
    printf("testspline: rho2=%f\n",(double)rho2);
    double flodeiv, flspline;
    double params[2];
    const gsl_odeiv2_step_type * T = gsl_odeiv2_step_rk8pd;
    float4 *coefs;

    gsl_odeiv2_step * s  = gsl_odeiv2_step_alloc (T, 1);
    gsl_odeiv2_control * c  = gsl_odeiv2_control_y_new (1e-9,1e-7);
    gsl_odeiv2_evolve * e  = gsl_odeiv2_evolve_alloc (1);

    params[0] = lambda2;
    params[1] = a;


    gsl_odeiv2_system sys = {func, jac, 1, params};

    double rho = 0.0;
    double rhof = sqrt(rho2);
    double h = rhof/((double)Npoints);
    double y[1] = { 0.00 };
    double *flambda;
    int i=0;
    flambda=(double *)malloc(Npoints*sizeof(*flambda));

    while(rho<rhof)
    {
    	
    	int status = gsl_odeiv2_evolve_apply (e, c, s,
                                        &sys, 
                                        &rho, rhof, &h, y);
    	if(status == GSL_FAILURE)
    	{
		printf("ERROR:GSL_FAILURE\n");
    	}
    }
    flodeiv=y[0];
    //rho+=h;
    flspline=(double)interp(rho2, rho2sp, flcoefs);
    return (flodeiv-flspline)/flodeiv;
}



