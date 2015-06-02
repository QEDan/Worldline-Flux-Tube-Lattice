//=========================================
// Kernel code for computing Wilson loops on the device
//=========================================
#include <builtin_types.h>
#include <cuda.h>
#include <stdio.h>
#include "intT.h"
#include <math.h>


extern "C" void EVMean(double *EV, float4 *Wsscal_h, float4 *Wsferm_h, int n, int *WL, double T, int fermion);
#if profile == fixflux
 #define THREADS_PER_BLOCK 128
 #define MY_KERNEL_MAX_THREADS THREADS_PER_BLOCK
 #define MY_KERNEL_MIN_BLOCKS 4
#else
 #define THREADS_PER_BLOCK 256
 #if __CUDA_ARCH__ >= 200
  #define MY_KERNEL_MAX_THREADS (2 * THREADS_PER_BLOCK)
  #define MY_KERNEL_MIN_BLOCKS 3
 #else
  #define MY_KERNEL_MAX_THREADS (2 * THREADS_PER_BLOCK)
  #define MY_KERNEL_MIN_BLOCKS 2
 #endif
#endif



extern "C"
__device__ float expint(const float x)
//Evaluates the exponential integral Ei(x)=-E1(-x) assuming x<-1.
//This algorithm is an abbreviated version of Numerical Recipes expint().
//See Chapter 6 on Special Functions. We assume (x > 1.0) and include only 
//the relevant code. This bit of code is Lentz's algorithm (section 5.2 of NR).
{
  const int MAXIT = 400;
  const float EPS = 1.0e-6;
  const float BIG = 1.0e10;
  int i;
  float a, b, c, d, del, h, ans;
  b = -x + 1.0f;
  c = BIG;
  d = 1.0f/b;
  h = d;
  for (i = 1; i <= MAXIT; i++)
  {
	a = -(float)i*(float)i;
	b += 2.0f;
	d = 1.0f/(a*d+b); //Denominators cannot be zero
	c = b + a/c;
	del = c*d;
	h *= del;
	if (fabsf(del-1.0f) <= EPS)
	{
		ans = -h*__expf(x);
		return ans;
	}
  }
  return 0.0f;

}


extern "C"
__host__ __device__ float interp(float rho2, float *rho2sp, float4 *coefs)
//interpolation function for periodic spline profile
{
  int j;
  int upperi = Nspline-1, loweri=0;
  float rho2diff = 0.0f;
  float flambda;

  #ifndef __CUDA_ARCH__
        #warning __CUDA_ARCH__ Undefined!
	//printf("Printf error\n");
  #else
	#warning __CUDA_ARCH__ defined!
  #endif


    //Discover which interval to look in using a binary search
  if(rho2 < rho2sp[Nspline-1] && rho2 > rho2sp[0])
//  if(0)
  {
   	while(upperi-loweri > 1)
	{
		if(rho2 >= rho2sp[(upperi+loweri)/2]) loweri=(upperi+loweri)/2;
		else upperi = (upperi+loweri)/2;
  	}
  	//interpolate using the jth interval
  	j = loweri;
	rho2diff = rho2-rho2sp[j];
	//rho2diff=0.0;
  	flambda = coefs[j].x+rho2diff*(coefs[j].y+rho2diff*(coefs[j].z+rho2diff*coefs[j].w));
  }
  else
  {
	flambda=0.0f;
  }
  
  //*flambda= coefs[j].x+rho2diff*(coefs[j].y+rho2diff*(coefs[j].z+rho2diff*coefs[j].w));
  //*fprime = coefs[j].y+rho2diff*(2.0f*coefs[j].z + rho2diff*(3.0f*coefs[j].w));
  //*flambda=1.0f-exp(-1.0f*rho2);
  //*fprime=1.0f*exp(-1.0f*rho2);
  return flambda;

}



extern "C"
__host__ double EV (double T, void * p, int* WLlist) {
//Function for calling the Kernel, then computing the 
//Expectation value from the results of each worldline
  struct Wparams params = *(struct Wparams *) p;
  double EV;
  const int groupsize = 128;
  double rtT = sqrt((double)T);
  cudaError_t errorcode;
  // call to integrate the function func
  if( verbosity >= 5)
  	printf("call to CUDA device\n");
  ExpectValue<<<params.nBlocks, params.nThreads>>>(params.Wsscal_d, params.Wsferm_d,
	params.worldlines,params.xcm, (float)params.F, (float)params.l2, (float)rtT, params.Nl, 
	params.Nppl, params.flcoefs, params.rho2sp, params.fermion);
  errorcode = cudaGetLastError();
  if ( errorcode>0) printf("cuda getLastError EV(): %s\n", cudaGetErrorString(errorcode));
  if (verbosity >= 6)
  	printf("return from CUDA\n");
  //Copy device memory back to host
  errorcode = cudaMemcpy(params.Wsscal_h, params.Wsscal_d,
	params.Nl*sizeof(params.Wsscal_h[0]), cudaMemcpyDeviceToHost);
  if(errorcode > 0) printf("cuda memcpy scal Error EV(): %s\n", cudaGetErrorString(errorcode));
  if(params.fermion == 1)
  {
  	errorcode = cudaMemcpy(params.Wsferm_h, params.Wsferm_d,
		params.Nl*sizeof(params.Wsferm_h[0]), cudaMemcpyDeviceToHost);
	if(errorcode > 0) printf("cuda memcpy ferm Error EV(): %s\n",cudaGetErrorString(errorcode));
  }
  //Compute the expectation value from the Wilson Loop data
  EV = 0.0;
  EVMean(&EV, params.Wsscal_h, params.Wsferm_h, groupsize,WLlist, T, params.fermion);
  //printf("EV=%f\n",EV);
  //result=exp(-m2*T)/(T*T*T)*(EV-1.0+1.0/6.0*TB*TB);
  //*SEout=(double) exp(-m2*T)/(T*T*T)*SE;
  //printf("%f %f \n",T,result);
  return EV;
}


extern "C"
__device__ float bump(const float x)
//Device version of the bump function
{
  //the 0.999 makes no numerical difference compared to 1.0
  //but seems to prevent some unpredictable, unspecified CUDA errors
  if(x*x < 0.999f) return __expf( -1.0f/(1.0f-x*x) );
  else return 0.0f;
}

extern "C"
__device__ float phi(const float x)
//Computes the \Phi function which is defined in the thesis
{
  const float onemx2 = 1.0f-x*x;
  //the 0.999 makes no numerical difference compared to 1.0
  //but seems to prevent some unpredictable, unspecified CUDA errors
  if( x < 0.999f )
  {
	return 1.0f - ( 0.5f/q2 )*( onemx2*__expf(-1.0f/onemx2) + expint(-1.0f/onemx2) );
  }
  else
	return 1.0f;
}

extern "C"
__device__ float chi(const float x, const float n, const float lambda)
//Computes the \Chi function which is defined in the thesis
{
  float ans;
  const float onemx2 = 1.0f-x*x;
  const float x2 = x*x;
  const float x4 = x2*x2;
  if(x <= -0.999f) ans = 0.0f;
  else if(x2 < 0.999f)
  {
	ans = 0.5f*(-onemx2*__expf(-1.0f/onemx2) - expint(-1.0f/onemx2));
	ans += 2.0f*n*tubedist/lambda*
		0.444f/( 1.0f + __expf(-(3.31f*x +
			5.25f*x*x2*sin(x)*cos( -0.907f*x2 - 1.29f*x4*x4 ))/cos(x)) );
	//an older, less precise approximation of the integral
	//ans += 2.0f*n*tubedist/lambda*
	//	(0.218f+0.393f*x*coshf(0.806f*x4-0.696f*x2+0.0902f)/coshf(0.825f*x2-0.0234f*x+0.375f));
  }
  else //for x >= 1.0f
  	ans = 2.0f*n*tubedist*q1/lambda;
  return ans;
}


extern "C"
__device__ float flperi(const float rho2, const float lambda2)
//f_lambda for periodic flux tube profile
{
  float expr, exprp1, exprap1, rho, lambda, flraf, f;
  int N = 2;
  int i;
  rho = sqrtf(rho2);
  lambda = sqrtf(lambda2);

  flraf = floorf(rho/tubedist);
  f = 0.0f;
  if((int)flraf-N > 0)
  {
	f=0.5f*(flraf-(float)N-1.0f)*(flraf-(float)N);
  }
  for(i = (int)flraf-N; i <= (int)flraf+N; i++)
  {
	if(i == 0)
	{
		expr = __expf(-rho/lambda);
		exprp1 = 1.0f+expr;
		f = -rho/tubedist*(expr/exprp1) - lambda/tubedist*log(exprp1);
	}
        else if(i > 0)
	{
		exprap1 = 1.0f+__expf(-(rho-(float)i*tubedist)/lambda);
		f = f+(float)i/exprap1;
	}
  }
  return tubedist/(lambda*log(2.0f))*f + 1.0f;
}



extern "C"
__device__ float ffixflux(const float rho2, const float lambda2)
//f_lambda for the fixflux profile
{
  const float lam = sqrtf(lambda2);
  const float lmlmin = (lam-lmin)/(tubedist-lmin);
  const float aml = (tubedist-lam)/(tubedist-lmin);
  const float n = floorf((sqrt(rho2)+tubedist/2.0f)/tubedist);
  float ans;
  if(rho2 <= tubedist*tubedist/4.0f)
  {
	ans = (1.0f-0.75f*lmlmin)*phi(2.0f*sqrt(rho2/lambda2)) + 3.0f*rho2/(tubedist*tubedist)*lmlmin;

  }
  else
  {
	ans = 1.0f+0.75f*(4.0f*rho2/(tubedist*tubedist)-1.0f)*lmlmin + 3.0f*n*(n-1.0f)*aml; 
	ans += 3.0f*lam/(q1*tubedist)*aml*chi(2.0f*(sqrt(rho2)-n*tubedist)/lam, n, lam);
	//ans = 1.01f;
  }
  return ans;
}


extern "C"
__device__ float flambda(float rho2, float lambda2, float4 *flcoefs, float *rho2sp)
//flambda(rho^2,lambda^2) defines the magnetic vector potential in cylindrical coordinates
{
  float f;

  switch(profile){
	case step: 
	  if(rho2 < lambda2)
		f = rho2/lambda2;
	  else
		f = 1.0f;
	  break;
	case smooth:
  	  f = rho2/(lambda2+rho2);
	  break;
	case quadratic:
	  if(rho2 < lambda2)
	  	f = rho2/lambda2*(2.0f-rho2/lambda2);
	  else
		f = 1.0f;
	  break;
	case gaussian:
	  f = 1.0f-exp(-rho2/lambda2);
	  break;
	case periodic:
	  f = flperi(rho2, lambda2);
	  break;
	case spline:
	  f = interp(rho2, rho2sp, flcoefs);
	  break;
	case fixflux:
	  f = ffixflux(rho2, lambda2);
	  break;
  }
  return f;
  //return rho2/lambda2;
} 


extern "C"
__device__ float fplperi(float rho2, float lambda2)
//f'_lambda for periodic profile
{
  float expr,exprp1,rho,lambda,flraf,f;
  float expral, expralp1;
  int N = 2;
  int i;
  rho = sqrtf(rho2);
  lambda = sqrtf(lambda2); 
  flraf = floorf(rho/tubedist);
  f=0.0f;
  for(i = (int)flraf-N; i <= (int)flraf+N; i++)
  {
	if(i == 0)
	{
		expr = __expf(-rho/lambda);
		exprp1 = 1.0f+expr;
		f = expr/(tubedist*exprp1*exprp1);
	}
	else if(i > 0)
	{
		expral = __expf(-(rho-(float)i*tubedist)/lambda);
		expralp1 = 1.0f+expral;
		f = f +((float)i/rho)*expral/(expralp1*expralp1);
	}
  }
  return tubedist/(2.0f*lambda2*log(2.0f))*f;
}




extern "C"
__device__ float fpfixflux(const float rho2, const float lambda2)
//f'_lambda(rho^2) for the fixflux field profile
{
  const float lam = sqrtf(lambda2);
  const float lmlmin = (lam-lmin)/(tubedist-lmin);
  const float rho = sqrtf(rho2);
  const float n = floorf((rho+tubedist/2.0f)/tubedist);
  const float a2 = tubedist*tubedist;
  const float aml = (tubedist - lam)/(tubedist - lmin);
  float ans;
  if(rho <= tubedist/2.0f)
  {
	ans = 2.0f/(lambda2*q2)*(1.0f-0.75f*lmlmin)*bump(2.0f*rho/lam);
	ans += 3.0f/a2*lmlmin;
  }
  else
  {
	ans = 3.0f/a2*lmlmin + 6.0f/(q1*lam*tubedist)*aml*bump(2.0f*(rho-n*tubedist)/lam);
	//ans = 3.0f/(tubedist*tubedist)*lmlmin + 6.0f/(q1*lam*tubedist)*
	//	(tubedist-lam)/(tubedist-lmin)*bump(2.0f*(rho-n*tubedist)/lam);
  }
  //ans = 1.0f;
  return ans;
}


extern "C"
__device__ float fplambda(const float rho2, const float lambda2, float4 *flcoefs, float *rho2sp)
//fplambda(rho^2,lambda^2) defines the magnetic vector potential derivative 
// wrt rho^2 in cylindrical coordinates
{
  float f, fjunk;
  f=1.0f;
  switch(profile){
	case step: 
	  if(rho2 < lambda2)
		f = 1.0f/lambda2;
	  else
		f = 0.0f;
	  break;
	case smooth:
  	  f = lambda2/((lambda2+rho2)*(lambda2+rho2));
	  break;
	case quadratic:
	  if(rho2<lambda2)
	  	f = 2.0f/lambda2*(1.0f-rho2/lambda2);
	  else
		f = 0.0f;
	  break;
	case gaussian:
	  f = 1.0f/lambda2*exp(-1.0f*rho2/lambda2);
	  break;
	case periodic:
	  f = fplperi(rho2, lambda2);
	  break;
	case spline:
	  fjunk = __sinf(sqrtf(rho2)*pi);
	  if(fjunk>100.0f)
		f = 0.0f;
	  else
	  {
		f = __expf(-1.0f*fjunk*fjunk/lambda2)/lambda2;
		//f=1.0f/lambda2*exp(-1.0f*rho2/lambda2);
	  }
	  break;
	case fixflux:
	  f = fpfixflux(rho2, lambda2);
	  break;
  }
  return f;
  //return 1.0f/lambda2;
  //return lambda2/((lambda2+rho2)*(lambda2+rho2));
} 

extern "C"
__device__ void Idt(float *scalI, float *fermI, float4 Ai, 
	const float l2, float4 *flcoefs, float *rho2sp, int fermion)
//Computes the integral over t from 0 to 1 in the scalar and fermion Wilson loop factors
{
	int i;
	const int n = 50;        //number of points in point-to-point proper time integral
	float t,  rhoi2; //proper time and rho squared
	const float h = 1.0f/((float) n);  //distance between points in integral
	float4 xiscal, xiferm;     //scalar and fermi integrands
	if (Ai.x<1.0e-8) Ai.x = 1.0e-8;
	if (Ai.y<1.0e-8) Ai.y = 1.0e-8;
	if (Ai.z<1.0e-8) Ai.z = 1.0e-8;
	float Aip1 = Ai.x+2.0f*Ai.y+Ai.z;  //rho^2 for the final point
	if(Aip1<1.0e-8) Aip1 = 1.0e-8;
	//if(profile == periodic && Aip1 > 10.0f*tubedist) Aip1 = 1.0e-8;
	//Begin the Simpson's method algorithm
	xiscal.x = flambda(Ai.x,l2,flcoefs,rho2sp)/Ai.x + flambda(Aip1,l2,flcoefs,rho2sp)/Aip1;
	xiscal.y = 0.0f;
	xiscal.z = 0.0f;
	if (fermion == 1)
	{
	  xiferm.x = fplambda(Ai.x, l2, flcoefs, rho2sp) + fplambda(Aip1, l2, flcoefs, rho2sp);
	  xiferm.y = 0.0f;
	  xiferm.z = 0.0f;
	}
	for(i = 1; i < n; i++)
	{
		t = (float)i*h;
		//rho2 at the point
		rhoi2 = Ai.x + 2.0f*Ai.y*t + Ai.z*t*t;
		if(rhoi2 < 1.0e-10) rhoi2 = 1.0e-10;
		//if(profile == periodic && rhoi2 > 10.0f*tubedist) rhoi2 = 1.0e-8;
		if(i%2==0) 
		{
			xiscal.z += flambda(rhoi2, l2, flcoefs, rho2sp)/rhoi2;
			if(fermion == 1)
			  xiferm.z += fplambda(rhoi2, l2, flcoefs, rho2sp);
		}
		else 
		{
			xiscal.y += flambda(rhoi2, l2, flcoefs, rho2sp)/rhoi2;
			if(fermion == 1)
			  xiferm.y += fplambda(rhoi2, l2, flcoefs, rho2sp);
		}
	
	}
	*scalI = (xiscal.x + 2.0f*xiscal.z + 4.0f*xiscal.y)*h/3.0f;
	if(fermion == 1)
	  *fermI = (xiferm.x + 2.0f*xiferm.z + 4.0f*xiferm.y)*h/3.0f;
	//*fermI=1.0f/l2;
}

extern "C"
__device__ void getzp1(float4 *zip1, float4 *worldlines, 
	float rtT, float4 xcm, int i, int inx, int Nppl)
//Function for determining the next point on the 
//worldline loop for each of the sub loops
{
  int inxp1;
  //get the next worldline index for the N/2 group
  if(i%2 == 1){
	if(i == Nppl-1)
	{
		inxp1 = inx*Nppl+1;
	}
	else
	{
		inxp1 = inx*Nppl+i+2;
	}
  }
  //get the next worldline index for the first N/4 group
  else if(i%4 == 0){
	if(i == Nppl-4)
	{
		inxp1 = inx*Nppl;
	}
	else
	{
		inxp1 = inx*Nppl+i+4;
	}
  }
  //get the next worldline index for the second N/4 group
  else if((i-2)%2 == 0){
	if(i == Nppl-2)
	{
		inxp1 = inx*Nppl+2;
	}
	else
	{
		inxp1 = inx*Nppl+i+4;
	}
  }
  //compute the next point
  zip1->x = xcm.x + rtT*worldlines[inxp1].x;
  zip1->y = xcm.y + rtT*worldlines[inxp1].y;
  zip1->z = xcm.z + rtT*worldlines[inxp1].z;

}

extern "C"
__device__ void WilsonLoop(float4 *worldlines, float4 *Wsscal, float4 *Wsferm, 
	float4 xcm, int inx, float F, float l2, float rtT, 
	int Nppl, float4 *flcoefs, float *rho2sp, int fermion)
//Returns the Wilson loop value
{
	int i;
	//const float e = 1.0;
        float4 WLstemp, WLftemp;
	float4 zi, zip1;
	float4 Ai;
	float xyyx;
	float scalI, fermI;
	//Compute the scalar contribution
	WLstemp.x = 0.0f; WLstemp.y = 0.0f; WLstemp.z = 0.0f;
	WLftemp.x = 0.0f; WLftemp.y = 0.0f; WLftemp.z = 0.0f;
	for(i = 0; i < Nppl; i++){
		//Compute the scaled, shifted coordinate
		zi.x = xcm.x + rtT*worldlines[inx*Nppl+i].x;
		zi.y = xcm.y + rtT*worldlines[inx*Nppl+i].y;
		getzp1(&zip1, worldlines, rtT, xcm, i, inx, Nppl);
		//Ai Bi and Ci coefficients for the rho2 polynomial
		Ai.x = zi.x*zi.x + zi.y*zi.y;
		Ai.y = zi.x*(zip1.x-zi.x)+zi.y*(zip1.y-zi.y);
		Ai.z = (zip1.x-zi.x)*(zip1.x-zi.x) 
			+ (zip1.y-zi.y)*(zip1.y-zi.y);
		Idt(&scalI, &fermI, Ai, l2, flcoefs, rho2sp, fermion);
		//scalI=1.0f/l2;
		//Compute the contribution to the N/2 integral
		xyyx = (zi.x*zip1.y-zi.y*zip1.x);
		if(i%2 == 1){
			WLstemp.x += xyyx*scalI;
			WLftemp.x += fermI;
		}
		//Compute the contribution to the first N/4 integral
		else if(i%4 == 0){
			WLstemp.z += xyyx*scalI;
			WLftemp.z += fermI;
		}
		//Compute the contribution to the second N/4 integral
		else if((i-2)%2 == 0){
			WLstemp.y += xyyx*scalI;
			WLftemp.y += fermI;
		}
	}
	Wsscal[inx].x = F*WLstemp.x;
	Wsscal[inx].y = F*WLstemp.y;
	Wsscal[inx].z = F*WLstemp.z;
	if( fermion == 1)
	{
	  Wsferm[inx].x = 2.0f*F*WLftemp.x*rtT*rtT/(Nppl/2.0f);
	  Wsferm[inx].y = 2.0f*F*WLftemp.y*rtT*rtT/(Nppl/4.0f);
	  Wsferm[inx].z = 2.0f*F*WLftemp.z*rtT*rtT/(Nppl/4.0f);
	}
	//Wsferm[inx].x=2.0f*F/l2*rtT*rtT;
	//Wsferm[inx].y=2.0f*F/l2*rtT*rtT;
	//Wsferm[inx].z=2.0f*F/l2*rtT*rtT;

	//Wsferm[inx].x=1.0f;
	//Wsferm[inx].y=1.0f;
	//Wsferm[inx].z=1.0f;

	
}

__global__ void 
__launch_bounds__(MY_KERNEL_MAX_THREADS, MY_KERNEL_MIN_BLOCKS)
ExpectValue(float4 *Wsscal, float4 *Wsferm, float4 *worldlines, 
	float4 xcm, float F, float l2, float rtT, int Nl, int Nppl, float4 *flcoefs, float *rho2sp, int fermion)
//Each thread computes the Wilson loop value for a single 
//worldline
{
        int inx = blockIdx.x * blockDim.x + threadIdx.x;       
        WilsonLoop(worldlines, Wsscal, Wsferm, xcm, inx, F, l2, rtT, Nppl, flcoefs, rho2sp, fermion);     
}




