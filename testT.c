#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <builtin_types.h> 
#include "intT.h" 

int main (int argc, char *argv[])
{
  //Number of worldlines, points per line, and dimensions in the worldline file
  const int Nl=5120, Nppl=1024, Dim=2;
  const int fermion=1;
  int i, nBlocks, nThreads=512;
  char *filename="worldlines.5120.1024.2.dat";
  //Arrays to hold values of the Wilson loops on the device and host
  float4 *Wsscal_d, *Wsferm_d;
  float4 *Wsscal_h, *Wsferm_h;
  float4 xcm;
  double exactX;
  //double F=2.0*0.5*100.0, l2=100.0;
  double F=0.5*100.0, l2=100.00;
  double texe, T,X;
  float4 *worldlines_h;
  float4 *worldlines_d;
  cudaError_t errorcode;
  double result, resultSE;
  double* intout;

  xcm.x=0.0;

  if(Nl%nThreads==0)
	nBlocks=Nl/nThreads; 
  else 
	nBlocks=Nl/nThreads+1;

  //parameters for the integration function
  struct Wparams params;

  //allocate memory for worldlines on host and device
  worldlines_h=(float4*)malloc(Nppl*nBlocks*nThreads*sizeof(*worldlines_h));
  cudaMalloc((void **)&worldlines_d,Nppl*nBlocks*nThreads*sizeof(*worldlines_d));
  printf("allocating memory \n");
  if(worldlines_h==NULL | worldlines_d==NULL)
  {
	fprintf(stderr,"out of memory 1\n");
	return(1);
  }
  //read worldlines in from file
  getwl(worldlines_h, filename,Nl,Nppl,Dim);
  printf("Copying worldlines to GPU device \n");
  //Copy worldlines to device
  errorcode=cudaMemcpy(worldlines_d,worldlines_h,
		nThreads*nBlocks*Nppl*sizeof(*worldlines_h),cudaMemcpyHostToDevice);
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
  params.Nl=Nl;
  params.Nppl=Nppl;
  params.F=F; 
  params.l2=l2;
  params.xcm = xcm;
  params.Wsscal_h=Wsscal_h;
  params.Wsscal_d=Wsscal_d;
  params.Wsferm_h=Wsferm_h;
  params.Wsferm_d=Wsferm_d;
  params.worldlines=worldlines_d; 
  params.SE=0.0;
  params.nBlocks=nBlocks;
  params.nThreads=nThreads;
  params.ng=40;
  params.fermion=fermion;
  //intout=(double *)malloc(params.ng*sizeof(intout[0]));
  //if(!intout) printf("main(): malloc failed\n");

  intout=integrateT(20,&params);
  
  result=0.0;
  for(i=0;i<params.ng;i++)
  {
	result+=intout[i];
	printf("rho: intout[%d]=%f\n",i,intout[i]);	
  }
  result/=params.ng;
  resultSE=0.0;
  for(i=0;i<params.ng;i++)
  {
	resultSE+=(intout[i]-result)*(intout[i]-result);
  }
  resultSE=sqrt(resultSE/((double)params.ng*((double)params.ng-1.0)));
  printf("rho integral: %f +/- %f\n",result,resultSE);

  //Free the device memory
  cudaFree(Wsscal_d);
  cudaFree(Wsferm_d);
  cudaFree(worldlines_d);
  //Free the host memory
  free(Wsscal_h);
  free(Wsferm_h);
  free(worldlines_h);
  free(intout);
  
     
  return 0;
}
