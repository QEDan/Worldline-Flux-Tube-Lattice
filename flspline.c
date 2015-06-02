float4* flspline(double* fli, int N)
{
  int i;
  double h;
  double* alpha;
  float4* flscoefs;
  alpha = (double*)malloc(N*sizeof(*alpha));
  flscoefs=(float4 *)malloc(N*sizeof(*flscoefs));

  h=1.0/(double)N;

  alpha[0]=3.0*(fli[1]-fli[0])/h;
  alpha[N-1]=-3.0*(fli[N-1]-fli[N-2])/h;

  for(i=1;i++;i<N-1)
  {
	alpha[i]=3.0/h*(fli[i+1]-2.0*fli[i]+fli[i-1]);
  }

  

}
