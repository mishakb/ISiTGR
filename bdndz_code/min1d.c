/* 1D minimization of a quartic polynomial using 5 evaluations */
#if 0
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#endif

double minquartic(double *y, double xleft, double xright) {
  double alpha, beta;
  double a,b,c,d,e;
  double d0, d1, d2, d3;
  double x, yprime;
  int i;

  /* Auxiliary parameters */
  alpha = (y[1]+y[3])/2. - y[2];
  beta  = (y[0]+y[4])/2. - y[2];

  /* Quartic coefficients */
  a = y[2];
  b = (8*(y[3]-y[1])-y[4]+y[0])/12.;
  c = (16.*alpha-beta)/12.;
  d = (y[3]-y[1])/2.-b;
  e = (beta-4.*alpha)/12.;

  /* Derivatives coefficients */
  d0 = b;
  d1 = c*2;
  d2 = d*3;
  d3 = e*4;

  for(i=0; i<53; i++) {
    x = 0.5*(xleft+xright);
    yprime = d0 + x*(d1 + x*(d2 + x*d3));
    if (yprime>0) {
      xright=x;
    } else {
      xleft=x;
    }
  }

  return(x);
}

#if 0
int main(void) {
  double xmin;
  double y[]={81,16,1,0,1};

  xmin=minquartic(y, -5., 5.);
  printf("xmin=%9.6lf\n", xmin);

  return(0);
}
#endif
