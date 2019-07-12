/* ====================================================== *
 * MK 1/2005						  *
 * Nonlinear density (3-d) and convergence (2-d) power    *
 * spectra for LambdaCDM models.			  *
 * Based on Peacock&Dodds (1996) and Smith et al. (2002). *
 * The latter nonlinear model incorporates an improved    *
 * and faster version of halofit.f.			  *
 * ====================================================== */

#include <math.h>
#include <stdlib.h>
#include <malloc.h>
#include <stdio.h>
#include <assert.h>
#include "smith2.h"

typedef real fr();    /* function of real */

extern real Omega_m, Gamma, sigma_8, n_spec, z0, beta_p;
real s1glob, s2glob;


/* ============================================================ *
 * Functions available for calls from fortran.			*
 * ============================================================ */

void pl_(real *a, real *k, real *res)
{
   *res = P_L(*a, *k);
}

void pnl_(real *a, real *k, real *res)
{
   *res = P_NL(*a, *k);
}

void pkappa_(real *s, real *res)
{
   *res = Pkappa(*s);
}

void setparameters_(real *OMEGAM, real *OMEGAV, real *GAMMA, real *SIGMA8,
		   real *NSPEC, real *BETAP, real *Z0, int *NONLINEAR)
{
   Omega_m   = *OMEGAM;
   Omega_v   = *OMEGAV;
   Gamma     = *GAMMA;
   sigma_8   = *SIGMA8;
   n_spec    = *NSPEC;
   beta_p    = *BETAP;
   z0        = *Z0;
   nonlinear = *NONLINEAR;
}

void setparameters(real *OMEGAM, real *OMEGAV, real *GAMMA, real *SIGMA8,
		   real *NSPEC, real *BETAP, real *Z0, int *NONLINEAR)
{
}

void set_cosmological_parameters_to_default()
{
   Omega_m   = 0.3;
   Omega_v   = 0.7;
   Gamma     = 0.21;
   sigma_8   = 0.9;
   n_spec    = 1.0;
   beta_p    = 1.5;
   z0        = 1.0;
   nonlinear = 2;
}

void dump_param(FILE *F)
{
   if (!F) F = stderr;
   fprintf(F,
	   "#Om Ov G s8 n_s beta_p z0 nonlinear = %.4f %.4f %.4f %.4f %.4f %.4f %.4f %d\n",
	   Omega_m, Omega_v, Gamma, sigma_8, n_spec, beta_p, z0, nonlinear);
}

void error(const char *s)
{
   (void)fprintf(stderr, "error: ");
   (void)fprintf(stderr,s);
   (void)fprintf(stderr,"\n");
   exit(1);
}

real **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a real matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	real **m;

	/* allocate pointers to rows */
	m=(real **) malloc((size_t)((nrow+NR_END)*sizeof(real*)));
	if (!m) error("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(real *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(real)));
	if (!m[nrl]) error("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_matrix(real **m, long nrl, long nrh, long ncl, long nch)
/* free a real matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

/* ============================================================ *
 * dfridr.c							*
 * NR page 188. Returns derivate of func at x, initial step is  *
 * h. Error estimate in err.					*
 * Modified! func depends on two real! (like P_L)		*
 * ============================================================ */

#define CON 1.4
#define CON2 (CON*CON)
#define BIG 1.0e30
#define NTAB 10
#define SAFE 2.0

real *vector(long nl, long nh)
/* allocate a real vector with subscript range v[nl..nh] */
{
	real *v;

	v=(real *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(real)));
	if (!v) error("allocation failure in vector()");
	return v-nl+NR_END;
}

void free_vector(real *v, long nl, long nh)
/* free a real vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

real dfridr(real (*func)(real,real), real x, real h, real *err, real aa)
{
	int i,j;
	real errt,fac,hh,**a,ans;

	ans = 1e30; /* dummy initialization */
	if (h == 0.0) error("h must be nonzero in dfridr.");
	a=matrix(1,NTAB,1,NTAB);
	hh=h;
	a[1][1]=((*func)(aa,x+hh)-(*func)(aa,x-hh))/(2.0*hh);
	*err=BIG;
	for (i=2;i<=NTAB;i++) {
		hh /= CON;
		a[1][i]=((*func)(aa,x+hh)-(*func)(aa,x-hh))/(2.0*hh);
		fac=CON2;
		for (j=2;j<=i;j++) {
			a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
			fac=CON2*fac;
			errt=FMAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
			if (errt <= *err) {
				*err=errt;
				ans=a[j][i];
			}
		}
		if (fabs(a[i][i]-a[i-1][i-1]) >= SAFE*(*err)) break;
	}
	free_matrix(a,1,NTAB,1,NTAB);
	return ans;
}
#undef CON
#undef CON2
#undef BIG
#undef NTAB
#undef SAFE

void polint(real xa[], real ya[], int n, real x, real *y, real *dy)
{
	int i,m,ns=1;
	real den,dif,dift,ho,hp,w;
	real *c,*d;

	dif=fabs(x-xa[1]);
	c=vector(1,n);
	d=vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0)
			  error("Error in routine polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free_vector(d,1,n);
	free_vector(c,1,n);
}


void spline(real x[], real y[], int n, real yp1, real ypn, real y2[])
{
	int i,k;
	real p,qn,sig,un,*u;

	u=vector(1,n-1);
	if (yp1 > 0.99e30)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	for (i=2;i<=n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
	for (k=n-1;k>=1;k--)
		y2[k]=y2[k]*y2[k+1]+u[k];
	free_vector(u,1,n-1);
}

void splint(real xa[], real ya[], real y2a[], int n, real x, real *y)
{
	int klo,khi,k;
	real h,b,a;

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) {
	   error("Bad xa input to routine splint");
	}
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
	*y=a*ya[klo]+b*ya[khi]+
	  ((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
}

#define FUNC(x) ((*func)(x))

real trapzd(real (*func)(real), real a, real b, int n)
{
	real x,tnm,sum,del;
	static real s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) {
		   sum += FUNC(x);
		}
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC

#define FUNC(x) ((*func)(x))

real trapzd1(real (*func)(real), real a, real b, int n)
{
	real x,tnm,sum,del;
	static real s;
	int it,j;

	if (n == 1) {
		return (s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) {
		   sum += FUNC(x);
		}
		s=0.5*(s+(b-a)*sum/tnm);
		return s;
	}
}
#undef FUNC

/* ============================================================ *
 * Interpolates f at the value x, where f is a real[n] array,	*
 * representing a function between a and b, stepwidth dx.	*
 * 'lower' and 'upper' are powers of a logarithmic power law	*
 * extrapolation. If no	extrapolation desired, set these to 0	*
 * ============================================================ */
real interpol(real *f, int n, real a, real b, real dx, real x,
	      real lower, real upper)
{
   real r;
   int  i;
   if (x < a) {
      if (lower==0.) {
	 error("value too small in interpol");
	 return 0.0;
      }
      return f[0] + lower*(x - a);
   }
   r = (x - a)/dx;
   i = (int)(floor(r));
   if (i+1 >= n) {
      if (upper==0.0) {
	 if (i+1==n) {
	    return f[i];  /* constant extrapolation */
	 } else {
	    error("value too big in interpol");
	    return 0.0;
	 }
      } else {
	 return f[n-1] + upper*(x-b); /* linear extrapolation */
      }
   } else {
      return (r - i)*(f[i+1] - f[i]) + f[i]; /* interpolation */
   }
}


/* ============================================================ *
 * like interpol, but f beeing a 2d-function			*
 * 'lower' and 'upper' are the powers of a power law extra-	*
 * polation in the first argument				*
 * ============================================================ */
real interpol2d(real **f,
		int nx, real ax, real bx, real dx, real x,
		int ny, real ay, real by, real dy, real y,
		real lower, real upper)
{
   real t, dt, s, ds;
   int i, j;
   if (x < ax) {
      error("value too small in interpol2d");
   }
   if (x > bx) {
      error("value too big in interpol2d");
   }
   t = (x - ax)/dx;
   i = (int)(floor(t));
   if (i+1 > nx || i < 0) error("index out of range in interpol");
   dt = t - i;
   if (y < ay) {
      return ((1.-dt)*f[i][0] + dt*f[i+1][0]) + (y-ay)*lower;
   } else if (y > by) {
      return ((1.-dt)*f[i][ny-1] + dt*f[i+1][ny-1]) + (y-by)*upper;
   }
   s = (y - ay)/dy;
   j = (int)(floor(s));
   ds = s - j;
   return (1.-dt)*(1.-ds)*f[i][j] +
     (1.-dt)*ds*f[i][j+1] +
     dt*(1.-ds)*f[i+1][j] +
     dt*ds*f[i+1][j+1];
}

real gammln(real xx)
{
	double x,y,tmp,ser;
	static double cof[6]={76.18009172947146,-86.50532032941677,
		24.01409824083091,-1.231739572450155,
		0.1208650973866179e-2,-0.5395239384953e-5};
	int j;

	y=x=xx;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0;j<=5;j++) ser += cof[j]/++y;
	return (real)(-tmp+log(2.5066282746310005*ser/x));
}



/* ============================================================ *
 * qromb.c							*
 * Romberg Integration. Uses trapzd. NR p. 140			*
 * ============================================================ */

#define EPS 1.0e-6
#define JMAX 35
#define JMAXP (JMAX+1)
#define K 5

real qromb(real (*func)(real), real a, real b)
{
	real ss,dss;
	real s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=0.25*h[j];
	}
	error("Too many steps in routine qromb");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K

#define EPS 1.0e-6
#define JMAX 40
#define JMAXP (JMAX+1)
#define K 5

real qromb1(real (*func)(real), real a, real b)
{
	real ss,dss;
	real s[JMAXP],h[JMAXP+1];
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=trapzd1(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=0.25*h[j];
	}
	error("Too many steps in routine qromb");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K


#define EPS 1.0e-7
#define JMAX 35
#define JMAXP (JMAX+1)
#define K 5

/* ============================================================ *
 * Romberg integration of an open intervall [a,b]. choose is a  *
 * pointer to a routine using an open quadrature formula.	*
 * Uses polint (Neville-Aitken) to extrapolate. NR p. 143	*
 * ============================================================ */

real qromo(real (*func)(real), real a, real b,
	real (*choose)(real(*)(real), real, real, int))
{
	int j;
	real ss,dss,h[JMAXP+1],s[JMAXP];

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=(*choose)(func,a,b,j);
		if (j >= K) {
			polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=h[j]/9.0;
	}
	error("Too many steps in routing qromo");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K

#define FUNC(x) ((*func)(x))

real midpnt(real (*func)(real), real a, real b, int n)
{
	real x,tnm,sum,del,ddel;
	static real s;
	int it,j;

	if (n == 1) {
		return (s=(b-a)*FUNC(0.5*(a+b)));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += FUNC(x);
			x += ddel;
			sum += FUNC(x);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}
#undef FUNC




/* ============================================================ *
 * Global variables needed in integrand functions int_for_p_2	*
 * and int_for_g						*
 * ============================================================ */

real sglob, aglob, rglob;

/* CPT 9 */
real da_dtau(real a)
{
   real res, det;
   if (a < epsilon) error("Division by 0 in da_dtau");
   det = 1. + Omega_m*(1./a - 1.) + Omega_v*(DSQR(a) - 1.);
   res = sqrt(det);
   return res;
}


real da_dtau_m3(real a)
{
   real res;
   if (a < epsilon) return 0.;
   res = da_dtau(a);
   if (res < epsilon) error("Division by 0 in da_dtau_m3");
   return 1./(res*res*res);
}

real D_plus(real a)
{
   static real OMEGA_M = -42.;
   static real OMEGA_V = -42.;
   static real da = 0.;
   static real table[N_a];
   real delta, aa, delta0;
   int  i;

   if (fabs(Omega_m-OMEGA_M)>epsilon || fabs(Omega_v-OMEGA_V)>epsilon) {
      delta0 = qromb(da_dtau_m3, 0., 1.);
      da = (1. - a_min)/(N_a-1.);
      aa = a_min;
      for (i = 0; i<N_a; i++, aa += da) {
	 delta = qromb(da_dtau_m3, 0., aa);
	 table[i] = 1./aa * da_dtau(aa) * delta / delta0;
      }
      OMEGA_M = Omega_m;
      OMEGA_V = Omega_v;
   }
   return interpol(table, N_a, a_min, 1., da, a, .0, .0);
}

/* PD 96 - (15,16) */
void Omega_a(real a, real *omega_m, real *omega_v)
{
   real f, a3;
   a3 = a*a*a;
   f = a + Omega_m*(1.-a) + Omega_v*(a3-a);
   *omega_m = Omega_m/f;
   *omega_v = Omega_v*a3/f;
}

/* CPT 29, PD 15+16 */
real g(real a)
{
   real omega_m, omega_v;
   Omega_a(a, &omega_m, &omega_v);

   return 2.5*omega_m/(pow(omega_m,4./7.) - omega_v + (1. + .5*omega_m)*(1. + omega_v/70.));
}

/* BBKS G3, see also PD2 51, 52. Returns k^n * T^2(k) */
real Tsqr(real k)
{
   real q, f1, f2;
   q = k/(2998.*Gamma);
   f1 = log(1 + 2.34*q)/(2.34*q);
   f2 = 1 + q*(3.89 + q*(259.21 + q*(162.771336 + q*2027.16958081)));
   /* polynomial: 1 + 3.89q + (16.1q)^2 + (5.46q)^3 + (6.71q)^4 */

   assert(finite(f1)&&finite(f2));

   return pow(k, n_spec)*DSQR(f1)/sqrt(f2);
}

real int_for_sigma_8(real k)
{
   real kR, res, x;
   kR = k/375.;
   x = (sin(kR) - kR*cos(kR))/(kR*kR*kR);
   res = DSQR(k)*Tsqr(k)*x*x;
   return res;
}

/* PD2 42, for a=1, so D+=1 */
real sigma_8_sqr()
{
   static real N_SPEC = -42.;
   static real GAMMA  = -42.;
   static real res = -42.;
   real integral;

   if (fabs(N_SPEC-n_spec)>epsilon || fabs(GAMMA-Gamma)>epsilon) {
      integral = qromb(int_for_sigma_8, k_min, k_max);
      res = 4.5/pi_sqr*integral;
      N_SPEC = n_spec;
      GAMMA = Gamma;
   }
   assert(res>0.0);

   return res;
}

/* lens efficiency */
real G(real a)
{
   return 1.5*Omega_m/a*g_source(a);
}

real P_L(real a, real k)
{
   real d, t, s;
   d = D_plus(a);
   t = Tsqr(k);
   s = sigma_8_sqr();
   return DSQR(sigma_8)*DSQR(d)*t/s;
}


/* PD 22 */
real n_L(real a, real k)
{
   real diff, hh, err, n;

   hh   = k/20.;
   diff = dfridr(P_L, 0.5*k, hh, &err, a);
   n    = .5*k/P_L(a,.5*k)*diff;
   return n;
}

/* PD 21, 23-27 */
real f_NL(real x, real a, real k)
{
   real A, B, alpha, beta, V;
   real c, gg, f0, f1, f2, f3, f4;

   c = 1.+n_L(a,k)/3.;
   if (c<0) error("spectral index too small for PD fitting formula");
   A = .482/pow(c,.947);
   B = .226/pow(c,1.778);
   alpha = 3.31/pow(c,.244);
   beta = .862/pow(c,.287);
   V = 11.55/pow(c,.423);
   gg = g(a);

   f0 = pow(A*x, alpha);
   f1 = 1. + B*beta*x + pow(f0, beta);
   f2 = f0*gg*gg*gg/(V*sqrt(x));
   f3 = 1. + pow(f2, beta);
   f4 = x*pow(f1/f3, 1./beta);

   assert(finite(f4));
   return f4;
}

/* ============================================================ *
 * Bond&Efstathiou 1984 approximation to Delta_L		*
 * ============================================================ */

real Delta_L_BE2(real k)
{
#define N_kk (500)
   const real kmin = 1.e-5;     /* in units of H_0/c */
   const real kmax = 6000.;     /* dt., from wint    */
   static real GAMMA = -42.;
   static real SIGMA = -42.;
   static real NSPEC = -42.;
   static real keff, q8, table[N_kk], logkmin, logkmax, dlogk, tk8;
   real q, tk;
   real klog, kk, res;
   int j;

   if (GAMMA!=Gamma || SIGMA!=sigma_8 || NSPEC!=n_spec) {

      if (GAMMA!=Gamma) {
	 keff = 0.172 + 0.011*DSQR(log(Gamma/0.36));
	 q8   = 1e-20 + keff/Gamma;
	 tk8  = 1./pow(1+pow(6.4*q8+pow(3.0*q8,1.5)+(1.7*q8)*(1.7*q8),1.13),(1/1.13));
      }

      logkmin = log(kmin);
      logkmax = log(kmax);
      dlogk   = (logkmax - logkmin)/(N_kk-1.);
      for (j=0,klog=logkmin; j<N_kk; j++,klog+=dlogk) {

	 kk = exp(klog);
	 q  = 1e-20 + kk/Gamma;
	 tk = 1/pow(1+pow(6.4*q+pow(3.0*q,1.5)+(1.7*q)*(1.7*q),1.13),(1/1.13));
	 table[j] = (3+n_spec)*log(q/q8) + 2*log(sigma_8*tk/tk8);

      }

      GAMMA = Gamma;
      SIGMA = sigma_8;
      NSPEC = n_spec;

   }

   kk = log(k);
   if (kk<logkmin || kk>logkmax) {
      q  = 1e-20 + k/Gamma;
      tk = 1/pow(1+pow(6.4*q+pow(3.0*q,1.5)+(1.7*q)*(1.7*q),1.13),(1/1.13));
      res = sigma_8*sigma_8*pow(q/q8,3+n_spec)*tk*tk/tk8/tk8;
   } else {
      res = exp(interpol(table, N_kk, logkmin, logkmax, dlogk, kk, 0.0, 0.0));
   }

   return res;
#undef N_kk
}

real Delta_L_BE(real k)
{
   real keff, q, q8, tk, tk8;
 
   keff=0.172+0.011*log(Gamma/0.36)*log(Gamma/0.36);
   q=1e-20+k/Gamma;
   q8=1e-20+keff/Gamma;
   tk=1/pow(1+pow(6.4*q+pow(3.0*q,1.5)+(1.7*q)*(1.7*q),1.13),(1/1.13));
   tk8=1/pow(1+pow(6.4*q8+pow(3.0*q8,1.5)+(1.7*q8)*(1.7*q8),1.13),(1/1.13));
 
   return sigma_8*sigma_8*pow(q/q8,3+n_spec)*tk*tk/tk8/tk8;
}

/* ============================================================ *
 * Calculates k_NL, n_eff, n_cur				
 * ============================================================ */ 

real int_for_wint2_knl(real logk)
{
   real krsqr, k;

   k = exp(logk);
   krsqr = DSQR(k*rglob);
   return Delta_L_BE2(k)*exp(-krsqr);
}

real int_for_wint2_neff(real logk)
{
   real krsqr, k;

   k = exp(logk);
   krsqr = DSQR(k*rglob);
   return Delta_L_BE2(k)*2.0*krsqr*exp(-krsqr);
}

real int_for_wint2_ncur(real logk)
{
   real krsqr, k;

   k = exp(logk);
   krsqr = DSQR(k*rglob);
   return Delta_L_BE2(k)*4.0*krsqr*(1.0-krsqr)*exp(-krsqr);
}

void wint(real r, real *sig, real *d1, real *d2, real amp)
{
   real sum1, sum2, sum3, t, x, w1, w2, w3, k, xsqr, d2overktsqr;
   int	i, nint;

   nint = 10000.;
   sum1 = sum2 = sum3 = 0.0;

   for(i=1; i<=nint; i++)
   {
      t     = (i-0.5)/(real)nint;
      k     = -1.0 + 1.0/t;
      *d2   = Delta_L_BE2(k);
      x     = k*r;
      xsqr  = x*x;
      w1    = exp(-xsqr);
      w2    = 2.*xsqr*w1;
      w3    = 2.*(1.0-xsqr)*w2;
      d2overktsqr  = *d2/(k*t*t);
      sum1 += w1*d2overktsqr;
      sum2 += w2*d2overktsqr;
      sum3 += w3*d2overktsqr;
   }

   sum1 = sum1/(real)nint;
   sum2 = sum2/(real)nint;
   sum3 = sum3/(real)nint;
   *sig = amp*sqrt(sum1);
   *d1  = -sum2/sum1;
   *d2  = -DSQR(*d1) - sum3/sum1;

}

void wint2(real r, real *sig, real *d1, real *d2, real amp, int onlysig)
{
   const real kmin = 1.e-2;
   real kmax, logkmin, logkmax, s1, s2, s3;

   /* choose upper integration limit to where filter function dropped
    substantially */
   kmax  = sqrt(5.*log(10.))/r;
   if (kmax<8000.0) kmax = 8000.0;

   logkmin = log(kmin);
   logkmax = log(kmax);
   rglob = r;

   if (onlysig==1) {
      s1   = qromb1(int_for_wint2_knl, logkmin, logkmax);
      *sig = amp*sqrt(s1);
   } else s1 = DSQR(1/amp);   /* sigma = 1 */

   if (onlysig==0) {
      s2  = qromb1(int_for_wint2_neff, logkmin, logkmax);
      s3  = qromb1(int_for_wint2_ncur, logkmin, logkmax);
      *d1 = -s2/s1;
      *d2 = -DSQR(*d1) - s3/s1;
   }
}

/* slope in the highly nonlinear regime, c.f. Smith et al (2002) eq. (61) */
real slope_NL(real rn, real rncur, real om_m, real om_v)
{
   real gam, f1a, f1b, frac, f1;

   gam = 0.86485 + 0.2989*rn + 0.1631*rncur;
   if(fabs(1-om_m)>0.01) {
      f1a = pow(om_m,(-0.0732));
      f1b = pow(om_m,(-0.0307));
      frac = om_v/(1.-om_m);  
      f1 = frac*f1b + (1-frac)*f1a;
   } else {
      f1 = 1.0;
   }

   return 3.0*(f1-1.0) + gam - 3.0;
}

void halofit(real rk, real rn, real rncur, real rknl, real plin, 
	     real om_m, real om_v, real *pnl)
{
   real gam,a,b,c,xmu,xnu,alpha,beta,f1,f2,f3;
   real y, ysqr;
   real f1a,f2a,f3a,f1b,f2b,f3b,frac,pq,ph;
   real nsqr;

   nsqr = rn*rn;
   gam = 0.86485 + 0.2989*rn + 0.1631*rncur;
   a = 1.4861 + 1.83693*rn + 1.67618*nsqr + 0.7940*rn*nsqr
     + 0.1670756*nsqr*nsqr - 0.620695*rncur;
   a = pow(10,a);
   b = pow(10,(0.9463+0.9466*rn+0.3084*nsqr-0.940*rncur));
   c = pow(10,(-0.2807+0.6669*rn+0.3214*nsqr-0.0793*rncur));
   xmu = pow(10,(-3.54419+0.19086*rn));
   xnu = pow(10,(0.95897+1.2857*rn));
   alpha = 1.38848+0.3701*rn-0.1452*nsqr;
   beta = 0.8291+0.9854*rn+0.3400*nsqr;

   if(fabs(1-om_m)>0.01) {
      f1a = pow(om_m,(-0.0732));
      f2a = pow(om_m,(-0.1423));
      f3a = pow(om_m,(0.0725));
      f1b = pow(om_m,(-0.0307));
      f2b = pow(om_m,(-0.0585));
      f3b = pow(om_m,(0.0743));
  
      frac = om_v/(1.-om_m);  
      f1 = frac*f1b + (1-frac)*f1a;
      f2 = frac*f2b + (1-frac)*f2a;
      f3 = frac*f3b + (1-frac)*f3a;
   } else {      /* EdS Universe */
      f1 = f2 = f3 = 1.0;
   }

   y = rk/rknl;
   ysqr = y*y;
   ph = a*pow(y,f1*3)/(1+b*pow(y,f2)+pow(f3*c*y,3-gam));
   ph = ph/(1+xmu/y+xnu/ysqr);
   pq = plin*pow(1+plin,beta)/(1+plin*alpha)*exp(-y/4.0-ysqr/8.0);
   *pnl = pq + ph;

   assert(finite(*pnl));
}
 
real dlog(real x)
{
   return log(x)/log(10.0);
}

real P_NL(real a, real k_NL)
{
   const real kNLstern = 1.e6;       /* h/Mpc */

   static real OMEGA_M   = -42.;
   static real OMEGA_V   = -42.;
   static real N_SPEC    = -42.;
   static real GAMMA     = -42.;
   static real SIGMA_8   = -42.;
   static int  NONLINEAR = -42;
   static real upper;
   static real table_k[N_k], table_P[N_k], y2[N_k], table_slope[N_a];
   static real **table_P_NL = 0;
   static real logkmin = 0., logkmax = 0., dk = 0., da = 0.;
   real Delta_NL, Delta_L, k_L, lnk_NL;
   
   real	omm, omv, amp;
   real logr1, logr2 ,diff, rmid, sig, d1, d2;
   real	rknl, rneff, rncur;
      
   real aa, klog, val, logrmidtmp, logrmid, logr1start, logr2start;
   int i,j, iter, golinear;
   const real logstep = 5.0;
   const int itermax  = 20;

   if (OMEGA_M != Omega_m || OMEGA_V != Omega_v || N_SPEC != n_spec ||
       GAMMA != Gamma || SIGMA_8 != sigma_8 || NONLINEAR != nonlinear) {
 
      if (!table_P_NL) table_P_NL = matrix(0, N_a-1, 0, N_k-1);

      /* upper = (dlnP/dlnk)_{k=kmax}, for splines & extrapolation */
      /* Note that the in the range considered here the linear power
       * spectrum is still a bit shallower than k^(n-4) since T(k) has
       * not yet reached its asymptotic limit. 
       */
      if (nonlinear==0) upper = n_spec-4.0;
      else if (nonlinear==1) upper = -2.5;

      da = (1. - a_min)/(N_a-1.);
      aa = a_min;
      logkmin = log(k_min);
      logkmax = log(k_max);
      dk = (logkmax - logkmin)/(N_k-1.);
  
      for (i=0; i<N_a; i++, aa +=da) {
	 /* printf("#i=%d\n", i);  */
	 klog = logkmin;

	 if(nonlinear==2) {           /* Smith et al. (2002) */
	    amp = aa*g(aa)/g(1.0);
	    Omega_a(aa, &omm, &omv);
	    golinear = 0;

	    /* find non-linear scale with iterative bisection */
	    logr1 = -2.0;
	    logr2 =  3.5;

	 iterstart:

	    logr1start = logr1;
	    logr2start = logr2;

	    iter = 0;
	    do {
	       logrmid = (logr2+logr1)/2.0;
	       rmid    = pow(10,logrmid);
	       /* wint(rmid, &sig, &d1, &d2, amp); */
	       wint2(rmid, &sig, 0x0, 0x0, amp, 1);

	       diff = sig - 1.0;

	       /* printf("rmid=%e, sig=%f diff=%f iter=%d\n", rmid, sig, diff, iter); */
	       /* if (i==55) printf("%e %f\n", rmid, sig); */

	       if(diff>0.001)
		 logr1 = dlog(rmid);
	       if(diff<-0.001)
		 logr2 = dlog(rmid);

	    } while (fabs(diff)>=0.001 && ++iter<itermax);

	    if (iter>=itermax) {
	       logrmidtmp = (logr2start+logr1start)/2.0;
	       /* printf("1 m m 2 = %f %f %f %f\n", logr1start, logrmid, logrmidtmp, logr2start); */
	       if (logrmid<logrmidtmp) {
		  logr1 = logr1start-logstep;
		  logr2 = logrmid;
	       } else if (logrmid>=logrmidtmp) {
		  logr1 = logrmid;
		  logr2 = logr2start+logstep;
	       }

	       /* printf("i=%d, logr1=% .e logr2=% .e k2=%e\n", i, logr1, logr2, 1/pow(10, logr2)); */

	       /* non-linear scale far beyond maximum scale: set flag golinear */
	       if (1/pow(10, logr2)>kNLstern) {
		  golinear = 1;
		  upper = table_slope[i] = n_spec-4.0;
		  goto after_wint;
	       } else {
		  goto iterstart;
	       }
	    }

	    /* spectral index & curvature at non-linear scale */
	    wint2(rmid, &sig, &d1, &d2, amp, 0);
	    rknl  = 1./rmid;
	    rneff = -3-d1;
	    rncur = -d2;

	    upper = table_slope[i] = slope_NL(rneff, rncur, omm, omv);

	 }

      after_wint:

	 for (j=0; j<N_k; j++, klog+=dk) {
	    k_L = exp(klog);

	    if (nonlinear==0) {         /* linear power spectrum */

	       Delta_NL = P_L(aa,k_L)*k_L*k_L*k_L/(2.*pi_sqr);
	       lnk_NL   = klog;

	    } else if (nonlinear==1) {                 /* PD 1996 */

	       Delta_L  = P_L(aa,k_L)*k_L*k_L*k_L/(2.*pi_sqr);
	       Delta_NL = f_NL(Delta_L, aa, k_L);
	       lnk_NL   = klog + 1./3.*log(1 + Delta_NL);  /* PD (5) */

	    } else if (nonlinear==2) {  /* Smith et al. 2002 */

	       Delta_L = amp*amp*Delta_L_BE2(k_L/2998.);
	       if (golinear==0) {
		  halofit(k_L/2998.,rneff,rncur,rknl,Delta_L,omm,omv,&Delta_NL);
	       } else {
		  Delta_NL = Delta_L;
	       }
	       lnk_NL = klog;

	    }

	    table_k[j] = lnk_NL;
	    table_P[j] = log(2*pi_sqr*Delta_NL) - 3.*lnk_NL; /* PD (3) */
	 }

	 spline(table_k-1, table_P-1, N_k, n_spec, upper, y2-1);
	 klog = logkmin;
	 for (j=0; j<N_k; j++, klog += dk) {
	    splint(table_k-1, table_P-1, y2-1, N_k, klog, &val);
	    table_P_NL[i][j] = val;
	 }
      }

      OMEGA_M = Omega_m;
      OMEGA_V = Omega_v;
      N_SPEC  = n_spec;
      GAMMA   = Gamma;
      SIGMA_8 = sigma_8;
      NONLINEAR = nonlinear;
      assert(nonlinear==0 || nonlinear==1 || nonlinear==2);
   }

   klog = log(k_NL);

   if (nonlinear==2) {
      upper = interpol(table_slope, N_a, a_min, 1.0, da, a, 1e31, 1e31);
   }

   val = interpol2d(table_P_NL, N_a, a_min, 1., da, a, N_k, logkmin, logkmax, 
		    dk, klog, n_spec, upper);
   return exp(val);
}

real int_for_w(real a)
{
   real asqr, d;
   asqr = DSQR(a);
   d = a*Omega_m + asqr*(1.-Omega_m-Omega_v) + DSQR(asqr)*Omega_v;

   if (d<0.0) error("H^2 is negative. You are trying to probe a closed Universe at larger than maximum redshift");

   return 1.0/sqrt(d);
}

/* BS01 2.41, with a(z_2) = 0, a(z_1) = z, in units of Hubble radius c/H_0 = 2998/h Mpc */
real w(real a)
{
   static real OMEGA_M = -42.;
   static real OMEGA_V = -42.;
   static real table[N_a];
   static real da = 0.0;
   real   aa;
   int    i;

   if (OMEGA_M != Omega_m || OMEGA_V != Omega_v) {
      da = (1.-a_min)/(N_a-1.);
      aa = a_min;
      for (i=0; i<N_a-1; i++, aa+=da) {
	 table[i] = qromb(int_for_w, aa, 1.);
      }
      table[N_a-1] = .0;
      OMEGA_M = Omega_m;
      OMEGA_V = Omega_v;
   }
   return interpol(table, N_a, a_min, 1., da, a, .0, .0);
}

/* BS01 2.4, 2.30 */
real f_K(real w)
{
   real K, K_h, f;

   K = Omega_m + Omega_v - 1;
   if (K > epsilon) {           /* open */
      K_h = sqrt(K);
      f = 1./K_h*sin(K_h*w);
   } else if (K < -epsilon) {   /* closed */
      K_h = sqrt(-K);
      f = 1./K_h*sinh(K_h*w);
   } else {                     /* flat */
      f = w;
   }

   if (f<0.0) error("angular distance is negative. You are trying to probe a closed Universe at larger than maximum redshift");

   return f;
}

/* S98 2.11. Note that a single redshift sheet (beta_p=0) is implemented
 * only in g_source. */
real prob(real z)
{
   static real norm = 0.;
   static real BETA_P = -42.;
   static real Z0     = -42.;
   real x, f;

   if (fabs(beta_p)<epsilon) error("beta_p=0 in prob");

   if (BETA_P != beta_p || Z0 != z0) {
      norm = beta_p/(z0*exp(gammln(3./beta_p)));
      BETA_P = beta_p;
      Z0     = z0;
   }
   if (z<0 && z<-epsilon) error("negative z in prob");
   if (z<0) x = 0;
   else x = z/z0;
   if (fabs(beta_p - 1.) < epsilon) {
      f = x;
   } else {
      f = pow(x,beta_p);
   }
   f=exp(-f);
   return norm*DSQR(x)*f;
}

/* ! global variable 'aglob' is needed here ! */
real int_for_g(real aprime)
{
   real ww, wprime;
   ww = w(aglob);
   wprime = w(aprime);
   return prob(1./aprime-1.)*f_K(wprime-ww)/f_K(wprime)/DSQR(aprime);
}

/* S98 2.9 */
real g_source(real a)
{
   static real OMEGA_M = -42.;
   static real OMEGA_V = -42.;
   static real BETA_P  = -42.;
   static real Z0      = -42.;
   static real table[N_a];
   static real da = 0.0;
   real   aa, wdelta;
   int    i;

   /* single source redshift at z0, p(w) = delta(w-w0) */
   if (beta_p == 0.0) {
      aa = 1./(1.+z0);
      if (a<=aa) return 0.0;
      wdelta = w(aa);
      return f_K(wdelta-w(a))/f_K(wdelta);
   }
   if (OMEGA_M != Omega_m || OMEGA_V != Omega_v || BETA_P != beta_p || Z0 != z0) {
      da = (1.-a_min)/(N_a-1.);
      table[0] = 0.0;
      aa = a_min+da;
      for (i=1;i<N_a-1;i++,aa+=da) {
	 aglob = aa;
	 table[i] = qromb(int_for_g, a_min, aa);
      }
      table[N_a-1] = 1.;
      OMEGA_M = Omega_m;
      OMEGA_V = Omega_v;
      BETA_P  = beta_p;
      Z0      = z0;
   }
   return interpol(table, N_a, a_min, 1., da, a, .0, .0);
}

/* ! global variable sglob is needed here ! */
real int_for_p_2(real a)
{
   real hoverh0, asqr, s, fKw, f, res, wa;

   if (a >= 1.0) error("a>=1 in int_for_p_2");
   s       = sglob;
   asqr    = DSQR(a);
   wa      = w(a);
   fKw     = f_K(wa);
   assert(fKw>0);
   f       = s/fKw;

   hoverh0 = sqrt(Omega_m/(a*asqr) + (1.-Omega_m-Omega_v)/asqr + Omega_v);
   res = DSQR(g_source(a)/asqr)/hoverh0;
   if (!finite(res)) assert(0);

   if (nonlinear==1 || nonlinear==2) res *= P_NL(a, f);
   else if (nonlinear==0) res *= P_L(a, f);
   else assert(0);

   if (!finite(res)) assert(0);

   return res;
}

/* S98  3.4 */
real Pkappa(real s)
{
   static real OMEGA_M   = -42.;
   static real OMEGA_V   = -42.;
   static real N_SPEC    = -42.;
   static real GAMMA     = -42.;
   static real BETA_P    = -42.;
   static real Z0        = -42.;
   static real NONLINEAR = -42;
   static real SIGMA_8   = -42.;
   static real table[N_s];
   static real ds = .0, logsmin = .0, logsmax = .0;
   real   ss, slog, f1, f2;
   int    i;

   if (OMEGA_M != Omega_m || OMEGA_V != Omega_v || N_SPEC != n_spec ||
       GAMMA != Gamma || BETA_P != beta_p || Z0 != z0 ||
       NONLINEAR != nonlinear || SIGMA_8 != sigma_8) {
      logsmin = log(s_min);
      logsmax = log(s_max);
      ds = (logsmax - logsmin)/(N_s - 1.);
      slog = logsmin;
      for (i=0; i<N_s; i++, slog+=ds) {
	 ss = exp(slog);
	 sglob = ss;
	 f1 = qromb(int_for_p_2, a_min, 0.7);
	 f2 = qromo(int_for_p_2, .7, 1.0, midpnt);
	 table[i] = log(9./4.*DSQR(Omega_m)*(f1 + f2));
      }
      OMEGA_M =   Omega_m;
      OMEGA_V =   Omega_v;
      N_SPEC  =   n_spec;
      GAMMA   =   Gamma;
      BETA_P  =   beta_p;
      Z0      =   z0;
      NONLINEAR = nonlinear;
      SIGMA_8 =  sigma_8;
   }
   slog = log(s);
   f1 = interpol(table, N_s, logsmin, logsmax, ds, slog, 0.0, 0.0);
   return exp(f1);
}
