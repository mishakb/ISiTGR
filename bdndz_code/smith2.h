/* ============================================================ *
 * smith2.h							*
 * MK 1/2006							*
 * ============================================================ */


#ifndef __SMITH2_H
#define __SMITH2_H

#define N_a     (100)
#define N_k     (100)
#define N_s     (200)
#define N_theta (100)

#define a_min (0.01)
#define k_min (1.e-2)
#define k_max (1.e6)
#define s_min (1.e-2)
#define s_max (1.e6)

#define pi     (3.14159265358979323846)
#define pi_sqr (9.86960440108935861883)
#define twopi  (6.28318530717958647693)

#define epsilon0 (1E-2)
#define epsilon  (1E-5)
#define epsilon1 (1E-10)
#define epsilon2 (1E-18)

typedef double real;

static real darg __attribute__((unused)), maxarg1 __attribute__((unused)), maxarg2 __attribute__((unused)),
  darg2 __attribute__((unused));
static int iminarg1 __attribute__((unused)), iminarg2 __attribute__((unused)), iarg __attribute__((unused));
static long long Lminarg1 __attribute__((unused)), Lminarg2 __attribute__((unused));

#define IMIN(a,b) (iminarg1=(a), iminarg2=(b), (iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))
#define LMIN(a,b) (Lminarg1=(a), Lminarg2=(b), (Lminarg1) < (Lminarg2) ? (Lminarg1) : (Lminarg2))
#define IMAX(a,b) (iminarg1=(a), iminarg2=(b), (iminarg1) > (iminarg2) ? (iminarg1) : (iminarg2))
#define FMAX(a,b) (maxarg1=(a), maxarg2=(b), (maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))
#define ISQR(a) ( (iarg=(a))==0.0 ? 0.0 : iarg*iarg )
#define DSQR(a) ( (darg=(a))==0.0 ? 0.0 : darg*darg )

#define NR_END 1
#define FREE_ARG char*


/* ============================================================ *
 * Functions available for calls from fortran.			*
 * ============================================================ */

void pl_(real *a, real *k, real *res);
void pnl_(real *a, real *k, real *res);
void pkappa_(real *s, real *res);
void setparameters_(real *OMEGAM, real *OMEGAV, real *GAMMA, real *SIGMA8,
		     real *NSPEC, real *BETAP, real *Z0, int *NONLINEAR);

/* ============================================================ */


void set_cosmological_parameters_to_default();
void setparameters(real *, real *, real *, real *,
		    real *, real *, real *, int *);
void dump_param(FILE *F);

void error(const char *s);

real **matrix(long nrl, long nrh, long ncl, long nch);
real dfridr(real (*func)(real,real), real x, real h, real *err, real aa);
void spline(real x[], real y[], int n, real yp1, real ypn, real y2[]);
void splint(real xa[], real ya[], real y2a[], int n, real x, real *y);
void polint(real xa[], real ya[], int n, real x, real *y, real *dy);
real trapzd(real (*func)(real), real a, real b, int n);
real trapzd1(real (*func)(real), real a, real b, int n);
real interpol(real *f, int n, real a, real b, real dx, real x,
	      real lower, real upper);
real interpol2d(real **f,
		int nx, real ax, real bx, real dx, real x,
		int ny, real ay, real by, real dy, real y,
		real lower, real upper);
real gammln(real xx);
real qromb(real (*func)(real), real a, real b);
real qromb1(real (*func)(real), real a, real b);
real qromo(real (*func)(real), real a, real b,
	   real (*choose)(real(*)(real), real, real, int));
real midpnt(real (*func)(real), real a, real b, int n);


real D_plus(real a);
real Tsqr(real k);
real sigma_8_sqr();
real sigma_R_sqr(real);
real P_L(real, real);
real P_NL(real, real);

void Omega_a(real, real*, real*);
real g(real a);
real Delta_L_BE2(real k);
real Delta_L_BE(real k);
real dlog(real x);
void halofit(real rk,real rn,real rncur,real rknl,real plin,real om_m, real om_v,real *pnl);
void wint2(real r,real *sig,real *d1,real *d2, real amp, int onlysig);

real g_source(real a);
real f_K(real w);
real G(real a);
real int_for_w(real);
real w(real a);
real prob(real z);
real Pkappa(real);

/* ============================================================ *
 * cosmological parameter					*
 * ============================================================ */

real Omega_m;		/* matter density parameter */
real Omega_v;		/* cosmogical constant parameter */

real Gamma;		/* shape parameter, see PD2 51 */
real sigma_8;		/* power spectrum normalization */

real n_spec;		/* spectral index of initial power spectrum */

real z0;		/* redshift distribution scale */
real beta_p;		/* redshift distribution power index. If beta_b = 0, a distribution
			 * with a single source redshift at z0 is assumed */
real sigma_epsilon;	/* ellipticity dispersion */
int  nonlinear;         /* 0: linear power spectrum                    */

#endif
