#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "smith2.h"

#define HL_HIMPC 2997.92458

typedef struct {
  int cosmotype;
  double TCMB, omega_B, omega_nu, omega_M, omega_de, w0, wa, h;
  double sigma8, ns, alphas;
  /*JD Parameters for testing GR*/
  double TGR_Q0, TGR_DR0, TGR_s;
  int TGRTDEP, TGR_RFUNC;
  
  double A; /* amplitude used for PS normalization -- call DNnormalize first */

  unsigned long zflag; /* bits for z-determination flag */
} COSMOPARAM;

/* Reads cosmological parameters structure */
void DNread_cosmo(COSMOPARAM *p, char FileName[]) {
  FILE *fp;
  fp = fopen(FileName, "r");
  fscanf(fp,  "%d", &(p->cosmotype));
  fscanf(fp, "%lg", &(p->TCMB));
  fscanf(fp, "%lg", &(p->omega_B));
  fscanf(fp, "%lg", &(p->omega_nu));
  fscanf(fp, "%lg", &(p->omega_M));
  fscanf(fp, "%lg", &(p->omega_de));
  fscanf(fp, "%lg", &(p->w0));
  fscanf(fp, "%lg", &(p->wa));
  fscanf(fp, "%lg", &(p->h));
  fscanf(fp, "%lg", &(p->sigma8));
  fscanf(fp, "%lg", &(p->ns));
  fscanf(fp, "%lg", &(p->alphas));
  /*JD*/
  fscanf(fp, "%lg", &(p->TGR_Q0));
  fscanf(fp, "%lg", &(p->TGR_DR0));
  fscanf(fp, "%lg", &(p->TGR_s));
  fscanf(fp,  "%d", &(p->TGRTDEP));
  fscanf(fp,  "%d", &(p->TGR_RFUNC));
  fclose(fp);
}
/* JD for testing GR this function is in the source term of the growth equation*/
double geff(COSMOPARAM *p, double z){
  double TGR_q,TGR_dr,a;
  a = 1.0/(1.0+z);
  
  if(p->TGRTDEP==0){
    TGR_q = p->TGR_Q0;
    TGR_dr = p->TGR_DR0;
  }
  else{
    TGR_q=(p->TGR_Q0 - 1.0)*pow(a,p->TGR_s)+1.0;
    TGR_dr=(p->TGR_DR0 - 1.0)*pow(a,p->TGR_s)+1.0;
  }
  if(p->TGR_RFUNC ==1){
    return(TGR_q*TGR_dr);
  }
  else{
    return(2.0*TGR_dr - TGR_q);
  }
}

/* Obtains Hubble rate at z, H(z)/H(0) */
double DNhubble(COSMOPARAM *p, double z) {
  double y = 1.+z;
  double E;
  double omega_K = 1.-p->omega_M-p->omega_de;
  E = pow(y,3*(1+p->w0+p->wa)) * exp(-3*p->wa*z/y);
  return(sqrt(y*y*(omega_K+y*p->omega_M)+p->omega_de*E));
}

/* Obtains the radial distance from z1 .. z2 in h^-1 Mpc */
double DNget_chi1(COSMOPARAM *p, double z1, double z2) {
  int i, ii;
  int NGauss = 10;
  double x[]= {0.148874338981631, 0.433395934129247, 0.679409568299024, 0.865063366688985, 0.973906528517172};
  double w[]= {0.295524224714753, 0.269266719309996, 0.219086362515982, 0.149451349150581, 0.066671344308688};
  double distance = 0.;
  double z;

  for(i=0; i<10; i++) {
    ii = i%(NGauss/2);
    z = 0.5*(z1+z2) + 0.5*(z2-z1)*(i>=NGauss/2? x[ii]: -x[ii]);
    distance += w[ii]/DNhubble(p,z);
  }
  distance *= (z2-z1)/2.;
  return(distance * HL_HIMPC);
}

/* Obtains the radial distance from z1 .. z2 in h^-1 Mpc */
double DNget_chi11(COSMOPARAM *p, double z1, double z2) {
  int i,step;
  double za,zb,ratio,lratio;
  double result=0.;

  step = 1+(int)floor(0.5*fabs(log((1.+z2)/(1.+z1))));
  lratio = log((1.+z2)/(1.+z1))/(double)step;
  ratio = exp(lratio);

  for(i=0; i<step; i++) {
    za = (1.+z1)*exp(lratio*i) - 1.;
    zb = (1.+za)*ratio - 1.;
    result += DNget_chi1(p,za,zb);
  }
  return(result);
}

/* Obtains the radial distance from z1 .. z2 in h^-1 Mpc.
 * This routine stores values from DNget_chi11 in order
 * to speed up the calculation.
 */
double DNget_chi(COSMOPARAM *p, double z1, double z2) {
  static int is_initialized = 0;
  static double distances[1001];
  long i, i1, i2;
  double ii1, ii2;

#ifdef SAVE_TABLES
  if (!is_initialized) {
    is_initialized=1;
    for(i=0;i<=1000;i++)
      distances[i] = DNget_chi11(p,0.,0.005*i);
  }

  ii1 = z1/0.005;
  ii2 = z2/0.005;
  i1 = (long)floor(ii1 + 0.5);
  i2 = (long)floor(ii2 + 0.5);
  if (ii1>=0 && ii1<=1000 && ii2>=0 && ii2<=1000 && fabs(ii1-i1)<1e-7 && fabs(ii2-i2)<1e-7)
    return(distances[i2]-distances[i1]);
#endif

  return(DNget_chi11(p,z1,z2));
}

/* Obtains the angular diameter distance from z1 .. z2 in h^-1 Mpc */
double DNget_r(COSMOPARAM *p, double z1, double z2) {
  double K, sqrtaK, chi, corr;
  double omega_K = 1.-p->omega_M-p->omega_de;

  chi = DNget_chi(p, z1, z2);
  K = -omega_K/(HL_HIMPC*HL_HIMPC);
  corr = K*chi*chi;
  if (fabs(corr)<1.0e-6) {
    return(chi*(1.-corr/6.+corr*corr/120.));
  }
  sqrtaK = sqrt(fabs(K));
  if (K>0)
    return(sin(sqrtaK*chi)/sqrtaK);
  if (K<0)
    return(sinh(sqrtaK*chi)/sqrtaK);
  return(0.);
}

/* Growth function derivative */
void DNgrowth_derivative(COSMOPARAM *p, double lny, double *X, double *dXdy, int flag) {
  static double z, H, y, omH;

  if (flag) {
    y = exp(lny);
    z = y-1;
    H = DNhubble(p,z);
    omH = p->omega_M*y*y*y/H;
  }
  dXdy[0] = X[1]/H+X[0];
  dXdy[1] = 3*X[1]+1.5*geff(p,z)*omH*X[0];  /*JD added geff for testing GR*/
}

void DNgrowth_step(COSMOPARAM *p, double lny, double *Xin, double *Xout, double dlny, int flag) {
  double k1[2], k2[2], k3[2], k4[2], X[2];
  int i;
  double half_dlny = 0.5*dlny;
  double lnyh = lny+0.5*dlny;

  /* Get derivatives */
  DNgrowth_derivative(p, lny, Xin, k1, flag);

  for(i=0;i<2;i++) X[i] = Xin[i]+half_dlny*k1[i];
  DNgrowth_derivative(p, lnyh, X, k2, 1);

  for(i=0;i<2;i++) X[i] = Xin[i]+half_dlny*k2[i];
  DNgrowth_derivative(p, lnyh, X, k3, 0);

  for(i=0;i<2;i++) X[i] = Xin[i]+dlny*k3[i];
  DNgrowth_derivative(p, lny+dlny, X, k4, 1);

  for(i=0;i<2;i++) Xout[i] = Xin[i] + dlny/6.*(k1[i]+2*(k2[i]+k3[i])+k4[i]);
}

/* Extrapolated growth function.  DNget_growth_normhiz1 is the actual function;
 * DNget_growth_normhiz is used only to store a table of commonly used values.
 *
 * In its current version, stores NZ1 values from ZSTEP1*0.5 .. ZSTEP1(-0.5+NZ1).
 */
double DNget_growth_normhiz(COSMOPARAM *p, double z) {
#define NZ1 500
#define ZSTEP1 0.01
  static int is_initialized=0;
  static double growth[NZ1];
  long i;
  double ii;

  double DNget_growth_normhiz1(COSMOPARAM*, double);

#ifdef SAVE_TABLES
  /* Make table */
  if (!is_initialized) {
    is_initialized = 1;
    for(i=0; i<NZ1; i++)
      growth[i] = DNget_growth_normhiz1(p, (0.5+i)*ZSTEP1);
  }

  /* Return growth function from table if we're on one of those values */
  ii = z/ZSTEP1 - 0.5;
  i = (long)floor(ii+0.5);
  if (fabs(ii-i)<1e-8 && i>=0 && i<NZ1)
    return (growth[i]);
#endif

  /* Otherwise recalculate */
  return(DNget_growth_normhiz1(p,z));  
#undef NZ1
#undef ZSTEP1
}

double DNget_growth_normhiz1(COSMOPARAM *p, double z) {
  double Xold[2], Xnew[2];
  double lny, dlny=-0.3;
  long N=20;
  long i;

  lny = log(1.+z) - N*dlny;
  Xold[0] = 1.;
  Xold[1] = -DNhubble(p,exp(lny)-1);

  for(i=0; i<N; i++) {
    DNgrowth_step(p, lny, Xold, Xnew, dlny, i? 0: 1);
    lny += dlny;
    Xold[0] = Xnew[0]; Xold[1]=Xnew[1];
  }
  return(Xnew[0]/(1.+z));
}

/* Growth function normalized to today */
double DNget_growth_normtoday(COSMOPARAM *p, double z) {
  return(DNget_growth_normhiz(p,z)/DNget_growth_normhiz(p,0));
}

/* Transfer fcn normalized to 1 at k=0 (Eisenstein/Hu)
 * Returns: T(k) in matter era
 */
double DNeh_trf(COSMOPARAM *p, double k) {

  static double th27, kmpc, omh2, obh2;
  static double q, alphac, betac, a1, a2, b1, b2, bm, zd, zeq, keq, s, req, rd;
  static double f, C0, T0_k1b, T0_kab, Tc;
  static double ksilk, y, alphab, betab, betanode__ks, betab__ks, Tb, stilde, T;
  static double old_h=-1, old_T=-1, old_om_M=-1, old_om_nu=-1, old_om_B=-1;
  static double xsupp, xnumer, kmpc0;

  /* Re-initialize constants if cosmology has changed */
  if (old_h!=p->h || old_T!=p->TCMB || old_om_M!=p->omega_M || old_om_nu!=p->omega_nu || old_om_B!=p->omega_B) {
    /* Update variables */
    old_h=p->h;
    old_T=p->TCMB;
    old_om_M=p->omega_M;
    old_om_nu=p->omega_nu;
    old_om_B=p->omega_B;

    /* Numbers to re-compute */
    th27 = p->TCMB / 2.7;
    omh2 = (p->omega_M-p->omega_nu) * p->h * p->h;
    obh2 = p->omega_B * p->h * p->h;

    /* redshift at decoupling, equality, sound horizon */
    b1 = 0.313 * pow(omh2, -0.419) * ( 1 + 0.607*pow(omh2, 0.674) );
    b2 = 0.238 * pow(omh2, 0.223);
    zd = 1291. * pow(omh2, 0.251) / ( 1 + 0.659*pow(omh2, 0.828) ) * ( 1 + b1*pow(obh2,b2) );
    zeq = 25000. * omh2 / pow(th27,4.);
    keq = 0.0746 * omh2 / th27 / th27; /* in Mpc^-1 */
    rd = 31500.*obh2 / pow(th27,4.) / zd;
    req = zd*rd/zeq;
    s = 1.632993161855/keq/sqrt(req) * log((sqrt(1+rd)+sqrt(rd+req))/(1+sqrt(req)));

    /* EH parameters */
    a1 = pow(46.9*omh2, 0.670) * ( 1 + pow(32.1*omh2, -0.532) );
    a2 = pow(12.0*omh2, 0.424) * ( 1 + pow(45.0*omh2, -0.582) );
    b1 = 0.944 / ( 1 + pow(458*omh2, -0.708) );
    b2 = pow( 0.395*omh2, -0.0266);
    bm = obh2/omh2;
    alphac = pow(a1, -bm) * pow(a2, -bm*bm*bm);
    betac = 1./( 1 + b1*(pow(1-bm, b2)-1) );

    /* k-independent baryon parameters */
    ksilk = 1.6 * pow(obh2, 0.52) * pow(omh2, 0.73) * ( 1 + pow(10.4*omh2, -0.95) );
    y = (1.+zeq)/(1.+zd);
    alphab = 2.07*keq*s*pow(1+rd,-0.75) * y * ( -6.*sqrt(1+y) + (2.+3.*y)*log((sqrt(1.+y)+1.)/(sqrt(1.+y)-1.)) );
    betab = 0.5 + bm + (3.-2.*bm)*sqrt(295.84*omh2*omh2+1.);

    /* More parameters */
    xnumer = 8.41*pow(omh2,0.435)/s;
    kmpc0 = omh2/(th27*th27);
  }

  /* wavenumber in Mpc^-1 comoving, no h */
  kmpc = k * p->h;

  /* CDM piece */
  q = kmpc/kmpc0;
  f = 1./(1. + pow(kmpc*s/5.4, 4.));
  C0 = 386./(1+69.9*pow(q,1.08));
  xsupp = q*q/log( M_E + 1.8*betac*q );
  T0_kab = 1./(1. + (C0+14.2/alphac)*xsupp);
  T0_k1b = 1./(1. + (C0+14.2)*xsupp);
  Tc = f*T0_k1b + (1-f)*T0_kab;

  /* Baryonic piece */
  betanode__ks = xnumer/kmpc;
  betab__ks = betab/kmpc/s;
  stilde = s*pow(1.+betanode__ks*betanode__ks*betanode__ks, -0.33333333);
  Tb = 1./(1. + (C0+14.2)*q*q/log( M_E + 1.8*q ))/(1+kmpc*kmpc*s*s/27.04)
       + alphab/(1+betab__ks*betab__ks*betab__ks)*exp(-pow(kmpc/ksilk,1.4));
  Tb *= sin(kmpc*stilde)/(kmpc*stilde);

  T = bm*Tb + (1.-bm)*Tc;

  return(T);
}

/* Neutrino suppression factor in transfer function, if we have massive neutrinos */
double DNneutrino_suppression(COSMOPARAM *p, double k, double D) {
  double th27, fnu, pcb, omh2, zeq;
  double q, yfs, qnf, Tfactor;

  fnu = p->omega_nu/p->omega_M;
//  if (fnu<1e-4)
    return(1.);
  pcb = 0.25*(5.-sqrt(25.-24.*fnu));

  th27 = p->TCMB / 2.7;
  omh2 = (p->omega_M-p->omega_nu) * p->h * p->h;
  zeq = 25000. * omh2 / pow(th27,4.);
  q = k*th27*th27/p->omega_M;
  qnf = 3*q/fnu;
  yfs = 17.2*fnu*(1.+0.488*pow(fnu,-1.167))*qnf*qnf;

  Tfactor = pow( pow(1.-fnu,0.7/pcb)*pow(D*(1.+zeq),-0.7) + pow(1.+yfs,-0.7), pcb/0.7);

  /* Diagnostics */
#if 0
printf("fnu=%lf, k=%le, q=%le, yfs=%le, pcb=%lf T=%lf (%lf)\n", fnu, k, q, yfs, pcb, Tfactor, pow(1.+yfs,-pcb));
exit(0);
#endif

  return(Tfactor);
}

/* Power spectrum, not normalized */
double DNpk_unnorm(COSMOPARAM *p, double k, double z) {
  double D, T, k0;
  k0 = 0.05/p->h;
  D = DNget_growth_normhiz(p,z);
  T = DNeh_trf(p,k)*DNneutrino_suppression(p,k,D);
  return(D*D*T*T*pow(k,p->ns+0.5*p->alphas*log(k/k0)));
}

/* sigma_8 normalization integral */
double DNsig8int(COSMOPARAM *p) {
  double k, dlnk, expdlnk, T, y, var;
  dlnk = 0.1;
  expdlnk = exp(dlnk);
  var=0;
  for(k=1e-3; k<1e2; k*=expdlnk) {
    y = 8.*k;
    T = y>1e-5? 3./(y*y)*(sin(y)/y-cos(y)): 1.-y*y/10.*(1.-y*y/28.);
    var += DNpk_unnorm(p,k,0.)*T*T*k*k*k;
  }
  var *= dlnk/(2.*M_PI*M_PI);
  return(var);
}

/* Normalization constant for cosmology */
void DNnormalize(COSMOPARAM *p) {
  p->A = p->sigma8*p->sigma8/DNsig8int(p);
}

/* Normalized linear power spectrum */
double DNpk_lin(COSMOPARAM *p, double k, double z) {
  return(p->A*DNpk_unnorm(p,k,z));
}

/* Normalized linear power spectrum */
double DNpk_nonlin(COSMOPARAM *p, double k, double z) {
  return(P_NL(1./(1.+z), k*HL_HIMPC)*HL_HIMPC*HL_HIMPC*HL_HIMPC);
}

