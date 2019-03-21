#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "smith2.h"
#include "utils_func_scale.c"

#define ZSTEP 0.01

/* NZ *must* have NZ equiv 2 mod 3 */
#define NZ 500
#define NSAMPLE 9

#include "read.c"
#include "linalg.c"
#include "min1d.c"

/*JD functions for testing General Relativity*/

/*JD for Testing GR this is to be used when calculating the magnification bias
 where modification to the convergence is considered only at large k*/
double TGR_D_HK(COSMOPARAM *p, double z){
	double TGR_q,TGR_dr,a;
	a = 1.0/(1.0+z);
	if(p->TGRTDEP==0){
		TGR_q = p->TGR_Qinf;
		TGR_dr = p->TGR_DRinf;
	}
	else{
       	TGR_q=(p->TGR_Qinf - 1.0)*pow(a,p->TGR_s)+1.0;
       	TGR_dr=(p->TGR_DRinf - 1.0)*pow(a,p->TGR_s)+1.0;
	}
	
	if(p->TGR_RFUNC ==1){
		return(TGR_q*(1.0+TGR_dr)/2.0);
	}
	else{
		return(TGR_dr);
	}
}
/* Gets lensing kernel W: kappa = int W delta dzl */

/* JD edited for Testing GR: this function gets the lensing kernel
   for calculation of the magnification bias.  As discussed in 
   arxiv:astro-ph/0004392v1 this kernel should be calculated at large k values.
   When testing convergence has a factor of TGR_D in it so we put
   that here in the lensing kernel.  
   For scale dependent case, we calculate TGR_D as valid at large k values.*/
double lensingkernel_mag_bias(double zl, double zs, COSMOPARAM *p) {
  double W;

  if (zl>=zs) return(0.);
  /*JD multiply by TGR_D_HK  (large k)*/ 
  W = TGR_D_HK(p,zl) *1.5 * p->omega_M * DNget_r(p,0,zl) * DNget_r(p,zl,zs) * (1.+zl) / DNget_r(p,0,zs) / DNhubble(p,zl) / HL_HIMPC;
  return(W);
}

/*JD end new functions for testing GR

/* "Repairs" a redshift distribution to include magnification bias, given alpha = -dlogN/dkappa
 * Keeps only the lensing curve if flag=1.
 */
void repairdndz(double *dist_in, double *dist_out, double alpha, COSMOPARAM *p, int flag) {
  int i,j;
  double dist_temp[NZ];

  alpha *= ZSTEP;
  for(i=0; i<NZ; i++) {
    if (flag!=1) {
      dist_temp[i] = dist_in[i];
    } else {
      dist_temp[i] = 0.0;
    }
    for(j=i+1; j<NZ; j++)
      if (dist_in[j]!=0)     
      /*JD for testing GR modified next line lensingkernel->lensingkernel_mag_bias*/
        dist_temp[i] += alpha*lensingkernel_mag_bias((i+0.5)*ZSTEP,(j+0.5)*ZSTEP,p)*dist_in[j];
  }
  for(i=0; i<NZ; i++) dist_out[i]=dist_temp[i];
}

/* If out1 exists and is !=out0:
 *   Gets power spectrum @ L, adds b=1,Q=0 ps to out0, and b=1,(Q=1)-(Q=0) ps to out1
 * If out1==NULL:
 *   nonlinear PS -> out0, ignores out1
 * If out1==out0:
 *   linear PS -> out0, ignores out1
 */
void angularpower(double L, double *dist1, double *dist2, COSMOPARAM *p, double *out0, double *out1) {
  int i;
  double z, chi, r, dchidz;
  double k, Pk, prod;
  L += 0.5;
  for(i=0; i<NZ; i++) {
    prod = dist1[i]*dist2[i];
#if 1
    prod = i<=49? dist1[i]*dist2[i]:
           i%3==0? dist1[i-1]*dist2[i-1]+dist1[i]*dist2[i]+dist1[i+1]*dist2[i+1]: 0;
#endif
    if (prod!=0) {
      z = (i+0.5)*ZSTEP;
      chi = DNget_chi(p,0,z);
      r = DNget_r(p,0,z);
      dchidz = HL_HIMPC/DNhubble(p,z);
      k = L/r;
      Pk = out1==NULL? DNpk_nonlin(p,L/r,z): DNpk_lin(p,L/r,i);//JD For Testing GR changed the z input to an integer
      if (out1!=NULL && out1!=out0) {						   //indicating z table position for interpolation
        Pk /= (1.+1.7*k);
        *out1 += Pk/(r*r)*ZSTEP*prod/dchidz * k*k;
      }
      *out0 += Pk/(r*r)*ZSTEP*prod/dchidz;
    }
  }
}

/* Gets chi^2 for comparing the 2 distributions given against the cross-power in
 * the file.  The file format is:
 *
 * nbin, nl  (nl = #l's used to specify window function)
 * ps, sigma, norm, L1...Lnl (nbin lines)
 * covariance matrix (nbin lines)
 *
 * out = (e, ex, ey, exx, exy, eyy)
 * where e = chi^2(0,0), ex,ey = grad, and grad^2
 * parameters are b^2 and b^2*Q.
 */
void getchi2(char FileName[], double *dist1, double *dist2, COSMOPARAM *p, double *out) {

  FILE *fp;
  int i, j, k, nbin, nl;
  double *ps, *sigma, *predict0, *predict1, **cov, *norm, **L;

  /* Open file, setup */
  fp = fopen(FileName, "r");
  fscanf(fp, "%d %d", &nbin, &nl);
  ps = dvector(0,nbin-1);
  sigma = dvector(0,nbin-1);
  predict0 = dvector(0,nbin-1);
  predict1 = dvector(0,nbin-1);
  norm = dvector(0,nbin-1);
  cov = dmatrix(0,nbin-1,0,nbin-1);
  L = dmatrix(0,nbin-1,0,nl-1);

  /* Get power spectra from data */
  for(i=0; i<nbin; i++) {
    fscanf(fp, "%lg %lg %lg", ps+i, sigma+i, norm+i);
#ifndef USENORM
    norm[i]=1.;
#endif
    for(k=0; k<nl; k++) fscanf(fp, "%lg", L[i]+k);
  }
  for(i=0; i<nbin; i++) for(j=0; j<nbin; j++) fscanf(fp, "%lg", cov[i]+j);

  /* Invert covariance matrix */
  gaussjinv(cov,nbin);

  /* Get prediction */
  for(i=0; i<nbin; i++) {
    predict0[i] = 0.;
    predict1[i] = 0.;
    for(k=0; k<nl; k++) {
      angularpower(L[i][k],dist1,dist2,p, predict0+i, predict1+i);
    }
    predict0[i] *= norm[i]/(double)nl;
    predict1[i] *= norm[i]/(double)nl;
  }

  /* Get chi^2 */
  out[0]=out[1]=out[2]=out[3]=out[4]=out[5]=0.;
  for(i=0; i<nbin; i++) for(j=0; j<nbin; j++) {
    out[0] += cov[i][j]*ps[i]*ps[j];
    out[1] += -2*cov[i][j]*ps[i]*predict0[j];
    out[2] += -2*cov[i][j]*ps[i]*predict1[j];
    out[3] += 2*cov[i][j]*predict0[i]*predict0[j];
    out[4] += 2*cov[i][j]*predict0[i]*predict1[j];
    out[5] += 2*cov[i][j]*predict1[i]*predict1[j];
  }

  /* Cleanup */
  fclose(fp);
  free_dvector(ps,0,nbin-1);
  free_dvector(sigma,0,nbin-1);
  free_dvector(predict0,0,nbin-1);
  free_dvector(predict1,0,nbin-1);
  free_dvector(norm,0,nbin-1);
  free_dmatrix(cov,0,nbin-1,0,nbin-1);
  free_dmatrix(L,0,nbin-1,0,nl-1);
}

/* Gets chi^2 & its coefficients for quasars or NVSS sample:
 * chi^2 = out[0] + out[1]*x + out[2]*x^2/2
 *
 * where the redshift distribution is multiplied by x.
 *
 * We subtract "poisson" from the C_l's before using them to account for objects matched to both catalogs.
 *
 * Returns chi^2 for x=1.
 */
double getchi2_1p(char FileName[], double *dist1, double *dist2, COSMOPARAM *p, double *out, double poisson) {

  FILE *fp, *fq;
  int i, j, k, nbin, nl;
  double *ps, *sigma, *predict0, **cov, *norm, **L;

  /* Open file, setup */
  fp = fopen(FileName, "r");
  fscanf(fp, "%d %d", &nbin, &nl);
  ps = dvector(0,nbin-1);
  sigma = dvector(0,nbin-1);
  predict0 = dvector(0,nbin-1);
  norm = dvector(0,nbin-1);
  cov = dmatrix(0,nbin-1,0,nbin-1);
  L = dmatrix(0,nbin-1,0,nl-1);

#ifdef VERBOSE_OUT
  fq = fopen("chi2.log", "a");
  fprintf(fq, "%s\n", FileName);
#endif

  /* Get power spectra from data */
  for(i=0; i<nbin; i++) {
    fscanf(fp, "%lg %lg %lg", ps+i, sigma+i, norm+i);
#ifndef USENORM
    norm[i]=1.;
#endif
    ps[i] -= norm[i]*poisson;
    for(k=0; k<nl; k++) fscanf(fp, "%lg", L[i]+k);
  }
  for(i=0; i<nbin; i++) for(j=0; j<nbin; j++) fscanf(fp, "%lg", cov[i]+j);

  /* Invert covariance matrix */
  gaussjinv(cov,nbin);

  /* Get prediction */
  for(i=0; i<nbin; i++) {
    predict0[i] = 0.;
    for(k=0; k<nl; k++) {
      angularpower(L[i][k],dist1,dist2,p, predict0+i, predict0+i);
    }
    predict0[i] *= norm[i]/(double)nl;
#ifdef VERBOSE_OUT
    fprintf(fq, "  %6.1lf %12.5le %12.5le %12.5le %7.2lf\n", L[i][nl/2], ps[i], sigma[i], predict0[i], (ps[i]-predict0[i])/sigma[i]);
#endif
  }

  /* Get chi^2 */
  out[0]=out[1]=out[2]=0.;
  for(i=0; i<nbin; i++) for(j=0; j<nbin; j++) {
    out[0] += cov[i][j]*ps[i]*ps[j];
    out[1] += -2*cov[i][j]*ps[i]*predict0[j];
    out[2] += 2*cov[i][j]*predict0[i]*predict0[j];
  }

  /* Cleanup */
  fclose(fp);
  free_dvector(ps,0,nbin-1);
  free_dvector(sigma,0,nbin-1);
  free_dvector(predict0,0,nbin-1);
  free_dvector(norm,0,nbin-1);
  free_dmatrix(cov,0,nbin-1,0,nbin-1);
  free_dmatrix(L,0,nbin-1,0,nl-1);
#ifdef VERBOSE_OUT
  fclose(fq);
#endif
  return(out[0]+out[1]+0.5*out[2]);
}

void printalldist(char FileName[], double **dist) {
  FILE *fp;
  int i,iz;

  fp = fopen(FileName, "w");
  for(iz=0; iz<NZ; iz++) {
    fprintf(fp, "%6.4lf", (iz+0.5)*ZSTEP);
    for(i=0; i<NSAMPLE; i++) fprintf(fp, " %7.4lf", dist[i][iz]);
    fprintf(fp, "\n");
  }
  fclose(fp);
}

/* Obtains total chi^2 given distribution, low-z samples.
 * ip = sample number
 */
void getallchi2(double **dist, COSMOPARAM *p, int ip, double *out) {

  /* The auto-powers */
  if (ip==0)
    getchi2("data/2mass00.dat", dist[0], dist[0], p, out);
  if (ip==1)
    getchi2("data/2mass11.dat", dist[1], dist[1], p, out);
  if (ip==2)
    getchi2("data/2mass22.dat", dist[2], dist[2], p, out);
  if (ip==3)
    getchi2("data/2mass33.dat", dist[3], dist[3], p, out);
  if (ip==4)
    getchi2("data/auto-lrg0.dat", dist[4], dist[4], p, out);
  if (ip==5)
    getchi2("data/auto-lrg1.dat", dist[5], dist[5], p, out);
}

void construct_dist_loz(COSMOPARAM *p, double **dist) {
  int i,iz;
  double out[6];
  double e, e1, e2, e11, e12, e22;
  double det, inv11, inv12, inv22;
  double x, y, chi2min;

  /* Read sample data */
  get_2mass(dist);
  get_lrg(dist+4);

  /* For each sample, get chi^2 surface and find minimum */
  for(i=0; i<6; i++) {

    getallchi2(dist, p, i, out);
    e   = out[0];
    e1  = out[1];
    e2  = out[2];
    e11 = out[3];
    e12 = out[4];
    e22 = out[5];

    /* Inverse of 2x2 matrix */
    det = e11*e22-e12*e12;
    inv11 =  e22/det;
    inv12 = -e12/det;
    inv22 =  e11/det;

    /* chi^2 minimum */
    x = -(inv11*e1+inv12*e2);
    y = -(inv12*e1+inv22*e2);
    chi2min = e + 0.5*(e1*x+e2*y);

#if 1
    /* Report bias */
    /*fprintf(stderr, "<%1d> bias=%9.6lf +/- %9.6lf\n", i, sqrt(x), sqrt(inv11/2./x));*/
    /* Report Q */
    /*fprintf(stderr, "    Q = %9.6lf +/- %9.6lf\n", y/x, y/x*sqrt(2.)*sqrt(inv11/x/x-2*inv12/x/y+inv22/y/y));*/
#endif

    /* Multiply the dN/dz by the bias */
    x=sqrt(x);

    for(iz=0; iz<NZ; iz++) dist[i][iz] *= x;
  }
}

/* Constructs QSO redshift distribution and returns chi^2.
 * Input is zp[0..1], which controls redshift distributions.
 *
 * Currently
 * zp[0] = bias * growth factor, main peak
 * zp[1] = bias * growth factor, secondary peak
 *
 * Critical redshift is 1.18 (QSO1)
 */
double construct_test_dist_qso(COSMOPARAM *p, double **alldist, double *dist, double *magdist, double *refdist, int slice, double *zp) {
  int iz, ilrg;
  double out[3];
  char FileName[1024];
  double chi2;

  /* Fixes the redshift distribution for magnification bias */
  if (slice==0)
    for(iz=0; iz<NZ; iz++)
      dist[iz] = zp[0]*refdist[iz] -0.18*magdist[iz];
  if (slice==1)
    for(iz=0; iz<NZ; iz++)
      dist[iz] = (iz>=118? zp[0]: zp[1])*refdist[iz] -0.10*magdist[iz];

  /* Get the chi^2 */
  sprintf(FileName, "data/qso%1d-auto.dat", slice);
  chi2 = getchi2_1p(FileName, dist, dist, p, out, 0);

  for(ilrg=0; ilrg<2; ilrg++) {
    sprintf(FileName, "data/qso%1d-lrg%1d.dat", slice, ilrg);
    chi2 += getchi2_1p(FileName, dist, alldist[4+ilrg], p, out, 0);
  }

  return(chi2);
}

#define NZP_QSOMAX 2
void construct_dist_qso(COSMOPARAM *p, double **dist) {
  int offset=6;
  int NZP_QSO;
  int iz, doffset, ip, id;
  long istep=0;
  double chi2, norm, renorm, z, growth;
  double y[5], xmin;
  int ix;
  double zp[NZP_QSOMAX], zp_old[NZP_QSOMAX];
  double delta[NZP_QSOMAX];
  double **magdist, **refdist, **dirs;
  double c2b, chi2_old;
  double dchi2, c2_start, c2_end;
  int dirbest;
  double norm0;
  double t,dt,zptest[NZP_QSOMAX],chi2test,signtest;

  /* Allocations */
  magdist = dmatrix(0, 1, 0, NZ-1);
  refdist = dmatrix(0, 1, 0, NZ-1);
  dirs = dmatrix(0, NZP_QSOMAX-1, 0, NZP_QSOMAX-1);

  /* Get spectro-confirmed dN/dz */
  get_qso(dist+offset);

  /* Get magnification bias distribution */
  for(doffset=0; doffset<2; doffset++)
    repairdndz(dist[offset+doffset], magdist[doffset], 1., p, 1);

  /* Normalize to constant clustering amplitude instead of constant bias */
  norm0 = DNget_growth_normhiz(p,0.,0);
  for(iz=0; iz<NZ; iz++) {
    z = (iz+0.5)*ZSTEP;
    growth = DNget_growth_normhiz(p,0.,z)/norm0;
    for(doffset=0; doffset<2; doffset++)
      dist[offset+doffset][iz] /= growth;
  }

  /* Save distribution as reference */
  for(doffset=0; doffset<2; doffset++)
    for(iz=0; iz<NZ; iz++)
      refdist[doffset][iz] = dist[offset+doffset][iz];

  for(doffset=0; doffset<2; doffset++) {
    /* Intiialize redshift distribution parameters */
    zp[0] = 1.5; zp[1] = 1.5; zp[2] = 1.5;
    istep=0;

    /* Number of parameters actually used */
    NZP_QSO = doffset? 2: 1;

    chi2_old = 1e99;
    chi2 = 1e49;
    while(chi2<chi2_old-1e-7) {

      /* Initialize search directions */
      if (istep==0)
        for(id=0; id<NZP_QSO; id++)
          for(ip=0; ip<NZP_QSO; ip++)
            dirs[id][ip] = id==ip? 1.: 0.;

      /* Save old location */
      for(ip=0; ip<NZP_QSO; ip++)
        zp_old[ip] = zp[ip];

      chi2_old=chi2;
      c2b=chi2=construct_test_dist_qso(p, dist, dist[offset+doffset], magdist[doffset], refdist[doffset], doffset, zp);
      dchi2=-1;
      dirbest=0;

      for(id=0; id<NZP_QSO; id++) {
        c2_start = c2b;

        /* Make constrained step */
        norm=0;
        for(ip=0; ip<NZP_QSO; ip++) {
          delta[ip] = dirs[id][ip];
          norm += delta[ip]*delta[ip];
        }
        renorm = sqrt(1.0/norm);
        for(ip=0; ip<NZP_QSO; ip++)
          delta[ip]*=renorm;

        /* 1D minimization */
        for(ip=0; ip<NZP_QSO; ip++)
          zp[ip] -= 2*delta[ip];
        for(ix=0; ix<5; ix++) {
          y[ix] = ix==2? c2_start: construct_test_dist_qso(p, dist, dist[offset+doffset], magdist[doffset], refdist[doffset], doffset, zp);
          for(ip=0; ip<NZP_QSO; ip++)
            zp[ip] += delta[ip];
        }
        for(ip=0; ip<NZP_QSO; ip++)
          zp[ip] -= 3*delta[ip];
        xmin = minquartic(y,-1.25,1.25);
        for(ip=0; ip<NZP_QSO; ip++)
          zp[ip] += xmin*delta[ip];
        c2b = construct_test_dist_qso(p, dist, dist[offset+doffset], magdist[doffset], refdist[doffset], doffset, zp);

        c2_end = c2b;
        if (dchi2<c2_start-c2_end) {
          dchi2 = c2_start-c2_end;
          dirbest=id;
        }
      }

      /* Get new chi^2 */
      chi2=c2b;
      /*fprintf(stderr, "%1d %8.5lf %8.5lf  %16.12lf %1d\n", doffset, zp[0], NZP_QSO>1? zp[1]: 0, chi2, dirbest);*/
      istep++;

      /* Fix old directions */
      for(ip=0; ip<NZP_QSO; ip++)
        dirs[dirbest][ip] = zp[ip]-zp_old[ip];

      /* If 1 parameter the minimization is 1D so there's no point in going through again */
      if (NZP_QSO==1) break;
    }

    construct_test_dist_qso(p, dist, dist[offset+doffset], magdist[doffset], refdist[doffset], doffset, zp);

  }

  free_dmatrix(magdist, 0, 1, 0, NZ-1);
  free_dmatrix(refdist, 0, 1, 0, NZ-1);
  free_dmatrix(dirs, 0, NZP_QSOMAX-1, 0, NZP_QSOMAX-1);
}
#undef NZP_QSO

/* Constructs NVSS redshift distribution and returns chi^2.
 * Currently uses a Gamma distribution as we can't do significantly
 * better with the current data.
 *
 * f(z) = e^A (z/C)^B e^(-z/C) where zp=(A,B,C).
 */
double construct_test_dist_nvss(COSMOPARAM *p, double **alldist, double *dist, double *zp) {
  int n, iz, ilrg, iqso;
  double out[3], my_out[3];
  char FileName[1024];
  double z, chi2, amp;

  for(iz=0; iz<NZ; iz++) {
    z = (iz+0.5)*ZSTEP;
    dist[iz] = pow(z/zp[1],zp[0])*exp(-zp[0]*(z/zp[1]-1));
  }

  /* Get the chi^2 */
//  chi2 = getchi2_1p("data/nvss-nvss-trunc.dat", dist, dist, p, out, 0);
  chi2=0;
  for(n=0; n<3; n++) my_out[n] = 0.;

  chi2 += getchi2_1p("data/nvss-lrg0.dat", dist, alldist[4], p, out, 3.398E-07);
  for(n=0; n<3; n++) my_out[n] += out[n];
  chi2 += getchi2_1p("data/nvss-lrg1.dat", dist, alldist[5], p, out, 1.986E-07);
  for(n=0; n<3; n++) my_out[n] += out[n];
  chi2 += getchi2_1p("data/nvss-qso0.dat", dist, alldist[6], p, out, 1.359E-07);
  for(n=0; n<3; n++) my_out[n] += out[n];
  chi2 += getchi2_1p("data/nvss-qso1.dat", dist, alldist[7], p, out, 1.427E-07);
  for(n=0; n<3; n++) my_out[n] += out[n];
  chi2 += getchi2_1p("data/nvss-2mass0.dat", dist, alldist[0], p, out, 9.946E-07);
  for(n=0; n<3; n++) my_out[n] += out[n];
  chi2 += getchi2_1p("data/nvss-2mass1.dat", dist, alldist[1], p, out, 6.057E-07);
  for(n=0; n<3; n++) my_out[n] += out[n];
  chi2 += getchi2_1p("data/nvss-2mass2.dat", dist, alldist[2], p, out, 3.572E-07);
  for(n=0; n<3; n++) my_out[n] += out[n];
  chi2 += getchi2_1p("data/nvss-2mass3.dat", dist, alldist[3], p, out, 2.205E-07);
  for(n=0; n<3; n++) my_out[n] += out[n];

  /* Get amplitude; include +/- 1 sigma if flagged */
  amp = -my_out[1]/my_out[2];
  chi2 = my_out[0] + my_out[1]*amp/2.;
  for(iz=0; iz<NZ; iz++) dist[iz] *= amp;

  return(chi2);
}

#define NZP_NVSS 2
void construct_dist_nvss(COSMOPARAM *p, double **dist) {
  int offset=8;
  int ip, id;
  long istep=0;
  double chi2, norm, renorm, bias;
  double y[5], xmin;
  int ix;
  double zp_old[NZP_NVSS];
  double delta[NZP_NVSS];
  double **dirs;
  double c2b, chi2_old;
  double dchi2, c2_start, c2_end;
  double stepsize=0.05;
  int dirbest;
  double zp[]={1.1, 0.8};

  /* Allocations */
  dirs = dmatrix(0, NZP_NVSS-1, 0, NZP_NVSS-1);

  chi2_old = 1e99;
  chi2 = 1e49;
  while(chi2<chi2_old-1e-7) {

    /* Initialize search directions */
    if (istep==0)
      for(id=0; id<NZP_NVSS; id++)
        for(ip=0; ip<NZP_NVSS; ip++) {
          dirs[id][ip] = id==ip? 1.: 0.;
          dirs[id][ip] = cos(M_PI*id*(ip+0.5)/(double)NZP_NVSS);
        }

    /* Save old location */
    for(ip=0; ip<NZP_NVSS; ip++)
      zp_old[ip] = zp[ip];

    chi2_old=chi2;
    c2b=chi2=construct_test_dist_nvss(p, dist, dist[offset], zp);
    dchi2=-1;
    dirbest=0;

    for(id=0; id<NZP_NVSS; id++) {
      c2_start = c2b;

      /* Make constrained step */
      norm=0;
      for(ip=0; ip<NZP_NVSS; ip++) {
        delta[ip] = dirs[id][ip];
        norm += delta[ip]*delta[ip];
      }
      renorm = sqrt(1.0/norm);
      for(ip=0; ip<NZP_NVSS; ip++)
        delta[ip]*=stepsize*renorm;

      /* 1D minimization -- uses quartic interpolation method, not exact for
       * general functions, but ought to be good enough?
       */
      for(ip=0; ip<NZP_NVSS; ip++)
        zp[ip] -= 2*delta[ip];
      for(ix=0; ix<5; ix++) {
        y[ix] = ix==2? c2_start: construct_test_dist_nvss(p, dist, dist[offset], zp);
        for(ip=0; ip<NZP_NVSS; ip++)
          zp[ip] += delta[ip];
      }
      for(ip=0; ip<NZP_NVSS; ip++)
        zp[ip] -= 3*delta[ip];
      xmin = minquartic(y,-2,2);
      for(ip=0; ip<NZP_NVSS; ip++)
        zp[ip] += xmin*delta[ip];
      c2b = construct_test_dist_nvss(p, dist, dist[offset], zp);

      c2_end = c2b;
      if (dchi2<c2_start-c2_end) {
        dchi2 = c2_start-c2_end;
        dirbest=id;
      }
    }

    /* Get new chi^2 */
    chi2=c2b;
    /*fprintf(stderr, "N %8.5lf %8.5lf  %16.12lf %1d\n", zp[0], zp[1], chi2, dirbest);*/
    istep++;

#if 1
    /* Fix old directions */
    for(ip=0; ip<NZP_NVSS; ip++)
      dirs[dirbest][ip] = zp[ip]-zp_old[ip];
#endif

    /* If step size is not small enough, go back and re-minimize */
#if 0
    if ((!(chi2<chi2_old-1e-7)) && (stepsize>0.01)) {
      stepsize/=10.;
      chi2_old = 1e99;
    }
#endif

    /* If 1 parameter the minimization is 1D so there's no point in going through again */
    if (NZP_NVSS==1) break;
  }

  construct_test_dist_nvss(p, dist, dist[offset], zp);

  free_dmatrix(dirs, 0, NZP_NVSS-1, 0, NZP_NVSS-1);

  /* Get effective bias */
  bias=0;
  for(ix=0;ix<NZ;ix++) bias+=dist[offset][ix]*ZSTEP;
  /*fprintf(stderr, "NVSS effective bias = %9.6lf\n", bias);*/
}
#undef NZP_NVSS


/* utilization:
 *
 * ./dndz_func_scale.x <outfile> <cosmoparams>
 *
 */
int main(int argc, char **argv) {
  COSMOPARAM p;
  double xGamma, beta, z0;
  int nl;
  double **dist;
  int i,iz;
  double z;
  FILE *fp;

  /* Read cosmological parameter file; default=cosmoparam if 3rd arg doesn't exist */
  if (argc>2) {
    DNread_cosmo(&p, argv[2]);
  } else {
    DNread_cosmo(&p, "cosmoparam");
  }
  
  /* Get cosmology, setup nonlinear spectrum */
  p.is_initialized=0;
  p.pk_initialized=0;
  DNnormalize(&p);
  xGamma = p.omega_M * p.h * exp(-p.omega_B - sqrt(2*p.h)*p.omega_B/p.omega_M);
  beta=1.5; z0=1.0;
  nl = 2;
  setparameters_(&(p.omega_M), &(p.omega_de), &xGamma, &(p.sigma8), &(p.ns), &beta, &z0, &nl);

  /* Allocate & null out b*dN/dz */
  dist = (double**)malloc((size_t)(NSAMPLE*sizeof(double*)));
  for(i=0;i<NSAMPLE;i++) {
    dist[i] = (double*)malloc((size_t)(NZ*sizeof(double)));
    for(iz=0;iz<NZ;iz++) dist[i][iz]=0.;
  }
  
  construct_dist_loz(&p, dist);
  construct_dist_qso(&p, dist);
  construct_dist_nvss(&p, dist);
  printalldist(argv[1], dist);
  return(0);
}
