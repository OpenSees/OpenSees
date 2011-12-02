/* PD1.f -- translated by f2c (version 20020621).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/


/* Table of constant values */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int c0 = 0;
static int c1 = 1;
static int c9 = 9;
static int c5 = 5;

int reloading_(double *, double *, 
	       double *, double *, double *, double *, 
	       double *, double *, double *, double *, int *);

int  unloading_(double *, double *, double *, double 
		*, double *, double *, double *, int *);

int yield1_(int *, double *, double *,
	    double *, double *);

int vdtan1_(int *, int *, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, 
	    double *, double *);

int crstr1_(double *, double *, 
	    double *, double *, double *, double *);
  
int degrad1_(int *, int *, double *, 
	     double *, double *, double *, double *, 
	     double *, double *);

int plasto1_(double *, double *, 
	     int *, double *, double *, double *, double *,
	     double *, double *, double *, double *, 
	     double *, double *, double *, double *, 
	     double *, int *);

int setpara_(double *, double *);

int algotan1_(double *, double *, double *, double *, 
	      double *, double *, double *, double *);

int elastan1_(double *, double *, double *, double *);


int pd_(double *d__, double *hstvp, double *hstv,
	double *epsp, double *sigp, double *deps, double *
	str, double *dd, int *ist)
{
  /* System generated locals */
  double d1, d2;
  
  /* Builtin functions */
  double sqrt(double);
  
  /* Local variables */
  static double phibound, e;
  static double ck, kp, mu, xl[4]	/* was [2][2] */, deg, lam, fkp, eps, tol, tol2;
  static int flag__;
  static double kapa[2], vdeg, cohn[2], delt, cmax, resf, ieps, sign, temp, peps, estr, fstr, dplas;
  static int index;
  static double delps, d2_eps__, trstr, dplas1;
  static double chleng;
  static int crmode;
  static double deleps, oldeps, degstr, ktcrit, viscom;
  static int maxitr;
  static double viscot;
  static double matpara[4], fenergy;
  
  /*       One-D version of (ElastoVisco) Plastic-Damage Model */
  
  /*    % Modified for OpenSees: June 2006 */
  /*    % Original Coded: July 1998 */
  /*    % By Jeeho Lee (jeeholee@berkeley.edu) */
  
  /*    # Input:  d[*],delt,eps,stif_typ */
  /*    # Output: str,dd */
  /*    # History Variables (size=11): hstvp, hstv */
  /*     ->ieps,peps,kapa[2],cohn[2],deg,vdeg,dplas,oldeps,phibound */
  
  /*    ! ist = 1: dd is tangent stiffness in any case */
  
  /*  Local Variables ---------------------------------------------- */
  /* set numerical parameters ------------------------------------------ */
  /* Parameter adjustments */
  --hstv;
  --hstvp;
  --d__;
  
  /* Function Body */
  tol = 1e-11;
  tol2 = 1e-7;
  maxitr = 10;
  /* temperary setting ++++++++++----------------------- */
  eps = *epsp + *deps;
  xl[0] = 0.f;
  xl[2] = 1.f;
  xl[1] = 0.f;
  xl[3] = 0.f;
  delt = 1.;
  /* initialize the material properties & variables -------------------- */
  e = d__[1];
  ktcrit = d__[7];
  mu = d__[8] * delt;
  /*  retrieve history variables --------------------------------------- */
  ieps = hstvp[1];
  peps = hstvp[2];
  kapa[0] = hstvp[3];
  kapa[1] = hstvp[4];
  cohn[0] = hstvp[5];
  cohn[1] = hstvp[6];
  deg = hstvp[7];
  vdeg = hstvp[8];
  dplas = hstvp[9];
  oldeps = hstvp[10];
  phibound = hstvp[11];
  /*  ------------------------------------------------------------------ */
  crmode = 0;
  /*      if (delt .le. 0.d0) then */
  delt = 1.;
  /*      endif */
  viscom = mu / (mu + delt);
  viscot = delt / (mu + delt);
  if (kapa[0] + kapa[1] <= 0.) {
    cohn[0] = d__[4];
    cohn[1] = -d__[5];
  }
  /*  compute trial effective stresses --------------------------------- */
  deleps = eps - peps;
  trstr = e * deleps;
  /*  check whether tensile or compressive state */
  if (trstr >= 0.) {
    sign = 1.;
    index = 1;
    kp = kapa[0];
    dplas1 = 1. - dplas;
  } else {
    sign = -1.;
    index = 2;
    kp = kapa[1];
    dplas1 = 1.;
  }
  flag__ = 0;
  setpara_(&d__[1], matpara);
  yield1_(&index, cohn, &trstr, &resf, &temp);
  cmax = sqrt(d__[4] * d__[4] + d__[5] * d__[5]);
  if (resf < tol * cmax) {
    if (kp > 0. && index == 2) {
      *deps = eps - oldeps;
      if (*deps < 0.) {
	/* +++++++++++++++++++++++++ Reloading +++++++++++++++++++++++++++++++ */
	/* Computing 2nd power */
	d1 = xl[0] - xl[2];
	/* Computing 2nd power */
	d2 = xl[1] - xl[3];
	chleng = sqrt(d1 * d1 + d2 * d2);
	reloading_(&chleng, &kp, &d__[1], matpara, &eps, deps, &peps, 
		   &phibound, cohn, &tol, &maxitr);
	kapa[1] = kp;
      } else {
	unloading_(&d__[1], &eps, deps, &kp, cohn, &peps, &tol, &
		   maxitr);
	phibound = resf / cohn[1] + 1.;
      }
    }
    degrad1_(&c0, &index, &d__[1], matpara, kapa, &kapa[1], &ck, &deg, &
	     degstr);
    ieps = viscot * peps + viscom * ieps;
    deleps = eps - ieps;
    vdeg = viscom * vdeg + viscot * deg;
    elastan1_(&e, &vdeg, &dplas1, dd);
    *str = (1. - vdeg) * dplas1 * e * deleps;
  } else {
    /*  plastic loading state  --------------------------------------------- */
    /* Computing 2nd power */
    d1 = xl[0] - xl[2];
    /* Computing 2nd power */
    d2 = xl[1] - xl[3];
    chleng = sqrt(d1 * d1 + d2 * d2);
    if (kp > ktcrit && index == 1) {
      delps = 0.;
      crstr1_(&e, cohn, &trstr, &dplas, &dplas1, &d2_eps__);
      crmode = 1;
      /*          write(*,*) '!!! crmode ON' */
    } else {
      plasto1_(&d__[1], matpara, &index, &sign, &chleng, &eps, &trstr, &
	       lam, &kp, cohn, &fenergy, &fstr, &fkp, &ck, &dplas1, &tol,
	       &maxitr);
      delps = lam * sign;
      algotan1_(&e, &sign, &lam, &fenergy, &fstr, &fkp, &ck, dd);
    }
    peps += delps;
    estr = e * (eps - peps);
    if (index == 1) {
      kapa[0] = kp;
    } else {
      kapa[1] = kp;
    }
    degrad1_(&c1, &index, &d__[1], matpara, kapa, &kapa[1], &ck, &deg, &
	     degstr);
    trstr = dplas1 * e * (eps - ieps);
    ieps = viscot * peps + viscom * ieps;
    deleps = eps - ieps;
    vdeg = viscom * vdeg + viscot * deg;
    *str = (1. - vdeg) * dplas1 * e * deleps;
    vdtan1_(&crmode, &index, &d__[1], &dplas1, &mu, &delt, &trstr, &estr, 
	    &vdeg, &degstr, &d2_eps__, dd);
  }
  /*  Compute secant tangent stiffness (only for elastic part) */
  /*  Post-processing part */
  /*  save history variables ------------------------------------ */
  hstv[1] = ieps;
  hstv[2] = peps;
  hstv[3] = kapa[0];
  hstv[4] = kapa[1];
  hstv[5] = cohn[0];
  hstv[6] = cohn[1];
  hstv[7] = deg;
  hstv[8] = vdeg;
  hstv[9] = dplas;
  hstv[10] = eps;
  hstv[11] = phibound;
  /*  ------------------------------------------------------------ */
  return 0;
} /* pd_ */

 int setpara_(double *d__, double *matpara)
{
  /* System generated locals */
  double d1, d2;
  
  
  /* Local variables */
  static double deg_para1__, deg_para2__, deg_para3__;
  static double temp;
  
  /* Parameter adjustments */
  --matpara;
  --d__;
  
  /* Function Body */
  deg_para1__ = .7f;
  deg_para2__ = .5f;
  deg_para3__ = 1.f;
  /* ----- compute the parameters for the material model */
  temp = d__[6] / d__[5];
  matpara[1] = temp * 2. - 1. + sqrt(temp * temp - temp) * 2.;
  d1 = 1. - deg_para1__;
  d2 = 1. / deg_para3__;
  temp = pow(d1, d2);
  matpara[2] = deg_para3__;
  matpara[3] = log(1. - deg_para2__) / log(.5 / matpara[1] + .5);
  matpara[4] = (temp - .5) / (temp * temp - temp);
  return 0;
} /* setpara_ */

 int crstr1_(double *e, double *cohn, double *
			     trstr, double *dplas, double *dplas1, double *d2_eps__)
{
  /* Parameter adjustments */
  --cohn;
  
  /* Function Body */
  *dplas1 = *dplas1 * cohn[1] / *trstr;
  *dplas = 1. - *dplas1;
  *d2_eps__ = -(*e) * *dplas / *trstr;
  return 0;
} /* crstr1_ */



int damg1_(int *, int *, double *, 
	   double *, double *, double *, double *, 
	   double *, double *);
int coml1_(int *, double *, 
	   double *, double *, double *, double *, 
	   double *);

int plasto1_(double *d__, double *matpara, int *
	     index, double *sign, double *chleng, double *eps, 
	     double *trstr, double *lam, double *kp, double *cohn, 
	     double *fenergy, double *fstr, double *fkp, double *
	     ck, double *dplas1, double *toler, int *maxitr)
{
  /* Builtin functions */
  
  /* Local variables */
  static double e, kpn;
  static int iter;
  static double temp, resq;
  
  static double q_lam__, lamkp, error, q_kapa__;
  static int switch__;
  
   /* Fortran I/O blocks */
  
  /* plastic or viscoplastic loading case --------------------------------- */
  /* Parameter adjustments */
  --cohn;
    --matpara;
    --d__;
    
    /* Function Body */
    e = *dplas1 * d__[1];
    /* initial setting -------------------------------------------------- */
    iter = 0;
    switch__ = 1;
    *lam = 0.;
    if (*index == 1) {
      *fenergy = d__[2] / *chleng;
    } else {
      *fenergy = d__[3] / *chleng;
    }
    kpn = *kp;
    /* iteration for damage evolution ------------------------------------ */
 L100:
    ++iter;
    damg1_(&c1, index, &d__[1], &matpara[1], kp, &cohn[1], fstr, fkp, ck);
    coml1_(index, &e, trstr, &cohn[1], ck, lam, &lamkp);
    /*   check the convergence of the damage evolution eqn */
    resq = kpn - *kp + *sign * *lam * *fstr / *fenergy;
    error = fabs(resq);
    /*      write(*,*) '#################################' */
    /*      write(*,*) 'toler =',toler */
    /*      write(*,*) 'error =', error */
    /*      write(*,*) 'kp =',kp */
    if (error > *toler) {
      if (iter > *maxitr) {
	fprintf(stderr, "toler = %f\n error = %f\n kp = %f\n",*toler, error, *kp);
	fprintf(stderr, "VEPD_@D: exceed the maximum iteration(iter)!\n");
	exit(-1);
      }
      q_lam__ = *sign * *fstr / *fenergy;
      q_kapa__ = *sign * *lam * *fkp / *fenergy - 1.;
      q_kapa__ += q_lam__ * lamkp;
      resq /= q_kapa__;
      *kp -= resq;
      temp = 1. - *toler;
      if (*kp < kpn) {
	*kp = kpn;
      } else if (*kp > temp) {
	*kp = temp;
	switch__ = -1;
      }
      goto L100;
      /* ---> go back to 100 */
    }
    /*     end of iteration loop ---------------------------------------- */
    return 0;
} /* plasto1_ */

int algotan1_(double *e, double *sign, double *
	      lam, double *fenergy, double *fstr, double *fkp, 
	      double *ck, double *dd)
{
  static double omega;
  
  /*   construct the algorithmic tangent stiffness: dd */
  omega = *fstr * *ck / (*fenergy - *sign * *lam * *fkp);
  *dd = *e * omega / (omega - *e);
  return 0;
} /* algotan1_ */

int vdtan1_(int *crmode, int *index, double *d__,
	    double *dplas1, double *mu, double *delt, double *
	    tstr, double *estr, double *vdeg, double *degstr, 
	    double *d2_eps__, double *dd)
{
  /* Local variables */
  static double e, temp, vdeg1, temp2;
  
  /* Parameter adjustments */
  --d__;
  
  /* Function Body */
  e = d__[1];
  vdeg1 = 1. - *vdeg;
  /* --- Modify Inviscid Algorithmic Tangent for Large Crack Mode --- */
  if (*crmode >= 1) {
    *dd = vdeg1 * (e * *dplas1 - *estr * *d2_eps__);
  } else {
    temp = *mu / (*mu + *delt);
    temp2 = *delt / (*mu + *delt);
    *tstr = temp * *tstr + pow(temp2, *dplas1) * *estr;
    *dd = temp * vdeg1 * *dplas1 * e + temp2 * (vdeg1 - *tstr * *degstr) *
      *dd;
  }
  return 0;
} /* vdtan1_ */

int elastan1_(double *e, double *vdeg, double *
	      dplas1, double *dd)
{
  /*    compute the elastic tangent stiffness */
  *dd = (1. - *vdeg) * *dplas1 * *e;
  return 0;
} /* elastan1_ */

int coml1_(int *index, double *e, double *trstr, 
	   double *cohn, double *ck, double *lam, double *lamkp)
{
  /* Parameter adjustments */
  --cohn;
  
  /* Function Body */
  if (*index == 1) {
    *lam = (*trstr - cohn[1]) / *e;
    *lamkp = -(*ck) / *e;
  } else {
    *lam = -(*trstr + cohn[2]) / *e;
    *lamkp = -(*ck) / *e;
  }
  return 0;
} /* coml1_ */

int yield1_(int *index, double *cohn, double *
	    str, double *resf, double *fb)
{
  /* Parameter adjustments */
  --cohn;
  
  /* Function Body */
  if (*index == 1) {
    *fb = cohn[2] * *str / cohn[1];
  } else {
    *fb = -(*str);
  }
  *resf = *fb - cohn[2];
  return 0;
} /* yield1_ */

int damg1_(int *flag__, int *index, double *d__, 
	   double *matpara, double *kapa, double *cohn, double *
	   fstr, double *fkp, double *ck)
{
  /* System generated locals */
  double d1, d2;
  
  /* Local variables */
  static double ac, at, rc, cy, rt, ty, phc, pht, rphc, temp, rpht;
  
  /* Parameter adjustments */
  --cohn;
  --matpara;
  --d__;
  
  /* Function Body */
  ac = matpara[1];
  rt = matpara[2];
  rc = matpara[3];
  at = matpara[4];
  if (*index == 1) {
    /* --- Tensile Damage */
    ty = d__[4];
    pht = at * (at + 2.) * *kapa + 1.;
    rpht = sqrt(pht);
    *fstr = ty * ((at + 1.) * sqrt(pht) - pht) / at;
    /* ----------- compute cohesions ----------------------------- */
    d1 = (at + 1. - rpht) / at;
    d2 = 1. - rt;
    cohn[1] = ty * rpht * pow(d1, d2);
    /* ----------- compute derivatives --------------------------- */
    *fkp = ty * (at + 2.) * ((at + 1.) / (sqrt(pht) * 2.) - 1.);
    temp = (at + 1. - rpht) / at;
    d1 = 1. - rt;
    d2 = -rt;
    *ck = ty * .5 * at * (at + 2.) * (pow(temp, d1) / rpht - (1. - 
								   rt) * pow(temp, d2) / at);
    /* --- Compressive Damage */
  } else {
    cy = d__[5];
    phc = ac * (ac + 2.) * *kapa + 1.;
    rphc = sqrt(phc);
    *fstr = cy * ((ac + 1.) * sqrt(phc) - phc) / ac;
    /* ----------- compute cohesions ----------------------------- */
    d1 = (ac + 1. - rphc) / ac;
    d2 = 1. - rc;
    cohn[2] = -cy * rphc * pow(d1, d2);
    /* ----------- compute derivatives --------------------------- */
    *fkp = cy * (ac + 2.) * ((ac + 1.) / (sqrt(phc) * 2.) - 1.);
    temp = (ac + 1. - rphc) / ac;
    d1 = 1. - rc;
    d2 = -rc;
    *ck = cy * -.5 * ac * (ac + 2.) * (pow(temp, d1) / rphc - (1. 
								    - rc) * pow(temp, d2) / ac);
  }
  return 0;
} /* damg1_ */

int degrad1_(int *flag__, int *index, double *
	     d__, double *matpara, double *kapa1, double *kapa2, 
	     double *ck, double *deg, double *degstr)
{
  /* System generated locals */
  double d1;
  
  /* Local variables */
  static double s, ac, dc, rc, at, dt, rt, phi, phi2, dfac, temp, temp2,
    degkp;
  
  /* Parameter adjustments */
  --matpara;
  --d__;
  
  /* Function Body */
  ac = matpara[1];
  rt = matpara[2];
  rc = matpara[3];
  at = matpara[4];
  dfac = 1.;
  phi = sqrt(at * (at + 2.) * *kapa1 + 1.);
  temp = (at + 1. - phi) / at;
  dt = 1. - pow(temp, rt);
  phi2 = sqrt(ac * (ac + 2.) * *kapa2 + 1.);
  temp2 = (ac + 1. - phi2) / ac;
  dc = 1. - pow(temp2, rc);
  if (*index == 1) {
    s = 1.;
  } else {
    s = 1. - dfac;
  }
  *deg = 1. - (1. - s * dt) * (1. - dc);
  if (*flag__ >= 1) {
    if (*index == 1) {
      d1 = rt - 1.;
      degkp = s * .5 * (at + 2.) * (1. - dc) * rt * pow(temp, d1)
	/ phi;
      degkp /= *ck;
      *degstr = (1. - dc) * degkp;
    } else {
      d1 = rc - 1.;
      degkp = (ac + 2.) * .5 * (1. - s * dt) * rc * pow(temp2, d1) / phi2;
      degkp = -degkp / *ck;
      *degstr = (1. - s * dt) * degkp;
    }
  }
  return 0;
} /* degrad1_ */

int reloading_(double *chleng, double *kp, 
	       double *d__, double *matpara, double *eps, double *
	       deps, double *peps, double *phib, double *cohn, 
	       double *toler, int *maxitr)
{
  /* Local variables */
  static double e, ck, ek, ar, br, cr, fkp, phi, kpn, resf;
  static int iter;
  static double temp, resq, estr, fstr;
  extern  int damg1_(int *, int *, double *, 
		     double *, double *, double *, double *, 
		     double *, double *);
  static int index;
  static double pepsn, error;
  extern  int yield1_(int *, double *, double *,
		      double *, double *);
  static double q_kapa__;
  static int switch__;
  static double fenergy;
  
  /* plastic or viscoplastic loading case --------------------------------- */
  /* Parameter adjustments */
  --cohn;
  --matpara;
  --d__;
  
  /* Function Body */
  e = d__[1];
  ar = 0.;
  br = 1.;
  cr = 1.;
  if (br > 1.) {
    br = 1.;
  }
  /* initial setting -------------------------------------------------- */
  index = 2;
  iter = 0;
  switch__ = 1;
  fenergy = d__[3] / *chleng;
  fenergy /= cr;
  pepsn = *peps;
  kpn = *kp;
  /* iteration for damage evolution ------------------------------------ */
 L100:
  ++iter;
  estr = e * (*eps - *peps);
  damg1_(&c1, &index, &d__[1], &matpara[1], kp, &cohn[1], &fstr, &fkp, &
	 ck);
  yield1_(&index, &cohn[1], &estr, &resf, &temp);
  /*   check the convergence of the damage evolution eqn */
  phi = resf / cohn[2] + 1.;
  if (phi < 0.) {
    fprintf(stderr,"RELOADING: Negative phi!");
    exit(-1);
  }
  resq = -(*peps) + pepsn + ar * (phi - *phib) * *deps;
  error = fabs(resq);
  if (error > *toler) {
    if (iter > *maxitr) {
      fprintf(stderr,"toler = %f\n error = %f\n kp = %f\n", *toler, error, *kp);
      fprintf(stderr,"RELOADING: exceed the maximum iteration (iter)!\n");
      exit(-1);
    }
    ek = (fenergy / (1. - br) + (*peps - pepsn) * fkp) / fstr;
    q_kapa__ = -ek + ar * *deps * (e * ek + estr * ck / cohn[2]) / cohn[2]
      ;
    resq /= q_kapa__;
    *kp -= resq;
    temp = 1. - *toler;
    if (*kp < kpn) {
      *kp = kpn;
    } else if (*kp > temp) {
      *kp = temp;
      switch__ = -1;
    }
    *peps = pepsn + (*kp - kpn) * fenergy / (fstr * (1. - br));
    goto L100;
    /* ---> go back to 100 */
  }
  /*     end of iteration loop ---------------------------------------- */
  return 0;
} /* reloading_ */

int yield1_(int *, double *, double *, double *, double *);

int unloading_(double *d__, double *eps, double *
	       deps, double *kp, double *cohn, double *peps, double *
	       toler, int *maxitr)
{
  /* Local variables */
  static double e, ar, br, resf;
  static int iter;
  static double temp, resq, estr;
  static int index;
  static double pepsn, error;
  
  static double q_epsr__;
  static int switch__;
  
  /* plastic or viscoplastic loading case --------------------------------- */
  /* Parameter adjustments */
  --cohn;
  --d__;
  
  /* Function Body */
  e = d__[1];
  ar = 0.;
  br = 1.;
  /* initial setting --------------------------------------------------- */
  iter = 0;
  switch__ = 1;
  index = 2;
  pepsn = *peps;
  /* computation for damage evolution ------------------------------------ */
  br = br;
 L100:
  ++iter;
  estr = e * (*eps - *peps);
  yield1_(&index, &cohn[1], &estr, &resf, &temp);
  /*   check the convergence of the damage evolution eqn */
  resq = -(*peps) + pepsn - ar * *deps * resf / cohn[2];
  error = fabs(resq);
  if (error > *toler) {
    if (iter > *maxitr) {
      fprintf(stderr, "toler = %f\nerror = %f\n",*toler, error);
      fprintf(stderr, "UNLOADING: exceed the maximum iteration (iter)!\n");
      exit(-1);
    }
    q_epsr__ = -1. - ar * e * *deps / cohn[2];
    resq /= q_epsr__;
    *peps -= resq;
    goto L100;
    /* ---> go back to 100 */
  }
  return 0;
} /* unloading_ */

