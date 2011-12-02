/*
 * lanczos.c
 *
 * This file implements the Lanczos method for fielder vector computations.
 * It is an adaptation from lanczos2.f from Pothen Spectral Nested Disection code.
 *
 * Start 1/17/94
 * George
 *
 * Multilevel Adoption
 * Started 12/8/94
 *
 * $Id: lanczos.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include "multilevel.h"


/**************************************************************************
* This function implements Lanczos eigenvalue solver
* Parameters:
*       graph	The weighted graph to partition, assumed connected
*        evec	The eigenvector
***************************************************************************/
int lanczos(CoarseGraphType *graph, double *evec)   
{
  int ierr, i, ic, iseed, j, jmax, iflag;
  double *u, *r, *alf, *bet, *betc, *z;
  double alfa, beto, be2o, disc, res, maceps, rr, s, dif, du, dl,
         xc, xu, xl, bd, eig, xa, w, xb, da, db, sumi, tiny, tol, oldeps;
  double eval;
  double eps = LANCZOSEPS;
  int lanstp;
  int m = graph->nvtxs;

  /*
   * Memory allocation
   */
  u = (double *)GKmalloc(sizeof(double)*m, "lanczos: u");
  r = (double *)GKmalloc(sizeof(double)*m, "lanczos: r");
  alf = (double *)GKmalloc(sizeof(double)*MAXLANCZOS, "lanczos: alf");
  bet = (double *)GKmalloc(sizeof(double)*MAXLANCZOS, "lanczos: bet");
  betc = (double *)GKmalloc(sizeof(double)*MAXLANCZOS, "lanczos: betc");
  z = (double *)GKmalloc(sizeof(double)*MAXLANCZOS, "lanczos: z");

  ierr = 0;
  iflag = 0;
  maceps = 1.0e-12;
  tiny = maceps*maceps*maceps;
  eig = 0.0;


  /*
   * Initialize Lanczos. u is the starting vector
   */
  iseed = 1;
  bet[0] = 0.0e0;
  betc[0] = 0.0e0;
  for (i=0; i<m; i++) {
    r[i] = 0.0e0;
    u[i] = 0.5e0*((double) (-m-1+2*(i+1)));
  }
  res = dnrm2(m, u);
  res = 1.0e0/res;
  if (eps > res) {
    oldeps = eps;
    eps = res/2.0;
  }

  dscal(m, res, u);
  dcopy(m, u, evec);

  /*
   * Lanczos loop
   */
  lanstp = 0;
  sumi = 0.0e0;
  for (j=0; j<MAXLANCZOS-1; j++) {
    dscal(m, -bet[j], r);
    matvec(graph, u, r);
    alf[j] = ddot(m, u, r);
    daxpy(m, -alf[j], u, r);

    /*
     * Orthogonalize against 1, 1, ..
     */
    bet[j+1] = dnrm2(m, r);
    if (bet[j+1] <= maceps) 
      iflag = 1;
    else
      dscal(m, 1.0e0/bet[j+1], r);
    dswap(m, u, r);

    /*
     * End of Lanczos recurrence
     * u contains q_{j+1}, r contains q_{j}
     * Initialize the eigen analysis
     */
    alfa = alf[j];
    beto = bet[j];
    be2o = beto*beto;
    betc[j] = be2o;
    tol = 0.5*eps*fabs(eig);

    if (j == 0) 
      goto L4500;

    if (j == 1) {
      dif = (alf[0]-alfa)*0.5;
      xc = alfa+dif + sqrt(dif*dif + be2o);
      xu = xc + eps*fabs(xc); 
      du = alfa - xu - be2o/(alf[0]-xu);
      xl = 2.0*xc - xu; 
      eig = xc;
    }
    else {  /* j >= 2 */
      du = alfa - xu - be2o/du;
      ic = 0;
      if (du > 0.0) { 
        xl = xu;
        dl = du;
        xa = xl + bd;
        pivot(alf, &betc[1], j, xa, &da);
        ic++;
        eig = xc;
        xb = xl + 0.5*(xa - xl);
        xu = beto + amax(xl, alfa);

        for (i=0; i<50; i++) {
          pivot(alf, &betc[1], j, xb, &db);
          ic++;
          if (db > 0.0)
            xl = xb;
          else {
            xu = xb;
            du = db;
          }
          w = xb - xa;
          if (w <= 0.0e0)
            xc = xb;
          else {
            s = -w - (db*(xb-eig)-da*(xa-eig))/w;
            rr = db*(xb-eig);
            disc = s*s - 4.0*rr;
            if (disc <= tiny)
              xc = xb;
            else
              xc = xb + 2.0*rr/(s+sign(sqrt(disc),s));
          }
          if ((xc < (xl+0.5*tol)) || (xc > (xu-0.5*tol)))
            xc = xl+0.5*(xu-xl);
          xa = xb;
          da = db;
          xb = xc;

          if (xu-xl <= eps)
            break;
        } /* for */
        if (i==50) {
          printf("\nMore than 50 iterations were performed with no convergance");
          return -2;
        }
      } /* if (du > 0.0) */
    } /* else (j >= 2) */

    if (du == 0.0)
      du = eps;

    /*
     * use givens recursion to get eigenvector and bound
     */
    givens(alf, bet, j, xc, z);
    res = fabs(z[j])*bet[j+1];
    bd = res;
    sumi = sumi + (float) ic;

L4500:
    jmax = j;
    if (res < eps || iflag == 1)
      break;  /* Done with Lanczos iterations */
  } /* for (j=0; j<MAXLANCZOS-1; j++) */
  if (j == MAXLANCZOS-1) {
    ierr = 0;
  }

  eig = xc;
  sumi = sumi/(float) (1+jmax);
  if (jmax == 0)
    z[0] = 1.0;


  /*
   * Compute eigenvector after convergence
   */
  dcopy(m, evec, u);
  for (i=0; i<m; i++) {
    evec[i] = z[0]*u[i];
    r[i] = 0.0;
  }
  eval = xc;
  bet[0] = 0.0;
  betc[0] = 0.0;

  /*
   * run Lanczos loop again for computation of vector
   */
  for (j=0; j<jmax-1; j++) {
    dscal(m, -bet[j], r);
    matvec(graph, u, r);

    daxpy(m, -alf[j], u, r);

    dscal(m, 1.0/bet[j+1], r);
    dswap(m, u, r);
    daxpy(m, z[j+1], u, evec);
  }

  if (eps < oldeps)
    eps = oldeps;
  lanstp = jmax;

  GKfree(u, r, alf, bet, betc, z, -1);

  return -ierr;
}



/*************************************************************************
* This uses Given's recurence for computing the eigenvector of a 
* tridiagonal matrix, when the eigenvalue is known
**************************************************************************/
void givens(double *a, double *b, int n, double eig, double *evec)
{
  int i;
  double sum, bd;

  evec[n] = 1.0;
  evec[n-1] = (eig - a[n])*evec[n]/b[n];
  sum = 1.0 + evec[n-1]*evec[n-1];
  for (i=n-2; i>=0; i--) {
    evec[i] = ((eig - a[i+1]) * evec[i+1] - b[i+2]*evec[i+2])/b[i+1];
    sum = sum + evec[i]*evec[i];
  }
  bd = 1.0/sqrt(sum);

  dscal(n+1, bd, evec);

}

/*************************************************************************
* This function performs pivoting???
**************************************************************************/
void pivot(double *a, double *b, int n, double x, double *d)
{
  int i;

  *d = 1.0;
  if (n < 0)
    return;

  *d = a[0] - x;
  if (*d == 0.0)
    *d = b[0]*1.0e-6;
  if (n == 0)
    return;

  for (i=1; i<=n; i++) {
    *d = a[i] - x - b[i-1] / *d;
    if (*d == 0.0)
      *d = b[i-1]*1.0e-6;
  }

}


/*************************************************************************
* This function multiplies a matrix by a vector
* b = b + Aa  where A = adj - diag
**************************************************************************/
void matvec(CoarseGraphType *graph, double *a, double *b)
{
  VertexType *vtx;
  int i, j;
  EdgeType *edges;

  for (i=0; i<graph->nvtxs; i++) {
    vtx = graph->vtxs[i];
    edges = vtx->edges;

    b[i] -= a[i]*vtx->ewgtsum;
    for (j=0; j<vtx->nedges; j++) 
      b[i] += edges[j].ewgt * a[edges[j].edge];
  }

}


/*************************************************************************
* This function prints a vector
**************************************************************************/
void printvector(int n, double *x)
{
  int i;

  for (i=0; i<n; i++)
    printf("%3.2e ",x[i]);
  printf("\n");

}

