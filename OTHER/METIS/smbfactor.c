/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * smbfactor.c
 *
 * This file performs the symbolic factorization of a matrix
 *
 * Started 11/1/94
 * George
 *
 * $Id: smbfactor.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include "multilevel.h"


/*************************************************************************
* This function sets up data structures for fill-in computations
**************************************************************************/
void ComputeFillIn(char *filename, int *perm, int n, int maxpool, int *smbnz)
{
  int *xadj, *adjncy;
  int maxlnz;

  printf("  Symbolic factorization...\n");

  ResetPools();
  xadj = imalloc(n+1, "ComputeFillIn: xadj");
  adjncy = (int *)GetEdgePool();

  ReadGraphSKIT(xadj, adjncy, filename);

  smbfactor(xadj, adjncy, n, perm, maxpool, &maxlnz, smbnz);

  free(xadj);
}


/*************************************************************************
* This function performs a symbolic factorization and returns fill-in
* and operation count
**************************************************************************/
void smbfactor(int *xadj, int *adjncy, int m, int *perm, int maxsub, int *maxlnz, int *smbnz)
{
  int i;
  int *iperm;
  int *xlnz, *xnzsub, *nzsub;

  iperm = imalloc(m+1, "smbfactor: iperm");
  xlnz = imalloc(m+1, "smbfactor: xlnz");
  xnzsub = imalloc(m+1, "smbfactor: xnzsub");
  nzsub = ((int *)GetEdgePool()) + xadj[m]+1000; 

  perm--; iperm--;
  for (i=1; i<=m; i++) 
    iperm[perm[i]] = i;
  perm++; iperm++;

  /*
   * Call sparspak routine.
   * Note the vectors perm and iperm have opposite meaning in SPARSPAK and in
   * this program.
   */
  if (smbfct(m, xadj, adjncy, iperm, perm, xlnz, maxlnz, xnzsub, nzsub, &maxsub))
    errexit("MAXSUB too small!");

  for (i=0; i<m; i++)
    smbnz[i] = xlnz[i+1] - xlnz[i] + 1;

  GKfree(xnzsub, xlnz, iperm, -1);

}

/*************************************************************************
* This function computes the elimination tree
**************************************************************************/
void ComputeElTree(int m, int *xlnz, int *nzsub, int *xnzsub, int *parent)
{
  int i;

  /* Calculate the elimination tree */
  xlnz--;
  nzsub--;
  xnzsub--;
  parent--;

  for (i=1; i<=m; i++) {
    if (xlnz[i] < xlnz[i+1])
      parent[i] = nzsub[xnzsub[i]];
    else
      parent[i] = -1;
  }

  xlnz++;
  nzsub++;
  xnzsub++;
  parent++;
}

