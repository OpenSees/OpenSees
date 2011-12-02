/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * stat.c
 *
 * This file contains code that prints various statistics about the partitioning
 * algorithm at various points.
 *
 * Started 9/15/94
 * George
 *
 * $Id: stat.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include "multilevel.h"

/*************************************************************************
* External variables 
**************************************************************************/
extern CtrlType *__Ctrl;		/* mlevelpart.c */

/*************************************************************************
* This function prints min/max/mean info on the vertex weights and
* vertex edge sums of a graph
**************************************************************************/
void PrintGraphMMM(CoarseGraphType *graph)
{
  int i, j;
  int vwgtmin, vwgtmax, vwgtsum, cewgtsum;
  int esumwgtmin, esumwgtmax, esumwgtsum;
  int nedges, ewgt;
  int tmp;
  VertexType *vtx;
  int maxewgt = 0;
  int ndegrees[1024];
  int maxdeg, mindeg;
  float evratio=0.0;
  float epv=0.0;

  maxdeg = 0;
  mindeg = 20000;

  vtx = graph->vtxs[0];
  vwgtsum = vwgtmin = vwgtmax = vtx->vwgt;
  esumwgtsum = 0;
  for (i=0; i<vtx->nedges; i++)
    esumwgtsum += vtx->edges[i].ewgt;
  esumwgtmin = esumwgtmax = esumwgtsum;

  nedges = ewgt = 0;
  cewgtsum = vtx->cewgt;

  for (i=0; i<1024; i++)
    ndegrees[i] = 0;

  for (i=1; i<graph->nvtxs; i++) {
    vtx = graph->vtxs[i];
    ndegrees[vtx->nedges]++;

    mindeg = (vtx->nedges < mindeg ? vtx->nedges : mindeg);
    maxdeg = (vtx->nedges > maxdeg ? vtx->nedges : maxdeg);

    if (vtx->vwgt < vwgtmin)
      vwgtmin = vtx->vwgt;
    if (vtx->vwgt > vwgtmax)
      vwgtmax = vtx->vwgt;
    vwgtsum += vtx->vwgt;

    tmp = 0;
    for (j=0; j<vtx->nedges; j++)
      tmp += vtx->edges[j].ewgt;
    if (tmp < esumwgtmin)
      esumwgtmin = tmp;
    if (tmp > esumwgtmax)
      esumwgtmax = tmp;
    esumwgtsum += tmp;

    nedges += vtx->nedges;
    ewgt += tmp;
    cewgtsum += vtx->cewgt;

    for (j=0; j<vtx->nedges; j++) {
      if (maxewgt < vtx->edges[j].ewgt)
        maxewgt = vtx->edges[j].ewgt;
    }

    evratio += (.5*vtx->cewgt) / (1.0*vtx->vwgt);

    epv += 1.0*vtx->ewgtsum/(1.0*vtx->vwgt);
  }

  evratio = evratio / (1.0*graph->nvtxs);
  epv = epv/ (1.0*graph->nvtxs);

  printf("Vtcs: %4d %5d %6.1f Edgs: %4d %5d %6.1f AvgEwgt: %6.1f [%6.2f] AvCewgt: %7.1f AvEdgs: %6.1f [%d,%d,%d,%d,%4.2f]\n",
    vwgtmin, vwgtmax, ((float)vwgtsum)/((float)graph->nvtxs),
    esumwgtmin, esumwgtmax, ((float)esumwgtsum)/((float)graph->nvtxs),
    ((float)ewgt)/((float)nedges), epv,
    cewgtsum/((float) graph->nvtxs), ((float)nedges)/((float)graph->nvtxs), mindeg, maxdeg, esumwgtsum, nedges, evratio);
/*
  for (i=0; i<1024; i++) {
    if (ndegrees[i] > 0)
      printf("%d:%d ", i, ndegrees[i]);
  }
  printf("\n");
*/
}



/*************************************************************************
* This function prints information about the partition cuts
**************************************************************************/
void PrintPartResults(char *title, int nparts, int *cuts)
{
  int i, j, k, m;
  int sum=0;
  int n = log2(nparts);
  int flag = 0;

  printf("%s:\n ", title);
  j = 1;
  for (k=1, i=0; i<n; i++, k*=2) {
    printf(" %3d-way",2*k);
  }

  printf("\n  ");
  j = 1;
  for (k=1, i=0; i<n; i++, k*=2) {
    m = j + k;
    for (; j<m; j++) {
      if (cuts[j] == -1)
        flag = 1;
      else
        sum += cuts[j];
    }
    if (!flag)
      printf("%7d ",sum);
    else
      printf("%7d ",-sum);
  }
  printf("\n");

}


/*************************************************************************
* This function prints information about the ordering such as OPC and
* balance figures.
**************************************************************************/
void PrintOrderResults(int nvtxs, SepNodeType *stree, int *smbnz)
{
  int i, j, k, sum, gsum;
  double seropc, paropc;
  int npes;

  /*
   * Compute the opc for each node in stree
   */
  CalcNodeOpc(stree, smbnz, 1);
  seropc = stree[1].opc + stree[1].subopc;

  printf("  Nonzeros: %-12d OPC: %-6.4e ", iasum(nvtxs, smbnz), seropc);

  if (__Ctrl->dbglvl&DBG_ORDERBAL) {
    printf("\n  Balance: ");
    for (npes=2; npes<=256; npes*=2) {
      printf("%3d-PE ",npes);
    }
    printf("\n           ");
    for (npes=2; npes<=256; npes*=2) {
      paropc = CalcParOpc(stree, 1, npes);
      printf("%6.3f ",seropc/(paropc*npes));
    }
  }
  printf("\n");

/* The following prints the # of vertices on the separators */
/*
  for (i=1; i<amin(nvtxs, 32); i++) {
    if (stree[i].nvtxs > 0) {
      printf("%8d (%8d) [%8d-%8d] %12.0f %12.0f\n",i, stree[i].nvtxs, stree[i].li, 
                stree[i].hi, stree[i].opc, stree[i].subopc);
    }
  }

  j = 1;
  k = 1;
  gsum = 0;
  for (i=1; i<8; i++) {
    sum = 0;
    for (;j<=k; j++) {
      if (stree[j].nvtxs == -1) {
        sum = 0;
        break;
      }
      else
        sum += stree[j].nvtxs-stree[2*j].nvtxs-stree[2*j+1].nvtxs;
    }
    if (sum == 0)
      break;
    printf("%4d ", sum);
    k += (1<<i);
    gsum += sum;
  }
  printf(". Total: %5d\n", gsum);
*/
}



/*************************************************************************
* This function computes the opc for each node in the stree
**************************************************************************/
void CalcNodeOpc(SepNodeType *stree, int *smbnz, int root)
{
  int i;

  stree[root].opc = 0.0;
  stree[root].subopc = 0.0;
  if (stree[root].nvtxs == 0)
    return;


  for (i=stree[root].li; i<stree[root].hi; i++) {
    stree[root].opc += smbnz[i] + (smbnz[i]-1)*(smbnz[i]-1);
    /* printf("\n%8d %4d",i, smbnz[i]); */
  }

  if (!stree[root].isleaf) {
    CalcNodeOpc(stree, smbnz, 2*root);
    CalcNodeOpc(stree, smbnz, 2*root+1);

    stree[root].subopc  = stree[2*root].opc + stree[2*root].subopc +
                          stree[2*root+1].opc + stree[2*root+1].subopc;
  }

}

/*************************************************************************
* This function computes the parallel opc for factoring the stree
**************************************************************************/
double CalcParOpc(SepNodeType *stree, int root, int npes)
{
  double opc, opc1, opc2;

  opc = stree[root].opc / (1.0*npes);

  if (npes == 1) {
    opc += stree[root].subopc;
    return opc;
  }

  if (!stree[root].isleaf) {
    opc1 = CalcParOpc(stree, 2*root, npes/2);
    opc2 = CalcParOpc(stree, 2*root+1, npes/2);
    opc += amax(opc1, opc2);
    return opc;
  }

  return opc;
}


/*************************************************************************
* This function computes the balance of the partitions
**************************************************************************/
float ComputePartBalance(CoarseGraphType *graph, int nparts, int *kpwgts)
{
  int maxpwgt;

  maxpwgt = kpwgts[iamax(nparts, kpwgts)];

  return (float)(1.0*maxpwgt*nparts/(1.0*graph->tvwgt));
}
