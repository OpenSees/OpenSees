/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * initpart.c
 *
 * This file contains code that implements the BFS based initial partition
 * algorithm
 *
 * Started 11/22/94
 * George
 *
 * $Id: initpart.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 */

#include "multilevel.h"

/*************************************************************************
* External Variables
**************************************************************************/
extern CtrlType *__Ctrl;        /* mlevelpart.c */
#ifndef METISLIB
extern timer InitPartTmr;	/* main.c */
#endif


/*************************************************************************
* This function randomly selects a node among the large vertex-weight nodes
* and grows a region around it, by doing a gready BFS that minimizes the 
* edges crossed.
**************************************************************************/
void InitPartition(CoarseGraphType *graph, int zeropwgt)
{

  starttimer(&InitPartTmr);

  switch (__Ctrl->InitPartType) {
    case INITPART_GGP:
      GGPPartition(graph, zeropwgt, 0);
      break;
    case INITPART_GGGP:
      if (graph->nvtxs < 2000)
        GGGPPartition(graph, zeropwgt); 
      else
        GGPPartition(graph, zeropwgt, 0); 
      break;
    case INITPART_EIG:
      EigPartition(graph, zeropwgt);
      break;
    case INITPART_GGPKL:
      GGPPartition(graph, zeropwgt, 1);
      break;
    default:
      errexit("Unsupported Initial Partition Type: %d", __Ctrl->InitPartType);
  }
  ComputePartitionParams(graph);

  if (__Ctrl->dbglvl&DBG_INITCUT) 
    printf("[%5d, %5d, %5d, %5d]\n", graph->pwgts[0], graph->pwgts[1], zeropwgt, graph->mincut);

  stoptimer(&InitPartTmr);

}



/*************************************************************************
* This function randomly selects a node among the large vertex-weight nodes
* and grows a region around it, by doing a gready BFS that minimizes the 
* edges crossed.
**************************************************************************/
void GGPPartition(CoarseGraphType *graph, int zeropwgt, int usekl)
{
  int ii, j, k, ntrials;
  int *touched;
  int higain;
  int imaxvwgt;
  int partwgt, partcut;
  VertexType *vtx;
  int *where, *bestwhere, *tmpptr;
  int bestcut, bestwgt;
  int *queue;
  int qhead, qtail;
  int orgvwgt, growp, shrinkp;

  ntrials = amin((usekl ? GGPKL_NTRIALS : GGP_NTRIALS), graph->nvtxs);

  if (graph->tvwgt/graph->nvtxs < 6)
    ntrials = amax(1, ntrials/2);

  if (__Ctrl->dbglvl&DBG_INITCUT) {
    if (!usekl)
      printf("GGP-%d (%5d %5d) ", ntrials, graph->tvwgt, zeropwgt);
    else
      printf("GGPKL-%d (%5d %5d) ", ntrials, graph->tvwgt, zeropwgt);
    fflush(stdout);
  }

  bestwhere = imalloc(graph->nvtxs, "BFSPartition: bestwhere");
  where = imalloc(graph->nvtxs, "BFSPartition: graph->where");
  touched = icoremalloc(graph->nvtxs, "BFSPartition: touched", 0);
  queue = icoremalloc(graph->nvtxs, "TrueBFSPartition: queue", 0);

  growp = 1; shrinkp = 0;
  orgvwgt = graph->tvwgt - zeropwgt;
  bestcut = orgvwgt*2000;

  imaxvwgt = RandomInRange(graph->nvtxs);
  for (ii=0; ii<ntrials; ii++) {
    iset(graph->nvtxs, shrinkp, where);
    iset(graph->nvtxs, 0, touched);

    queue[0] = imaxvwgt;
    qhead = 0;
    qtail = 1;

    partwgt = 0;
    touched[imaxvwgt] = 1;

    while (partwgt <= orgvwgt) {
      if (qtail == qhead) { /* Disconnected graph */
        for (;;) {
          higain = RandomInRange(graph->nvtxs);
          if (where[higain] == shrinkp)
            break;
        }
      }
      else
        higain = queue[qhead++];

      vtx = graph->vtxs[higain];

      if (orgvwgt < partwgt + vtx->vwgt)
        break;

      partwgt += vtx->vwgt;
      where[higain] = growp;

      for (j=0; j<vtx->nedges; j++) {
        k = vtx->edges[j].edge;
        if (touched[k] == 0) {  /* Look only the untouched nodes */
          queue[qtail++] = k;
          touched[k] = 1;
        }
      }
    }

    graph->where = where;
    if (usekl) {
      ComputePartitionParams(graph);
      if (!__Ctrl->IsWeighted)
        BFMR_Refine(graph, zeropwgt, SMART, 3);
      else
        BFMR_Refine_Weighted(graph, zeropwgt, SMART, 3);
      GKfree(graph->id, graph->ed, graph->htable.ht, -1);
      partcut = graph->mincut;
    }
    else
      partcut = ComputeCut(graph);

    if (__Ctrl->dbglvl&DBG_INITCUT) {
      printf("(%4d, %5d %5d) ",queue[0], partcut, graph->pwgts[0]);
      fflush(stdout);
    }

    if (partcut < bestcut) {
      SWAP(bestwhere, where, tmpptr);
      bestcut = partcut;
      bestwgt = partwgt;
    }

    if (bestcut == 0)
      break;

    tmpptr = (ii == 0 ? bestwhere : where);
    for (;;) {
      imaxvwgt = RandomInRange(graph->nvtxs);
      if (tmpptr[imaxvwgt] == shrinkp)
        break;
    }
  }

  graph->mincut = bestcut;
  graph->where = bestwhere;

  free(where);
  icorefree(2*graph->nvtxs);
}






/*************************************************************************
* This function randomly selects a node among the large vertex-weight nodes
* and grows a region around it, by doing a gready BFS that minimizes the 
* edges crossed.
**************************************************************************/
void GGGPPartition(CoarseGraphType *graph, int zeropwgt)
{
  int ii, i, j, k, ntrials;
  int *touched;
  int oldgain;
  int higain;
  int *id, *ed;
  int partwgt, partcut;
  BucketListType part;
  VertexType *vtx;
  int *where, *tmpptr, *bestwhere;
  int kwgt, bestcut, bestwgt;
  int orgvwgt, growp, shrinkp;

  ntrials = amin(GGGP_NTRIALS, graph->nvtxs);

  if (graph->tvwgt/graph->nvtxs < 6)
    ntrials = amax(1, ntrials/2);

  if (__Ctrl->dbglvl&DBG_INITCUT) {
    printf("GGGP-%d ", ntrials);
    fflush(stdout);
  }

  bestwhere = imalloc(graph->nvtxs, "BFSPartition: bestwhere");
  where = imalloc(graph->nvtxs, "BFSPartition: graph->where");
  touched = icoremalloc(graph->nvtxs, "BFSPartition: touched", 0);
  id = icoremalloc(graph->nvtxs, "BFSPartition: graph->id", 0);
  ed = icoremalloc(graph->nvtxs, "BFSPartition: graph->ed", 0);

  growp = 1; shrinkp = 0;
  orgvwgt = graph->tvwgt - zeropwgt;
  bestcut = orgvwgt*2000;

  initbucket(&part, 2*orgvwgt, graph->nvtxs, graph->nvtxs, -1);

  higain = RandomInRange(graph->nvtxs);
  for (ii=0; ii<ntrials; ii++) {
    if (__Ctrl->dbglvl&DBG_INITCUT) {
      printf("(%4d, ",graph->vtxs[higain]->vwgt);
      fflush(stdout);
    }

    for (i=0; i<graph->nvtxs; i++) 
      id[i] = graph->vtxs[i]->ewgtsum;
    iset(graph->nvtxs, shrinkp, where);
    iset(graph->nvtxs, 0, ed);
    iset(graph->nvtxs, 0, touched);

    resetbucket(&part);

    partwgt = 0;
    partcut = 0;
    touched[higain] = 1;
    Add2Part(&part, higain, ed[higain] - id[higain]);

    while (partwgt <= orgvwgt) {
      higain = GetMaxGainVtx(&part);
      if (higain == -1)  { /* Disconnected graph */
        for (;;) {
          higain = RandomInRange(graph->nvtxs);
          if (where[higain] == shrinkp)
            break;
        }
        touched[higain] = 1;
      }

      vtx = graph->vtxs[higain];

      if (orgvwgt < partwgt+vtx->vwgt)
        break;

      partwgt += vtx->vwgt;
      partcut -= (ed[higain] - id[higain]);
      where[higain] = growp;

      for (j=0; j<vtx->nedges; j++) {
        k = vtx->edges[j].edge;
        kwgt = vtx->edges[j].ewgt;
        if (where[k] == shrinkp) {  /* Look only the un inserted nodes */
          oldgain = ed[k] - id[k];
          id[k] -= kwgt;
          ed[k] += kwgt;

          if (touched[k] == 0) 
            Add2Part(&part, k, ed[k]-id[k]);
          else
            UpdatePart(&part, k, oldgain, ed[k]-id[k]);
          touched[k] = 1;
        }
      }
    }

    if (__Ctrl->dbglvl&DBG_INITCUT) {
      printf("%5d) ",partcut);
      fflush(stdout);
    }

    if (partcut < bestcut) {
      SWAP(bestwhere, where, tmpptr);
      bestcut = partcut;
      bestwgt = partwgt;
    }

    if (bestcut == 0)
      break;

    tmpptr = (ii == 0 ? bestwhere : where);
    for (;;) {
      higain = RandomInRange(graph->nvtxs);
      if (tmpptr[higain] == shrinkp)
        break;
    }
  }

  graph->mincut = bestcut;
  graph->where = bestwhere;

  freebucket(&part);
  icorefree(3*graph->nvtxs);
  free(where);
}



/*************************************************************************
* This function bisects a graph using the 2nd evector
**************************************************************************/
void EigPartition(CoarseGraphType *graph, int zeropwgt)
{
  int i;
  double *evec;
  int errstatus, sum;
  EvecWgtType *y;
  int *where;

  evec = (double *)GKmalloc(graph->nvtxs*sizeof(double), "EigPartition: evec");
  errstatus = lanczos(graph, evec);

  if (errstatus != 0) {
    printf("Something wrong with Lanczos... Switching to BFSPartition\n");
    GGGPPartition(graph, zeropwgt);
    return;
  }

  y = (EvecWgtType *)GKmalloc(graph->nvtxs*sizeof(EvecWgtType), "EigPartition: y");
  for (i=0; i<graph->nvtxs; i++) {
    y[i].ei = evec[i];
    y[i].wgt = graph->vtxs[i]->vwgt;
    y[i].index = i;
  }

  qsort((void *)y, (size_t)graph->nvtxs, (size_t)sizeof(EvecWgtType), inccompeinz);

  where = graph->where = imalloc(graph->nvtxs, "EigPartition: graph->where");
  sum = 0;
  for (i=0; i<graph->nvtxs; i++) {
    sum += y[i].wgt;
    if (sum > zeropwgt) {
      if (abs(sum-zeropwgt) > abs(sum-y[i].wgt-zeropwgt))
        i--;
      break;
    }
    where[y[i].index] = 0;
  }

  for (; i<graph->nvtxs; i++)
    where[y[i].index] = 1;

  GKfree(evec, y, -1);
}


/*************************************************************************
* This function is used by qsort to sort doubles in increasing order 
**************************************************************************/
int inccompeinz(const void *v1, const void *v2)
{
  if (((EvecWgtType *)v1)->ei < ((EvecWgtType *)v2)->ei)
    return -1;
  else if (((EvecWgtType *)v1)->ei > ((EvecWgtType *)v2)->ei)
    return 1;
  else  {
    if (((EvecWgtType *)v1)->wgt < ((EvecWgtType *)v2)->wgt)
      return -1;
    else if (((EvecWgtType *)v1)->wgt > ((EvecWgtType *)v2)->wgt)
      return 1;
    else
      return 0;
  }
}


/*************************************************************************
* This function computes the cut of a given partition
**************************************************************************/
int ComputeCut(CoarseGraphType *graph)
{
  int cut;
  int i, j;
  VertexType *vtx;
  int *where;

  where = graph->where;
  cut = 0;

  for (i=0; i<graph->nvtxs; i++) {
    vtx = graph->vtxs[i];
    for (j=0; j<vtx->nedges; j++) 
      if (where[vtx->edges[j].edge] != where[i]) 
        cut += vtx->edges[j].ewgt;
  }

  return cut/2;
}


