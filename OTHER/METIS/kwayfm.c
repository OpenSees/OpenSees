/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * kwayfm.c
 *
 * This file contains functions that implement the Kernighan-Lin partition 
 * refinment algorithm.
 * This scheme uses an array of buckets to store gains for efficiency.
 * It also implements the FM algorithm instead of KL
 *
 * Started 6/30/95
 * George
 *
 * $Id: kwayfm.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include "multilevel.h"

/*************************************************************************
* External Variables
**************************************************************************/
extern CtrlType *__Ctrl;	/* mlevelpart.c */
#ifndef METISLIB
extern timer GreedyTmr;		/* main.c */
extern timer GreedyInitTmr;	/* main.c */
extern timer GreedyIterTmr;	/* main.c */
extern timer GreedyWrapUpTmr;	/* main.c */
#endif



/*************************************************************************
* This function performs Kernighan-Lin refinement by swapping a vertex 
* at a time. Essentially implements the FM algorithm.
**************************************************************************/
void KWay_RefineFM(CoarseGraphType *graph, int nparts, int npasses)
{
  int i, j, k, imax;
  BucketListType buckets; 
  int *moved, *swaps, nswaps;
  int higain;			/* High gain node from a partition */
  int mincut, mincutorder, newcut;	/* The minimum cut so far, the order at which it occured */
  int *where;
  RInfoType *rinfo;
  EdgeType *degrees;
  int initcut;
  int pass=0;
  int *pwgts;                 /* The weights of the partitions */
  int from, to, limit, vwgt;
  int minpwgt, maxpwgt, badmaxpwgt;

  starttimer(&GreedyInitTmr);

  initbucket(&buckets, graph->tvwgt, graph->nvtxs, graph->nvtxs, graph->level);

  moved = icoremalloc(graph->nvtxs, "KWay_RefineFM: moved", 0);
  swaps = icoremalloc(graph->nvtxs, "KWay_RefineFM: swaps", 0);

  rinfo = graph->rinfo;
  where = graph->where;
  pwgts = graph->kpwgts;

  badmaxpwgt = UNBALANCE_FRACTION*graph->tvwgt/nparts;

  stoptimer(&GreedyInitTmr);

  for (pass=0; pass<npasses; pass++) {
    starttimer(&GreedyInitTmr);

    resetbucket(&buckets);

    iset(graph->nvtxs, -2, moved);

    mincutorder = 0;
    initcut = mincut = newcut = graph->mincut;
    limit = amax(0.01*graph->nvtxs, 200);
    minpwgt = pwgts[iamin(nparts, pwgts)];
    maxpwgt = pwgts[iamax(nparts, pwgts)];

    for (i=0; i<graph->nvtxs; i++) {
      if (rinfo[i].ed >= rinfo[i].id)  { /* Add only boundary vertices into the partition */
        imax = GetMaxEwgtI(rinfo[i].ndegrees, rinfo[i].degrees);
        if (imax != -1) {
          Add2Part(&buckets, i, rinfo[i].degrees[imax].ewgt-rinfo[i].id);
          moved[i] = -1;
        }
      }
    }
    stoptimer(&GreedyInitTmr);

    if (__Ctrl->dbglvl&DBG_KWAYREF && pass == 0)
        printf("Partitions: Nvtxs: %6d, Npart: %6d, InitCut: %8d, [%6d, %6d]\n",
           graph->nvtxs, buckets.nnodes, newcut, minpwgt, maxpwgt);


    /******************************************************
    * Get into the FM loop
    *******************************************************/
    starttimer(&GreedyIterTmr);
    nswaps = 0;
    for (;;) {
      if ((higain = GetMaxGainVtx(&buckets)) == -1) 
        break;

      /* Select where you are going to put this vertex */
      degrees = rinfo[higain].degrees;
      from = where[higain];
      vwgt = graph->vtxs[higain]->vwgt;
      moved[higain] = from;	/* Mark it as moved, anyway */

      for (k=0; k<rinfo[higain].ndegrees; k++) {
        if (pwgts[from]-vwgt >= minpwgt && pwgts[degrees[k].edge]+vwgt <= badmaxpwgt)
          break;
      }

      if (k < rinfo[higain].ndegrees) {
        for (j=k+1; j<rinfo[higain].ndegrees; j++) {
          if ((degrees[j].ewgt > degrees[k].ewgt && pwgts[degrees[j].edge]+vwgt <= badmaxpwgt) ||
              (degrees[j].ewgt == degrees[k].ewgt && pwgts[degrees[j].edge] < pwgts[degrees[k].edge]))
          k = j;
        }
        

        to = degrees[k].edge;
        newcut -= (degrees[k].ewgt - rinfo[higain].id);
        swaps[nswaps++] = higain;

        KWayFMUpdateDegreesI(graph, higain, to, &buckets, moved, 1);

        if (newcut <= mincut) {
          mincut = newcut;
          mincutorder = nswaps;
        }
        else {
          if (nswaps - mincutorder > limit)
            break; /* No further improvement, break out */
        }

        if (__Ctrl->dbglvl&DBG_ITERCUT) 
          printf("Moving %8d from %3d to %3d, Newcut: %6d [%6d], [%5d, %5d]\n", higain, from, to,
                  newcut, ComputeCut(graph), pwgts[iamin(nparts, pwgts)], pwgts[iamax(nparts, pwgts)]);
      }
      else 
        moved[higain] = -3;
    }
    stoptimer(&GreedyIterTmr);

    if (__Ctrl->dbglvl&DBG_KWAYREF)
      printf("Minimum Cut: %8d at %5d [%5d], [%5d, %5d] (%d)\n", mincut, mincutorder, nswaps, pwgts[iamin(nparts, pwgts)], pwgts[iamax(nparts, pwgts)], higain);

    /****************************************************************
    * Roll back computation 
    *****************************************************************/
    starttimer(&GreedyWrapUpTmr);
    for (nswaps--; nswaps>=mincutorder; nswaps--) {
      higain = swaps[nswaps];
      to = moved[higain];

      KWayFMUpdateDegreesI(graph, higain, to, &buckets, moved, 0);
    }
    stoptimer(&GreedyWrapUpTmr);

    ASSERT(KWayCheckDegrees(graph));

    graph->mincut = mincut;
    if (mincutorder == 0 || mincut >= initcut)
      break;
  }

  freebucket(&buckets);

  icorefree(2*graph->nvtxs);

}



/*************************************************************************
* This function updates the degrees of a k-way partition
**************************************************************************/
void KWayFMUpdateDegreesI(CoarseGraphType *graph, int v, int to, BucketListType *buckets, 
                          int *moved, int flag)
{
  int i, j, k, oldgain, newgain, from, imax;
  RInfoType *myrinfo;
  EdgeType *degrees;
  VertexType *vtx;

  myrinfo = graph->rinfo + v;
  degrees = myrinfo->degrees;

  if (myrinfo->ndegrees > 0) {
    for (k=0; k<myrinfo->ndegrees; k++)
      if (degrees[k].edge == to)
        break;

    /* Take care the special but BAD case */
    if (k == myrinfo->ndegrees) {  
      degrees[k].edge = to;
      degrees[k].ewgt = 0;
      myrinfo->ndegrees++;
    }
  }
  else {
    if (myrinfo->degrees == NULL) 
      myrinfo->degrees = GetnExtDegrees(graph->vtxs[v]->nedges);

    myrinfo->ndegrees = 1;
    degrees[0].ewgt = 0;
    k = 0;
  }

  myrinfo->ed += (myrinfo->id - degrees[k].ewgt);
  SWAP(myrinfo->id, degrees[k].ewgt, oldgain);

  from = graph->where[v];
  graph->where[v] = to;
  degrees[k].edge = from;

  if (degrees[k].ewgt == 0) {
    for (i=k+1; i<myrinfo->ndegrees; i++)
      degrees[i-1] = degrees[i];
    myrinfo->ndegrees--;
  }

  vtx = graph->vtxs[v];
  for (i=0; i<vtx->nedges; i++) {
    j = vtx->edges[i].edge;
    if (flag && moved[j] == -1) {
      imax = GetMaxEwgtI(graph->rinfo[j].ndegrees, graph->rinfo[j].degrees);
      oldgain = (imax == -1 ? -graph->rinfo[j].id : graph->rinfo[j].degrees[imax].ewgt - graph->rinfo[j].id);
    }

    KWayUpdateVtxDegrees(graph, j, vtx->edges[i].ewgt, from, to);

    if (flag) {
      if (moved[j] == -1) {
        imax = GetMaxEwgtI(graph->rinfo[j].ndegrees, graph->rinfo[j].degrees);
        newgain = (imax == -1 ? -graph->rinfo[j].id : graph->rinfo[j].degrees[imax].ewgt - graph->rinfo[j].id);

        UpdatePart(buckets, j, oldgain, newgain);
      }
      else if (moved[j] == -3 && graph->rinfo[j].ed - graph->rinfo[j].id > 0) {
        imax = GetMaxEwgtI(graph->rinfo[j].ndegrees, graph->rinfo[j].degrees);
        if (imax != -1) {
          Add2Part(buckets, j, graph->rinfo[j].degrees[imax].ewgt - graph->rinfo[j].id);
          moved[j] = -1;
        }
      }
    }
  }

  INC_DEC(graph->kpwgts[to], graph->kpwgts[from], vtx->vwgt);
}




/*************************************************************************
* This function performs Kernighan-Lin refinement by swapping a vertex 
* at a time. Essentially implements the FM algorithm.
**************************************************************************/
void KWay_BalanceFM(CoarseGraphType *graph, int nparts, int npasses)
{
  int i, j, k;
  BucketListType buckets; 
  int *moved, nswaps;
  int higain;			/* High gain node from a partition */
  int *where;
  RInfoType *rinfo;
  EdgeType *degrees;
  int pass=0;
  int *pwgts;                 /* The weights of the partitions */
  int from, to, vwgt;
  int oldmaxpwgt, maxpwgt, goodsize, lastmax;

  starttimer(&GreedyInitTmr);

  rinfo = graph->rinfo;
  where = graph->where;
  pwgts = graph->kpwgts;

  goodsize = UNBALANCE_FRACTION*graph->tvwgt;
  maxpwgt = pwgts[iamax(nparts, pwgts)];

  stoptimer(&GreedyInitTmr);

  if (maxpwgt*nparts <= goodsize)
    return;

  initbucket(&buckets, graph->tvwgt, graph->nvtxs, graph->nvtxs, graph->level);
  moved = icoremalloc(graph->nvtxs, "KWay_BalanceFM: moved", 0);

  for (pass=0; pass<npasses; pass++) {
    starttimer(&GreedyInitTmr);

    resetbucket(&buckets);

    iset(graph->nvtxs, -2, moved);

    for (i=0; i<graph->nvtxs; i++) {
      if (rinfo[i].ed-rinfo[i].id > -graph->tvwgt/graph->nvtxs)  { /* Add only boundary vertices into the partition */
        Add2Part(&buckets, i, rinfo[i].ed-rinfo[i].id);
        moved[i] = -1;
      }
    }
    stoptimer(&GreedyInitTmr);

    if (__Ctrl->dbglvl&DBG_KWAYREF) {
      if (pass == 0)
        printf("Balance: Nvtxs: %6d, InitCut: %8d\n", graph->nvtxs, graph->mincut);
      printf("\t\tNpart: %5d, InitMax: %4d,", buckets.nnodes, maxpwgt);
    }

    /******************************************************
    * Get into the FM loop
    *******************************************************/
    starttimer(&GreedyIterTmr);
    nswaps = 0;
    oldmaxpwgt = maxpwgt;
    while (maxpwgt*nparts > goodsize) {
      if ((higain = GetMaxGainVtx(&buckets)) == -1)
        break;

      /* Select where you are going to put this vertex */
      degrees = rinfo[higain].degrees;
      from = where[higain];
      vwgt = graph->vtxs[higain]->vwgt;
      moved[higain] = from;	/* Mark it as moved, anyway */

      for (k=0; k<rinfo[higain].ndegrees; k++) {
        if (pwgts[from] - pwgts[degrees[k].edge] > 2*vwgt) 
          break;
      }

      if (k < rinfo[higain].ndegrees) {
        nswaps = 0;
        for (j=k+1; j<rinfo[higain].ndegrees; j++) {
          if (pwgts[degrees[j].edge] < pwgts[degrees[k].edge])
            k = j;
        }
        
        to = degrees[k].edge;
        graph->mincut -= (degrees[k].ewgt - rinfo[higain].id);

        if (pwgts[from] == maxpwgt)
          lastmax = 1;

        ASSERT(from != to);

        KWayFMUpdateDegreesBal(graph, higain, to, &buckets, moved);

        if (lastmax) {
          maxpwgt = pwgts[iamax(nparts, pwgts)];
          lastmax = 0;
        }

        if (__Ctrl->dbglvl&DBG_ITERCUT) 
          printf("Moving %8d from %3d to %3d, NewMax: %4d, [%4d]\n", higain, from, to,
                  maxpwgt, pwgts[iamin(nparts, pwgts)]);
      }
      else {
        if (++nswaps >= 2000)
          break;
      }
    }
    stoptimer(&GreedyIterTmr);

    if (__Ctrl->dbglvl&DBG_KWAYREF)
      printf("Cut: %6d, Maxpwgt: %4d, [%4d] [%4d] (%d)\n", graph->mincut, maxpwgt, pwgts[iamin(nparts, pwgts)], buckets.nnodes, higain);

    if (oldmaxpwgt == maxpwgt || maxpwgt*nparts <= goodsize)
      break;

  }

  freebucket(&buckets);

  icorefree(graph->nvtxs);
}


/*************************************************************************
* This function updates the degrees of a k-way partition
**************************************************************************/
void KWayFMUpdateDegreesBal(CoarseGraphType *graph, int v, int to, BucketListType *buckets, int *moved)
{
  int i, j, k, tmp, from;
  RInfoType *myrinfo;
  EdgeType *degrees;
  VertexType *vtx;

  myrinfo = graph->rinfo + v;
  degrees = myrinfo->degrees;

  if (myrinfo->ndegrees > 0) {
    for (k=0; k<myrinfo->ndegrees; k++)
      if (degrees[k].edge == to)
        break;
  }
  else {
    if (myrinfo->degrees == NULL) 
      myrinfo->degrees = GetnExtDegrees(graph->vtxs[v]->nedges);

    myrinfo->ndegrees = 1;
    degrees[0].ewgt = 0;
    k = 0;
  }

  myrinfo->ed += (myrinfo->id - degrees[k].ewgt);
  SWAP(myrinfo->id, degrees[k].ewgt, tmp);

  from = graph->where[v];
  graph->where[v] = to;
  degrees[k].edge = from;

  if (degrees[k].ewgt == 0) {
    for (i=k+1; i<myrinfo->ndegrees; i++)
      degrees[i-1] = degrees[i];
    myrinfo->ndegrees--;
  }

  vtx = graph->vtxs[v];
  for (i=0; i<vtx->nedges; i++) {
    j = vtx->edges[i].edge;
    tmp = graph->rinfo[j].ed - graph->rinfo[j].id;

    KWayUpdateVtxDegrees(graph, j, vtx->edges[i].ewgt, from, to);

    if (moved[j] == -1)
      UpdatePart(buckets, j, tmp, graph->rinfo[j].ed - graph->rinfo[j].id);
  }

  INC_DEC(graph->kpwgts[to], graph->kpwgts[from], vtx->vwgt);
}



/*************************************************************************
* This function returns the maximum edge weight
**************************************************************************/
int GetMaxEwgtI(int nedges, EdgeType *edges)
{
  int i, max;

  if (nedges == 0)
    return -1;

  max = 0;
  for (i=1; i<nedges; i++) 
    if (edges[i].ewgt > edges[max].ewgt)
      max = i;

  return max;
}

