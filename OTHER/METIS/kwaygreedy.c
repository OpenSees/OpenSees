/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * kwaygreedy.c
 *
 * This file contains functions that implement a k-way greedy refinement algorithm
 *
 * Started 6/7/95
 * George
 *
 * $Id: kwaygreedy.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include "multilevel.h"

/*************************************************************************
* External Variables
**************************************************************************/
extern CtrlType *__Ctrl;        /* mlevelpart.c */
#ifndef METISLIB
extern timer GreedyTmr;         /* main.c */
extern timer GreedyInitTmr;     /* main.c */
extern timer GreedyIterTmr;     /* main.c */
extern timer GreedyWrapUpTmr;   /* main.c */
#endif


/*************************************************************************
* This function performs a k-way FM refinement
**************************************************************************/
void KWay_RefineGreedy(CoarseGraphType *graph, int nparts, int npasses)
{
  int i, j, k, l;			/* The node and its adjacent edge weight */
  int *where;
  int *pwgts;                 /* The weights of the partitions */
  int from, vwgt;
  RInfoType *rinfo;
  EdgeType *degrees;
  int swaps;
  int minpwgt, maxpwgt, badmaxpwgt, oldcut;

  rinfo = graph->rinfo;
  where = graph->where;
  pwgts = graph->kpwgts;

  starttimer(&GreedyInitTmr);

  minpwgt = pwgts[iamin(nparts, pwgts)];
  maxpwgt = pwgts[iamax(nparts, pwgts)];
  badmaxpwgt = UNBALANCE_FRACTION*graph->tvwgt/nparts;

  if (__Ctrl->dbglvl&DBG_KWAYREF)
    printf("Partitions: \t[%6d %6d] %6d", minpwgt, maxpwgt, graph->mincut);

  for (l=0; l<npasses; l++) {
    swaps = 0;
    oldcut = graph->mincut;
    minpwgt = pwgts[iamin(nparts, pwgts)];
    maxpwgt = pwgts[iamax(nparts, pwgts)];

    if (l > 0 && __Ctrl->dbglvl&DBG_KWAYREF)
      printf("\t\t[%4d %4d] [%4d] %6d", minpwgt, maxpwgt, badmaxpwgt, graph->mincut);

    for (i=0; i<graph->nvtxs; i++) {
      if (rinfo[i].ed >= rinfo[i].id) {
        degrees = rinfo[i].degrees;
        from = where[i];
        vwgt = graph->vtxs[i]->vwgt;

        for (k=0; k<rinfo[i].ndegrees; k++) {
          if (pwgts[from]-vwgt >= minpwgt && pwgts[degrees[k].edge]+vwgt <= badmaxpwgt) 
            break;
        }

        if (k < rinfo[i].ndegrees) {
          for (j=k+1; j<rinfo[i].ndegrees; j++) {
           if ((degrees[j].ewgt > degrees[k].ewgt && pwgts[degrees[j].edge]+vwgt <= badmaxpwgt) ||
               (degrees[j].ewgt == degrees[k].ewgt && pwgts[degrees[j].edge] < pwgts[degrees[k].edge]))
              k = j;
          }

          if (degrees[k].ewgt > rinfo[i].id || (degrees[k].ewgt == rinfo[i].id && (pwgts[from] - pwgts[degrees[k].edge] >= vwgt))) {
            graph->mincut -= (degrees[k].ewgt-rinfo[i].id);

            KWayUpdateDegrees(graph, i, k);
            swaps++;
          }
        }
      }
    }

    if (__Ctrl->dbglvl&DBG_KWAYREF)
      printf(" %d\n", swaps);

    if (graph->mincut == oldcut && minpwgt == pwgts[iamin(nparts, pwgts)] && maxpwgt == pwgts[iamax(nparts, pwgts)])
      break;
  }

  stoptimer(&GreedyInitTmr);

}






/*************************************************************************
* This function updates the degrees of a k-way partition
**************************************************************************/
void KWayUpdateDegrees(CoarseGraphType *graph, int v, int k)
{
  int i, tmp, from, to;
  RInfoType *myrinfo;
  EdgeType *degrees;
  VertexType *vtx;

  myrinfo = graph->rinfo + v;
  degrees = myrinfo->degrees;

  myrinfo->ed += (myrinfo->id - degrees[k].ewgt);

  SWAP(myrinfo->id, degrees[k].ewgt, tmp);

  from = graph->where[v];
  graph->where[v] = to = degrees[k].edge;
  vtx = graph->vtxs[v];

  degrees[k].edge = from;

  if (degrees[k].ewgt == 0) {
    for (i=k+1; i<myrinfo->ndegrees; i++)
      degrees[i-1] = degrees[i];
    myrinfo->ndegrees--;
  }


  for (i=0; i<vtx->nedges; i++) 
    KWayUpdateVtxDegrees(graph, vtx->edges[i].edge, vtx->edges[i].ewgt, from, to);

  INC_DEC(graph->kpwgts[to], graph->kpwgts[from], vtx->vwgt);
}


/*************************************************************************
* This function recomputes the ext/int degrees of a vertex
**************************************************************************/
void KWayUpdateVtxDegrees(CoarseGraphType *graph, int v, int ewgt, int from, int to)
{
  int i, j, mypart, *where;
  RInfoType *myrinfo;
  EdgeType *degrees;

  ASSERT(from != to);

  where = graph->where;
  myrinfo = graph->rinfo + v;
  degrees = myrinfo->degrees;
  mypart = where[v];

  if (mypart == from) {
    INC_DEC(myrinfo->ed, myrinfo->id, ewgt);

    for (i=0; i<myrinfo->ndegrees; i++) {
      if (degrees[i].edge == to) {
        degrees[i].ewgt += ewgt;
        break;
      }
    }
    if (i == myrinfo->ndegrees) {
      if (degrees == NULL) 
        degrees = myrinfo->degrees = GetnExtDegrees(graph->vtxs[v]->nedges);

      degrees[i].ewgt = ewgt;
      degrees[i].edge = to;
      myrinfo->ndegrees++;
    }
  }
  else if (mypart == to) {
    INC_DEC(myrinfo->id, myrinfo->ed, ewgt);

    for (i=0; i<myrinfo->ndegrees; i++) {
      if (degrees[i].edge == from) {
        degrees[i].ewgt -= ewgt;
        break;
      }
    }
    ASSERT(i<myrinfo->ndegrees);

    if (degrees[i].ewgt == 0) {
      for (j=i+1; j<myrinfo->ndegrees; j++)
        degrees[j-1] = degrees[j];
      myrinfo->ndegrees--;
    }
  }
  else {
    for (i=0; i<myrinfo->ndegrees; i++) {
      if (degrees[i].edge == from) {
        degrees[i].ewgt -= ewgt;
        break;
      }
    }
    ASSERT(i<myrinfo->ndegrees);

    if (degrees[i].ewgt == 0) {
      for (j=i+1; j<myrinfo->ndegrees; j++)
        degrees[j-1] = degrees[j];
      myrinfo->ndegrees--;
    }


    for (i=0; i<myrinfo->ndegrees; i++) {
      if (degrees[i].edge == to) {
        degrees[i].ewgt += ewgt;
        break;
      }
    }
    if (i == myrinfo->ndegrees) {
      if (degrees == NULL) 
        degrees = myrinfo->degrees = GetnExtDegrees(graph->vtxs[v]->nedges);

      degrees[i].ewgt = ewgt;
      degrees[i].edge = to;
      myrinfo->ndegrees++;
    }
  }

  ASSERT(myrinfo->ndegrees <= graph->vtxs[v]->nedges);

}



