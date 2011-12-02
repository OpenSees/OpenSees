/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * mlevelpart.c
 *
 * This file contains the driving code for the multilevel graph partitioning
 * algorithm.
 *
 * Started 10/5/94
 * George
 *
 * $Id: mlevelpart.c,v 1.2 2001-10-05 00:51:05 fmk Exp $
 *
 */

#include "multilevel.h"

/*************************************************************************
* Global Variables
**************************************************************************/
CtrlType *__Ctrl; 	/* Stores Control Information */

/*************************************************************************
* External Global Variables
**************************************************************************/
#ifndef METISLIB
extern timer SplitTmr;		/* main.c */
#endif


/*************************************************************************
* This function is called recursively to partition the graph
* The array kpwgts stores the target weights of the partitions
**************************************************************************/
int MultiLevelPart(CoarseGraphType *graph, int nparts, int CoarsenTo, int MatchType, 
                   int InitPartType, int RefineType, int dbglvl, int IsWeighted, 
                   int *part, int *cuts, int *kpwgts)
{
  int i, j;
  int totalcut;
  CtrlType *ctrl, *oldctrl;

  ctrl = (CtrlType *)GKmalloc(sizeof(CtrlType), "MultilevelPart: ctrl");

  ctrl->CoarsenTo = CoarsenTo;
  ctrl->MatchType = MatchType;
  ctrl->InitPartType = InitPartType;
  ctrl->RefineType = RefineType;
  ctrl->IsWeighted = IsWeighted;
  ctrl->OpType = OP_RMLB;
  ctrl->dbglvl = dbglvl;
  ctrl->nparts = 2;
  ctrl->cfrac = COARSEN_FRACTION;
  
  if (graph->nedges/graph->nvtxs > 20)
    ctrl->maxedges = 1.8*graph->nedges; /* Set the maxpool to something large */
  else
    ctrl->maxedges = 2.4*graph->nedges; /* Set the maxpool to something large */

  ctrl->maxicore = 2*graph->nvtxs;
  ctrl->maxgain = 2*graph->nvtxs;
  ctrl->maxbucket = 2*graph->nvtxs;

  AllocatePools(ctrl);

  /* Set current Ctrl context */
  oldctrl = __Ctrl;
  __Ctrl = ctrl;

  /* 
   * Set the kpwgts array. At a later point you may want the user to be
   * able to do that.
   */
#ifdef MODULE_WEIGHTS
  j = (graph->nvtxs+nparts-1)/nparts;
  if (!IsWeighted && j*(nparts-1) < graph->nvtxs) { /* Let the last guy have a module of the work */
    j = (graph->nvtxs+nparts-1)/nparts;
    for (i=0; i<nparts-1; i++)
      kpwgts[i] = j;
    kpwgts[i] = graph->nvtxs - j*(nparts-1);
  }
  else {
    j = graph->nvtxs;
    for (i=0; i<nparts; i++) {
      kpwgts[i] = j/(nparts-i);
      j -= kpwgts[i];
    }
  }
#endif

  totalcut = 0;

  RMLB(graph, nparts, part, 0, 1, cuts, kpwgts, &totalcut);

  FreePools(ctrl);
  free(ctrl);

  /* Restore old Ctrl context */
  __Ctrl = oldctrl;

  return totalcut;
}
  

/*************************************************************************
* This function is called recursively to partition the graph
**************************************************************************/
int RMLB(CoarseGraphType *graph, int nparts, int *part, int fpart, int tnode, int *cuts, int *kpwgts, int *totalcut)
{
  int i, j;
  CoarseGraphType *coarsegraph, lgraph, rgraph;
  int zeropwgt;

  if (graph->nvtxs <= 1) 
    return 0;

  /* Determine the weight of the first part */
#ifdef MODULE_WEIGHTS
  j = fpart+nparts;
  for (i=fpart, zeropwgt=0; i<j; i++)
    zeropwgt += kpwgts[i];
  if (zeropwgt == graph->tvwgt) { /* Sizes were fine so far */
    j = fpart+nparts/2;
    for (i=fpart, zeropwgt=0; i<j; i++)
      zeropwgt += kpwgts[i];
  }
  else {
    zeropwgt = ((graph->tvwgt+1)*(nparts>>1))/nparts;
  }
#else
  zeropwgt = ((graph->tvwgt+1)*(nparts>>1))/nparts;
#endif

  ASSERT(__Ctrl->lastedge == 0);

  coarsegraph = Coarsen(graph, amin(__Ctrl->CoarsenTo, amax(20, 0.5*graph->nvtxs)));

  Refine(graph, coarsegraph, zeropwgt);

  *totalcut += graph->mincut;
  cuts[tnode] = graph->mincut;

  if (nparts == 3) {
    kpwgts[fpart] = graph->pwgts[0];

    /* Label the 0-partition */
    for (i=0; i<graph->nvtxs; i++) 
      if (graph->where[i] == 0)
        part[graph->label[i]] = fpart;
  }
  if (nparts == 2) {
    kpwgts[fpart] = graph->pwgts[0];
    kpwgts[fpart+1] = graph->pwgts[1];

    /* Label both partitions */
    for (i=0; i<graph->nvtxs; i++) 
      part[graph->label[i]] = fpart + graph->where[i];
  }


  if (nparts > 3) 
    SplitGraphPart(graph, &lgraph, &rgraph);

  if (nparts == 3)
    SplitGraphPart1_2(graph, &rgraph);


  FreeRootGraph(graph);

  if (nparts > 3) {
    RMLB(&rgraph, (nparts+1)/2, part, fpart+nparts/2, tnode<<1, cuts, kpwgts, totalcut);
    RMLB(&lgraph, nparts/2, part, fpart, (tnode<<1)+1, cuts, kpwgts, totalcut);
  }
  if (nparts == 3)
    RMLB(&rgraph, (nparts+1)/2, part, fpart+nparts/2, tnode<<1, cuts, kpwgts, totalcut);

}



/*************************************************************************
* This function splits a graph acording to graph->where[] array.
* The two graphs created are stored in lgraph (where==0) and 
* rgraph (where = 1), however, the actual space for each vertex is resused.
**************************************************************************/
void SplitGraphPart(CoarseGraphType *graph, CoarseGraphType *lgraph, CoarseGraphType *rgraph)
{
  int i, j;
  int *where;
  int *map;		/* Points to graph->cmap. Now it is used for different purpose */
   			/* map[i] = j, v_i in graph is mapped in v_j in lgraph or rgraph */
  EdgeType *ledges;	/* Pointers to edges and weigths of left graph */
  EdgeType *redges;	/* Pointers to edges and weigths of right graph */
  int li, ri;		/* lgraph and rgraph indices */
  VertexType *oldvtx;
  int nedges, ewgtsum;

  starttimer(&SplitTmr);

  InitGraph(lgraph);
  InitGraph(rgraph);

  if (graph->cmap == NULL)
    graph->cmap = imalloc(graph->nvtxs, "SplitGraphPart: graph->cmap");

  where = graph->where;
  map = graph->cmap; 	/* Use cmap as an aux vector */

  lgraph->nvtxs = 0;
  lgraph->nedges = 0;
  lgraph->tvwgt = graph->pwgts[0];

  rgraph->nvtxs = 0;
  rgraph->nedges = 0;
  rgraph->tvwgt = graph->pwgts[1];

  for (i=0; i<graph->nvtxs; i++) {
    if (where[i] == 0) {
      map[i] = lgraph->nvtxs;
      lgraph->nvtxs++;
      lgraph->nedges += graph->vtxs[i]->nedges;
    }
    else {
      map[i] = rgraph->nvtxs;
      rgraph->nvtxs++;
      rgraph->nedges += graph->vtxs[i]->nedges;
    }
  }

  lgraph->vtxs = (VertexType **)GKmalloc(sizeof(VertexType *)*lgraph->nvtxs, "SplitGraphPart: lgraph->vtxs");
  lgraph->label = imalloc(lgraph->nvtxs, "SplitGraphPart: lgraph->label");
  li = 0;

  rgraph->vtxs = (VertexType **)GKmalloc(sizeof(VertexType *)*rgraph->nvtxs, "SplitGraphPart: rgraph->vtxs");
  rgraph->label = imalloc(rgraph->nvtxs, "SplitGraphPart: rgraph->label");
  ri = 0;


  for (i=0; i<graph->nvtxs; i++) {
    oldvtx = graph->vtxs[i];

    if (where[i] == 0) {  /* Put this vertex in the lgraph */
      lgraph->vtxs[li] = oldvtx;
      ledges = oldvtx->edges;
      nedges = 0;
      ewgtsum = 0;

      for (j=0; j<oldvtx->nedges; j++) {
        if (where[oldvtx->edges[j].edge] == 0) {  /* I keep this edge */
          ledges[nedges].edge = map[oldvtx->edges[j].edge];
          ledges[nedges].ewgt = oldvtx->edges[j].ewgt;
          ewgtsum += ledges[nedges].ewgt;
          nedges++;
        }
      }

      lgraph->vtxs[li]->nedges = nedges;
      lgraph->vtxs[li]->ewgtsum = ewgtsum;
      lgraph->label[li] = graph->label[i];
      li++;
    }
    else { /* Put this vertex in the rgraph */
      rgraph->vtxs[ri] = oldvtx;
      redges = oldvtx->edges;
      nedges = 0;
      ewgtsum = 0;

      for (j=0; j<oldvtx->nedges; j++) {
        if (where[oldvtx->edges[j].edge] == 1) {  /* I keep this edge */
          redges[nedges].edge = map[oldvtx->edges[j].edge];
          redges[nedges].ewgt = oldvtx->edges[j].ewgt;
          ewgtsum += redges[nedges].ewgt;
          nedges++;
        }
      }

      rgraph->vtxs[ri]->nedges = nedges;
      rgraph->vtxs[ri]->ewgtsum = ewgtsum;
      rgraph->label[ri] = graph->label[i];
      ri++;
    }
  }

  stoptimer(&SplitTmr);
}



/*************************************************************************
* This function splits a graph acording to graph->where[] array.
* The two graphs created are stored in lgraph (where==0) and 
* rgraph (where = 1), however, the actual space for each vertex is resused.
**************************************************************************/
void SplitGraphPart1_2(CoarseGraphType *graph, CoarseGraphType *rgraph)
{
  int i, j;
  int *where;
  int *map;		/* Points to graph->cmap. Now it is used for different purpose */
   			/* map[i] = j, v_i in graph is mapped in v_j in lgraph or rgraph */
  EdgeType *redges;	/* Pointers to edges and weigths of right graph */
  int ri;		/* rgraph indices */
  VertexType *oldvtx;
  int nedges, ewgtsum;

  starttimer(&SplitTmr);

  InitGraph(rgraph);

  if (graph->cmap == NULL)
    graph->cmap = imalloc(graph->nvtxs, "SplitGraphPart: graph->cmap");

  where = graph->where;
  map = graph->cmap; 	/* Use cmap as an aux vector */

  rgraph->nvtxs = 0;
  rgraph->nedges = 0;
  rgraph->tvwgt = graph->pwgts[1];

  for (i=0; i<graph->nvtxs; i++) {
    if (where[i] == 1) {
      map[i] = rgraph->nvtxs;
      rgraph->nvtxs++;
      rgraph->nedges += graph->vtxs[i]->nedges;
    }
  }

  rgraph->vtxs = (VertexType **)GKmalloc(sizeof(VertexType *)*rgraph->nvtxs, "SplitGraphPart: rgraph->vtxs");
  rgraph->label = imalloc(rgraph->nvtxs, "SplitGraphPart: rgraph->label");
  ri = 0;

  for (i=0; i<graph->nvtxs; i++) {
    oldvtx = graph->vtxs[i];

    if (where[i] == 1) {  /* Put this vertex in the lgraph */
      rgraph->vtxs[ri] = oldvtx;
      redges = oldvtx->edges;
      nedges = 0;
      ewgtsum = 0;

      for (j=0; j<oldvtx->nedges; j++) {
        if (where[oldvtx->edges[j].edge] == 1) {  /* I keep this edge */
          redges[nedges].edge = map[oldvtx->edges[j].edge];
          redges[nedges].ewgt = oldvtx->edges[j].ewgt;
          ewgtsum += redges[nedges].ewgt;
          nedges++;
        }
      }

      rgraph->vtxs[ri]->nedges = nedges;
      rgraph->vtxs[ri]->ewgtsum = ewgtsum;
      rgraph->label[ri] = graph->label[i];
      ri++;
    }
  }

  stoptimer(&SplitTmr);
}


