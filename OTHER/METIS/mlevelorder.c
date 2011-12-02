/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * mlevelorder.c
 *
 * This file contains the driving code for the multilevel graph partitioning
 * algorithm.
 *
 * Started 10/5/94
 * George
 *
 * $Id: mlevelorder.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include "multilevel.h"

/*************************************************************************
* External Global Variables
**************************************************************************/
extern CtrlType *__Ctrl;	/* mlevelpart.c */
#ifndef METISLIB
extern timer SplitTmr;		/* main.c */
#endif

/*************************************************************************
* This function is the entry point of Recursive nested dissection
**************************************************************************/
void MultiLevelOrder(CoarseGraphType *graph, int CoarsenTo, int MatchType,
                     int InitPartType, int RefineType, int dbglvl, int *perm, SepNodeType *stree)
{
  CtrlType *ctrl, *oldctrl;


  ctrl = GKmalloc(sizeof(CtrlType), "MultilevelPart: ctrl");

  ctrl->CoarsenTo = CoarsenTo;
  ctrl->MatchType = MatchType;
  ctrl->InitPartType = InitPartType;
  ctrl->RefineType = RefineType;
  ctrl->dbglvl = dbglvl;
  ctrl->IsWeighted = 0;
  ctrl->OpType = OP_MLND;
  ctrl->nparts = 2;
  ctrl->cfrac = COARSEN_FRACTION;
  
  if (graph->nedges/graph->nvtxs > 20)
    ctrl->maxedges = 1.8*graph->nedges; /* Set the maxpool to something large */
  else
    ctrl->maxedges = 2.4*graph->nedges; /* Set the maxpool to something large */

  ctrl->maxdegrees = 0;
  ctrl->maxicore = 2*graph->nvtxs;
  ctrl->maxgain = 2*graph->nvtxs;
  ctrl->maxbucket = 2*graph->nvtxs;

  AllocatePools(ctrl);

  /* Set current Ctrl context */
  oldctrl = __Ctrl;
  __Ctrl = ctrl;

  MLND(graph, perm, graph->nvtxs, stree, 1);

#ifdef METISLIB
  FreePools(ctrl);
  free(ctrl);

  /* Restore old Ctrl context */
  __Ctrl = oldctrl;
#endif

}


/*************************************************************************
* This function is called recursively to partition the graph
* Parameters:
*   graph:	The graph itself
*    perm:	The permutation vector
*   order:	The vertices if graph will be numbered between
*                  (order-graph->nvtxs, order]
*   stree:	The separator tree
* mysnode:	The node into the separator tree of graph
**************************************************************************/
void MLND(CoarseGraphType *graph, int *perm, int order, SepNodeType *stree, int mysnode)
{
  CoarseGraphType *coarsegraph, lgraph, rgraph;

  stree[mysnode].nvtxs = graph->nvtxs;

  ResetPools();

  coarsegraph = Coarsen(graph, amax(50, __Ctrl->CoarsenTo));

  Refine(graph, coarsegraph, (graph->tvwgt+1)/2);

  SplitGraphOrder(graph, &lgraph, &rgraph, perm, order);

  stree[mysnode].hi = order;
  order -= (graph->nvtxs - lgraph.nvtxs - rgraph.nvtxs);
  stree[mysnode].li = order;
  stree[mysnode].isleaf = 0;

  FreeRootGraph(graph);

  if (rgraph.nvtxs <= MMDSWITCH) {
    MDOrder(&rgraph, perm, order, stree, 2*mysnode);
    FreeRootGraph(&rgraph);
  }
  else
    MLND(&rgraph, perm, order, stree, 2*mysnode);

  if (lgraph.nvtxs <= MMDSWITCH) {
    MDOrder(&lgraph, perm, order-rgraph.nvtxs, stree, 2*mysnode+1);
    FreeRootGraph(&lgraph);
  }
  else
    MLND(&lgraph, perm, order-rgraph.nvtxs, stree, 2*mysnode+1);

}



/*************************************************************************
* This function splits a graph acording to graph->where[] array.
* The two graphs created are stored in lgraph (where==0) and 
* rgraph (where = 1), however, the actual space for each vertex is resused.
**************************************************************************/
void SplitGraphOrder(CoarseGraphType *graph, CoarseGraphType *lgraph, CoarseGraphType *rgraph, 
                     int *perm, int order)
{
  int i, j;
  int *where;
  int *map;		/* Points to graph->cmap. Now it is used for different purpose */
   			/* map[i] = j, v_i in graph is mapped in v_j in lgraph or rgraph */
  EdgeType *ledges;	/* Pointers to edges and weigths of left graph */
  EdgeType *redges;	/* Pointers to edges and weigths of right graph */
  EdgeType *oedges;	/* Pinters to edges of the old graph */
  int li, ri;		/* lgraph and rgraph indices */
  VertexType *oldvtx;
  int nedges, ewgtsum;


  starttimer(&SplitTmr);

  /* Take care the separator first */
  FindMinCovNodeSeparator(graph, perm, order);

  InitGraph(lgraph);
  InitGraph(rgraph);

  if (graph->cmap == NULL)
    graph->cmap = imalloc(graph->nvtxs, "SplitGraphOrder: graph->cmap");

  where = graph->where;
  map = graph->cmap; 

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
    else if (where[i] == 1) {
      map[i] = rgraph->nvtxs;
      rgraph->nvtxs++;
      rgraph->nedges += graph->vtxs[i]->nedges;
    }
    else if (where[i] == 2) {
      lgraph->tvwgt -= graph->vtxs[i]->vwgt;
    }
    else if (where[i] == 3) {
      rgraph->tvwgt -= graph->vtxs[i]->vwgt;
    }
  }

  lgraph->vtxs = (VertexType **)GKmalloc(sizeof(VertexType *)*lgraph->nvtxs, "SplitGraph: lgraph->vtxs");
  lgraph->label = imalloc(lgraph->nvtxs, "SplitGraphOrder: lgraph->label");
  li = 0;

  rgraph->vtxs = (VertexType **)GKmalloc(sizeof(VertexType *)*rgraph->nvtxs, "SplitGraph: rgraph->vtxs");
  rgraph->label = imalloc(rgraph->nvtxs, "SplitGraphOrder: rgraph->label");
  ri = 0;

  for (i=0; i<graph->nvtxs; i++) {
    oldvtx = graph->vtxs[i];

    if (where[i] == 0) {  /* Put this vertex in the lgraph */
      lgraph->vtxs[li] = oldvtx;
      oedges = ledges = oldvtx->edges;
      nedges = ewgtsum = 0;

      for (j=0; j<oldvtx->nedges; j++) {
        if (where[oedges[j].edge] == 0) {  /* I keep this edge */
          ledges[nedges].edge = map[oedges[j].edge];
          ledges[nedges].ewgt = oedges[j].ewgt;
          ewgtsum += ledges[nedges].ewgt;
          nedges++;
        }
      }

      lgraph->vtxs[li]->nedges = nedges;
      lgraph->vtxs[li]->ewgtsum = ewgtsum;
      lgraph->label[li] = graph->label[i];
      li++;
    }
    else if (where[i] == 1) { /* Put this vertex in the rgraph */
      rgraph->vtxs[ri] = oldvtx;
      oedges = redges = oldvtx->edges;
      nedges = ewgtsum = 0;

      for (j=0; j<oldvtx->nedges; j++) {
        if (where[oedges[j].edge] == 1) {  /* I keep this edge */
          redges[nedges].edge = map[oedges[j].edge];
          redges[nedges].ewgt = oedges[j].ewgt;
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
* This function simply enumerates the nodes in the graph
**************************************************************************/
void SimpleOrder(CoarseGraphType *graph, int *perm, int order)
{
  int i;

  for (i=0; i<graph->nvtxs; i++)
    perm[graph->label[i]] = order--;

}


/*************************************************************************
* This function simply enumerates the nodes in the graph
**************************************************************************/
void MDOrder(CoarseGraphType *graph, int *perm, int order, SepNodeType *stree, int mysnode)
{
  int i, j;
  int *xadj, *adjncy;
  int *p, *ip, *head, *qsize, *list, *marker; 
  int nofsub;

  stree[mysnode].nvtxs = graph->nvtxs;
  stree[mysnode].hi = order;
  stree[mysnode].li = order - graph->nvtxs;
  stree[mysnode].isleaf = 1;

  if (graph->nvtxs == 0)
    return;

  adjncy = icoremalloc(graph->nedges, "MDOrder: adjncy", 0);
  xadj = icoremalloc(graph->nvtxs+5, "MDOrder: xadj", 0);
  p = icoremalloc(graph->nvtxs+5, "MDOrder: p", 0);
  ip = icoremalloc(graph->nvtxs+5, "MDOrder: ip", 0);
  head = icoremalloc(graph->nvtxs+5, "MDOrder: head", 0);
  qsize = icoremalloc(graph->nvtxs+5, "MDOrder: qsize", 0);
  list = icoremalloc(graph->nvtxs+5, "MDOrder: list", 0);
  marker = icoremalloc(graph->nvtxs+5, "MDOrder: marker", 0);

  xadj[0] = 1;
  for (i=0; i<graph->nvtxs; i++) {
    for (j=0; j<graph->vtxs[i]->nedges; j++) 
      adjncy[xadj[i]+j-1] = graph->vtxs[i]->edges[j].edge+1;
    xadj[i+1] = xadj[i]+graph->vtxs[i]->nedges;
  }

  for (i=0; i<graph->nvtxs; i++)
    p[i] = ip[i] = i+1;

  genmmd(graph->nvtxs, xadj, adjncy, p, ip, 1, head, qsize, list, marker, 1000000, &nofsub);

  order = order - graph->nvtxs + 1;
  for (i=0; i<graph->nvtxs; i++)
    perm[graph->label[i]] = order + p[i]-1;

  icorefree(7*(graph->nvtxs+5) + graph->nedges);
}
