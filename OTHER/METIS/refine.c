/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * refine.c
 *
 * This file contains code for the incremental refinement of the coarsergraphs
 *
 * Started 8/31/94
 * George
 *
 * $Id: refine.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
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
extern timer InitPartTmr;	/* main.c */
extern timer ProjectTmr;	/* main.c */
extern timer UncrsTmr;		/* main.c */
#endif


/*************************************************************************
* This function is the driver of the iterative refinment
**************************************************************************/
void Refine(CoarseGraphType *orggraph, CoarseGraphType *graph, int zeropwgt)
{
  CoarseGraphType *fgraph;
  int limit = 0.02*orggraph->nvtxs;

  if (limit < 400)
    limit = 400;

  InitPartition(graph, zeropwgt);

  starttimer(&UncrsTmr);

  if (__Ctrl->OpType != OP_MLND)
    FastInitBalance(graph, zeropwgt);

  for (fgraph = graph; ; fgraph = fgraph->finer) {
    ASSERT(CheckBndSize(fgraph));

    if (__Ctrl->OpType != OP_MLND && __Ctrl->RefineType != REFINE_NONE) 
      if (orggraph->nvtxs/fgraph->nvtxs < 15) 
        FastBalance(fgraph, zeropwgt, 5*fgraph->tvwgt/fgraph->nvtxs);


    ASSERT(CheckBndSize(fgraph));

    starttimer(&GreedyTmr);
    switch (__Ctrl->RefineType) {
      case REFINE_GR:
        FMR_Refine(fgraph, zeropwgt, 1);
        break;
      case REFINE_KLR:
        FMR_Refine(fgraph, zeropwgt, 20);
        break;
      case REFINE_GKLR:
        if (fgraph->nvtxs < 2*limit)  /* Small graph, do 2*greedy */
          FMR_Refine(fgraph, zeropwgt, 20);
        else  /* All other graphs do single Greedy */
          FMR_Refine(fgraph, zeropwgt, 1);
        break;
      case REFINE_BGR:
        if (!__Ctrl->IsWeighted)
          BFMR_Refine(fgraph, zeropwgt, SMART, 2);
        else
          BFMR_Refine_Weighted(fgraph, zeropwgt, SMART, 2);
        break;
      case REFINE_BKLR:
        if (!__Ctrl->IsWeighted)
          BFMR_Refine(fgraph, zeropwgt, SMART, 20);
        else
          BFMR_Refine_Weighted(fgraph, zeropwgt, SMART, 20);
        break;
      case REFINE_BGKLR:
        if (fgraph->nbnd < limit) 
          if (!__Ctrl->IsWeighted)
            BFMR_Refine(fgraph, zeropwgt, SMART, 20);
          else
            BFMR_Refine_Weighted(fgraph, zeropwgt, SMART, 20);
        else 
          if (!__Ctrl->IsWeighted)
            BFMR_Refine(fgraph, zeropwgt, SMART, 2);
          else
            BFMR_Refine_Weighted(fgraph, zeropwgt, SMART, 2);
        break;
      case 14:
        if (fgraph->nbnd < limit) 
          BFMR_Refine(fgraph, zeropwgt, STUPID, 20);
        else
          BFMR_Refine(fgraph, zeropwgt, STUPID, 2);
        break;
      case 15:
        Greedy_Refine(fgraph, 10);
        break;
      case REFINE_NONE:
        break;
      default:
        errexit("Unsupported Refine Type: %d", __Ctrl->RefineType);
    }

    if (fgraph->finer == NULL && __Ctrl->OpType != OP_MLND && __Ctrl->RefineType >= 10) 
      BFMR_Refine_EqWgt(fgraph, zeropwgt);

    stoptimer(&GreedyTmr);

    if (fgraph != orggraph)
      ProjectPartition(fgraph, limit);
    else
      break;
  }

  stoptimer(&UncrsTmr);

}




/*************************************************************************
* This function computes the id/ed, the partition weights and the boundary
**************************************************************************/
void ComputePartitionParams(CoarseGraphType *graph)
{
  int i, j, nedges, me;
  int *id, *ed, *where, *pwgts;
  HTableType *htable;
  int cut;
  VertexType *vtx;
  EdgeType *edges;

  pwgts = graph->pwgts;
  where = graph->where;
  id = graph->id = ismalloc(graph->nvtxs, 0, "ComputePartitionParams: id");
  ed = graph->ed = ismalloc(graph->nvtxs, 0, "ComputePartitionParams: id");

  htable = &(graph->htable);
  CreateHTable(htable, graph->nvtxs);

  pwgts[0] = pwgts[1] = 0;
  cut = 0;
  for (i=0; i<graph->nvtxs; i++) {
    me = where[i];
    vtx = graph->vtxs[i];
    pwgts[me] += vtx->vwgt;
    nedges = vtx->nedges;
    edges = vtx->edges;
    for (j=0; j<nedges; j++) {
      if (me == where[edges[j].edge]) 
        id[i] += edges[j].ewgt;
      else
        ed[i] += edges[j].ewgt;
    }
    if (ed[i] > 0) {
      AddHTable(htable, i);
      cut += ed[i];
    }
  }

  graph->mincut = cut/2;
  graph->nbnd = htable->nelem;

  ASSERT(CheckBndSize(graph));
}



/*************************************************************************
* This function takes a graph with a partition and projects the partition
* to the finer graph
**************************************************************************/
void ProjectPartition(CoarseGraphType *cgraph, int limit)
{
  int i, v, u, mapped;
  CoarseGraphType *fgraph;
  int *fid, *fed, *cid, *ced, *fwhere, *cwhere;
  HTableType *fhtable;
  int me;
  int *cmap;
  VertexType *vtx;

  if (cgraph->finer == NULL)
    return;

  starttimer(&ProjectTmr);

  cwhere = cgraph->where;
  cid = cgraph->id;
  ced = cgraph->ed;

  fgraph = cgraph->finer;
  cmap = fgraph->cmap;
  fwhere = fgraph->where = imalloc(fgraph->nvtxs, "ProjectPartition: fgraph->where");
  fid = fgraph->id = imalloc(fgraph->nvtxs, "ProjectPartition1: fgraph->id");
  fed = fgraph->ed = imalloc(fgraph->nvtxs, "ProjectPartition1: fgraph->ed");

  fhtable = &(fgraph->htable);
  CreateHTable(fhtable, SelHTSize(2*cgraph->nbnd));

  for (v=0; v<fgraph->nvtxs; v++) {
    if (v <= fgraph->match[v]) {
      u = fgraph->match[v];

      mapped = cmap[v];
      fwhere[u] = fwhere[v] = cwhere[mapped];

      if (u == v) {
        fid[v] = cid[mapped];
        fed[v] = ced[mapped];
        if (fed[v] > 0) 
          AddHTable(fhtable, v);
      }
      else {
        if (ced[mapped] == 0) {
          fed[u] = fed[v] = 0;
          fid[u] = fgraph->vtxs[u]->ewgtsum;
          fid[v] = fgraph->vtxs[v]->ewgtsum;
        }
        else if (cid[mapped] == 0) {
          fid[u] = fid[v] = (fgraph->vtxs[u]->ewgtsum + fgraph->vtxs[v]->ewgtsum - cgraph->vtxs[mapped]->ewgtsum)>>1;
          fed[u] = fgraph->vtxs[u]->ewgtsum - fid[u];
          fed[v] = fgraph->vtxs[v]->ewgtsum - fid[v];
          if (fed[v] > 0) 
            AddHTable(fhtable, v);
          if (fed[u] > 0) 
            AddHTable(fhtable, u);
        }
        else { /* Time to do some work */
          fid[u] = fed[u] = 0;
          me = fwhere[u];

          vtx = fgraph->vtxs[u];
          for (i=0; i<vtx->nedges; i++) {
            if (me == cwhere[cmap[vtx->edges[i].edge]])
              fid[u] += vtx->edges[i].ewgt;
            else
              fed[u] += vtx->edges[i].ewgt;
          }

          fid[v] = cid[mapped] - fid[u] + 
                     (vtx->ewgtsum + fgraph->vtxs[v]->ewgtsum - cgraph->vtxs[mapped]->ewgtsum);
          fed[v] = ced[mapped] - fed[u];

          if (fed[v] > 0) 
            AddHTable(fhtable, v);
          if (fed[u] > 0) 
            AddHTable(fhtable, u);
        }
      }
    }
  }

  fgraph->nbnd = fhtable->nelem;

  fgraph->pwgts[0] = cgraph->pwgts[0];
  fgraph->pwgts[1] = cgraph->pwgts[1];
  fgraph->mincut = cgraph->mincut;
  fgraph->coarser = NULL;

  FreeGraph(cgraph);

  stoptimer(&ProjectTmr);

  ASSERT(CheckBndSize(fgraph));
}


