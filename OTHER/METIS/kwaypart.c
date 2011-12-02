/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * kwaypart.c
 *
 * This file contains code for the k-way partition
 *
 * Started 6/6/95
 * George
 *
 * $Id: kwaypart.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include "multilevel.h"

/*************************************************************************
* External Global Variables
**************************************************************************/
extern CtrlType *__Ctrl;   	/* mlevelpart.c */
#ifndef METISLIB
extern timer CoarsenTmr;	/* main.c */
extern timer GreedyTmr;		/* main.c */
extern timer GreedyInitTmr;	/* main.c */
extern timer GreedyIterTmr;	/* main.c */
extern timer GreedyWrapUpTmr;	/* main.c */
extern timer InitPartTmr;	/* main.c */
extern timer ProjectTmr;	/* main.c */
extern timer UncrsTmr;		/* main.c */
#endif


/*************************************************************************
* This function is the entry point of the multilevel graph partition
**************************************************************************/
int KWayPart(CoarseGraphType *graph, int nparts, int CoarsenTo, int MatchType, 
             int RefineType, int dbglvl, int IsWeighted, int *part, int *kpwgts)
{
  int i;
  CtrlType *ctrl, *oldctrl;
  CoarseGraphType *cgraph, *fgraph;
  int *cpart, *ccuts, *ckpwgts;
  timer Tmr1, Tmr2;

  ctrl = GKmalloc(sizeof(CtrlType), "MultilevelPart: ctrl");

  ctrl->CoarsenTo = CoarsenTo;
  ctrl->MatchType = MatchType;
  ctrl->InitPartType = -1;
  ctrl->RefineType = RefineType;
  ctrl->IsWeighted = IsWeighted;
  ctrl->OpType = OP_MLKP;
  ctrl->dbglvl = dbglvl;
  ctrl->nparts = nparts;
  ctrl->cfrac = KWAY_COARSEN_FRACTION;

  if (graph->nedges/graph->nvtxs > 20)
    ctrl->maxedges = 2.0*graph->nedges; /* Set the maxpool to something large */
  else
    ctrl->maxedges = 2.6*graph->nedges; /* Set the maxpool to something large */

  ctrl->maxicore = 2*graph->nvtxs;
  ctrl->maxgain = 2*graph->nvtxs;
  ctrl->maxbucket = 2*graph->nvtxs;

  AllocatePools(ctrl);

  /* Set current Ctrl context */
  oldctrl = __Ctrl;
  __Ctrl = ctrl;

  cleartimer(&Tmr1);
  starttimer(&Tmr1);
  cgraph = KwayCoarsen(graph, amax(CoarsenTo, NPARTS_FACTOR*nparts));
  stoptimer(&Tmr1);

  if (cgraph->finer != NULL) {
    i = amax(CoarsenTo, NPARTS_FACTOR*nparts);
    if (abs(cgraph->nvtxs-i) > abs(cgraph->finer->nvtxs-i)) {
      cgraph = cgraph->finer;
      FreeGraph(cgraph->coarser);
      cgraph->coarser = NULL;  
      GKfree(cgraph->cmap, cgraph->match, -1);
      cgraph->cmap = cgraph->match = NULL;
    }
  }


  cpart = ismalloc(cgraph->nvtxs, 0, "KWayPart: cpart");
  ckpwgts = ismalloc(nparts, 0, "KWayPart: ckpwgts");
  ccuts = ismalloc(cgraph->nvtxs, -1, "KWayPart: ccuts");

  cleartimer(&Tmr2);
  starttimer(&Tmr2);

  if (cgraph->label == NULL) {
    cgraph->label = imalloc(cgraph->nvtxs, "KWayPart: cgraph->label");
    for (i=0; i<cgraph->nvtxs; i++)
      cgraph->label[i] = i;
  }

  cgraph->mincut = MultiLevelPart(cgraph, nparts, 100, MATCH_HEM, INITPART_GGPKL, REFINE_BGKLR, dbglvl, 1, cpart, ccuts, ckpwgts);
  cgraph->where = cpart;	/* No need to free cpart anymore */

  stoptimer(&Tmr2);

  GKfree(ckpwgts, ccuts, -1);

  cleartimer(&CoarsenTmr);
  cleartimer(&GreedyTmr);
  cleartimer(&GreedyInitTmr);
  cleartimer(&GreedyIterTmr);
  cleartimer(&InitPartTmr);
  cleartimer(&ProjectTmr);
  cleartimer(&UncrsTmr);

#ifndef METISLIB
  CoarsenTmr = Tmr1;
  InitPartTmr = Tmr2;
#endif
  
  starttimer(&UncrsTmr);
  starttimer(&ProjectTmr);
  fgraph = cgraph->finer;
  if (fgraph != NULL) { /* Take care the no k-way refinement */
    KWayProjectPartition(cgraph);
    KWayComputePartitionParams(fgraph, nparts, 1);
  }
  stoptimer(&ProjectTmr);


  if (fgraph != NULL)   /* If fgraph = NULL, no k-way refinement */
    KWayRefine(graph, fgraph, nparts, kpwgts);

  icopy(graph->nvtxs, graph->where, part);

  stoptimer(&UncrsTmr);

  FreeRootGraph(graph);
  FreePools(ctrl);
  free(ctrl);

  __Ctrl = oldctrl;

  return graph->mincut;
}



/*************************************************************************
* This function is the driver of the k-way refinment
**************************************************************************/
void KWayRefine(CoarseGraphType *orggraph, CoarseGraphType *graph, int nparts, int *kpwgts)
{
  CoarseGraphType *fgraph;

  if (__Ctrl->dbglvl&DBG_PARTSIZES)
    printf("\n");

  fgraph = graph;
  while (1) {
    ASSERT(KWayCheckDegrees(fgraph));
    KWay_BalanceFM(fgraph, nparts, 1);
    ASSERT(KWayCheckDegrees(fgraph));

    starttimer(&GreedyTmr);
    switch (__Ctrl->RefineType) {
      case REFINE_BGR:
        KWay_RefineGreedy(fgraph, nparts, KWAY_REF_GREEDY_NITER);
        break;
      case REFINE_BKLR:
        KWay_RefineFM(fgraph, nparts, KWAY_REF_FM_NITER);
        break;
      case REFINE_NONE:
        break;
      default:
        errexit("Unsupported Refine Type: %d", __Ctrl->RefineType);
    }
    stoptimer(&GreedyTmr);

    ASSERT(KWayCheckDegrees(fgraph));

    fgraph = fgraph->finer;

    if (fgraph == NULL) 
      break;

    starttimer(&ProjectTmr);
    KWayProjectPartition(fgraph->coarser); 
    KWayComputePartitionParams(fgraph, nparts, 0);
    stoptimer(&ProjectTmr);
  }

  icopy(nparts, orggraph->kpwgts, kpwgts);

}




/*************************************************************************
* This function computes the id/ed, the partition weights and the boundary
**************************************************************************/
void KWayComputePartitionParams(CoarseGraphType *graph, int nparts, int flag)
{
  int i, j, k, nedges, me, other;
  register int id, ed;
  int *where, *kpwgts;
  RInfoType *rinfo, *myrinfo;
  VertexType *vtx;
  EdgeType *edges, *edegrees;

  ASSERT(graph->rinfo == NULL);

  graph->rinfo = (RInfoType *)GKmalloc(graph->nvtxs*sizeof(RInfoType), "KWayComputePartitionParams: rinfo");

  if (flag) {
    ASSERT(graph->kpwgts == NULL);
    graph->kpwgts = ismalloc(nparts, 0, "KWayComputePartitionParams: kpwgts");
  }

  kpwgts = graph->kpwgts;
  rinfo = graph->rinfo;
  where = graph->where;

  ResetExtDegrees();

  for (i=0; i<graph->nvtxs; i++) {
    me = where[i];
    vtx = graph->vtxs[i];
    myrinfo = rinfo+i;
    nedges = vtx->nedges;
    edges = vtx->edges;

    id = ed = 0;
    for (j=0; j<nedges; j++) {
      if (me == where[edges[j].edge]) 
        id += edges[j].ewgt;
      else
        ed += edges[j].ewgt;
    }

    myrinfo->id = id;
    myrinfo->ed = ed;
    myrinfo->ndegrees = 0;
    myrinfo->degrees = NULL;

    if (flag) 
      kpwgts[me] += vtx->vwgt;

    if (ed > 0) {  /* Time to do some serious work */
      edegrees = myrinfo->degrees = GetnExtDegrees(nedges);

      for (j=0; j<nedges; j++) {
        other = where[edges[j].edge];
        if (me != other) {
          for (k=0; k<myrinfo->ndegrees; k++) {
            if (edegrees[k].edge == other) {
              edegrees[k].ewgt += edges[j].ewgt;
              break;
            }
          }
          if (k == myrinfo->ndegrees) {
            edegrees[k].edge = other;
            edegrees[k].ewgt = edges[j].ewgt;
            myrinfo->ndegrees++;
          }
          ASSERT(myrinfo->ndegrees <= nedges);
        }
      }
    }
  }

}



/*************************************************************************
* This function takes a graph with a partition and projects the partition
* to the finer graph
**************************************************************************/
void KWayProjectPartition(CoarseGraphType *cgraph)
{
  int v, u, mapped;
  CoarseGraphType *fgraph;
  int *fwhere, *cwhere;
  int *cmap;

  cwhere = cgraph->where;

  fgraph = cgraph->finer;
  cmap = fgraph->cmap;

  ASSERT(fgraph->where == NULL);
  fwhere = fgraph->where = imalloc(fgraph->nvtxs, "ProjectPartition: fgraph->where");

  for (v=0; v<fgraph->nvtxs; v++) {
    if (v <= fgraph->match[v]) {
      u = fgraph->match[v];

      mapped = cmap[v];
      fwhere[u] = fwhere[v] = cwhere[mapped];
    }
  }

  fgraph->mincut = cgraph->mincut;
  fgraph->kpwgts = cgraph->kpwgts;
  cgraph->kpwgts = NULL;
  fgraph->coarser = NULL;

  FreeGraph(cgraph);
}



/*************************************************************************
* This function checks the id/ed degrees of a Kway graph
**************************************************************************/
int KWayCheckDegrees(CoarseGraphType *graph)
{
  int v, i, ed, id, *where, err=0;
  VertexType *vtx;
  RInfoType *rinfo;

  where = graph->where;
  rinfo = graph->rinfo;

  for (v=0; v<graph->nvtxs; v++) {
    vtx = graph->vtxs[v];
    for (i=0, id=0, ed=0; i<vtx->nedges; i++) {
      if (where[v] == where[vtx->edges[i].edge])
        id += vtx->edges[i].ewgt;
      else
        ed += vtx->edges[i].ewgt;
    }
    if (id != rinfo[v].id || ed != rinfo[v].ed) {
      printf("%d: %d %d %d %d\n", v, id, ed, rinfo[v].id, rinfo[v].ed);
      err = 1;
    }
  }

  if (err)
    return 0;
  else
    return 1;
}




