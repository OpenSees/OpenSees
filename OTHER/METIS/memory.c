/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * memory.c
 *
 * This file contains routines that deal with memory allocation
 *
 * Started 8/27/94
 * George
 *
 * $Id: memory.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include "multilevel.h"

/*************************************************************************
* External Global Variables
**************************************************************************/
extern CtrlType *__Ctrl;		/* mlevelpart.c */


/*************************************************************************
* This function allocate various pools of memory
**************************************************************************/
void AllocatePools(CtrlType *ctrl)
{

  /* Edge pool + icore pool */
  ctrl->edgepool = (EdgeType *)GKmalloc(sizeof(EdgeType)*ctrl->maxedges + sizeof(int)*ctrl->maxicore, "AllocatePools: ctrl->edgepool");
  ctrl->lastedge = 0;

  /* External Degree pool */
  ctrl->degrees = NULL;
  ctrl->lastdegree = 0;
  ctrl->maxdegrees = 0;

  /* Integer pool starts at the end of the edge pool and moves backwards */
  ctrl->icore = (int *)ctrl->edgepool;
  ctrl->icore += (sizeof(EdgeType)/sizeof(int))*ctrl->maxedges + ctrl->maxicore;
  ctrl->cicore = 0;

  /* Gain pool */
  ctrl->gaincore = (ListNodeType *)GKmalloc(ctrl->maxgain*sizeof(ListNodeType), "InitBucketPools: ctrl->gaincore");
  ctrl->cgain = 0;

  /* bucket pool */
  ctrl->gbcore = (GainBucketType *)GKmalloc(ctrl->maxbucket*sizeof(GainBucketType), "InitBucketPools: ctrl->gbcore");
  ctrl->cbucket = 0;

}

/*************************************************************************
* This function de-allocate various pools of memory
**************************************************************************/
void FreePools(CtrlType *ctrl)
{
  GKfree(ctrl->edgepool, ctrl->gaincore, ctrl->gbcore, -1);
  ctrl->maxedges = ctrl->lastedge = ctrl->cicore = ctrl->maxicore = 0;
  ctrl->cgain = ctrl->maxgain = ctrl->cbucket = ctrl->maxbucket = 0;
  ctrl->maxdegrees = ctrl->lastdegree = 0;
}


/*************************************************************************
* This function creates a CoarseGraphType data structure and initializes
* the various fields
**************************************************************************/
CoarseGraphType *CreateGraph(void)
{
  CoarseGraphType *graph;

  graph = (CoarseGraphType *)GKmalloc(sizeof(CoarseGraphType), "CreateCoarseGraph: graph");

  InitGraph(graph);

  return graph;
}


/*************************************************************************
* This function takes the head of a coarsened graph list and frees memory
* stored in its arrays. It also frees the memory for the edges/ewgts of
* the head of the list.
**************************************************************************/
void FreeRootGraph(CoarseGraphType *graph)
{
  GKfree(graph->vtxs, graph->cmap, graph->match,
         graph->where, graph->id, graph->ed,
         graph->htable.ht, 
         graph->label,
         graph->rinfo, graph->kpwgts,
         -1);

  graph->vtxs = NULL;
  graph->cmap = NULL;
  graph->match = NULL;
  graph->where = NULL;
  graph->id = NULL;
  graph->ed = NULL;
  graph->htable.ht = NULL;
  graph->label = NULL;
  graph->rinfo = NULL;
  graph->kpwgts = NULL;

  ASSERT(graph->coarser == NULL);
}


/*************************************************************************
* This function frees a graph minus its edges
**************************************************************************/
void FreeGraph(CoarseGraphType *graph)
{
  int i;

  FreeEdgePool(graph->nedges);

  GKfree(graph->vtxs, graph->allvtxs, graph->cmap, graph->match,
         graph->where, graph->id, graph->ed,
         graph->htable.ht, 
         graph->label,
         graph->rinfo, graph->kpwgts,
         -1);

  free(graph);
}


/*************************************************************************
* This function creates a CoarseGraphType data structure and initializes
* the various fields
**************************************************************************/
void InitGraph(CoarseGraphType *graph) 
{
  graph->tvwgt = -1;
  graph->nvtxs = -1;
  graph->nedges = -1;
  graph->allvtxs = NULL;
  graph->vtxs = NULL;
  graph->cmap = NULL;
  graph->match = NULL;
  graph->coarser = NULL;
  graph->finer = NULL;
  graph->level = -1;

  graph->where = NULL;
  graph->id = NULL;
  graph->ed = NULL;
  graph->pwgts[0] = graph->pwgts[0] = graph->mincut = -1;

  graph->kpwgts = NULL;
  graph->rinfo = NULL;

  graph->nbnd = -1;
  graph->htable.size = -1;
  graph->htable.nelem = -1;
  graph->htable.ht = NULL;

  graph->label = NULL;

}



/*************************************************************************
* This function initializes the space for edges and ewgts
**************************************************************************/
void ResetPools(void)
{
  __Ctrl->lastedge = 0;
}


/*************************************************************************
* This function returns a pointer to the current edge in edges
**************************************************************************/
EdgeType *GetEdgePool(void)
{
  return __Ctrl->edgepool + __Ctrl->lastedge;
}


/*************************************************************************
* This function sets the end of the edge pool
**************************************************************************/
int SetEdgePool(int added)
{
  __Ctrl->lastedge += added;

  if (__Ctrl->lastedge > __Ctrl->maxedges) { /* An error condition */
    if (2*(__Ctrl->lastedge-__Ctrl->maxedges) < (__Ctrl->maxicore-__Ctrl->cicore)) {
      /* We are at the free memory of icore, no need to exit */
      return 0;
    }
    errexit("\nedgepool and ewgtpool have been exhusted! Used %d Had %d (%d %d)", 
       __Ctrl->lastedge, __Ctrl->maxedges, __Ctrl->cicore, __Ctrl->maxicore);
  }
  return 1;
}

/*************************************************************************
* This function sets the end of the edge pool
**************************************************************************/
void FreeEdgePool(int freed)
{
  __Ctrl->lastedge -= freed;

  if (__Ctrl->lastedge < 0)
    errexit("You are freeing the edgepool faster than you fill it up!\n");
}


/*************************************************************************
* This function returns the number of edges left in the pool
**************************************************************************/
int EdgePoolSizeLeft(void)
{
  return __Ctrl->maxedges - __Ctrl->lastedge;
}


/*************************************************************************
* This function returns a pointer to k ext degrees
**************************************************************************/
EdgeType *GetnExtDegrees(int n)
{
  __Ctrl->lastdegree += n;
  if (__Ctrl->lastdegree > __Ctrl->maxdegrees)
    errexit("Run out of External Degrees: %d %d\n", __Ctrl->lastdegree, __Ctrl->maxdegrees);

  return __Ctrl->degrees + (__Ctrl->lastdegree-n);
}

/*************************************************************************
* This function initializes the space for edges and ewgts
**************************************************************************/
void ResetExtDegrees()
{
  __Ctrl->maxdegrees = __Ctrl->maxedges - __Ctrl->lastedge;
  __Ctrl->lastdegree = 0;
  __Ctrl->degrees = __Ctrl->edgepool + __Ctrl->lastedge;
}


/*************************************************************************
* This function allocates n words from icore
**************************************************************************/
int *icoremalloc(int n, char *msg, int flag)
{
  if (n > 2*__Ctrl->maxedges+__Ctrl->maxicore - (2*__Ctrl->lastedge+__Ctrl->cicore)) {
    if (flag)
      return NULL;
    else 
      errexit("Memory allocation failed for %s, requested %d, had %d", msg, n, __Ctrl->maxicore-__Ctrl->cicore);
  }

  __Ctrl->cicore += n;
  return (__Ctrl->icore - __Ctrl->cicore);
}

/*************************************************************************
* This function frees n words from icore
**************************************************************************/
void icorefree(int n)
{
  if (__Ctrl->cicore - n < 0)
    errexit("You are freeing icore faster than you are filling it up!");

  __Ctrl->cicore -= n;
}
