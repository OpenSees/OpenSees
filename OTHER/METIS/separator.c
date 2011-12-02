/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * separator.c
 *
 * This file contains code that computes the node separator from a graph
 *
 * Started 10/28/94
 * George
 *
 */

#include "multilevel.h"

/*************************************************************************
* External Variables
**************************************************************************/
extern CtrlType *__Ctrl;		/* mlevelpart.c */




/*************************************************************************
* This function takes a graph, computes the smallest node separator
* based on mincover and modifies the graph->where to != (0,1) indicating 
* that this node is in the separator.
**************************************************************************/
void FindMinCovNodeSeparator(CoarseGraphType *graph, int *perm, int order)
{
  int i, j, k, l, nedges;
  int asize, bsize, ia, ib;
  VertexType **vtxs;
  EdgeType *edges;
  int *where, *ed;
  int *bxadj, *badjncy;
  int *vmap, *ivmap, *cover, csize, *bnd;


  if (graph->cmap == NULL)
    graph->cmap = imalloc(graph->nvtxs, "FindMinCovNodeSeparator: graph->cmap");
  if (graph->match == NULL)
    graph->match = imalloc(graph->nvtxs, "FindMinCovNodeSeparator: graph->match");

  vmap = graph->cmap;  		/* Use graph->cmap and graph->match as aux vectors */
  ivmap = graph->match;
  vtxs = graph->vtxs;
  where = graph->where;
  ed = graph->ed;
  bnd = graph->id;

  asize = bsize = nedges = 0;
  for (i=0; i<graph->nvtxs; i++) {
    if (ed[i] > 0) {
      bnd[bsize++] = i;
      if (where[i] == 0)
        asize++;
      nedges += ed[i];
    }
  }

  bxadj = imalloc(bsize+1, "FindMinCovSeparator: bxadj");
  badjncy = imalloc(nedges, "FindMinCovSeparator: badjncy");
  cover = imalloc(bsize, "FindMinCovSeparator: cover");

  ia = 0;
  ib = asize;
  for (i=0; i<bsize; i++) {
    j = bnd[i];
    k = (where[j] == 0 ? ia++ : ib++);
    vmap[j] = k;
    ivmap[k] = j;
    bxadj[k+1] = ed[j];
  }

  bxadj[0] = 0;
  for (i=2; i<bsize+1; i++)
    bxadj[i] += bxadj[i-1];

  for (i=0; i<bsize; i++) {
    k = bnd[i];
    l = bxadj[vmap[k]];
    edges = vtxs[k]->edges;
    for (j=0; j<vtxs[k]->nedges; j++) {
      if (where[edges[j].edge] != where[k]) 
        badjncy[l++] = vmap[edges[j].edge];
    }
  }

  if (__Ctrl->dbglvl&DBG_ORDERNODE)
    printf("Nvtxs: %6d, Cut: %6d, SS: [%4d %4d], ", graph->nvtxs, graph->mincut, asize, bsize-asize);

  MinCover(bxadj, badjncy, asize, bsize, cover, &csize);

  if (__Ctrl->dbglvl&DBG_ORDERNODE)
    printf("Cover: %5d\n", csize);

  for (i=0; i<csize; i++) {
    j = ivmap[cover[i]]; 
    where[j] += 2;
    perm[graph->label[j]] = order--;
  }

  GKfree(bxadj, badjncy, cover, -1);
}
