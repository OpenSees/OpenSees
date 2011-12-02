/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * util.c
 *
 * This file contains utility functions
 *
 * Started 8/28/94
 * George
 *
 * $Id: util.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include "multilevel.h"


/*************************************************************************
* This function initializes the random number generator
**************************************************************************/
void InitRandom(void)
{
  srand48(7654321L);  
}


/*************************************************************************
* This function takes a graph and randomly permutes the adjacency list
* of each vertex.
**************************************************************************/
void PermuteGraphRandom(CoarseGraphType *graph)
{
  int i, j, k;
  VertexType *vtx;
  int u, v;
  EdgeType tmp, *edges;

  edges = graph->vtxs[0]->edges;
  for (i=0; i<graph->nvtxs; i++) {
    vtx = graph->vtxs[i];
    k = vtx->nedges;

    for (j=k/2; j>0; j--) {
      u = RandomInRange(k);
      v = RandomInRange(k);
      SWAP(edges[u], edges[v], tmp);
    }
    edges += k;
  }
}



/*************************************************************************
* This file randomly permutes the contents of an array.
* flag == 0, don't initialize perm
* flag == 1, set p[i] = i 
**************************************************************************/
void RandomPermute(int *p, int n, int flag)
{
  int i, u, v, tmp, m;

  if (flag == 1)
    for (i=0; i<n; i++)
      p[i] = i;

  m = n-8;
  for (i=0; i<n/16; i++) {
    v = RandomInRange(m);
    u = RandomInRange(m);
    SWAP(p[v], p[u], tmp);
    SWAP(p[v+1], p[u+1], tmp);
    SWAP(p[v+2], p[u+2], tmp);
    SWAP(p[v+3], p[u+3], tmp);
  }
}



/*************************************************************************
* This function returns the relative difference of two ints
**************************************************************************/
float RelDiff(int a, int b)
{
  double da = (double)a;
  double db = (double)b;
  float rdiff;

  if (da > db)
    rdiff = fabs(da-db)/da;
  else
    rdiff = fabs(da-db)/db;

  return rdiff;

}


/*************************************************************************
* This function checks the id/ed of a graph and the ewgtsum
**************************************************************************/
int CheckDegrees(CoarseGraphType *graph)
{
  int i, j, l;
  int id, ed;
  VertexType *vtx;

  for (i=0; i<graph->nvtxs; i++) {
    id = ed = 0;
    vtx = graph->vtxs[i];
    for (j=0; j<vtx->nedges; j++) {
      if (graph->where[i] == graph->where[vtx->edges[j].edge])
        id += vtx->edges[j].ewgt;
      else
        ed += vtx->edges[j].ewgt;
    }
    if (id != graph->id[i])
      printf("\nID for %d failed: %d %d",i, graph->id[i], id);
    if (ed != graph->ed[i])
      printf("\nED for %d failed: %d %d",i, graph->ed[i], ed);

    l = 0;
    for (j=0; j<vtx->nedges; j++)
      l += vtx->edges[j].ewgt;
    if (vtx->ewgtsum != l)
      printf("\nEWGTSUM for %d failed: %d %d",i, vtx->ewgtsum, l);
  }

  return 1;
}



/*************************************************************************
* This function sorts an array of type KeyValueType in increasing order
**************************************************************************/
void SortKeyValueNodesDec(KeyValueType *nodes, int n)
{
  qsort((void *)nodes, (size_t)n, (size_t)sizeof(KeyValueType), DecKeyValueCmp);
}


/*************************************************************************
* This function compares 2 KeyValueType variables for sorting in inc order
**************************************************************************/
int DecKeyValueCmp(const void *v1, const void *v2)
{
  KeyValueType *n1, *n2;

  n1 = (KeyValueType *)v1;
  n2 = (KeyValueType *)v2;

  return n2->key - n1->key;

}


/*
printmemstat()
{
  struct mallinfo mstat;

  mstat = mallinfo();
  printf("%d %d %d \n", mstat.arena, mstat.usmblks/4+mstat.uordblks/4, mstat.fsmblks/4+mstat.fordblks/4);

}
*/
