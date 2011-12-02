/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * match.c
 *
 * This file contains routines that perform the maximal graph matching
 *
 * Started 8/27/94
 * George
 *
 * $Id: match.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include "multilevel.h"

/*************************************************************************
* External Global variables 
**************************************************************************/
extern CtrlType *__Ctrl;		/* mlevelpart.c */

/*************************************************************************
* Global variables within match.c
**************************************************************************/
static int marker;	/* Used during vertex merging */
static KeyValueType *htable;	/* Index hash table used during vertex merging */
static int lastedge;	/* Used to mark the last edge position in free array */
static EdgeType *edges;	/* Pointer to the edge memory array */


/*************************************************************************
* This function finds a maximal matching for a graph randomly.
**************************************************************************/
int RM_Match(CoarseGraphType *graph)
{
  int i, j, k, ii;
  VertexType *cvtx;		/* The current vertex */
  int nvtxs;			/* The number of vertices in coersed graph */
  int *perm;			/* Random permutation array */
  int *cmap, *match;
  int idx;

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "RandomMatch: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "RandomMatch: graph->cmap");
  perm = icoremalloc(graph->nvtxs, "RandomMatch: perm", 0);

  RandomPermute(perm, graph->nvtxs, 1);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      idx = -1;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        if (cmap[k] == UNMATCHED) {
          idx = k;
          break;
        }
      }
      if (idx == -1)
        idx = i;

      cmap[i] = cmap[idx] = nvtxs++;  
      match[i] = idx;
      match[idx] = i;
    }
  }

  icorefree(graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}


/*************************************************************************
* This function finds a maximal matching for a graph randomly.
**************************************************************************/
int RM_Match_W(CoarseGraphType *graph)
{
  int i, j, k, ii;
  VertexType *cvtx;		/* The current vertex */
  int nvtxs;			/* The number of vertices in coersed graph */
  int *perm;			/* Random permutation array */
  int *cmap, *match;
  int idx;
  int maxvwgt;

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "RandomMatch: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "RandomMatch: graph->cmap");
  perm = icoremalloc(graph->nvtxs, "RandomMatch: perm", 0);

  maxvwgt = graph->tvwgt/(NPARTS_FACTOR*__Ctrl->nparts);

  RandomPermute(perm, graph->nvtxs, 1);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      idx = -1;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        if (cmap[k] == UNMATCHED && cvtx->vwgt+graph->vtxs[k]->vwgt < maxvwgt) {
          idx = k;
          break;
        }
      }
      if (idx == -1)
        idx = i;

      cmap[i] = cmap[idx] = nvtxs++;  
      match[i] = idx;
      match[idx] = i;
    }
  }

  icorefree(graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}


/*************************************************************************
* This function finds a maximal matching for a graph by combining the
* random and heaby edge heuristics
**************************************************************************/
int HEM_Match(CoarseGraphType *graph)
{
  int i, j, k, kwgt, ii;
  VertexType *cvtx;		/* The current vertex */
  int nvtxs;			/* The number of vertices in coersed graph */
  int *perm;			/* Random permutation array */
  int maxewgt;		        /* The maximum edge-weight and the index */
  int maxidx;
  int *cmap, *match;

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "RandomHeavyEdgeMatch: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "RandomHeavyEdgeMatch: graph->match");
  perm = icoremalloc(graph->nvtxs, "RandomHeavyEdgeMatch: perm", 0);

  RandomPermute(perm, graph->nvtxs, 1);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      maxewgt = -1;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;
        if (cmap[k] == UNMATCHED) {
          if (maxewgt < kwgt) {
            maxidx = k;
            maxewgt = kwgt;
          }
        }
      }
      if (maxewgt == -1)
        maxidx = i;

      cmap[i] = cmap[maxidx] = nvtxs++;  
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  icorefree(graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}


/*************************************************************************
* This function finds a maximal matching for a graph by combining the
* random and heaby edge heuristics
**************************************************************************/
int HEM_Match_W(CoarseGraphType *graph)
{
  int i, j, k, kwgt, ii;
  VertexType *cvtx;		/* The current vertex */
  int nvtxs;			/* The number of vertices in coersed graph */
  int *perm;			/* Random permutation array */
  int maxewgt;		        /* The maximum edge-weight and the index */
  int maxidx;
  int *cmap, *match;
  int maxvwgt;

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "RandomHeavyEdgeMatch: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "RandomHeavyEdgeMatch: graph->match");
  perm = icoremalloc(graph->nvtxs, "RandomHeavyEdgeMatch: perm", 0);

  maxvwgt = graph->tvwgt/(NPARTS_FACTOR*__Ctrl->nparts);

  RandomPermute(perm, graph->nvtxs, 1);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      maxewgt = -1;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;
        if (cmap[k] == UNMATCHED && cvtx->vwgt+graph->vtxs[k]->vwgt < maxvwgt) {
          if (maxewgt < kwgt) {
            maxidx = k;
            maxewgt = kwgt;
          }
        }
      }
      if (maxewgt == -1)
        maxidx = i;

      cmap[i] = cmap[maxidx] = nvtxs++;  
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  icorefree(graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}




/*************************************************************************
* This function finds a maximal matching for a graph by combining the
* random and light edge heuristics
**************************************************************************/
int LEM_Match(CoarseGraphType *graph)
{
  int i, j, k, kwgt, ii;
  VertexType *cvtx;		/* The current vertex */
  int nvtxs;			/* The number of vertices in coersed graph */
  int *perm;			/* Random permutation array */
  int minewgt;		        /* The maximum edge-weight and the index */
  int minidx;
  int *cmap, *match;

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "RandomLightEdgeMatch: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "RandomLightEdgeMatch: graph->match");
  perm = icoremalloc(graph->nvtxs, "RandomHeavyEdgeMatch: perm", 0);

  RandomPermute(perm, graph->nvtxs, 1);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      minewgt = 1000000;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;
        if (cmap[k] == UNMATCHED) {
          if (minewgt > kwgt) {
            minidx = k;
            minewgt = kwgt;
          }
        }
      }
      if (minewgt == 1000000)
        minidx = i;

      cmap[i] = cmap[minidx] = nvtxs++;  
      match[minidx] = i;
      match[i] = minidx;
    }
  }

  icorefree(graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}



/*************************************************************************
* This function finds a maximal matching for a graph by combining the
* random and light edge heuristics
**************************************************************************/
int LEM_Match_W(CoarseGraphType *graph)
{
  int i, j, k, kwgt, ii;
  VertexType *cvtx;		/* The current vertex */
  int nvtxs;			/* The number of vertices in coersed graph */
  int *perm;			/* Random permutation array */
  int minewgt;		        /* The maximum edge-weight and the index */
  int minidx, maxvwgt;
  int *cmap, *match;

  maxvwgt = graph->tvwgt/(NPARTS_FACTOR*__Ctrl->nparts);

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "RandomLightEdgeMatch: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "RandomLightEdgeMatch: graph->match");
  perm = icoremalloc(graph->nvtxs, "RandomHeavyEdgeMatch: perm", 0);

  RandomPermute(perm, graph->nvtxs, 1);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      minewgt = 1000000;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;
        if (cmap[k] == UNMATCHED && cvtx->vwgt+graph->vtxs[k]->vwgt < maxvwgt) {
          if (minewgt > kwgt) {
            minidx = k;
            minewgt = kwgt;
          }
        }
      }
      if (minewgt == 1000000)
        minidx = i;

      cmap[i] = cmap[minidx] = nvtxs++;  
      match[minidx] = i;
      match[i] = minidx;
    }
  }

  icorefree(graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}




/*************************************************************************
* This function finds a maximal matching for a graph by combining the
* random and maximu clique edge heuristics
**************************************************************************/
int HCM_Match(CoarseGraphType *graph)
{
  int i, j, k, kwgt, ii;
  VertexType *cvtx;		/* The current vertex */
  int nvtxs;			/* The number of vertices in coersed graph */
  int *perm;			/* Random permutation array */
  float maxclique;	        /* The maximum edge-weight and the index */
  int maxidx;
  float tmp;
  int *cmap, *match;

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "RandomHeavyEdgeMatch: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "RandomHeavyEdgeMatch: graph->match");
  perm = icoremalloc(graph->nvtxs, "RandomHeavyEdgeMatch: perm", 0);

  RandomPermute(perm, graph->nvtxs, 1);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      maxclique = -1;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;
        if (cmap[k] == UNMATCHED) {
          tmp = ((float)(cvtx->cewgt + graph->vtxs[k]->cewgt + 2.0*kwgt)) /
                ((float)(cvtx->vwgt + graph->vtxs[k]->vwgt));
          if (maxclique < tmp) {
            maxidx = k;
            maxclique = tmp;
          }
        }
      }
      if (maxclique == -1)
        maxidx = i;

      cmap[i] = cmap[maxidx] = nvtxs++;  
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  icorefree(graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}



/*************************************************************************
* This function finds a maximal matching for a graph by combining the
* random and maximu clique edge heuristics
**************************************************************************/
int HCM_Match_W(CoarseGraphType *graph)
{
  int i, j, k, kwgt, ii;
  VertexType *cvtx;		/* The current vertex */
  int nvtxs;			/* The number of vertices in coersed graph */
  int *perm;			/* Random permutation array */
  float maxclique;	        /* The maximum edge-weight and the index */
  int maxidx, maxvwgt;
  float tmp;
  int *cmap, *match;

  maxvwgt = graph->tvwgt/(NPARTS_FACTOR*__Ctrl->nparts);

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "RandomHeavyEdgeMatch: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "RandomHeavyEdgeMatch: graph->match");
  perm = icoremalloc(graph->nvtxs, "RandomHeavyEdgeMatch: perm", 0);

  RandomPermute(perm, graph->nvtxs, 1);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      maxclique = -1;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;
        if (cmap[k] == UNMATCHED && cvtx->vwgt+graph->vtxs[k]->vwgt < maxvwgt) {
          tmp = ((float)(cvtx->cewgt + graph->vtxs[k]->cewgt + 2.0*kwgt)) /
                ((float)(cvtx->vwgt + graph->vtxs[k]->vwgt));
          if (maxclique < tmp) {
            maxidx = k;
            maxclique = tmp;
          }
        }
      }
      if (maxclique == -1)
        maxidx = i;

      cmap[i] = cmap[maxidx] = nvtxs++;  
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  icorefree(graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}




/*************************************************************************
* This function implements the fast version of Modified HEM
**************************************************************************/
int MHEM_Match(CoarseGraphType *graph)
{
  int i, ii, j, k, l, cewgt, kwgt, nedges;
  VertexType *cvtx;             /* The current vertex */
  EdgeType *edges;
  int nvtxs;                    /* The number of vertices in coersed graph */
  int *perm;                    /* Random permutation array */
  int maxewgt, maxaewgt;        /* The maximum edge-weight and the index */
  int maxidx;
  int *cmap, *match, *lookup;
  int *same, nsame;

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "MHEM_Match: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "MHEM_Match: graph->match");
  perm = icoremalloc(graph->nvtxs, "MHEM_Match: perm", 0);
  lookup = icoremalloc(graph->nvtxs, "MHEM_Match: lookup", 0);
  same = icoremalloc(1000, "MHEM_Match: same", 1);
  if (same == NULL) {
    icorefree(2*graph->nvtxs);
    return 0;
  }

  iset(graph->nvtxs, -1, lookup);

  RandomPermute(perm, graph->nvtxs, 1);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      nsame = maxewgt = maxaewgt = 0;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;

        if (cmap[k] == UNMATCHED) {
          if (maxewgt < kwgt) {
            maxidx = k;
            maxewgt = kwgt;
            same[0] = k;
            nsame = 1;
          }
          else { 
            if (maxewgt == kwgt) 
              same[nsame++] = k;
          }
        }
      }

      if (nsame > 2) {
        for (j=0; j<cvtx->nedges; j++)
          lookup[cvtx->edges[j].edge] = ii;

        for (j=0; j<nsame; j+=2) {
          k = same[j];
          edges = graph->vtxs[k]->edges;
          nedges = graph->vtxs[k]->nedges;

          cewgt = 0;
          for (l=0; l<nedges; l++) {
            if (lookup[edges[l].edge] == ii)
              cewgt += edges[l].ewgt; 
          }
          if (maxaewgt < cewgt) {
            maxidx = k;
            maxaewgt = cewgt;
          }
        }
      }

      if (maxewgt == 0)
        maxidx = i;

      cmap[i] = cmap[maxidx] = nvtxs++;  
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  icorefree(1000+2*graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}



/*************************************************************************
* This function implements the fast version of Modified HEM
**************************************************************************/
int MHEM_Match_W(CoarseGraphType *graph)
{
  int i, ii, j, k, l, cewgt, kwgt, nedges;
  VertexType *cvtx;             /* The current vertex */
  EdgeType *edges;
  int nvtxs;                    /* The number of vertices in coersed graph */
  int *perm;                    /* Random permutation array */
  int maxewgt, maxaewgt;        /* The maximum edge-weight and the index */
  int maxidx;
  int *cmap, *match, *lookup;
  int maxvwgt;
  int *same, nsame;

  maxvwgt = graph->tvwgt/(NPARTS_FACTOR*__Ctrl->nparts);

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "MHEM_Match: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "MHEM_Match: graph->match");
  perm = icoremalloc(graph->nvtxs, "MHEM_Match: perm", 0);
  lookup = icoremalloc(graph->nvtxs, "MHEM_Match: lookup", 0);
  same = icoremalloc(1000, "MHEM_Match: same", 1);
  if (same == NULL) {
    icorefree(2*graph->nvtxs);
    return 0;
  }

  iset(graph->nvtxs, -1, lookup);

  RandomPermute(perm, graph->nvtxs, 1);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      nsame = maxewgt = maxaewgt = 0;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;

        if (cmap[k] == UNMATCHED) {
          if (maxewgt < kwgt && (maxvwgt > cvtx->vwgt+graph->vtxs[k]->vwgt)) {
            maxidx = k;
            maxewgt = kwgt;
            same[0] = k;
            nsame = 1;
          }
          else { 
            if (maxewgt == kwgt && (maxvwgt > cvtx->vwgt+graph->vtxs[k]->vwgt)) 
              same[nsame++] = k;
          }
        }
      }

      if (nsame > 2) {
        for (j=0; j<cvtx->nedges; j++)
          lookup[cvtx->edges[j].edge] = ii;

        for (j=0; j<nsame; j+=2) {
          k = same[j];
          edges = graph->vtxs[k]->edges;
          nedges = graph->vtxs[k]->nedges;

          cewgt = 0;
          for (l=0; l<nedges; l++) {
            if (lookup[edges[l].edge] == ii)
              cewgt += edges[l].ewgt; 
          }
          if (maxaewgt < cewgt) {
            maxidx = k;
            maxaewgt = cewgt;
          }
        }
      }

      if (maxewgt == 0)
        maxidx = i;

      cmap[i] = cmap[maxidx] = nvtxs++;  
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  icorefree(1000+2*graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}



/*************************************************************************
* This function implements the Sorted RM matching
**************************************************************************/
int SRM_Match(CoarseGraphType *graph)
{
  int i, j, k, kwgt, ii;
  VertexType *cvtx;		/* The current vertex */
  int nvtxs;			/* The number of vertices in coersed graph */
  int maxidx, maxvwgt;
  int *cmap, *match, *perm;
  KeyValueType *cand;

  maxvwgt = graph->tvwgt/(NPARTS_FACTOR*__Ctrl->nparts);

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "SRM_Match: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "SRM_Match: graph->match");
  perm = icoremalloc(graph->nvtxs, "SRM_Match: perm", 0);
  cand = (KeyValueType *)icoremalloc(2*graph->nvtxs, "SRM_Match: cand", 1);
  if (cand == NULL) {
    icorefree(graph->nvtxs);
    return 0;
  }

  RandomPermute(perm, graph->nvtxs, 1);
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    cand[i].key = -graph->vtxs[i]->nedges;
    cand[i].val = i;
  }
  SortKeyValueNodesDec(cand, graph->nvtxs);


  for (i=0; i<graph->nvtxs; i++) {
    cand[i].key = -graph->vtxs[i]->nedges;
    cand[i].val = i;
  }
  SortKeyValueNodesDec(cand, graph->nvtxs);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = cand[ii].val;
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      maxidx = -1;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;
        if (cmap[k] == UNMATCHED && cvtx->vwgt+graph->vtxs[k]->vwgt < maxvwgt) {
          maxidx = k;
          break;
        }
      }
      if (maxidx == -1)
        maxidx = i;

      cmap[i] = cmap[maxidx] = nvtxs++;  
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  icorefree(3*graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}


/*************************************************************************
* This function implements the Sorted HEM matching scheme
**************************************************************************/
int SHEM_Match(CoarseGraphType *graph)
{
  int i, j, k, kwgt, ii;
  VertexType *cvtx;		/* The current vertex */
  int nvtxs;			/* The number of vertices in coersed graph */
  int maxewgt;		        /* The maximum edge-weight and the index */
  int maxidx, maxvwgt;
  int *cmap, *match, *perm;
  KeyValueType *cand;

  maxvwgt = graph->tvwgt/(NPARTS_FACTOR*__Ctrl->nparts);

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "SHEM_Match: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "SHEM_Match: graph->match");
  perm = icoremalloc(graph->nvtxs, "SHEM_Match: perm", 0);
  cand = (KeyValueType *)icoremalloc(2*graph->nvtxs, "SHEM_Match: cand", 1);
  if (cand == NULL) {
    icorefree(graph->nvtxs);
    return 0;
  }

  RandomPermute(perm, graph->nvtxs, 1);
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    cand[i].key = -graph->vtxs[i]->nedges;
    cand[i].val = i;
  }
  SortKeyValueNodesDec(cand, graph->nvtxs);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = cand[ii].val;
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      maxewgt = -1;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;
        if (cmap[k] == UNMATCHED && cvtx->vwgt+graph->vtxs[k]->vwgt < maxvwgt) {
          if (maxewgt < kwgt) {
            maxidx = k;
            maxewgt = kwgt;
          }
        }
      }
      if (maxewgt == -1)
        maxidx = i;

      cmap[i] = cmap[maxidx] = nvtxs++;  
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  icorefree(3*graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}




/*************************************************************************
* This function implements the fast version of Sorted Modified HEM
**************************************************************************/
int SMHEM_Match(CoarseGraphType *graph)
{
  int i, ii, j, k, l, cewgt, kwgt, nedges;
  VertexType *cvtx;             /* The current vertex */
  EdgeType *edges;
  int nvtxs;                    /* The number of vertices in coersed graph */
  int *perm;                    /* Random permutation array */
  int maxewgt, maxaewgt;        /* The maximum edge-weight and the index */
  int maxidx;
  int *cmap, *match, *lookup;
  int maxvwgt;
  int *same, nsame;
  KeyValueType *cand;

  maxvwgt = graph->tvwgt/(NPARTS_FACTOR*__Ctrl->nparts);

  cmap = graph->cmap = ismalloc(graph->nvtxs, UNMATCHED, "SMHEM_Match: graph->cmap");
  match = graph->match = imalloc(graph->nvtxs, "SMHEM_Match: graph->match");
  perm = icoremalloc(graph->nvtxs, "SMHEM_Match: perm", 0);
  lookup = icoremalloc(graph->nvtxs, "SMHEM_Match: lookup", 0);
  same = icoremalloc(1000, "SMHEM_Match: same", 1);
  cand = (KeyValueType *)icoremalloc(2*graph->nvtxs, "SMHEM_Match: cand", 1);
  if (same == NULL || cand == NULL) {
    icorefree(2*graph->nvtxs);
    return 0;
  }

  iset(graph->nvtxs, -1, lookup);

  RandomPermute(perm, graph->nvtxs, 1);
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = perm[ii];
    cand[i].key = -graph->vtxs[i]->nedges;
    cand[i].val = i;
  }
  SortKeyValueNodesDec(cand, graph->nvtxs);

  nvtxs = 0;
  for (ii=0; ii<graph->nvtxs; ii++) {
    i = cand[ii].val;
    if (cmap[i] == UNMATCHED) {
      cvtx = graph->vtxs[i];
      nsame = maxewgt = maxaewgt = 0;

      for (j=0; j<cvtx->nedges; j++) {
        k = cvtx->edges[j].edge;
        kwgt = cvtx->edges[j].ewgt;

        if (cmap[k] == UNMATCHED) {
          if (maxewgt < kwgt && (maxvwgt > cvtx->vwgt+graph->vtxs[k]->vwgt)) {
            maxidx = k;
            maxewgt = kwgt;
            same[0] = k;
            nsame = 1;
          }
          else { 
            if (maxewgt == kwgt && (maxvwgt > cvtx->vwgt+graph->vtxs[k]->vwgt)) 
              same[nsame++] = k;
          }
        }
      }

      if (nsame > 2) {
        for (j=0; j<cvtx->nedges; j++)
          lookup[cvtx->edges[j].edge] = ii;

        for (j=0; j<nsame; j+=2) {
          k = same[j];
          edges = graph->vtxs[k]->edges;
          nedges = graph->vtxs[k]->nedges;

          cewgt = 0;
          for (l=0; l<nedges; l++) {
            if (lookup[edges[l].edge] == ii)
              cewgt += edges[l].ewgt; 
          }
          if (maxaewgt < cewgt) {
            maxidx = k;
            maxaewgt = cewgt;
          }
        }
      }

      if (maxewgt == 0)
        maxidx = i;

      cmap[i] = cmap[maxidx] = nvtxs++;  
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  icorefree(1000+4*graph->nvtxs);

  return CreateCoarseGraph(graph, nvtxs);

}




/*************************************************************************
* This function creates the coarser graph
**************************************************************************/
int CreateCoarseGraph(CoarseGraphType *graph, int nvtxs)
{
  int i;
  CoarseGraphType *coarser;
  int *cmap, *match;
  int status;
  
  coarser = CreateGraph();
  coarser->nvtxs = nvtxs;

  coarser->vtxs = (VertexType **)GKmalloc(sizeof(VertexType *)*nvtxs, "CreateCoarseGraph: coarser->vtxs");
  coarser->allvtxs = (VertexType *)GKmalloc(sizeof(VertexType)*nvtxs, "CreateCoarseGraph: coarser->vtxs");
  for (i=0; i<nvtxs; i++)
    coarser->vtxs[i] = coarser->allvtxs+i;

  coarser->finer = graph;
  graph->coarser = coarser;
  cmap = graph->cmap;
  match = graph->match;

  htable = (KeyValueType *)GKmalloc(sizeof(KeyValueType)*graph->nvtxs, "RandomHeavyEdgeMatch: htable");
  for (i=0; i<graph->nvtxs; i++)
    htable[i].key = htable[i].val = 0;

  marker = 0;
  lastedge = 0;
  edges = GetEdgePool();

  for (i=0; i<graph->nvtxs; i++) {
    if (i <= match[i]) 
      mergevertices(coarser->vtxs, cmap, i, graph->vtxs[i], match[i], graph->vtxs[match[i]]);
  }

  free(htable);

  status = SetEdgePool(lastedge);
  if (status == 0) {
    FreeGraph(coarser);
    graph->coarser = NULL;
    GKfree(graph->cmap, graph->match, -1);
    graph->cmap = graph->match = NULL;
  }
  else {
    coarser->nedges = lastedge;
  }

  return status;
}



/*************************************************************************
* This function merges two vertices into a single vertex in the coarser
* graph.
**************************************************************************/
void mergevertices(VertexType **vtxs, int *cmap, int v, VertexType *vvtx, int u, VertexType *uvtx)
{
  int i, l, m;
  int vtx;	/* The single coarsed vertex */
  VertexType *newvtx;
  EdgeType *oldedges;

  vtx = cmap[v];
  newvtx = vtxs[vtx];

  newvtx->vwgt = vvtx->vwgt;
  newvtx->edges = edges + lastedge;

  m = lastedge;
  marker++;

  newvtx->cewgt = vvtx->cewgt;
  newvtx->ewgtsum = vvtx->ewgtsum;
  oldedges = vvtx->edges;
  for (i=0; i<vvtx->nedges; i++) {
    l = cmap[oldedges[i].edge];

    if (l == vtx) 
      continue;

    if (htable[l].key != marker) { /* I'm seeing this vertex for the first time */
      htable[l].key = marker;
      htable[l].val = m;
      edges[m].edge = l;
      edges[m].ewgt = oldedges[i].ewgt;
      m++;
    }
    else {
      edges[htable[l].val].ewgt += oldedges[i].ewgt;
    }
  }


  if (v != u) {
    newvtx->vwgt += uvtx->vwgt;
    newvtx->cewgt += uvtx->cewgt;
    newvtx->ewgtsum += uvtx->ewgtsum;

    oldedges = uvtx->edges;
    for (i=0; i<uvtx->nedges; i++) {
      l = cmap[oldedges[i].edge];

      if (l == vtx) { /* Add/Remove the compacted edge weight */
        INC_DEC(newvtx->cewgt, newvtx->ewgtsum, 2*oldedges[i].ewgt);
        continue;
      }

      if (htable[l].key == marker) {
        edges[htable[l].val].ewgt += oldedges[i].ewgt;
      }
      else { /* I'm seeing this vertex for the first time */
        htable[l].key = marker;
        htable[l].val = m;
        edges[m].edge = l;
        edges[m].ewgt = oldedges[i].ewgt;
        m++;
      }
    }
  }

  newvtx->nedges = m-lastedge;
  lastedge = m;
}

