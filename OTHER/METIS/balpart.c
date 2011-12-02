/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * balpart.c
 *
 * This file contains code that does a fast balancing .
 *
 * Started 1/10/95
 * $Id: balpart.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include "multilevel.h"


/*************************************************************************
* External Variables
**************************************************************************/
extern CtrlType *__Ctrl;        /* mlevelpart.c */
#ifndef METISLIB
extern timer BalanceTmr;        /* main.c */
#endif


/*************************************************************************
* This function performs a fast refinement by swapping a group of vertices
* from one partition to the other
**************************************************************************/
void FastInitBalance(CoarseGraphType *graph, int zeropwgt)
{
  int i, j, k, kwgt;
  KeyValueType *cand[2];
  VertexType *vtx;
  HTableType *htable;
  int n[2], idx[2];
  int *id, *ed, *where, *pwgts;
  int part;
  int from, me, olded;
  int higain, halfsplit;
  int idealwgts[2];

  pwgts = graph->pwgts;

  idealwgts[0] = zeropwgt;
  idealwgts[1] = graph->tvwgt - idealwgts[0];

  halfsplit = (graph->tvwgt > 2*(zeropwgt-1) && graph->tvwgt < 2*(zeropwgt+1) ? 1 : 0);

  if (halfsplit) {
    if (((float)abs(pwgts[0]-pwgts[1]))/((float)graph->tvwgt) < 0.10)
      return;
  }
  else {
    if (((float)abs(pwgts[0]-idealwgts[0]) + abs(pwgts[1]-idealwgts[1]))/((float)graph->tvwgt) < 0.10)
      return;
  }

  starttimer(&BalanceTmr);

  id = graph->id;
  ed = graph->ed;
  htable = &(graph->htable);
  where = graph->where;

  cand[0] = (KeyValueType *)GKmalloc(sizeof(KeyValueType)*graph->nvtxs, "FastGroupSwap: cand1");
  cand[1] = (KeyValueType *)GKmalloc(sizeof(KeyValueType)*graph->nvtxs, "FastGroupSwap: cand1");

  n[0] = n[1] = 0;

  for (i=0; i<graph->nvtxs; i++) {
    part = where[i];
    cand[part][n[part]].key = ed[i] - id[i];
    cand[part][n[part]].val = i;
    n[part]++;
  }

  SortKeyValueNodesDec(cand[0], n[0]);
  SortKeyValueNodesDec(cand[1], n[1]);

  if (__Ctrl->dbglvl&DBG_BALCUT) {
    printf("[%7d %5d], [%7d %5d], ORG-CUT: %5d ",pwgts[0], n[0], pwgts[1], n[1], graph->mincut);
    fflush(stdout);
  }

  idx[0] = idx[1] = 0;
  for (;;) {
    if (halfsplit) {
      if (fabs(((double)pwgts[0]-pwgts[1])/((double)graph->tvwgt)) < 0.10)
        break;
      from = (pwgts[0] < pwgts[1] ? 1 : 0);
    }
    else {
      if (((float)abs(pwgts[0]-idealwgts[0]) + abs(pwgts[1]-idealwgts[1]))/((float)graph->tvwgt) < 0.10)
        break;
      from = (pwgts[0] < idealwgts[0] ? 1 : 0);
    }

    if (idx[from] == n[from])
      break;

    higain = cand[from][idx[from]].val;
    idx[from]++;

    INC_DEC(pwgts[(from+1)%2], pwgts[from], graph->vtxs[higain]->vwgt);
  }

  /* Perform the combined swap */
  for (from=0; from<2; from++) {
    for (i=0; i<idx[from]; i++) {
      higain = cand[from][i].val;

      me = where[higain] = (where[higain] + 1)%2;
      if (ed[higain] == 0 && id[higain] > 0) 
        AddHTable(htable, higain);
      else if (ed[higain] > 0 && id[higain] == 0) 
        DelHTable(htable, higain);
      SWAP(id[higain], ed[higain], j);

      graph->mincut += ed[higain] - id[higain];

      vtx = graph->vtxs[higain];

      for (j=0; j<vtx->nedges; j++) {
        k = vtx->edges[j].edge;
        olded = ed[k];

        kwgt = (me == where[k] ? vtx->edges[j].ewgt : -vtx->edges[j].ewgt);
        INC_DEC(id[k], ed[k], kwgt);

        if (olded == 0 && ed[k] > 0) 
          AddHTable(htable, k);
        else if (olded > 0 && ed[k] == 0) 
          DelHTable(htable, k);
      }
    }
  }
  stoptimer(&BalanceTmr);

  graph->nbnd = htable->nelem;

  if (__Ctrl->dbglvl&DBG_BALCUT) 
    printf(" [%7d %5d], [%7d %5d], FNL-CUT: %5d\n",pwgts[0], idx[0], pwgts[1], idx[1], graph->mincut);

  GKfree(cand[0], cand[1], -1);

}


/*************************************************************************
* This function performs a fast refinement by swapping a group of vertices
* from one partition to the other
**************************************************************************/
void FastBalance(CoarseGraphType *graph, int zeropwgt, int frac)
{
  int i, j, k, kwgt;
  KeyValueType *cand[2];
  VertexType *vtx;
  HTableType *htable;
  int n[2], idx[2];
  int *id, *ed, *where, *pwgts;
  int part;
  int from, me, olded;
  int higain, halfsplit;
  int diff, idealwgts[2];

  pwgts = graph->pwgts;

  idealwgts[0] = zeropwgt;
  idealwgts[1] = graph->tvwgt - idealwgts[0];

  halfsplit = (graph->tvwgt > 2*(zeropwgt-1) && graph->tvwgt < 2*(zeropwgt+1) ? 1 : 0);

  if (halfsplit)
    diff = abs(pwgts[0] - pwgts[1]);
  else
    diff = abs(pwgts[0]-idealwgts[0]) + abs(pwgts[1]-idealwgts[1]);

  if (((float)diff)/((float)graph->tvwgt) > 0.10)
    frac = frac/2;
  if (diff <= frac)
    return;

  starttimer(&BalanceTmr);

  id = graph->id;
  ed = graph->ed;
  htable = &(graph->htable);
  where = graph->where;

  cand[0] = (KeyValueType *)GKmalloc(sizeof(KeyValueType)*graph->nvtxs, "FastGroupSwap: cand1");
  cand[1] = (KeyValueType *)GKmalloc(sizeof(KeyValueType)*graph->nvtxs, "FastGroupSwap: cand1");

  n[0] = n[1] = 0;

  for (i=0; i<graph->nvtxs; i++) {
    if (ed[i] - id[i] > -20*frac) { 
      part = where[i];
      cand[part][n[part]].key = ed[i] - id[i];
      cand[part][n[part]].val = i;
      n[part]++;
    }
  }

  SortKeyValueNodesDec(cand[0], n[0]);
  SortKeyValueNodesDec(cand[1], n[1]);

  if (__Ctrl->dbglvl&DBG_BALCUT) {
    printf("[%7d %5d], [%7d %5d], OC: %5d ",pwgts[0], n[0], pwgts[1], n[1], graph->mincut);
    fflush(stdout);
  }

  idx[0] = idx[1] = 0;

  for (;;) {
    if (halfsplit) {
      if (abs(pwgts[0]-pwgts[1]) <= frac)
        break;
      from = (pwgts[0] < pwgts[1] ? 1 : 0);
    }
    else {
      if (abs(pwgts[0]-idealwgts[0]) + abs(pwgts[1]-idealwgts[1]) <= frac)
        break;
      from = (pwgts[0] < idealwgts[0] ? 1 : 0);
    }

    if (idx[from] == n[from])
      break;

    higain = cand[from][idx[from]].val;
    idx[from]++;

    INC_DEC(pwgts[(from+1)%2], pwgts[from], graph->vtxs[higain]->vwgt);
  }

  /* Perform the combined swap */
  for (from=0; from<2; from++) {
    for (i=0; i<idx[from]; i++) {
      higain = cand[from][i].val;

      me = where[higain] = (where[higain] + 1)%2;
      if (ed[higain] == 0 && id[higain] > 0) 
        AddHTable(htable, higain);
      else if (ed[higain] > 0 && id[higain] == 0) 
        DelHTable(htable, higain);
      SWAP(id[higain], ed[higain], j);

      graph->mincut += ed[higain] - id[higain];

      vtx = graph->vtxs[higain];

      for (j=0; j<vtx->nedges; j++) {
        k = vtx->edges[j].edge;
        olded = ed[k];

        kwgt = (me == where[k] ? vtx->edges[j].ewgt : -vtx->edges[j].ewgt);
        INC_DEC(id[k], ed[k], kwgt);

        if (olded == 0 && ed[k] > 0) 
          AddHTable(htable, k);
        else if (olded > 0 && ed[k] == 0) 
          DelHTable(htable, k);
      }
    }
  }
  stoptimer(&BalanceTmr);

  graph->nbnd = htable->nelem;

  if (__Ctrl->dbglvl&DBG_BALCUT) 
    printf(" [%7d %5d], [%7d %5d], FC: %5d [B]\n",pwgts[0], idx[0], pwgts[1], idx[1], graph->mincut);

  free(cand[0]);
  free(cand[1]);

}


/*************************************************************************
* This function performs a fast refinement by swapping a group of vertices
* from one partition to the other
**************************************************************************/
void FastBalance2(CoarseGraphType *graph, int zeropwgt)
{
  int i, j, k, kwgt;
  KeyValueType *cand;
  VertexType *vtx;
  int ncand;
  int *id, *ed, *where, *pwgts;
  int from, to, olded;
  int higain, nswaps;


  if (graph->pwgts[0] == zeropwgt)
    return;

  starttimer(&BalanceTmr);

  pwgts = graph->pwgts;

  id = graph->id;
  ed = graph->ed;
  where = graph->where;

  if (pwgts[0] < zeropwgt) {
    from = 1;
    to = 0;
    nswaps = zeropwgt - pwgts[0];
  }
  else {
    from = 0;
    to = 1;
    nswaps = pwgts[0] - zeropwgt;
  }

  cand = (KeyValueType *)GKmalloc(sizeof(KeyValueType)*graph->nvtxs, "FastGroupSwap: cand");

  ncand = 0;
  for (i=0; i<graph->nvtxs; i++) {
    if (where[i] == from && ed[i] > 0) {
      cand[ncand].key = ed[i] - id[i];
      cand[ncand].val = i;
      ncand++;
    }
  }
  if (ncand < nswaps) {  /* Ugly but I have to do it */
    printf("**** Complicated Balancing Starting ****\n");
    ncand = 0;
    j = graph->nedges/(1.5*graph->nvtxs);
    for (i=0; i<graph->nvtxs; i++) {
      if (where[i] == from && ed[i]-id[i] > -j) {
        cand[ncand].key = ed[i] - id[i];
        cand[ncand].val = i;
        ncand++;
      }
    }
  }

  SortKeyValueNodesDec(cand, ncand);

  if (__Ctrl->dbglvl&DBG_BALCUT) {
    printf("[%7d %7d, %5d], OC: %5d ",pwgts[from], pwgts[to], ncand, graph->mincut);
    fflush(stdout);
  }

  /* Perform nswaps */
  nswaps = amin(nswaps, ncand);
  INC_DEC(pwgts[to], pwgts[from], nswaps);
  for (i=0; i<nswaps; i++) {
    higain = cand[i].val;

    where[higain] = to;
    SWAP(id[higain], ed[higain], j);

    graph->mincut += ed[higain] - id[higain];

    vtx = graph->vtxs[higain];

    for (j=0; j<vtx->nedges; j++) {
      k = vtx->edges[j].edge;
      olded = ed[k];

      kwgt = (to == where[k] ? vtx->edges[j].ewgt : -vtx->edges[j].ewgt);
      INC_DEC(id[k], ed[k], kwgt);
    }
  }
  stoptimer(&BalanceTmr);

  if (__Ctrl->dbglvl&DBG_BALCUT) 
    printf(" [%7d %7d %5d] FC: %5d [B]\n",pwgts[0], pwgts[1], nswaps, graph->mincut);

  free(cand);
}


