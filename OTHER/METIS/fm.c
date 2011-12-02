/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * fm.c
 *
 * This file contains functions that implement the Kernighan-Lin partition 
 * refinment algorithm.
 * This scheme uses an array of buckets to store gains for efficiency.
 * It also implements the FM algorithm instead of KL
 *
 * Started 9/10/94
 * George
 *
 * $Id: fm.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include "multilevel.h"

/*************************************************************************
* External Variables
**************************************************************************/
extern CtrlType *__Ctrl;	/* mlevelpart.c.c */
#ifndef METISLIB
extern timer GreedyTmr;		/* main.c */
extern timer GreedyInitTmr;	/* main.c */
extern timer GreedyIterTmr;	/* main.c */
extern timer GreedyWrapUpTmr;	/* main.c */
#endif


/*************************************************************************
* This function performs Kernighan-Lin refinement by swapping a vertex 
* at a time. Essentially implements the FM algorithm.
**************************************************************************/
void FMR_Refine(CoarseGraphType *graph, int zeropwgt, int npasses)
{
  int i, j;
  int k, kwgt;			/* The node and its adjacent edge weight */
  BucketListType parts[2]; 
  int *moved, *swaps, nswaps;
  int higain;			/* High gain node from a partition */
  int order;			/* The order at which nodes got moved */
  int mincut, newcut, mincutorder;	/* The minimu cut so far, the order at which it occured */
  int *id, *ed, *where;
  VertexType *vtx;
  int oldgain;
  int initcut;
  int pass=0;
  int *pwgts;                 /* The weights of the partitions */
  int from, me, limit, status;

  starttimer(&GreedyInitTmr);

  limit = amax(0.001*graph->nvtxs, 15);
  limit = amin(limit, 50);

  initbucket(&parts[0], graph->tvwgt, graph->nvtxs, graph->nvtxs, graph->level);
  initbucket(&parts[1], graph->tvwgt, graph->nvtxs, graph->nvtxs, graph->level);

  moved = icoremalloc(graph->nvtxs, "RefineKernLinEqualSingleVertex: moved", 0);
  swaps = icoremalloc(graph->nvtxs, "RefineKernLinEqualSingleVertex: moved", 0);

  id = graph->id;
  ed = graph->ed;
  where = graph->where;
  pwgts = graph->pwgts;

  stoptimer(&GreedyInitTmr);

  for (pass=0; pass<npasses; pass++) {
    starttimer(&GreedyInitTmr);

    resetbucket(&parts[0]);
    resetbucket(&parts[1]);

    status = 0;
    mincutorder = 0;
    initcut = mincut = newcut = graph->mincut;
    iset(graph->nvtxs, 0, moved);

    for (i=0; i<graph->nvtxs; i++) 
      Add2Part(&parts[where[i]], i, ed[i]-id[i]);

    stoptimer(&GreedyInitTmr);

    if (__Ctrl->dbglvl&DBG_FFCUT)
      if (pass == 0)
        printf("Partitions: [%6d,%8d]  [%6d,%8d], Initial Cut: %8d [%d]\n",
              parts[0].nnodes, pwgts[0], parts[1].nnodes, pwgts[1], initcut, ComputeCut(graph));


    /******************************************************
    * Get into the FM loop
    *******************************************************/
    starttimer(&GreedyIterTmr);
    order = nswaps = 0;
    for (;;) {
      order++;

      from = (pwgts[0] < zeropwgt ? 1 : 0);

      higain = GetMaxGainVtx(&parts[from]);

      if (higain == -1)
        break;

      newcut -= (ed[higain] - id[higain]);
      if (newcut <= mincut) {
        mincut = newcut;
        mincutorder = order;
      }
      else {
        if (order - mincutorder > limit)
          break; /* No further improvement, break out */
      }

      moved[higain] = order;
      swaps[nswaps++] = higain;  

      me = where[higain] = (from+1)%2;
      vtx = graph->vtxs[higain];
      INC_DEC(pwgts[me], pwgts[from], vtx->vwgt);

    
      if (__Ctrl->dbglvl&DBG_ITERCUT) {
        printf("Gains from moving %8d from %2d: %4d, \tNewcut: %d\n", higain, from, 
                ed[higain]-id[higain], newcut);
      }

      /**********************************************************
      * Update the id[i]/ed[i] values of the affected nodes
      ***********************************************************/
      SWAP(id[higain], ed[higain], j);
      for (j=0; j<vtx->nedges; j++) {
        k = vtx->edges[j].edge;
        oldgain = ed[k]-id[k];

        kwgt = (me == where[k] ? vtx->edges[j].ewgt : -vtx->edges[j].ewgt);
        INC_DEC(id[k], ed[k], kwgt);

        if (moved[k] == 0 && status == 0)
          status = UpdatePart(&parts[where[k]], k, oldgain, ed[k]-id[k]);
      }
      if (status == -1) {  /* Exit if there where problems in the UpdatePart */
        order++;
        break;
      }
    }
    stoptimer(&GreedyIterTmr);

    if (__Ctrl->dbglvl&DBG_FFCUT)
      printf("\tMinimum Cut: %8d at %5d [%6d %6d]\n",mincut, mincutorder, pwgts[0], pwgts[1]);


    /****************************************************************
    * Roll back computation 
    *****************************************************************/
    starttimer(&GreedyWrapUpTmr);
    for (nswaps--; nswaps>=0; nswaps--) {
      higain = swaps[nswaps];
      if (moved[higain] > mincutorder) {
        me = where[higain] = (where[higain]+1)%2;

        SWAP(id[higain], ed[higain], j);

        vtx = graph->vtxs[higain];
        INC_DEC(pwgts[me], pwgts[(me+1)%2], vtx->vwgt);
        for (j=0; j<vtx->nedges; j++) {
          k = vtx->edges[j].edge;

          kwgt = (me == where[k] ? vtx->edges[j].ewgt : -vtx->edges[j].ewgt);
          INC_DEC(id[k], ed[k], kwgt);
        }
      }
      else
        break;
    }
    stoptimer(&GreedyWrapUpTmr);

    graph->mincut = mincut;

    if (mincutorder == 0 || mincut >= initcut)
      break;
  }

  freebucket(&parts[0]);
  freebucket(&parts[1]);

  icorefree(2*graph->nvtxs);

}




/*************************************************************************
* This function performs Kernighan-Lin refinement by swapping a vertex 
* at a time. Essentially implements the FM algorithm.
**************************************************************************/
void BFMR_Refine(CoarseGraphType *graph, int zeropwgt, int smart, int npasses)
{
  int i, j;
  int k, kwgt;			/* The node and its adjacent edge weight */
  BucketListType parts[2]; 
  int *moved, *swaps, nswaps;
  int higain;			/* High gain node from a partition */
  int order;			/* The order at which nodes got moved */
  int mincut, mincutorder, newcut;	/* The minimum cut so far, the order at which it occured */
  int *id, *ed, *where;
  HTableType *htable;
  VertexType *vtx;
  int oldgain;
  int initcut;
  int pass=0;
  int *pwgts;                 /* The weights of the partitions */
  int from, me, olded, limit, status, halfsplit;

  ASSERT(CheckBndSize(graph));

  smart = (smart == SMART);

  TIMELVL(starttimer(&GreedyInitTmr));

  initbucket(&parts[0], graph->tvwgt, graph->nvtxs, graph->nvtxs, graph->level);
  initbucket(&parts[1], graph->tvwgt, graph->nvtxs, graph->nvtxs, graph->level);

  moved = icoremalloc(graph->nvtxs, "RefineKernLinEqualSingleVertex: moved", 0);
  swaps = icoremalloc(graph->nvtxs, "RefineKernLinEqualSingleVertex: swaps", 0);

  id = graph->id;
  ed = graph->ed;
  where = graph->where;
  htable = &(graph->htable);
  pwgts = graph->pwgts;

  halfsplit = (graph->tvwgt > 2*(zeropwgt-1) && graph->tvwgt < 2*(zeropwgt+1) ? 1 : 0);

  TIMELVL(stoptimer(&GreedyInitTmr));

  for (pass=0; pass<npasses; pass++) {
    TIMELVL(starttimer(&GreedyInitTmr));

    resetbucket(&parts[0]);
    resetbucket(&parts[1]);

    iset(graph->nvtxs, 0, moved);

    status = 0;
    mincutorder = 0;
    initcut = mincut = newcut = graph->mincut;
    limit = amax(0.01*htable->nelem, 15);
    limit = amin(limit, 50);

    for (j=0; j<htable->size; j++) {
      i = htable->ht[j];
      if (i != HT_EMPTY) {
        if (ed[i] > 0)  { /* Add only boundary vertices into the partition */
          Add2Part(&parts[where[i]], i, ed[i]-id[i]);
          moved[i] = -1;
        }
      }
    }
    TIMELVL(stoptimer(&GreedyInitTmr));

    if (__Ctrl->dbglvl&DBG_FFCUT)
      if (pass == 0)
        printf("Partitions: [%6d,%8d]  [%6d,%8d], BND: %5d, Initial Cut: %8d [%d] %d\n",
              parts[0].nnodes, pwgts[0], parts[1].nnodes, pwgts[1], htable->nelem, newcut, ComputeCut(graph), halfsplit);


    /******************************************************
    * Get into the FM loop
    *******************************************************/
    TIMELVL(starttimer(&GreedyIterTmr));
    order = nswaps = 0;
    for (;;) {
      order++;

      if (halfsplit)
        from = (pwgts[0] < pwgts[1] ? 1 : 0);
      else
        from = (pwgts[0] < zeropwgt ? 1 : 0);

      higain = GetMaxGainVtx(&parts[from]);

      if (higain == -1)
        break;

      newcut -= (ed[higain] - id[higain]);
      if (newcut <= mincut) {
        mincut = newcut;
        mincutorder = order;
      }
      else {
        if (order - mincutorder > limit)
          break; /* No further improvement, break out */
      }

      moved[higain] = order;
      swaps[nswaps++] = higain;

      me = where[higain] = (from+1)%2;
      vtx = graph->vtxs[higain];
      INC_DEC(pwgts[me], pwgts[from], vtx->vwgt);

      if (__Ctrl->dbglvl&DBG_ITERCUT) {
        printf("Gains from moving %8d from %2d: %4d, \tNewcut: %d [%d] %d %d %d\n", higain, from, 
                ed[higain]-id[higain], newcut, ComputeCut(graph), pwgts[0], pwgts[1], vtx->vwgt);
      }

      /**********************************************************
      * Update the id[i]/ed[i] values of the affected nodes
      ***********************************************************/
      if (ed[higain] == 0 && id[higain] > 0) 
        AddHTable(htable, higain);
      else if (ed[higain] > 0 && id[higain] == 0) 
        DelHTable(htable, higain);
      SWAP(id[higain], ed[higain], j);

      for (j=0; j<vtx->nedges; j++) {
        k = vtx->edges[j].edge;
        oldgain = ed[k]-id[k];
        olded = ed[k];

        kwgt = (me == where[k] ? vtx->edges[j].ewgt : -vtx->edges[j].ewgt); 
        INC_DEC(id[k], ed[k], kwgt);

        if (moved[k] == -1 && status == 0)  /* Update only boundary vertices */
          status = UpdatePart(&parts[where[k]], k, oldgain, ed[k]-id[k]);
        else if (smart && status == 0) {
          if (moved[k] == 0 && ed[k]-id[k] > 0) {
            status = Add2Part(&parts[where[k]], k, ed[k]-id[k]);
            moved[k] = -1;
          }
        }

        if (olded == 0 && ed[k] > 0) 
          AddHTable(htable, k);
        else if (olded > 0 && ed[k] == 0) 
          DelHTable(htable, k);
      }
      if (status == -1) {  /* Exit if there where problems in the UpdatePart */
        order++;
        break;
      }
    }
    TIMELVL(stoptimer(&GreedyIterTmr));

    if (__Ctrl->dbglvl&DBG_FFCUT)
      printf("\tMinimum Cut: %8d at %5d (%d) [%6d %6d]\n",mincut, mincutorder, order, pwgts[0], pwgts[1]);

    /****************************************************************
    * Roll back computation 
    *****************************************************************/
    TIMELVL(starttimer(&GreedyWrapUpTmr));
    for (nswaps--; nswaps>=0; nswaps--) {
      higain = swaps[nswaps];
      if (moved[higain] > mincutorder) {
        me = where[higain] = (where[higain]+1)%2;

        if (ed[higain] == 0 && id[higain] > 0) 
          AddHTable(htable, higain);
        else if (ed[higain] > 0 && id[higain] == 0) 
          DelHTable(htable, higain);
        SWAP(id[higain], ed[higain], j);
  
        vtx = graph->vtxs[higain];
        INC_DEC(pwgts[me], pwgts[(me+1)%2], vtx->vwgt);
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
      else
        break;
    }
    TIMELVL(stoptimer(&GreedyWrapUpTmr));

    graph->mincut = mincut;
    graph->nbnd = htable->nelem;

    ASSERT(CheckBndSize(graph));

    if (mincutorder == 0 || mincut >= initcut)
      break;
  }

  freebucket(&parts[0]);
  freebucket(&parts[1]);

  icorefree(2*graph->nvtxs);
}




/*************************************************************************
* This function performs Kernighan-Lin refinement by swapping a vertex 
* at a time. Essentially implements the FM algorithm.
**************************************************************************/
void BFMR_Refine_Weighted(CoarseGraphType *graph, int zeropwgt, int smart, int npasses)
{
  int i, j;
  int k, kwgt;			/* The node and its adjacent edge weight */
  BucketListType parts[2]; 
  int *moved, *swaps, nswaps;
  int higain;			/* High gain node from a partition */
  int order;			/* The order at which nodes got moved */
  int mincut, mincutorder, newcut;	/* The minimum cut so far, the order at which it occured */
  int *id, *ed, *where;
  HTableType *htable;
  VertexType *vtx;
  int oldgain;
  int initcut;
  int pass=0;
  int *pwgts;                 /* The weights of the partitions */
  int from, me, olded, limit, status, halfsplit;
  int maxpwgt[2];

  ASSERT(CheckBndSize(graph));

  maxpwgt[0] = 1.15*zeropwgt;
  maxpwgt[1] = 1.15*(graph->tvwgt - zeropwgt);

  smart = (smart == SMART);

  starttimer(&GreedyInitTmr);

  initbucket(&parts[0], graph->tvwgt, graph->nvtxs, graph->nvtxs, graph->level);
  initbucket(&parts[1], graph->tvwgt, graph->nvtxs, graph->nvtxs, graph->level);

  moved = icoremalloc(graph->nvtxs, "RefineKernLinEqualSingleVertex: moved", 0);
  swaps = icoremalloc(graph->nvtxs, "RefineKernLinEqualSingleVertex: swaps", 0);

  id = graph->id;
  ed = graph->ed;
  where = graph->where;
  htable = &(graph->htable);
  pwgts = graph->pwgts;

  halfsplit = (graph->tvwgt > 2*(zeropwgt-1) && graph->tvwgt < 2*(zeropwgt+1) ? 1 : 0);

  stoptimer(&GreedyInitTmr);

  for (pass=0; pass<npasses; pass++) {
    starttimer(&GreedyInitTmr);

    resetbucket(&parts[0]);
    resetbucket(&parts[1]);

    iset(graph->nvtxs, 0, moved);

    status = 0;
    mincutorder = 0;
    initcut = mincut = newcut = graph->mincut;
    limit = amax(0.01*htable->nelem, 15);
    limit = amin(limit, 50);

    for (j=0; j<htable->size; j++) {
      i = htable->ht[j];
      if (i != HT_EMPTY) {
        if (ed[i] > 0)  { /* Add only boundary vertices into the partition */
          Add2Part(&parts[where[i]], i, ed[i]-id[i]);
          moved[i] = -1;
        }
      }
    }
    stoptimer(&GreedyInitTmr);

    if (__Ctrl->dbglvl&DBG_FFCUT)
      if (pass == 0)
        printf("[W]Partitions: [%6d,%8d]  [%6d,%8d], BND: %5d, Initial Cut: %8d [%d] %d\n",
              parts[0].nnodes, pwgts[0], parts[1].nnodes, pwgts[1], htable->nelem, newcut, ComputeCut(graph), smart);


    /******************************************************
    * Get into the FM loop
    *******************************************************/
    starttimer(&GreedyIterTmr);
    order = nswaps = 0;
    for (;;) {
      order++;

      if (halfsplit)
        from = (pwgts[0] < pwgts[1] ? 1 : 0);
      else
        from = (pwgts[0] < zeropwgt ? 1 : 0);

      higain = GetMaxGainVtx(&parts[from]);

      if (higain == -1)
        break;

      if (pwgts[(from+1)%2] + graph->vtxs[higain]->vwgt > maxpwgt[(from+1)%2]) {
        moved[higain] = 0;
        order--;
        continue;
      }

      newcut -= (ed[higain] - id[higain]);
      if (newcut <= mincut) {
        mincut = newcut;
        mincutorder = order;
      }
      else {
        if (order - mincutorder > limit)
          break; /* No further improvement, break out */
      }

      moved[higain] = order;
      swaps[nswaps++] = higain;

      me = where[higain] = (from+1)%2;
      vtx = graph->vtxs[higain];
      INC_DEC(pwgts[me], pwgts[from], vtx->vwgt);

      if (__Ctrl->dbglvl&DBG_ITERCUT) {
        printf("Gains from moving %8d from %2d: %4d, \tNewcut: %d [%d] %d %d %d\n", higain, from, 
                ed[higain]-id[higain], newcut, ComputeCut(graph), pwgts[0], pwgts[1], vtx->vwgt);
      }

      /**********************************************************
      * Update the id[i]/ed[i] values of the affected nodes
      ***********************************************************/
      if (ed[higain] == 0 && id[higain] > 0) 
        AddHTable(htable, higain);
      else if (ed[higain] > 0 && id[higain] == 0) 
        DelHTable(htable, higain);
      SWAP(id[higain], ed[higain], j);

      for (j=0; j<vtx->nedges; j++) {
        k = vtx->edges[j].edge;
        oldgain = ed[k]-id[k];
        olded = ed[k];

        kwgt = (me == where[k] ? vtx->edges[j].ewgt : -vtx->edges[j].ewgt); 
        INC_DEC(id[k], ed[k], kwgt);

        if (moved[k] == -1 && status == 0)  /* Update only boundary vertices */
          status = UpdatePart(&parts[where[k]], k, oldgain, ed[k]-id[k]);
        else if (smart && status == 0) {
          if (moved[k] == 0 && ed[k]-id[k] > 0) {
            status = Add2Part(&parts[where[k]], k, ed[k]-id[k]);
            moved[k] = -1;
          }
        }

        if (olded == 0 && ed[k] > 0) 
          AddHTable(htable, k);
        else if (olded > 0 && ed[k] == 0) 
          DelHTable(htable, k);
      }
      if (status == -1) {  /* Exit if there where problems in the UpdatePart */
        order++;
        break;
      }
    }
    stoptimer(&GreedyIterTmr);

    if (__Ctrl->dbglvl&DBG_FFCUT)
      printf("\tMinimum Cut: %8d at %5d (%d) [%6d %6d]\n",mincut, mincutorder, order, pwgts[0], pwgts[1]);

    /****************************************************************
    * Roll back computation 
    *****************************************************************/
    starttimer(&GreedyWrapUpTmr);
    for (nswaps--; nswaps>=0; nswaps--) {
      higain = swaps[nswaps];
      if (moved[higain] > mincutorder) {
        me = where[higain] = (where[higain]+1)%2;

        if (ed[higain] == 0 && id[higain] > 0) 
          AddHTable(htable, higain);
        else if (ed[higain] > 0 && id[higain] == 0) 
          DelHTable(htable, higain);
        SWAP(id[higain], ed[higain], j);
  
        vtx = graph->vtxs[higain];
        INC_DEC(pwgts[me], pwgts[(me+1)%2], vtx->vwgt);
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
      else
        break;
    }
    stoptimer(&GreedyWrapUpTmr);

    graph->mincut = mincut;
    graph->nbnd = htable->nelem;

    ASSERT(CheckBndSize(graph));

    if (mincutorder == 0 || mincut >= initcut)
      break;
  }

  freebucket(&parts[0]);
  freebucket(&parts[1]);

  icorefree(2*graph->nvtxs);
}



/*************************************************************************
* This function performs Kernighan-Lin refinement by swapping a vertex 
* at a time. Essentially implements the FM algorithm.
**************************************************************************/
void BFMR_Refine_EqWgt(CoarseGraphType *graph, int zeropwgt)
{
  int i, j;
  int k, kwgt;			/* The node and its adjacent edge weight */
  BucketListType parts[2]; 
  int *moved, *swaps, nswaps;
  int higain;			/* High gain node from a partition */
  int order;			/* The order at which nodes got moved */
  int mindiff, mindifforder, newdiff;	/* The minimum cut so far, the order at which it occured */
  int mincut, newcut;
  int *id, *ed, *where;
  HTableType *htable;
  VertexType *vtx;
  int oldgain;
  int *pwgts;                 /* The weights of the partitions */
  int from, to, olded, limit, status, maxcut;

  if (graph->pwgts[0] == zeropwgt)
    return;

  ASSERT(CheckBndSize(graph));

  TIMELVL(starttimer(&GreedyInitTmr));

  initbucket(&parts[0], graph->tvwgt, graph->nvtxs, graph->nvtxs, graph->level);
  initbucket(&parts[1], graph->tvwgt, graph->nvtxs, graph->nvtxs, graph->level);

  moved = icoremalloc(graph->nvtxs, "RefineKernLinEqualSingleVertex: moved", 0);
  swaps = icoremalloc(graph->nvtxs, "RefineKernLinEqualSingleVertex: swaps", 0);

  id = graph->id;
  ed = graph->ed;
  where = graph->where;
  htable = &(graph->htable);
  pwgts = graph->pwgts;

  iset(graph->nvtxs, 0, moved);

  status = 0;
  mindifforder = 0;
  mindiff = newdiff = abs(pwgts[0] - zeropwgt); 
  newcut = mincut = graph->mincut;
  limit = amax(0.01*htable->nelem, 30);
  limit = amin(limit, 80);

#ifdef MODULE_WEIGHTS
  if (__Ctrl->IsWeighted == 0)
    maxcut = mincut + graph->nedges/(2*graph->nvtxs);
  else
    maxcut = mincut;
#else
  maxcut = mincut;
#endif

  for (j=0; j<htable->size; j++) {
    i = htable->ht[j];
    if (i != HT_EMPTY) {
      if (ed[i] > 0)  { /* Add only boundary vertices into the partition */
        Add2Part(&parts[where[i]], i, ed[i]-id[i]);
        moved[i] = -1;
      }
    }
  }
  TIMELVL(stoptimer(&GreedyInitTmr));

  if (__Ctrl->dbglvl&DBG_FFCUT)
    printf("[D]Partitions: [%6d,%8d]  [%6d,%8d], BND: %5d, ICut: %8d, IDiff: %d\n",
          parts[0].nnodes, pwgts[0], parts[1].nnodes, pwgts[1], htable->nelem, newcut, newdiff);

  /******************************************************
  * Get into the FM loop
  *******************************************************/
  TIMELVL(starttimer(&GreedyIterTmr));
  order = nswaps = 0;
  for (;;) {
    order++;

    from = (pwgts[0] < zeropwgt ? 1 : 0);

    higain = GetMaxGainVtx(&parts[from]);

    if (higain == -1)
      break;

    to = (from+1)%2;
    vtx = graph->vtxs[higain];
    INC_DEC(pwgts[to], pwgts[from], vtx->vwgt);

    newdiff = abs(pwgts[0] - zeropwgt);
    newcut -= (ed[higain] - id[higain]);

    if (newcut <= maxcut && (newdiff < mindiff || (newdiff == mindiff && newcut <= mincut))) {
      mindiff = newdiff;
      mincut = newcut;
      mindifforder = order;
    }
    else {
      if (order - mindifforder > limit) {
        INC_DEC(pwgts[from], pwgts[to], vtx->vwgt);
        break; /* No further improvement, break out */
      }
    }

    moved[higain] = order;
    swaps[nswaps++] = higain;

    where[higain] = to;

    if (__Ctrl->dbglvl&DBG_ITERCUT) {
      printf("Gains from moving %8d from %2d: %4d, \tNewcut: %d \tNewdiff: %d\n", higain, from, 
              ed[higain]-id[higain], newcut, newdiff);
    }

    /**********************************************************
    * Update the id[i]/ed[i] values of the affected nodes
    ***********************************************************/
    if (ed[higain] == 0 && id[higain] > 0) 
      AddHTable(htable, higain);
    else if (ed[higain] > 0 && id[higain] == 0) 
      DelHTable(htable, higain);
    SWAP(id[higain], ed[higain], j);

    for (j=0; j<vtx->nedges; j++) {
      k = vtx->edges[j].edge;
      oldgain = ed[k]-id[k];
      olded = ed[k];

      kwgt = (to == where[k] ? vtx->edges[j].ewgt : -vtx->edges[j].ewgt); 
      INC_DEC(id[k], ed[k], kwgt);

      if (moved[k] == -1 && status == 0)  /* Update only boundary vertices */
        status = UpdatePart(&parts[where[k]], k, oldgain, ed[k]-id[k]);
      else if (status == 0) {
        if (moved[k] == 0 && ed[k]-id[k] > 0) {
          status = Add2Part(&parts[where[k]], k, ed[k]-id[k]);
          moved[k] = -1;
        }
      }

      if (olded == 0 && ed[k] > 0) 
        AddHTable(htable, k);
      else if (olded > 0 && ed[k] == 0) 
        DelHTable(htable, k);
    }
    if (status == -1) {  /* Exit if there where problems in the UpdatePart */
      order++;
      break;
    }
  }
  TIMELVL(stoptimer(&GreedyIterTmr));

  if (__Ctrl->dbglvl&DBG_FFCUT)
    printf("\tMinimum Cut: %8d at %5d, MinDiff: %4d\n",mincut, mindifforder, mindiff);

  /****************************************************************
  * Roll back computation 
  *****************************************************************/
  TIMELVL(starttimer(&GreedyWrapUpTmr));
  for (nswaps--; nswaps>=0; nswaps--) {
    higain = swaps[nswaps];
    if (moved[higain] > mindifforder) {
      to = where[higain] = (where[higain]+1)%2;

      if (ed[higain] == 0 && id[higain] > 0) 
        AddHTable(htable, higain);
      else if (ed[higain] > 0 && id[higain] == 0) 
        DelHTable(htable, higain);
      SWAP(id[higain], ed[higain], j);

      vtx = graph->vtxs[higain];
      INC_DEC(pwgts[to], pwgts[(to+1)%2], vtx->vwgt);
      for (j=0; j<vtx->nedges; j++) {
        k = vtx->edges[j].edge;
        olded = ed[k];

        kwgt = (to == where[k] ? vtx->edges[j].ewgt : -vtx->edges[j].ewgt);
        INC_DEC(id[k], ed[k], kwgt);

        if (olded == 0 && ed[k] > 0) 
          AddHTable(htable, k);
        else if (olded > 0 && ed[k] == 0) 
          DelHTable(htable, k);
      }
    }
    else
      break;
  }
  TIMELVL(stoptimer(&GreedyWrapUpTmr));

  graph->mincut = mincut;
  graph->nbnd = htable->nelem;

  ASSERT(CheckBndSize(graph));

  freebucket(&parts[0]);
  freebucket(&parts[1]);

  icorefree(2*graph->nvtxs);
}



/*************************************************************************
* This function performs Kernighan-Lin refinement by swapping a vertex 
* at a time. Essentially implements the FM algorithm.
**************************************************************************/
void Greedy_Refine(CoarseGraphType *graph, int npasses)
{
  int ii, i, j;
  int k, kwgt;			/* The node and its adjacent edge weight */
  int newcut;			/* The minimum cut so far, the order at which it occured */
  int *id, *ed, *where;
  HTableType *htable;
  VertexType *vtx;
  int pass=0;
  int *pwgts;                 /* The weights of the partitions */
  int from, me, olded, movewgt, halfwgt;

  ASSERT(CheckBndSize(graph));

  TIMELVL(starttimer(&GreedyInitTmr));

  id = graph->id;
  ed = graph->ed;
  where = graph->where;
  htable = &(graph->htable);
  pwgts = graph->pwgts;

  halfwgt = graph->tvwgt/2;
  movewgt = halfwgt*0.90;

  TIMELVL(stoptimer(&GreedyInitTmr));

  if (__Ctrl->dbglvl&DBG_FFCUT)
    printf("Initial Cut: %8d [%6d %6d] [%d]\n", graph->mincut, pwgts[0], pwgts[1], ComputeCut(graph));

  starttimer(&GreedyIterTmr);
  for (pass=0; pass<npasses; pass++) {
    newcut = graph->mincut;

    for (ii=0; ii<htable->size; ii++) {
      if ((i = htable->ht[ii]) == HT_EMPTY)
        continue;

      if ((ed[i] > id[i] && pwgts[where[i]] > movewgt) || (ed[i] == id[i] && pwgts[where[i]] > halfwgt)) {
        vtx = graph->vtxs[i];
        newcut -= (ed[i] - id[i]);
        from = where[i];
        me = where[i] = (from+1)%2;
        INC_DEC(pwgts[me], pwgts[from], vtx->vwgt);

        if (ed[i] == 0 && id[i] > 0) 
          AddHTable(htable, i);
        else if (ed[i] > 0 && id[i] == 0) 
          DelHTable(htable, i);
        SWAP(id[i], ed[i], j);

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

    if (__Ctrl->dbglvl&DBG_FFCUT)
      printf("\t%8d [%6d %6d] [%d]\n", newcut, pwgts[0], pwgts[1], ComputeCut(graph));

    if (newcut >= graph->mincut)
      break;

    graph->mincut = newcut;
  }
  stoptimer(&GreedyIterTmr);

}





/*************************************************************************
* This function prints a partition
**************************************************************************/
void printwhere(CoarseGraphType *graph)
{
  int i, j, *where;

  where = graph->where;

  printf("Partition\n");
  for (i=0, j=0; i<graph->nvtxs; i++) {
    if (where[i] == 0) {
      printf("%8d",i);
      j++;
      if (j%10 == 0)
        printf("\n");
    }
  }

  printf("\n");
  for (i=0, j=0; i<graph->nvtxs; i++) {
    if (where[i] == 1) {
      printf("%8d",i);
      j++;
      if (j%10 == 0)
        printf("\n");
    }
  }

  printf("\n");
}


/*************************************************************************
* This function checks the size of the boundary against graph->nbnd
**************************************************************************/
int CheckBndSize(CoarseGraphType *graph)
{
  int i, j=0;

  for (i=0; i<graph->nvtxs; i++)
    if (graph->ed[i] > 0)
      j++;

  if (j !=graph->nbnd)
    printf("Computed Bnd: %d, Stored Bnd: %d\n", j, graph->nbnd);

  return j==graph->nbnd;
}
