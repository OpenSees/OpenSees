/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * coarsen.c
 *
 * This file contains functions that perform the graph coarsening
 *
 * Started 8/28/94
 * George
 *
 * $Id: coarsen.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include "multilevel.h"

/*************************************************************************
* External Variables
**************************************************************************/
extern CtrlType *__Ctrl;	/* mlevelpart.c.c */
#ifndef METISLIB
extern timer CoarsenTmr;	/* main.c */
#endif


/*************************************************************************
* This function takes a graph and creates the link list of coarser graphs
**************************************************************************/
CoarseGraphType *Coarsen(CoarseGraphType *graph, int CoarsenTo)
{
  CoarseGraphType *cgraph;
  int i;
  int switchlvl, maxvwgt;

  starttimer(&CoarsenTmr);

  cgraph = graph;
  cgraph->level = 0;

  if (graph->nvtxs / CoarsenTo < 4)
    switchlvl = 0;
  else {
    if (__Ctrl->IsWeighted == 0) {
      maxvwgt = graph->tvwgt/CoarsenTo;
      maxvwgt = maxvwgt/8;
      for (switchlvl=0, i=1; ; switchlvl++, i*=2) {
        if (i > maxvwgt)
          break;
      }
      if (switchlvl > 4)
        switchlvl = 4;
    }
    else
      switchlvl = 0;
  }

  i = 0;
  while (cgraph->nvtxs > CoarsenTo) {
    if (EdgePoolSizeLeft() < 0.70*cgraph->nedges) { /* Let's pretend that we run out of edges */
      printf("\n***Warning*** Coarsening aborted prematurily [%d %d]\n",graph->nvtxs, cgraph->nvtxs);
      break;
    }

    if (__Ctrl->dbglvl&DBG_PROGRESS) {
      printf("%7d ",cgraph->nvtxs);  
      fflush(stdout); 
    }
    if (__Ctrl->dbglvl&DBG_CRSSTAT) {
      PrintGraphMMM(cgraph);
      fflush(stdout); 
    }

    if (SelectMatching(cgraph, switchlvl) == 0) {  /* No memory for further coarsening */
      printf("\n***Warning*** Coarsening aborted prematurily [%d %d]\n",graph->nvtxs, cgraph->nvtxs);
      break;
    }

    cgraph->coarser->tvwgt = cgraph->tvwgt;
    cgraph = cgraph->coarser;
    cgraph->level = i+1;

    /* Make sure that we exit if things do not coarsen!  */
    if (cgraph->finer->nvtxs*__Ctrl->cfrac < cgraph->nvtxs) {
      if (__Ctrl->dbglvl&DBG_PROGRESS)
        printf("\nSlow progress in coarsening. Aborted at %d nodes\n", cgraph->nvtxs);
      break;
    }

    i++;

#ifndef METISLIB
    if (__Ctrl->dbglvl&DBG_GRAPHOUT)
      WriteGraph(cgraph, "cgr");
#endif
  }

  if (__Ctrl->dbglvl&DBG_PROGRESS) {
    printf("%7d [%d]",cgraph->nvtxs, EdgePoolSizeLeft());  
    fflush(stdout); 
  }
  if (__Ctrl->dbglvl&DBG_CRSSTAT) {
    PrintGraphMMM(cgraph);
    fflush(stdout); 
  }
  if (__Ctrl->dbglvl&DBG_PROGRESS) 
    printf("\n");

  stoptimer(&CoarsenTmr);

  return cgraph;
}


/*************************************************************************
* This function takes a graph and creates the link list of coarser graphs
**************************************************************************/
CoarseGraphType *KwayCoarsen(CoarseGraphType *graph, int CoarsenTo)
{
  CoarseGraphType *cgraph;
  int i;
  int switchlvl, maxvwgt;

  starttimer(&CoarsenTmr);

  cgraph = graph;
  cgraph->level = 0;

  if (graph->nvtxs / CoarsenTo < 4)
    switchlvl = 0;
  else {
    if (__Ctrl->IsWeighted == 0) {
      maxvwgt = graph->tvwgt/(NPARTS_FACTOR*__Ctrl->nparts);
      maxvwgt = maxvwgt/8;
      for (switchlvl=0, i=1; ; switchlvl++, i*=2) {
        if (i > maxvwgt)
          break;
      }
      if (switchlvl > 4)
        switchlvl = 4;
    }
    else
      switchlvl = 0;
  }

  i = 0;
  while (cgraph->nvtxs > CoarsenTo) {
    if (EdgePoolSizeLeft() < cgraph->nedges) { /* Let's pretend that we run out of edges */
      printf("\n***Warning*** Coarsening aborted prematurily [%d %d]\n",graph->nvtxs, cgraph->finer->nvtxs);
      cgraph = cgraph->finer;
      FreeGraph(cgraph->coarser);
      cgraph->coarser = NULL;
      GKfree(cgraph->cmap, cgraph->match, -1);
      cgraph->cmap = cgraph->match = NULL;
      break;
    }

    if (__Ctrl->dbglvl&DBG_PROGRESS) {
      printf("%7d ",cgraph->nvtxs);  
      fflush(stdout); 
    }
    if (__Ctrl->dbglvl&DBG_CRSSTAT) {
      PrintGraphMMM(cgraph);
      fflush(stdout); 
    }

    KwaySelectMatching(cgraph, switchlvl);

    cgraph->coarser->tvwgt = cgraph->tvwgt;
    cgraph = cgraph->coarser;
    cgraph->level = i+1;

    /* Make sure that we exit if things do not coarsen!  */
    if (cgraph->finer->nvtxs*__Ctrl->cfrac < cgraph->nvtxs) {
      if (__Ctrl->dbglvl&DBG_PROGRESS)
        printf("\nSlow progress in coarsening. Aborted at %d nodes\n", cgraph->nvtxs);
      break;
    }

    i++;

#ifndef METISLIB
    if (__Ctrl->dbglvl&DBG_GRAPHOUT)
      WriteGraph(cgraph, "cgr");
#endif
  }

  if (__Ctrl->dbglvl&DBG_PROGRESS) {
    printf("%7d [%d]",cgraph->nvtxs, EdgePoolSizeLeft());  
    fflush(stdout); 
  }
  if (__Ctrl->dbglvl&DBG_CRSSTAT) {
    PrintGraphMMM(cgraph);
    fflush(stdout); 
  }
  if (__Ctrl->dbglvl&DBG_PROGRESS) 
    printf("\n");

  stoptimer(&CoarsenTmr);

  return cgraph;
}




/*************************************************************************
* This function selects the appropriate function depending on MatchType
**************************************************************************/
int SelectMatching(CoarseGraphType *cgraph, int switchlvl)
{
  switch (__Ctrl->MatchType) {
    case MATCH_RM:
      if (__Ctrl->IsWeighted == 0 || cgraph->level < switchlvl)
        return RM_Match(cgraph);
      else
        return RM_Match_W(cgraph);
      break;
    case MATCH_HEM:
      if (__Ctrl->IsWeighted == 0 && cgraph->level == 0) {
        return RM_Match(cgraph);
      }
      else {
        if (__Ctrl->IsWeighted == 0 || cgraph->level < switchlvl)
          return HEM_Match(cgraph);
        else
          return HEM_Match_W(cgraph);
      }
      break;
    case MATCH_LEM:
      if (__Ctrl->IsWeighted == 0 && cgraph->level == 0) {
        return RM_Match(cgraph);
      }
      else {
        if (__Ctrl->IsWeighted == 0 || cgraph->level < switchlvl)
          return LEM_Match(cgraph);
        else
          return LEM_Match_W(cgraph);
      }
      break;
    case MATCH_HCM:
      if (__Ctrl->IsWeighted == 0 && cgraph->level == 0) {
        return RM_Match(cgraph);
      }
      else {
        if (__Ctrl->IsWeighted == 0 || cgraph->level < switchlvl)
          return HCM_Match(cgraph);
        else
          return HCM_Match_W(cgraph);
      }
      break;
    case MATCH_MHEM:
      if (cgraph->level < switchlvl)
        return MHEM_Match(cgraph);
      else
        return MHEM_Match_W(cgraph);
      break;
    case MATCH_SRM:
      if (cgraph->level < switchlvl)
        return RM_Match(cgraph);
      else
        return SRM_Match(cgraph);
      break;
    case MATCH_SHEM:
      if (__Ctrl->IsWeighted == 0 && cgraph->level == 0) {
        return RM_Match(cgraph);
      }
      else {
        if (cgraph->level < switchlvl)
          return HEM_Match(cgraph);
        else
          return SHEM_Match(cgraph);
      }
      break;
    case MATCH_SMHEM:
      if (cgraph->level < switchlvl)
        return MHEM_Match(cgraph);
      else
        return SMHEM_Match(cgraph);
      break;
    default:
      errexit("Unsupported Match Type: %d", __Ctrl->MatchType);
  }

  return 0;
}



/*************************************************************************
* This function selects the appropriate function depending on MatchType
**************************************************************************/
void KwaySelectMatching(CoarseGraphType *cgraph, int switchlvl)
{
  switch (__Ctrl->MatchType) {
    case MATCH_RM:
      if (cgraph->level < switchlvl)
        RM_Match(cgraph);
      else
        RM_Match_W(cgraph);
      break;
    case MATCH_HEM:
      if (__Ctrl->IsWeighted == 0 && cgraph->level == 0) {
        RM_Match(cgraph);
      }
      else {
        if (cgraph->level < switchlvl)
          HEM_Match(cgraph);
        else
          HEM_Match_W(cgraph);
      }
      break;
    case MATCH_LEM:
      if (__Ctrl->IsWeighted == 0 && cgraph->level == 0) {
        RM_Match(cgraph);
      }
      else {
        if (cgraph->level < switchlvl)
          LEM_Match(cgraph);
        else
          LEM_Match_W(cgraph);
      }
      break;
    case MATCH_HCM:
      if (__Ctrl->IsWeighted == 0 && cgraph->level == 0) {
        RM_Match(cgraph);
      }
      else {
        if (cgraph->level < switchlvl)
          HCM_Match(cgraph);
        else
          HCM_Match_W(cgraph);
      }
      break;
    case MATCH_MHEM:
      if (cgraph->level < switchlvl)
        MHEM_Match(cgraph);
      else
        MHEM_Match_W(cgraph);
      break;
    case MATCH_SRM:
      if (cgraph->level < switchlvl)
        RM_Match(cgraph);
      else
        SRM_Match(cgraph);
      break;
    case MATCH_SHEM:
      if (__Ctrl->IsWeighted == 0 && cgraph->level == 0) {
        RM_Match(cgraph);
      }
      else {
        if (cgraph->level < switchlvl)
          HEM_Match(cgraph);
        else
          SHEM_Match(cgraph);
      }
      break;
    case MATCH_SMHEM:
      if (cgraph->level < switchlvl)
        MHEM_Match(cgraph);
      else
        SMHEM_Match(cgraph);
      break;
    default:
      errexit("Unsupported Match Type: %d", __Ctrl->MatchType);
  }

}

