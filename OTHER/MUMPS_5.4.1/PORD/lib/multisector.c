/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: ms.c
/
/ author        J"urgen Schulze, University of Paderborn
/ created       01jan04
/
/ This file contains functions dealing with the multisector object
/
******************************************************************************

Data type:  struct multisector
              graph_t *G;         pointer to original graph
              int     *stage;     stage[u]=i => node u will be elim. in stage i
              int     nstages;    number of stages
              int     nnodes;     number of nodes in multisector
              int     totmswght;  weigth of nodes in multisector
Comments:
  o Structure does not own graph object G => it will not be freed
    Note: G is the original graph
Methods in lib/multisector.c:
- ms = newMultisector(graph_t *G);
    o Initial: nstages = nnodes = totmswght = 0;
- void freeMultisector(ms_t *ms); 
- ms = trivialMultisector(graph_t *G);
    o allocates memory for the multisector object by a call to newMultisector
      and sets stage[u] = 0 for all vertices u and nstages = 1; the trivial
      multisector can be used for pure bottom-up orderings
- ms = constructMultisector(graph_t *G, options_t* options, timings_t *cpus);
    o MASTER_FUNCTION: computes a multisector for G according to the specified
      ordtype e { MINIMUM_PRIORITY, INCOMPLETE_ND, MULTISECTION,
                  TRISTAGE_MULTISECTION }
      MINIMUM_PRIORTY:
         return the multisector obtained by a call to trivialMultisector
      INCOMPLETE_ND, MULTISECTION, TRISTAGE_MULTISECTION:
         build separator tree by calling buildNDtree and extract multisector
         by calling extractMS2stage (MULTISECTION) or extractMSmultistage
         (INCOMPLETE_ND, TRISTAGE_MULTISECTION)
    o used options: (also see buildNDtree)
       OPTION_ORDTYPE, OPTION_DOMAIN_SIZE, OPTION_MSGLVL, OPTION_NODE_SELECTION3
    o returned timings: (see buildNDtree)
       TIME_INITDOMDEC, TIME_COARSEDOMDEC, TIME_INITSEP, TIME_REFINESEP
       TIME_MULTILEVEL, TIME_SMOOTH
- ms = extractMS2stage(nestdiss_t *ndroot);
    o extracts a 2-stage multisector from the nested dissection tree with root
      ndroot: stage[u] = 0 => u belongs to a domain
              stage[u] = 1 => u belongs to the multisector 
      and nstages = 2; the 2-stage multisector can be used for classical
      multisection orderings
- ms = extractMSmultistage(nestdiss_t *ndroot);
    o extracts a multi-stage multisector from the nested dissection tree at
      ndroot: stage[u] = 0 => u belongs to a domain
              stage[u] = i, i > 0 => u belongs to the multisector, i.e.:
                   stage[u] = 1 => u belongs to a leaf separator
                      :
                   stage[u] = nstages-1 => u belongs to the root separator
      the multisector can be used for incomplete nested dissection orderings
      or for three-stage multisection orderings

******************************************************************************/

#include <space.h>


/*****************************************************************************
******************************************************************************/
multisector_t*
newMultisector(graph_t *G)
{ multisector_t *ms;

  mymalloc(ms, 1, multisector_t);
  mymalloc(ms->stage, G->nvtx, PORD_INT);

  ms->G = G;
  ms->nstages = 0;
  ms->nnodes = 0;
  ms->totmswght = 0;

  return(ms);
}


/*****************************************************************************
******************************************************************************/
void
freeMultisector(multisector_t *ms)
{
  free(ms->stage);
  free(ms);
}


/*****************************************************************************
******************************************************************************/
multisector_t*
trivialMultisector(graph_t *G)
{ multisector_t *ms;
  PORD_INT           *stage, nvtx, u;
  
  /* ----------------------------------------------------------------- 
     allocate memory for the multisector object and init. stage vector
     ----------------------------------------------------------------- */
  nvtx = G->nvtx;
  ms = newMultisector(G);
  stage = ms->stage;

  for (u = 0; u < nvtx; u++)
    stage[u] = 0;                      /* no vertex belongs to a separator */

  /* -------------------------------
     finalize the multisector object
     ------------------------------- */
  ms->nstages = 1;
  ms->nnodes = 0;
  ms->totmswght = 0;

  return(ms);
}


/*****************************************************************************
******************************************************************************/
multisector_t*
constructMultisector(graph_t *G, options_t* options, timings_t *cpus)
{ multisector_t *ms;
  nestdiss_t    *ndroot;
  PORD_INT           *map, nvtx, ordtype;

  nvtx = G->nvtx;

  /* ------------------------------
     check number of nodes in graph
     ------------------------------ */
  /* -----------------------------------
     JY: inserted the condition
     "&& (options[OPTION_MSGLVL] > 0)"
     below, to avoid systematic printing
     ----------------------------------- */
    if ((nvtx <= MIN_NODES) && (options[OPTION_ORDTYPE] != MINIMUM_PRIORITY)
        && (options[OPTION_MSGLVL] > 0))
   { printf("\nWarning in constructMultisector\n"
         "  graph has less than %d nodes, skipping separator construction\n\n",
         MIN_NODES);
     options[OPTION_ORDTYPE] = MINIMUM_PRIORITY;
   }
  /* --------------------------------------------------------
     determine the multisector according to the ordering type
     -------------------------------------------------------- */
  ordtype = options[OPTION_ORDTYPE];
   switch(ordtype)
   { case MINIMUM_PRIORITY:
       ms = trivialMultisector(G);
       break;

     case INCOMPLETE_ND:
     case MULTISECTION:
     case TRISTAGE_MULTISECTION:
       mymalloc(map, nvtx, PORD_INT);
       ndroot = setupNDroot(G, map);
       buildNDtree(ndroot, options, cpus);
       if (ordtype == MULTISECTION)
         ms = extractMS2stage(ndroot);
       else
         ms = extractMSmultistage(ndroot);
       freeNDtree(ndroot);
       freeNDnode(ndroot);
       free(map);
       break;

     default:
       fprintf(stderr, "\nError in function constructMultisector\n"
            "  unrecognized ordering type %d\n", ordtype);
       quit();
   }
  return(ms);
}

  
/*****************************************************************************
******************************************************************************/
multisector_t*
extractMS2stage(nestdiss_t *ndroot)
{ multisector_t *ms;
  nestdiss_t    *nd, *parent;
  PORD_INT           *stage, *intvertex, *intcolor;
  PORD_INT           nvint, nnodes, totmswght, i;

  /* ----------------------------------------------------------------- 
     allocate memory for the multisector object and init. stage vector
     ----------------------------------------------------------------- */
  ms = trivialMultisector(ndroot->G);
  stage = ms->stage;

  /* ------------------------------------------------------------
     extract the stages of the separator vertices:
     stage[u] = 1, iff u belongs to a separator
     ------------------------------------------------------------ */
  nnodes = totmswght = 0;
  for (nd = ndroot; nd->childB != NULL; nd = nd->childB);
  while (nd != ndroot)
   { parent = nd->parent;
     if ((parent == NULL) || (parent->childB == NULL)
        || (parent->childW == NULL))
       { fprintf(stderr, "\nError in function extractMS2stage\n"
              "  nested dissection tree corrupted\n");
         quit();
       }
      if (parent->childB == nd)        /* left subtree of parent visited */
        for (nd = parent->childW; nd->childB != NULL; nd = nd->childB);
      else                             /* right subtree of parent visited */
       { nd = parent;                  /* extract the separator of parent */
         totmswght += nd->cwght[GRAY];
         nvint = nd->nvint;
         intvertex = nd->intvertex;
         intcolor = nd->intcolor;
         for (i = 0; i < nvint; i++)
           if (intcolor[i] == GRAY)
            { nnodes++;
              stage[intvertex[i]] = 1;
            }
       }
   }

  /* ------------------------------------------
     finalize the multisector object and return
     ------------------------------------------ */
  ms->nstages = 2;
  ms->nnodes = nnodes;
  ms->totmswght = totmswght;

  return(ms);
}

  
/*****************************************************************************
******************************************************************************/
multisector_t*
extractMSmultistage(nestdiss_t *ndroot)
{ multisector_t *ms;
  nestdiss_t    *nd, *parent;
  PORD_INT           *stage, *intvertex, *intcolor;
  PORD_INT           nvtx, nvint, maxstage, istage, nnodes, totmswght, i, u;

  /* -----------------------------------------------------------------
     allocate memory for the multisector object and init. stage vector
     ----------------------------------------------------------------- */
  ms = trivialMultisector(ndroot->G);
  stage = ms->stage;

  /* ------------------------------------------------------------
     extract the stages of the separator vertices:
     stage[u] = i, i>0, iff u belongs to a separator in depth i-1
     ------------------------------------------------------------ */
  maxstage = nnodes = totmswght = 0;
  for (nd = ndroot; nd->childB != NULL; nd = nd->childB);
  while (nd != ndroot)
   { parent = nd->parent;
     if ((parent == NULL) || (parent->childB == NULL)
        || (parent->childW == NULL))
       { fprintf(stderr, "\nError in function extractMSmultistage\n"
              "  nested dissection tree corrupted\n");
         quit();
       }
      if (parent->childB == nd)        /* left subtree of parent visited */
        for (nd = parent->childW; nd->childB != NULL; nd = nd->childB);
      else                             /* right subtree of parent visited */
       { nd = parent;                  /* extract the separator of parent */
         istage = nd->depth + 1;       /* sep. vertices belong to this stage */
         maxstage = max(maxstage, istage);
         totmswght += nd->cwght[GRAY];
         nvint = nd->nvint;
         intvertex = nd->intvertex;
         intcolor = nd->intcolor;
         for (i = 0; i < nvint; i++)
           if (intcolor[i] == GRAY)
            { nnodes++;
              stage[intvertex[i]] = istage;
            }
       }
   }

  /* --------------------------------------------------------------------
     we have: stage[u] = 0 => u belongs to a domain
              stage[u] = 1 => u belongs to the root separator (depth = 0)
                 :
              stage[u] = maxstage => u belongs to a leaf separator
     but we must eliminate the separators in a bottom-up fashion; we like
     to have: stage[u] = 0 => u belongs to a domain
              stage[u] = 1 => u belongs to a leaf separator
                 :
              stage[u] = maxstage => u belongs to the root separator
     -------------------------------------------------------------------- */
  nvtx = ndroot->G->nvtx;
  for (u = 0; u < nvtx; u++)
    if (stage[u] > 0)
      stage[u] = maxstage - stage[u] + 1;

  /* ------------------------------------------
     finalize the multisector object and return
     ------------------------------------------ */
  ms->nstages = maxstage + 1;
  ms->nnodes = nnodes;
  ms->totmswght = totmswght;

  return(ms);
}
