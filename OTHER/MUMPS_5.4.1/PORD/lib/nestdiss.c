/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: nestdiss.c
/
/ author        J"urgen Schulze, University of Paderborn
/ created       00dec29
/
/ This file contains functions dealing with the rec. nested dissection object
/
******************************************************************************

Data type:  struct nestdiss
              graph_t *G;               pointer to original graph
              int     *map;             maps nodes of G to constructed subgraph
              int     depth;            depth in nested dissection tree
              int     nvint;            number of vertices in subgraph
              int     *intvertex;       internal vertices of subgraph
              int     *intcolor;        color of vertices in intvertex
              int     cwght[3];         weights of bisection
              struct nestdiss *parent;  pointer to parent nd node
              struct nestdiss *childB;  pointer to black descendant nd node
              struct nestdiss *childW;  pointer to white descendand nd node
Comments:
  o Structure used to build the nested dissection tree. Vector intvertex
    holds the vertices of the subgraph to be partitioned. Once a separator
    has been computed, the coloring of vertex u = intvertex[i] is stored in
    vector intcolor[i] and the partition weights are stored in cwght[GRAY],
    cwght[BLACK], and cwght[WHITE].
  o Structure does not own graph object G => it will not be freed
    Note: G is the original graph
  o Structure does not own map array => it will not be freed
    Note: map is a "global" array that is used when constructing the subgraph
          induced by the vertices in intvertex. The array maps the vertices
          of the original graph G to the vertices of the subgraph.
Methods in lib/nestdiss.c:
- nd = newNDnode(graph_t *G, int *map, int nvint);
    o Initial: depth = 0, cwght[GRAY] = cwght[BLACK] = cwght[WHITE] = 0,
               and parent = childB = childW = NULL;
- void freeNDnode(nestdiss_t *nd);
- ndroot = setupNDroot(graph_t *G, int *map);
    o sets up the root of the nested dissection tree; the function first
      calls newNDnode to allocate memory for ndroot and, then, sets
      intvertex[i] = i for all 0 <= i < G->nvtx
- void splitNDnode(nestdiss_t *nd, options_t *options, timings_t *cpus);
    o constructs the subgraph induced by nd->intvertex and computes a
      bisection for it by calling constructSeparator and smoothSeparator. 
      Then, the nd object is splitted in a black one that holds the black 
      partition and a white one that holds the white partition.
    o used options: (see constructSeparator and smoothSeparator)
       OPTION_MSGLVL, OPTION_NODE_SELECTION3
    o returned timings: (also see constructSeparator)
       TIME_INITDOMDEC, TIME_COARSEDOMDEC, TIME_INITSEP, TIME_REFINESEP
       TIME_MULTILEVEL, TIME_SMOOTH
- void buildNDtree(nestdiss_t *ndroot, options_t *options, timings_t *cpus);
    o builds the nested dissection tree under root ndroot, i.e. it applies
      the nested dissection process to the (sub)graph induced by
      ndroot->intvertex by iteratively calling function splitNDnode.
    o used options: (also see splitNDnode)
       OPTION_DOMAIN_SIZE, OPTION_MSGLVL, OPTION_NODE_SELECTION3
    o returned timings: (see splitNDnode)
       TIME_INITDOMDEC, TIME_COARSEDOMDEC, TIME_INITSEP, TIME_REFINESEP
       TIME_MULTILEVEL, TIME_SMOOTH
- void freeNDtree(nestdiss_t *ndroot);
    o removes the nested dissection tree under root ndroot
      Note: ndroot is not freed

******************************************************************************/

#include <space.h>


/*****************************************************************************
******************************************************************************/
nestdiss_t*
newNDnode(graph_t *G, PORD_INT *map, PORD_INT nvint)
{ nestdiss_t *nd;

  mymalloc(nd, 1, nestdiss_t);
  mymalloc(nd->intvertex, nvint, PORD_INT);
  mymalloc(nd->intcolor, nvint, PORD_INT);

  nd->G = G;
  nd->map = map;
  nd->depth = 0;
  nd->nvint = nvint;
  nd->cwght[GRAY] = nd->cwght[BLACK] = nd->cwght[WHITE] = 0;
  nd->parent = nd->childB = nd->childW = NULL;

  return(nd);
}


/*****************************************************************************
******************************************************************************/
void
freeNDnode(nestdiss_t *nd)
{
  free(nd->intvertex);
  free(nd->intcolor);
  free(nd);
}


/*****************************************************************************
******************************************************************************/
nestdiss_t*
setupNDroot(graph_t *G, PORD_INT *map)
{ nestdiss_t *ndroot;
  PORD_INT        *intvertex, nvtx, i;

  nvtx = G->nvtx;
  ndroot = newNDnode(G, map, nvtx);
  intvertex = ndroot->intvertex;

  for (i = 0; i < nvtx; i++)
    intvertex[i] = i;

  return(ndroot);
}


/*****************************************************************************
******************************************************************************/
void
splitNDnode(nestdiss_t *nd, options_t *options, timings_t *cpus)
{ nestdiss_t *b_nd, *w_nd;
  graph_t    *Gsub;
  gbisect_t  *Gbisect;
  PORD_INT        *map, *intvertex, *intcolor, *b_intvertex, *w_intvertex;
  PORD_INT        nvint, b_nvint, w_nvint, u, i;

  map = nd->map;
  nvint = nd->nvint;
  intvertex = nd->intvertex;
  intcolor = nd->intcolor;

  /* -------------------------------------------------------------
     extract the subgraph for which a bisection has to be computed
     ------------------------------------------------------------- */
  if (nd->G->nvtx == nd->nvint)
   { Gsub = nd->G;                    /* a hack to save time and space */
     for (u = 0; u < nd->nvint; u++)  /* but do not forget the map vector */
       map[u] = u;
   }
  else
    Gsub = setupSubgraph(nd->G, intvertex, nvint, map);
  Gbisect = newGbisect(Gsub);

  /* ---------------------------------
     compute the bisection for Gbisect
     --------------------------------- */
  pord_starttimer(cpus[TIME_MULTILEVEL]);
  constructSeparator(Gbisect, options, cpus);
  pord_stoptimer(cpus[TIME_MULTILEVEL]);

  pord_starttimer(cpus[TIME_SMOOTH]);
  if (Gbisect->cwght[GRAY] > 0)
    smoothSeparator(Gbisect, options);
  pord_stoptimer(cpus[TIME_SMOOTH]);

  /* ----------------------------------------
     copy the bisection back to the nd object
     ---------------------------------------- */
  b_nvint = w_nvint = 0;
  nd->cwght[GRAY] = Gbisect->cwght[GRAY];
  nd->cwght[BLACK] = Gbisect->cwght[BLACK];
  nd->cwght[WHITE] = Gbisect->cwght[WHITE];
  for (i = 0; i < nvint; i++)
   { u = intvertex[i];
     intcolor[i] = Gbisect->color[map[u]];
     switch(intcolor[i])
      { case GRAY: break;
        case BLACK: b_nvint++; break;
        case WHITE: w_nvint++; break;
        default:
          fprintf(stderr, "\nError in function splitNDnode\n"
               "  node %d has unrecognized color %d\n", u, intcolor[i]);
          quit();
      }
   }

  /* ------------------------------------------------------
     and now split the nd object according to the bisection
     ------------------------------------------------------ */
  b_nd = newNDnode(nd->G, map, b_nvint);
  b_intvertex = b_nd->intvertex;
  w_nd = newNDnode(nd->G, map, w_nvint);
  w_intvertex = w_nd->intvertex;

  b_nvint = w_nvint = 0;
  for (i = 0; i < nvint; i++)
   { u = intvertex[i];
     if (intcolor[i] == BLACK) b_intvertex[b_nvint++] = u;
     if (intcolor[i] == WHITE) w_intvertex[w_nvint++] = u;
   }
  nd->childB = b_nd; b_nd->parent = nd;
  nd->childW = w_nd; w_nd->parent = nd;
  b_nd->depth = nd->depth + 1;
  w_nd->depth = nd->depth + 1;

  /* -----------------
     free the subgraph
     ----------------- */
  if (Gsub != nd->G)
    freeGraph(Gsub);
  freeGbisect(Gbisect);
}


/*****************************************************************************
******************************************************************************/
void
buildNDtree(nestdiss_t *ndroot, options_t *options, timings_t *cpus)
{ nestdiss_t *nd;
  nestdiss_t *queue[2*MAX_SEPS+1];
  PORD_INT        maxseps, seps, domainsize, qhead, qtail;

  maxseps = MAX_SEPS;
  domainsize = options[OPTION_DOMAIN_SIZE];
  if (domainsize == 1) maxseps = DEFAULT_SEPS;   /* secret switch */

  /* --------------------------------------------------
     build the nested dissection tree under root ndroot
     -------------------------------------------------- */
  queue[0] = ndroot;
  qhead = 0; qtail = 1; seps = 0;
  while ((qhead != qtail) && (seps < maxseps))
   { seps++;
     nd = queue[qhead++];

     splitNDnode(nd, options, cpus);
     if ((nd->childB == NULL) || (nd->childW == NULL))
      { fprintf(stderr, "\nError in function buildNDtree\n"
             "  recursive nested dissection process failed\n");
        quit();
      }

     if (options[OPTION_MSGLVL] > 1)
       printf("%4d. S %6d, B %6d, W %6d [bal %4.2f, rel %6.4f, cost %7.2f]\n",
             seps, nd->cwght[GRAY], nd->cwght[BLACK], nd->cwght[WHITE],
             (FLOAT)min(nd->cwght[BLACK], nd->cwght[WHITE])
                     / max(nd->cwght[BLACK], nd->cwght[WHITE]),
             (FLOAT)nd->cwght[GRAY]
                     / (nd->cwght[GRAY] + nd->cwght[BLACK] + nd->cwght[WHITE]),
             F(nd->cwght[GRAY], nd->cwght[BLACK], nd->cwght[WHITE]));

     if ((nd->childB->nvint > MIN_NODES)
        && ((nd->cwght[BLACK] > domainsize) || (qtail < DEFAULT_SEPS)))
       queue[qtail++] = nd->childB;
     if ((nd->childW->nvint > MIN_NODES)
        && ((nd->cwght[WHITE] > domainsize) || (qtail < DEFAULT_SEPS)))
       queue[qtail++] = nd->childW;
   }
}


/*****************************************************************************
******************************************************************************/
void
freeNDtree(nestdiss_t *ndroot)
{ nestdiss_t *nd, *parent;

  /* ------------------------------------------------------
     to remove the nested dissection tree under root ndroot
     visit the nodes in post-order
     ------------------------------------------------------ */
  for (nd = ndroot; nd->childB != NULL; nd = nd->childB);
  while (nd != ndroot)
   { parent = nd->parent;
     if ((parent == NULL) || (parent->childB == NULL)
        || (parent->childW == NULL))
       { fprintf(stderr, "\nError in function removeNDtree\n"
              "  nested dissection tree corrupted\n");
         quit();
       }
     if (parent->childB == nd)  /* left subtree of parent visited */
      { freeNDnode(nd);         /* free root of left subtree and goto right */
        for (nd = parent->childW; nd->childB != NULL; nd = nd->childB);
      }
     else                       /* right subtree of parent visited */
      { freeNDnode(nd);         /* free root of right subtree and goto parent */
        nd = parent;
      }
   }
}
