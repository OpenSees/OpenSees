/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: gbisect.c
/
/ author        J"urgen Schulze, University of Paderborn
/ created       00dec29
/
/ This file contains functions dealing with the graph bisection object
/
******************************************************************************

Data type:  struct gbisect
              graph_t *G;         pointer to graph that will be partitioned
              int     *color;     color of node (GRAY, BLACK, or WHITE)
              int     cwght[3];   weights of GRAY, BLACK, WHITE partitions
Comments:
  o Structure used to compute the bisection of a graph. Structure does not
    own graph object => it will not be freed.
Methods in lib/gbisect.c:
- Gbisect = newGbisect(graph_t *G);
    o Initial: cwght[GRAY] = cwght[BLACK] = cwght[WHITE] = 0
- void freeGbisect(gbisect_t *Gbisect);
- void printGbisect(gbisect_t *Gbisect);
- void checkSeparator(gbisect_t *Gbisect);
- void constructSeparator(gbisect_t *Gbisect, options_t *options,
                          timings_t *cpus);
    o constructs a vertex separator by applying the new multilevel approach;
      it first constructs an initial domain decomposition for Gbisect->G
      by calling constructDomainDecomposition; the dd is then coarsed by
      several calls to shrinkDomainDecomposition; the last dd is colored
      by a call to initialDDSep; this coloring is refined during the
      uncoarsening phase by several calls to improveDDSep
    o used options:
       OPTION_MSGLVL, OPTION_NODE_SELECTION3
      returned timings: 
       TIME_INITDOMDEC, TIME_COARSEDOMDEC, TIME_INITSEP, TIME_REFINESEP
- int smoothBy2Layers(gbisect_t *Gbisect, int *bipartvertex, int *pnX,
                      int black, int white);
    o on start, bipartvertex contains the nodes of the separator; the
      separator is then paired with eiter the black or the white partition
      so that the nodes in bipartvertex induce a bipartite graph; this
      graph is constructed by setupBipartiteGraph; a Dulmage-Mendelsohn
      decomposition is computed and the separator is smoothed; the
      vertices of the smoothed separator are returned in bipartvertex
- void smoothSeparator(gbisect_t *Gbisect, options_t *options);
    o smoothes a given separator by repeatedly calling smoothBy2Layers
    o used options: OPTION_MSGLVL

******************************************************************************/

#include <space.h>
/* #define DEBUG */
/* #define BE_CAUTIOUS */


/*****************************************************************************
******************************************************************************/
gbisect_t*
newGbisect(graph_t *G)
{ gbisect_t *Gbisect;

  mymalloc(Gbisect, 1, gbisect_t);
  mymalloc(Gbisect->color, G->nvtx, PORD_INT);

  Gbisect->G = G;
  Gbisect->cwght[GRAY] = 0;
  Gbisect->cwght[BLACK] = 0;
  Gbisect->cwght[WHITE] = 0;

  return(Gbisect);
}


/*****************************************************************************
******************************************************************************/
void
freeGbisect(gbisect_t *Gbisect)
{
  free(Gbisect->color);
  free(Gbisect);
}


/*****************************************************************************
******************************************************************************/
void
printGbisect(gbisect_t *Gbisect)
{ graph_t *G;
  PORD_INT     count, u, v, i, istart, istop;

  G = Gbisect->G;
  printf("\n#nodes %d, #edges %d, totvwght %d\n", G->nvtx, G->nedges >> 1,
         G->totvwght);
  printf("partition weights: S %d, B %d, W %d\n", Gbisect->cwght[GRAY],
         Gbisect->cwght[BLACK], Gbisect->cwght[WHITE]);
  for (u = 0; u < G->nvtx; u++)
   { count = 0;
     printf("--- adjacency list of node %d (weight %d, color %d)\n", u,
            G->vwght[u], Gbisect->color[u]);
     istart = G->xadj[u];
     istop = G->xadj[u+1];
     for (i = istart; i < istop; i++)
      { v = G->adjncy[i];
        printf("%5d (color %2d)", v, Gbisect->color[v]);
        if ((++count % 4) == 0)
          printf("\n");
      }
     if ((count % 4) != 0)
       printf("\n");
   }
}


/*****************************************************************************
******************************************************************************/
void
checkSeparator(gbisect_t *Gbisect)
{ PORD_INT *xadj, *adjncy, *vwght, *color, *cwght;
  PORD_INT nvtx, err, checkS, checkB, checkW, a, b, u, v, i, istart, istop;

  nvtx = Gbisect->G->nvtx;
  xadj = Gbisect->G->xadj;
  adjncy = Gbisect->G->adjncy;
  vwght = Gbisect->G->vwght;
  color = Gbisect->color;
  cwght = Gbisect->cwght;

  err = FALSE;
  printf("checking separator of induced subgraph (S %d, B %d, W %d)\n",
         cwght[GRAY], cwght[BLACK], cwght[WHITE]);

  checkS = checkB = checkW = 0;
  for (u = 0; u < nvtx; u++)
   { istart = xadj[u];
     istop = xadj[u+1];
     switch(color[u])
      { case GRAY:                /* is it a minimal separator? */
          checkS += vwght[u];
          a = b = FALSE;
          for (i = istart; i < istop; i++)
           { v = adjncy[i];
             if (color[v] == WHITE) a = TRUE;
             if (color[v] == BLACK) b = TRUE;
           }
          if (!((a) && (b)))
            printf("WARNING: not a minimal separator (node %d)\n", u);
          break;
        case BLACK:               /* is it realy a separator? */
          checkB += vwght[u];
          for (i = istart; i < istop; i++)
           { v = adjncy[i];
             if (color[v] == WHITE)
              { printf("ERROR: white node %d adjacent to black node %d\n", u,v);
                err = TRUE;
              }
           }
          break;
        case WHITE:
          checkW += vwght[u];
          break;
        default:
          printf("ERROR: node %d has unrecognized color %d\n", u, color[u]);
          err = TRUE;
      }
   }

  /* check cwght[GRAY], cwght[BLACK], cwght[WHITE] */
  if ((checkS != cwght[GRAY]) || (checkB != cwght[BLACK])
     || (checkW != cwght[WHITE]))
   { printf("ERROR in partitioning: checkS %d (S %d), checkB %d (B %d), "
            "checkW %d (W %d)\n", checkS, cwght[GRAY], checkB, cwght[BLACK],
             checkW, cwght[WHITE]);
     err = TRUE;
   }
  if (err) quit();
}


/*****************************************************************************
******************************************************************************/
void
constructSeparator(gbisect_t *Gbisect, options_t *options, timings_t *cpus)
{ domdec_t *dd, *dd2;
  PORD_INT      *color, *cwght, *map, nvtx, u, i;

  nvtx = Gbisect->G->nvtx;
  color = Gbisect->color;
  cwght = Gbisect->cwght;

  /* --------------------------------------------------------------
     map vector identifies vertices of Gbisect->G in domain decomp.
     -------------------------------------------------------------- */
  mymalloc(map, nvtx, PORD_INT);

  /* --------------------------------------
     construct initial domain decomposition
     -------------------------------------- */
  pord_starttimer(cpus[TIME_INITDOMDEC]);
  dd = constructDomainDecomposition(Gbisect->G, map);

#ifdef BE_CAUTIOUS
  checkDomainDecomposition(dd);
#endif

  if (options[OPTION_MSGLVL] > 2)
    printf("\t  0. dom.dec.: #nodes %d (#domains %d, weight %d), #edges %d\n",
           dd->G->nvtx, dd->ndom, dd->domwght, dd->G->nedges >> 1);
  pord_stoptimer(cpus[TIME_INITDOMDEC]);

  /* ---------------------------------------------------
     construct sequence of coarser domain decompositions
     --------------------------------------------------- */
  pord_starttimer(cpus[TIME_COARSEDOMDEC]);
  i = 0;
  while ((dd->ndom > MIN_DOMAINS) && (i < MAX_COARSENING_STEPS)
         && ((dd->G->nedges >> 1) > dd->G->nvtx))
   { shrinkDomainDecomposition(dd, options[OPTION_NODE_SELECTION3]);
     dd = dd->next; i++;

#ifdef BE_CAUTIOUS
     checkDomainDecomposition(dd);
#endif

     if (options[OPTION_MSGLVL] > 2)
       printf("\t %2d. dom.dec.: #nodes %d (#domains %d, weight %d), #edges %d"
              "\n", i, dd->G->nvtx, dd->ndom, dd->domwght, dd->G->nedges >> 1);
   }
  pord_stoptimer(cpus[TIME_COARSEDOMDEC]);

  /* -----------------------------------------------
     determine coloring of last domain decomposition
     ------------------------------------------------ */
  pord_starttimer(cpus[TIME_INITSEP]);
  initialDDSep(dd);
  if (dd->cwght[GRAY] > 0)
    improveDDSep(dd);

#ifdef BE_CAUTIOUS
  checkDDSep(dd);
#endif

  if (options[OPTION_MSGLVL] > 2)
    printf("\t %2d. dom.dec. sep.: S %d, B %d, W %d [cost %7.2f]\n",
           i, dd->cwght[GRAY], dd->cwght[BLACK], dd->cwght[WHITE],
           F(dd->cwght[GRAY], dd->cwght[BLACK], dd->cwght[WHITE]));
  pord_stoptimer(cpus[TIME_INITSEP]);

  /* --------------
     refine coloring
     --------------- */

  pord_starttimer(cpus[TIME_REFINESEP]);
  while (dd->prev != NULL)
   { dd2 = dd->prev;
     dd2->cwght[GRAY] = dd->cwght[GRAY];
     dd2->cwght[BLACK] = dd->cwght[BLACK];
     dd2->cwght[WHITE] = dd->cwght[WHITE];
     for (u = 0; u < dd2->G->nvtx; u++)
       dd2->color[u] = dd->color[dd2->map[u]];
     freeDomainDecomposition(dd);
     if (dd2->cwght[GRAY] > 0)
       improveDDSep(dd2);

#ifdef BE_CAUTIOUS
     checkDDSep(dd2);
#endif

     dd = dd2;
     i--;
     if (options[OPTION_MSGLVL] > 2)
       printf("\t %2d. dom.dec. sep.: S %d, B %d, W %d [cost %7.2f]\n",
              i, dd->cwght[GRAY], dd->cwght[BLACK], dd->cwght[WHITE],
              F(dd->cwght[GRAY], dd->cwght[BLACK], dd->cwght[WHITE]));
   }
  pord_stoptimer(cpus[TIME_REFINESEP]);

  /* ---------------------------------
     copy coloring to subgraph Gbisect
     --------------------------------- */
  cwght[GRAY] = dd->cwght[GRAY];
  cwght[BLACK] = dd->cwght[BLACK];
  cwght[WHITE] = dd->cwght[WHITE];
  for (u = 0; u < nvtx; u++)
    color[u] = dd->color[map[u]];
  freeDomainDecomposition(dd);
  free(map);
}


/*****************************************************************************
******************************************************************************/
PORD_INT
smoothBy2Layers(gbisect_t *Gbisect, PORD_INT *bipartvertex, PORD_INT *pnX,
                PORD_INT black, PORD_INT white)
{ gbipart_t *Gbipart;
  PORD_INT       *xadj, *adjncy, *color, *cwght, *map;
  PORD_INT       *flow, *rc, *matching, *dmflag, dmwght[6];
  PORD_INT       nvtx, smoothed, nX, nX2, nY, x, y, u, i, j, jstart, jstop;

  nvtx = Gbisect->G->nvtx;
  xadj = Gbisect->G->xadj;
  adjncy = Gbisect->G->adjncy;
  color = Gbisect->color;
  cwght = Gbisect->cwght;
  nX = *pnX;

  /* ----------------------------------------------------
     map vector identifies vertices of Gbisect in Gbipart
     ---------------------------------------------------- */
  mymalloc(map, nvtx, PORD_INT);

  /* ----------------------------------
     construct set Y of bipartite graph
     ---------------------------------- */
  nY = 0;
  for (i = 0; i < nX; i++)
   { x = bipartvertex[i];
     jstart = xadj[x];
     jstop = xadj[x+1];
     for (j = jstart; j < jstop; j++)
      { y = adjncy[j];
        if (color[y] == black)
         { bipartvertex[nX+nY++] = y;
           color[y] = GRAY;
         }
      }
   }
  for (i = nX; i < nX+nY; i++)
   { y = bipartvertex[i];
     color[y] = black;
   }

  /* --------------------------------------------
     compute the Dulmage-Mendelsohn decomposition
     -------------------------------------------- */
  Gbipart = setupBipartiteGraph(Gbisect->G, bipartvertex, nX, nY, map);

  mymalloc(dmflag, (nX+nY), PORD_INT);
  switch(Gbipart->G->type)
   { case UNWEIGHTED:
       mymalloc(matching, (nX+nY), PORD_INT);
       maximumMatching(Gbipart, matching);
       DMviaMatching(Gbipart, matching, dmflag, dmwght);
       free(matching);
       break;
     case WEIGHTED:
       mymalloc(flow, Gbipart->G->nedges, PORD_INT);
       mymalloc(rc, (nX+nY), PORD_INT);
       maximumFlow(Gbipart, flow, rc);
       DMviaFlow(Gbipart, flow, rc, dmflag, dmwght);
       free(flow);
       free(rc);
       break;
     default:
       fprintf(stderr, "\nError in function smoothSeparator\n"
            "  unrecognized bipartite graph type %d\n", Gbipart->G->type);
       quit();
   }

#ifdef DEBUG
  printf("Dulmage-Mendelsohn decomp. computed\n"
         "SI %d, SX %d, SR %d, BI %d, BX %d, BR %d\n", dmwght[SI], dmwght[SX],
         dmwght[SR], dmwght[BI], dmwght[BX], dmwght[BR]);
#endif

  /* -----------------------------------------------------------------------
     1st TEST: try to exchange SI with BX, i.e. nodes in SI are moved from
     the separator into white (white grows), and nodes in BX are moved from
     black into the separator (black shrinks)
     ----------------------------------------------------------------------- */
  smoothed = FALSE;
  if (F(cwght[GRAY]-dmwght[SI]+dmwght[BX], cwght[black]-dmwght[BX],
        cwght[white]+dmwght[SI]) + EPS < F(cwght[GRAY], cwght[black],
                                           cwght[white]))
   { smoothed = TRUE;

#ifdef DEBUG
     printf("exchange SI with BX\n");
#endif

     cwght[white] += dmwght[SI]; cwght[GRAY] -= dmwght[SI];
     cwght[black] -= dmwght[BX]; cwght[GRAY] += dmwght[BX];
     for (i = 0; i < nX+nY; i++)
      { u = bipartvertex[i];
        if (dmflag[map[u]] == SI)
          color[u] = white;
        if (dmflag[map[u]] == BX)
          color[u] = GRAY;
      }
   }

  /* -----------------------------------------------------------------------
     2nd TEST: try to exchange SR with BR, i.e. nodes in SR are moved from
     the separator into white (white grows), and nodes in BR are moved from
     black into the separator (black shrinks)
     NOTE: SR is allowed to be exchanged with BR only if SI = BX = 0 or if
           SI has been exchanged with BX (Adj(SR) is a subset of BX u BR)
     ----------------------------------------------------------------------- */
  if ((F(cwght[GRAY]-dmwght[SR]+dmwght[BR], cwght[black]-dmwght[BR],
         cwght[white]+dmwght[SR]) + EPS < F(cwght[GRAY], cwght[black],
                                            cwght[white]))
      && ((smoothed) || (dmwght[SI] == 0)))
   { smoothed = TRUE;

#ifdef DEBUG
     printf("exchange SR with BR\n");
#endif

     cwght[white] += dmwght[SR]; cwght[GRAY] -= dmwght[SR];
     cwght[black] -= dmwght[BR]; cwght[GRAY] += dmwght[BR];
     for (i = 0; i < nX+nY; i++)
      { u = bipartvertex[i];
        if (dmflag[map[u]] == SR)
          color[u] = white;
        if (dmflag[map[u]] == BR)
          color[u] = GRAY;
      }
   }

  /* -----------------------------------------------------
     fill bipartvertex with the nodes of the new separator
     ----------------------------------------------------- */
  nX2 = 0;
  for (i = 0; i < nX+nY; i++)
   { u = bipartvertex[i];
     if (color[u] == GRAY)
       bipartvertex[nX2++] = u;
   }
  *pnX = nX2;

  /* -------------------------------
     free working storage and return
     ------------------------------- */
  free(map); free(dmflag);
  freeBipartiteGraph(Gbipart);
  return(smoothed);
}


/*****************************************************************************
******************************************************************************/
void
smoothSeparator(gbisect_t *Gbisect, options_t *options)
{ PORD_INT *xadj, *adjncy, *vwght, *color, *cwght, *bipartvertex;
  PORD_INT nvtx, nX, nX2, u, x, y, a, b, i, j, jstart, jstop;

  nvtx = Gbisect->G->nvtx;
  xadj = Gbisect->G->xadj;
  adjncy = Gbisect->G->adjncy;
  vwght = Gbisect->G->vwght;
  color = Gbisect->color;
  cwght = Gbisect->cwght;

  mymalloc(bipartvertex, nvtx, PORD_INT);

  /* ----------------------------------------------------------
     extract the separator (store its vertices in bipartvertex)
     ---------------------------------------------------------- */
  nX = 0;
  for (u = 0; u < nvtx; u++)
    if (color[u] == GRAY)
      bipartvertex[nX++] = u;

  do
   { /* ---------------------------------------------------------------
        minimize the separator (i.e. minimize set X of bipartite graph)
        --------------------------------------------------------------- */
     cwght[GRAY] = nX2 = 0;
     for (i = 0; i < nX; i++)
      { x = bipartvertex[i];
        a = b = FALSE;
        jstart = xadj[x];
        jstop = xadj[x+1];
        for (j = jstart; j < jstop; j++)
         { y = adjncy[j];
           if (color[y] == WHITE) a = TRUE;
           if (color[y] == BLACK) b = TRUE;
         }
        if ((a) && (!b))
         { color[x] = WHITE; cwght[WHITE] += vwght[x]; }
        else if ((!a) && (b))
         { color[x] = BLACK; cwght[BLACK] += vwght[x]; }
        else
         { bipartvertex[nX2++] = x; cwght[GRAY] += vwght[x]; }
      }
     nX = nX2;

#ifdef BE_CAUTIOUS
     checkSeparator(Gbisect);
#endif

     /* ------------------------------------------------------------------
        smooth the unweighted/weighted separator
        first pair it with the larger set; if unsuccessful try the smaller
        ------------------------------------------------------------------ */
     if (cwght[BLACK] >= cwght[WHITE])
      { a = smoothBy2Layers(Gbisect, bipartvertex, &nX, BLACK, WHITE);
        if (!a)
          a = smoothBy2Layers(Gbisect, bipartvertex, &nX, WHITE, BLACK);
      }
     else
      { a = smoothBy2Layers(Gbisect, bipartvertex, &nX, WHITE, BLACK);
        if (!a)
          a = smoothBy2Layers(Gbisect, bipartvertex, &nX, BLACK, WHITE);
      }
     if ((options[OPTION_MSGLVL] > 2) && (a))
       printf("\t separator smoothed: S %d, B %d, W %d [cost %7.2f]\n",
              cwght[GRAY], cwght[BLACK], cwght[WHITE],
              F(cwght[GRAY], cwght[BLACK], cwght[WHITE])); 
   } while (a);
     
  free(bipartvertex);
}

