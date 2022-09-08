/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: ddbisect.c
/
/ author        J"urgen Schulze, University of Paderborn
/ created       00mar09
/
/ This file contains code for the construction/improvement of a vertex
/ separator for a domain decomposition
/
******************************************************************************

Data type:  struct domdec
              graph_t *G;            pointer to graph object
              int     ndom;          number of domains
              int     domwght;       total weight of domains
              int     *vtype;        type of node (DOMAIN or MULTISEC)
              int     *color;        color of node (GRAY, BLACK, or WHITE)
              int     cwght[3];      weights of GRAY, BLACK, WHITE partitions
              int     *map;          maps nodes to next coarser domain decomp.
              struct domdec *prev;   pointer to previous finer domain decomp.
              struct domdec *next;   pointer to next coarser domain decomp.
Comments:
  o Structure holds the domain decompositions constructed by the
    coarsening process; it also holds the colorings of the domain decomp.
    computed by the refinement process
  o vtype[v]: represents the status of a node in the domain decomposition
              0, iff status of v is unknown
              1, iff v is a domain vertex
              2, iff v is a multisector vertex
              3, iff multisec v is eliminated and now forms a domain
              4, iff multisec v is absorbed by another multisec/domain
Methods in lib/ddbisect.c:
- void checkDDSep(domdec_t *dd);
- int findPseudoPeripheralDomain(domdec_t *dd, int domain);
    o returns a domain with maximal excentricity by repeated breadth first
      search; first bfs starts at node domain
- void constructLevelSep(domdec_t *dd, int domain);
    o determines a vertex separator by breadth first search starting at node
      domain;
- void initialDDSep(domdec_t *dd);
    o computes an initial separator for the domain decomposition dd;
      initially, all domains/multisecs are colored black; the function scans
      over all connected components of dd; it first calls findPseudoPeripheral-
      Domain to obtain a domain with maximal excentricity and then it calls
      constructLevelSep for that domain.
- void updateB2W(bucket_t *w_bucket, bucket_t *b_bucket, domdec_t *dd,
            int domain, int *tmp_color, int *deltaW, int *deltaB, int *deltaS);
    o if domain flips its color from BLACK to WHITE, all neighboring domains
      that share a common variable have to be updated (see my PhD thesis)
- void updateW2B(bucket_t *w_bucket, bucket_t *b_bucket, domdec_t *dd,
            int domain, int *tmp_color, int *deltaW, int *deltaB, int *deltaS);
    o if domain flips its color from WHITE to BLACK, all neighboring domains
      that share a common variable have to be updated (see my PhD thesis)
- void improveDDSep(domdec_t *dd);
    o Fiducia-Mattheyses variant to improve the coloring/separator of a
      domain decomposition (see my PhD thesis)

******************************************************************************/

#include <space.h>
/* #define DEBUG */


/******************************************************************************
******************************************************************************/
void
checkDDSep(domdec_t *dd)
{ PORD_INT *xadj, *adjncy, *vwght, *vtype, *color, *cwght;
  PORD_INT nvtx, err, u, v, i, istart, istop, nBdom, nWdom;
  PORD_INT checkS, checkB, checkW;

  nvtx = dd->G->nvtx;
  xadj = dd->G->xadj;
  adjncy = dd->G->adjncy;
  vwght = dd->G->vwght;
  vtype = dd->vtype;
  color = dd->color;
  cwght = dd->cwght;

  err = FALSE;
  printf("checking separator of domain decomposition (S %d, B %d, W %d)\n",
         cwght[GRAY], cwght[BLACK], cwght[WHITE]);

  checkS = checkB = checkW = 0;
  for (u = 0; u < nvtx; u++)
    /* check neighborhood of multisector nodes */
    if (vtype[u] == 2)
     { nBdom = nWdom = 0;
       istart = xadj[u];
       istop = xadj[u+1];
       for (i = istart; i < istop; i++)
        { v = adjncy[i];
          if (color[v] == BLACK) nBdom++;
          if (color[v] == WHITE) nWdom++;
        }
       switch(color[u])
        { case GRAY:
            checkS += vwght[u];
            if ((nBdom == 0) || (nWdom == 0))
              printf("WARNING: multisec %d belongs to S, but nBdom = %d and "
                     "nWdom = %d\n", u, nBdom, nWdom);
            break;
          case BLACK: 
            checkB += vwght[u];
            if (nWdom > 0)
             { printf("ERROR: black multisec %d adjacent to white domain\n", u);
               err = TRUE;
             }
            break;
          case WHITE:
            checkW += vwght[u];
            if (nBdom > 0)
             { printf("ERROR: white multisec %d adjacent to black domain\n", u);
               err = TRUE;
             }
            break;
          default:
            printf("ERROR: multisec %d has unrecognized color %d\n", u,
                   color[u]);
            err = TRUE;
        }
     }

    /* sum up size of white/black domains */
    else /* if (vtype[u] == 1) */
      switch(color[u])
       { case BLACK:
           checkB += vwght[u]; break;
         case WHITE:
           checkW += vwght[u]; break;
         default:
           printf("ERROR: domain %d has unrecognized color %d\n", u, color[u]);
           err = TRUE;
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
PORD_INT
findPseudoPeripheralDomain(domdec_t* dd, PORD_INT domain)
{ PORD_INT *xadj, *adjncy, *vtype, *level, *queue;
  PORD_INT nvtx, qhead, qtail, nlev, lastdomain, u, v, i, istart, istop;

  nvtx = dd->G->nvtx;
  xadj = dd->G->xadj;
  adjncy = dd->G->adjncy;
  vtype = dd->vtype;

  /* ------------------------
     allocate working storage
     ------------------------ */
  mymalloc(level, nvtx, PORD_INT);
  mymalloc(queue, nvtx, PORD_INT);

  /* ---------------------------------------
     find a domain with maximal excentricity
     --------------------------------------- */
  nlev = 0; lastdomain = domain;
  while (TRUE)
   { for (u = 0; u < nvtx; u++)
       level[u] = -1;
     queue[0] = domain; level[domain] = 0;
     qhead = 0; qtail = 1;
     while (qhead != qtail)
      { u = queue[qhead++];
        if (vtype[u] == 1)     /* remember last domain */
          lastdomain = u;
        istart = xadj[u];
        istop = xadj[u+1];
        for (i = istart; i < istop; i++)
         { v = adjncy[i];
           if (level[v] == -1)
            { queue[qtail++] = v;
              level[v] = level[u] + 1;
            }
         }
      }
     if (level[lastdomain] > nlev)
      { nlev = level[lastdomain];
        domain = lastdomain;
      }
     else break;
   }

  /* -------------------------------
     free working storage and return
     ------------------------------- */
  free(level); free(queue);
  return(domain);
}


/*****************************************************************************
*****************************************************************************/
void
constructLevelSep(domdec_t* dd, PORD_INT domain)
{ PORD_INT *xadj, *adjncy, *vwght, *vtype, *color, *cwght;
  PORD_INT *queue, *deltaS, *deltaB, *deltaW;
  PORD_INT nvtx, bestvalue, weight, qhead, qtail, qopt, q, dS, dB, dW;
  PORD_INT u, v, w, i, istart, istop, j, jstart, jstop;

  /* ======================================================================
     vtype[u]: (u domain)
         1 => domain u has not been touched yet (not in queue, no color flip)
        -1 => domain u is in queue and its deltaS, deltaB, deltaW values
              have to be updated
        -2 => domain u is in queue and no update necessary
        -3 => domain u has flipped its color to black
     deltaS[u], deltaB[u], deltaW[u]:
        u domain: denotes the change in partition size, if u flips its color
        u multisec: deltaB/deltaW denote number of adj. black/white domains
     ====================================================================== */

  nvtx = dd->G->nvtx;
  xadj = dd->G->xadj;
  adjncy = dd->G->adjncy;
  vwght = dd->G->vwght;
  vtype = dd->vtype;
  color = dd->color;
  cwght = dd->cwght;

  /* ------------------------------------------
     allocate working storage + initializations
     ------------------------------------------ */
  mymalloc(queue, nvtx, PORD_INT);
  mymalloc(deltaS, nvtx, PORD_INT);
  mymalloc(deltaB, nvtx, PORD_INT);
  mymalloc(deltaW, nvtx, PORD_INT);
  for (u = 0; u < nvtx; u++)
   { deltaS[u] = deltaB[u] = deltaW[u] = 0;
     if (vtype[u] == 2)
       deltaW[u] = xadj[u+1] - xadj[u];
   }

  /* ---------------------------------------------
     build a BFS tree rooted at domain
     the separator is given by the level structure
     --------------------------------------------- */
  queue[0] = domain;
  qhead = 0; qtail = 1;
  vtype[domain] = -1;
  while ((cwght[BLACK] < cwght[WHITE]) && (qhead != qtail))
   { qopt = 0;
     bestvalue = MAX_INT;
     
     /* --------------------------------------------------------------------
        run through queue, update domains if necessary, and find best domain
        -------------------------------------------------------------------- */
     for (q = qhead; q < qtail; q++)
      { u = queue[q];
        if (vtype[u] == -1)
         { dB = vwght[u]; dW = -dB; dS = 0;
           istart = xadj[u];
           istop = xadj[u+1];
           for (i = istart; i < istop; i++)
            { v = adjncy[i];                    /* color of multisec v */
              weight = vwght[v];                /* is GRAY or WHITE */
              if (color[v] == WHITE)
               { dW -= weight; dS += weight; }  /* multisec will move to S */
              else if (deltaW[v] == 1)
               { dB += weight; dS -= weight; }  /* multisec will move to B */
            }
           deltaS[u] = dS; deltaB[u] = dB; deltaW[u] = dW;
           vtype[u] = -2;
         }
        if (cwght[GRAY] + deltaS[u] < bestvalue)
         { bestvalue = cwght[GRAY] + deltaS[u];
           qopt = q;
         }
      }

     /* ----------------------------------------------------
        move best domain to head of queue and color it black
        ---------------------------------------------------- */
     u = queue[qopt];
     swap(queue[qopt], queue[qhead], v);
     qhead++;
     color[u] = BLACK;
     cwght[GRAY] += deltaS[u];
     cwght[BLACK] += deltaB[u];
     cwght[WHITE] += deltaW[u];
     vtype[u] = -3;

     /* ------------------------------------------------------------
        update all multisecs that are adjacent to domain u and check
        domains adjacent to the multisecs
        ------------------------------------------------------------ */
     istart = xadj[u];
     istop = xadj[u+1];
     for (i = istart; i < istop; i++)
      { v = adjncy[i];
        deltaB[v]++; deltaW[v]--;
        if (deltaW[v] == 0)       /* color of multisec v changed to BLACK */
          color[v] = BLACK;
        else if (deltaB[v] == 1)  /* color of multisec v changed to GRAY */
         { color[v] = GRAY;
           jstart = xadj[v];
           jstop = xadj[v+1];
           for (j = jstart; j < jstop; j++)
            { w = adjncy[j];
              if (vtype[w] == 1)           /* a new domain enters the queue */
               { queue[qtail++] = w;
                 vtype[w] = -1;
               }
              else if (vtype[w] == -2)     /* update (old) domain in queue */
                vtype[w] = -1;
            }
         }
        else if (deltaW[v] == 1)  /* color of multisec v remains GRAY for */
         { jstart = xadj[v];      /* the last time */
           jstop = xadj[v+1];
           for (j = jstart; j < jstop; j++)
            { w = adjncy[j];
              if (vtype[w] == -2)
                vtype[w] = -1;
            }
         }
      }
   }

  /* ---------------------------
     reset vtype and free memory
     --------------------------- */
  for (i = 0; i < qtail; i++)
   { u = queue[i];
     vtype[u] = 1;
   }
  free(queue);
  free(deltaS); free(deltaB); free(deltaW);
}


/*****************************************************************************
******************************************************************************/
void
initialDDSep(domdec_t *dd)
{  PORD_INT *vtype, *color, *cwght;
   PORD_INT nvtx, totvwght, domain, u;

   nvtx = dd->G->nvtx;
   totvwght = dd->G->totvwght;
   vtype = dd->vtype;
   color = dd->color;
   cwght = dd->cwght;

  /* --------------------------------------------------------
     initializations (all nodes are colored white by default)
     -------------------------------------------------------- */
  cwght[GRAY] = 0;
  cwght[BLACK] = 0;
  cwght[WHITE] = totvwght;
  for (u = 0; u < nvtx; u++)
    color[u] = WHITE;

  /* ----------------------------------------------------------------------
     scan over connected components and create level based vertex separator
     ---------------------------------------------------------------------- */
  for (u = 0; u < nvtx; u++)
    if ((vtype[u] == 1) && (color[u] == WHITE))
     { domain = findPseudoPeripheralDomain(dd, u);
       constructLevelSep(dd, domain);
       if (cwght[BLACK] >= cwght[WHITE])
         break;
     }
}


/*****************************************************************************
*****************************************************************************/
void
updateB2W(bucket_t *w_bucket, bucket_t *b_bucket, domdec_t *dd, PORD_INT domain,
          PORD_INT *tmp_color, PORD_INT *deltaW, PORD_INT *deltaB, PORD_INT *deltaS)
{ PORD_INT *xadj, *adjncy, *vwght, *vtype;
  PORD_INT weight, u, v, i, istart, istop, j, jstart, jstop;

  xadj = dd->G->xadj;
  adjncy = dd->G->adjncy;
  vwght = dd->G->vwght;
  vtype = dd->vtype;
 
  istart = xadj[domain];
  istop = xadj[domain+1];
  for (i = istart; i < istop; i++)
   { u = adjncy[i];
     weight = vwght[u];
     jstart = xadj[u];
     jstop = xadj[u+1];

     /* ---------------------------------------------------------------
        subcase (1): before flipping domain to WHITE there was only one
          other WHITE domain v. update deltaB[v] and deltaS[v]
        --------------------------------------------------------------- */
     if (deltaW[u] < 0)
      { v = -(deltaW[u]+1);
        deltaW[u] = 1;

#ifdef DEBUG
        printf(" B2W case (1): (via multisec %d) removing domain %d from "
               "w_bucket\n", u, v);
#endif

        removeBucket(w_bucket, v);
        deltaB[v] -= weight; deltaS[v] += weight;
        insertBucket(w_bucket, deltaS[v], v);
      }

     /* ---------------------------------------------------------------
        subcase (2): all other domains are BLACK. update deltaB, deltaS
          of these BLACK domains. NOTE: subcase (3) may directly follow
        --------------------------------------------------------------- */
     if (deltaW[u] == 0)
      { tmp_color[u] = GRAY;
        for (j = jstart; j < jstop; j++)
         { v = adjncy[j];
           if (vtype[v] == 1)
            {
#ifdef DEBUG
              printf(" B2W case (2): (via multisec %d) removing domain %d from "
                     "b_bucket\n", u, v);
#endif

              removeBucket(b_bucket, v);
              deltaB[v] += weight; deltaS[v] -= weight;
              insertBucket(b_bucket, deltaS[v], v);
            }
         }
      }

     if (deltaB[u] < 0) deltaB[u] = 1;   /* the unique BLACK dom. flipped */
     deltaB[u]--; deltaW[u]++;

     /* -------------------------------------------------------------
        subcase (3): after flipping domain to WHITE there is only one
          remaining BLACK domain. search it and update deltaW, deltaS
          furthermore, store the remaining BLACK domain in deltaB[u]
        ------------------------------------------------------------- */
     if (deltaB[u] == 1)
      { for (j = jstart; j < jstop; j++)
         { v = adjncy[j];
           if ((tmp_color[v] == BLACK) && (vtype[v] == 1))
            {
#ifdef DEBUG
              printf(" B2W case (3): (via multisec %d) removing domain %d from "
                     "b_bucket\n", u, v);
#endif

              removeBucket(b_bucket, v);
              deltaW[v] += weight; deltaS[v] -= weight;
              deltaB[u] = -(v+1);
              insertBucket(b_bucket, deltaS[v], v);
            }
         }
      }

     /* -------------------------------------------------------------
        subcase (4): after flipping domain to WHITE there is no other
          BLACK domain. update deltaW, deltaS of the WHITE domains
        ------------------------------------------------------------- */
     if (deltaB[u] == 0)
      { tmp_color[u] = WHITE;
        for (j = jstart; j < jstop; j++)
         { v = adjncy[j];
           if (vtype[v] == 1)
            {
#ifdef DEBUG
              printf(" B2W case (4): (via multisec %d) removing domain %d from "
                     "w_bucket\n", u, v);
#endif

              removeBucket(w_bucket, v);
              deltaW[v] -= weight; deltaS[v] += weight;
              insertBucket(w_bucket, deltaS[v], v);
            }
         }
      }
   }
}


/*****************************************************************************
*****************************************************************************/
void
updateW2B(bucket_t *w_bucket, bucket_t *b_bucket, domdec_t *dd, PORD_INT domain,
          PORD_INT *tmp_color, PORD_INT *deltaW, PORD_INT *deltaB, PORD_INT *deltaS)
{ PORD_INT *xadj, *adjncy, *vwght, *vtype;
  PORD_INT weight, u, v, i, istart, istop, j, jstart, jstop;

  xadj = dd->G->xadj;
  adjncy = dd->G->adjncy;
  vwght = dd->G->vwght;
  vtype = dd->vtype;

  istart = xadj[domain];
  istop = xadj[domain+1];
  for (i = istart; i < istop; i++)
   { u = adjncy[i];
     weight = vwght[u];
     jstart = xadj[u];
     jstop = xadj[u+1];

     /* ---------------------------------------------------------------
        subcase (1): before flipping domain to BLACK there was only one
          other BLACK domain v. update deltaW[v] and deltaS[v]
        --------------------------------------------------------------- */
     if (deltaB[u] < 0)
      { v = -(deltaB[u]+1);
        deltaB[u] = 1;

#ifdef DEBUG
        printf(" W2B case (1): (via multisec %d) removing domain %d from "
               "b_bucket\n", u, v);
#endif

        removeBucket(b_bucket, v);
        deltaW[v] -= weight; deltaS[v] += weight;
        insertBucket(b_bucket, deltaS[v], v);
      }

     /* ---------------------------------------------------------------
        subcase (2): all other domains are WHITE. update deltaW, deltaS
          of these WHITE domains. NOTE: subcase (3) may directly follow
        --------------------------------------------------------------- */
     if (deltaB[u] == 0)
      { tmp_color[u] = GRAY;
        for (j = jstart; j < jstop; j++)
         { v = adjncy[j];
           if (vtype[v] == 1)
            {
#ifdef DEBUG
              printf(" W2B case (2): (via multisec %d) removing domain %d from "
                     "w_bucket\n", u, v);
#endif

              removeBucket(w_bucket, v);
              deltaW[v] += weight; deltaS[v] -= weight;
              insertBucket(w_bucket, deltaS[v], v);
            }
         }
      }

     if (deltaW[u] < 0) deltaW[u] = 1;   /* the unique WHITE dom. flipped */
     deltaB[u]++; deltaW[u]--;

     /* -------------------------------------------------------------
        subcase (3): after flipping domain to BLACK there is only one
          remaining WHITE domain. search it and update deltaB, deltaS
          furthermore, store the remaining WHITE domain in deltaW[u]
        ------------------------------------------------------------- */
     if (deltaW[u] == 1)
      { for (j = jstart; j < jstop; j++)
         { v = adjncy[j];
           if ((tmp_color[v] == WHITE) && (vtype[v] == 1))
            {
#ifdef DEBUG
              printf(" W2B case (3): (via multisec %d) removing domain %d from "
                     "w_bucket\n", u, v);
#endif

              removeBucket(w_bucket, v);
              deltaB[v] += weight; deltaS[v] -= weight;
              deltaW[u] = -(v+1);
              insertBucket(w_bucket, deltaS[v], v);
            }
         }
      }

     /* ---------------------------------------------------------------
        subcase (4): after flipping domain to BLACK there is no other
          WHITE domain. update deltaB, deltaS of the BLACK domains
        --------------------------------------------------------------- */
     if (deltaW[u] == 0)
      { tmp_color[u] = BLACK;
        for (j = jstart; j < jstop; j++)
         { v = adjncy[j];
           if (vtype[v] == 1)
            {
#ifdef DEBUG
              printf(" W2B case (4): (via multisec %d) removing domain %d from "
                     "b_bucket\n", u, v);
#endif

              removeBucket(b_bucket, v);
              deltaB[v] -= weight; deltaS[v] += weight;
              insertBucket(b_bucket, deltaS[v], v);
            }
         }
      }
   }
}


/*****************************************************************************
******************************************************************************/
void
improveDDSep(domdec_t *dd)
{ bucket_t *b_bucket, *w_bucket;
  PORD_INT      *xadj, *adjncy, *vwght, *vtype, *color, *cwght;
  PORD_INT      *tmp_color, *deltaS, *deltaB, *deltaW;
  PORD_INT      nvtx, weight, tmp_S, tmp_B, tmp_W;
  PORD_INT      pos, bestglobalpos, badflips, b_domain, w_domain, domain, nxtdomain;
  PORD_INT      fhead, ftail, u, v, i, istart, istop;
  FLOAT    bestglobalvalue, b_value, w_value, value;

  /* ======================================================================
     vtype[u]: (u domain)
         1 => color of domain u has not been changed
       < 0 => points to next domain in flipping list
              (fhead points to first, ftail points to last domain in list)
       = 0 => domain is last domain in flipping list
     ====================================================================== */

  nvtx = dd->G->nvtx;
  xadj = dd->G->xadj;
  adjncy = dd->G->adjncy;
  vwght = dd->G->vwght;
  vtype = dd->vtype;
  color = dd->color;
  cwght = dd->cwght;

  mymalloc(tmp_color, nvtx, PORD_INT);
  mymalloc(deltaS, nvtx, PORD_INT);
  mymalloc(deltaB, nvtx, PORD_INT);
  mymalloc(deltaW, nvtx, PORD_INT);

OUTER_LOOP_START:

  /* ----------------------------------------------------------------------
     copy data of actual bisection and initialize buckets and flipping list
     ---------------------------------------------------------------------- */
  tmp_S = cwght[GRAY];
  tmp_B = cwght[BLACK];
  tmp_W = cwght[WHITE];
  bestglobalpos = badflips = 0;
  bestglobalvalue = F(tmp_S, tmp_B, tmp_W);

  b_bucket = setupBucket(nvtx, nvtx, (nvtx >> 1));
  w_bucket = setupBucket(nvtx, nvtx, (nvtx >> 1));

  fhead = 0; ftail = -1;
  pos = 0;

  /* ----------------------------------------------------------
     initialize tmp_color, deltaB, and deltaW for all multisecs
     ---------------------------------------------------------- */
  for (u = 0; u < nvtx; u++)
    if (vtype[u] == 2)
     { deltaB[u] = deltaW[u] = 0;
       istart = xadj[u];
       istop = xadj[u+1];
       for (i = istart; i < istop; i++)
        { v = adjncy[i];
          if (color[v] == BLACK) deltaB[u]++;
          else deltaW[u]++;
        }
       if ((deltaB[u] > 0) && (deltaW[u] > 0))  /* update multisec coloring */
         tmp_color[u] = GRAY;
       else if (deltaB[u] > 0) tmp_color[u] = BLACK;
       else tmp_color[u] = WHITE;
       color[u] = tmp_color[u];
     }

  /* -----------------------------------------------------------------
     initialize tmp_color, deltaS,B,W for all domains and fill buckets
     ----------------------------------------------------------------- */
  for (u = 0; u < nvtx; u++)
    if (vtype[u] == 1)
     { tmp_color[u] = color[u];
       if (tmp_color[u] == BLACK)          /* domain may be flipped to WHITE */
        { deltaW[u] = vwght[u]; deltaB[u] = -deltaW[u]; deltaS[u] = 0;
          istart = xadj[u];
          istop = xadj[u+1];
          for (i = istart; i < istop; i++)
           { v = adjncy[i];                /* tmp_color[v] e {GRAY, BLACK} */
             weight = vwght[v];
             if (tmp_color[v] == BLACK)    /* multisec v will move into S */
              { deltaB[u] -= weight;
                deltaS[u] += weight;
              }
             else if (deltaB[v] == 1)      /* multisec v will move into W */
              { deltaW[u] += weight;
                deltaS[u] -= weight;
                deltaB[v] = -(u+1);
              }
           }
          insertBucket(b_bucket, deltaS[u], u);
        }
       if (tmp_color[u] == WHITE)          /* domain may be flipped to BLACK */
        { deltaB[u] = vwght[u]; deltaW[u] = -deltaB[u]; deltaS[u] = 0;
          istart = xadj[u];
          istop = xadj[u+1];
          for (i = istart; i < istop; i++)
           { v = adjncy[i];                /* tmp_color[v] e {GRAY, WHITE} */
             weight = vwght[v];
             if (tmp_color[v] == WHITE)    /* multisec v will move into S */
              { deltaW[u] -= weight;
                deltaS[u] += weight;
              }
             else if (deltaW[v] == 1)      /* multisec v will move into B */
              { deltaB[u] += weight;
                deltaS[u] -= weight;
                deltaW[v] = -(u+1);
              }
           }
          insertBucket(w_bucket, deltaS[u], u);
        }
     }

#ifdef DEBUG
  printf("starting inner loop: b_bucket->nobj %d, w_bucket->nobj %d\n",
         b_bucket->nobj, w_bucket->nobj);
  waitkey();
#endif

INNER_LOOP_START:

  /* -------------------------------------------
     extract best domain from b_bucket, w_bucket
     ------------------------------------------- */
  b_value = w_value = MAX_FLOAT;
  if ((b_domain = minBucket(b_bucket)) != -1)
   { b_value = F((tmp_S+deltaS[b_domain]), (tmp_B+deltaB[b_domain]),
                 (tmp_W+deltaW[b_domain]));

#ifdef DEBUG
     printf("best black domain: %d, deltaS %d, deltaB %d, deltaW %d, "
            "cost %7.2f\n", b_domain, deltaS[b_domain], deltaB[b_domain],
            deltaW[b_domain], b_value);
#endif
   }
  if ((w_domain = minBucket(w_bucket)) != -1)
   { w_value = F((tmp_S+deltaS[w_domain]), (tmp_B+deltaB[w_domain]),
                 (tmp_W+deltaW[w_domain]));

#ifdef DEBUG
     printf("best white domain: %d, deltaS %d, deltaB %d, deltaW %d, "
            "cost %7.2f\n", w_domain, deltaS[w_domain], deltaB[w_domain],
            deltaW[w_domain], w_value);
#endif
   }

  if ((b_domain == ERR) && (w_domain == ERR)) goto INNER_LOOP_END;

  if (b_value + EPS < w_value)
   { domain = b_domain; value = b_value;
     removeBucket(b_bucket, domain);
   }
  else
   { domain = w_domain; value = w_value;
     removeBucket(w_bucket, domain);
   }

#ifdef DEBUG
  printf(" domain %d removed from bucket\n", domain);
#endif

  /* -------------------------------------------------------------------
     flip the color of domain and put it in list of log. flipped domains
     ------------------------------------------------------------------- */
  if (ftail != -1)
    vtype[ftail] = -(domain+1);   /* append domain */
  else fhead = -(domain+1);       /* list starts with domain */
  vtype[domain] = 0;              /* mark end of list */
  ftail = domain;                 /* domain is last element in list */

  if (tmp_color[domain] == BLACK)
   { tmp_color[domain] = WHITE;
     updateB2W(w_bucket,b_bucket,dd,domain,tmp_color,deltaW,deltaB,deltaS);
   }
  else if (tmp_color[domain] == WHITE)
   { tmp_color[domain] = BLACK;
     updateW2B(w_bucket,b_bucket,dd,domain,tmp_color,deltaW,deltaB,deltaS);
   }
  tmp_S += deltaS[domain];
  tmp_B += deltaB[domain];
  tmp_W += deltaW[domain];

  pos++;
  if (value + EPS < bestglobalvalue)
   { bestglobalvalue = value;
     bestglobalpos = pos;
     badflips = 0;
   }
  else badflips++;
  if (badflips < MAX_BAD_FLIPS) goto INNER_LOOP_START;

INNER_LOOP_END:

  /* --------------------------------------------
     end of inner loop: now do the physical flips
     -------------------------------------------- */
  pos = 0;
  nxtdomain = fhead;
  while (nxtdomain != 0)
   { domain = -nxtdomain - 1;
     if (pos < bestglobalpos)
      { if (color[domain] == BLACK) color[domain] = WHITE;
        else color[domain] = BLACK;
        cwght[GRAY] += deltaS[domain];
        cwght[BLACK] += deltaB[domain];
        cwght[WHITE] += deltaW[domain];
        pos++;
      }
     nxtdomain = vtype[domain];
     vtype[domain] = 1;
   }

  /* ----------------------------------------------
     partition improved => re-start the whole stuff
     ---------------------------------------------- */
#ifdef DEBUG
  printf(" INNER_LOOP_END (#pyhs. flips %d): S %d, B %d, W %d (%7.2f)\n",
         bestglobalpos, cwght[GRAY], cwght[BLACK], cwght[WHITE],
         bestglobalvalue);
  waitkey();
#endif

  /* JY: moved next instruction after the two
   *     freeBucket instructions because
   *     this was the cause of a memory leak.
   * if (bestglobalpos > 0) goto OUTER_LOOP_START;
   */

  freeBucket(b_bucket);
  freeBucket(w_bucket);

  if (bestglobalpos > 0) goto OUTER_LOOP_START;
  free(tmp_color); free(deltaS); free(deltaB); free(deltaW);
}

