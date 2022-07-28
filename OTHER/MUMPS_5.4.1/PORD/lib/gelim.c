/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: gelim.c
/
/ author        J"urgen Schulze, University of Paderborn
/ created       01jan10
/
/ This file contains functions dealing with the elimination graph object
/
******************************************************************************

Data type: struct gelim
             graph_t *G;        pointer to graph object
             int     maxedges;  max number of edges that can be stored
             int     *len;      length of v's adjacency list
             int     *elen;     number of elements adjacent to v
             int     *parent;   parent in front tree / representative of v
             int     *degree;   boundary size / (approximate) degree
             int     *score;    holds the score of uneliminated vertex v
Comments:
  o Structure used to hold the elimination graphs of a bottom-up ordering
  o G->totvwght: total weight of all uneliminated vertices
  o G->xadj[v] = -1 => there is no adjacency list for variable/element v
              => variable v has degree 0 (in this case G->vwght[v] > 0)
              => variable v istinguishable/removed by mass elimination
                 or element v has been absorbed (in this case G->vwght[v] = 0)
  o G->vwght[v]: weight of the princial variable v; if v becomes an element,
                 weight[v] remains unchanged for the rest of the elim. process
                 = 0 => variable v is nonprincipal/removed by mass elimination
  o len[v], elen[v]: the adjacency list of vertex/element v contains len[v]
                     entries; the first elen[v] entries are elements
                     (if v is an element, then elen[v] = 0 will hold)
  o parent[v]: for an (absorbed) element, parent[v] points to the parent of
               element v in the front tree; for an indistinguishable vertex,
               parent[v] points to its representative vertex (which may have
               also found to be indistinguishable to another one)
  o degree[v]: for an uneliminated vertex, the (approximate) degree in Gelim;
               for an element, the weight of its boundary (i.e. degree[v]
               gives the exakt degree of v at the time of its elimination)
  o score[v]: vertices are eliminated according to their score value >= 0;
              additionally, the score vector is used to represent the status
              of a node in the actual stage:
              -1, iff variable v will be eliminated in an upcomming stage
              -2, iff variable v is nonprincipal/removed by mass elim.
              -3, iff variable v has been eliminated and now forms an element
              -4, iff element v has been absorbed
Methods in lib/gelim.c
- Gelim = newElimGraph(int nvtx, int nedges);
- void freeElimGraph(gelim_t *Gelim);
- void printElimGraph(gelim_t *Gelim);
- Gelim = setupElimGraph(graph_t *G);
    o allocates memory for the elimination graph by calling newElimGraph and 
      initializes the vectors, i.e. len[u] = xadj[u+1]-xadj[u]; elen[u] = 0; 
      parent[u] = -1; degree[u] = exact (external) degree of vertex u;
      score[u] = -1; xadj[u] = -1, if len[u] = 0
- int crunchElimGraph(gelim_t *Gelim);
    o tries to compress the adjacency vector
      on success the function return TRUE, otherwise FALSE
- void buildElement(gelim_t *Gelim, int me);
    o turns variable me into an element; if me is an leaf, the element is
      constructed in-place, otherwise its adjacency list is appended to G
    o all relevant vectors are updated, i.e.
      vwght[me] = 0, degree[me] = |Lme|, score[me] = -3
      for all neighboring elements: parent[e] = me, score[e] = -4
- void updateAdjncy(gelim_t *Gelim, int *reachset, int nreach, int *tmp,
                    int *pflag);
    o updates the adjacency structure of all vertices in reachset
      IMPORTANT REQUIREMENTS:
      (1) all values stored in tmp[u] are smaller than *pflag
- void findIndNodes(gelim_t *Gelim, int *reachset, int nreach, int *bin,
                    int *next, int *tmp, int *pflag);
    o searches reachset for indistinguishable vertices
      IMPORTANT REQUIREMENTS:
       (1) the adjacency lists of all vertices in reachset have been updated
           by a call to updateAdjncy
       (2) bin[i] = -1 for all 0 <= i < G->nvtx
       (3) all values stored in tmp[u] are smaller than *pflag
    o on return bin[i] = -1 holds again
- void updateDegree(gelim_t *Gelim, int *reachset, int nreach, int *bin);
    o computes new approximate degrees for all vertices in reachset
      IMPORTANT REQUIREMENTS:
       (1) the adjacency lists of all vertices in reachset have been updated
           by a call to updateAdjncy
       (2) the boundary size of each newly formed element has been computed
       (3) bin[i] = -1 for all 0 <= i < G->nvtx
    o on return bin[i] = -1 holds again
- void updateScore(gelim_t *Gelim, int *reachset, int nreach, int scoretype,
                   int *bin);
    o updates the score of all vertices in reachset
      IMPORTANT REQUIREMENTS:
       (1) the approximate degrees are correctly computed (by updateDegree)
       (2) bin[i] = -1 for all 0 <= i < G->nvtx
    o on return bin[i] = -1 holds again
- T = extractElimTree(gelim_t *Gelim);
    o uses the status of the nodes (stored in the score vector) and the
      parent vector to set up the elimination tree T; vectors T->ncolfactor
      and T->ncolupdate are initialized using vectors G->vwght and degree

******************************************************************************/

#include <space.h>
/* #define DEBUG */


/*****************************************************************************
******************************************************************************/
gelim_t*
newElimGraph(PORD_INT nvtx, PORD_INT nedges)
{ gelim_t *Gelim;

  mymalloc(Gelim, 1, gelim_t);
  Gelim->G = newGraph(nvtx, nedges);
  Gelim->maxedges = nedges;

  mymalloc(Gelim->len, nvtx, PORD_INT);
  mymalloc(Gelim->elen, nvtx, PORD_INT);
  mymalloc(Gelim->parent, nvtx, PORD_INT);
  mymalloc(Gelim->degree, nvtx, PORD_INT);
  mymalloc(Gelim->score, nvtx, PORD_INT);

  return(Gelim);
}


/*****************************************************************************
******************************************************************************/
void
freeElimGraph(gelim_t *Gelim)
{
  freeGraph(Gelim->G);
  free(Gelim->len);
  free(Gelim->elen);
  free(Gelim->parent);
  free(Gelim->degree);
  free(Gelim->score);
  free(Gelim);
}


/*****************************************************************************
******************************************************************************/
void
printElimGraph(gelim_t *Gelim)
{ graph_t *G;
  PORD_INT     count, u, v, i, istart;

  G = Gelim->G;
  for (u = 0; u < G->nvtx; u++)
   { istart = G->xadj[u];

     /* ---------------------------------------------------------------
        case 1: u is a principal variable
          => vwght[u]: weight of all mapped indistinguishable variables
          => degree[u]: approximate degree
        ---------------------------------------------------------------- */
     if ((Gelim->score[u] == -1) || (Gelim->score[u] >= 0))
      { printf("--- adjacency list of variable %d (weight %d, degree %d, "
               "score %d):\n", u, G->vwght[u], Gelim->degree[u],
               Gelim->score[u]);
        printf("elements:\n");
        count = 0;
        for (i = istart; i < istart + Gelim->elen[u]; i++)
         { printf("%5d", G->adjncy[i]);
           if ((++count % 16) == 0)
             printf("\n");
         }
        if ((count % 16) != 0)
          printf("\n");
        printf("variables:\n");
        count = 0;
        for (i = istart + Gelim->elen[u]; i < istart + Gelim->len[u]; i++)
         { printf("%5d", G->adjncy[i]);
           if ((++count % 16) == 0)
             printf("\n");
         }
        if ((count % 16) != 0)
          printf("\n");
      }

     /* ---------------------------------------------------------------
        case 2: u is nonprincipal/removed by mass elimination
        ---------------------------------------------------------------- */
     else if (Gelim->score[u] == -2)
       printf("--- variable %d is nonprincipal/removed by mass elim. "
              "(parent %d)\n", u, Gelim->parent[u]);

     /* -----------------------------------------------
        case 3: u is an element:
          => degree[u]: weight of boundary
        ----------------------------------------------- */
     else if (Gelim->score[u] == -3)
      { printf("--- boundary of element %d (degree %d, score %d):"
               "\n", u, Gelim->degree[u], Gelim->score[u]);
        count = 0;
        for (i = istart; i < istart + Gelim->len[u]; i++)
         { v = G->adjncy[i];
           if (G->vwght[v] > 0)
            { printf("%5d", G->adjncy[i]);
              if ((++count % 16) == 0)
                printf("\n");
            }
         }
        if ((count % 16) != 0)
          printf("\n");
      }

     /* --------------------------------
        case 4: u is an absorbed element
        -------------------------------- */
     else if (Gelim->score[u] == -4)
       printf("--- element %d has been absorbed (parent %d)\n", u,
              Gelim->parent[u]);

     /* ----------------------------------------
        none of the above cases is true => error
        ---------------------------------------- */
     else
      { fprintf(stderr, "\nError in function printElimGraph\n"
             "  node %d has invalid score %d\n", u, Gelim->score[u]);
        quit();
      }
   }
}


/*****************************************************************************
******************************************************************************/
gelim_t*
setupElimGraph(graph_t *G)
{ gelim_t *Gelim;
  PORD_INT     *xadj, *adjncy, *vwght, *xadjGelim, *adjncyGelim, *vwghtGelim;
  PORD_INT     *len, *elen, *parent, *degree, *score;
  PORD_INT     nvtx, nedges, deg, u, i, istart, istop;

  nvtx = G->nvtx;
  nedges = G->nedges;
  xadj = G->xadj;
  adjncy = G->adjncy;
  vwght = G->vwght;

  Gelim = newElimGraph(nvtx, nedges+nvtx);
  xadjGelim = Gelim->G->xadj;
  adjncyGelim = Gelim->G->adjncy;
  vwghtGelim = Gelim->G->vwght;
  len = Gelim->len;
  elen = Gelim->elen;
  parent = Gelim->parent;
  degree = Gelim->degree;
  score = Gelim->score;

  /* --------------
     copy the graph
     -------------- */
  Gelim->G->type = G->type;
  Gelim->G->totvwght = G->totvwght;
  for (u = 0; u < nvtx; u++)
   { xadjGelim[u] = xadj[u];
     vwghtGelim[u] = vwght[u];
   }
  xadjGelim[nvtx] = xadj[nvtx];
  for (i = 0; i < nedges; i++)
    adjncyGelim[i] = adjncy[i];
  Gelim->G->nedges = nedges;

  /* ----------------------
     initialize all vectors
     ---------------------- */
  for (u = 0; u < nvtx; u++)
   { istart = xadj[u];
     istop = xadj[u+1];
     len[u] = istop - istart;
     elen[u] = 0;
     parent[u] = -1;
     deg = 0;

     switch(Gelim->G->type)  /* compute the external degree of u */
      { case UNWEIGHTED:
          deg = len[u];
          break;
        case WEIGHTED:
          for (i = istart; i < istop; i++)
            deg += vwght[adjncy[i]];
          break;
        default:
          fprintf(stderr, "\nError in function setupElimGraph\n"
               "  unrecognized graph type %d\n", Gelim->G->type);
      }
     degree[u] = deg;

     if (len[u] == 0)        /* len(u) = 0 => adjncy list of u not in use */
       xadjGelim[u] = -1;    /* mark with -1, otherwise crunchElimGraph fails */
     score[u] = -1;
   }

  return(Gelim);
}


/*****************************************************************************
******************************************************************************/
PORD_INT
crunchElimGraph(gelim_t *Gelim)
{ PORD_INT *xadj, *adjncy, *len;
  PORD_INT nvtx, nedges, u, i, isrc, idest;

  nvtx = Gelim->G->nvtx;
  nedges = Gelim->G->nedges;
  xadj = Gelim->G->xadj;
  adjncy = Gelim->G->adjncy;
  len = Gelim->len;

  /* ---------------------------------------------
     mark begining of u's adjacency list by -(u+1)
     --------------------------------------------- */
  for (u = 0; u < nvtx; u++)
   { i = xadj[u];                /* is adjacency list of u still in use? */
     if (i != -1)                /* verify that list is non-empty */
      { if (len[u] == 0)
         { fprintf(stderr, "\nError in function crunchElimGraph\n"
                "  adjacency list of node %d is empty\n", u);
           quit();
         }
        xadj[u] = adjncy[i];     /* if so, move first item to xadj[u] */
        adjncy[i] = -(u+1);      /* u's adjacency list is headed by -(u+1) */
        if (len[u] == 0)
          printf("error: u %d, len %d\n", u, len[u]);
      }
   }

  /* --------------------------
     crunch all adjacency lists
     -------------------------- */
  idest = isrc = 0;
  while (isrc < Gelim->G->nedges)
   { u = adjncy[isrc++];
     if (u < 0)                        /* a new adjacency list starts here */
      { u = -u - 1;                    /* it's the adjacency list of u */
        adjncy[idest] = xadj[u];       /* first item was stored in xadj[u] */
        xadj[u] = idest++;
        for (i = 1; i < len[u]; i++)
          adjncy[idest++] = adjncy[isrc++];
      }
   }
  Gelim->G->nedges = idest;

  /* ------------------
     was it successful?
     ------------------ */
  if (idest < nedges) return(TRUE);
  else return (FALSE);
}


/*****************************************************************************
******************************************************************************/
void
buildElement(gelim_t *Gelim, PORD_INT me)
{ graph_t *G;
  PORD_INT     *xadj, *adjncy, *vwght, *len, *elen, *parent, *degree, *score;
  PORD_INT     degme, elenme, vlenme, mesrcptr, medeststart, medeststart2;
  PORD_INT     medestptr, ln, p, i, j, v, e;

  G = Gelim->G;
  xadj = G->xadj;
  adjncy = G->adjncy;
  vwght = G->vwght;
  len = Gelim->len;
  elen = Gelim->elen;
  parent = Gelim->parent;
  degree = Gelim->degree;
  score = Gelim->score;

  /* ---------------------------------
     construct boundary of element Lme
     --------------------------------- */
  degme = 0;
  G->totvwght -= vwght[me];  /* me eliminated => reduce weight of Gelim */
  vwght[me] = -vwght[me];
  score[me] = -3;            /* variable me becomes an element */

  elenme = elen[me];
  vlenme = len[me] - elenme;
  mesrcptr = xadj[me];

  /* -----------------------------------------------------------
     if me is a leaf => its boundary can be constructed in-place
     ----------------------------------------------------------- */
  if (elenme == 0)
   { medeststart = xadj[me];         /* Lme overwrites old variable */
     medestptr = medeststart;        /* boundary of Lme starts here */
     for (i = 0; i < vlenme; i++)
      { v = adjncy[mesrcptr++];
        if (vwght[v] > 0)            /* v not yet placed in boundary */
         { degme += vwght[v];        /* increase size of Lme */
           vwght[v] = -vwght[v];     /* flag v as being in Lme */
           adjncy[medestptr++] = v;
         }
      }
   }

  /* -------------------------------------------------------------------
     me is not a leaf => its boundary must be constructed in empty space
     ------------------------------------------------------------------- */
  else
   { medeststart = G->nedges;        /* Lme appended to graph */
     medestptr = medeststart;        /* boundary of Lme starts here */
     for (i = 0; i <= elenme; i++)
      { if (i < elenme)              /* working on elements */
         { len[me]--;
           e = adjncy[mesrcptr++];   /* merge boundary of element e with Lme */
           p = xadj[e];              /* adjacency list of e starts here */
           ln = len[e];
         }
        else
         { e = me;                   /* merge uncovered variables with Lme */
           p = mesrcptr;             /* variables start here */
           ln = vlenme;
         }
        for (j = 0; j < ln; j++)
         { len[e]--;                 /* pick next variable, decrease length */
           v = adjncy[p++];
           if (vwght[v] > 0)
            { degme += vwght[v];     /* increase size of Lme */
              vwght[v] = -vwght[v];  /* flag v as being in Lme */

              /* ------------------------------------------
                 add v to Lme, compress adjncy if necessary
                 ------------------------------------------ */
              if (medestptr == Gelim->maxedges)
               { if (len[me] == 0) xadj[me] = -1;
                 else xadj[me] = mesrcptr;
                 if (len[e] == 0) xadj[e] = -1;
                 else xadj[e] = p;

                 /* crunch adjacency list -- !!!we need more memory!!! */
                 if (!crunchElimGraph(Gelim))
                  { fprintf(stderr, "\nError in function buildElement\n"
                         "  unable to construct element (not enough memory)\n");
                    quit();
                  }

                 /* crunch partially constructed element me */
                 medeststart2 = G->nedges;
                 for (p = medeststart; p < medestptr; p++)
                   adjncy[G->nedges++] = adjncy[p];
                 medeststart = medeststart2;
                 medestptr = G->nedges;

                 mesrcptr = xadj[me];
                 p = xadj[e];
               }
              adjncy[medestptr++] = v;
            }
         }

        /* ----------------------
           mark absorbed elements
           ---------------------- */
        if (e != me)
         { xadj[e] = -1;
           parent[e] = me;
           score[e] = -4;
         }
      }

     G->nedges = medestptr;          /* new element Lme ends here */
   }

  /* -----------------------------------
     element me successfully constructed
     ----------------------------------- */
  degree[me] = degme;
  xadj[me] = medeststart;
  vwght[me] = -vwght[me];
  elen[me] = 0;
  len[me] = medestptr - medeststart;
  if (len[me] == 0)
    xadj[me] = -1;

  /* ---------------------------
     unmark all variables in Lme
     --------------------------- */
  mesrcptr = xadj[me];
  vlenme = len[me];
  for (i = 0; i < vlenme; i++)
   { v = adjncy[mesrcptr++];
     vwght[v] = -vwght[v];
   }
}


/*****************************************************************************
******************************************************************************/
void
updateAdjncy(gelim_t *Gelim, PORD_INT *reachset, PORD_INT nreach, PORD_INT *tmp, PORD_INT *pflag)
{ PORD_INT *xadj, *adjncy, *vwght, *len, *elen, *parent, *score;
  PORD_INT u, v, e, me, i, j, jj, jdest, jfirstolde, jfirstv, jstart, jstop;
  PORD_INT covered, marku;

  xadj = Gelim->G->xadj;
  adjncy = Gelim->G->adjncy;
  vwght = Gelim->G->vwght;
  len = Gelim->len;
  elen = Gelim->elen;
  parent = Gelim->parent;
  score = Gelim->score;

  /* -----------------------------------------------------------------
     build the new element/variable list for each variable in reachset
     ----------------------------------------------------------------- */
  for (i = 0; i < nreach; i++)
   { u = reachset[i];
     vwght[u] = -vwght[u];        /* mark all variables in reachset */
     jstart = xadj[u];
     jstop = xadj[u] + len[u];
     jdest = jfirstolde = jstart;

#ifdef DEBUG
     printf("Updating adjacency list of node %d\n", u);
#endif

     /* --------------------------------------------------------
        scan the list of elements associated with variable u
        place newly formed elements at the beginning of the list
        -------------------------------------------------------- */
     for (j = jstart; j < jstart + elen[u]; j++)
      { e = adjncy[j];

#ifdef DEBUG
        printf("  >> element %d (score %d, parent %d)\n", e,score[e],parent[e]);
#endif

        if (score[e] == -4)       /* e has been absorbed in this elim. step */
         { me = parent[e];        /* me is the newly formed element */
           if (tmp[me] < *pflag)
            { adjncy[jdest++] = adjncy[jfirstolde];  /* move 1st old e to end */
              adjncy[jfirstolde++] = me;             /* append me at the beg. */
              tmp[me] = *pflag;
            }
         }
        else                      /* e has not been absorbed, i.e. it is */
          if (tmp[e] < *pflag)    /* an old element */
           { adjncy[jdest++] = e;
             tmp[e] = *pflag;
           }
      }
     jfirstv = jdest;             /* list of variables starts here */

     /* -------------------------------------------------------
        scan the list of variables associated with variable u
        place newly formed elements at the begining of the list
        ------------------------------------------------------- */
     for (j = jstart + elen[u]; j < jstop; j++)
      { v = adjncy[j];

#ifdef DEBUG
        printf("  >> variable %d (score %d)\n", v, score[v]);
#endif

        if (score[v] == -3)       /* v has been eliminated in this step */
         { if (tmp[v] < *pflag)   /* and, thus, forms a newly created elem. */
           { adjncy[jdest++] = adjncy[jfirstv];      /* move 1st var. to end  */
             adjncy[jfirstv++] = adjncy[jfirstolde]; /* move 1st old e to end */
             adjncy[jfirstolde++] = v;               /* append v at the beg.  */
             tmp[v] = *pflag;
           }
         }
        else
          adjncy[jdest++] = v;    /* v is still a variable */
      }
     elen[u] = jfirstv - jstart;
     len[u] = jdest - jstart;
     (*pflag)++;                  /* clear tmp for next round */

#ifdef DEBUG
     printf(" node %d: neighboring elements:\n", u);
     for (j = jstart; j < jstart + elen[u]; j++)
       printf("%5d", adjncy[j]);
     printf("\n node %d: neighboring variables:\n", u);
     for (j = jstart + elen[u]; j < jstart + len[u]; j++)
       printf("%5d", adjncy[j]);
     printf("\n");
#endif
   }

  /* ---------------------------------------------------------
     remove from each list all covered edges between variables
     --------------------------------------------------------- */
  for (i = 0; i < nreach; i++)
   { u = reachset[i];
     jstart = xadj[u];
     jstop = jstart + len[u];
     marku = FALSE;

     for (jdest = j = jstart + elen[u]; j < jstop; j++)
      { v = adjncy[j];
        if (vwght[v] > 0)         /* v does not belong to reachset */
          adjncy[jdest++] = v;    /* edge (u,v) not covered */
        if (vwght[v] < 0)         /* both vertices belong to reachset */
         { covered = FALSE;       /* check for a common element */
           if (!marku)
            { for (jj = jstart; jj < jstart + elen[u]; jj++)  /* mark elem. */
                tmp[adjncy[jj]] = *pflag;                     /* of u       */
              marku = TRUE;
            }
           for (jj = xadj[v]; jj < xadj[v] + elen[v]; jj++)   /* check elem. */
             if (tmp[adjncy[jj]] == *pflag)                   /* of v        */
              { covered = TRUE;
                break;
              }
           if (!covered)
             adjncy[jdest++] = v;
         }
      }
     len[u] = jdest - jstart;
     (*pflag)++;                  /* clear tmp for next round */

#ifdef DEBUG
     printf(" node %d: neighboring uncovered variables:\n", u);
     for (j = jstart + elen[u]; j < jstart + len[u]; j++)
       printf("%5d", adjncy[j]);
     printf("\n");
#endif
   }

  /* --------------------------------
     unmark all variables in reachset
     -------------------------------- */
  for (i = 0; i < nreach; i++)
   { u = reachset[i]; 
     vwght[u] = -vwght[u];
   }
}


/*****************************************************************************
******************************************************************************/
void
findIndNodes(gelim_t *Gelim, PORD_INT *reachset, PORD_INT nreach, PORD_INT *bin, PORD_INT *next,
             PORD_INT *tmp, PORD_INT *pflag)
{ PORD_INT *xadj, *adjncy, *vwght, *len, *elen, *parent, *score;
  PORD_INT nvtx, chk, keepon, u, v, w, wlast, i, j, jstart, jstop, jstep, jj, jjstop;
  nvtx = Gelim->G->nvtx;
  xadj = Gelim->G->xadj;
  adjncy = Gelim->G->adjncy;
  vwght = Gelim->G->vwght;
  len = Gelim->len;
  elen = Gelim->elen;
  parent = Gelim->parent;
  score = Gelim->score;

#ifdef DEBUG
  printf("Checking reachset for indistinguishable variables\n");
#endif

  /* -----------------------------------------------------------------------
     compute checksums for all principal variables on reachset and fill bins
     NOTE: checksums are stored in parent vector
     ----------------------------------------------------------------------- */
  for (i = 0; i < nreach; i++)
   { u = reachset[i];
     chk = 0;
     jstart = xadj[u];
     jstop = jstart + len[u];
     /* Modified by JYL: 16 march 2005:
      * This code was failing in case of
      * overflow.
     for (j = jstart; j < jstop; j++)
         chk += adjncy[j];
     chk = chk % nvtx;
     */
     jstep=max(1000000000/nvtx,1);
     for (j = jstart; j < jstop; j+=jstep)
     {
       jjstop = min(jstop, j+jstep);
       for (jj = j; jj < jjstop; jj++)
         chk += adjncy[jj];
       chk = chk % nvtx;
     }

     parent[u] = chk;
     /* JYL: temporary:
        if (parent[u] < - 10)
        printf("Probleme %d \n",chk);*/
     next[u] = bin[chk];
     bin[chk] = u;
   }

  /* -----------------------
     supervariable detection
     ----------------------- */
  for (i = 0; i < nreach; i++)
   { u = reachset[i];
     if (vwght[u] > 0)         /* u is a principal variable */
      { chk = parent[u];       /* search bin[chk] for ind. nodes */
        v = bin[chk];          /* okay, v is the first node in this bin */
        bin[chk] = -1;         /* no further examinations of this bin */
        while (v != -1)
         { jstart = xadj[v];
           jstop = xadj[v] + len[v];
           for (j = jstart; j < jstop; j++)
             tmp[adjncy[j]] = *pflag;
           w = next[v];        /* v is principal and w is a potential */
           wlast = v;          /* nonprincipal variable */
           while (w != -1)
            { keepon = TRUE;
              if ((len[w] != len[v]) || (elen[w] != elen[v])
                 || ((score[w] < 0) && (score[v] >= 0))
                 || ((score[w] >= 0) && (score[v] < 0)))
                keepon = FALSE;
              if (keepon)
               { for (jj = xadj[w]; jj < xadj[w] + len[w]; jj++)
                   if (tmp[adjncy[jj]] < *pflag)
                    { keepon = FALSE;
                      break;
                    }
               }
              if (keepon)                /* found it! mark w as nonprincipal */
               { parent[w] = v;          /* representative of w is v */
                 /* Temporary JY
		    if (parent[w] < - 10)
	               printf("Probleme\n");
                  */
#ifdef DEBUG
                 printf(" non-principal variable %d (score %d) mapped onto "
                        "%d (score %d)\n", w, score[w], v, score[v]); 
#endif

                 vwght[v] += vwght[w];   /* add weight of w */
                 vwght[w] = 0;
                 xadj[w] = -1;           /* w's adjacency list can be over- */
                 score[w] = -2;          /* written during next crunch */
                 w = next[w];
                 next[wlast] = w;        /* remove w from bin */
               }
              else                       /* failed */
               { wlast = w;
                 w = next[w];
               }
            }
           v = next[v];        /* no more variables can be absorbed by v */
           (*pflag)++;         /* clear tmp vector for next round */
         }
      }
   }

  /* -------------------------------------------------------
     re-initialize parent vector for all principal variables
     ------------------------------------------------------- */
  for (i = 0; i < nreach; i++)
   { u = reachset[i];
     if (vwght[u] > 0)
       parent[u] = -1;
   }
}   


/*****************************************************************************
******************************************************************************/
void
updateDegree(gelim_t *Gelim, PORD_INT *reachset, PORD_INT nreach, PORD_INT *bin)
{ PORD_INT *xadj, *adjncy, *vwght, *len, *elen, *degree;
  PORD_INT totvwght, deg, vwghtv, u, v, w, e, me, r, i, istart, istop;
  PORD_INT j, jstart, jstop;

  totvwght = Gelim->G->totvwght;
  xadj = Gelim->G->xadj;
  adjncy = Gelim->G->adjncy;
  vwght = Gelim->G->vwght;
  len = Gelim->len;
  elen = Gelim->elen;
  degree = Gelim->degree;

  /* -------------------------------------------------------------------
     degree update only for those vertices in reachset that are adjacent
     to an element
     ------------------------------------------------------------------- */
  for (r = 0; r < nreach; r++)
   { u = reachset[r];
     if (elen[u] > 0)
       bin[u] = 1;
   }

  /* -----------------------------------------
     and now do the approximate degree updates
     ----------------------------------------- */
  for (r = 0; r < nreach; r++)
   { u = reachset[r];
     if (bin[u] == 1)          /* me is the most recently formed element */
      { me = adjncy[xadj[u]];  /* in the neighborhood of u */

#ifdef DEBUG
        printf("Updating degree of all variables in L(%d) (initiated by %d)\n",
               me, u);
#endif

        /* ----------------------------------------------------------------
           compute in bin[e] the size of Le\Lme for all unabsorbed elements
           ---------------------------------------------------------------- */
        istart = xadj[me];
        istop = istart + len[me];         /* compute in bin[e] the size */
        for (i = istart; i < istop; i++)  /* of Le/Lme for all elements */
         { v = adjncy[i];                 /* e != me that are adjacent  */
           vwghtv = vwght[v];             /* to a principal var. e Lme  */
           if (vwghtv > 0)
            { jstart = xadj[v];
              jstop = jstart + elen[v];
              for (j = jstart; j < jstop; j++)
               { e = adjncy[j];
                 if (e != me)
                  { if (bin[e] > 0) bin[e] -= vwghtv;
                    else bin[e] = degree[e] - vwghtv;
                  }
               }
            }
         }

#ifdef DEBUG
        for (i = istart; i < istop; i++)
         { v = adjncy[i];
           if (vwght[v] > 0)
             for (j = xadj[v]; j < xadj[v] + elen[v]; j++)
              { e = adjncy[j];
                if (e != me)
                  printf("  >> element %d: degree %d, outer degree %d\n", e,
                         degree[e], bin[e]);
              }
         }
#endif

        /* ------------------------------------------------------
           update approx. degree for all v in Lme with bin[v] = 1
           ------------------------------------------------------ */
        for (i = istart; i < istop; i++)
         { v = adjncy[i];                  /* update the upper bound deg. */
           vwghtv = vwght[v];              /* of all principal variables  */
           deg = 0;                        /* in Lme that have not been   */
           if (bin[v] == 1)                /* updated yet                 */
            { jstart = xadj[v];
              jstop = jstart + len[v];

              /* scan the element list associated with principal v */
              for (j = jstart; j < jstart + elen[v]; j++)
               { e = adjncy[j];            
                 if (e != me) deg += bin[e];
               }

              /* scan the supervariables in the list associated with v */
              for (j = jstart + elen[v]; j < jstop; j++)
               { w = adjncy[j];
                 deg += vwght[w];
               }

              /* compute the external degree of v (add size of Lme) */
              deg = min(degree[v], deg);
              degree[v] = max(1, min(deg+degree[me]-vwghtv, totvwght-vwghtv));
              bin[v] = -1;

#ifdef DEBUG
              printf("  >> variable %d (totvwght %d, vwght %d): deg %d, "
                     "degme %d, approx degree %d\n", v, totvwght, vwghtv, deg,
                     degree[me], degree[v]);
#endif
            }
         }

        /* ------------------------------------
           clear bin[e] of all elements e != me
           ------------------------------------ */
        for (i = istart; i < istop; i++)
         { v = adjncy[i];
           vwghtv = vwght[v];
           if (vwghtv > 0)
            { jstart = xadj[v];
              jstop = jstart + elen[v];
              for (j = jstart; j < jstop; j++)
               { e = adjncy[j];
                 if (e != me) bin[e] = -1;
               }
            }
         }
      }
   }
}


/*****************************************************************************
******************************************************************************/
void
updateScore(gelim_t *Gelim, PORD_INT *reachset, PORD_INT nreach, PORD_INT scoretype, PORD_INT *bin)
{ PORD_INT *xadj, *adjncy, *vwght, *len, *elen, *degree, *score;
  PORD_INT vwghtv, deg, degme, u, v, me, r, i, istart, istop;
  /* Modified by JYL, 16 march 2005.
   * scr could overflow for quasi dense rows.
   * Use a double instead for large degrees
   * aset it near to MAX_INT in case of problem.
   */
  double scr_dbl;
  PORD_INT scr;

  xadj = Gelim->G->xadj;
  adjncy = Gelim->G->adjncy;
  vwght = Gelim->G->vwght;
  len = Gelim->len;
  elen = Gelim->elen;
  degree = Gelim->degree;
  score = Gelim->score;

  /* ------------------------------------------------------------------
     score update only for those vertices in reachset that are adjacent
     to an element
     ------------------------------------------------------------------ */
  for (r = 0; r < nreach; r++)
   { u = reachset[r];
     if (elen[u] > 0)
       bin[u] = 1;
   }

  /* ----------------------------
     and now do the score updates
     ---------------------------- */
  scoretype = scoretype % 10;
  for (r = 0; r < nreach; r++)
   { u = reachset[r];
     if (bin[u] == 1)          /* me is the most recently formed element */
      { me = adjncy[xadj[u]];  /* in the neighborhood of u */

#ifdef DEBUG
        printf("Updating score of all variables in L(%d) (initiated by %d)\n",
               me, u);
#endif

        istart = xadj[me];
        istop = xadj[me] + len[me];
        for (i = istart; i < istop; i++)
         { v = adjncy[i];                  /* update score of all principal */
           if (bin[v] == 1)                /* variables in Lme that have not */
            { vwghtv = vwght[v];           /* been updated yet */
              deg = degree[v];
              degme = degree[me] - vwghtv;
	      if (deg > 40000 || degme > 40000)
	      {
              switch(scoretype)
               { case AMD:
                   scr_dbl = (double)deg;
                   break;
                 case AMF:
                   scr_dbl = (double)deg*(double)(deg-1)/2 - (double)degme*(double)(degme-1)/2;
                   break;
                 case AMMF:
                   scr_dbl = ((double)deg*(double)(deg-1)/2 - (double)degme*(double)(degme-1)/2) / (double)vwghtv;
                   break;
                 case AMIND:
                   scr_dbl = max(0, ((double)deg*(double)(deg-1)/2 - (double)degme*(double)(degme-1)/2)
                             - (double)deg*(double)vwghtv);
                   break;
                 default:
                   fprintf(stderr, "\nError in function updateScore\n"
                        "  unrecognized selection strategy %d\n", scoretype);
                   quit();
               }
              /* Some buckets have offset nvtx / 2.
	       * Using MAX_INT - nvtx should then be safe */
              score[v] = (PORD_INT) (min(scr_dbl,MAX_INT-Gelim->G->nvtx));
	      }
	      else
	      {
              switch(scoretype)
               { case AMD:
                   scr = deg;
                   break;
                 case AMF:
                   scr = deg*(deg-1)/2 - degme*(degme-1)/2;
                   break;
                 case AMMF:
                   scr = (deg*(deg-1)/2 - degme*(degme-1)/2) / vwghtv;
                   break;
                 case AMIND:
                   scr = max(0, (deg*(deg-1)/2 - degme*(degme-1)/2)
                             - deg*vwghtv);
                   break;
                 default:
                   fprintf(stderr, "\nError in function updateScore\n"
                        "  unrecognized selection strategy %d\n", scoretype);
                   quit();
                 }
               score[v] = scr;
	      }
              bin[v] = -1;

#ifdef DEBUG
              printf("  >> variable %d (me %d): weight %d, (ext)degme %d, "
                     "degree %d, score %d\n", u, me, vwghtv, degme, degree[v],
                     score[v]);
#endif
         
              if (score[v] < 0)
               { fprintf(stderr, "\nError in function updateScore\n"
                      " score[%d] = %d is negative\n", v, score[v]);
                 quit();
               }
            }
         }
      }
   }
}


/*****************************************************************************)
******************************************************************************/
elimtree_t*
extractElimTree(gelim_t *Gelim)
{ elimtree_t *T;
  PORD_INT        *vwght, *par, *degree, *score, *sib, *fch;
  PORD_INT        *ncolfactor, *ncolupdate, *parent, *vtx2front;
  PORD_INT        nvtx, nfronts, root, u, v, front;

  nvtx = Gelim->G->nvtx;
  vwght = Gelim->G->vwght;
  par = Gelim->parent;
  degree = Gelim->degree;
  score = Gelim->score;

  /* ------------------------
     allocate working storage
     ------------------------ */
  mymalloc(sib, nvtx, PORD_INT);
  mymalloc(fch, nvtx, PORD_INT);
  for (u = 0; u < nvtx; u++)
    sib[u] = fch[u] = -1;

  /* --------------------------------------------------------------
     count fronts and create top-down view of the tree given by par
     -------------------------------------------------------------- */
  nfronts = 0;
  root = -1;
  for (u = 0; u < nvtx; u++)
    switch(score[u])
     { case -2:          /* variable u is nonprincipal */
         break;
       case -3:          /* variable u has been elim. and now forms an elem. */
         sib[u] = root;
         root = u;
         nfronts++;
         break;
       case -4:          /* element u has been absorbed by par[u] */
         v = par[u];
         sib[u] = fch[v];
         fch[v] = u;
         nfronts++;
         break;
       default:
         fprintf(stderr, "\nError in function extractElimTree\n"
              "  ordering not complete (score[%d] = %d)\n", u, score[u]);
         quit();
     }

#ifdef DEBUG
  for (u = 0; u < nvtx; u++)
    printf("node %d: score %d, par %d, fch %d, sib %d\n", u, score[u],
           par[u], fch[u], sib[u]);
#endif

  /* --------------------------------------
     allocate space for the elimtree object
     -------------------------------------- */
  T = newElimTree(nvtx, nfronts);
  ncolfactor = T->ncolfactor;
  ncolupdate = T->ncolupdate;
  parent = T->parent;
  vtx2front = T->vtx2front;

  /* -------------------------------------------------------------
     fill the vtx2front vector so that representative vertices are
     mapped in a post-order traversal
     ------------------------------------------------------------- */
  nfronts = 0;
  u = root;
  while (u != -1)
   { while (fch[u] != -1)
       u = fch[u];
     vtx2front[u] = nfronts++;
     while ((sib[u] == -1) && (par[u] != -1))
      { u = par[u];
        vtx2front[u] = nfronts++;
      }
     u = sib[u];
   }

  /* ---------------------------------------------------
     fill in the vtx2front map for nonprincipal vertices
     --------------------------------------------------- */
  for (u = 0; u < nvtx; u++)
    if (score[u] == -2)
     { v = u;
       while ((par[v] != -1) && (score[v] == -2))
         v = par[v];
       vtx2front[u] = vtx2front[v];
     }

  /* -------------------------------------------------------------
     set up the parent vector of T and fill ncolfactor, ncolupdate
     ------------------------------------------------------------- */
  for (u = 0; u < nvtx; u++)
   { front = vtx2front[u];
     if (score[u] == -3)
      { parent[front] = -1;
        ncolfactor[front] = vwght[u];
        ncolupdate[front] = degree[u];
      }
     if (score[u] == -4)
      { parent[front] = vtx2front[par[u]];
        ncolfactor[front] = vwght[u];
        ncolupdate[front] = degree[u];
      }
   }

  /* ----------------------------
     set up all other arrays of T
     ---------------------------- */
  initFchSilbRoot(T);

    /* ----------------------
     free memory and return
     ---------------------- */
  free(sib); free(fch);
  return(T);
}

