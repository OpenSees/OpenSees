/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: minpriority.c
/
/ author        J"urgen Schulze, University of Paderborn
/ created       01jan15
/
/ This file contains functions dealing with the minimum priority object
/
******************************************************************************

Data type:  struct minprior
              gelim_t       *Gelim;      the elimination graph of G
              multisector_t *ms;         the multisector for G
              bucket_t      *bucket;     holds unelim. vert. of actual stage
              stageinfo_t   *stageinfo;  contains statistics for each stage
              int           *reachset;   holds boundary vert. in each step
              int           nreach;      number of vertices in reachset
              int           *auxaux;     general purpose auxiliary vector
              int           *auxbin;     special auxiliary vector
              int           *auxtmp;     special auxiliary vector
              int           flag;        flag for vector auxtmp (see below)
            struct stageinfo
              int          nstep;        # of elim. steps in each stage
              int          welim;        weight of elim. vert. in each stage
              int          nzf;          # of factor entries in each stage
              FLOAT        ops;          # of factor ops. in each stage
Comments:
  o Structure used to compute a minimum priority ordering for a graph G
    with multisector ms. The elimination process is organized in stages.
    The stages are given by the multisector (i.e. ms->stage). The vertices
    of a stage are eliminated in steps. In each elimination step a maximal
    independent set of vertices with minimum priority is eliminated
  o Structure does not own multisector object => it will not be freed
  o Three auxiliary vectors can be used by functions working on minprior
    IMPORTANT INVARIANTS for vectors auxbin, auxtmp
    auxbin[i] = -1 holds at start and at end of each function
    auxtmp[i] < flag holds at start and at end of each function
Methods in lib/minpriority.c:
- minprior = newMinPriority(int nvtx, int nstages);
    o Initial: Gelim = ms = bucket = NULL,
               nreach = 0, flag = 1;
- void freeMinPriority(minprior_t *minprior);
- minprior = setupMinPriority(multisector_t *ms);
    o allocates memory for the minprior object by calling newMinPriority and
      sets up the elimination graph by a call to setupElimGraph and the bucket
      by a call to setupBucket; finally, it initializes the vectors, i.e.
      auxbin[u] = -1, auxtmp[u] = 0 for all 0 <= u <= nvtx, and
      nstep = welim = nzf = ops = 0 for all stages
- T = orderMinPriority(minprior_t *minprior options_t *options,timings_t *cpus);
    o MASTER_FUNCTION: computes a bottom-up ordering according to the specified
      ordtype e { MINIMUM_PRIORITY, INCOMPLETE_ND, MULTISECTION,
                  TRISTAGE_MULTISECTION }
    o used options:
       OPTION_ORDTYPE, OPTION_NODE_SELECTION1, OPTION_NODE_SELECTION2
    o returned timings: (see eliminateStage)
       TIME_UPDSCORE, TIME_UPDADJNCY, TIME_FINDINODES
- void eliminateStage(minprior_t *minprior, int istage, int scoretype,
                      timings_t *cpus);
    o eliminates all principal variables u with stage[u] <= istage using
      the score function given by scoretype
    o returned timings:
       TIME_UPDSCORE, TIME_UPDADJNCY, TIME_FINDINODES
- int eliminateStep(minprior_t *minprior, int istage, int scoretype);
    o the variables u with stage[u] <= istage are eliminated in steps;
      in each step a maximal independet set of variables with minimum score
      is eliminated
    o the function returns the size of the independent set, i.e. the number
      of variables that have been eliminated in the actual step

******************************************************************************/

#include <space.h>
/* #define DEBUG */
/* #define BE_CAUTIOUS */


/*****************************************************************************
******************************************************************************/
minprior_t*
newMinPriority(PORD_INT nvtx, PORD_INT nstages)
{ minprior_t  *minprior;
  stageinfo_t *stageinfo;

  mymalloc(stageinfo, nstages, stageinfo_t);
  mymalloc(minprior, 1, minprior_t);
  minprior->Gelim = NULL;
  minprior->ms = NULL;
  minprior->bucket = NULL;
  minprior->stageinfo = stageinfo;

  mymalloc(minprior->reachset, nvtx, PORD_INT);
  mymalloc(minprior->auxaux, nvtx, PORD_INT);
  mymalloc(minprior->auxbin, nvtx, PORD_INT);
  mymalloc(minprior->auxtmp, nvtx, PORD_INT);

  minprior->nreach = 0;
  minprior->flag = 1;

  return(minprior);
}


/*****************************************************************************
******************************************************************************/
void
freeMinPriority(minprior_t *minprior)
{
  freeElimGraph(minprior->Gelim);
  freeBucket(minprior->bucket);
  free(minprior->stageinfo);
  free(minprior->reachset);
  free(minprior->auxaux);
  free(minprior->auxbin);
  free(minprior->auxtmp);
  free(minprior);
}


/*****************************************************************************
******************************************************************************/
minprior_t*
setupMinPriority(multisector_t *ms)
{ minprior_t  *minprior;
  stageinfo_t *stageinfo;
  PORD_INT         *auxbin, *auxtmp;
  PORD_INT         nvtx, nstages, istage, u;

  nvtx = ms->G->nvtx;
  nstages = ms->nstages;

  minprior = newMinPriority(nvtx, nstages);
  minprior->ms = ms;
  minprior->Gelim = setupElimGraph(ms->G);
  minprior->bucket = setupBucket(nvtx, nvtx, 0);

  auxbin = minprior->auxbin;
  auxtmp = minprior->auxtmp;
  for (u = 0; u < nvtx; u++)
   { auxbin[u] = -1;
     auxtmp[u] = 0;
   }

  for (istage = 0; istage < nstages; istage++)
   { stageinfo = minprior->stageinfo + istage;
     stageinfo->nstep = 0;
     stageinfo->welim = 0;
     stageinfo->nzf = 0;
     stageinfo->ops = 0.0;
   }

  return(minprior);
}


/*****************************************************************************
******************************************************************************/
elimtree_t*
orderMinPriority(minprior_t *minprior, options_t *options, timings_t *cpus)
{ elimtree_t *T;
  PORD_INT        nvtx, nstages, istage, scoretype, ordtype;

  nvtx = minprior->Gelim->G->nvtx;
  nstages = minprior->ms->nstages;

  ordtype = options[OPTION_ORDTYPE];
  scoretype = options[OPTION_NODE_SELECTION2];

  /* ------------------------------
     check whether nstages is valid
     ------------------------------ */
  if ((nstages < 1) || (nstages > nvtx))
   { fprintf(stderr, "\nError in function orderMinPriority\n"
          "  no valid number of stages in multisector (#stages = %d)\n",
          nstages);
     quit();
   }

  if ((nstages < 2) && (ordtype != MINIMUM_PRIORITY))
   { fprintf(stderr, "\nError in function orderMinPriority\n"
          "  not enough stages in multisector (#stages = %d)\n", nstages);
     quit();
   }

  /* --------------------------------------------------------------
     first stage: eliminate all vertices in the remaining subgraphs
     -------------------------------------------------------------- */
  scoretype = options[OPTION_NODE_SELECTION1];
  eliminateStage(minprior, 0, scoretype, cpus);

  /* -------------------------------------------------------
     other stages: eliminate all vertices in the multisector
     ------------------------------------------------------- */
  switch(ordtype)
   { case MINIMUM_PRIORITY:
       break;
     case INCOMPLETE_ND:
       for (istage = 1; istage < nstages; istage++)
         eliminateStage(minprior, istage, scoretype, cpus);
       break;
     case MULTISECTION:
       eliminateStage(minprior, nstages-1, scoretype, cpus);
       break;
     default:
       fprintf(stderr, "\nError in function orderMinPriority\n"
            "  unrecognized ordering type %d\n", ordtype);
       quit();
   }

  /* -------------------------------------------
     print statistics for the elimination stages
     ------------------------------------------- */
  if ((ordtype != MINIMUM_PRIORITY) && (options[OPTION_MSGLVL] > 1))
    for (istage = 0; istage < nstages; istage++)
      printf("%4d. stage: #steps %6d, weight %6d, nzl %8d, ops %e\n", istage,
             minprior->stageinfo[istage].nstep,
             minprior->stageinfo[istage].welim,
             minprior->stageinfo[istage].nzf,
             minprior->stageinfo[istage].ops);

  /* -----------------------------------
     extract elimination tree and return
     ----------------------------------- */
  T = extractElimTree(minprior->Gelim);
  return(T);
}
       

/*****************************************************************************
******************************************************************************/
void
eliminateStage(minprior_t *minprior, PORD_INT istage, PORD_INT scoretype, timings_t *cpus)
{ gelim_t     *Gelim;
  bucket_t    *bucket;
  stageinfo_t *stageinfo;
  PORD_INT         *stage, *reachset, *auxbin, *auxtmp, *auxaux;
  PORD_INT         *degree, *score;
  PORD_INT         *pflag, nreach, nvtx, r, u, i;

  Gelim = minprior->Gelim;
  bucket = minprior->bucket;
  stage = minprior->ms->stage;
  stageinfo = minprior->stageinfo + istage;
  reachset = minprior->reachset;
  auxaux = minprior->auxaux;
  auxbin = minprior->auxbin;
  auxtmp = minprior->auxtmp;
  pflag = &(minprior->flag);
  
  nvtx = Gelim->G->nvtx;
  degree = Gelim->degree;
  score = Gelim->score;

#ifdef DEBUG
  printf("\nSTARTING NEW ELIMINATION STAGE (nedges %d, maxedges %d)\n\n",
         Gelim->G->nedges, Gelim->maxedges);
  if (istage> 0) printElimGraph(Gelim);
  /* waitkey(); */
#endif

  /* -------------------------------------------------------------
     load reachset with all principal variables in stage <= istage
     ------------------------------------------------------------- */
  nreach = 0;
  for (u = 0; u < nvtx; u++)
    if ((score[u] == -1) && (stage[u] <= istage))
     { reachset[nreach++] = u;
       score[u] = degree[u];
       /* score[u] = degree[u]*(degree[u]-1)/2; */
     }

  /* ----------------------------------------------------------------
     do an initial update of the vertices in reachset and fill bucket
     ---------------------------------------------------------------- */
  pord_starttimer(cpus[TIME_UPDSCORE]);
  updateDegree(Gelim, reachset, nreach, auxbin);
  updateScore(Gelim, reachset, nreach, scoretype, auxbin);
  pord_stoptimer(cpus[TIME_UPDSCORE]);
  for (i = 0; i < nreach; i++)
   { u = reachset[i];
     insertBucket(bucket, score[u], u);
   }

  /* -------------------------------------
     and now start the elimination process
     ------------------------------------- */
  while (TRUE)
   { if (eliminateStep(minprior, istage, scoretype) == 0)
       break;
     nreach = minprior->nreach;

#ifdef BE_CAUTIOUS
     printf("checking arrays auxtmp and auxbin\n");
     for (u = 0; u < nvtx; u++)
       if ((auxtmp[u] >= *pflag) || (auxbin[u] != -1))
        { printf("ERROR: flag = %d, auxtmp[%d] = %d, auxbin[%d] = %d\n",
                 *pflag, u, auxtmp[u], u, auxbin[u]);
          quit();
        }
#endif

     /* ----------------------------------------------------------
        update the adjacency structure of all vertices in reachset
        ---------------------------------------------------------- */
     pord_starttimer(cpus[TIME_UPDADJNCY]);
     updateAdjncy(Gelim, reachset, nreach, auxtmp, pflag);
     pord_stoptimer(cpus[TIME_UPDADJNCY]);

     /* ----------------------------------------
        find indistinguishable nodes in reachset
        ---------------------------------------- */
     pord_starttimer(cpus[TIME_FINDINODES]);
     findIndNodes(Gelim, reachset, nreach, auxbin, auxaux, auxtmp, pflag);
     pord_stoptimer(cpus[TIME_FINDINODES]);

#ifdef BE_CAUTIOUS
     printf("checking arrays auxtmp and auxbin\n");
     for (u = 0; u < nvtx; u++)
       if ((auxtmp[u] >= *pflag) || (auxbin[u] != -1))
        { printf("ERROR: flag = %d, auxtmp[%d] = %d, auxbin[%d] = %d\n",
                 *pflag, u, auxtmp[u], u, auxbin[u]);
          quit();
        }
#endif

     /* ----------------------------------------------------------------
        clean reachset of nonprincipal nodes and nodes not in this stage
        ---------------------------------------------------------------- */
     r = 0;
     for (i = 0; i < nreach; i++)
      { u = reachset[i];
        if (score[u] >= 0)
          reachset[r++] = u;
      }
     nreach = r;

     /* ---------------------------------------------------
        update the degree/score of all vertices in reachset
        --------------------------------------------------- */
     pord_starttimer(cpus[TIME_UPDSCORE]);
     updateDegree(Gelim, reachset, nreach, auxbin);
     updateScore(Gelim, reachset, nreach, scoretype, auxbin);
     pord_stoptimer(cpus[TIME_UPDSCORE]);

     /* ----------------------------
        re-insert vertices in bucket
        ---------------------------- */
     for (i = 0; i < nreach; i++)
      { u = reachset[i];
        insertBucket(bucket, score[u], u);
      }

     stageinfo->nstep++;
   }
}


/*****************************************************************************
******************************************************************************/
PORD_INT
eliminateStep(minprior_t *minprior, PORD_INT istage, PORD_INT scoretype)
{ gelim_t     *Gelim;
  bucket_t    *bucket;
  stageinfo_t *stageinfo;
  PORD_INT         *stage, *reachset, *auxtmp;
  PORD_INT         *xadj, *adjncy, *vwght, *len, *degree, *score;
  PORD_INT         *pflag, *pnreach, nelim, minscr, vwghtu, u, v, i, istart, istop;
  FLOAT       tri, rec;

  Gelim = minprior->Gelim;
  bucket = minprior->bucket;
  stage = minprior->ms->stage;
  stageinfo = minprior->stageinfo + istage;
  reachset = minprior->reachset;
  pnreach = &(minprior->nreach);
  auxtmp = minprior->auxtmp;
  pflag = &(minprior->flag);

  xadj = Gelim->G->xadj;
  adjncy = Gelim->G->adjncy;
  vwght = Gelim->G->vwght;
  len = Gelim->len;
  degree = Gelim->degree;
  score = Gelim->score;

#ifdef DEBUG
  printf("\nStarting new elimination step (nedges %d, maxedges %d)\n",
         Gelim->G->nedges, Gelim->maxedges);
  /* waitkey(); */
#endif

  /* ----------------------
     check for empty bucket
     ---------------------- */
  if ((u = minBucket(bucket)) == -1)
    return(0);
  minscr = score[u];

  /* ----------------------------------------
     loop while nodes of minimum score remain
     ---------------------------------------- */
  nelim = 0;
  *pnreach = 0;
  while (TRUE)
   { vwghtu = vwght[u];

     /* --------------------------------------------------
        increment welim and nelim and remove u from bucket
        -------------------------------------------------- */
     removeBucket(bucket, u);
     stageinfo->welim += vwghtu;
     nelim++;

     /* -----------------------------------------------------------------
        call buildElement to create element u and merge u's boundary with 
        the nodes in reachset; remove any vertex from bucket that belongs
        to u's boundary and to the actual stage
        ----------------------------------------------------------------- */
     buildElement(Gelim, u);
     istart = xadj[u];
     istop = istart + len[u];
     for (i = istart; i < istop; i++)
      { v = adjncy[i];                 /* v belongs to u's boundary */
        if (auxtmp[v] < *pflag)        /* v not yet in reachset */
         { auxtmp[v] = *pflag;
           if (stage[v] <= istage)     /* v belongs to actual stage */
             removeBucket(bucket, v);
           reachset[(*pnreach)++] = v;
         }
      }

#ifdef DEBUG
     printf("Node %d (weight %d, score %d) eliminated: (boundary weight %d)\n",
            u, vwghtu, minscr, degree[u]);
     for (i = istart; i < istop; i++)
       printf("%4d (degree %2d)", adjncy[i], degree[adjncy[i]]);
     printf("\n");
#endif

      /* ---------------------------------------------------------------
         increment the storage and operation counts for this elim. stage
         --------------------------------------------------------------- */
     tri = vwghtu;
     rec = degree[u];
     stageinfo->nzf += (PORD_INT)((tri * (tri+1)) / 2);
     stageinfo->nzf += (PORD_INT)(tri * rec);
     stageinfo->ops += (tri*tri*tri) / 3.0 + (tri*tri) / 2.0 - (5*tri) / 6.0;
     stageinfo->ops += (tri*tri*rec) + (rec*(rec+1)*tri);

     /* ---------------------------------------------------------------
        end this elim. step, if one of the following conditions is true
         (1) no multiple elimination
         (2) bucket empty
         (3) no further variable with minimum score
        ---------------------------------------------------------------- */
     if (scoretype / 10 == 0)
       break;
     if ((u = minBucket(bucket)) == -1)
       break;
     if (score[u] > minscr)
       break;
   }

  /* -----------------------
     clear auxtmp and return
     ----------------------- */
  (*pflag)++;
  return(nelim);
}

