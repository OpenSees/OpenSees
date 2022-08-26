/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: interface.c
/
/ author        J"urgen Schulze, University of Paderborn
/ created       01jan26
/
/ This file contains some high level interface functions (only these
/ functions should be called by a user).
/
******************************************************************************/


#include <space.h>

/*****************************************************************************
    o Input:
        undirected graph G
        options -- if NULL, default options are used
           option[0] holds OPTION_ORDTYPE
           option[1] holds OPTION_NODE_SELECTION1
           option[2] holds OPTION_NODE_SELECTION2
           option[3] holds OPTION_NODE_SELECTION3
           option[4] holds OPTION_DOMAIN_SIZE
           option[5] holds OPTION_MSGLVL
    o Output:
        elimination/front tree T reflecting the ordering of G
        cpus -- if NULL, no timing information is pulled back
           cpus[0]  holds TIME_COMPRESS
           cpus[1]  holds TIME_MS
           cpus[2]  holds TIME_MULTILEVEL
           cpus[3]  holds TIME_INITDOMDEC
           cpus[4]  holds TIME_COARSEDOMDEC
           cpus[5]  holds TIME_INITSEP
           cpus[6]  holds TIME_REFINESEP
           cpus[7]  holds TIME_SMOOTH
           cpus[8]  holds TIME_BOTTOMUP
           cpus[9]  holds TIME_UPDADJNCY
           cpus[10] holds TIME_FINDINODES
           cpus[11] holds TIME_UPDSCORE
    o Comments:
        this function computes an ordering for G; it returns an elimination
        tree T; permutation vectors perm, invp can be extracted from T by
        calling function permFromElimTree(T, perm, invp)
******************************************************************************/
elimtree_t*
SPACE_ordering(graph_t *G, options_t *options, timings_t *cpus)
{ graph_t       *Gc;
  multisector_t *ms;
  minprior_t    *minprior;
  elimtree_t    *T, *T2;
  timings_t	cpusOrd[ORD_TIME_SLOTS];
  options_t     default_options[] = { SPACE_ORDTYPE, SPACE_NODE_SELECTION1,
                       SPACE_NODE_SELECTION2, SPACE_NODE_SELECTION3,
                       SPACE_DOMAIN_SIZE, SPACE_MSGLVL };
  PORD_INT           *vtxmap, istage, totnstep, totnzf;
  FLOAT         totops;

  /* --------------------------------------------------
     set default options, if no other options specified
     -------------------------------------------------- */
  if (options == NULL)
    options = default_options;

  /* ----------------
     reset all timers
     ---------------- */
  pord_resettimer(cpusOrd[TIME_COMPRESS]);
  pord_resettimer(cpusOrd[TIME_MS]);
  pord_resettimer(cpusOrd[TIME_MULTILEVEL]);
  pord_resettimer(cpusOrd[TIME_INITDOMDEC]);
  pord_resettimer(cpusOrd[TIME_COARSEDOMDEC]);
  pord_resettimer(cpusOrd[TIME_INITSEP]);
  pord_resettimer(cpusOrd[TIME_REFINESEP]);
  pord_resettimer(cpusOrd[TIME_SMOOTH]);
  pord_resettimer(cpusOrd[TIME_BOTTOMUP]);
  pord_resettimer(cpusOrd[TIME_UPDADJNCY]);
  pord_resettimer(cpusOrd[TIME_FINDINODES]);
  pord_resettimer(cpusOrd[TIME_UPDSCORE]);

  /* ------------------
     compress the graph
     ------------------ */
  pord_starttimer(cpusOrd[TIME_COMPRESS]);
  mymalloc(vtxmap, G->nvtx, PORD_INT);
  Gc = compressGraph(G, vtxmap);
  pord_stoptimer(cpusOrd[TIME_COMPRESS]);

  if (Gc != NULL)
   { if (options[OPTION_MSGLVL] > 0)
       printf("compressed graph constructed (#nodes %d, #edges %d)\n",
              Gc->nvtx, Gc->nedges >> 1);
   }
  else
   { Gc = G;
     free(vtxmap);
     if (options[OPTION_MSGLVL] > 0)
       printf("no compressed graph constructed\n");
   }

  /* -------------------
     compute multisector
     ------------------- */
  
  
  pord_starttimer(cpusOrd[TIME_MS]);
  ms = constructMultisector(Gc, options, cpusOrd);
  pord_stoptimer(cpusOrd[TIME_MS]);
	

  if (options[OPTION_MSGLVL] > 0)
    printf("quality of multisector: #stages %d, #nodes %d, weight %d\n",
           ms->nstages, ms->nnodes, ms->totmswght);

  /* ---------------------------------
     compute minimum priority ordering
     --------------------------------- */
  pord_starttimer(cpusOrd[TIME_BOTTOMUP])
  minprior = setupMinPriority(ms);
  T = orderMinPriority(minprior, options, cpusOrd);
  pord_stoptimer(cpusOrd[TIME_BOTTOMUP]);

  if (options[OPTION_MSGLVL] > 0)
   { totnstep = totnzf = 0;
     totops = 0.0;
     for (istage = 0; istage < ms->nstages; istage++)
      { totnstep += minprior->stageinfo[istage].nstep;
        totnzf   += minprior->stageinfo[istage].nzf;
        totops   += minprior->stageinfo[istage].ops;
      }
     printf("quality of ordering: #steps %d, nzl %d, ops %e\n", totnstep,
            totnzf, totops);
   }

  /* -----------------------
     expand elimination tree
     ----------------------- */
  if (Gc != G)
   { T2 = expandElimTree(T, vtxmap, G->nvtx);
     freeElimTree(T);
     freeGraph(Gc);
     free(vtxmap);
   }
  else T2 = T;

  /* --------------------------------------------------
     pull back timing results, if vector cpus available
     -------------------------------------------------- */
  if (cpus != NULL)
   { cpus[0]  = cpusOrd[TIME_COMPRESS];
     cpus[1]  = cpusOrd[TIME_MS];
     cpus[2]  = cpusOrd[TIME_MULTILEVEL];
     cpus[3]  = cpusOrd[TIME_INITDOMDEC];
     cpus[4]  = cpusOrd[TIME_COARSEDOMDEC];
     cpus[5]  = cpusOrd[TIME_INITSEP];
     cpus[6]  = cpusOrd[TIME_REFINESEP];
     cpus[7]  = cpusOrd[TIME_SMOOTH];
     cpus[8]  = cpusOrd[TIME_BOTTOMUP];
     cpus[9]  = cpusOrd[TIME_UPDADJNCY];
     cpus[10] = cpusOrd[TIME_FINDINODES];
     cpus[11] = cpusOrd[TIME_UPDSCORE];
   }

  /* ----------------------
     free memory and return
     ---------------------- */
  freeMultisector(ms);
  freeMinPriority(minprior);
  return(T2);
}


#if defined(cleaned_version)
/*****************************************************************************
    o Input:
        elimination/front tree T
        max. number of zeros that is allowed to be introduced in front
    o Output:
        transformed elimination/front tree T'
    o Comments:
        the goal is to make T (obtained by orderMinPriority or
        setupElimTree) more appropiate for the multifrontal algorithm
******************************************************************************/
elimtree_t*
SPACE_transformElimTree(elimtree_t *T, PORD_INT maxzeros)
{ elimtree_t *T2, *T3;

  /* -----------------------------------------------------
     1st: determine the fundamental fronts
          this step significantly improves the cache reuse
     ----------------------------------------------------- */
  T2 = fundamentalFronts(T);

  /* -----------------------------------------------------------------
     2nd: group together small subtrees into one front
          this step reduces the number of fronts and thus the overhead
          associated with them; the expense is added storage for the
          logically zero entries and the factor operations on them
     ------------------------------------------------------------------ */
  T3 = mergeFronts(T2, maxzeros);
  freeElimTree(T2);

  /* --------------------------------------------------------------
     3rd: order the children of a front so that the working storage
          in the multifrontal algorithm is minimized
     -------------------------------------------------------------- */
  (void)justifyFronts(T3);
  return(T3);
}

/*****************************************************************************
    o Input:
        transformed elimination/front tree T, input matrix A
    o Output:
        initial factor matrix L of the permuted input matrix PAP
    o Comments: L contains nonzeros of PAP; all other entries are set to 0.0
******************************************************************************/
factorMtx_t*
SPACE_symbFac(elimtree_t *T, inputMtx_t *A)
{ factorMtx_t *L;
  frontsub_t  *frontsub;
  css_t       *css;
  inputMtx_t  *PAP;
  elimtree_t  *PTP;
  PORD_INT         *perm, neqs, nelem;

  /* ------------------------------------------------------
     extract permutation vectors from T and permute T and A
     ------------------------------------------------------ */
  neqs = A->neqs;
  mymalloc(perm, neqs, PORD_INT);
  permFromElimTree(T, perm);
  PTP = permuteElimTree(T, perm);
  PAP = permuteInputMtx(A, perm);

  /* -------------------------------------------------------------------
     create factor matrix L of PAP, i.e.
      (1) create the subscript structure of the fronts, i.e. frontsub
      (2) use frontsub to create the compressed subscript structure of L
      (3) allocate memory for L and the nonzeros of L, i.e. L->nzl
      (4) init. L with the nonzeros of PAP
     ------------------------------------------------------------------- */
  frontsub = setupFrontSubscripts(PTP, PAP);
  css = setupCSSFromFrontSubscripts(frontsub);

  nelem = css->xnzl[neqs];
  L = newFactorMtx(nelem);
  L->perm = perm;
  L->frontsub = frontsub;
  L->css = css;

  initFactorMtx(L, PAP);

  /* -----------------------------------------------------
     free permuted input matrix and return
     note: PTP and perm have been inherited by frontsub, L
     ----------------------------------------------------- */
  freeInputMtx(PAP);
  return(L);
}


/*****************************************************************************
    o Input:
        transformed elimination/front tree
        initial factor matrix L of the permuted input matrix PAP
    o Output:
        factor matrix L of the permuted input matrix PAP
        cpus -- if NULL no timing information is pulled back
           cpus[0] holds TIME_INITFRONT
           cpus[1] holds TIME_EXPAND
           cpus[2] holds TIME_KERNEL
           cpus[3] holds TIME_INITUPD
    o Comments:
        this function does the actual numerical factorization; to
        improve register and cache reuse it uses a kernel of size 3x3
******************************************************************************/
void
SPACE_numFac(factorMtx_t *L, timings_t *cpus)
{ timings_t  cpusFactor[NUMFAC_TIME_SLOTS];

  /* ----------------
     reset all timers
     ---------------- */
  pord_resettimer(cpusFactor[TIME_INITFRONT]);
  pord_resettimer(cpusFactor[TIME_EXADD]);
  pord_resettimer(cpusFactor[TIME_KERNEL]);
  pord_resettimer(cpusFactor[TIME_INITUPD]);

  /* -------------------------
     compute Cholesky factor L
     ------------------------- */
  numfac(L, cpusFactor);

 /* --------------------------------------------------
     pull back timing results, if vector cpus available
     -------------------------------------------------- */
  if (cpus != NULL)
   { cpus[0] = cpusFactor[TIME_INITFRONT];
     cpus[1] = cpusFactor[TIME_EXADD];
     cpus[2] = cpusFactor[TIME_KERNEL];
     cpus[3] = cpusFactor[TIME_INITUPD];
   }
}


/*****************************************************************************
    o Input:
        transformed elimination/front tree
        factor matrix L of the permuted input matrix PAP
        right hand side vector rhs of the original system Ax = b
    o Output:
        solution vector xvec of the original system Ax = b
    o Comments:
        this function solves the remaining triangular systems;
******************************************************************************/
void
SPACE_solveTriangular(factorMtx_t *L, FLOAT *rhs, FLOAT *xvec)
{ FLOAT *yvec;
  PORD_INT   *perm;
  PORD_INT   neqs, k;

  perm = L->perm;
  neqs = L->css->neqs;

  /* -------------------------------------------
     set up permuted right hand side vector yvec
     ------------------------------------------- */
  mymalloc(yvec, neqs, FLOAT);
  for (k = 0; k < neqs; k++)
    yvec[perm[k]] = rhs[k];

  /* -------------------------
     solve Ly = b and L^Tz = y
     ------------------------- */
  forwardSubst1x1(L, yvec);
  backwardSubst1x1(L, yvec);

  /* ---------------------------------------------------------------
     extract from yvec the solution vector of the un-permuted system
     --------------------------------------------------------------- */
  for (k = 0; k < neqs; k++)
    xvec[k] = yvec[perm[k]];
  free(yvec);
}


/*****************************************************************************
    o Input:
        sparse matrix A, right hand side vector rhs
        options -- if NULL, default options are used
           option[0] holds OPTION_ORDTYPE
           option[1] holds OPTION_NODE_SELECTION1
           option[2] holds OPTION_NODE_SELECTION2
           option[3] holds OPTION_NODE_SELECTION3
           option[4] holds OPTION_DOMAIN_SIZE
           option[5] holds OPTION_MSGLVL
           option[6] holds OPTION_ETREE_NONZ
    o Output:
        solution vector xvec of the original system Ax = b
        cpus -- if NULL, no timing information is pulled back
           cpus[0]  holds time to construct the graph
           cpus[1]  holds time to compute the ordering
           cpus[2]  holds TIME_COMPRESS
           cpus[3]  holds TIME_MS
           cpus[4]  holds TIME_MULTILEVEL
           cpus[5]  holds TIME_INITDOMDEC
           cpus[6]  holds TIME_COARSEDOMDEC
           cpus[7]  holds TIME_INITSEP
           cpus[8]  holds TIME_REFINESEP
           cpus[9]  holds TIME_SMOOTH
           cpus[10] holds TIME_BOTTOMUP
           cpus[11] holds TIME_UPDADJNCY;
           cpus[12] holds TIME_FINDINODES
           cpus[13] holds TIME_UPDSCORE
           cpus[14] holds time to transform the elimination tree
           cpus[15] holds time to compute the symbolical factorization
           cpus[16] holds time to compute the numerical factorization
           cpus[17] holds TIME_INITFRONT
           cpus[18] holds TIME_EXADD
           cpus[19] holds TIME_KERNEL
           cpus[20] holds TIME_INITUPD
           cpus[21] holds time to solve the triangular systems
    o Comments:
        this is the final topmost function that can be used as a black
        box in other algorithm; it provides a general purpose direct
        solver for large sparse positive definite systems
******************************************************************************/
void
SPACE_solve(inputMtx_t *A, FLOAT *rhs, FLOAT *xvec, options_t *options,
            timings_t *cpus)
{ graph_t     *G;
  elimtree_t  *T, *T2;
  factorMtx_t *L;
  timings_t   cpusOrd[ORD_TIME_SLOTS], cpusFactor[NUMFAC_TIME_SLOTS];
  timings_t   t_graph, t_ord, t_etree, t_symb, t_num, t_solvetri;
  options_t   default_options[] = { SPACE_ORDTYPE, SPACE_NODE_SELECTION1,
                     SPACE_NODE_SELECTION2, SPACE_NODE_SELECTION3,
                     SPACE_DOMAIN_SIZE, SPACE_MSGLVL, SPACE_ETREE_NONZ };

  /* --------------------------------------------------
     set default options, if no other options specified
     -------------------------------------------------- */
  if (options == NULL)
    options = default_options;

  /* ----------------
     reset all timers
     ---------------- */
  pord_resettimer(t_graph);
  pord_resettimer(t_ord);
  pord_resettimer(t_etree);
  pord_resettimer(t_symb);
  pord_resettimer(t_num);
  pord_resettimer(t_solvetri);

  /* -----------------
     set up graph G(A)
     ----------------- */
  pord_starttimer(t_graph);
  G = setupGraphFromMtx(A);
  pord_stoptimer(t_graph);

  if (options[OPTION_MSGLVL] > 0)
    printf("\ninduced graph constructed: #vertices %d, #edges %d, #components "
           "%d\n", G->nvtx, G->nedges >> 1, connectedComponents(G));

  /* --------------------------------------------
     construct ordering/elimination tree for G(A)
     -------------------------------------------- */
  pord_starttimer(t_ord);
  T = SPACE_ordering(G, options, cpusOrd);
  pord_stoptimer(t_ord);
  freeGraph(G);

  if (options[OPTION_MSGLVL] > 0)
    printf("quality of initial elim. tree: #fronts %d, #indices %d\n\t"
           "nzl %d, ops %e, wspace %d\n", T->nfronts, nFactorIndices(T),
           nFactorEntries(T), nFactorOps(T), nWorkspace(T));

  /* -------------------------------
     elimination tree transformation
     ------------------------------- */
  pord_starttimer(t_etree);
  T2 = SPACE_transformElimTree(T, options[OPTION_ETREE_NONZ]);
  pord_stoptimer(t_etree);
  freeElimTree(T);

  if (options[OPTION_MSGLVL] > 0)
    printf("quality of transformed elim. tree: #fronts %d, #indices %d\n\t"
           "nzl %d, ops %e, wspace %d\n", T2->nfronts, nFactorIndices(T2),
           nFactorEntries(T2), nFactorOps(T2), nWorkspace(T2));

  /* ------------------------
     symbolical factorization
     ------------------------ */
  pord_starttimer(t_symb);
  L = SPACE_symbFac(T2, A);
  pord_stoptimer(t_symb);

  if (options[OPTION_MSGLVL] > 0)
    printf("quality of factor matrix:\n\tneqs %d, #indices %d, nzl %d\n",
           L->css->neqs, L->css->nind, L->nelem);

  /* -----------------------
     numerical factorization
     ----------------------- */
  pord_starttimer(t_num);
  SPACE_numFac(L, cpusFactor);
  pord_stoptimer(t_num);

  if (options[OPTION_MSGLVL] > 0)
    printf("performance of numerical factorization: %6.2f mflops\n",
            (double)nFactorOps(T2) / t_num / 1000000);

  /* ------------------------------
     solution of triangular systems
     ------------------------------ */
  pord_starttimer(t_solvetri);
  SPACE_solveTriangular(L, rhs, xvec);
  pord_stoptimer(t_solvetri);

  if (options[OPTION_MSGLVL] > 0)
    printf("performance of forward/backward solve:  %6.2f mflops\n",
            (double)nTriangularOps(T2) / t_solvetri / 1000000);

  freeElimTree(T2);
  freeFactorMtx(L);

  /* --------------------------------------------------
     pull back timing results, if vector cpus available
     -------------------------------------------------- */
  if (cpus != NULL)
   { cpus[0]  = t_graph;
     cpus[1]  = t_ord;
     cpus[2]  = cpusOrd[TIME_COMPRESS];
     cpus[3]  = cpusOrd[TIME_MS];
     cpus[4]  = cpusOrd[TIME_MULTILEVEL];
     cpus[5]  = cpusOrd[TIME_INITDOMDEC];
     cpus[6]  = cpusOrd[TIME_COARSEDOMDEC];
     cpus[7]  = cpusOrd[TIME_INITSEP];
     cpus[8]  = cpusOrd[TIME_REFINESEP];
     cpus[9]  = cpusOrd[TIME_SMOOTH];
     cpus[10] = cpusOrd[TIME_BOTTOMUP];
     cpus[11] = cpusOrd[TIME_UPDADJNCY];
     cpus[12] = cpusOrd[TIME_FINDINODES];
     cpus[13] = cpusOrd[TIME_UPDSCORE];
     cpus[14] = t_etree;
     cpus[15] = t_symb;
     cpus[16] = t_num;
     cpus[17] = cpusFactor[TIME_INITFRONT];
     cpus[18] = cpusFactor[TIME_EXADD];
     cpus[19] = cpusFactor[TIME_KERNEL];
     cpus[20] = cpusFactor[TIME_INITUPD];
     cpus[21] = t_solvetri;
   }
}


/*****************************************************************************
    o Input:
        sparse matrix A with permutation vector perm
        right hand side vector rhs
        options -- if NULL, default options are used
          option[0] holds OPTION_MSGLVL
          option[1] holds OPTION_ETREE_NONZ
    o Output:
        solution vector xvec of the original system Ax = b
        cpus -- if NULL, no timing information is pulled back
           cpus[0] holds time to construct the graph
           cpus[1] holds time to construct the elimination tree
           cpus[2] holds time to transform the elimination tree
           cpus[3] holds time to compute the symbolical factorization
           cpus[4] holds time to compute the numerical factorization
           cpus[5] holds TIME_INITFRONT
           cpus[6] holds TIME_EXADD
           cpus[7] holds TIME_KERNEL
           cpus[8] holds TIME_INITUPD
           cpus[9] holds time to solve the triangular systems
    o Comments: 
        this function can be used to solve an equation system
        using an externally computed permutation vector
******************************************************************************/
void
SPACE_solveWithPerm(inputMtx_t *A, PORD_INT *perm, FLOAT *rhs, FLOAT *xvec,
                    options_t *options, timings_t *cpus)
{ graph_t     *G;
  elimtree_t  *T, *T2;
  factorMtx_t *L;
  timings_t   cpusFactor[NUMFAC_TIME_SLOTS], t_graph, t_etree_construct;
  timings_t   t_etree_merge, t_symb, t_num, t_solvetri;
  options_t   default_options[] = { SPACE_MSGLVL, SPACE_ETREE_NONZ };
  PORD_INT         *invp, i, msglvl, maxzeros;

  /* --------------------------------------------------
     set default options, if no other options specified
     -------------------------------------------------- */
  if (options == NULL)
    options = default_options;
  msglvl = options[0];
  maxzeros = options[1];

  /* ----------------
     reset all timers
     ---------------- */
  pord_resettimer(t_graph);
  pord_resettimer(t_etree_construct);
  pord_resettimer(t_etree_merge);
  pord_resettimer(t_symb);
  pord_resettimer(t_num);
  pord_resettimer(t_solvetri);

  /* -----------------
     set up graph G(A)
     ----------------- */
  pord_starttimer(t_graph);
  G = setupGraphFromMtx(A);
  pord_stoptimer(t_graph);

  if (msglvl > 0)
    printf("\ninduced graph constructed: #vertices %d, #edges %d, #components "
           "%d\n", G->nvtx, G->nedges >> 1, connectedComponents(G));

  /* ---------------------------------------------------
     construct inital elimination tree according to perm
     --------------------------------------------------- */
  pord_starttimer(t_etree_construct);
  mymalloc(invp, G->nvtx, PORD_INT);
  for (i = 0; i < G->nvtx; i++)
    invp[perm[i]] = i;
  T = setupElimTree(G, perm, invp);
  pord_stoptimer(t_etree_construct);
  freeGraph(G);
  free(invp);

  if (msglvl > 0)
    printf("quality of initial elim. tree: #fronts %d, #indices %d\n\t"
           "nzl %d, ops %e, wspace %d\n", T->nfronts, nFactorIndices(T),
           nFactorEntries(T), nFactorOps(T), nWorkspace(T));

  /* -------------------------------
     elimination tree transformation
     ------------------------------- */
  pord_starttimer(t_etree_merge);
  T2 = SPACE_transformElimTree(T, maxzeros);
  pord_stoptimer(t_etree_merge);
  freeElimTree(T);

  if (msglvl > 0)
    printf("quality of transformed elim. tree: #fronts %d, #indices %d\n\t"
           "nzl %d, ops %e, wspace %d\n", T2->nfronts, nFactorIndices(T2),
           nFactorEntries(T2), nFactorOps(T2), nWorkspace(T2));

  /* ------------------------
     symbolical factorization
     ------------------------ */
  pord_starttimer(t_symb);
  L = SPACE_symbFac(T2, A);
  pord_stoptimer(t_symb);

  if (msglvl > 0)
    printf("quality of factor matrix:\n\tneqs %d, #indices %d, nzl %d\n",
           L->css->neqs, L->css->nind, L->nelem);

  /* -----------------------
     numerical factorization
     ----------------------- */
  pord_starttimer(t_num);
  SPACE_numFac(L, cpusFactor);
  pord_stoptimer(t_num);

  if (msglvl > 0)
    printf("performance of numerical factorization: %6.2f mflops\n",
            (double)nFactorOps(T2) / t_num / 1000000);

  /* ------------------------------
     solution of triangular systems
     ------------------------------ */
  pord_starttimer(t_solvetri);
  SPACE_solveTriangular(L, rhs, xvec);
  pord_stoptimer(t_solvetri);

  if (msglvl > 0)
    printf("performance of forward/backward solve:  %6.2f mflops\n",
            (double)nTriangularOps(T2) / t_solvetri / 1000000);

  freeElimTree(T2);
  freeFactorMtx(L);

  /* --------------------------------------------------
     pull back timing results, if vector cpus available
     -------------------------------------------------- */
  if (cpus != NULL)
   { cpus[0] = t_graph;
     cpus[1] = t_etree_construct;
     cpus[2] = t_etree_merge;
     cpus[3] = t_symb;
     cpus[4] = t_num;
     cpus[5] = cpusFactor[TIME_INITFRONT];
     cpus[6] = cpusFactor[TIME_EXADD];
     cpus[7] = cpusFactor[TIME_KERNEL];
     cpus[8] = cpusFactor[TIME_INITUPD];
     cpus[9] = t_solvetri;
   }
}


/*****************************************************************************
    o Input:
        graph G with permutation vector perm
        options -- if NULL, default options are used
          option[0] holds OPTION_MSGLVL
          option[1] holds OPTION_ETREE_NONZ
          option[2] holds OPTION_ETREE_BAL
          option[3] holds dimension of hypercube
    o Output:
        mapping object map
        cpus -- if NULL, no timing information is pulled back
           cpus[0] holds time to construct the elimination tree
           cpus[1] holds time to transform the elimination tree
           cpus[2] holds time to compute the mapping
    o Comments:
        this function can be used to obtain a mapping object for the
        parallel factorization
******************************************************************************/
mapping_t*
SPACE_mapping(graph_t *G, PORD_INT *perm, options_t *options, timings_t *cpus)
{ mapping_t *map;
  elimtree_t *T, *T2;
  timings_t  t_etree_construct, t_etree_merge, t_map;
  options_t  default_options[] = { SPACE_MSGLVL, SPACE_ETREE_NONZ,
                                   SPACE_ETREE_BAL, 2 };
  PORD_INT        *invp, i, msglvl, maxzeros, bal, dimQ;

  /* --------------------------------------------------
     set default options, if no other options specified
     -------------------------------------------------- */
  if (options == NULL)
    options = default_options;
  msglvl = options[0];
  maxzeros = options[1];
  bal = options[2];
  dimQ = options[3];

  /* ----------------
     reset all timers
     ---------------- */
  pord_resettimer(t_etree_construct);
  pord_resettimer(t_etree_merge);
  pord_resettimer(t_map);

  /* ---------------------------------------------------
     construct inital elimination tree according to perm
     --------------------------------------------------- */
  pord_starttimer(t_etree_construct);
  mymalloc(invp, G->nvtx, PORD_INT);
  for (i = 0; i < G->nvtx; i++)
    invp[perm[i]] = i;
  T = setupElimTree(G, perm, invp);
  pord_stoptimer(t_etree_construct);
  free(invp);

  if (msglvl > 0)
    printf("quality of initial elim. tree: #fronts %d, #indices %d\n\t"
           "nzl %d, ops %e, wspace %d\n", T->nfronts, nFactorIndices(T),
           nFactorEntries(T), nFactorOps(T), nWorkspace(T));

  /* -------------------------------
     elimination tree transformation
     ------------------------------- */
  pord_starttimer(t_etree_merge);
  T2 = SPACE_transformElimTree(T, maxzeros);
  pord_stoptimer(t_etree_merge);
  freeElimTree(T);

  if (msglvl > 0)
    printf("quality of transformed elim. tree: #fronts %d, #indices %d\n\t"
           "nzl %d, ops %e, wspace %d\n", T2->nfronts, nFactorIndices(T2),
           nFactorEntries(T2), nFactorOps(T2), nWorkspace(T2));

  /* -------------------
     compute the mapping
     ------------------- */
  pord_starttimer(t_map);
  map = setupMapping(T2, dimQ, bal);
  pord_stoptimer(t_map);

  /* --------------------------------------------------
     pull back timing results, if vector cpus available
     -------------------------------------------------- */
  if (cpus != NULL)
   { cpus[0] = t_etree_construct;
     cpus[1] = t_etree_merge;
     cpus[2] = t_map;
   }

  /* --------------------------------------------------------------
     return mapping object (don't free T2, since it belongs to map)
     -------------------------------------------------------------- */
  return(map);
}
#endif
