/*! @file
 * \brief Performs static scheduling for the look-ahead factorization algorithm.
 *
 * <pre>
 * -- Distributed SuperLU routine (version 4.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * August 15, 2014
 *
 * Modified:
 *
 * Reference:
 * 
 * </pre>
 */

#include "superlu_ddefs.h"

#define ISORT                   /* Note: qsort() has bug on Mac */

#ifdef ISORT
extern void isort (int_t N, int_t * ARRAY1, int_t * ARRAY2);
extern void isort1 (int_t N, int_t * ARRAY);

#else

int
superlu_sort_perm (const void *arg1, const void *arg2)
{
    const int_t *val1 = (const int_t *) arg1;
    const int_t *val2 = (const int_t *) arg2;
    return (*val2 < *val1);
}
#endif

int
static_schedule(superlu_options_t * options, int m, int n, 
		LUstruct_t * LUstruct, gridinfo_t * grid, SuperLUStat_t * stat,
		int_t *perm_c_supno, int_t *iperm_c_supno, int *info)
{
    int_t *xsup;
    int_t  i, ib, jb, lb, nlb, il, iu;
    int_t Pc, Pr;
    int iam, krow, yourcol, mycol, myrow; 
    int j, k, nsupers;  /* k - current panel to work on */
    int_t *index;
    Glu_persist_t *Glu_persist = LUstruct->Glu_persist;
    LocalLU_t *Llu = LUstruct->Llu;
    int ncb, nrb, p, pr, pc, nblocks;
    int_t *etree_supno_l, *etree_supno, *blocks, *blockr, *Ublock, *Urows,
        *Lblock, *Lrows, *sf_block, *sf_block_l, *nnodes_l,
        *nnodes_u, *edag_supno_l, *recvbuf, **edag_supno;
    float edag_supno_l_bytes;
    int nnodes, *sendcnts, *sdispls, *recvcnts, *rdispls, *srows, *rrows;
    etree_node *head, *tail, *ptr;
    int *num_child;

    int iword = sizeof (int_t);

    /* Test the input parameters. */
    *info = 0;
    if (m < 0) *info = -2;
    else if (n < 0) *info = -3;
    if (*info) {
        pxerr_dist ("static_schedule", grid, -*info);
        return (-1);
    }

    /* Quick return if possible. */
    if (m == 0 || n == 0) return 0;
 
    /* 
     * Initialization.  
     */
    iam = grid->iam;
    Pc = grid->npcol; 
    Pr = grid->nprow;
    myrow = MYROW (iam, grid);
    mycol = MYCOL (iam, grid);
    nsupers = Glu_persist->supno[n - 1] + 1;
    xsup = Glu_persist->xsup;
    nblocks = 0;
    ncb = nsupers / Pc;
    nrb = nsupers / Pr;

#if ( DEBUGlevel >= 1 ) 
    print_memorylog(stat, "before static schedule");
#endif

    /* ================================================== *
     * static scheduling of j-th step of LU-factorization *
     * ================================================== */
    if (options->lookahead_etree == YES &&  /* use e-tree of symmetrized matrix and */
        (options->ParSymbFact == NO ||  /* 1) symmetric fact with serial symbolic, or */
         (options->SymPattern == YES && /* 2) symmetric pattern, and                  */
          options->RowPerm == NOROWPERM))) { /* no rowperm to destroy symmetry */

        /* if symmetric pattern or using e-tree of |A^T|+|A|,
           then we can use a simple tree structure for static schduling */

        if (options->ParSymbFact == NO) {
            /* Use the etree computed from serial symb. fact., and turn it
               into supernodal tree.  */
            int_t *etree = LUstruct->etree;
#if ( PRNTlevel>=1 )
            if (grid->iam == 0)
                printf (" === using column e-tree ===\n");
#endif

            /* look for the first off-diagonal blocks */
            etree_supno = SUPERLU_MALLOC (nsupers * sizeof (int_t));
	    log_memory(nsupers * iword, stat);

            for (i = 0; i < nsupers; i++) etree_supno[i] = nsupers;

            for (j = 0, lb = 0; lb < nsupers; lb++) {
                for (k = 0; k < SuperSize (lb); k++) {
                    jb = Glu_persist->supno[etree[j + k]];
                    if (jb != lb)
                        etree_supno[lb] = SUPERLU_MIN (etree_supno[lb], jb);
                }
                j += SuperSize (lb);
            }
        } else { /* ParSymbFACT==YES and SymPattern==YES and RowPerm == NOROWPERM */
            /* Compute an "etree" based on struct(L),
               assuming struct(U) = struct(L').   */
#if ( PRNTlevel>=1 )
            if (grid->iam == 0)
                printf (" === using supernodal e-tree ===\n");
#endif

            /* find the first block in each supernodal-column of local L-factor */
            etree_supno_l = SUPERLU_MALLOC (nsupers * sizeof (int_t));
	    log_memory(nsupers * iword, stat);

            for (i = 0; i < nsupers; i++) etree_supno_l[i] = nsupers;
            for (lb = 0; lb < ncb; lb++) {
                jb = lb * grid->npcol + mycol;
                index = Llu->Lrowind_bc_ptr[lb];
                if (index) {   /* Not an empty column */
                    i = index[0];
                    k = BC_HEADER;
                    krow = PROW (jb, grid);
                    if (krow == myrow) {  /* skip the diagonal block */
                        k += LB_DESCRIPTOR + index[k + 1];
                        i--;
                    }
                    if (i > 0)
                    {
                        etree_supno_l[jb] = index[k];
                        k += LB_DESCRIPTOR + index[k + 1];
                        i--;
                    }

                    for (j = 0; j < i; j++)
                    {
                        etree_supno_l[jb] =
                            SUPERLU_MIN (etree_supno_l[jb], index[k]);
                        k += LB_DESCRIPTOR + index[k + 1];
                    }
                }
            }
            if (mycol < nsupers % grid->npcol) {
                jb = ncb * grid->npcol + mycol;
                index = Llu->Lrowind_bc_ptr[ncb];
                if (index) {     /* Not an empty column */
                    i = index[0];
                    k = BC_HEADER;
                    krow = PROW (jb, grid);
                    if (krow == myrow) { /* skip the diagonal block */
                        k += LB_DESCRIPTOR + index[k + 1];
                        i--;
                    }
                    if (i > 0) {
                        etree_supno_l[jb] = index[k];
                        k += LB_DESCRIPTOR + index[k + 1];
                        i--;
                    }
                    for (j = 0; j < i; j++) {
                        etree_supno_l[jb] =
                            SUPERLU_MIN (etree_supno_l[jb], index[k]);
                        k += LB_DESCRIPTOR + index[k + 1];
                    }
                }
            }

            /* form global e-tree */
            etree_supno = SUPERLU_MALLOC (nsupers * sizeof (int_t));

            MPI_Allreduce (etree_supno_l, etree_supno, nsupers, mpi_int_t,
                           MPI_MIN, grid->comm);

            SUPERLU_FREE (etree_supno_l);
        }

        /* initialize number of children for each node */
        num_child = SUPERLU_MALLOC (nsupers * sizeof (int_t));
        for (i = 0; i < nsupers; i++) num_child[i] = 0;
        for (i = 0; i < nsupers; i++)
            if (etree_supno[i] != nsupers)  num_child[etree_supno[i]]++;

        /* push initial leaves to the fifo queue */
        nnodes = 0;
        for (i = 0; i < nsupers; i++) {
            if (num_child[i] == 0) {
                ptr = SUPERLU_MALLOC (sizeof (etree_node));
                ptr->id = i;
                ptr->next = NULL;
                /*printf( " == push leaf %d (%d) ==\n",i,nnodes ); */
                nnodes++;

                if (nnodes == 1) {
                    head = ptr;
                    tail = ptr;
                } else {
                    tail->next = ptr;
                    tail = ptr;
                }
            }
        }

        /* process fifo queue, and compute the ordering */
        i = 0;

        while (nnodes > 0) {
            ptr = head;
            j = ptr->id;
            head = ptr->next;
            perm_c_supno[i] = j;
            SUPERLU_FREE (ptr);
            i++;
            nnodes--;

            if (etree_supno[j] != nsupers) {
                num_child[etree_supno[j]]--;
                if (num_child[etree_supno[j]] == 0) {
                    nnodes++;

                    ptr = SUPERLU_MALLOC (sizeof (etree_node));
                    ptr->id = etree_supno[j];
                    ptr->next = NULL;

                    /*printf( "=== push %d ===\n",ptr->id ); */
                    if (nnodes == 1) {
                        head = ptr;
                        tail = ptr;
                    } else {
                        tail->next = ptr;
                        tail = ptr;
                    }
                }
            }
            /*printf( "\n" ); */
        }
        SUPERLU_FREE (num_child);
        SUPERLU_FREE (etree_supno);
	log_memory(-2 * nsupers * iword, stat);

    } else {         /* Unsymmetric pattern */

        /* Need to process both L- and U-factors, use the symmetrically
           pruned graph of L & U instead of tree (very naive implementation) */
        int nrbp1 = nrb + 1;
	float Ublock_bytes, Urows_bytes, Lblock_bytes, Lrows_bytes;

        /* allocate some workspace */
        if (! (sendcnts = SUPERLU_MALLOC ((4 + 2 * nrbp1) * Pr * Pc * sizeof (int))))
            ABORT ("Malloc fails for sendcnts[].");
	log_memory((4 + 2 * nrbp1) * Pr * Pc * sizeof (int), stat);

        sdispls = &sendcnts[Pr * Pc];
        recvcnts = &sdispls[Pr * Pc];
        rdispls = &recvcnts[Pr * Pc];
        srows = &rdispls[Pr * Pc];
        rrows = &srows[Pr * Pc * nrbp1];

        myrow = MYROW (iam, grid);
#if ( PRNTlevel>=1 )
        if (grid->iam == 0)
            printf (" === using DAG ===\n");
#endif

        /* send supno block of local U-factor to a processor *
         * who owns the corresponding block of L-factor      */

        /* srows   : # of block to send to a processor from each supno row */
        /* sendcnts: total # of blocks to send to a processor              */
        for (p = 0; p < Pr * Pc * nrbp1; p++) srows[p] = 0;
        for (p = 0; p < Pr * Pc; p++) sendcnts[p] = 0;

        /* sending blocks of U-factors corresponding to L-factors */
        /* count the number of blocks to send */
        for (lb = 0; lb < nrb; ++lb) {
            jb = lb * Pr + myrow;
            pc = jb % Pc;
            index = Llu->Ufstnz_br_ptr[lb];

            if (index) {         /* Not an empty row */
                k = BR_HEADER;
                nblocks += index[0];
                for (j = 0; j < index[0]; ++j) {
                    ib = index[k];
                    pr = ib % Pr;
                    p = pr * Pc + pc;
                    sendcnts[p]++;
                    srows[p * nrbp1 + lb]++;

                    k += UB_DESCRIPTOR + SuperSize (index[k]);
                }
            }
        }

        if (myrow < nsupers % grid->nprow) {
            jb = nrb * Pr + myrow;
            pc = jb % Pc;
            index = Llu->Ufstnz_br_ptr[nrb];

            if (index) {         /* Not an empty row */
                k = BR_HEADER;
                nblocks += index[0];
                for (j = 0; j < index[0]; ++j) {
                    ib = index[k];
                    pr = ib % Pr;
                    p = pr * Pc + pc;
                    sendcnts[p]++;
                    srows[p * nrbp1 + nrb]++;
                    k += UB_DESCRIPTOR + SuperSize (index[k]);
                }
            }
        }

        /* insert blocks to send */
        sdispls[0] = 0;
        for (p = 1; p < Pr * Pc; p++) sdispls[p] = sdispls[p - 1] + sendcnts[p - 1];
        if (!(blocks = intMalloc_dist (nblocks)))
            ABORT ("Malloc fails for blocks[].");
	log_memory( nblocks * iword, stat );

        for (lb = 0; lb < nrb; ++lb) {
            jb = lb * Pr + myrow;
            pc = jb % Pc;
            index = Llu->Ufstnz_br_ptr[lb];

            if (index) {       /* Not an empty row */
                k = BR_HEADER;
                for (j = 0; j < index[0]; ++j) {
                    ib = index[k];
                    pr = ib % Pr;
                    p = pr * Pc + pc;
                    blocks[sdispls[p]] = ib;
                    sdispls[p]++;

                    k += UB_DESCRIPTOR + SuperSize (index[k]);
                }
            }
        }

        if (myrow < nsupers % grid->nprow) {
            jb = nrb * Pr + myrow;
            pc = jb % Pc;
            index = Llu->Ufstnz_br_ptr[nrb];

            if (index) {       /* Not an empty row */
                k = BR_HEADER;
                for (j = 0; j < index[0]; ++j) {
                    ib = index[k];
                    pr = ib % Pr;
                    p = pr * Pc + pc;
                    blocks[sdispls[p]] = ib;
                    sdispls[p]++;

                    k += UB_DESCRIPTOR + SuperSize (index[k]);
                }
            }
        }

        /* communication */
        MPI_Alltoall (sendcnts, 1, MPI_INT, recvcnts, 1, MPI_INT, grid->comm);
        MPI_Alltoall (srows, nrbp1, MPI_INT, rrows, nrbp1, MPI_INT, grid->comm);

	log_memory( -(nblocks * iword), stat );  /* blocks[] to be freed soon */

        nblocks = recvcnts[0];
        rdispls[0] = sdispls[0] = 0;
        for (p = 1; p < Pr * Pc; p++) {
            rdispls[p] = rdispls[p - 1] + recvcnts[p - 1];
            sdispls[p] = sdispls[p - 1] + sendcnts[p - 1];
            nblocks += recvcnts[p];
        }

        if (!(blockr = intMalloc_dist (nblocks))) ABORT ("Malloc fails for blockr[].");
	log_memory( nblocks * iword, stat );

        MPI_Alltoallv (blocks, sendcnts, sdispls, mpi_int_t, blockr, recvcnts,
                       rdispls, mpi_int_t, grid->comm);

        SUPERLU_FREE (blocks); /* memory logged before */

	
        /* store the received U-blocks by rows */
        nlb = nsupers / Pc;
        if (!(Ublock = intMalloc_dist (nblocks))) ABORT ("Malloc fails for Ublock[].");
        if (!(Urows = intMalloc_dist (1 + nlb))) ABORT ("Malloc fails for Urows[].");

	Ublock_bytes = nblocks * iword;
	Urows_bytes = (1 + nlb) * iword;
	log_memory( Ublock_bytes + Urows_bytes, stat );

        k = 0;
        for (jb = 0; jb < nlb; jb++) {
            j = jb * Pc + mycol;
            pr = j % Pr;
            lb = j / Pr;
            Urows[jb] = 0;

            for (pc = 0; pc < Pc; pc++) {
                p = pr * Pc + pc; /* the processor owning this block of U-factor */

                for (i = rdispls[p]; i < rdispls[p] + rrows[p * nrbp1 + lb];
                     i++) {
                    Ublock[k] = blockr[i];
                    k++;
                    Urows[jb]++;
                }
                rdispls[p] += rrows[p * nrbp1 + lb];
            }
            /* sort by the column indices to make things easier for later on */

#ifdef ISORT
            isort1 (Urows[jb], &(Ublock[k - Urows[jb]]));
#else
            qsort (&(Ublock[k - Urows[jb]]), (size_t) (Urows[jb]),
                   sizeof (int_t), &superlu_sort_perm);
#endif
        }
        if (mycol < nsupers % grid->npcol) {
            j = nlb * Pc + mycol;
            pr = j % Pr;
            lb = j / Pr;
            Urows[nlb] = 0;

            for (pc = 0; pc < Pc; pc++) {
                p = pr * Pc + pc;
                for (i = rdispls[p]; i < rdispls[p] + rrows[p * nrbp1 + lb];
                     i++) {
                    Ublock[k] = blockr[i];
                    k++;
                    Urows[nlb]++;
                }
                rdispls[p] += rrows[p * nrb + lb];
            }
#ifdef ISORT
            isort1 (Urows[nlb], &(Ublock[k - Urows[nlb]]));
#else
            qsort (&(Ublock[k - Urows[nlb]]), (size_t) (Urows[nlb]),
                   sizeof (int_t), &superlu_sort_perm);
#endif
        }
        SUPERLU_FREE (blockr);
	log_memory( -nblocks * iword, stat );

        /* sort the block in L-factor */
        nblocks = 0;
        for (lb = 0; lb < ncb; lb++) {
            jb = lb * Pc + mycol;
            index = Llu->Lrowind_bc_ptr[lb];
            if (index) {        /* Not an empty column */
                nblocks += index[0];
            }
        }
        if (mycol < nsupers % grid->npcol) {
            jb = ncb * Pc + mycol;
            index = Llu->Lrowind_bc_ptr[ncb];
            if (index) {       /* Not an empty column */
                nblocks += index[0];
            }
        }

        if (!(Lblock = intMalloc_dist (nblocks))) ABORT ("Malloc fails for Lblock[].");
        if (!(Lrows = intMalloc_dist (1 + ncb))) ABORT ("Malloc fails for Lrows[].");

	Lblock_bytes = nblocks * iword;
	Lrows_bytes = (1 + ncb) * iword;
	log_memory(Lblock_bytes + Lrows_bytes, stat);

        for (lb = 0; lb <= ncb; lb++) Lrows[lb] = 0;
        nblocks = 0;
        for (lb = 0; lb < ncb; lb++) {
            Lrows[lb] = 0;

            jb = lb * Pc + mycol;
            index = Llu->Lrowind_bc_ptr[lb];
            if (index) {      /* Not an empty column */
                i = index[0];
                k = BC_HEADER;
                krow = PROW (jb, grid);
                if (krow == myrow)  /* skip the diagonal block */
                {
                    k += LB_DESCRIPTOR + index[k + 1];
                    i--;
                }

                for (j = 0; j < i; j++) {
                    Lblock[nblocks] = index[k];
                    Lrows[lb]++;
                    nblocks++;

                    k += LB_DESCRIPTOR + index[k + 1];
                }
            }
#ifdef ISORT
            isort1 (Lrows[lb], &(Lblock[nblocks - Lrows[lb]]));
#else
            qsort (&(Lblock[nblocks - Lrows[lb]]), (size_t) (Lrows[lb]),
                   sizeof (int_t), &superlu_sort_perm);
#endif
        }
        if (mycol < nsupers % grid->npcol) {
            Lrows[ncb] = 0;
            jb = ncb * Pc + mycol;
            index = Llu->Lrowind_bc_ptr[ncb];
            if (index) {       /* Not an empty column */
                i = index[0];
                k = BC_HEADER;
                krow = PROW (jb, grid);
                if (krow == myrow) { /* skip the diagonal block */
                    k += LB_DESCRIPTOR + index[k + 1];
                    i--;
                }
                for (j = 0; j < i; j++) {
                    Lblock[nblocks] = index[k];
                    Lrows[ncb]++;
                    nblocks++;
                    k += LB_DESCRIPTOR + index[k + 1];
                }
#ifdef ISORT
                isort1 (Lrows[ncb], &(Lblock[nblocks - Lrows[ncb]]));
#else
                qsort (&(Lblock[nblocks - Lrows[ncb]]), (size_t) (Lrows[ncb]),
                       sizeof (int_t), &superlu_sort_perm);
#endif
            }
        }

        /* look for the first local symmetric nonzero block match */
        if (!(sf_block = intMalloc_dist (nsupers))) ABORT ("Malloc fails for sf_block[].");
        if (!(sf_block_l = intMalloc_dist (nsupers))) ABORT ("Malloc fails for sf_block_l[].");

	log_memory( 2 * nsupers * iword, stat );

        for (lb = 0; lb < nsupers; lb++)
            sf_block_l[lb] = nsupers;
        i = 0;
        j = 0;
        for (jb = 0; jb < nlb; jb++) {
            if (Urows[jb] > 0) {
                ib = i + Urows[jb];
                lb = jb * Pc + mycol;
                for (k = 0; k < Lrows[jb]; k++) {
                    while (Ublock[i] < Lblock[j] && i + 1 < ib)
                        i++;

                    if (Ublock[i] == Lblock[j]) {
                        sf_block_l[lb] = Lblock[j];
                        j += (Lrows[jb] - k);
                        k = Lrows[jb];
                    } else {
                        j++;
                    }
                }
                i = ib;
            } else {
                j += Lrows[jb];
            }
        }
        if (mycol < nsupers % grid->npcol) {
            if (Urows[nlb] > 0) {
                ib = i + Urows[nlb];
                lb = nlb * Pc + mycol;
                for (k = 0; k < Lrows[nlb]; k++) {
                    while (Ublock[i] < Lblock[j] && i + 1 < ib)
                        i++;

                    if (Ublock[i] == Lblock[j])
                    {
                        sf_block_l[lb] = Lblock[j];
                        j += (Lrows[nlb] - k);
                        k = Lrows[nlb];
                    }
                    else
                    {
                        j++;
                    }
                }
                i = ib;
            } else {
                j += Lrows[nlb];
            }
        }

        /* compute the first global symmetric matchs */
        MPI_Allreduce (sf_block_l, sf_block, nsupers, mpi_int_t, MPI_MIN,
                       grid->comm);
        SUPERLU_FREE (sf_block_l);
	log_memory( -nsupers * iword, stat );

        /* count number of nodes in DAG (i.e., the number of blocks on and above the first match) */
        if (!(nnodes_l = intMalloc_dist (nsupers))) ABORT ("Malloc fails for nnodes_l[].");
        if (!(nnodes_u = intMalloc_dist (nsupers))) ABORT ("Malloc fails for nnodes_u[].");
	log_memory( 2 * nsupers * iword, stat );

        for (lb = 0; lb < nsupers; lb++)  nnodes_l[lb] = 0;
        for (lb = 0; lb < nsupers; lb++)  nnodes_u[lb] = 0;

        nblocks = 0;
        /* from U-factor */
        for (i = 0, jb = 0; jb < nlb; jb++) {
            lb = jb * Pc + mycol;
            ib = i + Urows[jb];
            while (i < ib) {
                if (Ublock[i] <= sf_block[lb]) {
                    nnodes_u[lb]++;
                    i++;
                    nblocks++;
                } else {     /* get out */
                    i = ib;
                }
            }
            i = ib;
        }
        if (mycol < nsupers % grid->npcol) {
            lb = nlb * Pc + mycol;
            ib = i + Urows[nlb];
            while (i < ib) {
                if (Ublock[i] <= sf_block[lb]) {
                    nnodes_u[lb]++;
                    i++;
                    nblocks++;
                } else {         /* get out */
                    i = ib;
                }
            }
            i = ib;
        }

        /* from L-factor */
        for (i = 0, jb = 0; jb < nlb; jb++) {
            lb = jb * Pc + mycol;
            ib = i + Lrows[jb];
            while (i < ib) {
                if (Lblock[i] < sf_block[lb]) {
                    nnodes_l[lb]++;
                    i++;
                    nblocks++;
                } else {
                    i = ib;
                }
            }
            i = ib;
        }
        if (mycol < nsupers % grid->npcol) {
            lb = nlb * Pc + mycol;
            ib = i + Lrows[nlb];
            while (i < ib) {
                if (Lblock[i] < sf_block[lb]) {
                    nnodes_l[lb]++;
                    i++;
                    nblocks++;
                } else {
                    i = ib;
                }
            }
            i = ib;
        }

#ifdef USE_ALLGATHER
        /* insert local nodes in DAG */
        if (!(edag_supno_l = intMalloc_dist (nsupers + nblocks)))
            ABORT ("Malloc fails for edag_supno_l[].");
	edag_supno_l_bytes = (nsupers + nblocks) * iword;
	log_memory(edag_supno_l_bytes, stat);

        iu = il = nblocks = 0;
        for (lb = 0; lb < nsupers; lb++) {
            j = lb / Pc;
            pc = lb % Pc;

            edag_supno_l[nblocks] = nnodes_l[lb] + nnodes_u[lb];
            nblocks++;
            if (mycol == pc) {
                /* from U-factor */
                ib = iu + Urows[j];
                for (jb = 0; jb < nnodes_u[lb]; jb++) {
                    edag_supno_l[nblocks] = Ublock[iu];
                    iu++;
                    nblocks++;
                }
                iu = ib;

                /* from L-factor */
                ib = il + Lrows[j];
                for (jb = 0; jb < nnodes_l[lb]; jb++) {
                    edag_supno_l[nblocks] = Lblock[il];
                    il++;
                    nblocks++;
                }
                il = ib;
            }
        }
        SUPERLU_FREE (nnodes_u);
	log_memory(-nsupers * iword, stat);

        /* form global DAG on each processor */
        MPI_Allgather (&nblocks, 1, MPI_INT, recvcnts, 1, MPI_INT,
                       grid->comm);
        nblocks = recvcnts[0];
        rdispls[0] = 0;
        for (lb = 1; lb < Pc * Pr; lb++) {
            rdispls[lb] = nblocks;
            nblocks += recvcnts[lb];
        }
        if (!(recvbuf = intMalloc_dist (nblocks))) ABORT ("Malloc fails for recvbuf[].");
	log_memory(nblocks * iword, stat);

        MPI_Allgatherv (edag_supno_l, recvcnts[iam], mpi_int_t,
                        recvbuf, recvcnts, rdispls, mpi_int_t, grid->comm);
        SUPERLU_FREE (edag_supno_l);
	log_memory(-edag_supno_l_bytes, stat);

        if (!(edag_supno = SUPERLU_MALLOC (nsupers * sizeof (int_t *))))
            ABORT ("Malloc fails for edag_supno[].");
	log_memory(nsupers * iword, stat);

        k = 0;
        for (lb = 0; lb < nsupers; lb++) nnodes_l[lb] = 0;
        for (p = 0; p < Pc * Pr; p++) {
            for (lb = 0; lb < nsupers; lb++) {
                nnodes_l[lb] += recvbuf[k];
                k += (1 + recvbuf[k]);
            }
        }
        for (lb = 0; lb < nsupers; lb++) {
            if (nnodes_l[lb] > 0)
                if (!(edag_supno[lb] = intMalloc_dist (nnodes_l[lb])))
                    ABORT ("Malloc fails for edag_supno[lb].");
            nnodes_l[lb] = 0;
        }
        k = 0;
        for (p = 0; p < Pc * Pr; p++) {
            for (lb = 0; lb < nsupers; lb++) {
                jb = k + recvbuf[k] + 1;
                k++;
                for (; k < jb; k++) {
                    edag_supno[lb][nnodes_l[lb]] = recvbuf[k];
                    nnodes_l[lb]++;
                }
            }
        }
        SUPERLU_FREE (recvbuf);
	log_memory(-nblocks * iword, stat);

#else   /* not USE_ALLGATHER */
        int nlsupers = nsupers / Pc;
        if (mycol < nsupers % Pc) nlsupers++;

        /* insert local nodes in DAG */
        if (!(edag_supno_l = intMalloc_dist (nlsupers + nblocks)))
            ABORT ("Malloc fails for edag_supno_l[].");
	edag_supno_l_bytes = (nlsupers + nblocks) * iword;
	log_memory(edag_supno_l_bytes, stat);

        iu = il = nblocks = 0;
        for (lb = 0; lb < nsupers; lb++) {
            j = lb / Pc;
            pc = lb % Pc;
            if (mycol == pc) {
                edag_supno_l[nblocks] = nnodes_l[lb] + nnodes_u[lb];
                nblocks++;
                /* from U-factor */
                ib = iu + Urows[j];
                for (jb = 0; jb < nnodes_u[lb]; jb++) {
                    edag_supno_l[nblocks] = Ublock[iu];
                    iu++;
                    nblocks++;
                }
                iu = ib;

                /* from L-factor */
                ib = il + Lrows[j];
                for (jb = 0; jb < nnodes_l[lb]; jb++) {
                    edag_supno_l[nblocks] = Lblock[il];
                    il++;
                    nblocks++;
                }
                il = ib;
            } else if (nnodes_l[lb] + nnodes_u[lb] != 0)
                printf (" # %d: nnodes[" IFMT "]=" IFMT "+" IFMT "\n",
			grid->iam, lb, nnodes_l[lb], nnodes_u[lb]);
        }
        SUPERLU_FREE (nnodes_u);
	log_memory(-nsupers * iword, stat);

        /* form global DAG on each processor */  
        MPI_Allgather (&nblocks, 1, MPI_INT, recvcnts, 1, MPI_INT, grid->comm);
        nblocks = recvcnts[0];
        rdispls[0] = 0;
        for (lb = 1; lb < Pc * Pr; lb++) {
            rdispls[lb] = nblocks;
            nblocks += recvcnts[lb];
        }
        if (!(recvbuf = intMalloc_dist (nblocks))) ABORT ("Malloc fails for recvbuf[].");
	log_memory(nblocks * iword, stat);

        MPI_Allgatherv (edag_supno_l, recvcnts[iam], mpi_int_t,
                        recvbuf, recvcnts, rdispls, mpi_int_t, grid->comm);

        SUPERLU_FREE (edag_supno_l);
	log_memory(-edag_supno_l_bytes, stat);

        if (!(edag_supno = SUPERLU_MALLOC (nsupers * sizeof (int_t *))))
            ABORT ("Malloc fails for edag_supno[].");
	log_memory(nsupers * sizeof(int_t *), stat);

        k = 0;
        for (lb = 0; lb < nsupers; lb++) nnodes_l[lb] = 0;
        for (p = 0; p < Pc * Pr; p++) {
            yourcol = MYCOL (p, grid);

            for (lb = 0; lb < nsupers; lb++) {
                j = lb / Pc;
                pc = lb % Pc;
                if (yourcol == pc) {
                    nnodes_l[lb] += recvbuf[k];
                    k += (1 + recvbuf[k]);
                }
            }
        }
        for (lb = 0; lb < nsupers; lb++) {
            if (nnodes_l[lb] > 0)
                if (!(edag_supno[lb] = intMalloc_dist (nnodes_l[lb])))
                    ABORT ("Malloc fails for edag_supno[lb].");
            nnodes_l[lb] = 0;
        }
        k = 0;
        for (p = 0; p < Pc * Pr; p++) {
            yourcol = MYCOL (p, grid);

            for (lb = 0; lb < nsupers; lb++) {
                j = lb / Pc;
                pc = lb % Pc;
                if (yourcol == pc)
                {
                    jb = k + recvbuf[k] + 1;
                    k++;
                    for (; k < jb; k++)
                    {
                        edag_supno[lb][nnodes_l[lb]] = recvbuf[k];
                        nnodes_l[lb]++;
                    }
                }
            }
        }
        SUPERLU_FREE (recvbuf);
	log_memory( -nblocks * iword , stat);

#endif  /* end USE_ALL_GATHER */

        /* initialize the num of child for each node */
        num_child = SUPERLU_MALLOC (nsupers * sizeof (int_t));
        for (i = 0; i < nsupers; i++) num_child[i] = 0;
        for (i = 0; i < nsupers; i++) {
            for (jb = 0; jb < nnodes_l[i]; jb++) {
                num_child[edag_supno[i][jb]]++;
            }
        }

        /* push initial leaves to the fifo queue */
        nnodes = 0;
        for (i = 0; i < nsupers; i++) {
            if (num_child[i] == 0) {
                ptr = SUPERLU_MALLOC (sizeof (etree_node));
                ptr->id = i;
                ptr->next = NULL;
                /*printf( " == push leaf %d (%d) ==\n",i,nnodes ); */
                nnodes++;

                if (nnodes == 1) {
                    head = ptr;
                    tail = ptr;
                } else {
                    tail->next = ptr;
                    tail = ptr;
                }
            }
        }

        /* process fifo queue, and compute the ordering */
        i = 0;

        while (nnodes > 0) {
            /*printf( "=== pop %d (%d) ===\n",head->id,i ); */
            ptr = head;
            j = ptr->id;
            head = ptr->next;

            perm_c_supno[i] = j;
            SUPERLU_FREE (ptr);
            i++;
            nnodes--;

            for (jb = 0; jb < nnodes_l[j]; jb++) {
                num_child[edag_supno[j][jb]]--;
                if (num_child[edag_supno[j][jb]] == 0) {
                    nnodes++;

                    ptr = SUPERLU_MALLOC (sizeof (etree_node));
                    ptr->id = edag_supno[j][jb];
                    ptr->next = NULL;

                    /*printf( "=== push %d ===\n",ptr->id ); */
                    if (nnodes == 1) {
                        head = ptr;
                        tail = ptr;
                    } else {
                        tail->next = ptr;
                        tail = ptr;
                    }
                }
            }
            /*printf( "\n" ); */
        }
        for (lb = 0; lb < nsupers; lb++)
            if (nnodes_l[lb] > 0)  SUPERLU_FREE (edag_supno[lb]);

        SUPERLU_FREE (num_child);
        SUPERLU_FREE (edag_supno);
        SUPERLU_FREE (nnodes_l);
        SUPERLU_FREE (sf_block);
        SUPERLU_FREE (sendcnts);

	log_memory(-(4 * nsupers + (4 + 2 * nrbp1)*Pr*Pc) * iword, stat);

        SUPERLU_FREE (Ublock);
        SUPERLU_FREE (Urows);
        SUPERLU_FREE (Lblock);
        SUPERLU_FREE (Lrows);
	log_memory(-(Ublock_bytes + Urows_bytes + Lblock_bytes + Lrows_bytes), stat);
    }
    /* ======================== *
     * end of static scheduling *
     * ======================== */

    for (lb = 0; lb < nsupers; lb++) iperm_c_supno[perm_c_supno[lb]] = lb;

#if ( DEBUGlevel >= 1 )
    print_memorylog(stat, "after static schedule");
#endif

    return 0;
} /* STATIC_SCHEDULE */

