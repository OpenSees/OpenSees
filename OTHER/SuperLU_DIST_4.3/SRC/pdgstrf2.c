

/*! @file 
 * \brief Performs panel LU factorization.
 *
 * <pre>
 * -- Distributed SuperLU routine (version 4.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * August 15, 2014
 *
 * <pre>
 * Purpose
 * =======
 *   Panel factorization -- block column k
 *
 *   Factor diagonal and subdiagonal blocks and test for exact singularity.
 *   Only the column processes that own block column *k* participate
 *   in the work.
 *
 * Arguments
 * =========
 * options (input) superlu_options_t* (global)
 *         The structure defines the input parameters to control
 *         how the LU decomposition will be performed.
 *
 * k0     (input) int (global)
 *        Counter of the next supernode to be factorized.
 *
 * k      (input) int (global)
 *        The column number of the block column to be factorized.
 *
 * thresh (input) double (global)
 *        The threshold value = s_eps * anorm.
 *
 * Glu_persist (input) Glu_persist_t*
 *        Global data structures (xsup, supno) replicated on all processes.
 *
 * grid   (input) gridinfo_t*
 *        The 2D process mesh.
 *
 * Llu    (input/output) LocalLU_t*
 *        Local data structures to store distributed L and U matrices.
 *
 * U_diag_blk_send_req (input/output) MPI_Request*
 *        List of send requests to send down the diagonal block of U.
 *
 * tag_ub (input) int
 *        Upper bound of MPI tag values.
 *
 * stat   (output) SuperLUStat_t*
 *        Record the statistics about the factorization.
 *        See SuperLUStat_t structure defined in util.h.
 *
 * info   (output) int*
 *        = 0: successful exit
 *        < 0: if info = -i, the i-th argument had an illegal value
 *        > 0: if info = i, U(i,i) is exactly zero. The factorization has
 *             been completed, but the factor U is exactly singular,
 *             and division by zero will occur if it is used to solve a
 *             system of equations.
 * </pre>
 */

#include <math.h>
#include "superlu_ddefs.h"

/* This pdgstrf2 is based on TRSM function */
void
pdgstrf2_trsm
    (superlu_options_t * options, int_t k0, int_t k, double thresh,
     Glu_persist_t * Glu_persist, gridinfo_t * grid, LocalLU_t * Llu,
     MPI_Request * U_diag_blk_send_req, int tag_ub,
     SuperLUStat_t * stat, int *info)
{
    /* printf("entering pdgstrf2 %d \n", grid->iam); */
    int cols_left, iam, l, pkk, pr;
    int incx = 1, incy = 1;

    int nsupr;                  /* number of rows in the block (LDA) */
    int nsupc;                /* number of columns in the block */
    int luptr;
    int_t i, myrow, krow, j, jfst, jlst, u_diag_cnt;
    int_t *xsup = Glu_persist->xsup;
    double *lusup, temp;
    double *ujrow, *ublk_ptr;   /* pointer to the U block */
    double alpha = -1, zero = 0.0;
    int_t Pr;
    MPI_Status status;
    MPI_Comm comm = (grid->cscp).comm;

    /* Initialization. */
    iam = grid->iam;
    Pr = grid->nprow;
    myrow = MYROW (iam, grid);
    krow = PROW (k, grid);
    pkk = PNUM (PROW (k, grid), PCOL (k, grid), grid);
    j = LBj (k, grid);          /* Local block number */
    jfst = FstBlockC (k);
    jlst = FstBlockC (k + 1);
    lusup = Llu->Lnzval_bc_ptr[j];
    nsupc = SuperSize (k);
    if (Llu->Lrowind_bc_ptr[j])
        nsupr = Llu->Lrowind_bc_ptr[j][1];
    else
        nsupr = 0;
#ifdef PI_DEBUG
    printf ("rank %d  Iter %d  k=%d \t dtrsm nsuper %d \n",
            iam, k0, k, nsupr);
#endif
    ublk_ptr = ujrow = Llu->ujrow;

    luptr = 0;                  /* Point to the diagonal entries. */
    cols_left = nsupc;          /* supernode size */
    int ld_ujrow = nsupc;       /* leading dimension of ujrow */
    u_diag_cnt = 0;
    incy = ld_ujrow;

    if ( U_diag_blk_send_req && 
	 U_diag_blk_send_req[myrow] != MPI_REQUEST_NULL ) {
        /* There are pending sends - wait for all Isend to complete */
        for (pr = 0; pr < Pr; ++pr)
            if (pr != myrow) {
                MPI_Wait (U_diag_blk_send_req + pr, &status);
            }

	/* flag no more outstanding send request. */
	U_diag_blk_send_req[myrow] = MPI_REQUEST_NULL;
    }

    if (iam == pkk) {            /* diagonal process */
        for (j = 0; j < jlst - jfst; ++j) {  /* for each column in panel */
            /* Diagonal pivot */
            i = luptr;
            /* Not to replace zero pivot.  */
            if (options->ReplaceTinyPivot == YES && lusup[i] != 0.0 )  {
                if (fabs (lusup[i]) < thresh) {  /* Diagonal */

#if ( PRNTlevel>=2 )
                    printf ("(%d) .. col %d, tiny pivot %e  ",
                            iam, jfst + j, lusup[i]);
#endif
                    /* Keep the new diagonal entry with the same sign. */
                    if (lusup[i] < 0)  lusup[i] = -thresh;
                    else  lusup[i] = thresh;
#if ( PRNTlevel>=2 )
                    printf ("replaced by %e\n", lusup[i]);
#endif
                    ++(stat->TinyPivots);
                }
            }

#if 0
            for (l = 0; l < cols_left; ++l, i += nsupr, ++u_diag_cnt)
                 ublk_ptr[u_diag_cnt] = lusup[i]; /* copy one row of U */
#endif

            /* storing U in full form  */
            int st;
            for (l = 0; l < cols_left; ++l, i += nsupr, ++u_diag_cnt) {
                st = j * ld_ujrow + j;
                ublk_ptr[st + l * ld_ujrow] = lusup[i]; /* copy one row of U */
            }

            if ( ujrow[0] == zero ) { /* Test for singularity. */
                *info = j + jfst + 1;
            } else {              /* Scale the j-th column within diag. block. */
                temp = 1.0 / ujrow[0];
                for (i = luptr + 1; i < luptr - j + nsupc; ++i)
		    lusup[i] *= temp;
                stat->ops[FACT] += nsupc - j - 1;
            }

            /* Rank-1 update of the trailing submatrix within diag. block. */
            if (--cols_left) {
                /* l = nsupr - j - 1;  */
                l = nsupc - j - 1;  /* Piyush */
                dger_ (&l, &cols_left, &alpha, &lusup[luptr + 1], &incx,
                       &ujrow[ld_ujrow], &incy, &lusup[luptr + nsupr + 1],
                       &nsupr);
                stat->ops[FACT] += 2 * l * cols_left;
            }

            /* ujrow = ublk_ptr + u_diag_cnt;  */
            ujrow = ujrow + ld_ujrow + 1; /* move to next row of U */
            luptr += nsupr + 1; /* move to next column */

        }                       /* for column j ...  first loop */

	/* ++++++++++second step ====== */

        ublk_ptr = ujrow = Llu->ujrow;

        if (U_diag_blk_send_req && iam == pkk)  { /* Send the U block */
            /** ALWAYS SEND TO ALL OTHERS - TO FIX **/
            for (pr = 0; pr < Pr; ++pr)
                if (pr != krow) {
                    /* tag = ((k0<<2)+2) % tag_ub;        */
                    /* tag = (4*(nsupers+k0)+2) % tag_ub; */
                    MPI_Isend (ublk_ptr, nsupc * nsupc, MPI_DOUBLE, pr,
                               SLU_MPI_TAG (4, k0) /* tag */ ,
                               comm, U_diag_blk_send_req + pr);

                }

	    /* flag outstanding Isend */
            U_diag_blk_send_req[krow] = (MPI_Request) TRUE; /* Sherry */
        }

        /* pragma below would be changed by an MKL call */

        char uplo = 'u', side = 'r', transa = 'n', diag = 'n';

        l = nsupr - nsupc;
        // n = nsupc;
        double alpha = 1.0;
#ifdef PI_DEBUG
        printf ("calling dtrsm\n");
        printf ("dtrsm diagonal param 11:  %d \n", nsupr);
#endif

#if defined (USE_VENDOR_BLAS)
        dtrsm_ (&side, &uplo, &transa, &diag,
                &l, &nsupc,
                &alpha, ublk_ptr, &ld_ujrow, &lusup[nsupc], &nsupr,
		1, 1, 1, 1);
#else
        dtrsm_ (&side, &uplo, &transa, &diag,
                &l, &nsupc,
                &alpha, ublk_ptr, &ld_ujrow, &lusup[nsupc], &nsupr);
#endif

    } else {  /* non-diagonal process */
        /* ================================================ *
         * Receive the diagonal block of U                  *
         * for panel factorization of L(:,k)                *
         * note: we block for panel factorization of L(:,k) *
         * but panel factorization of U(:,k) don't          *
         * ================================================ */

        /* tag = ((k0<<2)+2) % tag_ub;        */
        /* tag = (4*(nsupers+k0)+2) % tag_ub; */
        // printf("hello message receiving%d %d\n",(nsupc*(nsupc+1))>>1,SLU_MPI_TAG(4,k0));
        MPI_Recv (ublk_ptr, (nsupc * nsupc), MPI_DOUBLE, krow,
                  SLU_MPI_TAG (4, k0) /* tag */ ,
                  comm, &status);
        if (nsupr > 0) {
            char uplo = 'u', side = 'r', transa = 'n', diag = 'n';
            double alpha = 1.0;

#ifdef PI_DEBUG
            printf ("dtrsm non diagonal param 11:  %d \n", nsupr);
            if (!lusup)
                printf (" Rank :%d \t Empty block column occured :\n", iam);
#endif
#if defined (USE_VENDOR_BLAS)
            dtrsm_ (&side, &uplo, &transa, &diag,
                    &nsupr, &nsupc,
                    &alpha, ublk_ptr, &ld_ujrow, lusup, &nsupr, 1, 1, 1, 1);
#else
            dtrsm_ (&side, &uplo, &transa, &diag,
                    &nsupr, &nsupc,
                    &alpha, ublk_ptr, &ld_ujrow, lusup, &nsupr);
#endif
        }

    }                           /* end if pkk ... */

    /* printf("exiting pdgstrf2 %d \n", grid->iam);  */

}  /* PDGSTRF2_trsm */


/************************************************************************/
void pdgstrs2_omp
/************************************************************************/
(int_t k0, int_t k, Glu_persist_t * Glu_persist,
 gridinfo_t * grid, LocalLU_t * Llu, SuperLUStat_t * stat)
{
#ifdef PI_DEBUG
    printf("====Entering pdgstrs2==== \n");
#endif
    int iam, pkk;
    int incx = 1;
    int nsupr;                /* number of rows in the block L(:,k) (LDA) */
    int segsize;
    int nsupc;                /* number of columns in the block */
    int_t luptr, iukp, rukp;
    int_t b, gb, j, klst, knsupc, lk, nb;
    int_t *xsup = Glu_persist->xsup;
    int_t *usub;
    double *lusup, *uval;

#ifdef _OPENMP
    int thread_id = omp_get_thread_num ();
    int num_thread = omp_get_num_threads ();
#else
    int thread_id = 0;
    int num_thread = 1;
#endif

    /* Quick return. */
    lk = LBi (k, grid);         /* Local block number */
    if (!Llu->Unzval_br_ptr[lk]) return;

    /* Initialization. */
    iam = grid->iam;
    pkk = PNUM (PROW (k, grid), PCOL (k, grid), grid);
    int k_row_cycle = k / grid->nprow;  // for know which cycle k exist; (to assign thread wise blocking) 
    int gb_col_cycle;
    klst = FstBlockC (k + 1);
    knsupc = SuperSize (k);
    usub = Llu->Ufstnz_br_ptr[lk];  /* index[] of block row U(k,:) */
    uval = Llu->Unzval_br_ptr[lk];
    nb = usub[0];
    iukp = BR_HEADER;
    rukp = 0;
    if (iam == pkk) {
        lk = LBj (k, grid);
        nsupr = Llu->Lrowind_bc_ptr[lk][1]; /* LDA of lusup[] */
        lusup = Llu->Lnzval_bc_ptr[lk];
    } else {
        nsupr = Llu->Lsub_buf_2[k0 % (1 + stat->num_look_aheads)][1];   /* LDA of lusup[] */
        lusup = Llu->Lval_buf_2[k0 % (1 + stat->num_look_aheads)];
    }

    /* Loop through all the row blocks. */
    for (b = 0; b < nb; ++b)  {
        /* assuming column cyclic distribution of data among threads */
        gb = usub[iukp];
        gb_col_cycle = gb / grid->npcol;
        nsupc = SuperSize (gb);
        iukp += UB_DESCRIPTOR;

        /* Loop through all the segments in the block. */
        for (j = 0; j < nsupc; ++j) {
#ifdef PI_DEBUG
            printf("segsize %d klst %d usub[%d] : %d",segsize,klst ,iukp,usub[iukp]);
#endif 
            segsize = klst - usub[iukp++];
            if (segsize) {    /* Nonzero segment. */
                luptr = (knsupc - segsize) * (nsupr + 1);

		/* if gb belongs to present thread then do the factorize */
                if ((gb_col_cycle + k_row_cycle + 1) % num_thread == thread_id) {
#ifdef PI_DEBUG
                    printf ("dtrsv param 4 %d param 6 %d\n", segsize, nsupr);
#endif
#if defined (USE_VENDOR_BLAS)
                    dtrsv_ ("L", "N", "U", &segsize, &lusup[luptr], &nsupr,
                            &uval[rukp], &incx, 1, 1, 1);
#else
                    dtrsv_ ("L", "N", "U", &segsize, &lusup[luptr], &nsupr,
                            &uval[rukp], &incx);
#endif
                }

                if (thread_id == 0)
                    stat->ops[FACT] += segsize * (segsize + 1); // master thread updated the stats
                rukp += segsize;
            }
        }
    }                           /* for b ... */

} /* PDGSTRS2_omp */

