

/*! @file 
 * \brief Scatter the computed blocks into LU destination.
 *
 * <pre>
 * -- Distributed SuperLU routine (version 4.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * October 1, 2014
 *
 */
#include <math.h>
#include "superlu_ddefs.h"

static void
dscatter_l_1 (int ib,
           int ljb,
           int nsupc,
           int_t iukp,
           int_t* xsup,
           int klst,
           int nbrow,
           int_t lptr,
           int temp_nbrow,
           int * usub,
           int * lsub,
           double *tempv,
           int * indirect_thread,
           int_t ** Lrowind_bc_ptr, double **Lnzval_bc_ptr,
	   gridinfo_t * grid)
{
    // TAU_STATIC_TIMER_START("SCATTER_LB");
    // printf("hello\n");
    int_t rel, i, segsize, jj;
    double *nzval;
    int_t *index = Lrowind_bc_ptr[ljb];
    int_t ldv = index[1];       /* LDA of the dest lusup. */
    int_t lptrj = BC_HEADER;
    int_t luptrj = 0;
    int_t ijb = index[lptrj];
    while (ijb != ib)
    {
        /* Search for dest block --
           blocks are not ordered! */
        luptrj += index[lptrj + 1];
        lptrj += LB_DESCRIPTOR + index[lptrj + 1];

        ijb = index[lptrj];
    }
    /*
     * Build indirect table. This is needed because the
     * indices are not sorted for the L blocks.
     */
    int_t fnz = FstBlockC (ib);
    lptrj += LB_DESCRIPTOR;
    for (i = 0; i < index[lptrj - 1]; ++i)
    {
        rel = index[lptrj + i] - fnz;
        indirect_thread[rel] = i;

    }

    nzval = Lnzval_bc_ptr[ljb] + luptrj;
    // tempv =bigV + (cum_nrow + cum_ncol*nbrow);
    for (jj = 0; jj < nsupc; ++jj)
    {
        segsize = klst - usub[iukp + jj];
        // printf("segsize %d \n",segsize);
        if (segsize) {
            /*#pragma _CRI cache_bypass nzval,tempv */
            for (i = 0; i < temp_nbrow; ++i) {
                rel = lsub[lptr + i] - fnz;
                nzval[indirect_thread[rel]] -= tempv[i];
                // printf("i (src) %d, perm (dest) %d  \n",i,indirect_thread[rel]);
#ifdef PI_DEBUG
                double zz = 0.0;
                // if(!(*(long*)&zz == *(long*)&tempv[i]) )
                printf ("(%d %d, %0.3e, %0.3e, %3e ) ", ljb,
                        nzval - Lnzval_bc_ptr[ljb] + indirect_thread[rel],
                        nzval[indirect_thread[rel]] + tempv[i],
                        nzval[indirect_thread[rel]],tempv[i]);
                //printing triplets (location??, old value, new value ) if none of them is zero
#endif
            }
            // printf("\n");
            tempv += nbrow;
#ifdef PI_DEBUG
            // printf("\n");
#endif
        }
        nzval += ldv;
        // printf("%d\n",nzval );
    }
    // TAU_STATIC_TIMER_STOP("SCATTER_LB");
} /* dscatter_l_1 */

static void
dscatter_l (int ib,
           int ljb,
           int nsupc,
           int_t iukp,
           int_t* xsup,
           int klst,
           int nbrow,
           int_t lptr,
           int temp_nbrow,
           int_t* usub,
           int_t* lsub,
           double *tempv,
           int* indirect_thread,int* indirect2,
           int_t ** Lrowind_bc_ptr, double **Lnzval_bc_ptr,
           gridinfo_t * grid)
{
    
    int_t rel, i, segsize, jj;
    double *nzval;
    int_t *index = Lrowind_bc_ptr[ljb];
    int_t ldv = index[1];       /* LDA of the dest lusup. */
    int_t lptrj = BC_HEADER;
    int_t luptrj = 0;
    int_t ijb = index[lptrj];
    
    while (ijb != ib)
    {
        luptrj += index[lptrj + 1];
        lptrj += LB_DESCRIPTOR + index[lptrj + 1];
        ijb = index[lptrj];
    }
    
    /*
     * Build indirect table. This is needed because the
     * indices are not sorted for the L blocks.
     */
    int_t fnz = FstBlockC (ib);
    int_t dest_nbrow; 
    lptrj += LB_DESCRIPTOR;
    dest_nbrow=index[lptrj - 1];
    
    for (i = 0; i < dest_nbrow; ++i)
    {
        rel = index[lptrj + i] - fnz;
        indirect_thread[rel] = i;

    }

    /* can be precalculated */
    for (i = 0; i < temp_nbrow; ++i)
    {
        rel = lsub[lptr + i] - fnz;
        indirect2[i] =indirect_thread[rel]; 
    }


    nzval = Lnzval_bc_ptr[ljb] + luptrj;
    for (jj = 0; jj < nsupc; ++jj)
    {
        segsize = klst - usub[iukp + jj];
        if (segsize)
        {
            for (i = 0; i < temp_nbrow; ++i)
            {
                nzval[indirect2[i]] -= tempv[i];
            }
            tempv += nbrow;
        }
        nzval += ldv;
    }
    
} /* dscatter_l */


static void
dscatter_u (int ib,
           int jb,
           int nsupc,
           int_t iukp,
           int_t * xsup,
           int klst,
           int nbrow,
           int_t lptr,
           int temp_nbrow,
           int_t* lsub,
           int_t* usub,
           double* tempv,
           int_t ** Ufstnz_br_ptr, double **Unzval_br_ptr,
           gridinfo_t * grid)
{
#ifdef PI_DEBUG
    printf ("A(%d,%d) goes to U block \n", ib, jb);
#endif
    // TAU_STATIC_TIMER_START("SCATTER_U");
    // TAU_STATIC_TIMER_START("SCATTER_UB");

    int_t jj, i, fnz, rel;
    int segsize;
    double *ucol;
    int_t ilst = FstBlockC (ib + 1);
    int_t lib = LBi (ib, grid);
    int_t *index = Ufstnz_br_ptr[lib];

    /* reinitialize the pointer to each row of U */
    int_t iuip_lib, ruip_lib;
    iuip_lib = BR_HEADER;
    ruip_lib = 0;

    int_t ijb = index[iuip_lib];
    while (ijb < jb)            /* Search for dest block. */
    {
        ruip_lib += index[iuip_lib + 1];
        // printf("supersize[%ld] \t:%ld \n",ijb,SuperSize( ijb ) );
        iuip_lib += UB_DESCRIPTOR + SuperSize (ijb);
        ijb = index[iuip_lib];
    }
    /* Skip descriptor.  Now point to fstnz index of
       block U(i,j). */
    iuip_lib += UB_DESCRIPTOR;

    // tempv = bigV + (cum_nrow + cum_ncol*nbrow);
    for (jj = 0; jj < nsupc; ++jj)
    {
        segsize = klst - usub[iukp + jj];
        fnz = index[iuip_lib++];
        if (segsize)            /* Nonzero segment in U(k.j). */
        {
            ucol = &Unzval_br_ptr[lib][ruip_lib];

            // printf("========Entering loop=========\n");
            for (i = 0; i < temp_nbrow; ++i)
            {

                rel = lsub[lptr + i] - fnz;
                // printf("%d %d %d %d %d \n",lptr,i,fnz,temp_nbrow,nbrow );
                // printf("hello   ucol[%d] %d %d : \n",rel,lsub[lptr + i],fnz);

                ucol[rel] -= tempv[i];

                // printf("hello\n");

#ifdef PI_DEBUG
                double zz = 0.0;
                if (!(*(long *) &zz == *(long *) &tempv[i]))
                    printf ("(%d, %0.3e, %0.3e ) ", rel, ucol[rel] + tempv[i],
                            ucol[rel]);
                //printing triplets (location??, old value, new value ) if none of them is zero
#endif
            }                   /* for i=0..temp_nbropw */
            tempv += nbrow;
#ifdef PI_DEBUG
            // printf("\n");
#endif
        }                       /*ig segsize */
        ruip_lib += ilst - fnz;

    }                           /*for jj=0:nsupc */
#ifdef PI_DEBUG
    // printf("\n");
#endif
    // TAU_STATIC_TIMER_STOP("SCATTER_UB");
} /* dscatter_u */


static void
arrive_at_ublock (int_t j,      //block number 
                  int_t * iukp, // output 
                  int_t * rukp, int_t * jb, /* Global block number of block U(k,j). */
                  int_t * ljb,  /* Local block number of U(k,j). */
                  int_t * nsupc,    /*supernode size of destination block */
                  int_t iukp0,  //input
                  int_t rukp0, int_t * usub,    /*usub scripts */
                  int_t * perm_u,   /*permutation matrix */
                  int_t * xsup, /*for SuperSize and LBj */
                  gridinfo_t * grid)
{
    int_t jj;
    *iukp = iukp0;
    *rukp = rukp0;

#ifdef ISORT
    for (jj = 0; jj < perm_u[j]; jj++)
#else
    for (jj = 0; jj < perm_u[2 * j + 1]; jj++)
#endif
    {
        /* reinitilize the pointers to the begining of the */
        /* kth column/row of L/U factors                   */
        // printf("iukp %d \n",*iukp );
        *jb = usub[*iukp];      /* Global block number of block U(k,j). */
        // printf("jb %d \n",*jb );
        *nsupc = SuperSize (*jb);
        // printf("nsupc %d \n",*nsupc );
        *iukp += UB_DESCRIPTOR; /* Start fstnz of block U(k,j). */

        *rukp += usub[*iukp - 1];   /* Move to block U(k,j+1) */
        *iukp += *nsupc;
    }

    /* reinitilize the pointers to the begining of the */
    /* kth column/row of L/U factors                   */
    *jb = usub[*iukp];          /* Global block number of block U(k,j). */
    *ljb = LBj (*jb, grid);     /* Local block number of U(k,j). */
    *nsupc = SuperSize (*jb);
    *iukp += UB_DESCRIPTOR;     /* Start fstnz of block U(k,j). */
}

typedef struct pair pair;

int_t
get_min (int_t * sums, int_t nprocs)
{
    int_t min_ind, min_val;
    min_ind = 0;
    min_val = 2147483647;
    for (int i = 0; i < nprocs; i++)
    {
        if (sums[i] < min_val)
        {
            min_val = sums[i];
            min_ind = i;
        }
    }

    return min_ind;
}

int
compare_pair (const void *a, const void *b)
{
    return (((pair *) a)->val - ((pair *) b)->val);
}

int_t
static_partition (pair * work_load, int_t nwl, int_t * partition, int_t ldp,
                  int_t * sums, int_t * counts, int nprocs)
{
    //initialization loop
    for (int i = 0; i < nprocs; ++i)
    {
        counts[i] = 0;
        sums[i] = 0;
    }
    qsort (work_load, nwl, sizeof (pair), compare_pair);
    // for(int i=0;i<nwl;i++)
    for (int i = nwl - 1; i >= 0; i--)
    {
        int_t ind = get_min (sums, nprocs);
        // printf("ind %d\n",ind );
        partition[ldp * ind + counts[ind]] = work_load[i].ind;
        counts[ind]++;
        sums[ind] += work_load[i].val;

    }

    return 0;
}

/*Divide CPU-GPU dgemm work here*/
#ifdef PI_DEBUG
int Ngem = 2;
// int_t Ngem = 0;
int min_gpu_col = 6;
#else

    // int_t Ngem = 0;

#endif


#ifdef GPU_ACC

void
gemm_division_cpu_gpu(
    int* num_streams_used,  /*number of streams that will be used */
    int* stream_end_col,    /*array holding last column blk for each partition */
    int * ncpu_blks,        /*Number of CPU dgemm blks */
    /*input */
    int nbrow,              /*number of row in A matrix */
    int ldu,                /*number of k in dgemm */
    int nstreams, 
    int* full_u_cols,       /*array containing prefix sum of work load */
    int num_blks            /*Number of work load */
)
{
    int Ngem = sp_ienv(7);  /*get_mnk_dgemm ();*/
    int min_gpu_col = get_cublas_nb ();

    // Ngem = 1000000000;
    /*
       cpu is to gpu dgemm should be ideally 0:1 ratios to hide the total cost
       However since there is gpu latency of around 20,000 ns implying about
       200000 floating point calculation be done in that time so ~200,000/(2*nbrow*ldu)
       should be done in cpu to hide the latency; we Ngem =200,000/2 
     */
    int i, j;

    // {
    //     *num_streams_used=0;
    //     *ncpu_blks = num_blks;
    //     return;
    // }

    for (i = 0; i < nstreams; ++i)
    {
        stream_end_col[i] = num_blks;
    }

    *ncpu_blks = 0;
    /*easy returns -1 when number of column are less than threshold */
    if (full_u_cols[num_blks - 1] < (Ngem / (nbrow * ldu)) || num_blks == 1 )
    {
        *num_streams_used = 0;
        *ncpu_blks = num_blks;
#ifdef PI_DEBUG
        printf ("full_u_cols[num_blks-1] %d  %d \n",
                full_u_cols[num_blks - 1], (Ngem / (nbrow * ldu)));
        printf ("Early return \n");
#endif
        return;

    }

    /* Easy return -2 when number of streams =0 */
    if (nstreams == 0)
    {
        *num_streams_used = 0;
        *ncpu_blks = num_blks;
        return;
        /* code */
    }
    /*find first block where count > Ngem */


    for (i = 0; i < num_blks - 1; ++i)  /*I can use binary search here */
    {
        if (full_u_cols[i + 1] > Ngem / (nbrow * ldu))
            break;
    }
    *ncpu_blks = i + 1;

    int_t cols_remain =
        full_u_cols[num_blks - 1] - full_u_cols[*ncpu_blks - 1];

#ifdef PI_DEBUG
    printf ("Remaining cols %d num_blks %d cpu_blks %d \n", cols_remain,
            num_blks, *ncpu_blks);
#endif
    if (cols_remain > 0)
    {
        *num_streams_used = 1;  /* now atleast one stream would be used */

#ifdef PI_DEBUG
        printf ("%d %d  %d %d \n", full_u_cols[num_blks - 1],
                full_u_cols[*ncpu_blks], *ncpu_blks, nstreams);
#endif
        int_t FP_MIN = 200000 / (nbrow * ldu);
        int_t cols_per_stream = SUPERLU_MAX (min_gpu_col, cols_remain / nstreams);
        cols_per_stream = SUPERLU_MAX (cols_per_stream, FP_MIN);
#ifdef PI_DEBUG
        printf ("cols_per_stream :\t%d\n", cols_per_stream);
#endif

        int_t cutoff = cols_per_stream + full_u_cols[*ncpu_blks - 1];
        for (int_t i = 0; i < nstreams; ++i)
        {
            stream_end_col[i] = num_blks;
        }
        j = *ncpu_blks;
        for (i = 0; i < nstreams - 1; ++i)
        {
            int_t st = (i == 0) ? (*ncpu_blks) : stream_end_col[i - 1];

            for (j = st; j < num_blks - 1; ++j)
            {
#ifdef PI_DEBUG
                printf ("i %d, j %d, %d  %d ", i, j, full_u_cols[j + 1],
                        cutoff);
#endif
                if (full_u_cols[j + 1] > cutoff)
                {
#ifdef PI_DEBUG
                    printf ("cutoff met \n");
#endif
                    cutoff = cols_per_stream + full_u_cols[j];
                    stream_end_col[i] = j + 1;
                    *num_streams_used += 1;
                    j++;
                    break;
                }
#ifdef PI_DEBUG
                printf ("\n");
#endif
            }

        }

    }
}

void
gemm_division_new (int * num_streams_used,   /*number of streams that will be used */
                   int * stream_end_col, /*array holding last column blk for each partition */
                   int * ncpu_blks,  /*Number of CPU dgemm blks */
                        /*input */
                   int nbrow,    /*number of row in A matrix */
                   int ldu,  /*number of k in dgemm */
                   int nstreams,
                   Ublock_info_t *Ublock_info,    /*array containing prefix sum of work load */
                   int num_blks  /*Number of work load */
    )
{
    int Ngem = sp_ienv(7); /*get_mnk_dgemm ();*/
    int min_gpu_col = get_cublas_nb ();

    Ngem = 1000000000;
    /*
       cpu is to gpu dgemm should be ideally 0:1 ratios to hide the total cost
       However since there is gpu latency of around 20,000 ns implying about
       200000 floating point calculation be done in that time so ~200,000/(2*nbrow*ldu)
       should be done in cpu to hide the latency; we Ngem =200,000/2 
     */
    int_t i, j;


    for (int i = 0; i < nstreams; ++i)
    {
        stream_end_col[i] = num_blks;
    }

    *ncpu_blks = 0;
    /*easy returns -1 when number of column are less than threshold */
    if (Ublock_info[num_blks - 1].full_u_cols < (Ngem / (nbrow * ldu)) || num_blks == 1)
    {
        *num_streams_used = 0;
        *ncpu_blks = num_blks;

        return;

    }

    /* Easy return -2 when number of streams =0 */
    if (nstreams == 0)
    {
        *num_streams_used = 0;
        *ncpu_blks = num_blks;
        return;
        /* code */
    }
    /*find first block where count > Ngem */


    for (i = 0; i < num_blks - 1; ++i)  /*I can use binary search here */
    {
        if (Ublock_info[i + 1].full_u_cols > Ngem / (nbrow * ldu))
            break;
    }
    *ncpu_blks = i + 1;

    int_t cols_remain =
       Ublock_info [num_blks - 1].full_u_cols - Ublock_info[*ncpu_blks - 1].full_u_cols;

    if (cols_remain > 0)
    {
        *num_streams_used = 1;  /* now atleast one stream would be used */

        int_t FP_MIN = 200000 / (nbrow * ldu);
        int_t cols_per_stream = SUPERLU_MAX (min_gpu_col, cols_remain / nstreams);
        cols_per_stream = SUPERLU_MAX (cols_per_stream, FP_MIN);

        int_t cutoff = cols_per_stream + Ublock_info[*ncpu_blks - 1].full_u_cols;
        for (int_t i = 0; i < nstreams; ++i)
        {
            stream_end_col[i] = num_blks;
        }
        j = *ncpu_blks;
        for (i = 0; i < nstreams - 1; ++i)
        {
            int_t st = (i == 0) ? (*ncpu_blks) : stream_end_col[i - 1];

            for (j = st; j < num_blks - 1; ++j)
            {
                if (Ublock_info[j + 1].full_u_cols > cutoff)
                {

                    cutoff = cols_per_stream + Ublock_info[j].full_u_cols;
                    stream_end_col[i] = j + 1;
                    *num_streams_used += 1;
                    j++;
                    break;
                }

            }

        }

    }
}

#endif  /* defined GPU_ACC */
