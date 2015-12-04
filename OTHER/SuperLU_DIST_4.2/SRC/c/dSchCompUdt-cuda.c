

/*! @file 
 * \brief This file contains the main loop of pdgstrf which involves
 *        rank k update of the Schur complement.
 *        Uses CUDA GPU.
 *
 * <pre>
 * -- Distributed SuperLU routine (version 4.0) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * October 1, 2014
 *
 */

#define SCHEDULE_STRATEGY dynamic

#define cublasCheckErrors(fn) \
    do { \
        cublasStatus_t __err = fn; \
        if (__err != CUBLAS_STATUS_SUCCESS) { \
            fprintf(stderr, "Fatal cublas error: %d (at %s:%d)\n", \
                (int)(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while(0);


if ( msg0 && msg2 ) {  /* L(:,k) and U(k,:) are not empty. */
    ldu   =0;
    full  =1;
    int cum_nrow;
    int temp_nbrow;

    lptr = lptr0;
    luptr = luptr0;
    
    nbrow= lsub[1];
    if (myrow==krow) nbrow = lsub[1]-lsub[3];

    if (nbrow>0) {
        
        int ncol_max = SUPERLU_MIN(buffer_size/nbrow,bigu_size/ldt);
        int num_streams_used,        /*number of streams that will be used*/
        ncpu_blks;                     /*Number of CPU dgemm blks*/

        int jjj, jjj_st,jjj_global;        
        for (j = jj0; j < nub; ++j) {
            arrive_at_ublock( j,&iukp,&rukp,&jb,&ljb,&nsupc,
	    		      iukp0,rukp0,usub,perm_u,xsup,grid );

            ncols =0 ;  //initialize at 0 
            jj = iukp;
            int temp_ldu=0; 
            for (; jj < iukp+nsupc; ++jj) {
                segsize = klst - usub[jj];
                if ( segsize ) {
		    ++ncols;
		}
                temp_ldu = SUPERLU_MAX(temp_ldu, segsize);
            }

            full_u_cols[j] = ncols;
            blk_ldu[j] = temp_ldu;
        } /* end for j = jj0..nub */

        jjj = jj0; /* initialization */
            
        // #pragma omp barrier 
        while ( jjj < nub ) {
            jjj_st=jjj;
#ifdef _OPENMP
#pragma omp single
#endif
            {
                ldu = blk_ldu[jjj_st];
                for (j = jjj_st; j < nub ; ++j) {
                    
                    /* prefix sum */
                    if (j != jjj_st) full_u_cols[j] += full_u_cols[j-1];

                    ldu = SUPERLU_MAX(ldu, blk_ldu[j]);   

                    /* break condition */
                    /* the number of columns that can be processed is limited by buffer size*/
                    if (full_u_cols[j]+((j+1==nub)?0:full_u_cols[j+1]) > ncol_max) {
                        break;
                    }
                } /* end for j=jjj_st to nub */  

                jjj_global = SUPERLU_MIN(nub, j+1); /* Maximum value of jjj will be nub */
                
                // TAU_STATIC_TIMER_START("work_divison");
                /* Divide CPU-GPU gemm here */
                gemm_division_cpu_gpu(
		       &num_streams_used, /*number of streams that will be used*/
		       stream_end_col,    /*array holding last column blk for each partition*/
		       &ncpu_blks,        /*Number of CPU gemm blks*/
		       			  /*input*/
		       nbrow,             /*number of row in A matrix*/
		       ldu,               /*number of k in dgemm*/
		       nstreams,
		       full_u_cols + jjj_st, /*array containing prefix sum of work load*/
		       jjj_global-jjj_st     /*Number of work load */
                );
                // TAU_STATIC_TIMER_STOP("work_divison");

            } /* pragma omp single */

            jjj = jjj_global;
            // printf("thread_id %d, jjj %d \n",thread_id,jjj );
            if (jjj == jjj_st+1 && full_u_cols[jjj_st] > ncol_max) {
                printf("allocate more memory for buffer !!!!\n");
                if(nbrow * full_u_cols[jjj_st] > buffer_size)
                    printf("%d buffer_size %d\n",nbrow*full_u_cols[jjj_st],buffer_size );
            }
            
            // #pragma omp barrier 
            /* gathering circuit */
            assert(jjj_st<nub);
            assert(jjj-1<nub);
            // TAU_STATIC_TIMER_START("GATHER_U");
#ifdef _OPENMP
#pragma omp for schedule( SCHEDULE_STRATEGY )
#endif
            for (j = jjj_st; j < jjj; ++j) {
                if (j==jjj_st) tempu = bigU;
                else tempu = bigU + ldu*full_u_cols[j-1];

                /* == processing each of the remaining columns == */
                arrive_at_ublock(j,&iukp,&rukp,&jb,&ljb,&nsupc,
				 iukp0,rukp0,usub,perm_u,xsup,grid);

                // tempu = tempU2d;
                for (jj = iukp; jj < iukp+nsupc; ++jj) {
                    segsize = klst - usub[jj];
                    if ( segsize ) {
                        lead_zero = ldu - segsize;
                        for (i = 0; i < lead_zero; ++i) tempu[i] = zero;
                        tempu += lead_zero;
                        for (i = 0; i < segsize; ++i)
                            tempu[i] = uval[rukp+i];
                        rukp += segsize;
                        tempu += segsize;
                    }
                }

                rukp -= usub[iukp - 1]; /* Return to start of U(k,j). */

            } /* end for j=jjj_st to jjj */  

	    if ( num_streams_used > 0 ) {
#ifdef PI_DEBUG
		printf("nbrow %d *ldu %d  =%d < ldt %d * max_row_size %d =%d \n",nbrow,ldu,nbrow*ldu,ldt,max_row_size,ldt*max_row_size );
		assert(nbrow*ldu<=ldt*max_row_size);
#endif 
		cudaMemcpy2DAsync(dA, nbrow*sizeof(double),
				  &lusup[luptr+(knsupc-ldu)*nsupr],
				  nsupr*sizeof(double), nbrow*sizeof(double),
				  ldu, cudaMemcpyHostToDevice, streams[0]);
	    }
                
	    for (int i = 0; i < num_streams_used; ++i) {
		int st = (i==0) ? ncpu_blks+jjj_st : jjj_st+stream_end_col[i-1]; 
		int st_col = full_u_cols[st-1];
		int num_col_stream = full_u_cols[jjj_st+stream_end_col[i]-1]-full_u_cols[st-1];
		tempu = bigU;
                    
		double *tempv1 = bigV + full_u_cols[st-1]*nbrow;

		/* Following is for testing purpose */
#ifdef GPU_ACC
		int stream_id = i;
		int b_offset  = ldu * st_col;
		int c_offset  = st_col * nbrow;
		size_t B_stream_size = ldu * num_col_stream * sizeof(double);
		size_t C_stream_size = nbrow * num_col_stream * sizeof(double);
		
		assert(ldu*(st_col+num_col_stream) < bigu_size);
		assert(nbrow*(st_col+num_col_stream) < buffer_size);
		
		cudaMemcpyAsync(dB+b_offset, tempu+b_offset, B_stream_size,
				cudaMemcpyHostToDevice, streams[stream_id]);
		
		cublasCheckErrors(
				  cublasSetStream(handle[stream_id],
						  streams[stream_id])
				  );
		
		cublasCheckErrors(
				  cublasDgemm(handle[stream_id],
					      CUBLAS_OP_N, CUBLAS_OP_N,
					      nbrow, num_col_stream, ldu,
                                              &alpha, dA, nbrow,
					      &dB[b_offset], ldu, 
					      &beta, &dC[c_offset],
                                              nbrow)
				  );
		
		checkCuda( cudaMemcpyAsync(tempv1, dC+c_offset,
					   C_stream_size,
					   cudaMemcpyDeviceToHost,
					   streams[stream_id]) );
#else 
		if ( num_col_stream > 0 ) {   
		    my_dgemm_("N", "N", &nbrow, &num_col_stream, &ldu,
			      &alpha, &lusup[luptr+(knsupc-ldu)*nsupr],
			      &nsupr, tempu+ldu*st_col, &ldu, &beta,
			      tempv1, &nbrow, 1, 1);
		}
		
#endif 
		
	    } /* end for i = 1 to num_streams used */
	    
	    int num_col = full_u_cols[jjj_st+ncpu_blks-1];
	    int st_col = 0;        /*special case for cpu */
	    tempv = bigV + nbrow * st_col;
	    tempu = bigU;
	    
	    double tstart = SuperLU_timer_();
#if defined (USE_VENDOR_BLAS)            
	    dgemm_("N", "N", &nbrow, &num_col, &ldu, &alpha,
		  &lusup[luptr+(knsupc-ldu)*nsupr], &nsupr,
		  tempu+ldu*st_col, &ldu, &beta, tempv, &nbrow, 1, 1);
#else
	    dgemm_("N", "N", &nbrow, &num_col, &ldu, &alpha,
		  &lusup[luptr+(knsupc-ldu)*nsupr], &nsupr,
		  tempu+ldu*st_col, &ldu, &beta, tempv, &nbrow);
#endif
	    gemm_timer += SuperLU_timer_() -tstart;
	    stat->ops[FACT] += 2 * nbrow * ldu * full_u_cols[jjj-1];
	    
	    // printf("after dgemm \n");
	    
            /* Now scattering blocks handled by cpu */
            int temp_ncol;
	    
            /* scatter first blocks which cpu has computated*/
            tstart = SuperLU_timer_();

#ifdef _OPENMP
#pragma omp parallel  \
    private(j,iukp,rukp, tempu, tempv, cum_nrow, jb, nsupc,ljb,	\
	    segsize,lead_zero,					\
	    ib, temp_nbrow,ilst,lib,index,			\
	    ijb,fnz,ucol,rel,ldv,lptrj,luptrj,			\
	    nzval,     lb ,                     jj, i)		\
    firstprivate(luptr,lptr) default (shared)
#endif
            {
                int thread_id = omp_get_thread_num();
        
                int* indirect_thread = indirect + ldt*thread_id;
                int* indirect2_thread = indirect2 + ldt*thread_id;
                double* tempv1;
                
                if (ncpu_blks< omp_get_num_threads()) {
                    // TAU_STATIC_TIMER_START("SPECIAL_CPU_SCATTER");
                    
                    for (j = jjj_st; j < jjj_st+ncpu_blks; ++j) {
                        /* code */
                        #ifdef PI_DEBUG
                            printf("scattering %d  block column\n",j);
                        #endif

                        /* == processing each of the remaining columns == */

                        if(j==jjj_st) tempv1 = bigV;
                        else tempv1 = bigV + full_u_cols[j-1]*nbrow;

                        arrive_at_ublock( j,&iukp,&rukp,&jb,&ljb,&nsupc,
					  iukp0,rukp0,usub,perm_u,xsup,grid );

                        cum_nrow =0 ;

                        /* do update with the kth column of L and (k,j)th block of U */
                        lptr = lptr0;
                        luptr = luptr0;

#ifdef _OPENMP
#pragma omp for schedule( SCHEDULE_STRATEGY ) nowait
#endif
                        for (lb = 0; lb < nlb; lb++ ) {
                            int cum_nrow = 0;
                            int temp_nbrow;
                            lptr = lptr0;
                            luptr = luptr0;
                            for (int i = 0; i < lb; ++i) {
                                ib = lsub[lptr];        /* Row block L(i,k). */
                                temp_nbrow = lsub[lptr+1];   /* Number of full rows. */
                                lptr += LB_DESCRIPTOR;  /* Skip descriptor. */
                                lptr += temp_nbrow;
                                luptr += temp_nbrow;
                                cum_nrow +=temp_nbrow;
                            }

                            ib = lsub[lptr];       /* Row block L(i,k). */
                            temp_nbrow = lsub[lptr+1];  /* Number of full rows. */
                            assert(temp_nbrow<=nbrow);

                            lptr += LB_DESCRIPTOR; /* Skip descriptor. */

                            /* Now gather the result into the destination block. */
                            if ( ib < jb ) {  /* A(i,j) is in U. */
                                #ifdef PI_DEBUG
                                    printf("cpu scatter \n");
                                    printf("A(%d,%d) goes to U block %d \n", ib,jb,ljb);
                                #endif

                                tempv = tempv1+cum_nrow;
                                dscatter_u (
						 ib,jb,
						 nsupc,iukp,xsup,
						 klst,nbrow,
						 lptr,temp_nbrow,lsub,
						 usub,tempv,
						 Ufstnz_br_ptr,
						 Unzval_br_ptr,
						 grid
						 );
                            } else {    /* A(i,j) is in L. */
#ifdef PI_DEBUG
                                printf("cpu scatter \n");
                                printf("A(%d,%d) goes to L block %d \n", ib,jb,ljb);
#endif
                                
                                tempv = tempv1+cum_nrow;

                                dscatter_l (
						 ib, ljb,nsupc,iukp,xsup,klst,nbrow,lptr,
						 temp_nbrow,usub,lsub,tempv,
						 indirect_thread,indirect2_thread,
						 Lrowind_bc_ptr,Lnzval_bc_ptr,grid
						 );
                            } /* if ib < jb ... */

                            lptr += temp_nbrow;
                            luptr += temp_nbrow;
                            cum_nrow += temp_nbrow;

                        } /* for lb ... */

                        luptr=luptr0;
                    } /* for j = jjj_st ... */

                    // TAU_STATIC_TIMER_STOP("SPECIAL_CPU_SCATTER");
                } else {
#ifdef _OPENMP
#pragma omp for schedule(SCHEDULE_STRATEGY) nowait
#endif
                    for (j = jjj_st; j < jjj_st+ncpu_blks; ++j) {
                        /* code */
                        #ifdef PI_DEBUG
                            printf("scattering %d  block column\n",j);
                        #endif 

                        /* == processing each of the remaining columns == */
                        if(j==jjj_st) tempv1 = bigV;
                        else tempv1 = bigV + full_u_cols[j-1]*nbrow;

                        arrive_at_ublock( j,&iukp,&rukp,&jb,&ljb,&nsupc,
					  iukp0,rukp0,usub,perm_u,xsup,grid );
                        cum_nrow =0 ;

                        /* do update with the kth column of L and (k,j)th block of U */
                        lptr = lptr0;
                        luptr = luptr0;

                        for (lb = 0; lb < nlb; lb++ ) {
                            ib = lsub[lptr];       /* Row block L(i,k). */
                            temp_nbrow = lsub[lptr+1];  /* Number of full rows. */
                            assert(temp_nbrow<=nbrow);

                            lptr += LB_DESCRIPTOR; /* Skip descriptor. */
#ifdef DGEMM_STAT
			    if(j==jjj_st) {
				temp_ncol = full_u_cols[j];
			    } else {
				temp_ncol = full_u_cols[j]- full_u_cols[j-1];  
			    }
			    printf("%d %d %d \n",temp_nbrow, temp_ncol,ldu);
#endif

			    /* Now gather the result into the destination block. */
			    if ( ib < jb ) {  /* A(i,j) is in U. */
#ifdef PI_DEBUG
				printf("cpu scatter \n");
				printf("A(%d,%d) goes to U block %d \n", ib,jb,ljb);
#endif

				tempv = tempv1+cum_nrow;
                                dscatter_u (
						 ib,jb,
						 nsupc,iukp,xsup,
						 klst,nbrow,
						 lptr,temp_nbrow,lsub,
						 usub,tempv,
						 Ufstnz_br_ptr,
						 Unzval_br_ptr,
						 grid
						 );
			    } else {    /* A(i,j) is in L. */
#ifdef PI_DEBUG
                                printf("cpu scatter \n");
                                printf("A(%d,%d) goes to L block %d \n", ib,jb,ljb);
#endif
                                tempv = tempv1+cum_nrow;

                                dscatter_l (
						 ib, ljb,nsupc,iukp,xsup,klst,nbrow,lptr,
						 temp_nbrow,usub,lsub,tempv,
						 indirect_thread,indirect2_thread,
						 Lrowind_bc_ptr,Lnzval_bc_ptr,grid
						 );
			    } /* if ib < jb ... */

			    lptr += temp_nbrow;
			    luptr += temp_nbrow;
			    cum_nrow += temp_nbrow;
			
			} /* for lb ... */

			luptr=luptr0;
		    } /* for j = jjj_st ... */
		}     /* else if (ncpu_blks >= omp_get_num_threads()) */
	    }         /* parallel region */

	    scatter_timer += SuperLU_timer_() - tstart; 
#ifdef _OPENMP
#pragma omp parallel							\
    private(j,iukp,rukp, tempu, tempv, cum_nrow, jb, nsupc,ljb,		\
	    segsize,lead_zero,						\
	    ib, temp_nbrow,ilst,lib,index,				\
	    ijb,fnz,ucol,rel,ldv,lptrj,luptrj,				\
	    nzval,     lb ,                     jj, i)			\
    firstprivate(luptr,lptr) default (shared)
#endif
            {
                int thread_id = omp_get_thread_num();
        
                int* indirect_thread = indirect + ldt*thread_id;
                int* indirect2_thread = indirect2 + ldt*thread_id;
                double* tempv1;
                for(i = 0; i < num_streams_used; i++) { /* i is private variable */
                    checkCuda(cudaStreamSynchronize (streams[i]));
                    int jjj_st1 = (i==0) ? jjj_st + ncpu_blks : jjj_st + stream_end_col[i-1];
                    int jjj_end = jjj_st + stream_end_col[i];
                    assert(jjj_end-1<nub);
                    assert(jjj_st1>jjj_st) ;

                    /* now scatter it */
#pragma omp for schedule( SCHEDULE_STRATEGY ) nowait 
                    for (j = jjj_st1; j < jjj_end; ++j) {
                        /* code */
#ifdef PI_DEBUG
			printf("scattering %d  block column\n",j);
#endif 
                        /* == processing each of the remaining columns == */

                        if(j==jjj_st) tempv1 = bigV;
                        else tempv1 = bigV + full_u_cols[j-1]*nbrow;

                        arrive_at_ublock( j,&iukp,&rukp,&jb,&ljb,&nsupc,
					  iukp0,rukp0,usub,perm_u,xsup,grid );
                        cum_nrow =0 ;

                        /* do update with the kth column of L and (k,j)th block of U */
                        lptr = lptr0;
                        luptr = luptr0;
                        for (lb = 0; lb < nlb; lb++) {
                            ib = lsub[lptr];       /* Row block L(i,k). */
                            temp_nbrow = lsub[lptr+1];  /* Number of full rows. */
                            assert(temp_nbrow<=nbrow);

                            lptr += LB_DESCRIPTOR; /* Skip descriptor. */
#ifdef DGEMM_STAT
			    if(j==jjj_st) {
				temp_ncol = full_u_cols[j];
			    } else {
				temp_ncol = full_u_cols[j]- full_u_cols[j-1];  
			    }
			    printf("%d %d %d \n",temp_nbrow, temp_ncol,ldu);
#endif

                            /* Now gather the result into the destination block. */
                            if ( ib < jb ) { /* A(i,j) is in U. */
#ifdef PI_DEBUG
				printf("gpu scatter \n");
				printf("A(%d,%d) goes to U block %d \n", ib,jb,ljb);
#endif
                                tempv = tempv1+cum_nrow;
                                dscatter_u (
						 ib,jb,
						 nsupc,iukp,xsup,
						 klst,nbrow,
						 lptr,temp_nbrow,lsub,
						 usub,tempv,
						 Ufstnz_br_ptr,
						 Unzval_br_ptr,
						 grid
						 );
                            } else {    /* A(i,j) is in L. */
#ifdef PI_DEBUG
                                printf("gpu scatter \n");
                                printf("A(%d,%d) goes to L block %d \n", ib,jb,ljb);
#endif
                                tempv = tempv1+cum_nrow;

                                dscatter_l (
						 ib, ljb,nsupc,iukp,xsup,klst,nbrow,lptr,
						 temp_nbrow,usub,lsub,tempv,
						 indirect_thread,indirect2_thread,
						 Lrowind_bc_ptr,Lnzval_bc_ptr,grid
						 );
                            } /* if ib < jb ... */

                            lptr += temp_nbrow;
                            luptr += temp_nbrow;
                            cum_nrow += temp_nbrow;
			    
                        } /* for lb ... */

                        luptr=luptr0;
                    } /* for j = jjj_st ... */
                    
                } /* end for i = 0 to nstreams */
                // TAU_STATIC_TIMER_STOP("GPU_SCATTER");
                // TAU_STATIC_TIMER_STOP("INSIDE_OMP");
            } /* end pragma omp parallel */
            // TAU_STATIC_TIMER_STOP("OUTSIDE_OMP");
        }  /* end while(jjj<nub) */
 
    } /* if nbrow>0 */

 }   /* if msg1 and msg 2 */



