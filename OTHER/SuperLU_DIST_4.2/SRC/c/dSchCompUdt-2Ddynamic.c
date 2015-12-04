

/*! @file 
 * \brief THis file contains the main loop of pdgstrf which involves rank k
 *        update of the Schur complement.
 *        Uses 2D partitioning for the scatter phase.
 *
 * <pre>
 * -- Distributed SuperLU routine (version 4.1) --
 * Lawrence Berkeley National Lab, Univ. of California Berkeley.
 * October 1, 2014
 *
 */

#define SCHEDULE_STRATEGY guided 

double tt_start;
double tt_end;
if ( msg0 && msg2 ) {  /* L(:,k) and U(k,:) are not empty. */
    int cum_nrow=0;
    int temp_nbrow;
    lptr = lptr0;
    luptr = luptr0;
    /**
     * seperating L blocks
     */
     int lookAheadBlk=0, RemainBlk=0;

     tt_start = SuperLU_timer_();

     for (int i = 0; i < nlb; ++i) {
	 ib = lsub[lptr];            /* Row block L(i,k). */
	 temp_nbrow = lsub[lptr+1];  /* Number of full rows. */
        
	 int look_up_flag=1;
	 for (int j = k0+1; j < SUPERLU_MIN (k0 + num_look_aheads+2, nsupers ); ++j)
	     {
		 if(ib == perm_c_supno[j]) look_up_flag=0;
	     }
	 
	 if(!look_up_flag) {
	     /* ib is within look up window */
	     if (lookAheadBlk==0) {
		 lookAheadFullRow[lookAheadBlk] = temp_nbrow;
	     } else {
		 lookAheadFullRow[lookAheadBlk] = temp_nbrow+lookAheadFullRow[lookAheadBlk-1];   
	     }
	     lookAheadStRow[lookAheadBlk] = cum_nrow;
	     lookAhead_lptr[lookAheadBlk] = lptr;
	     lookAhead_ib[lookAheadBlk] = ib; 
	     lookAheadBlk++;
	 } else { /* ib is not in look up window */

	     if (RemainBlk==0) {
		 Remain_info[RemainBlk].FullRow = temp_nbrow;
	     } else {
		 Remain_info[RemainBlk].FullRow = temp_nbrow+Remain_info[RemainBlk-1].FullRow;   
	     }

             RemainStRow[RemainBlk] = cum_nrow;
             // Remain_lptr[RemainBlk] = lptr;
	     Remain_info[RemainBlk].lptr = lptr;
	     // Remain_ib[RemainBlk] = ib; 
	     Remain_info[RemainBlk].ib = ib; 
	     RemainBlk++;
	 }
	 
         cum_nrow +=temp_nbrow;
	 
	 lptr += LB_DESCRIPTOR;  /* Skip descriptor. */
	 lptr += temp_nbrow;
	 luptr += temp_nbrow;
     }  /* for i ... */

     lptr = lptr0;
     luptr = luptr0;

     /* leading dimension of L buffer */
#if 0
     int LDlookAhead_LBuff = lookAheadFullRow[lookAheadBlk-1]; /* may go negative.*/
#else /* Piyush fix */
     int LDlookAhead_LBuff = lookAheadBlk==0? 0 :lookAheadFullRow[lookAheadBlk-1];
#endif

     /* #pragma omp parallel for  */
     for (int i = 0; i < lookAheadBlk; ++i) {
	 int StRowDest  = 0;
	 int temp_nbrow;
	 if (i==0) {
	     temp_nbrow = lookAheadFullRow[0];
	 } else {
	     StRowDest   = lookAheadFullRow[i-1];
	     temp_nbrow  = lookAheadFullRow[i]-lookAheadFullRow[i-1];
	 }
	 
	 int StRowSource=lookAheadStRow[i];
	 
	 /* Now copying the matrix*/
	 // #pragma omp parallel for (gives slow down)
	 for (int j = 0; j < knsupc; ++j) {
	     memcpy(&lookAhead_L_buff[StRowDest+j*LDlookAhead_LBuff],
		    &lusup[luptr+j*nsupr+StRowSource],
		    temp_nbrow * sizeof(double) );
	 }
     }

     int LDRemain_LBuff = RemainBlk==0 ? 0 : Remain_info[RemainBlk-1].FullRow;
#ifdef _OPENMP
#pragma omp parallel for 
#endif
     for (int i = 0; i < RemainBlk; ++i) {
	 int StRowDest  = 0;
	 int temp_nbrow;
         if (i==0)  {
	     temp_nbrow = Remain_info[0].FullRow;
	 } else  {
	     StRowDest   = Remain_info[i-1].FullRow;
	     temp_nbrow  = Remain_info[i].FullRow-Remain_info[i-1].FullRow;
	 }

	 int StRowSource=RemainStRow[i];

	 /* Now copying the matrix*/
	 // #pragma omp parallel for (gives slow down)
	 for (int j = 0; j < knsupc; ++j) {
	     // printf("StRowDest %d LDRemain_LBuff %d StRowSource %d \n", StRowDest ,LDRemain_LBuff ,StRowSource );
	     memcpy(&Remain_L_buff[StRowDest+j*LDRemain_LBuff],
		    &lusup[luptr+j*nsupr+StRowSource],
                    temp_nbrow * sizeof(double) );
	 }
     } /* parallel for i ... */

     tt_end = SuperLU_timer_();

     LookAheadRowSepTimer += tt_end-tt_start;
#if 0
     LookAheadRowSepMOP  +=  2*knsupc*(lookAheadFullRow[lookAheadBlk-1]+Remain_info[RemainBlk-1].FullRow );
#else
     int_t lnbrow, rnbrow;
     lnbrow = lookAheadBlk==0 ? 0  : lookAheadFullRow[lookAheadBlk-1];
     rnbrow = RemainBlk==0 ? 0 : Remain_info[RemainBlk-1].FullRow;
     nbrow = lnbrow + rnbrow;
     LookAheadRowSepMOP += 2*knsupc*(nbrow);
#endif     
     ldu   =0;
     full  =1;
     
     /*updating lookahead rows */
     
     tt_start = SuperLU_timer_();
#if 0     
     nbrow = lookAheadFullRow[lookAheadBlk-1]+Remain_info[RemainBlk-1].FullRow;
#endif

     if ( nbrow>0 ) {
	 /*
	  *   counting U blocks
	  */
	 ncols=0;
	 ldu=0;
	 full=1;
	 int temp_ncols=0;
	 for (j = jj0; j < nub; ++j) {
	     temp_ncols=0;
	     arrive_at_ublock(
			      j,&iukp,&rukp,&jb,&ljb,&nsupc,
			      iukp0,rukp0,usub,perm_u,xsup,grid
			      );
	     Ublock_info[j].iukp = iukp;
	     Ublock_info[j].rukp = rukp;
	     Ublock_info[j].jb = jb;
	     
	     /* Prepare to call GEMM. */
	     jj = iukp;
	     
	     for (; jj < iukp+nsupc; ++jj) {
		 segsize = klst - usub[jj];
		 if ( segsize ) {
                    ++temp_ncols;
                    if ( segsize != ldu ) full = 0;
                    if ( segsize > ldu ) ldu = segsize;
		 }
	     }

	     Ublock_info[j].full_u_cols = temp_ncols;
	     ncols += temp_ncols;
	 }

	 /* Now doing prefix sum on  on full_u_cols */
	 for ( j = jj0+1; j < nub; ++j) {
	     Ublock_info[j].full_u_cols += Ublock_info[j-1].full_u_cols;
	 }
            
	 tempu = bigU;

#ifdef _OPENMP        
#pragma omp parallel for private(j,iukp,rukp,tempu, jb, nsupc,ljb,segsize,\
	lead_zero, jj, i) \
        default (shared) schedule(SCHEDULE_STRATEGY)
#endif
        for (j = jj0; j < nub; ++j) {

            if(j==jj0) tempu = bigU;
            else tempu = bigU + ldu*Ublock_info[j-1].full_u_cols;

            /* == processing each of the remaining columns == */
            arrive_at_ublock(j,&iukp,&rukp,&jb,&ljb,&nsupc,
			     iukp0,rukp0,usub,perm_u,xsup,grid);
            
            for (jj = iukp; jj < iukp+nsupc; ++jj) {
                segsize = klst - usub[jj];
                if ( segsize ) {
                    lead_zero = ldu - segsize;
                    for (i = 0; i < lead_zero; ++i) tempu[i] = zero;
                    tempu += lead_zero;
                    for (i = 0; i < segsize; ++i) tempu[i] = uval[rukp+i];
                    rukp += segsize;
                    tempu += segsize;
                }
            }

            rukp -= usub[iukp - 1]; /* Return to start of U(k,j). */

        }   /* parallel for j:jjj_st..jjj */

        tempu = bigU;       //setting it to starting of the matrix

    }  /* if(nbrow>0) */

    tt_end = SuperLU_timer_();
    GatherTimer += tt_end-tt_start;
    GatherMOP += 2*ldu*ncols;

    int Lnbrow   = lookAheadBlk==0 ? 0 :lookAheadFullRow[lookAheadBlk-1];
    int Rnbrow   = RemainBlk==0 ? 0 : Remain_info[RemainBlk-1].FullRow;
    int jj_cpu=nub;       /*limit between CPU and GPU */
    tempv = bigV;

    if (Lnbrow>0 && ldu >0 && ncols>0) {
        ncols   = Ublock_info[nub-1].full_u_cols;
        schur_flop_counter  += 2 * (double)Lnbrow * (double)ldu * (double)ncols;
        stat->ops[FACT]     += 2 * (double)Lnbrow * (double)ldu * (double)ncols;

        tt_start = SuperLU_timer_();

#ifdef _OPENMP
#pragma omp parallel for default (shared) \
    private (j,i,lb,rukp,iukp,jb,nsupc,ljb,lptr,ib,temp_nbrow,cum_nrow)	\
    schedule(dynamic)
#endif
        for (int ij = 0; ij < lookAheadBlk*(nub-jj0); ++ij) {
            int j   = ij/lookAheadBlk + jj0; 
            int lb  = ij%lookAheadBlk;

#ifdef _OPENMP            
            int thread_id = omp_get_thread_num();
#else
            int thread_id = 0;
#endif
            int* indirect_thread    = indirect + ldt*thread_id;
            int* indirect2_thread   = indirect2 + ldt*thread_id;
            double* tempv1 = bigV + thread_id*ldt*ldt; 

            /* Getting U block information */
            /* unsigned long long ut_start, ut_end; */
            int_t rukp =  Ublock_info[j].rukp;
            int_t iukp =  Ublock_info[j].iukp;
            int jb   =  Ublock_info[j].jb;
            int nsupc = SuperSize(jb);
            int ljb = LBj (jb, grid);
            int st_col;
            int ncols;
            if (j>jj0) {
                ncols  = Ublock_info[j].full_u_cols-Ublock_info[j-1].full_u_cols;
                st_col = Ublock_info[j-1].full_u_cols;
            } else {
                ncols  = Ublock_info[j].full_u_cols;
                st_col = 0;   
            }

            /* Getting L block information */
            int_t lptr = lookAhead_lptr[lb];
            int ib   = lookAhead_ib[lb];
            int temp_nbrow = lsub[lptr+1];
            lptr += LB_DESCRIPTOR;
            int cum_nrow = (lb==0 ? 0 : lookAheadFullRow[lb-1]);

#if defined (USE_VENDOR_BLAS)            
            dgemm_("N", "N", &temp_nbrow, &ncols, &ldu, &alpha,
                  &lookAhead_L_buff[(knsupc-ldu)*Lnbrow+cum_nrow], &Lnbrow,
                  &tempu[st_col*ldu], &ldu, &beta, tempv1, &temp_nbrow, 1, 1);
#else
            dgemm_("N", "N", &temp_nbrow, &ncols, &ldu, &alpha,
                  &lookAhead_L_buff[(knsupc-ldu)*Lnbrow+cum_nrow], &Lnbrow,
                  &tempu[st_col*ldu], &ldu, &beta, tempv1, &temp_nbrow);
#endif
            if ( ib < jb ) {
                dscatter_u (
				 ib, jb,
				 nsupc, iukp,xsup,
				 klst, temp_nbrow,
				 lptr, temp_nbrow,lsub,
				 usub, tempv1,
				 Ufstnz_br_ptr,
				 Unzval_br_ptr,
				 grid
				 );
            } else {
                dscatter_l (
				 ib, ljb, nsupc,iukp,xsup,klst,temp_nbrow,lptr,
				 temp_nbrow,usub,lsub,tempv1,
				 indirect_thread, indirect2_thread,
				 Lrowind_bc_ptr,Lnzval_bc_ptr,grid
				 );
            }
        } /* for ij = ... */

        tt_end = SuperLU_timer_();
        LookAheadGEMMTimer += tt_end- tt_start;
        LookAheadGEMMFlOp  += 2 * (double ) Lnbrow * (double )ldu * (double )ncols;
        stat->ops[FACT]    += 2 * (double ) Lnbrow * (double )ldu * (double )ncols;
        LookAheadScatterTimer += tt_end-tt_start;
        LookAheadScatterMOP   += 3*Lnbrow*ncols;
    } /* if Lnbrow < ... */
    
    /***************************************************************
     *   Updating remaining rows and column on CPU
     ***************************************************************/
    Rnbrow  = RemainBlk==0 ? 0 : Remain_info[RemainBlk-1].FullRow;
    ncols   = jj_cpu==0 ? 0 : Ublock_info[jj_cpu-1].full_u_cols;

    schur_flop_counter  += 2 * (double)Rnbrow * (double)ldu * (double)ncols;
    stat->ops[FACT]     += 2 * (double)Rnbrow * (double)ldu * (double)ncols;

    tt_start = SuperLU_timer_();

#ifdef _OPENMP
#pragma omp parallel for default (shared) \
    private (j,i,lb,rukp,iukp,jb,nsupc,ljb,lptr,ib,temp_nbrow,cum_nrow)	\
    schedule(dynamic)
#endif
    for (int ij = 0; ij < RemainBlk*(jj_cpu-jj0); ++ij) {
	int j   = ij / RemainBlk + jj0; 
	int lb  = ij % RemainBlk;

#ifdef _OPENMP            
	int thread_id = omp_get_thread_num();
#else
	int thread_id = 0;
#endif
	int* indirect_thread = indirect + ldt*thread_id;
	int* indirect2_thread = indirect2 + ldt*thread_id;
	double* tempv1 = bigV + thread_id*ldt*ldt; 

	/* Getting U block information */
	/* unsigned long long ut_start, ut_end; */
	int_t rukp =  Ublock_info[j].rukp;
	int_t iukp =  Ublock_info[j].iukp;
	int jb   =  Ublock_info[j].jb;
	int nsupc = SuperSize(jb);
	int ljb = LBj (jb, grid);
	int st_col;
	int ncols;
	if (j>jj0) {
	    ncols  = Ublock_info[j].full_u_cols-Ublock_info[j-1].full_u_cols;
	    st_col = Ublock_info[j-1].full_u_cols;
	} else {
	    ncols  = Ublock_info[j].full_u_cols;
	    st_col = 0;   
	}

	/* Getting L block information */
	int_t lptr = Remain_info[lb].lptr;
	int ib   = Remain_info[lb].ib;
	int temp_nbrow = lsub[lptr+1];
	lptr += LB_DESCRIPTOR;
	int cum_nrow = (lb==0 ? 0 : Remain_info[lb-1].FullRow);

	/* calling GEMM */
#if defined (USE_VENDOR_BLAS)
	dgemm_("N", "N", &temp_nbrow, &ncols, &ldu, &alpha,
	      &Remain_L_buff[(knsupc-ldu)*Rnbrow+cum_nrow], &Rnbrow,
	      &tempu[st_col*ldu], &ldu, &beta, tempv1, &temp_nbrow, 1, 1);
#else
	dgemm_("N", "N", &temp_nbrow, &ncols, &ldu, &alpha,
	      &Remain_L_buff[(knsupc-ldu)*Rnbrow+cum_nrow], &Rnbrow,
	      &tempu[st_col*ldu], &ldu, &beta, tempv1, &temp_nbrow);
#endif

	/* Now scattering the block */
	if ( ib<jb ) {
	    dscatter_u (
			     ib, jb,
			     nsupc, iukp,xsup,
			     klst, temp_nbrow,
			     lptr, temp_nbrow,lsub,
			     usub, tempv1,
			     Ufstnz_br_ptr,
			     Unzval_br_ptr,
			     grid
			     );
	} else {
	    dscatter_l (
			     ib, ljb, nsupc,iukp,xsup,klst,temp_nbrow,lptr,
			     temp_nbrow,usub,lsub,tempv1,
			     indirect_thread, indirect2_thread,
			     Lrowind_bc_ptr,Lnzval_bc_ptr,grid
			     );
	}
    } /* for (int ij =... */
        
}  /* if  k L(:,k) and U(k,:) are not empty */
