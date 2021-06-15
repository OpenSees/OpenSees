#include "mex.h"

#define MUMPS_ARITH_d 2
#define MUMPS_ARITH_z 8


#if MUMPS_ARITH == MUMPS_ARITH_z

# include "zmumps_c.h"
# define dmumps_c       zmumps_c
# define dmumps_par     zmumps_par
# define DMUMPS_STRUC_C ZMUMPS_STRUC_C
# define DMUMPS_alloc   ZMUMPS_alloc     
# define DMUMPS_free    ZMUMPS_free
# define double2        mumps_double_complex
# define mxREAL2        mxCOMPLEX

#elif MUMPS_ARITH == MUMPS_ARITH_d

# include "dmumps_c.h"
# define double2 double
# define mxREAL2 mxREAL
# define EXTRACT_CMPLX_FROM_C_TO_MATLAB EXTRACT_FROM_C_TO_MATLAB
# define EXTRACT_CMPLX_FROM_MATLAB_TOPTR EXTRACT_FROM_MATLAB_TOPTR

#else

# error "Only d and z arithmetics are supported"

#endif

#define SYM        (prhs[0])
#define JOB        (prhs[1])
#define ICNTL_IN   (prhs[2])
#define CNTL_IN    (prhs[3])
#define PERM_IN    (prhs[4])
#define COLSCA_IN  (prhs[5])
#define ROWSCA_IN  (prhs[6])
#define RHS_IN     (prhs[7])
#define VAR_SCHUR  (prhs[8])
#define INST       (prhs[9])
#define REDRHS_IN  (prhs[10])
#define KEEP_IN    (prhs[11])
#define DKEEP_IN   (prhs[12])
#define A_IN       (prhs[13])

#define INFO_OUT   (plhs[0])
#define RINFO_OUT  (plhs[1])
#define RHS_OUT    (plhs[2])
#define INST_OUT   (plhs[3])
#define SCHUR_OUT  (plhs[4])
#define REDRHS_OUT (plhs[5])
#define PIVNUL_LIST (plhs[6])
#define PERM_OUT   (plhs[7])
#define UNS_PERM   (plhs[8])
#define ICNTL_OUT  (plhs[9])
#define CNTL_OUT   (plhs[10])
#define COLSCA_OUT (plhs[11])
#define ROWSCA_OUT (plhs[12])
#define KEEP_OUT   (plhs[13])
#define DKEEP_OUT  (plhs[14])

#define MYMALLOC(ptr,l,type)                      \
  if(!(ptr = (type *) malloc(l*sizeof(type)))){   \
    mexErrMsgTxt ("Malloc failed in mumpsmex.c"); \
  }                                               

#define MYFREE(ptr) \
    if(ptr){        \
        free(ptr);  \
        ptr = 0;    \
    }

#define EXTRACT_FROM_MATLAB_TOPTR(mxcomponent,mumpspointer,type,length)         \
  ptr_matlab = mxGetPr(mxcomponent);                                            \
  if(ptr_matlab[0] != -9999){                                                   \
    MYFREE(mumpspointer);                                                       \
    MYMALLOC(mumpspointer,length,type);                                         \
    for(i=0;i<length;i++){                                                      \
      mumpspointer[i] = ptr_matlab[i];                                          \
    }                                                                           \
  }


/* For scaling arrays, if they were previously allocated by MUMPS, touch nothing   */
/* This is not quite correct (user may want to modify MUMPS scaling and use given  */
/* scaling, or provide a new scaling vector on input after a previous call where   */
/* it was computed by MUMPS). One way to solve this might be to separate COLSCA_IN */
/* and COLSCA_OUT in the C interface (and possibly Fortran) too, but breaking      */
/* backward compatibility.                                                         */
#define EXTRACT_SCALING_FROM_MATLAB_TOPTR(mxcomponent,mumpspointer,is_a_pointer_from_mumps,length)   \
  ptr_matlab = mxGetPr(mxcomponent);                                            \
  if( ptr_matlab[0] != -9999 && ! (is_a_pointer_from_mumps) ) {                 \
    MYFREE(mumpspointer);                                                       \
    MYMALLOC(mumpspointer,length,double);                                       \
    for(i=0;i<length;i++){                                                      \
      mumpspointer[i] = ptr_matlab[i];                                          \
    }                                                                           \
  }

#define EXTRACT_FROM_MATLAB_TOARR(mxcomponent,mumpsarray,type,length)           \
  ptr_matlab = mxGetPr(mxcomponent);                                            \
  if(ptr_matlab[0] != -9999){                                                   \
    for(i=0;i<length;i++){                                                      \
      if(ptr_matlab[i] != -9998){                                               \
        mumpsarray[i] = ptr_matlab[i];                                          \
      }                                                                         \
    }                                                                           \
  }

#define EXTRACT_FROM_MATLAB_TOVAL(mxcomponent,mumpsvalue)                       \
  ptr_matlab = mxGetPr(mxcomponent);                                            \
  if(ptr_matlab[0] != -9999){                                                   \
      mumpsvalue = ptr_matlab[0];                                               \
  }

#define EXTRACT_FROM_C_TO_MATLAB(mxcomponent,mumpspointer,length)               \
  if(mumpspointer == 0){                                                        \
    mxcomponent = mxCreateDoubleMatrix (1, 1, mxREAL);                          \
    ptr_matlab = mxGetPr (mxcomponent);                                         \
    ptr_matlab[0] = -9999;                                                      \
  }else{                                                                        \
    mxcomponent = mxCreateDoubleMatrix (1,length,mxREAL);                       \
    ptr_matlab = mxGetPr (mxcomponent);                                         \
    for(i=0;i<length;i++){                                                      \
      ptr_matlab[i]=(double)(mumpspointer)[i];                                  \
    }                                                                           \
  }

#if MUMPS_ARITH == MUMPS_ARITH_z

#define EXTRACT_CMPLX_FROM_MATLAB_TOPTR(mxcomponent,mumpspointer,type,length)   \
  ptr_matlab = mxGetPr(mxcomponent);                                            \
  if(ptr_matlab[0] != -9999){                                                   \
    MYFREE(mumpspointer);                                                       \
    MYMALLOC(mumpspointer,length,double2);                                      \
    for(i=0;i<length;i++){                                                      \
      (mumpspointer[i]).r = ptr_matlab[i];                                      \
    }                                                                           \
    ptr_matlab = mxGetPi(mxcomponent);                                          \
    if(ptr_matlab){                                                             \
      for(i=0;i<length;i++){                                                    \
        (mumpspointer[i]).i = ptr_matlab[i];                                    \
      }                                                                         \
    }else{                                                                      \
      for(i=0;i<length;i++){                                                    \
        (mumpspointer[i]).i = 0.0;                                              \
      }                                                                         \
    }                                                                           \
  }


#define EXTRACT_CMPLX_FROM_C_TO_MATLAB(mxcomponent,mumpspointer,length)         \
  if(mumpspointer == 0){                                                        \
    mxcomponent = mxCreateDoubleMatrix (1, 1, mxCOMPLEX);                       \
    ptr_matlab = mxGetPr (mxcomponent);                                         \
    ptr_matlab[0] = -9999;                                                      \
    ptr_matlab = mxGetPi (mxcomponent);                                         \
    ptr_matlab[0] = -9999;                                                      \
  }else{                                                                        \
    mxcomponent = mxCreateDoubleMatrix (1,length,mxCOMPLEX);                    \
    ptr_matlab = mxGetPr (mxcomponent);                                         \
    for(i=0;i<length;i++){                                                      \
      ptr_matlab[i] = (mumpspointer[i]).r;                                      \
    }                                                                           \
    ptr_matlab = mxGetPi (mxcomponent);                                         \
    for(i=0;i<length;i++){                                                      \
      ptr_matlab[i] = (mumpspointer[i]).i;                                      \
    }                                                                           \
  }

#endif

void DMUMPS_free(DMUMPS_STRUC_C **dmumps_par){
  if(*dmumps_par){
  MYFREE( (*dmumps_par)->irn );
  MYFREE( (*dmumps_par)->jcn  );
  MYFREE( (*dmumps_par)->a );
  MYFREE( (*dmumps_par)->irn_loc );
  MYFREE( (*dmumps_par)->jcn_loc );
  MYFREE( (*dmumps_par)->a_loc );
  MYFREE( (*dmumps_par)->eltptr );
  MYFREE( (*dmumps_par)->eltvar );
  MYFREE( (*dmumps_par)->a_elt );
  MYFREE( (*dmumps_par)->perm_in );
  /* colsca/rowsca might have been allocated by
   * MUMPS but in that case the corresponding pointer
   * is already equal to 0 so line below will do nothing */
  MYFREE( (*dmumps_par)->colsca );
  MYFREE( (*dmumps_par)->rowsca  );
  MYFREE( (*dmumps_par)->pivnul_list );
  MYFREE( (*dmumps_par)->listvar_schur );
  MYFREE( (*dmumps_par)->sym_perm );
  MYFREE( (*dmumps_par)->uns_perm );
  MYFREE( (*dmumps_par)->irhs_ptr);
  MYFREE( (*dmumps_par)->irhs_sparse);
  MYFREE( (*dmumps_par)->rhs_sparse);
  MYFREE( (*dmumps_par)->rhs);
  MYFREE( (*dmumps_par)->redrhs);
  MYFREE(*dmumps_par);
  }
}

void DMUMPS_alloc(DMUMPS_STRUC_C **dmumps_par){

  MYMALLOC((*dmumps_par),1,DMUMPS_STRUC_C);
  (*dmumps_par)->irn  = NULL;
  (*dmumps_par)->jcn  = NULL;
  (*dmumps_par)->a  = NULL;
  (*dmumps_par)->irn_loc  = NULL;
  (*dmumps_par)->jcn_loc  = NULL;
  (*dmumps_par)->a_loc  = NULL;
  (*dmumps_par)->eltptr  = NULL;
  (*dmumps_par)->eltvar  = NULL;
  (*dmumps_par)->a_elt  = NULL;
  (*dmumps_par)->perm_in  = NULL;
  (*dmumps_par)->colsca  = NULL;
  (*dmumps_par)->rowsca  = NULL;
  (*dmumps_par)->rhs  = NULL;
  (*dmumps_par)->redrhs  = NULL;
  (*dmumps_par)->rhs_sparse = NULL;
  (*dmumps_par)->irhs_sparse = NULL;
  (*dmumps_par)->irhs_ptr = NULL;
  (*dmumps_par)->pivnul_list  = NULL;
  (*dmumps_par)->listvar_schur  = NULL;
  (*dmumps_par)->schur  = NULL;
  (*dmumps_par)->sym_perm  = NULL;
  (*dmumps_par)->uns_perm  = NULL;
}

void mexFunction(int nlhs, mxArray *plhs[ ],
                 int nrhs, const mxArray *prhs[ ]) { 
  
  int i,j,pos;
  int *ptr_int;
  double *ptr_matlab;
#if MUMPS_ARITH == MUMPS_ARITH_z
  double *ptri_matlab;
#endif
  mwSize tmp_m,tmp_n;

  /* C pointer for input parameters */
  size_t inst_address;
  mwSize n,m,ne, netrue ;
  int job;
  mwIndex *irn_in,*jcn_in;
  
  /* variable for multiple and sparse rhs */
  int posrhs;
          mwSize  nbrhs,ldrhs, nz_rhs;
  mwIndex *irhs_ptr, *irhs_sparse;
  double *rhs_sparse;
#if MUMPS_ARITH == MUMPS_ARITH_z
  double *im_rhs_sparse;
#endif

  DMUMPS_STRUC_C *dmumps_par;
  int dosolve = 0;
  int donullspace = 0;
  int doanalysis = 0;
  int dofactorize = 0;
  
  
  EXTRACT_FROM_MATLAB_TOVAL(JOB,job);

  doanalysis = (job == 1 || job == 4 || job == 6);
  dofactorize = (job == 2 || job == 4 || job == 5 || job == 6);
  dosolve = (job == 3 || job == 5 || job == 6);

  if(job == -1){
    DMUMPS_alloc(&dmumps_par);
    EXTRACT_FROM_MATLAB_TOVAL(SYM,dmumps_par->sym);
    dmumps_par->job = -1;
    dmumps_par->par = 1;
    dmumps_c(dmumps_par);
    dmumps_par->nz = -1;
    dmumps_par->nz_alloc = -1;
  }else{
    EXTRACT_FROM_MATLAB_TOVAL(INST,inst_address);
    ptr_int = (int *) inst_address;

    dmumps_par = (DMUMPS_STRUC_C *) ptr_int;

    if(job == -2){
      dmumps_par->job = -2;
      dmumps_c(dmumps_par);
      /* If colsca/rowsca were freed by MUMPS,
         dmumps_par->colsca/rowsca are now null.
         Application of MYFREE in call below thus ok */
      DMUMPS_free(&dmumps_par);
    }else{

      /* check of input arguments */
      n = mxGetN(A_IN);
      m = mxGetM(A_IN);

      if (!mxIsSparse(A_IN) || n != m )
          mexErrMsgTxt("Input matrix must be a sparse square matrix");
      
      jcn_in = mxGetJc(A_IN);
      ne = jcn_in[n];
      irn_in = mxGetIr(A_IN);
      dmumps_par->n = (int)n;
      if(dmumps_par->n != n)
          mexErrMsgTxt("Input is too big; will not work...barfing out\n");
      
      if(dmumps_par->sym != 0)
          netrue = (n+ne)/2;
      else
          netrue = ne;
      
      if(dmumps_par->nz_alloc < netrue || dmumps_par->nz_alloc >= 2*netrue){  
        MYFREE(dmumps_par->jcn);
        MYFREE(dmumps_par->irn);
        MYFREE(dmumps_par->a);
        MYMALLOC((dmumps_par->jcn),(int)netrue,int);
        MYMALLOC((dmumps_par->irn),(int)netrue,int);
        MYMALLOC((dmumps_par->a),(int)netrue,double2);
        dmumps_par->nz_alloc = (int)netrue;
    if (dmumps_par->nz_alloc != netrue)
        mexErrMsgTxt("Input is too big; will not work...barfing out\n");
      }


      if(dmumps_par->sym == 0){
        /* if analysis already performed then we only need to read
           numerical values
           Note that we suppose that matlab did not change the internal
           format of the matrix between the 2 calls */
        if(doanalysis){ 
          /* || dmumps_par->info[22] == 0 */
          for(i=0;i<dmumps_par->n;i++){
            for(j=jcn_in[i];j<jcn_in[i+1];j++){
              (dmumps_par->jcn)[j] = i+1;
              (dmumps_par->irn)[j] = ((int)irn_in[j])+1;
            }
          }
        }
    dmumps_par->nz = (int)ne;
    if( dmumps_par->nz != ne)
        mexErrMsgTxt("Input is too big; will not work...barfing out\n");
#if MUMPS_ARITH == MUMPS_ARITH_z
        ptr_matlab = mxGetPr(A_IN);
        for(i=0;i<dmumps_par->nz;i++){                                                   
          ((dmumps_par->a)[i]).r = ptr_matlab[i];
        }
        ptr_matlab = mxGetPi(A_IN);
        if(ptr_matlab){
          for(i=0;i<dmumps_par->nz;i++){                                                   
            ((dmumps_par->a)[i]).i = ptr_matlab[i];
          }
        }else{
          for(i=0;i<dmumps_par->nz;i++){                                                   
             ((dmumps_par->a)[i]).i = 0.0;
             }
        }
#else
        ptr_matlab = mxGetPr(A_IN);
        for(i=0;i<dmumps_par->nz;i++){                                                   
          (dmumps_par->a)[i] = ptr_matlab[i];
        }
#endif
      }else{
        /* in the symmetric case we do not need to check doanalysis */
        pos = 0;
        ptr_matlab = mxGetPr(A_IN);
#if MUMPS_ARITH == MUMPS_ARITH_z
        ptri_matlab = mxGetPi(A_IN);
#endif
        for(i=0;i<dmumps_par->n;i++){
          for(j=jcn_in[i];j<jcn_in[i+1];j++){
            if(irn_in[j] >= i){
              if(pos >= netrue)
              mexErrMsgTxt("Input matrix must be symmetric");
              (dmumps_par->jcn)[pos] = i+1;
              (dmumps_par->irn)[pos] = (int)irn_in[j]+1;
#if MUMPS_ARITH == MUMPS_ARITH_z
              ((dmumps_par->a)[pos]).r = ptr_matlab[j];
              if(ptri_matlab){
                ((dmumps_par->a)[pos]).i = ptri_matlab[j];
              }else{
                ((dmumps_par->a)[pos]).i = 0.0;
              }
#else
              (dmumps_par->a)[pos] = ptr_matlab[j];
#endif
              pos++;
             }
          }
        }
        dmumps_par->nz = pos;
      }
    

      EXTRACT_FROM_MATLAB_TOVAL(JOB,dmumps_par->job);
      EXTRACT_FROM_MATLAB_TOARR(ICNTL_IN,dmumps_par->icntl,int,60);
      EXTRACT_FROM_MATLAB_TOARR(CNTL_IN,dmumps_par->cntl,double,15);
      EXTRACT_FROM_MATLAB_TOPTR(PERM_IN,(dmumps_par->perm_in),int,((int)n));

      /* colsca and rowsca are treated differently: it may happen that
         dmumps_par-> colsca is nonzero because it was set to a nonzero
         value on output (COLSCA_OUT) from MUMPS. Unfortunately if scaling
         was on output, one cannot currently provide scaling on input
         afterwards without reinitializing the instance */

      EXTRACT_SCALING_FROM_MATLAB_TOPTR(COLSCA_IN,(dmumps_par->colsca),(dmumps_par->colsca_from_mumps),((int)n)); /* type always double */
      EXTRACT_SCALING_FROM_MATLAB_TOPTR(ROWSCA_IN,(dmumps_par->rowsca),(dmumps_par->rowsca_from_mumps),((int)n)); /* type always double */

      EXTRACT_FROM_MATLAB_TOARR(KEEP_IN,dmumps_par->keep,int,500);
      EXTRACT_FROM_MATLAB_TOARR(DKEEP_IN,dmumps_par->dkeep,double,230);

      dmumps_par->size_schur = (int)mxGetN(VAR_SCHUR);
      EXTRACT_FROM_MATLAB_TOPTR(VAR_SCHUR,(dmumps_par->listvar_schur),int,dmumps_par->size_schur);
      if(!dmumps_par->listvar_schur) dmumps_par->size_schur = 0;

      ptr_matlab = mxGetPr (RHS_IN);

/*
 * To follow the "spirit" of the Matlab/Scilab interfaces, treat case of null
 * space separately. In that case, we initialize lrhs and nrhs, automatically
 * allocate the space needed, and do not rely on what is provided by the user
 * in component RHS, that is not touched.
 *
 * Note that, at the moment, the user should not call the solution step combined
 * with the factorization step when he/she sets icntl[25-1] to a non-zero value.
 * Hence we suppose in the following that infog[28-1] is available and that we
 * can use it.
 * 
 * For users of scilab/matlab, it would still be nice to be able to set ICNTL(25)=-1,
 * and use JOB=6. If we want to make such a feature available, we should
 * call separately job=2 and job=3 even if job=5 or 6 and set nbrhs (and allocate
 * space correctly) between job=2 and job=3 calls to MUMPS.
 *
 */
      if ( dmumps_par->icntl[25-1] == -1 && dmumps_par->infog[28-1] > 0 ) {
          dmumps_par->nrhs=dmumps_par->infog[28-1];
          donullspace = dosolve;
         }
      else if ( dmumps_par->icntl[25-1] > 0 && dmumps_par->icntl[25-1] <= dmumps_par->infog[28-1] ) {
           dmumps_par->nrhs=1;
           donullspace = dosolve;
         }
      else {
           donullspace=0;
         }
      if (donullspace) {
        nbrhs=dmumps_par->nrhs; ldrhs=n;
        dmumps_par->lrhs=(int)n;
        MYMALLOC((dmumps_par->rhs),((dmumps_par->n)*(dmumps_par->nrhs)),double2);
         }
      else if((!dosolve) || ptr_matlab[0] == -9999 ) { /* rhs not already provided, or not used */
/*     Case where dosolve is true and ptr_matlab[0]=-9999, this could cause problems:
 *        1/ RHS was not initialized while it should have been
 *        2/ RHS was explicitely initialized to -9999 but is not allocated of the right size
 */
        EXTRACT_CMPLX_FROM_MATLAB_TOPTR(RHS_IN,(dmumps_par->rhs),double,1);
      }else{
        nbrhs = mxGetN(RHS_IN);
        ldrhs = mxGetM(RHS_IN);
        dmumps_par->nrhs = (int)nbrhs;
        dmumps_par->lrhs = (int)ldrhs;
        if(ldrhs != n){
          mexErrMsgTxt ("Incompatible number of rows in RHS");
        }
        if (!mxIsSparse(RHS_IN)){ /* full rhs */
          dmumps_par->icntl[20-1] = 0;
          EXTRACT_CMPLX_FROM_MATLAB_TOPTR(RHS_IN,(dmumps_par->rhs),double,(int)( dmumps_par->nrhs*ldrhs));
        }else{ /* sparse rhs */
          /* printf("sparse RHS ldrhs = %d nrhs = %d\n",ldrhs,nbrhs); */
          if (dmumps_par->icntl[30-1] == 0) {
            /* A-1 feature was not requested => we are in the standard
             * sparse RHS case and thus we set ICNTL(20) accordingly. */
            dmumps_par->icntl[20-1] = 1;
          }
          irhs_ptr = mxGetJc(RHS_IN);
          irhs_sparse = mxGetIr(RHS_IN);
          rhs_sparse = mxGetPr(RHS_IN);
#if MUMPS_ARITH == MUMPS_ARITH_z
          im_rhs_sparse = mxGetPi(RHS_IN);
#endif
          nz_rhs = irhs_ptr[nbrhs];
          dmumps_par->nz_rhs = (int)nz_rhs;

          MYMALLOC((dmumps_par->irhs_ptr),(dmumps_par->nrhs+1),int);
          MYMALLOC((dmumps_par->irhs_sparse), dmumps_par->nz_rhs,int);
          MYMALLOC((dmumps_par->rhs_sparse), dmumps_par->nz_rhs,double2);
          /* dmumps_par->rhs will store the solution*/
          MYMALLOC((dmumps_par->rhs),((dmumps_par->nrhs*dmumps_par->lrhs)),double2);

          for(i=0;i< dmumps_par->nrhs;i++){
            for(j=irhs_ptr[i];j<irhs_ptr[i+1];j++){
              (dmumps_par->irhs_sparse)[j] = irhs_sparse[j]+1;
            }
            (dmumps_par->irhs_ptr)[i] = irhs_ptr[i]+1;
          }
          (dmumps_par->irhs_ptr)[dmumps_par->nrhs] = dmumps_par->nz_rhs+1;
#if MUMPS_ARITH == MUMPS_ARITH_z
          if(im_rhs_sparse){
            for(i=0;i<dmumps_par->nz_rhs;i++){                                                   
              ((dmumps_par->rhs_sparse)[i]).r = rhs_sparse[i];
              ((dmumps_par->rhs_sparse)[i]).i = im_rhs_sparse[i];
            }
          }else{
            for(i=0;i<dmumps_par->nz_rhs;i++){                                                   
              ((dmumps_par->rhs_sparse)[i]).r = rhs_sparse[i];
              ((dmumps_par->rhs_sparse)[i]).i = 0.0;
            }
          }
#else
          for(i=0;i<dmumps_par->nz_rhs;i++){                                                   
            (dmumps_par->rhs_sparse)[i] = rhs_sparse[i];
          }
#endif
        }
      }

      if(dmumps_par->size_schur > 0){
        if (dofactorize) {
          MYMALLOC((dmumps_par->schur),((dmumps_par->size_schur)*(dmumps_par->size_schur)),double2);
        }
        dmumps_par->icntl[18] = 1;
      }else{
        dmumps_par->icntl[18] = 0;
      }
       /* Reduced RHS */
       if ( dmumps_par->size_schur > 0 && dosolve ) {
          if ( dmumps_par->icntl[26-1] == 2 ) {
            /* REDRHS is on input */
            tmp_m= mxGetM(REDRHS_IN);
            tmp_n= mxGetN(REDRHS_IN);
            if (tmp_m != dmumps_par->size_schur || tmp_n != dmumps_par->nrhs) {
              mexErrMsgTxt ("bad dimensions for REDRHS in mumpsmex.c");
            }
            EXTRACT_CMPLX_FROM_MATLAB_TOPTR(REDRHS_IN,(dmumps_par->redrhs),double,((int)tmp_m*tmp_n));
            dmumps_par->lredrhs=dmumps_par->size_schur;
          }
          if ( dmumps_par->icntl[26-1] == 1 ) {
            /* REDRHS on output. Must be allocated before the call */
            MYFREE(dmumps_par->redrhs);
            if(!(dmumps_par->redrhs=(double2 *)malloc((dmumps_par->size_schur*dmumps_par->nrhs)*sizeof(double2)))){
              mexErrMsgTxt("malloc redrhs failed in intmumpsc.c");
            }
          }
       }
      dmumps_c(dmumps_par);
    }
  }
  if(nlhs > 0){
    EXTRACT_FROM_C_TO_MATLAB( INFO_OUT  ,(dmumps_par->infog),80);
    EXTRACT_FROM_C_TO_MATLAB( RINFO_OUT ,(dmumps_par->rinfog),40);
    /* A-1 on output */
    if ( dmumps_par->icntl[30-1] != 0 && dosolve ) {
      RHS_OUT = mxCreateSparse(dmumps_par->n, dmumps_par->n,dmumps_par->nz_rhs,mxREAL2);

      irhs_ptr = mxGetJc(RHS_OUT);
      irhs_sparse = mxGetIr(RHS_OUT);
      for(j=0;j<dmumps_par->nrhs+1;j++){
         irhs_ptr[j] = (mwIndex) ((dmumps_par->irhs_ptr)[j]-1);
      }
      ptr_matlab = mxGetPr(RHS_OUT);
#if MUMPS_ARITH == MUMPS_ARITH_z
      ptri_matlab = mxGetPi(RHS_OUT);
#endif
      for(i=0;i<dmumps_par->nz_rhs;i++){
#if MUMPS_ARITH == MUMPS_ARITH_z
        /* complex arithmetic */
        ptr_matlab[i] = (dmumps_par->rhs_sparse)[i].r;
        ptri_matlab[i] = (dmumps_par->rhs_sparse)[i].i;
#else
        /* real arithmetic */
        ptr_matlab[i] = (dmumps_par->rhs_sparse)[i];
#endif
        irhs_sparse[i] = (mwIndex)((dmumps_par->irhs_sparse)[i]-1);
      }

    }
    else if(dmumps_par->rhs && dosolve){
      /* nbrhs may not have been set (case of null space) */
      nbrhs=dmumps_par->nrhs;
      RHS_OUT = mxCreateDoubleMatrix (dmumps_par->n,dmumps_par->nrhs,mxREAL2);
      ptr_matlab = mxGetPr (RHS_OUT);
#if MUMPS_ARITH == MUMPS_ARITH_z
      ptri_matlab = mxGetPi (RHS_OUT);
      for(j=0;j<dmumps_par->nrhs;j++){
        posrhs = j*(int)n;
        for(i=0;i<dmumps_par->n;i++){
          ptr_matlab[posrhs+i]= (dmumps_par->rhs)[posrhs+i].r;
          ptri_matlab[posrhs+i]= (dmumps_par->rhs)[posrhs+i].i;
        }              
      }
#else
      for(j=0;j<dmumps_par->nrhs;j++){
        posrhs = j*dmumps_par->n;
        for(i=0;i<dmumps_par->n;i++){
          ptr_matlab[posrhs+i]= (dmumps_par->rhs)[posrhs+i];
        }              
      }
#endif
    }else{
      EXTRACT_CMPLX_FROM_C_TO_MATLAB( RHS_OUT,(dmumps_par->rhs),1);
    }

    ptr_int = (int *)dmumps_par;
    inst_address = (size_t) ptr_int;
    EXTRACT_FROM_C_TO_MATLAB( INST_OUT   ,&inst_address,1); 
    EXTRACT_FROM_C_TO_MATLAB( PIVNUL_LIST,dmumps_par->pivnul_list,dmumps_par->infog[27]);
    EXTRACT_FROM_C_TO_MATLAB( PERM_OUT   ,dmumps_par->sym_perm,dmumps_par->n);
    EXTRACT_FROM_C_TO_MATLAB( UNS_PERM   ,dmumps_par->uns_perm,dmumps_par->n);
    EXTRACT_FROM_C_TO_MATLAB( ICNTL_OUT  ,dmumps_par->icntl,60);
    EXTRACT_FROM_C_TO_MATLAB( CNTL_OUT   ,dmumps_par->cntl,15);
    EXTRACT_FROM_C_TO_MATLAB( ROWSCA_OUT ,dmumps_par->rowsca,dmumps_par->n);
    EXTRACT_FROM_C_TO_MATLAB( COLSCA_OUT ,dmumps_par->colsca,dmumps_par->n);
    EXTRACT_FROM_C_TO_MATLAB( KEEP_OUT   ,dmumps_par->keep,500);
    EXTRACT_FROM_C_TO_MATLAB( DKEEP_OUT  ,dmumps_par->dkeep,230);

    if(dmumps_par->size_schur > 0 && dofactorize){
      SCHUR_OUT = mxCreateDoubleMatrix(dmumps_par->size_schur,dmumps_par->size_schur,mxREAL2);
      ptr_matlab = mxGetPr (SCHUR_OUT);
#if MUMPS_ARITH == MUMPS_ARITH_z
      ptri_matlab = mxGetPi (SCHUR_OUT);
      for(i=0;i<dmumps_par->size_schur;i++){
        pos = i*(dmumps_par->size_schur);
        for(j=0;j<dmumps_par->size_schur;j++){
          ptr_matlab[j+pos] = ((dmumps_par->schur)[j+pos]).r;
          ptri_matlab[j+pos] = ((dmumps_par->schur)[j+pos]).i;
        }
      }
#else
      for(i=0;i<dmumps_par->size_schur;i++){
        pos = i*(dmumps_par->size_schur);
        for(j=0;j<dmumps_par->size_schur;j++){
          ptr_matlab[j+pos] = (dmumps_par->schur)[j+pos];
        }
      }
#endif
    }else{
      SCHUR_OUT = mxCreateDoubleMatrix(1,1,mxREAL2);
      ptr_matlab = mxGetPr (SCHUR_OUT);
      ptr_matlab[0] = -9999; 
#if MUMPS_ARITH == MUMPS_ARITH_z
      ptr_matlab = mxGetPi (SCHUR_OUT);
      ptr_matlab[0] = -9999;
#endif 
    }
    /* REDRHS on output */
    if ( dmumps_par->icntl[26-1]==1 && dmumps_par->size_schur > 0 && dosolve ) {
      REDRHS_OUT = mxCreateDoubleMatrix(dmumps_par->size_schur,dmumps_par->nrhs,mxREAL2);
      ptr_matlab = mxGetPr(REDRHS_OUT);
#if MUMPS_ARITH == MUMPS_ARITH_z
      ptri_matlab = mxGetPi (REDRHS_OUT);
#endif
      for(i=0;i<dmumps_par->nrhs*dmumps_par->size_schur;i++){
#if MUMPS_ARITH == MUMPS_ARITH_z
        ptr_matlab[i] = ((dmumps_par->redrhs)[i]).r;
        ptri_matlab[i] = ((dmumps_par->redrhs)[i]).i;
#else
        ptr_matlab[i] = ((dmumps_par->redrhs)[i]);
#endif
      }
    }else{
      REDRHS_OUT = mxCreateDoubleMatrix(1,1,mxREAL2);
      ptr_matlab = mxGetPr (REDRHS_OUT);
      ptr_matlab[0] = -9999; 
#if MUMPS_ARITH == MUMPS_ARITH_z
      ptr_matlab = mxGetPi (REDRHS_OUT);
      ptr_matlab[0] = -9999;
#endif 
    }

    MYFREE(dmumps_par->redrhs);
    MYFREE(dmumps_par->schur);
    MYFREE(dmumps_par->irhs_ptr);
    MYFREE(dmumps_par->irhs_sparse);
    MYFREE(dmumps_par->rhs_sparse);
    MYFREE(dmumps_par->rhs);
  }
}
