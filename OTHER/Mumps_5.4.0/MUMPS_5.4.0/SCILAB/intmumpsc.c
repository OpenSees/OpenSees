#include "mex.h"
#include "stack-c.h"
#include "sci_gateway.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

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

#elif MUMPS_ARITH == MUMPS_ARITH_d

# include "dmumps_c.h"
# define double2 double
# define EXTRACT_CMPLX_FROM_C_TO_SCILAB EXTRACT_DOUBLE_FROM_C_TO_SCILAB
# define EXTRACT_CMPLX_FROM_SCILAB_TOPTR EXTRACT_FROM_SCILAB_TOPTR

#else

# error "Only d and z arithmetics are supported"

#endif


#define nb_RHS 12

#define MYFREE(ptr)\
if(ptr){	\
  free(ptr);  \
  ptr=0;} 	\


#define EXTRACT_FROM_SCILAB_TOPTR(it,ptr_scilab1,ptr_scilab2,mumpspointer,type,length)\
if(ptr_scilab1[0] != -9999){                                            	\
  free(mumpspointer);                                                   	\
  mumpspointer = (type *) malloc(length*sizeof(type));  			\
  for(i=0;i<length;i++){                                        		\
  mumpspointer[i] = ptr_scilab1[i];                             		\
  }                                                             		\
}

#define EXTRACT_FROM_SCILAB_TOARR(ptr_scilab,mumpsarray,type,length)    \
if(ptr_scilab[0] != -9999){                                            	\
  for(i=0;i<length;i++){                                                \
    if(ptr_scilab[i] != -9998){                                         \
      mumpsarray[i] = ptr_scilab[i];                                    \
      }                                                                 \
   }                                                                    \
}

#define EXTRACT_INT_FROM_C_TO_SCILAB(num,ptr_scilab,mumpsptr,length1,length2,one)\
if(mumpsptr == 0){       							\
    CreateVar(nb_RHS+num,"i",&one,&one,&ptr_scilab);                            \
    *istk(ptr_scilab)=-9999;                                                    \
}else{                                                  			\
CreateVar(nb_RHS+num,"i",&length1,&length2,&ptr_scilab);			\
int l=length1*length2;  							\
for (i=0;i<l;i++){ 								\
    *istk(ptr_scilab+i)=(mumpsptr)[i];}                                         \
 }                                                                              \
LhsVar(num)=nb_RHS+num;

#define EXTRACT_DOUBLE_FROM_C_TO_SCILAB(num,it,ptr_scilab1,ptr_scilab2,mumpsptr,length1,length2,one)\
if(mumpsptr == 0){                            					\
CreateVar(nb_RHS+num,"d",&one,&one,&ptr_scilab1);                		\
*stk(ptr_scilab1)=-9999;                                         		\
}else{                                                          		\
CreateVar(nb_RHS+num,"d",&length1,&length2,&ptr_scilab1); 			\
int l=length1*length2;  							\
for (i=0;i<l;i++){                                      			\
   *stk(ptr_scilab1+i)=(mumpsptr)[i];}                       			\
}                                                           			\
LhsVar(num)=nb_RHS+num;

#if MUMPS_ARITH == MUMPS_ARITH_z

#define EXTRACT_CMPLX_FROM_SCILAB_TOPTR(it,ptr_scilab1,ptr_scilab2,mumpsptr,type,length)\
if(ptr_scilab1[0] != -9999){                                                    \
  free(mumpsptr);	 			         			\
  mumpsptr=(double2 *) malloc(length*sizeof(double2));      			\
  for(i=0;i<length;i++){                                                	\
    (mumpsptr[i]).r = ptr_scilab1[i];} 						\
  if(it==1){ 									\
    for(i=0;i<length;i++){							\
      (mumpsptr[i]).i = ptr_scilab2[i];}					\
  }else{									\
    for(i=0;i<length;i++){							\
      (mumpsptr[i]).i = 0.0;}							\
  }                                                                             \
}

#define EXTRACT_CMPLX_FROM_C_TO_SCILAB(num,it,ptr_scilab1,ptr_scilab2,mumpsptr,length1,length2,one)\
  if(it!=1){									\
    Scierror(999,"Internal error, the variable built must be complex.\n");}  	\
  if(mumpsptr == 0){                                                        	\
    CreateCVar(nb_RHS+num,"d",&it, &one,&one,&ptr_scilab1,&ptr_scilab2);        \
    *stk(ptr_scilab1) = -9999;  						\
    *stk(ptr_scilab2) = -9999;  						\
  }else{									\
    int l=length1*length2;							\
    CreateCVar(nb_RHS+num,"d",&it,&length1,&length2,&ptr_scilab1,&ptr_scilab2); \
    for(i=0;i<l;i++){                                                      	\
      *stk(ptr_scilab1+i) = (mumpsptr[i]).r; 					\
    }                                                                           \
    for(i=0;i<l;i++){								\
      *stk(ptr_scilab2+i) = (mumpsptr[i]).i;					\
    }										\
  }										\
  LhsVar(num)=nb_RHS+num;

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

  *dmumps_par = (DMUMPS_STRUC_C *) malloc(sizeof(DMUMPS_STRUC_C));
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
  (*dmumps_par)->irhs_ptr  = NULL;
  (*dmumps_par)->irhs_sparse  = NULL;
  (*dmumps_par)->rhs_sparse  = NULL;
  (*dmumps_par)->pivnul_list  = NULL;
  (*dmumps_par)->listvar_schur  = NULL;
  (*dmumps_par)->schur  = NULL;
  (*dmumps_par)->sym_perm  = NULL;
  (*dmumps_par)->uns_perm  = NULL;
}



 static int dmumpsc(char *fname){


  /* RhsVar parameters */
  int njob, mjob, ljob, mint, nint, lint, nsym, msym, lsym, nA, mA, nRHS, nREDRHS, mRHS,lRHS, liRHS;
  int mREDRHS,lREDRHS,liREDRHS;
  int nicntl, micntl, licntl, ncntl, mcntl, lcntl, nperm, mperm, lperm;
  int ncols, mcols, lcols, licols, nrows, mrows, lrows, lirows, ns_schu , ms_schu, ls_schu;
  int nv_schu, mv_schu, lv_schu, nschu, mschu, lschu;
  int type_rhs, mtype_rhs, ntype_rhs, ltype_rhs;

  /* LhsVar parameters */
  int linfog, lrinfog, lrhsout,lrhsouti, linstout, lschurout, lschurouti, ldef;
  int lpivnul_list, lmapp, lsymperm, lunsperm;
  int one=1, temp1=80, temp2=40, temp3, temp4;
  int it, itRHS, itREDRHS; /* parameter for real/complex types */

  int i,j,k1,k2, nb_in_row,netrue;
  int *ptr_int;
  double *ptr_double;
  double *ptr_scilab;
#if MUMPS_ARITH == MUMPS_ARITH_z
  double * ptri_scilab;
#endif

  /* Temporary length variables */
  int len1, len2;
  /* Temporary pointers in stack */
  int stkptr, stkptri;

  /* C pointer for input parameters */
  int inst_address;
  int ne,inst;
  int *irn_in,*jcn_in;

  /* Variable for multiple and sparse RHS*/
  int posrhs, posschur, nz_RHS,col_ind,k;
  int *irhs_ptr;
  int *irhs_sparse;
  double *rhs_sparse;
#if MUMPS_ARITH == MUMPS_ARITH_z
  double *im_rhs_sparse;
  char * function_name="zmumpsc";
#else
  char * function_name="dmumpsc";
#endif

  SciSparse A;
  SciSparse RHS_SPARSE;
  DMUMPS_STRUC_C *dmumps_par;

  int dosolve=0;
  int donullspace=0;
  int doanal = 0;
  /* Check number of input parameters */
  CheckRhs(11,12);

  /* Get job value. njob/mjob are the dimensions of variable job. */
  GetRhsVar(2,"i",&mjob,&njob,&ljob);
  dosolve = (*istk(ljob) == 3 || *istk(ljob) == 5 ||*istk(ljob) == 6);
  doanal = (*istk(ljob) == 1 || *istk(ljob) == 4 || *istk(ljob) == 6);
  if(*istk(ljob) == -1){

    DMUMPS_alloc(&dmumps_par);
    GetRhsVar(1,"i",&msym,&nsym,&lsym);
    dmumps_par->sym=*istk(lsym);
    dmumps_par->job = -1;
    dmumps_par->par = 1;
    dmumps_c(dmumps_par);
    dmumps_par->nz = -1;
    dmumps_par->nz_alloc=-1;
    it=1;
  }else{
    /* Obtain pointer on instance */
    GetRhsVar(10,"i",&mint,&nint,&lint);
    inst_address=*istk(lint); /* EXTRACT_FROM_SCILAB_TOVAL(INST,inst_address); */
    ptr_int = (int *) inst_address;

    dmumps_par = (DMUMPS_STRUC_C *) ptr_int;
    if(*istk(ljob) == -2){
      dmumps_par->job = -2;
      dmumps_c(dmumps_par);
      DMUMPS_free(&dmumps_par);
    }else{
      /* Get the sparse matrix A */
      GetRhsVar(12,"s",&mA,&nA,&A);

      if (nA != mA || mA<1 ){
	Scierror(999,"%s: Bad dimensions for mat\n",function_name);
	return 0;
      }

      ne=A.nel;
      dmumps_par->n = nA;
	
      if(dmumps_par->sym != 0){
	netrue = (nA+ne)/2;
      }else{
	netrue = ne;
      }

      if(dmumps_par->nz_alloc < netrue ||dmumps_par->nz_alloc >= 2*netrue){
	MYFREE(dmumps_par->jcn);
	MYFREE(dmumps_par->irn);
	MYFREE(dmumps_par->a);
	
	dmumps_par->jcn = (int*)malloc(netrue*sizeof(int));
	dmumps_par->irn = (int*)malloc(netrue*sizeof(int));
	dmumps_par->a = (double2 *) malloc(netrue*sizeof(double2));
	dmumps_par->nz_alloc = netrue;
      }
      /* Check for symmetry in order to initialize only
       * lower triangle on entry to symmetric MUMPS code */
      if ((dmumps_par->sym)==0){
        /*
	 * Unsymmetric case:
	 * build irn from mnel for MUMPS format
	 * mA : number of rows
	 */

        if(doanal){
	  for(i=0;i<ne;i++){
	    (dmumps_par->jcn)[i]=(A.icol)[i];}
	  k1=0;
	  for (k2=1;k2<mA+1;k2++){
	    nb_in_row=0;
	    while(nb_in_row<(A.mnel)[k2-1]){
	      dmumps_par->irn[k1]=k2; /* matrix indices start at 1 */
	      k1=k1+1;
	      nb_in_row=nb_in_row+1;
	    }
	  }
        }

#if MUMPS_ARITH == MUMPS_ARITH_z
	for(i=0;i<ne;i++){
          ((dmumps_par->a)[i]).r = (A.R)[i];}
        if(A.it == 1){
          for(i=0;i<ne;i++){
            ((dmumps_par->a)[i]).i = (A.I)[i];}
        }else{
          for(i=0;i<ne;i++){
            ((dmumps_par->a)[i]).i = 0.0;}
        }
#else
	for(i=0;i<ne;i++){
          ((dmumps_par->a)[i]) = (A.R)[i];}
#endif
	dmumps_par->nz = ne;
	}
      else{
	/* symmetric case */
	k1=0;
        i=0;
	for (k2=1;k2<mA+1;k2++){
	  nb_in_row=0;
	  while(nb_in_row<(A.mnel)[k2-1]){
            if( k2 >= (A.icol)[i]){
	      if(k1>=netrue){	
	        Scierror(999,"%s: The matrix must be symmetric\n",function_name);
	  	return 0;
	      }
              (dmumps_par->jcn)[k1]=(A.icol)[i];
	      (dmumps_par->irn)[k1]=k2;
#if MUMPS_ARITH == MUMPS_ARITH_z
	      (dmumps_par->a)[k1].r=(A.R)[i];
              if(A.it == 1){
                ((dmumps_par->a)[k1]).i = (A.I)[i];}
              else{
                ((dmumps_par->a)[k1]).i = 0.0;}
#else
	      ((dmumps_par->a)[k1]) = (A.R)[i];
#endif
	      k1=k1+1;}
		
	      nb_in_row=nb_in_row+1;
	      i=i+1;
	     }
	  }
	dmumps_par->nz = k1;
      	}
	
        GetRhsVar(2,"i",&mjob,&njob,&ljob);
	dmumps_par->job=*istk(ljob);
	
	GetRhsVar(3,"i",&micntl,&nicntl,&licntl);
	EXTRACT_FROM_SCILAB_TOARR(istk(licntl),dmumps_par->icntl,int,60);

	GetRhsVar(4,"d",&mcntl,&ncntl,&lcntl);
	EXTRACT_FROM_SCILAB_TOARR(stk(lcntl),dmumps_par->cntl,double,15);
	
        GetRhsVar(5,"i",&mperm, &nperm, &lperm);
	EXTRACT_FROM_SCILAB_TOPTR(IT_NOT_USED,istk(lperm),istk(lperm),(dmumps_par->perm_in),int,nA);

	GetRhsCVar(6,"d",&it,&mcols,&ncols,&lcols,&licols);
        EXTRACT_FROM_SCILAB_TOPTR(it,stk(lcols),stk(licols),(dmumps_par->colsca),double2,nA);
	
        GetRhsCVar(7,"d",&it,&mrows,&nrows,&lrows,&lirows);
        EXTRACT_FROM_SCILAB_TOPTR(it,stk(lrows),stk(lirows),(dmumps_par->rowsca),double2,nA);


/*
 * To follow the "spirit" of the Matlab/Scilab interfaces, treat case of null
 * space separately. In that case, we initialize lrhs and nrhs automatically,
 * allocate the space needed, and do not rely on what is provided by the user
 * in component RHS, that is not touched.
 * At the moment the user should not call the solution step combined
 * with the factorization step when he/she sets icntl[25] to a non-zero value.
 * Hence we suppose infog[28-1] is available and we can use it.
 *
 * For users of scilab/matlab, it would still be nice to be able to set ICNTL(25)=-1,
 * and use JOB=6. If we want to make this functionality available, we should
 * call separately job=2 and job=3 even if job=5 or 6 and set nrhs (and allocate
 * space correctly) between job=2 and job=3 calls to MUMPS.
 *
 */
      if ( dmumps_par->icntl[25-1] == -1 && dmumps_par->infog[28-1] > 0) {
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
        nRHS=dmumps_par->nrhs;
        dmumps_par->lrhs=dmumps_par->n;
        dmumps_par->rhs=(double2 *)malloc((dmumps_par->n)*(dmumps_par->nrhs)*sizeof(double2));
        dmumps_par->icntl[19]=0;
         }

      else if(GetType(8)!=5){
/*        Dense RHS */
          GetRhsCVar(8,"d",&itRHS,&mRHS,&nRHS,&lRHS,&liRHS);

          if((!dosolve) || (stk(lRHS)[0]) == -9999){
          /* Could be dangerous ? See comment in Matlab interface */
            EXTRACT_CMPLX_FROM_SCILAB_TOPTR(itRHS,stk(lRHS),stk(liRHS),(dmumps_par->rhs),double2,one);
          }else{

            dmumps_par->nrhs = nRHS;
            dmumps_par->lrhs = mRHS;
            if(mRHS!=nA){
	      Scierror(999,"%s: Incompatible number of rows in RHS\n",function_name);
            }
            dmumps_par->icntl[19]=0;
            EXTRACT_CMPLX_FROM_SCILAB_TOPTR(itRHS,stk(lRHS),stk(liRHS),(dmumps_par->rhs),double2,(nRHS*mRHS));
          }
        }else{
/*        Sparse RHS */
          GetRhsVar(8,"s",&mRHS,&nRHS,&RHS_SPARSE);
          dmumps_par->icntl[19]=1;
          dmumps_par->nrhs = nRHS;
          dmumps_par->lrhs = mRHS;
          nz_RHS=RHS_SPARSE.nel;
          dmumps_par->nz_rhs=nz_RHS;

          irhs_ptr=(int*)malloc((nRHS+1)*sizeof(int));

          dmumps_par->irhs_ptr=(int*)malloc((nRHS+1)*sizeof(int));
          dmumps_par->irhs_sparse=(int*)malloc(nz_RHS*sizeof(int));
          dmumps_par->rhs_sparse=(double2*)malloc(nz_RHS*sizeof(double2));
          dmumps_par->rhs=(double2*)malloc((nRHS*mRHS)*sizeof(double2));
          /* transform row-oriented sparse multiple rhs (scilab)
	   * into column-oriented sparse multiple rhs (mumps) */
          k=0;
          for(i=0;i<nRHS+1;i++){
            irhs_ptr[i]=0;
            dmumps_par->irhs_ptr[i]=0;}
          for(i=1;i<mRHS+1;i++){
            for(j=0;j<(RHS_SPARSE.mnel)[i-1];j++){
              col_ind=(RHS_SPARSE.icol)[k];
              k++;
              ((dmumps_par->irhs_ptr)[col_ind])++;
            }
          }
          (dmumps_par->irhs_ptr)[0]=1;
          irhs_ptr[0]=(dmumps_par->irhs_ptr)[0];
          for(i=1;i<nRHS+1;i++){
            (dmumps_par->irhs_ptr)[i]=(dmumps_par->irhs_ptr)[i]+(dmumps_par->irhs_ptr)[i-1];
            irhs_ptr[i]= (dmumps_par->irhs_ptr)[i];
          }
          k=RHS_SPARSE.nel-1;
          for(i=mRHS;i>=1;i--){

            for(j=0;j<(RHS_SPARSE.mnel)[i-1];j++){
              col_ind=(RHS_SPARSE.icol)[k];
             (dmumps_par->irhs_sparse)[irhs_ptr[col_ind]-2]=i;
#if MUMPS_ARITH == MUMPS_ARITH_z
              if(RHS_SPARSE.it==1){
                ((dmumps_par->rhs_sparse)[irhs_ptr[col_ind]-2]).r=RHS_SPARSE.R[k];
                ((dmumps_par->rhs_sparse)[irhs_ptr[col_ind]-2]).i=RHS_SPARSE.I[k];
              }else{
                ((dmumps_par->rhs_sparse)[irhs_ptr[col_ind]-2]).r=RHS_SPARSE.R[k];
                ((dmumps_par->rhs_sparse)[irhs_ptr[col_ind]-2]).i=0.0;
              }
#else
              (dmumps_par->rhs_sparse)[irhs_ptr[col_ind]-2]=RHS_SPARSE.R[k];
#endif
              k--;
              irhs_ptr[col_ind]=irhs_ptr[col_ind]-1;
            }
          }
 	MYFREE(irhs_ptr);
	}
	
	GetRhsVar(9,"i",&nv_schu,&mv_schu,&lv_schu);
	dmumps_par-> size_schur=mv_schu;
	EXTRACT_FROM_SCILAB_TOPTR(IT_NOT_USED,istk(lv_schu),istk(lv_schu),(dmumps_par->listvar_schur),int,dmumps_par->size_schur);
	if(!dmumps_par->listvar_schur) dmumps_par->size_schur=0;

        if(dmumps_par->size_schur > 0){
	  MYFREE(dmumps_par->schur);	
          if(!(dmumps_par->schur=(double2 *)malloc((dmumps_par->size_schur*dmumps_par->size_schur)*sizeof(double2)))){
	    Scierror(999,"%s: malloc Schur failed in intmumpsc.c\n",function_name);
          }
          dmumps_par->icntl[18]=1;
        }else{
          dmumps_par->icntl[18]=0;
        }

        /* Reduced RHS */
        if ( dmumps_par->size_schur > 0 && dosolve ) {

          if ( dmumps_par->icntl[26-1] == 2 ) {
            /* REDRHS is on input */
            GetRhsCVar(11,"d",&itREDRHS,&mREDRHS,&nREDRHS,&lREDRHS,&liREDRHS);
            if (mREDRHS != dmumps_par->size_schur || nREDRHS != dmumps_par->nrhs ) {
             Scierror(999,"%s: bad dimensions for REDRHS\n");
            }
            /* Fill dmumps_par->redrhs */
            EXTRACT_CMPLX_FROM_SCILAB_TOPTR(itREDRHS,stk(lREDRHS),stk(liREDRHS),(dmumps_par->redrhs),double2,(nREDRHS*mREDRHS));
            dmumps_par->lrhs=mREDRHS;
          }

          if ( dmumps_par->icntl[26-1] == 1 ) {
            /* REDRHS on output. Must be allocated before the call */
            MYFREE(dmumps_par->redrhs);
            if(!(dmumps_par->redrhs=(double2 *)malloc((dmumps_par->size_schur*dmumps_par->nrhs)*sizeof(double2)))){
              Scierror(999,"%s: malloc redrhs failed in intmumpsc.c\n",function_name);
            }
          }
        }

        /* call C interface to MUMPS */
        dmumps_c(dmumps_par);

      }
    }

    if(*istk(ljob)==-2){
      return 0;
    }else{

      CheckLhs(11,11);

      EXTRACT_INT_FROM_C_TO_SCILAB(1,linfog,(dmumps_par->infog),one,temp1,one);

      EXTRACT_DOUBLE_FROM_C_TO_SCILAB(2,it,lrinfog,lrinfog,(dmumps_par->rinfog),one,temp2,one);

       if(dmumps_par->rhs && dosolve){ /* Just to know if solution step was called */
        it =1;
        EXTRACT_CMPLX_FROM_C_TO_SCILAB(3,it,lrhsout,lrhsouti,(dmumps_par->rhs),nA,nRHS,one);

      }else{
        it=1;
        EXTRACT_CMPLX_FROM_C_TO_SCILAB(3,it,lrhsout,lrhsouti,(dmumps_par->rhs),one,one,one);
      }

      ptr_int = (int *)dmumps_par;
      inst_address = (int) ptr_int;
      EXTRACT_INT_FROM_C_TO_SCILAB(4,linstout,&inst_address,one,one,one);


      temp4=dmumps_par->size_schur;
      if(temp4>0){
        it=1;
        EXTRACT_CMPLX_FROM_C_TO_SCILAB(5,it,lschurout,lschurouti,(dmumps_par->schur),temp4,temp4,one);
   }else{
        it=1;
        EXTRACT_CMPLX_FROM_C_TO_SCILAB(5,it,lschurout,lschurouti,(dmumps_par->schur),one,one,one);
      }

      /* REDRHS on output */
      it=1;
      if ( dmumps_par->icntl[26-1]==1 && dmumps_par->size_schur > 0 && dosolve ) {
        len1=dmumps_par->size_schur;
        len2=dmumps_par->nrhs;
      }
      else {
        len1=1;
        len2=1;
      }
      it=1;
      EXTRACT_CMPLX_FROM_C_TO_SCILAB(6,it,stkptr,stkptri,(dmumps_par->redrhs),len1,len2,one)


      MYFREE(dmumps_par->redrhs);
      MYFREE(dmumps_par->schur);
      MYFREE(dmumps_par->irhs_ptr);
      MYFREE(dmumps_par->irhs_sparse);
      MYFREE(dmumps_par->rhs_sparse);
      MYFREE(dmumps_par->rhs);

      /* temp3=dmumps_par->deficiency;*/
      temp3=dmumps_par->infog[27];
      EXTRACT_INT_FROM_C_TO_SCILAB(7,lpivnul_list,(dmumps_par->pivnul_list),one,temp3,one);

      EXTRACT_INT_FROM_C_TO_SCILAB(8,lsymperm,(dmumps_par->sym_perm),one,nA,one);

      EXTRACT_INT_FROM_C_TO_SCILAB(9,lunsperm,(dmumps_par->uns_perm),one,nA,one);

      nicntl=60;
      EXTRACT_INT_FROM_C_TO_SCILAB(10,licntl,(dmumps_par->icntl),one,nicntl,one);
      ncntl=15;
      EXTRACT_DOUBLE_FROM_C_TO_SCILAB(11,it,lcntl,lcntl,(dmumps_par->cntl),one,ncntl,one);
      return 0;

    }
}


static GenericTable Tab[]={
#if MUMPS_ARITH == MUMPS_ARITH_z
{(Myinterfun) sci_gateway, dmumpsc,"zmumpsc"}
#else
{(Myinterfun) sci_gateway, dmumpsc,"dmumpsc"}
#endif
};

#if MUMPS_ARITH == MUMPS_ARITH_z
int C2F(scizmumps)()
#else
int C2F(scidmumps)()
#endif
{Rhs = Max(0, Rhs);
(*(Tab[Fin-1].f))(Tab[Fin-1].name,Tab[Fin-1].F);
return 0;
}
