/*
 * File:  nmat.c
 * =============
 *
 * Changed by Jun Peng (junpeng@leland.stanford.edu)
 * on April 2000.
 * altered to improve data access. Use division instead of sqrt.
 */


#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include "FeStructs.h"
#include "utility.h"
#include "nmat.h"


#define MAX(x,y)   (((x) < (y)) ? (y) : (x))

/***************************************************************
 ******    pfsfct ..... general envelope symmetric fact   ******
 ***************************************************************
 
   purpose - this routine performs sparse cholesky factorization
        on a matrix stored in eneralized envelope format.
 
   input parameters -
        neqns - number of equations.
        penv  - index vector for block diagonal envelope
        plnz  - index vector for lnz.  plnz(i) points to the
                start of nonzeros in column i of factor l.
        (nblks,xblk) - partitioning blocks
   updated parameters -
        diag  - the diagonal of l overwrites that of a.
        env   - block diag envelope
        lnz   - on input, contains nonzeros of a, and on
                return, the nonzeros of l.
        iflag - the error flag.  it is set to 1 if a zero
                or negative square root occurs during the
                factorization.
   working parameters -
        nxtlnz - temporary vector pointing to the location of
                 the next segment under each block
        nxtsub - temporary vector pointin to the location of
                 the next row subscript under each block
        temp   - temporary vector to accumulate row modifications.

   program subroutines 
        pfefct, pflslv 
 
 ***************************************************************/
int pfsfct(int neqns, double *diag, double **penv, int nblks, 
	   int *xblk, OFFDBLK **begblk, OFFDBLK *first, int *rowblks)
/*************************************************************** 
 ***************************************************************/
{  
   int blk, nextblk, jbeg, iflag ;
   int iband, blkbeg, blkend, blksze ;
   int jrow, krow ;
   int jblk, jb, kb, pos ;
   OFFDBLK *marker, *ks, *js, *ls ;
   double *work;
   int ii;
   
   if  ( neqns <= 0 )  return(0) ;
   marker = begblk[0];

   marker = first;
   js = marker;
/* ----------------------------------------------------------
   for each block blk, do ...
   ----------------------------------------------------------*/
   for (blk = 0; blk < nblks; blk++)
   {  
      nextblk = blk + 1 ;
      blkbeg = xblk[blk] ;
      blkend = xblk[nextblk]  ;
      blksze = blkend - blkbeg ;
/*    --------------------------------------------------------
      update rows from row segments
      The function Dotrows();
      -------------------------------------------------------*/
      while( js->row < blkend)
      {
	 jrow = js->row;
	 jbeg = js->beg;
 
	 jblk = rowblks[jbeg];
         ls = begblk[blk] ;
         ks = js->bnext ;
/*    -------------------------------------------------------
      update the diagonals from the off diagonal row segments
      ------------------------------------------------------*/
	 
	 iband = xblk[jblk+1] - jbeg;
	 work = (double*) calloc(iband, sizeof(double)); 
	 for (ii = 0; ii < iband; ii++) {
	     work[ii] = js->nz[ii];
	     js->nz[ii] /= diag[ii + jbeg]; 	    
	 }
	 diag[jrow] -= dot_real(js->nz, work, iband);
	 if (diag[jrow] == 0) return (1);
	 free (work);
	 
	 if (ks->row < blkend )
	 {  /* part of envelop block*/
	    for ( ; ks->row < blkend ; ks = ks->bnext)
	    {
	       krow = ks->row ;
	       pos = MAX(jbeg, ks->beg) ;
	       iband = xblk[jblk+1] - pos;
	       jb = pos - jbeg ;
	       kb = pos - ks->beg ;
	       pos = jrow - krow + (penv[krow + 1] - penv[krow]) ;
	       *(penv[krow] + pos) -= 
		  dot_real(js->nz+jb, ks->nz+kb, iband);
            }
         }
	 for ( ; ks->beg < blkend ; ks = ks->bnext)
	 {
	    krow = ks->row ;
	    pos = MAX(jbeg, ks->beg);
	    iband = xblk[jblk+1] - pos;
	    jb = pos - jbeg ;
	    kb = pos - ks->beg ;
	    /* part of another row segment */
	    while ( ls->row != krow) ls = ls->bnext ;
	    pos = jrow - ls->beg ;
	    ls->nz[pos] -= 
	        dot_real(js->nz+jb, ks->nz+kb, iband) ;
         }

	 js = js->next ;
      }
/*    -------------------------------------------------------
      perform envelope fct on diag block blk.
      -------------------------------------------------------
*/
      iflag = pfefct(blksze, penv+blkbeg, diag+ blkbeg) ;
      if (iflag) return(nextblk);

/*    -------------------------------------------------------
      for each row "node" in this block, do
         update row segments under block blk with a backsolve
      -------------------------------------------------------
*/
      for (ks = begblk[blk]; ks->beg < blkend ; ks = ks->bnext )
      {  jbeg = ks->beg ;
         iband = blkend - jbeg ;
         pflslv(iband, (penv + jbeg), (diag + jbeg), ks->nz);
      }
   }

   return(0) ;
}

/***************************************************************
 ****     pfefct ..... envelope block diagonal env fact      ***
 ***************************************************************
 
   purpose - this routine factors a positive definite
         matrix a into l*l(transpose).  the matrix a is stored
         in the envelope format.  the algorithm used is the
         standard bordering method.
 
   input parameters -
        neqns - number of equations.
        penv - the envelope index vector.
   updated parameters -
        penv - the envelope of l overwrites that of a.
        diag - the diagonal of l overwrites that of a.
        iflag - the error flag.  it is set to 1 if a zero or
                negative square root is detected during the
                factorization.
   program subroutines
      pflslv, dot_real
 
 ***************************************************************/

int pfefct(int neqns, double **penv, double *diag)
{  
   double *ptenv ; 
   int iband, i, jj, ifirst ;
   double *work;
   
/*    -------------------------------------------------
      for each row i, ...
      -------------------------------------------------
*/
   for (i=1; i < neqns ; i++)
   {  
      ptenv = penv[i] ;
      iband = penv[i+1] - ptenv ;
      work = (double *)calloc(iband, sizeof(double));

      if ( iband > 0 )
      {  
	 ifirst = i - iband ;
/*       ---------------------------------------
          compute row i of the triangular factor.
         --------------------------------------- */

         pflslv( iband, penv+ifirst, diag+ifirst, ptenv );
	 for (jj = 0; jj < iband; jj++) {
	     work[jj] = ptenv[jj];
	     ptenv[jj] = ptenv[jj] / diag[i+jj-iband]; 
	 }
	 
         diag[i] -= dot_real(ptenv, work, iband ) ;
      }
      
      free (work);

      if ( fabs(diag[i]) < 1.0e-16)  return (1);  
   }

   return(0) ;
}

/***************************************************************
 *****  pfsslv ..... generalized envelope symmetric solve  *****
 ***************************************************************
 
   purpose - to solve a symmetric system usin gen envelope.
 
   input parameters -
        neqns  - number of equations.
        dia   - diagonal components of l.
        (penv,env) - block diag envelope of l.
        (plnz,lnz) - structure of nonzeros in l.
        (nblks, xblk) - partitioning blocks.
   updated parameters -
        rhs - on input, it contains the rhs vector, and on
                output, the solution vector.

   program subroutines
      pflslv, pfuslv
 
 ***************************************************************/
 
void pfsslv(int neqns, double *diag, double **penv, int nblks, 
	    int *xblk, double *rhs, OFFDBLK **begblk)

/***************************************************************
 
 ***************************************************************/
{  int ii, j, irow, blk ;
   int nextblk, blkbeg, blkend, blksze ;
   OFFDBLK *is ;
   double *ptr ;
 
   if  ( neqns <= 0 )  return ;

/*       --------------------------------------------
         forward substitution: for each block, do ...
         --------------------------------------------
*/
   for (blk = 0; blk < nblks ; blk++)
   {  nextblk = blk + 1 ;
      blkbeg = xblk[blk] ;
      blkend = xblk[nextblk] ;
      blksze = blkend - blkbeg ;

/*         ----------------------------------------------------------
           lower solve using diag block,
           ----------------------------------------------------------
*/
      pflslv( blksze, penv+blkbeg, diag+blkbeg, rhs+blkbeg ) ;

/*    ----------------------------------------------------------
      then update the rhs.
      ----------------------------------------------------------
*/
      for (is = begblk[blk] ; is->beg < blkend ; is = is->bnext ) 
      {  ptr = is->nz ;
         j = is->beg ;
         irow = is->row ;
         rhs[irow] -= dot_real(ptr, (rhs+j), (blkend-j)) ;
      }
   }
/* --------------------------------------------------------------
   backward substitution: for each block in reverse order, do ...
   --------------------------------------------------------------
*/
   for (blk = nblks-1; blk >= 0; blk--)
   {  
      nextblk = blk + 1 ;
      blkbeg = xblk[blk] ;
      blkend = xblk[nextblk] ;
      blksze = blkend - blkbeg ;
/*    ----------------------------------------------------------
      update the rhs,
      ----------------------------------------------------------
*/
      for (ii = blkbeg; ii < blkend; ii++) {
	 rhs[ii] /= diag[ii];
      }
      
      for (is=begblk[blk] ; is->beg < blkend ; is = is->bnext ) 
      {  ptr = is->nz ;
         j = is->beg ;
         irow = is->row ;
	 saxpy((rhs+j), ptr, -rhs[irow], (blkend-is->beg)) ;
      }
/*         ----------------------------------------------------------
           then upper solve using diag block.
           ----------------------------------------------------------
*/
      pfuslv( blksze, penv+blkbeg, diag+blkbeg, rhs+blkbeg ) ;
   }

   return ;
}

/***************************************************************
 **********     pflslv ..... gen envelope lower solve     ******
 ***************************************************************
 
   purpose - this routine solves a lower triangular system
        l x = rhs, where l is stored in the envelope format.
 
   input parameters -
        neqns  - number of equations.
        (penv) - array pair for the envelope of l.
        diag   - array for the diagonal of l.
   updated parameters -
        rhs    - on input, it contains the right hand vector.
                on return, it contains the solution vector.
 
 ***************************************************************/
 
void pflslv (int neqns, double **penv, double *diag, double *rhs)
{ 
  int i, iband ;

   if ( neqns <= 1 )  return ;
/*      -----------------------------------------------
        for each row i, do ...
        -----------------------------------------------
*/
   for (i = 1; i < neqns; i++)
   {  
      iband = penv[i+1] - penv[i] ;
      if (iband > i) iband = i ;
      if (iband > 0)
      {
         rhs[i] -= dot_real(penv[i+1] - iband, rhs+i-iband, iband) ;
      }
   }

   return ;
}

/***************************************************************
 **********     pfuslv ..... envelope upper solve     ******
 ***************************************************************
 
   purpose - this routine solves an upper triangular system
         u x = rhs, where u is stored in the envelope format.
 
   input parameters -
        neqns  - number of equations.
        (penv,env) - array pair for the envelope of u.
        diag   - array for the diagonal of u.
   updated parameters -
        rhs    - on input, it contains the right hand side.
                on output, it contains the solution vector.
 
 ***************************************************************/
 
void pfuslv(int neqns, double **penv, double *diag, double *rhs)
/***************************************************************/
{  int i, k ;
   double s, *ptr ;

   for (i=neqns-1; i >= 0 ;i--)
   {
/*    --------------------------------------------------
      for each row "i" in reverse order, do ...
      --------------------------------------------------
*/
      if  ( rhs[i] == 0.0e0 )  continue ;
/*    ------------------------------------------------
      compute the solution rhs(i) and
      ------------------------------------------------
*/
      /*      rhs[i] /= diag[i] ; */
      s = rhs[i] ;
/*    ------------------------------------------------
      use it to update the rhs entries.
      iband = penv[i+1] - penv[i] ;
      ------------------------------------------------
*/
      k = i-1 ;
      for (ptr = penv[i+1]-1; ptr >= penv[i] ; ptr--,k--)
         rhs[k] -= ( *ptr * s) ;
/*    or saxpy operation */
   }
   return ;
}
