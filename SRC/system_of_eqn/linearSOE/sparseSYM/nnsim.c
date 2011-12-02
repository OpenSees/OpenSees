/*
 * File:  nnsim.c
 * ==============
 * Subroutines in this program perform the symbolic factorization
 * and set up the data structure and memory allocation 
 * for the generalized profile matrix solver
 *
 * Altered to avoid eliminate the realloc calls in setting up pplnz 
 * March 8 1990 David Mackay
 */


#include <stdio.h>
#include <assert.h>
#include <malloc.h>
#include "FeStructs.h"

#define ABS(x)     (((x) < 0) ? -(x) : (x))


/***********************************************************
************************nodfac ******************************
*************************************************************

       Purpose:  to construct structure of lower triangular
                 factor by rows.  It is based on teh node
                 addition model with ordered tree.

       Input:
          (perm, invp) - permutation vector and its inverse
          padj         - adjancency structure of original ordering
          ancstr       - ancstr of each node in teh ordered tree

       Output:
          pnzbeg, pnzsub - structure of matrix factor l
          nonz           - number of nonzeros in l

       Working variables:
          list           - temporary list of subscripts in a row

*****************************************************************/

int nodfac(int *perm, int *invp, int **padj, int *ancstr , int *list, int neqns, 
	   int nblks, int *xblk, int *envlen, OFFDBLK **segfirst, 
	   OFFDBLK **first, int *rowblks )
{ 
   int i, node, nbr, qm, m, nnext ;
   int bcount, knz, cnz ;
   int  nbrblk ;
   int *pt ;
   int *len ;
   int count ;
   OFFDBLK **segprv ;
   OFFDBLK *p, *po, *nbeg ;
   OFFDBLK junk ;

   cnz &= 0 ;
   p = NULL ;
   *first = NULL ;
   po = &junk ;
   count = 0 ;
   bcount = 0 ;
   segprv = (OFFDBLK **) calloc((nblks+1),sizeof(OFFDBLK *)) ;
   len = (int *) calloc(nblks, sizeof(int)) ;
   assert (segprv && len != NULL) ;
   for (i=0;i<=nblks;i++)
   { 
      segfirst[i] = NULL ;
      segprv[i]  =  NULL ;
   }
   zeroi(nblks, len) ;
   /* initialize link */
   for (i = 0; i< neqns ; i++)
      list[i] = i ;
   zeroi(neqns,envlen) ;
   for (node = 1 ; node < neqns; node++)
   {
      knz = 0 ;
      i = perm[node] ;
/*    sort adjacency list */
      for (pt = padj[i] ; pt < padj[i+1] ; pt++)
      {  
	 nbr = invp[*pt] ;
         if ( nbr >= node) continue ;
         qm = node ;
         do
         {  m = qm ;
            qm = list[m] ;
         }  while(qm<=nbr) ;
         list[m] = nbr ;
         list[nbr] = qm ;
      }
      nbr = list[node] ;
      list[node] = node ;
      nbeg = NULL ;
      while (ancstr[nbr] <= node)
      {  
	 p = (OFFDBLK *)malloc( sizeof(OFFDBLK));
	 assert (p != NULL) ;
         p->row = node ;
         p->beg = nbr ;
	 po->next = p ;
	 po = p ;
         nbrblk = rowblks[nbr] ;
         knz += (xblk[nbrblk+1] - nbr) ;
         len[count - bcount] = xblk[nbrblk+1] - nbr ;
         count++;
	 if (*first == NULL) *first = p ;
	 if (nbeg == NULL) nbeg = p ;
	 if (segprv[nbrblk] != NULL) segprv[nbrblk]->bnext = p ;
	 segprv[nbrblk] = p ;
	 if (segfirst[nbrblk] == NULL) segfirst[nbrblk] = p ;
         qm = nbr ;
         do
         {  nnext = list[qm] ;
            list[qm] = qm ;
            qm = nnext ;
	 } while (nnext < xblk[nbrblk+1] ) ;
	 nbr = ancstr[nbr] ;
	 if (nbr >= nnext ) nbr = nnext ;
	 else list[nbr] = nnext ;
      }
      /* part of the diagonal envelop block */
      envlen[node] = node - nbr ;
/*    should now allocate space for row  and set up pointers */
      if (knz > 0) 
      {  
	 nbeg->nz = (double *)calloc(knz, sizeof(double)) ;
         assert(nbeg->nz != NULL) ;
         if ( bcount < count) bcount++ ;
         m = bcount ;
         while (bcount < count)
         {  (nbeg->next)->nz = nbeg->nz + len[bcount - m] ;
            nbeg = nbeg->next ;
            bcount++ ;
         }
      }
      cnz += knz ;   
   }  /* end for node */
   printf("nozeros in row segments %d\n",cnz) ;

/* -----------------------------------------------
   add on ending pieces for loops in factorization
   and backsolving to catch on 
   -----------------------------------------------*/
   node = neqns ;
   nbr = neqns ;
   p = (OFFDBLK *)calloc(1, sizeof(OFFDBLK));
   assert (p != NULL) ;
   p->row = neqns ;
   p->beg = neqns ;
   po->next = p ;
   p->next = p ;
   p->bnext = p ;
   for (i=0 ; i<=nblks ; i++)
   {
      if (segfirst[i] == NULL) segfirst[i] = p ;
      else segprv[i]->bnext = p ;
   }
   if (*first == NULL) *first = p ;

   free(len) ;
   free(segprv) ;
   return ;
}



/************************************************************************
 ************  setenv ..... set up envelope   ***************************
 ************************************************************************
 
    purpose - allocate space for the envelope structure
	      and setup pointers

    input parameters -
        neqns - no of equations
	penv - array of pointers to be filled
	envlen - an array with the lengths of each row.

    output parameters -
        penv - filled array of pointers pointing to allocated space
        knz  - number of nonzeros in the envelope structure .

 
 ************************************************************************/
     
int setenv(int neqns, double **penv, int *envlen)
{
   int i, knz ;

   knz = 0 ;
   for (i=1;i<neqns;i++)
      knz += envlen[i] ;
   
   penv[0] = (double *)calloc(knz+1,sizeof(double)) ;
   assert(penv[0] != NULL ) ;
   for (i=0;i<neqns;i++)
   {  penv[i+1] = penv[i] + envlen[i] ;
   }

   return(knz);
}
