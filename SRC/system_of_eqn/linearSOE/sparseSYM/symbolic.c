/* 
 * This part is for the ordering, post-ordering and symoblic
 * factorization.
 * ----------------------
 * The subrotines in this file form the adjacency structure of the matrix
 * The first subroutines also calls the subroutines to form the
 * elimination tree and symbolic factgorization

 * This version is made to work with the lm array of a finite element problem
 * where each node has multiple dof's.
 * ----------------------------------
 * Originally written by:  David R. Mackay
 *
 * Modified by:
 *  Jun Peng (junpeng@stanford.edu)
 *  Prof. Kincho H. Law
 *  Stanford University
 * --------------------
 */


#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <assert.h>
#include "utility.h"
#include "FeStructs.h"
/*#include "globalVars.h"*/


#ifdef _WIN32
extern int MYGENMMD(int *neq, int *fxadj, int *adjncy, int *winvp,
			     int *wperm, int *delta, int *fchild, int *parent,
			     int *sibling, int *marker, int *maxint, int *nofsub,
			     int *kdx);

#else
extern int mygenmmd_(int *neq, int *fxadj, int *adjncy, int *winvp,
		     int *wperm, int *delta, int *fchild, int *parent,
		     int *sibling, int *marker, int *maxint, int *nofsub,
		     int *kdx);
#endif

void gennd(int neqns, int **padj, int *mask, int *perm, 
	   int *xls, int *ls, int *work);
void forminv(int neqns, int *perm, int *invp);
int pfordr(int neqns, int **padj, int *perm, int *invp, int *parent, int *fchild, 
	   int *sibling, int *winvp, int *wperm, int *list, int *rowblks);
void genrcm(int neqns, int **padj, int *perm, int *mask, int *xls, int *work);
void pfblk (int nblks, int *xblk, int *list);

int nodfac(int *perm, int *invp, int **padj, int *ancstr , int *list, int neqns, 
	   int nblks, int *xblk, int *envlen, OFFDBLK **segfirst, 
	   OFFDBLK **first, int *rowblks );
int setenvlpe(int neqns, double **penv, int *envlen);



/* int symFactorization(int *fxadj, int *adjncy, int neq, int LSPARSE) */


int symFactorization(int *fxadj, int *adjncy, int neq, int LSPARSE, 
		     int **xblkMY,
		     int **invpMY, int **rowblksMY, OFFDBLK ***begblkMY,
		     OFFDBLK **firstMY, double ***penvMY, double **diagMY)

{
    int delta, maxint;
    int nofsub, kdx;
    int ndnz;
    int *marker;
    int *winvp, *wperm;
    int i;
    int *perm, *parent, *fchild, *sibling;
    int **padj;

    int nblks;
    int *xblk;
    int *invp;
    int *rowblks;
    OFFDBLK **begblk;
    OFFDBLK *first;
    double **penv;
    double *diag;


 /* set up storage space and pointers */ 

    perm = (int *)calloc(neq +1   , sizeof(int)) ;
    invp = (int *)calloc(neq +1   , sizeof(int)) ;
    parent = (int *)calloc(neq +1 , sizeof(int)) ;
    fchild = (int *)calloc(neq +1 , sizeof(int)) ;
    sibling = (int *)calloc(neq +1, sizeof(int)) ;
    marker = (int *) calloc(neq +1, sizeof(int)) ;
    winvp  = (int *) calloc(neq +1, sizeof(int)) ;
    wperm  = (int *) calloc(neq +1, sizeof(int)) ;
    assert( perm && invp && parent && fchild && sibling && marker
	    && winvp && wperm != NULL) ;

    kdx = 0;
    delta = 1;
    maxint = 99999999;
    nofsub = 99999999;

 /* Using (fxadj, adjncy) pair to form the padj  */

    for(i=0; i<=neq; i++) {
        fxadj[i]++;
    }
    padj = (int **)calloc(neq+1,sizeof(int *)) ;
    assert(padj != NULL) ;
    padj[0] = (int *)calloc(fxadj[neq]+1, sizeof(int)) ;
    assert(padj[0] != NULL) ;
    copyi(fxadj[neq], adjncy, padj[0]);
    for (i=1; i<=neq; i++)
       padj[i] = padj[0] + fxadj[i] - 1;
    for (i=0; i<fxadj[neq]-1; i++)
       adjncy[i]++ ;

 /* Choose different ordering schema */

    switch(LSPARSE)
    {
       case 1:
   /* Now call minimum degree ordering  ( a fortran subroutine) */
#ifdef WIN32 
	 MYGENMMD( &neq, fxadj, adjncy, winvp, wperm, &delta, fchild, parent,
		   sibling, marker, &maxint, &nofsub, &kdx ) ;
#else
	 mygenmmd_( &neq, fxadj, adjncy, winvp, wperm, &delta, fchild, parent,
		    sibling, marker, &maxint, &nofsub, &kdx ) ;
#endif
         /* reset subscripts for c rather than fortran */
         for (i=0;i<=neq;i++)
         {
            winvp[i]-- ;
            wperm[i]-- ;
         }
         break ;

      case 2:
	/* Now call the nested dissection ordering */

         gennd(neq,padj,marker,wperm,fchild,sibling,parent) ;
         forminv(neq,wperm,winvp) ;
         break ; 

      case 3:
	/* Now call the general reverse chuthill-mckee ordering */

         genrcm(neq, padj, wperm, marker, fchild, sibling ) ;
         forminv(neq,wperm, winvp) ;
         break ;
   }

   /* free up space used just for mygenmmd and the fortran program */
   /*
    free(fxadj);
    free(adjncy);
   */

   rowblks = (int *)calloc(neq+1,sizeof(int)) ;
   assert(rowblks != 0) ;

/* set up the elimination tree, perform postordering           */
   if (LSPARSE < 4) {
       nblks = pfordr( neq, padj, perm, invp, parent, fchild, sibling,
		       winvp, wperm, marker, rowblks ) ;
   } 
   else { 
      for (i=0;i<=neq;i++)
      { 
	 invp[i] = i ;   
	 perm[i] = i ;
	 parent[i] = neq ;
	 rowblks[i] = 0 ;
      }
      marker[0] = 0 ;
      marker[1] = neq ;
      nblks = 1 ;
   }
         
   free(winvp) ;
   free(wperm) ;
   free(sibling) ;

/*  set up xblk profile blocks  and space for numerical values */
   xblk = (int *)calloc(nblks+1, sizeof(int)) ;
   begblk = (OFFDBLK **)  calloc(nblks+1, sizeof(OFFDBLK *)) ;
   assert(xblk && begblk != NULL) ;
         
/* set up xblk index: the beginning row/column of each block */      
        
   pfblk( nblks, xblk, marker );
        
/*       -------------------------------------------------
         perform the symbolic factorization and obtain the
         number of nonzeros
         -------------------------------------------------
*/  
           
   nodfac(perm, invp, padj, parent, fchild , neq, nblks,
	  xblk, marker, begblk, &first, rowblks) ;

   free(perm) ;
   free(parent) ;
   free(fchild) ;
   free(padj[0]) ;
   free(padj);

   penv = (double **)calloc(neq + 1, sizeof(double *)) ;
   diag = (double *)calloc(neq + 1,sizeof(double )) ;
   assert ( penv && diag != NULL) ;
   ndnz = setenvlpe(neq, penv, marker) ;
        
   free(marker);

   *xblkMY = xblk;
   *invpMY = invp;
   *rowblksMY = rowblks;
   *begblkMY = begblk;
   *firstMY = first;
   *penvMY = penv;
   *diagMY = diag;


   for(i=0; i<=neq; i++) {
       fxadj[i]--;
   }
   for (i=0; i<fxadj[neq]; i++) {
       adjncy[i]-- ;
   }

   return(nblks);
}





