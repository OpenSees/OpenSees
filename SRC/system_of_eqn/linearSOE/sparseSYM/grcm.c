/*
 * File:  grcm.c
 * =============
 * Do the general reverse Cuthill-Mckee ordering
 *
 * Originally written by:  David R. Mackay
 *
 * Modified by:
 *  Jun Peng (junpeng@stanford.edu)
 *  Prof. Kincho H. Law
 *  Stanford University
 * --------------------
 */

#include <assert.h>

/* fmk - adding prototypes */
void zeroi(int, int *);
int rcm(int root, int **padj, int *mask, int *perm, int *deg, int *work);
int fnroot(int root, int **padj, int *mask, int *nlvl, int *xls, int *ls);
int ndegree(int root, int **padj, int *mask, int *deg, int *ls, int *work);
void revrse(int n, int *v);

/*************************************************************************
***********************genrcm . . . general reverse cuthill mckee ********
**************************************************************************

  purpose - genrcm finds the reverse cuthill-mckee
            ordering for a general graph. for each connected
            component in the graph genrcm obtains the ordering
            by calling the subroutine rcm.

  input parameters -
            neqns - number of equations
            padj - the adjacency structure
  output parameter -
            perm - vector that contains the rcm ordering
  working parameters
            mask - is used to mark variables that have been
                   numbered during the ordering process.  it is
                   initialized to 0 and set to -1 as each nod
		   is numbered.
            xls - the index vector for a level structure. the level
                  structure is sorted in the currently unused spaces
		  in the permutation vector perm.
            work - a working vector
  program routines -
            fnroot, rcm
*************************************************************************/

void genrcm(int neqns, int **padj, int *perm, int *mask, int *xls, int *work)
{
   int num, i, root, nlvl, ccsize;
   zeroi(neqns, work);
   zeroi(neqns, mask);

   num = 0;
   for (i=0;i<neqns ; i++)
   {
/*      ---------------------------------------------------
        for each masked connected component
        -------------------------------------------------*/
        if (mask[i] < 0) continue ;
        root = i ;
/*    	------------------------------------------------------------
	first find a pseudo-peripheral node root.  note that the level
	structure found by fnroot is stored starting at perm[num].
	then rcm is called to order the component using root as the
	starting node.
	----------------------------------------------------------*/
        root = fnroot(root, padj, mask, &nlvl, xls, perm + num) ;
	ccsize = rcm(root, padj, mask, perm+num,xls, work) ;
	num += ccsize ;
	if (num > neqns) return ;
   }
   return;
}

/*************************************************************************
 *****************subrcm . . . substructure reverse cuthill mckee ********
 *************************************************************************

  purpose - subrcm finds the reverse cuthill-mckee
            ordering for a subgraph. for each connected
	    component in the graph subrcm obtains the ordering
	    by calling the subroutine rcm.

  input parameters -
            neqns - number of equations
	    padj - the adjacency structure
	    root - first trial root ( any node that lies
			int the subgraph)
  output parameter -
            perm - vector that contains the rcm ordering
  working parameters
            mask - is used to mark variables that have been
                   numbered during the ordering process.  it is
		   initialized to 0 and set to -1 as each node
		   is numbered.
            xls - the index vector for a level structure. the level
                  structure is sorted in the currently unused spaces
		  in the permutation vector perm.
            work - a working vector
  program routines -
            fnroot, rcm
*************************************************************************/
void subrcm (int neqns, int root, int **padj, int *perm,
	     int *mask, int *xls, int *work)
{
   int num, nlvl, ccsize ;
   zeroi(neqns, work) ;

   num = 0 ;
   if (mask[root] <0) return ;
/* ------------------------------------------------------------
   first find a pseudo-peripheral node root.  note that the level
   structure found by fnroot is stored starting at perm[num].
   then rcm is called to order the component using root as the
   starting node.
   ----------------------------------------------------------*/
   root = fnroot(root, padj, mask, &nlvl, xls, perm + num) ;
   ccsize = rcm(root, padj, mask, perm+num,xls, work) ;

   num += ccsize ;
   if (num > neqns) return ;
   return ;
}


/***********************************************************************
********************  rcm . . . reverse cuthill mckee ******************
************************************************************************

  purpose - rcm numbers a connected component specified by
            mask and root, using the rcm algorithm.  the
	    numbering is to be started at the node root.
  input parameters -
            root - is the node that defines the  connected
	    component and it is used as the starting
	    node for the rcm ordering.
	    padj - the adjacency structure

  updated parameters -
            mask - only those nodes with nonnegative input mask
	    values are considered by the routine.
	    the nodes numbered by rcm will have their
	    mask values set to zero.

  output parameters -
            perm - will contain the rcm ordering.
	    ccsize - is the size of the connected component that
	    ahas been numbered by rcm.

  working parameter -
            work - must be set to zero's before calling!!!!!!!!!
	    deg - is a temporary vector used to hold the degree
	    of the nodes in the section graph specified
	    by mask and root.
  program routines
            degree
*************************************************************************/
int rcm(int root, int **padj, int *mask, int *perm, int *deg, int *work)
{
   int ccsize ;
   int i, lbegin, lvlend, lnbr, nbr, node, fnbr, k, l, lperm ;
   int *ptr ;

/* ---------------------------------------
   find the degrees of the nodes in the
   component specified by mask and root.
   --------------------------------------*/
   ccsize = ndegree(root, padj, mask, deg, perm, work) ;
   mask[root] |= -1 ;
   if (ccsize <= 1) return(ccsize) ;
   lvlend &= 0 ;
   lnbr = 1 ;
/* -----------------------------------------------------
   lbegin and lvlend point to the beginning and
   the end of the current level respectively.
   -----------------------------------------------------*/
   do
   {
      lbegin = lvlend ;
      lvlend = lnbr ;
      for (i = lbegin ; i < lvlend ; i++)
      {
/*       ------------------------------------------
         for each node in current level . . .
         ----------------------------------------*/
         node = perm[i] ;
/*       -----------------------------------------------------
         find the unnumbered neighbors of node.
         fnbr and lnbr point to the first and last
         unnumbered neighbors respectively of the current
         node in perm.
         -----------------------------------------------------*/
         fnbr = lnbr ;
         for (ptr = padj[node] ; ptr < padj[node + 1] ; ptr++)
         {  if (mask[*ptr] < 0 ) continue ;
            mask[*ptr] |= -1 ;
            perm[lnbr] = *ptr ;
            lnbr++ ;
         }
         if (fnbr < lnbr-1)
         {
/*          ---------------------------------------------------
            sort the neighbors of node in increasing
            order by degree. linear insertion i sused.
            --------------------------------------------------*/
            k = fnbr ;
            do
            {  l = k ;
               k++ ;
               nbr = perm[k] ;
               while( l >= fnbr )
               {  lperm = perm[l] ;
                  if (deg[lperm] <= deg[nbr]) break ;
                  perm[l+1] = lperm ;
                  l-- ;
                }
                perm[l+1] = nbr ;
            } while(k < lnbr-1) ;
         } /* endif */
      }  /* end for i =  */
   }  while(lnbr > lvlend ) ; /* end do */
/* ----------------------------------------------------------
   we now haave the cuthill mckee ordering
   now reverse it
   --------------------------------------------------------*/
   revrse(ccsize, perm) ;

   return(ccsize) ;
}


/**************************************************************************
******************** ndegree . . . degree in masked component **************
***************************************************************************
        purpose - This routine computes the degrees of the nodes in the
                connected component specified by mask and root.
                nodes for which mask is -1 are ignored.

        input parameter -
                root - is the input node that defines the component
                padj - the adjacency structure
                mask - specifies a section subgraph
                work - must be initialized to zero before degree is
                       called !!!!

        output parameters
                deg - array containing the degrees of the nodes in
                        the component.
                ccsize - size of the component specified by mask and root

        working parameter -
                ls - a temporary vector used to store the nodes of the
                        component level by level.
************************************************************************/
int ndegree(int root, int **padj, int *mask, int *deg, int *ls, int *work)
{
   int ccsize ;
   int i, lbegin, lvlend, node, ideg, lvsize ;
   int *ptr ;
   int minone = -1 ;

/* ------------------------------------------
   initialization . . ..
   the array work is used to mark which nodes
   have been considered so far.
   -----------------------------------------*/
   ls[0] = root ;
   work[root] |= minone ;
   lvlend &= 0 ;
   ccsize = 1 ;
/* -------------------------------------------------------------
   lbegin is the pointer to the beginning of the current
   level, and lvlend points to the end of this level.
   ------------------------------------------------------------*/
   do
   {  lbegin = lvlend ;
      lvlend = ccsize ;
/*    ---------------------------------------------------------
      find the degrees of nodes in the current level,
      and at the same time, generate the next level.
      --------------------------------------------------------*/
      for (i=lbegin ; i < lvlend ; i++)
      {  node = ls[i] ;
         ideg = 0 ;
         for (ptr = padj[node] ; ptr < padj[node+1] ; ptr++)
         {  if (mask[*ptr] < 0) continue ;
            ideg++ ;
            if (work[*ptr] < 0) continue ;
            work[*ptr] |= minone ;
            ls[ccsize] = *ptr ;
            ccsize++ ;
         }
         deg[node] = ideg ;
      }
/*    -----------------------------------------------------
      compute the current level width.
      if it is nonnegative. generate another level.
      ---------------------------------------------------*/
      lvsize = ccsize - lvlend ;
   }  while (lvsize > 0) ;
/* ------------------------------------------------------
   reset work to 0
   -----------------------------------------------------*/
   for (i=0;i<ccsize ; i++)
   {  node = ls[i] ;
      work[node] &= 0 ;
   }

   return(ccsize) ;
}
