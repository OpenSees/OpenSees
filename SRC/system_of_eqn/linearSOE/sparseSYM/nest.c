/*
 * File:  nest.c
 * =============
 * Do the nested dissection ordering.
 */


#include <assert.h>
int i_greater();

#define INCLEVEL 6
/******************************************************************************
*************************** gennd . . . . general nested dissection  **********
*******************************************************************************

        purpose - subroutine gennd finds a nested dissection ordering for
                a general graph.

        input parameters -
                neqns - number of equations.
                padj - adjacency structure
        output parameters -
                perm - the nested dissection ordering

        working parameters
                mask - is used to mask off variables that haave been
                        numbered during the ordering process
                xls, ls - this level strucrre pari is used as temporary storage
                        by fnroot.
        routines called - 
                fndsep, revrse
*******************************************************************************/

void gennd(int neqns, int **padj, int *mask, int *perm, 
	   int *xls, int *ls, int *work)
{ 
   int num, i, root, nsep ;

   zeroi(neqns, mask) ;
   num = 0 ;
/* -------------------------------
   for each masked component
   -----------------------------*/
/* modified to operate on equations rather than nodes*/

   for (i=0;i<neqns ; i++)
   {
      while (mask[i] >= 0)
      {
	 root = i ;
/*       -----------------------------------------------------------
         find a separator and number the nodes next.
         ---------------------------------------------------------*/
         nsep = fndsep(root, padj, mask,(perm + num), xls, ls, work, neqns);
         num += nsep ;
      }
      if (num >= neqns) printf("breaking out at i %d nums %d neqns %d\n",i,num, neqns);
      if (num >= neqns ) break ;
   }

/* -----------------------------------------------------------------
   since separators found first should be ordered last, 
   routine revrse is called to adjust the ordering vector.
   ---------------------------------------------------------------*/
   revrse(neqns, perm) ;

   return ;
}


/**********************************************************
*****************fndsep . . . .  find separator   *********
***********************************************************

        purpose - this routine is used to find a small
                separator for a connected component
                specified by mask in the given graph.

        input parameters
                root - is the node that determines the masked
                        component
                padj - the adjacency structure

        output parameters
                nsep - number of variablesin the separator
                sep - vector containing the separator nodes

        updated parameter
                mask - nodes in the separator have their mask
                        values set to zero.

        working parameters
                xls, ls - level structur pair for level structure
                        found by fnroot

        routines called -
                fnroot
******************************************************************/
int fndsep(int root, int **padj, int *mask, int *sep, int *xls,
		   int *ls, int *work, int neqns)
{
   int nlvl, nsep, i, node, midlvl, midbeg, mp1beg, mp1end ;
   int *ptr ;
   int minone = -1 ;

   zeroi(neqns, work) ;
   root = fnroot(root, padj, mask, &nlvl, xls, ls) ;

/* ----------------------------------------------------------------
   if the number of levels is less than 3 (or some other value I
   choose for incomplete nested dissection), return the whole
   component as the separator.
   ---------------------------------------------------------------*/

   if (nlvl < INCLEVEL)
   {
      nsep = xls[nlvl+1] ;
      node = ls[0] ;
   /* renumber the segment by rcm */
      subrcm(nsep, node, padj, sep, mask, xls, work) ;

      /*return(nsep) */
     
      for (i=0;i<nsep;i++)
      {
	     node = ls[i] ;

         sep[i] = node ;
         mask[node] |= minone ;
      }
      return(nsep) ;
      
   }
  
/* ---------------------------------------------------------------
   find the middle level of the rooted level structure.
   ---------------------------------------------------------*/
   midlvl = (nlvl + 2)/2 ;
   { 
     int j, k, l ; ;
     k = xls[nlvl] ;
     j = k/ 2 ;
     l = 0 ;
     for (k = 0 ; k < nlvl && l < j ; k++)
     {
       if (l < j) {
	 l += (xls[k+1] - xls[k]) ;
       }
     }
 
     k-- ;
     midlvl = k ;
     midbeg = xls[midlvl] ;
     mp1beg = xls[midlvl + 1] ;
     mp1end = xls[midlvl + 2] ;
   }	 

/* ----------------------------------------------------------------
   the separator is obtained by including only those middle-level
   nodes with neighbors in the middle+1 level. work is used
   temporarily to mark those nodes in the middle+1 level.
   ---------------------------------------------------------------*/
   for (i = mp1beg;i<mp1end;i++)
   {
      node = ls[i] ;
      work[node] |= minone ;
   }
   nsep = 0 ;
   for (i=midbeg;i< mp1beg ; i++)
   {
      node = ls[i] ;
      for (ptr = padj[node] ; ptr < padj[node+1] ; ptr++)
      {
	 if (work[*ptr] < 0)
         {
	    sep[nsep] = node ;
            nsep++ ;
            mask[node] |= minone ;
            ptr = padj[node + 1] ; /* node has been added get out of ptr
                                      loop and start a new node */
         }
      }
   }
/* -------------------
   reset work
   ------------------*/
   for (i=mp1beg; i< mp1end; i++)
   {
      node = ls[i] ;
      work[i] &= 0 ;
   }

   return( nsep) ;
}


/***************************************************************************
***********************  fnroot ...  find pseudo-peripheral node ***********
****************************************************************************
        purpose - fnroot implements a modified version of the scheme
                by Gibbs, Polle, and Stockmeyer to find pseudo-peripheral
                nodes.  It determines such a node for the section 
                subgraph specified by mask and root.

        input parameters
                padj- adjacency structure
                mask - specifies a section subgraph. nodes for which
                        mask is less than zero are ignored by fnroot.
                root - on input, it (along with mask defines the
                        component for which a pseudo-peripheral node is 
                        to be found. 

        updated parameters
                nlvl -  the number of levels in the level structure rooted at
                        the node root.
        output
                newroot - the pseudo-peripheral node.
                (xls,ls)        the level structure array pair containing the
                        the level structure found

        subroutines called
                rootls
******************************************************************************/
int fnroot(int root, int **padj, int *mask, int *nlvl, int *xls, int *ls)
{
   int ccsize, j, jstrt, mindeg, ndeg, node ;
   int *kptr, nabor, nunlvl ;
   int oldroot ;
   
/* --------------------------------------------
   determine the level structure rooted at root
   --------------------------------------------*/
   oldroot = root ;
   *nlvl = rootls(root, padj, mask, xls, ls) ;
   ccsize = xls[*nlvl + 1] ;
   if (*nlvl == 0 || *nlvl == ccsize-1) return(root) ;
/* --------------------------------------------------
   pick a node with minimum degree from the last level
   ---------------------------------------------------*/
   do
   {
      jstrt = xls[*nlvl] ;
      mindeg = ccsize ;
      root  = ls[jstrt] ;
      if (ccsize != jstrt)
      {
	 for (j=jstrt; j<ccsize; j++)
         {
	    node = ls[j] ;
            ndeg = 0 ;
            for (kptr = padj[node] ; kptr < padj[node+1] ; kptr++)
            {
	       nabor = *kptr ;
               if (mask[nabor] >= 0) ndeg++ ;
            }
            if (ndeg < mindeg)
            {
	       root = node ;
               mindeg = ndeg ;
            }
         }
      }
/*    -----------------------------------------------
      and generate its rooted level structure
      -----------------------------------------------*/
      nunlvl = rootls(root, padj, mask, xls, ls) ;
      if ( nunlvl < *nlvl)
      {
	 root = oldroot ;
         *nlvl = rootls(root, padj, mask, xls, ls) ;
      }
      if ( nunlvl <= *nlvl) return(root) ;
      *nlvl = nunlvl ;
      oldroot = root ;
   } while (*nlvl < ccsize-1) ;
   
   return(root) ;
}


/**************************************************************************
**********************rootls ....... rooted level structure  **************
***************************************************************************
        Purpos - rootls generates the level structer rooted
                 at the inp8ut node called root.  Only those
                 nodes for which mask is nonnegative will considered.
        Input parameters
                root - the node at which the level structetr is to be
                       rooted
                padj - adjacency structure
                mask - is used to specify a section subgraph. 
                       nodes with mask[i] < 0 are ignored.

        Output parameters
                nlvl - is the number of levels in the level structure.
                xls,ls - array pair for the rooted level structure.
************************************************************************/

int rootls(int root, int **padj, int *mask, int *xls, int *ls)
{  
   int i, lbegin ;
   int lvlend, ccsize, node ;
   int nlvl ;
   int *jptr ;
   int minone = -1 ;

/* ---------------------------------
   initialization
   ---------------------------------*/
   mask[root] |= minone ;
   ls[0] = root ;
   nlvl |= minone ;
   lvlend &= 0 ;
   ccsize = 1 ;
/* ----------------------------------------------------------
   lbegin is the pointer to the begining of the current level
   and lvlend points to the end of this level
   ---------------------------------------------------------*/
   do
   {
      lbegin = lvlend ;
      lvlend = ccsize ;
      nlvl++ ;
      xls[nlvl] = lbegin ;
/*    -------------------------------------------------------
      generate the next level by finding all the masked 
      neighbors of nodes in the current level.
      ------------------------------------------------------*/
      for ( i = lbegin;i<lvlend ; i++)
      {
	 node = ls[i] ;
         for ( jptr = padj[node] ; jptr < padj[node+1] ; jptr++)
         {
	    if (mask[*jptr] < 0) continue ;
            ls[ccsize] = *jptr ;
            ccsize++ ;
            mask[*jptr] |= minone ;
         }
      }
/*    ----------------------------------------------------------
      compute the current level width. if it is nonzero generate
      the next level.
      --------------------------------------------------------*/
      /*`lvsize = ccsize - lvlend ; */
   }  while ( ccsize > lvlend ) ;
/* --------------------------------------------------------
   reset mask to zero for the nodes in the level structure 
   --------------------------------------------------------*/
   xls[nlvl+1] = lvlend ;
   for (i=0;i<ccsize;i++)
   {  node = ls[i] ;
      mask[node] &= 0 ;
   }

   return(nlvl) ;
}


/******************************************************************************
************************** revrse  . . . reverse a vector *********************
*******************************************************************************
        purpose - to reverse the order of a vector

*******************************************************************************/

revrse(int n, int *v)
{  
   int  *ve ;
   int temp ;

   ve = v + n - 1 ;

   while (ve > v)
   {  temp = *v ;
      *v = *ve ;
      *ve = temp ;
      ve-- ;
      v++ ;
   }
   
   return ;
}
   

/*****************************************************************************
**************************** forminv . . .  form invp ************************
******************************************************************************
        purpose - form the inverse of the perm vector
        input - neqns
                perm
        updated - invp
            
*****************************************************************************/
void forminv(int neqns, int *perm, int *invp)
{  
   int i ;
   for (i=0;i<neqns;i++)
      invp[perm[i]] = i ;
   return ;
}

