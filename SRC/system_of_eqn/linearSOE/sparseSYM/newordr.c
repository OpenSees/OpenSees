/*
 * File:  newordr.c
 * ================
 * Subroutines in this file form the elimination tree, and postorder
 * the elimination tree.
 * March 1990 David Mackay
 */


/************************************************************************
 ************  pfordr ..... profile reordering  *************************
 ************************************************************************
 
    purpose - to form elimination  tree,
              and perform an equivalent post ordering 
              suitable for the sparse generalized profile method
 
    input parameters -
        neqns - no of equations
        padj - adjacency structure
        wperm - original perm vector
        winvp - original invp vector
    output parameters -
        (nperm,ninvp) - new permutation and inverse permutation vectors
        parent - parent vector of elimination tree
        nblks - number of blocks 
        list - a linked list of begining row/columns 
    working parameters -
        fchild - first son vector
        sibling - sibling vector
    subroutines - etree, bntree, postordr
 
 ************************************************************************/
     
int pfordr(int neqns, int **padj, int *perm, int *invp, int *parent, int *fchild, 
	   int *sibling, int *winvp, int *wperm, int *list, int *rowblks)
{  
   int nblks;
   int *i, j;

   if ( neqns <= 0 )  return(0);
/*       form elimination tree */
   etree ( neqns, padj, wperm, winvp, parent, fchild );
/*         -------------------------------------------------------
           obtain a binary representation of the elim tree, and
           then perform postordering 
           -------------------------------------------------------
*/
   bntree ( neqns, parent, fchild, sibling ) ;
   zeroi(neqns, list ) ;
   list[0] = neqns ;
   minoni(neqns,list) ;
   postordr ( neqns-1, parent, fchild, sibling, winvp, wperm,
      invp, perm, list, rowblks  ) ;


/* count number of blocks */
/*   temp1 = list ; */

/* this sets parent to be the ancestor array */
   nblks = 0 ;
   i = parent ;
   while (*list >=0 ) 
   {  j = parent[list[1] - 1 ] ;
      for ( ; i < parent + list[1]; i++)
	 *i = j ;
      nblks++ ;
      list++ ;
   }
   *list = neqns ;
   for ( ; i < parent + neqns ; i++)
      *i = neqns ;

   return(nblks) ;

}

/*
  revised and written in c by David Mackay Jan 1990
  
  acknowledgements:
    this routine is based on a fortran routine
    written and owned by dr. joseph liu,
    department of computer science, york university.
 
 ***********************************************************************
 ****************     etree ..... elimination tree     *****************
 ***********************************************************************
 
        purpose -
            this subroutine computes the elimination tree from a given
            ordering and adjacency structure.
 
        input parameters -
            neqns       - number of equations.
            padj        - the adjacency structure.
            (perm,invp) - the permutation and inverse permutation
                          vectors.
 
        output parameters -
            parent      - the parent vector of the elimination tree.
 
        Working parameters -
            ancstr      - the ancestor vector.
 
 ***********************************************************************/

etree(int neqns, int **padj, int *perm, int *invp, int *parent, int *ancstr)
{  
   int  i, nbr, next, node, mone;
   int *pt ;

   mone = -1 ;

   for (i = 0;i<neqns;i++)
   {  parent[i] |= mone ;
      ancstr[i] |= mone ;
      node = perm[i] ;
      for (pt = padj[node ] ; pt < padj[node+1] ; pt++)
      {  nbr = invp[*pt] ;
         if  ( nbr >=  i )  continue ;
         while(ancstr[nbr] >= 0 && ancstr[nbr] != i)
         {  next = ancstr[nbr] ;
            ancstr[nbr] = i  ;
            nbr = next ;
         }
         if ( ancstr[nbr] <0)
         {  parent[nbr] = i;
            ancstr[nbr] = i ;
         }
      }
   }
   parent[neqns-1] = neqns ;

   return;
}


/************************************************************************
 ************  bntree ....binary tree representation  *******************
 ************************************************************************
 
    purpose - to set up the binary tree representation 
 
    input parameters -
        neqns - no of equations
        parent - parent vector of elimination tree

    output parameters -
        fchild - first son vector (left children)
        sibling - siblings or right children

  modeled from a fortran program bntree by Kincho Law
 
************************************************************************/
     
bntree (int neqns, int *parent, int *fchild, int *sibling)
{
   int node, p ;

   minoni(neqns,fchild) ;
   minoni(neqns,sibling) ;

/* start processing */

   for (node = 0; node < neqns; node++)
   {  p = parent[node] ;
      if (p >= neqns ) continue ;
      if (fchild[p] == -1)
         fchild[p] = node ;
      else
      {  sibling[node] = fchild[p] ;
         fchild[p] = node ;
      }
   }

   return ;
}

/************************************************************************
 ************  postordr ....post order elimination tree *****************
 ************************************************************************
 
    purpose - to form an equivalent post ordering of the elimination
              tree.  It also forms a list of the begining of 
	      new blocks in the generalized profile matrix.
 
    input parameters -
        i - begin postordering at i from first call this should
            be the last equation number
        fchild - first child vector
        sibling - sibling vector
        oinvp -  old invp vector 
        operm - old perm vector 

    updated parameters -
        ninvp -  new invp vector
        nperm -  new perm vector
        parent - new parent vector 

modeled from a recursive postorder traversal in
Fundamentals of Data Structures by Horowitz and Sahni 
page 231
	additional info like forming the link list to form xblk 
        and the number of blocks has been tacked on

        list 
 
 ************************************************************************/
 static int count = 0 ;
 static int xcount = -1 ;
     
postordr(int i, int *parent, int *fchild, int *sibling, int *oinvp, int *operm, 
	 int *ninvp, int *nperm, int *list, int *rowblks)
{
   int t ;

/* ----------------------------------------------------------
      postordr all children
   ----------------------------------------------------------*/
   if (fchild[i] >= 0)  postordr(fchild[i],parent,fchild,
      sibling, oinvp, operm, ninvp, nperm, list, rowblks  ) ;
   else
   {  /* add count to list of nodes begining a new block */
      xcount++ ;
      list[xcount] = count;
   }
/* --------------------------------------------------------
      renumber node i
   --------------------------------------------------------*/
   t = operm[i] ;
   nperm[count] = t ;
   ninvp[t] = count ;
/* save new number in t to update the parent node vector after sibligns are
   renumbered and the value of the new parent number is determined */ 
   t = count ;
   rowblks[count] = xcount ;
   count++ ;
/* update subtrees of siblings */
   if (sibling[i] >= 0)
   {  postordr(sibling[i],parent,fchild, sibling, oinvp, operm, 
	 ninvp, nperm, list, rowblks ) ;
      /* add parent to list of nodes begining a new block */
      /* comment out below to cut down number of blocks and 
	 follow Liu's scheme */ 
      if (list[xcount] != count)
      {  xcount++ ;
         list[xcount] = count ;
      }
      /**/
/*    update the parent vector */ 
      parent[t] = count ;
      parent[count - 1] = count ;
   }
   
   return ;
}

/************************************************************************
 ************  pfblk ..... profile      blocks **************************
 ************************************************************************
 
    purpose -  to set up xblk - index of begining row/column of each block
 
    input parameters -
        neqns - no of equations
        list - a linked list of begining row/columns of each block
    output parameters -
        xblk - index vector for blocks
 
 ************************************************************************/
     
pfblk (int nblks, int *xblk, int *list)
{ 
   int *stop ;

   stop = list + nblks ;
   for ( ; list <= stop; list++, xblk++)
      *xblk = *list ;

   return ;
}
