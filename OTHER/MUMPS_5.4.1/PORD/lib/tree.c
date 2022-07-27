/*****************************************************************************
/
/ SPACE SPArse Cholesky Elimination) Library: tree.c
/
/ author        J"urgen Schulze, University of Paderborn
/ created       09/15/99
/
/ This file contains functions dealing with elimination/front tree object
/
******************************************************************************

Data type:  struct elimtree
              int  nvtx;          number of vertices in the tree
              int  nfronts;       number of fronts in the tree
              int  root;          root of the tree
              int  *ncolfactor;   number of factor columns in front
              int  *ncolupdate;   number of update columns for front
              int  *parent;       parent in front tree
              int  *firstchild;   first child in front tree
              int  *silbings;     silbings in front tree
              int  *vtx2front;    maps vertices to fronts
Comments:
  o Structure used to hold the elimination/front tree; the tree is used to
    guide the symbolical and numerical factorization; a "node" in the tree
    can be a single vertex (in the context of an elimination tree) or a group
    of vertices (as for a front tree)
  o NOTE: Also the ordering can be expressed in terms of front trees; the
    permutation vector perm is then obtained by a post order traversal
    of the tree (see method permFromElimTree below)
Methods in lib/tree.c:
- T = newElimTree(int nvtx, int nfronts);
    o Initial: root = -1
- void freeElimTree(elimtree_t *T);
- void printElimTree(elimtree_t *T);
- int firstPostorder(elimtree_t *T);
    o returns the first front in a post order traversal of T
- int firstPostorder2(elimtree_t *T, int root);
    o returns the first front in a post order traversal of T[root]
- int nextPostorder(elimtree_t *T, int J);
    o returns the front that follows J in a post order traversal of T
- int firstPreorder(elimtree_t *T);
    o returns the first front in a pre order traversal of T
- int nextPreorder(elimtree_t *T, int J);
    o returns the front that follows J in a pre order traversal of T
- T = setupElimTree(graph_t *G, int *perm, int *invp);
    o constructs an elimination tree for G with permutation vectors perm,
      invp; a union-find algorithm is used to set up the parent vector of T;
      T->root and vectors T->firstchild, T->silbings are initialized by
      calling initFchSilbRoot; vector T->ncolupdate is filled by calling
      function setupCSSFromGraph (see below)
- void initFchSilbRoot(elimtree_t *T);
    o uses vector T->parent to initialize T->firstchild, T->silbings, T->root
- void permFromElimTree(elimtree_t *T, int *perm);
    o fills vectors perm, invp according to a post order traversal of T
- T2 = expandElimTree(elimtree_t *T, int *vtxmap, int nvtxorg)
    o creates and returns an elimtree object for the uncompressed graph;
      the map from original vertices to compressed vertices is found in vector 
      vtxmap; the number of original vertices (i.e. the length of vector
      vtxmap) is nvtxorg
    o NOTE: the function only expands vector T->vtx2front and sets
      T2->nvtx to nvtxorg; all other vectors are copied from T to T2, i.e.
      the number of fronts and the tree structure are the same in T and T2
- PTP = permuteElimTree(elimtree_t *T, int *perm);
    o in T: vtx2front[u] points to front containing vertex u
      in PTP: vtx2front[k] points to front containing column k = perm[u]
    o NOTE: the function only permutes vector T->vtx2front; all other vectors
      are copied from T to PTP, i.e. the number of fronts and the tree
      structure are the same in T and PTP
- T2 = fundamentalFronts(elimtree_t *T);
    o compresses chains of fronts to a single front; once a map from original
      fronts to compressed fronts is known, the compressed elimtree object T2
      can be created by calling compressElimTree (see below)
- T2 = mergeFronts(elimtree_t *T, int maxzeros);
    o merges small subtrees together in one front; it returns an elimtree
      object T2 where a front has either been merged with none or all of its
      children; the maximal number of zero entries that is allowed to be
      introduced when merging the fronts together is given by maxzeros
- T2 = compressElimTree(elimtree_t *T, int *frontmap, int cnfronts);
    o creates a new front tree using frontmap; vector frontmap maps the
      original fronts of T to a smaller set of fronts; cnfronts is number of
      new fronts (i.e. the maximal entry in frontmap)
- int justifyFronts(elimtree_t *T);
    o orders the children of a front so that the working storage in the
      multifrontal algorithm is minimized; the function returns the amount
      of extra working storage for the justified tree
- int nWorkspace(elimtree_t *T);
    o returns the size of the working storage in the multifrontal algorithm
      (measured in terms of FLOATS, for BYTES multiply with sizeof(FLOAT))
- int nFactorIndices(elimtree_t *T);
    o returns the number of indices taken by the factor matrix represented by T
- int nFactorEntries(elimtree_t *T);
    o returns the number of entries taken by the factor matrix represented by T
- FLOAT nFactorOps(elimtree_t *T);
    o returns the number of operations required to compute the factor matrix
      represented by T
- void subtreeFactorOps(elimtree *T, FLOAT *ops)
    o returns in ops[K] the number of operations required to factor the fronts
      in tree T(K) (this includes front K)
- FLOAT nTriangularOps(elimtree_t *T);
    o returns the number of operations required to solve the triangular systems

******************************************************************************/

#include <space.h>


/*****************************************************************************
******************************************************************************/
elimtree_t*
newElimTree(PORD_INT nvtx, PORD_INT nfronts)
{ elimtree_t *T;

  mymalloc(T, 1, elimtree_t);
  mymalloc(T->ncolfactor, nfronts, PORD_INT);
  mymalloc(T->ncolupdate, nfronts, PORD_INT);
  mymalloc(T->parent, nfronts, PORD_INT);
  mymalloc(T->firstchild, nfronts, PORD_INT);
  mymalloc(T->silbings, nfronts, PORD_INT);
  mymalloc(T->vtx2front, nvtx, PORD_INT);

  T->nvtx = nvtx;
  T->nfronts = nfronts;
  T->root = -1;

  return(T);
}


/*****************************************************************************
******************************************************************************/
void
freeElimTree(elimtree_t *T)
{
  free(T->ncolfactor);
  free(T->ncolupdate);
  free(T->parent);
  free(T->firstchild);
  free(T->silbings);
  free(T->vtx2front);
  free(T);
}


/*****************************************************************************
******************************************************************************/
void
printElimTree(elimtree_t *T)
{ PORD_INT *ncolfactor, *ncolupdate, *parent, *firstchild, *silbings, *vtx2front;
  PORD_INT *first, *link, nvtx, nfronts, root, J, K, u, count, child;

  nvtx = T->nvtx;
  nfronts = T->nfronts;
  root = T->root;
  ncolfactor = T->ncolfactor;
  ncolupdate = T->ncolupdate;
  parent = T->parent;
  firstchild = T->firstchild;
  silbings = T->silbings;
  vtx2front = T->vtx2front;

  printf("#fronts %d, root %d\n", nfronts, root);

  /* -----------------------------------------------------------
     store the vertices/columns of a front in a bucket structure
     ----------------------------------------------------------- */
  mymalloc(first, nfronts, PORD_INT);
  mymalloc(link, nvtx, PORD_INT);

  for (J = 0; J < nfronts; J++)
    first[J] = -1;
  for (u = nvtx-1; u >= 0; u--)
   { J = vtx2front[u];
     link[u] = first[J];
     first[J] = u;
   }

  /* -----------------------------------------------------------
     print fronts according to a postorder traversal of the tree
     ----------------------------------------------------------- */
  for (K = firstPostorder(T); K != -1; K = nextPostorder(T, K))
   { printf("--- front %d, ncolfactor %d, ncolupdate %d, parent %d\n",
            K, ncolfactor[K], ncolupdate[K], parent[K]);
     count = 0;
     printf("children:\n");
     for (child = firstchild[K]; child != -1; child = silbings[child])
      { printf("%5d", child);
        if ((++count % 16) == 0)
          printf("\n");
      }
     if ((count % 16) != 0)
       printf("\n");
     count = 0;
     printf("vertices mapped to front:\n");
     for (u = first[K]; u != -1; u = link[u])
      { printf("%5d", u);
        if ((++count % 16) == 0)
          printf("\n");
      }
     if ((count % 16) != 0)
       printf("\n");
   }

  /* ----------------------
     free memory and return
     ---------------------- */
  free(first); free(link);
}


/*****************************************************************************
******************************************************************************/
PORD_INT
firstPostorder(elimtree_t *T)
{ PORD_INT *firstchild, J;

  firstchild = T->firstchild;

  if ((J = T->root) != -1)
    while (firstchild[J] != -1)
      J = firstchild[J];
  return(J);
}


/*****************************************************************************
******************************************************************************/
PORD_INT
firstPostorder2(elimtree_t *T, PORD_INT root)
{ PORD_INT *firstchild, J;

  firstchild = T->firstchild;

  if ((J = root) != -1)
    while (firstchild[J] != -1)
      J = firstchild[J];
  return(J);
}


/*****************************************************************************
******************************************************************************/
PORD_INT
nextPostorder(elimtree_t *T, PORD_INT J)
{ PORD_INT *parent, *firstchild, *silbings;

  parent = T->parent;
  firstchild = T->firstchild;
  silbings = T->silbings;

  if (silbings[J] != -1)
   { J = silbings[J];
     while (firstchild[J] != -1)
       J = firstchild[J];
   }
  else J = parent[J];
  return(J);
}


/*****************************************************************************
******************************************************************************/
PORD_INT
firstPreorder(elimtree_t *T)
{
  return(T->root);
}


/*****************************************************************************
******************************************************************************/
PORD_INT
nextPreorder(elimtree_t *T, PORD_INT J)
{ PORD_INT *parent, *firstchild, *silbings;

  parent = T->parent;
  firstchild = T->firstchild;
  silbings = T->silbings;

  if (firstchild[J] != -1)
    J = firstchild[J];
  else
   { while ((silbings[J] == -1) && (parent[J] != -1))
       J = parent[J];
     J = silbings[J];
   }
  return(J);
}


/*****************************************************************************
******************************************************************************/
elimtree_t*
setupElimTree(graph_t *G, PORD_INT *perm, PORD_INT *invp)
{ elimtree_t *T;
  css_t      *css;
  PORD_INT        *xadj, *adjncy, *vwght, *ncolfactor, *ncolupdate, *parent;
  PORD_INT        *vtx2front, *realroot, *uf_father, *uf_size;
  PORD_INT        *xnzl, *nzlsub, *xnzlsub;
  PORD_INT        nvtx, front, front2, froot, f, r, u, v, i, istart, istop;
  PORD_INT        prevlen, len, h, hsub;

  nvtx = G->nvtx;
  xadj = G->xadj;
  adjncy = G->adjncy;
  vwght = G->vwght;

  /* --------------------------
     set up the working storage
     -------------------------- */
  mymalloc(realroot, nvtx, PORD_INT);
  mymalloc(uf_father, nvtx, PORD_INT);
  mymalloc(uf_size, nvtx, PORD_INT);

  /* ------------------------------------------------
     allocate storage for the elimination tree object
     ------------------------------------------------ */
  T = newElimTree(nvtx, nvtx);
  ncolfactor = T->ncolfactor;
  ncolupdate = T->ncolupdate;
  parent = T->parent;
  vtx2front = T->vtx2front;

  /* -----------------------------
     set up the parent vector of T
     ----------------------------- */
  for (front = 0; front < nvtx; front++)
   { parent[front] = -1;
     u = invp[front];           /* only vertex u belongs to this front */
     uf_father[front] = front;  /* front forms a set in union-find structure */
     uf_size[front] = 1;        /* the set consists of a single front */
     realroot[front] = front;
     froot = front;

     /* run through the adjacency list of u */
     istart = xadj[u];
     istop = xadj[u+1];
     for (i = istart; i < istop; i++)
      { v = adjncy[i];
        front2 = perm[v];
        if (front2 < front)
         { r = front2;

           while (uf_father[r] != r)  /* find root of front2 in union-find */
             r = uf_father[r];
           while (front2 != r)        /* path compression */
            { f = front2;
              front2 = uf_father[front2];
              uf_father[f] = r;
            }

           f = realroot[r];           /* merge union-find sets */
           if ((parent[f] == -1) && (f != front))
            { parent[f] = front;
              if (uf_size[froot] < uf_size[r])
               { uf_father[froot] = r;
                 uf_size[r] += uf_size[froot];
                 froot = r;
               }
              else
               { uf_father[r] = froot;
                 uf_size[froot] += uf_size[r];
               }
              realroot[froot] = front;
            }
         }
      }
   }

  /* ---------------------------------------------
     set the vectors T->firstchild and T->silbings
     --------------------------------------------- */
  initFchSilbRoot(T);

  /* ----------------------------------------------------------
     set the vectors T->vtx2front, T->ncolfactor, T->ncolupdate
     ---------------------------------------------------------- */
  css = setupCSSFromGraph(G, perm, invp);
  xnzl = css->xnzl;
  nzlsub = css->nzlsub;
  xnzlsub = css->xnzlsub;

  prevlen = 0;
  for (front = 0; front < nvtx; front++)
   { u = invp[front];
     ncolfactor[front] = vwght[u];
     ncolupdate[front] = 0;
     vtx2front[u] = front;
     len = xnzl[front+1] - xnzl[front];
     if (prevlen - 1 == len)
       ncolupdate[front] = ncolupdate[front-1] - vwght[u];
     else
      { h = xnzlsub[front] + 1;
        for (i = 1; i < len; i++)
         { hsub = nzlsub[h++];
           v = invp[hsub];
           ncolupdate[front] += vwght[v];
         }
      }
     prevlen = len;
   }

  /* ----------------------
     free memory and return
     ---------------------- */
  free(css);
  free(realroot); free(uf_father); free(uf_size);
  return(T);
}


/*****************************************************************************
******************************************************************************/
void
initFchSilbRoot(elimtree_t *T)
{ PORD_INT *parent, *firstchild, *silbings, nfronts, J, pJ;

  nfronts = T->nfronts;
  parent = T->parent;
  firstchild = T->firstchild;
  silbings = T->silbings;

  for (J = 0; J < nfronts; J++)
    silbings[J] = firstchild[J] = -1;

  for (J = nfronts-1; J >= 0; J--)
    if ((pJ = parent[J]) != -1)
     { silbings[J] = firstchild[pJ];
       firstchild[pJ] = J;
     }
    else
     { silbings[J] = T->root;
       T->root = J;
     }
}


/*****************************************************************************
******************************************************************************/
void
permFromElimTree(elimtree_t *T, PORD_INT *perm)
{ PORD_INT *vtx2front, *first, *link;
  PORD_INT nvtx, nfronts, K, u, count;

  nvtx = T->nvtx;
  nfronts = T->nfronts;
  vtx2front = T->vtx2front;

  /* -----------------------------------------------------------
     store the vertices/columns of a front in a bucket structure
     ----------------------------------------------------------- */
  mymalloc(first, nfronts, PORD_INT);
  mymalloc(link, nvtx, PORD_INT);

  for (K = 0; K < nfronts; K++)
    first[K] = -1;
  for (u = nvtx-1; u >= 0; u--)
   { K = vtx2front[u];
     link[u] = first[K];
     first[K] = u;
   }

  /* -----------------------------------------------------
     postorder traversal of the elimination tree to obtain
     the permutation vectors perm, invp
     ----------------------------------------------------- */
  count = 0;
  for (K = firstPostorder(T); K != -1; K = nextPostorder(T, K))
    for (u = first[K]; u != -1; u = link[u])
     { perm[u] = count;
       count++;
     }

  /* ----------------------
     free memory and return
     ---------------------- */
  free(first); free(link);
}


/*****************************************************************************
******************************************************************************/
elimtree_t*
permuteElimTree(elimtree_t *T, PORD_INT *perm)
{ elimtree_t *PTP;
  PORD_INT        nvtx, nfronts, J, u;

  nvtx = T->nvtx;
  nfronts = T->nfronts;

  /* --------------------------------------------------------------
     allocate space for the new elimtree object and copy front data
     the permuted tree has the same number of fronts/tree structure
     -------------------------------------------------------------- */
  PTP = newElimTree(nvtx, nfronts);
  PTP->root = T->root;
  for (J = 0; J < nfronts; J++)
   { PTP->ncolfactor[J] = T->ncolfactor[J];
     PTP->ncolupdate[J] = T->ncolupdate[J];
     PTP->parent[J] = T->parent[J];
     PTP->firstchild[J] = T->firstchild[J];
     PTP->silbings[J] = T->silbings[J];
   }

  /* ---------------------------------------------------------------------
     set up the new vtx2front vector; the trees only differ in this vector
     --------------------------------------------------------------------- */
  for (u = 0; u < nvtx; u++)
    PTP->vtx2front[perm[u]] = T->vtx2front[u];

  return(PTP);
}


/*****************************************************************************
******************************************************************************/
elimtree_t*
expandElimTree(elimtree_t *T, PORD_INT *vtxmap, PORD_INT nvtxorg)
{ elimtree_t *T2;
  PORD_INT        *vtx2front, *vtx2front2;
  PORD_INT        nfronts, J, u;

  nfronts = T->nfronts;

  /* --------------------------------------------------------------
     allocate space for the new elimtree object and copy front data
     the expanded tree has the same number of fronts/tree structure
     -------------------------------------------------------------- */
  T2 = newElimTree(nvtxorg, nfronts);
  T2->root = T->root;
  for (J = 0; J < nfronts; J++)
   { T2->ncolfactor[J] = T->ncolfactor[J];
     T2->ncolupdate[J] = T->ncolupdate[J];
     T2->parent[J] = T->parent[J];
     T2->firstchild[J] = T->firstchild[J];
     T2->silbings[J] = T->silbings[J];
   }

  /* ---------------------------------------------------------------------
     set up the new vtx2front vector; the trees only differ in this vector
     --------------------------------------------------------------------- */
  vtx2front = T->vtx2front;
  vtx2front2 = T2->vtx2front;
  for (u = 0; u < nvtxorg; u++)
    vtx2front2[u] = vtx2front[vtxmap[u]];

  return(T2);
}


/*****************************************************************************
******************************************************************************/
elimtree_t*
fundamentalFronts(elimtree_t *T)
{ elimtree_t *T2;
  PORD_INT        *ncolfactor, *ncolupdate, *parent, *firstchild, *silbings;
  PORD_INT        *frontmap, nfronts, cnfronts, J, child;

  nfronts = T->nfronts;
  ncolfactor = T->ncolfactor;
  ncolupdate = T->ncolupdate;
  parent = T->parent;
  firstchild = T->firstchild;
  silbings = T->silbings;

  /* -------------------------
     set up the working arrays
     ------------------------- */
  mymalloc(frontmap, nfronts, PORD_INT);

  /* -----------------------------
     search the fundamental fronts
     ----------------------------- */
  cnfronts = 0;
  J = T->root;
  while (J != -1)
   { while (firstchild[J] != -1)
       J = firstchild[J];
     frontmap[J] = cnfronts++;
     while ((silbings[J] == -1) && (parent[J] != -1))
      { J = parent[J];
        child = firstchild[J];
        if ((silbings[child] != -1)
            || (ncolupdate[child] != ncolfactor[J] + ncolupdate[J]))
          frontmap[J] = cnfronts++;
        else
          frontmap[J] = frontmap[child];
      }
     J = silbings[J];
   }

  /* ------------------------------
     construct new elimination tree
     ------------------------------ */
  T2 = compressElimTree(T, frontmap, cnfronts);

  /* ----------------------
     free memory and return
     ---------------------- */
  free(frontmap);
  return(T2);
}


/*****************************************************************************
******************************************************************************/
elimtree_t*
mergeFronts(elimtree_t *T, PORD_INT maxzeros)
{ elimtree_t *T2;
  PORD_INT        *ncolfactor, *ncolupdate, *firstchild, *silbings;
  PORD_INT        *frontmap, *newncolfactor, *nzeros, *rep;
  PORD_INT        nfronts, cnfronts, K, ncolfrontK, J, Jall, cost;

  nfronts = T->nfronts;
  ncolfactor = T->ncolfactor;
  ncolupdate = T->ncolupdate;
  firstchild = T->firstchild;
  silbings = T->silbings;

  /* -------------------------
     set up the working arrays
     ------------------------- */
  mymalloc(frontmap, nfronts, PORD_INT);
  mymalloc(newncolfactor, nfronts, PORD_INT);
  mymalloc(nzeros, nfronts, PORD_INT);
  mymalloc(rep, nfronts, PORD_INT);
  for (K = 0; K < nfronts; K++)
   { newncolfactor[K] = ncolfactor[K];
     nzeros[K] = 0;
     rep[K] = K;
   }

  /* -----------------------------------------------------
     perform a postorder traversal of the elimination tree
     ----------------------------------------------------- */
  for (K = firstPostorder(T); K != -1; K = nextPostorder(T, K))
    if (firstchild[K] != -1)
     { ncolfrontK = newncolfactor[K] + ncolupdate[K];
       Jall = 0;
       cost = 0;
       for (J = firstchild[K]; J != -1; J = silbings[J])
        { Jall += newncolfactor[J];
          cost -= newncolfactor[J] * newncolfactor[J];
          cost += 2*newncolfactor[J] * (ncolfrontK - ncolupdate[J]);
          cost += 2*nzeros[J];
        }
       cost += Jall * Jall;
       cost = cost / 2;
       if (cost < maxzeros)
        { for (J = firstchild[K]; J != -1; J = silbings[J])
           { rep[J] = K;
             newncolfactor[K] += newncolfactor[J];
           }
          nzeros[K] = cost;
        }
     }

  /* ----------------------------------
     construct frontmap from vector rep
     ---------------------------------- */
  cnfronts = 0;
  for (K = 0; K < nfronts; K++)
    if (rep[K] == K)
      frontmap[K] = cnfronts++;
    else
     { for (J = K; rep[J] != J; J = rep[J]);
       rep[K] = J;
     }
  for (K = 0; K < nfronts; K++)
    if ((J = rep[K]) != K)
      frontmap[K] = frontmap[J];

  /* ------------------------------
     construct new elimination tree
     ------------------------------ */
  T2 = compressElimTree(T, frontmap, cnfronts);

  /* ----------------------
     free memory and return
     ---------------------- */
  free(frontmap); free(newncolfactor);
  free(nzeros); free(rep);
  return(T2);
}


/*****************************************************************************
******************************************************************************/
elimtree_t*
compressElimTree(elimtree_t *T, PORD_INT *frontmap, PORD_INT cnfronts)
{ elimtree_t *T2;
  PORD_INT        *ncolfactor, *ncolupdate, *parent, *vtx2front;
  PORD_INT        nvtx, nfronts, u, K, pK, newfront, pnewfront;

  nvtx = T->nvtx;
  nfronts = T->nfronts;
  ncolfactor = T->ncolfactor;
  ncolupdate = T->ncolupdate;
  parent = T->parent;
  vtx2front = T->vtx2front;

  /* --------------------------------------------
     allocate memory for the new elimtree T2
     and init. ncolfactor, ncolupdate, and parent
     -------------------------------------------- */
  T2 = newElimTree(nvtx, cnfronts);
  for (K = 0; K < cnfronts; K++)
   { T2->ncolfactor[K] = T2->ncolupdate[K] = 0;
     T2->parent[K] = -1;
   }

  /* --------------------------------------------------------------
     set the new vectors T2->ncolfactor, T2->ncolupdate, T2->parent
     -------------------------------------------------------------- */
  for (K = 0; K < nfronts; K++)
   { newfront = frontmap[K];
     T2->ncolfactor[newfront] += ncolfactor[K];
     if (((pK = parent[K]) != -1)
        && ((pnewfront = frontmap[pK]) != newfront))
      { T2->parent[newfront] = pnewfront;
        T2->ncolupdate[newfront] = ncolupdate[K];
      }
   }

  /* ---------------------------------------------------
     set the new vectors T2->firstchild and T2->silbings
     --------------------------------------------------- */
  initFchSilbRoot(T2);

  /* ------------------------------------
     set the the new vector T2->vtx2front
     ------------------------------------ */
  for (u = 0; u < nvtx; u++)
    T2->vtx2front[u] = frontmap[vtx2front[u]];
  return(T2);
}


/*****************************************************************************
******************************************************************************/
PORD_INT
justifyFronts(elimtree_t *T)
{ PORD_INT *ncolfactor, *ncolupdate, *firstchild, *silbings, *minWspace, *list;
  PORD_INT nfronts, K, ncolfrontK, frontsizeK, wspace, child, nxtchild;
  PORD_INT count, m, s, i;

  nfronts = T->nfronts;
  ncolfactor = T->ncolfactor;
  ncolupdate = T->ncolupdate;
  firstchild = T->firstchild;
  silbings = T->silbings;

  /* -------------------------
     set up the working arrays
     ------------------------- */
  mymalloc(minWspace, nfronts, PORD_INT);
  mymalloc(list, nfronts, PORD_INT);

  /* ---------------------------------------------------------
     postorder traversal of the elimination tree to obtain the
     optimal justification of the children of each front
     ---------------------------------------------------------- */
  wspace = 0;
  for (K = firstPostorder(T); K != -1; K = nextPostorder(T, K))
   { ncolfrontK = ncolfactor[K] + ncolupdate[K];
     frontsizeK = (ncolfrontK * (ncolfrontK + 1)) >> 1;

     if ((child = firstchild[K]) == -1)
       minWspace[K] = frontsizeK;
     else
      { count = 0;

        /* sort children according to their minWspace value */
        while (child != -1)
         { list[count++] = child;
           child = silbings[child];
         }
        insertUpIntsWithStaticIntKeys(count, list, minWspace);
        firstchild[K] = -1;
        for (i = 0; i < count; i++)
         { child = list[i];
           silbings[child] = firstchild[K];
           firstchild[K] = child;
         }

        /* compute minWspace[K] */
        child = firstchild[K];
        nxtchild = silbings[child];
        m = s = minWspace[child];
        while (nxtchild != -1)
         { s = s - minWspace[child]
               + ((ncolupdate[child] * (ncolupdate[child] + 1)) >> 1)
               + minWspace[nxtchild];
           m = max(m, s);
           child = nxtchild;
           nxtchild = silbings[nxtchild];
         }
        s = s - minWspace[child]
            + ((ncolupdate[child] * (ncolupdate[child] + 1)) >> 1)
            + frontsizeK;
        minWspace[K] = max(m, s);
      }

     wspace = max(wspace, minWspace[K]);
   }

  /* ----------------------
     free memory and return
     ---------------------- */
  free(minWspace); free(list);
  return(wspace);
}


/*****************************************************************************
******************************************************************************/
PORD_INT
nWorkspace(elimtree_t *T)
{ PORD_INT *ncolfactor, *ncolupdate, *firstchild, *silbings, *minWspace;
  PORD_INT nfronts, K, ncolfrontK, frontsizeK, wspace, child, nxtchild, m, s;

  nfronts = T->nfronts;
  ncolfactor = T->ncolfactor;
  ncolupdate = T->ncolupdate;
  firstchild = T->firstchild;
  silbings = T->silbings;

  /* -------------------------
     set up the working arrays
     ------------------------- */
  mymalloc(minWspace, nfronts, PORD_INT);

  /* -------------------------------------------
     postorder traversal of the elimination tree
     ------------------------------------------- */
  wspace = 0;
  for (K = firstPostorder(T); K != -1; K = nextPostorder(T, K))
   { ncolfrontK = ncolfactor[K] + ncolupdate[K];
     frontsizeK = (ncolfrontK * (ncolfrontK + 1)) >> 1;

     if ((child = firstchild[K]) == -1)
       minWspace[K] = frontsizeK;
     else
      { child = firstchild[K];
        nxtchild = silbings[child];
        m = s = minWspace[child];
        while (nxtchild != -1)
         { s = s - minWspace[child]
               + ((ncolupdate[child] * (ncolupdate[child] + 1)) >> 1)
               + minWspace[nxtchild];
           m = max(m, s);
           child = nxtchild;
           nxtchild = silbings[nxtchild];
         }
        s = s - minWspace[child]
            + ((ncolupdate[child] * (ncolupdate[child] + 1)) >> 1)
            + frontsizeK;
        minWspace[K] = max(m, s);
      }

     wspace = max(wspace, minWspace[K]);
   }

  /* ----------------------
     free memory and return
     ---------------------- */
  free(minWspace);
  return(wspace);
}


/*****************************************************************************
******************************************************************************/
PORD_INT
nFactorIndices(elimtree_t *T)
{ PORD_INT *ncolfactor, *ncolupdate;
  PORD_INT nfronts, ind, K;

  nfronts = T->nfronts;
  ncolfactor = T->ncolfactor;
  ncolupdate = T->ncolupdate;

  ind = 0;
  for (K = 0; K < nfronts; K++)
    ind += (ncolfactor[K] + ncolupdate[K]);
  return(ind);
}


/*****************************************************************************
******************************************************************************/
PORD_INT
nFactorEntries(elimtree_t *T)
{ PORD_INT *ncolfactor, *ncolupdate;
  PORD_INT ent, tri, rec, K;

  ncolfactor = T->ncolfactor;
  ncolupdate = T->ncolupdate;

  ent = 0;
  for (K = firstPostorder(T); K != -1; K = nextPostorder(T, K))
   { tri = ncolfactor[K];
     rec = ncolupdate[K];
     ent += (tri * (tri+1)) / 2;
     ent += (tri * rec);
   }
  return(ent);
}


/*****************************************************************************
******************************************************************************/
FLOAT
nFactorOps(elimtree_t *T)
{ PORD_INT   *ncolfactor, *ncolupdate;
  FLOAT ops, tri, rec;
  PORD_INT   K;

  ncolfactor = T->ncolfactor;
  ncolupdate = T->ncolupdate;

  ops = 0.0;
  for (K = firstPostorder(T); K != -1; K = nextPostorder(T, K))
   { tri = ncolfactor[K];
     rec = ncolupdate[K];
     ops += (tri*tri*tri) / 3.0 + (tri*tri) / 2.0 - (5*tri) / 6.0;
     ops += (tri*tri*rec) + (rec*(rec+1)*tri);
   }
  return(ops);
}


/*****************************************************************************
******************************************************************************/
void
subtreeFactorOps(elimtree_t *T, FLOAT *ops)
{ PORD_INT   *ncolfactor, *ncolupdate;
  FLOAT tri, rec;
  PORD_INT   J, K;

  ncolfactor = T->ncolfactor;
  ncolupdate = T->ncolupdate;

  for (K = firstPostorder(T); K != -1; K = nextPostorder(T, K))
   { tri = ncolfactor[K];
     rec = ncolupdate[K];
     ops[K] = (tri*tri*tri) / 3.0 + (tri*tri) / 2.0 - (5*tri) / 6.0;
     ops[K] += (tri*tri*rec) + (rec*(rec+1)*tri);
     for (J = T->firstchild[K]; J != -1; J = T->silbings[J])
       ops[K] += ops[J];
   }
}


/*****************************************************************************
******************************************************************************/
FLOAT
nTriangularOps(elimtree_t *T)
{ PORD_INT   *ncolfactor, *ncolupdate;
  FLOAT ops, tri, rec;
  PORD_INT   K;

  ncolfactor = T->ncolfactor;
  ncolupdate = T->ncolupdate;

  ops = 0.0;
  for (K = firstPostorder(T); K != -1; K = nextPostorder(T, K))
   { tri = ncolfactor[K];
     rec = ncolupdate[K];
     ops += (tri*tri) + 2.0*tri*rec;  /* forward ops */
     ops += (tri*tri) + 2.0*tri*rec;  /* backward ops */
   }
  return(ops);
}

