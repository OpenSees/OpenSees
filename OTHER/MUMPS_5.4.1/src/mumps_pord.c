/*
 *
 *  This file is part of MUMPS 5.4.1, released
 *  on Tue Aug  3 09:49:43 UTC 2021
 *
 *
 *  Copyright 1991-2021 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  Mumps Technologies, University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license 
 *  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
 *  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
 *
 */
/*
 * This file contains interfaces to external ordering packages.
 * At the moment, PORD (J. Schulze) and SCOTCH are interfaced.
 */
#include "mumps_pord.h"
void MUMPS_CALL MUMPS_PORD_INTSIZE(MUMPS_INT *pord_intsize)
{
#if defined(pord)
#   if defined(PORD_INTSIZE64) || defined(INTSIZE64)
    *pord_intsize=64;
#   else
    *pord_intsize=32;
#   endif
#else
    *pord_intsize=-99999;
#endif
}
#if defined(pord)
/* Interface to PORD */
#if defined(INTSIZE64) || defined(PORD_INTSIZE64)
void MUMPS_CALL
MUMPS_PORDF( MUMPS_INT8 *nvtx, MUMPS_INT8 *nedges,
             MUMPS_INT8 *xadj, MUMPS_INT8 *adjncy,
             MUMPS_INT8 *nv, MUMPS_INT *ncmpa )
#else
void MUMPS_CALL
MUMPS_PORDF( MUMPS_INT *nvtx, MUMPS_INT *nedges,
             MUMPS_INT *xadj, MUMPS_INT *adjncy,
             MUMPS_INT *nv, MUMPS_INT *ncmpa )
#endif
{
    *ncmpa = mumps_pord( *nvtx, *nedges, xadj, adjncy, nv );
}
/* Interface to PORD with weighted graph */
#if defined(INTSIZE64) || defined(PORD_INTSIZE64)
void MUMPS_CALL
MUMPS_PORDF_WND( MUMPS_INT8 *nvtx, MUMPS_INT8 *nedges,
                 MUMPS_INT8 *xadj, MUMPS_INT8 *adjncy,
                 MUMPS_INT8 *nv, MUMPS_INT *ncmpa, MUMPS_INT8 *totw )
#else
void MUMPS_CALL
MUMPS_PORDF_WND( MUMPS_INT *nvtx, MUMPS_INT *nedges,
                 MUMPS_INT *xadj, MUMPS_INT *adjncy,
                 MUMPS_INT *nv, MUMPS_INT *ncmpa, MUMPS_INT *totw )
#endif
{
    *ncmpa = mumps_pord_wnd( *nvtx, *nedges, xadj, adjncy, nv, totw );
}
/************************************************************
 mumps_pord is used in ana_aux.F
        permutation and inverse permutation not set on output,
        but may be printed in default file: "perm_pord" and "iperm_pord",
        if associated part uncommneted.
        But, if uncommetnted a bug occurs in psl_ma41_analysi.F
******************************************************************/
/*********************************************************/
MUMPS_INT mumps_pord
(
   PORD_INT nvtx,
   PORD_INT nedges,      /* NZ-like */
   PORD_INT *xadj_pe,    /* NZ-like */
   PORD_INT *adjncy,
   PORD_INT *nv
)
{
/**********************************
Arguments:
input:
-----
- nvtx          : dimension of the Problem (N)
- nedges        : number of entries (NZ)
- adjncy        : non-zeros entries (IW input)
input/output:
-------------
- xadj_pe       : in: pointer through beginning of column non-zeros entries
                  out: "father array" (PE)
ouput:
------
- nv            : "nfront array" (NV)
*************************************/
  graph_t    *G;
  elimtree_t *T;
  timings_t  cpus[12];
  options_t  options[] = { SPACE_ORDTYPE, SPACE_NODE_SELECTION1,
                    SPACE_NODE_SELECTION2, SPACE_NODE_SELECTION3,
                    SPACE_DOMAIN_SIZE, 0 };
  PORD_INT *ncolfactor, *ncolupdate, *parent, *vtx2front;
  PORD_INT *first, *link, nfronts, J, K, u, vertex, vertex_root, count;
 /* Explicit shifting of indices to be optimized */
  for (u = nvtx; u >= 0; u--)
   {
     xadj_pe[u] = xadj_pe[u] - 1;
   }
   for (K = nedges-1; K >= 0; K--)
   {
      adjncy[K] = adjncy[K] - 1;
   }
   /* initialization of the graph */
   mymalloc(G, 1, graph_t);
   G->xadj   = xadj_pe;
   G->adjncy = adjncy;
   G->nvtx = nvtx;
   G->nedges = nedges;
   /* FIXME: G->vwght and G->tocwght accessed if G->type==UNWEIGHTED? */
   mymalloc(G->vwght, nvtx, PORD_INT);
   G->type = UNWEIGHTED;
   G->totvwght = nvtx;
   for (u = 0; u < nvtx; u++)
     G->vwght[u] = 1;
   /* main function of the Ordering */
   T = SPACE_ordering(G, options, cpus);
   nfronts = T->nfronts;
   ncolfactor = T->ncolfactor;
   ncolupdate = T->ncolupdate;
   parent = T->parent;
 /*    firstchild = T->firstchild; */
   vtx2front = T->vtx2front;
    /* -----------------------------------------------------------
     store the vertices/columns of a front in a bucket structure
     ----------------------------------------------------------- */
   mymalloc(first, nfronts, PORD_INT);
   mymalloc(link, nvtx, PORD_INT);
   for (J = 0; J < nfronts; J++)
      first[J] = -1;
   for (u = nvtx-1; u >= 0; u--)
      {
        J = vtx2front[u];
        link[u] = first[J];
        first[J] = u;
      }
  /* -----------------------------------------------------------
     fill the two arrays corresponding to the MUMPS tree structure
     ----------------------------------------------------------- */
   count = 0;
   for (K = firstPostorder(T); K != -1; K = nextPostorder(T, K))
     {
       vertex_root = first[K];
       if (vertex_root == -1)
         {
           /* Should never happen */
#          if defined(PORD_INTSIZE64) || defined(INTSIZE64)
           printf(" Internal error in mumps_pord, %ld\n",K);
#          else
           printf(" Internal error in mumps_pord, %d\n",K);
#          endif
           exit(-1);
         }
       /* for the principal column of the supervariable */
       if (parent[K] == -1)
         xadj_pe[vertex_root] = 0; /* root of the tree */
       else
         xadj_pe[vertex_root] = - (first[parent[K]]+1);
       nv[vertex_root] = ncolfactor[K] + ncolupdate[K];
       count++;
       for (vertex = link[vertex_root]; vertex != -1; vertex = link[vertex])
        /* for the secondary columns of the supervariable */
         {
           xadj_pe[vertex] = - (vertex_root+1);
           nv[vertex] = 0;
          count++;
         }
     }
  /* ----------------------
     free memory and return
     ---------------------- */
  free(first); free(link);
  free(G->vwght);
  free(G);
  freeElimTree(T);
  return (0);
}
/*********************************************************/
MUMPS_INT mumps_pord_wnd
(
        PORD_INT nvtx,
        PORD_INT nedges,
        PORD_INT *xadj_pe,
        PORD_INT *adjncy,
        PORD_INT *nv,
        PORD_INT *totw
)
{
/**********************************
Arguments:
input:
-----
- nvtx          : dimension of the Problem (N)
- nedges        : number of entries (NZ)
- adjncy        : non-zeros entries (IW input)
- totw          : sum of the weigth of the vertices
input/output:
-------------
- xadj_pe       : in: pointer through beginning of column non-zeros entries
                  out: "father array" (PE)
- nv            : in: weight of the vertices
                  out: "nfront array" (NV)
*************************************/
        graph_t    *G;
        elimtree_t *T;
        timings_t  cpus[12];
        options_t  options[] = { SPACE_ORDTYPE, SPACE_NODE_SELECTION1,
                    SPACE_NODE_SELECTION2, SPACE_NODE_SELECTION3,
                    SPACE_DOMAIN_SIZE, 0 };
        PORD_INT *ncolfactor, *ncolupdate, *parent, *vtx2front;
        PORD_INT *first, *link, nfronts, J, K, u, vertex, vertex_root, count;
 /* Explicit shifting of indices to be optimized */
        for (u = nvtx; u >= 0; u--)
        {
          xadj_pe[u] = xadj_pe[u] - 1;
        }
        for (K = nedges-1; K >= 0; K--)
        {
          adjncy[K] = adjncy[K] - 1;
        }
        /* initialization of the graph */
        mymalloc(G, 1, graph_t);
        G->xadj  = xadj_pe;
        G->adjncy= adjncy;
        G->nvtx = nvtx;
        G->nedges = nedges;
        G->type = WEIGHTED;
        G->totvwght = (*totw);
        /* FIXME: avoid allocation and do: G->vwght=nv; instead? */
        mymalloc(G->vwght, nvtx, PORD_INT);
        for (u = 0; u < nvtx; u++)
          G->vwght[u] = nv[u];
        /* main function of the Ordering */
        T = SPACE_ordering(G, options, cpus);
        nfronts = T->nfronts;
        ncolfactor = T->ncolfactor;
        ncolupdate = T->ncolupdate;
        parent = T->parent;
  /*    firstchild = T->firstchild; */
        vtx2front = T->vtx2front;
    /* -----------------------------------------------------------
     store the vertices/columns of a front in a bucket structure
     ----------------------------------------------------------- */
        mymalloc(first, nfronts, PORD_INT);
        mymalloc(link, nvtx, PORD_INT);
        for (J = 0; J < nfronts; J++)
          first[J] = -1;
        for (u = nvtx-1; u >= 0; u--)
        {
          J = vtx2front[u];
          link[u] = first[J];
          first[J] = u;
        }
  /* -----------------------------------------------------------
     fill the two arrays corresponding to the MUMPS tree structure
     ----------------------------------------------------------- */
  count = 0;
  for (K = firstPostorder(T); K != -1; K = nextPostorder(T, K))
     {
       vertex_root = first[K];
       if (vertex_root == -1)
         {
           /* Should never happen */
#          if defined(PORD_INTSIZE64) || defined(INTSIZE64)
           printf(" Internal error in mumps_pord, %ld\n",K);
#          else
           printf(" Internal error in mumps_pord, %d\n",K);
#          endif
           exit(-1);
         }
         /* for the principal column of the supervariable */
       if (parent[K] == -1)
         xadj_pe[vertex_root] = 0; /* root of the tree */
       else
         xadj_pe[vertex_root] = - (first[parent[K]]+1);
       nv[vertex_root] = ncolfactor[K] + ncolupdate[K];
       count++;
       for (vertex = link[vertex_root]; vertex != -1; vertex = link[vertex])
         /* for the secondary columns of the supervariable */
         {
           xadj_pe[vertex] = - (vertex_root+1);
           nv[vertex] = 0;
           count++;
         }
     }
  /* ----------------------
     free memory and return
     ---------------------- */
  free(first); free(link);
  free(G->vwght);
  free(G);
  freeElimTree(T);
  return (0);
}
#endif /* pord */
