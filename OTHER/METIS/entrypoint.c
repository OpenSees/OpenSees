/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * entrypoint.c 
 *
 * This file provides the entry points for the libmetis.a
 * library.
 *
 * Started 8/22/95
 * George
 *
 * $Id: entrypoint.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include <multilevel.h>

static EdgeType *__edges;
static VertexType *__vtxs;


/*************************************************************************
* This function is the entry point of the multilevel recursive bisection
* algorithm.
* 
* Input Parameters:
*   n			The number of vertices
*   (xadj,adjncy)	The adjacency structure of the graph in compressed 
*                       column format. xadj[n+1], adjncy[m] (m = # of edges).
*   vwgts		A vector of size n, that stores the weights of each
*                       vertex. It can be NULL, see weightflag.
*   ewgts		A vector of size m, that stores the weights of each edge.
*                       It can be NULL, see weightflag
*   weightflag		Used to indicate if the graph is weighted depending on
*                       its value:
*                         0 	No weights (vwgts=NULL, and edgeweights=NULL)
*                         1	Weights on edges (vwgts=NULL)
*			  2	Weights on vertices (edgeweights=NULL)
*			  3     Weights both on vertices and edges
*  nparts		The number of parts to partition the graph
*  options		This is an array of size 5 that is used to pass
*                       parameters for the various phases of the recursive
*                       bisection algorithm. If options[0] = 0, then the
*                       defaults parameters are used, otherwise if
*                       options[0] = 1, the meaning of the rest 4 integers
*                       in options is as follows:
*                         options[1] = CoarsenTo
*                         options[2] = MType
*                         options[3] = IPType
*                         options[4] = RType
*  numbering		If numbering=0 it indicates that array indices start 
*                       from 0 (as in C), and if numbering=1 it indicates that
*                       array numbering starts from 1 (as in Fortran). 
*                          
* Output Parameters:
*   edgecut		The edge-cut of the nparts-partition
*   partition		An array of size n, that stores the resulting partition
*                       vector.
**************************************************************************/
int PMETIS(int *n, int *xadj, int *adjncy, int *vwgts, int *ewgts, 
           int *weightflag, int *nparts, int *options, int *numbering,
           int *edgecut, int *partition)
{
  int i, j;
  CoarseGraphType graph;
  int *kpwgts, *cuts;

  /* 
   * Do some sanity checks
   */
  if (*n <= 0) 
    errexit("PMETIS: n is less or equal to zero!");
  if (xadj == NULL)
    errexit("PMETIS: xadj array is NULL!");
  if (adjncy == NULL)
    errexit("PMETIS: adjncy array is NULL!");
  if (*nparts <= 1)
    errexit("PMETIS: nparts is less or equal to 1!");
  if (*numbering != 0 && *numbering != 1)
    errexit("PMETIS: numbering is not 0 or 1!");
  if (partition == NULL)
    errexit("PMETIS: partition is NULL!");

  ConvertGraph(&graph, *n, xadj, adjncy, vwgts, ewgts, *weightflag, *numbering);
  InitRandom();
  PermuteGraphRandom(&graph);

  iset(*n, 0, partition);
  kpwgts = ismalloc(*n, 0, "PMETIS: kpwgts");
  cuts = ismalloc(2*(*n), -1, "PMETIS: cuts");

  if (options[0] == 0)
    *edgecut = MultiLevelPart(&graph, *nparts, PMETIS_CTO, PMETIS_MTYPE, PMETIS_IPTYPE, PMETIS_RTYPE,
                              0, *weightflag, partition, cuts, kpwgts);
  else
    *edgecut = MultiLevelPart(&graph, *nparts, options[1], options[2], options[3], options[4],
                              0, *weightflag, partition, cuts, kpwgts);

  free(kpwgts);
  free(cuts);
  CleanUpRootGraph();
}



/*************************************************************************
* This function is the entry point of the multilevel k-way partition
* algorithm.
* 
* Input Parameters:
*   n			The number of vertices
*   (xadj,adjncy)	The adjacency structure of the graph in compressed 
*                       column format. xadj[n+1], adjncy[m] (m = # of edges).
*   vwgts		A vector of size n, that stores the weights of each
*                       vertex. It can be NULL, see weightflag.
*   ewgts		A vector of size m, that stores the weights of each edge.
*                       It can be NULL, see weightflag
*   weightflag		Used to indicate if the graph is weighted depending on
*                       its value:
*                         0 	No weights (vwgts=NULL, and edgeweights=NULL)
*                         1	Weights on edges (vwgts=NULL)
*			  2	Weights on vertices (edgeweights=NULL)
*			  3     Weights both on vertices and edges
*  nparts		The number of parts to partition the graph
*  options		This is an array of size 5 that is used to pass
*                       parameters for the various phases of the recursive
*                       bisection algorithm. If options[0] = 0, then the
*                       defaults parameters are used, otherwise if
*                       options[0] = 1, the meaning of the rest 4 integers
*                       in options is as follows:
*                         options[1] = CoarsenTo
*                         options[2] = MType
*                         options[3] = IPType
*                         options[4] = RType
*  numbering		If numbering=0 it indicates that array indices start 
*                       from 0 (as in C), and if numbering=1 it indicates that
*                       array numbering starts from 1 (as in Fortran). 
*                          
* Output Parameters:
*   edgecut		The edge-cut of the nparts-partition
*   partition		An array of size n, that stores the resulting partition
*                       vector.
**************************************************************************/
int KMETIS(int *n, int *xadj, int *adjncy, int *vwgts, int *ewgts, 
           int *weightflag, int *nparts, int *options, int *numbering,
           int *edgecut, int *partition)
{
  int i, j;
  CoarseGraphType graph;
  int *kpwgts;

  /* 
   * Do some sanity checks
   */
  if (*n <= 0) 
    errexit("KMETIS: n is less or equal to zero!");
  if (xadj == NULL)
    errexit("KMETIS: xadj array is NULL!");
  if (adjncy == NULL)
    errexit("KMETIS: adjncy array is NULL!");
  if (*nparts <= 1)
    errexit("KMETIS: nparts is less or equal to 1!");
  if (*numbering != 0 && *numbering != 1)
    errexit("KMETIS: numbering is not 0 or 1!");
  if (partition == NULL)
    errexit("KMETIS: partition is NULL!");

  ConvertGraph(&graph, *n, xadj, adjncy, vwgts, ewgts, *weightflag, *numbering);
  InitRandom();
  PermuteGraphRandom(&graph);

  iset(*n, 0, partition);
  kpwgts = ismalloc(*n, 0, "KMETIS: kpwgts");

  if (options[0] == 0)
    *edgecut = KWayPart(&graph, *nparts, KMETIS_CTO, KMETIS_MTYPE, KMETIS_RTYPE,
                              0, *weightflag, partition, kpwgts);
  else
    *edgecut = KWayPart(&graph, *nparts, options[1], options[2], options[4],
                              0, *weightflag, partition, kpwgts);

  free(kpwgts);
  CleanUpRootGraph();
}


/*************************************************************************
* This function is the entry point of the multilevel nested disection
* ordering algorithm.
* 
* Input Parameters:
*   n			The number of vertices
*   (xadj,adjncy)	The adjacency structure of the graph in compressed 
*                       column format. xadj[n+1], adjncy[m] (m = # of edges).
*  options		This is an array of size 5 that is used to pass
*                       parameters for the various phases of the recursive
*                       bisection algorithm. If options[0] = 0, then the
*                       defaults parameters are used, otherwise if
*                       options[0] = 1, the meaning of the rest 4 integers
*                       in options is as follows:
*                         options[1] = CoarsenTo
*                         options[2] = MType
*                         options[3] = IPType
*                         options[4] = RType
*  numbering		If numbering=0 it indicates that array indices start 
*                       from 0 (as in C), and if numbering=1 it indicates that
*                       array numbering starts from 1 (as in Fortran). 
*                          
* Output Parameters:
*   perm		The permutation vector of size n.
*   iperm		The inverse permutation vector of size n.
**************************************************************************/
int OMETIS(int *n, int *xadj, int *adjncy, int *options, int *numbering, int *perm, int *iperm)
{
  int i, j;
  CoarseGraphType graph;
  SepNodeType *stree;

  /* 
   * Do some sanity checks
   */
  if (*n <= 0) 
    errexit("OMETIS: n is less or equal to zero!");
  if (xadj == NULL)
    errexit("OMETIS: xadj array is NULL!");
  if (adjncy == NULL)
    errexit("OMETIS: adjncy array is NULL!");
  if (*numbering != 0 && *numbering != 1)
    errexit("OMETIS: numbering is not 0 or 1!");
  if (perm == NULL)
    errexit("OMETIS: perm is NULL!");
  if (iperm == NULL)
    errexit("OMETIS: iperm is NULL!");

  ConvertGraph(&graph, *n, xadj, adjncy, NULL, NULL, 0, *numbering);
  InitRandom();
  PermuteGraphRandom(&graph);

  stree = (SepNodeType *)GKmalloc(sizeof(SepNodeType)*graph.nvtxs, "main: stree");
  for (i=0; i<graph.nvtxs; i++)
    stree[i].nvtxs = -1;

  if (options[0] == 0)
    MultiLevelOrder(&graph, OMETIS_CTO, OMETIS_MTYPE, OMETIS_IPTYPE, OMETIS_RTYPE,
                              0, iperm, stree);
  else
    MultiLevelOrder(&graph, options[1], options[2], options[3], options[4],
                              0, iperm, stree);

  if (*numbering == 0) {
    for (i=0; i<*n; i++)
      iperm[i]--;
    for (i=0; i<*n; i++)
      perm[iperm[i]] = i;
  }
  else {
    for (i=0; i<*n; i++)
      perm[iperm[i]-1] = i+1;
  }

  free(stree);
  CleanUpRootGraph();
}



/*************************************************************************
* This function converts a graph in the adjancency format to the one used
* internaly by METIS
**************************************************************************/
void ConvertGraph(CoarseGraphType *graph, int n, int *xadj, int *adjncy, 
                  int *vwgts, int *ewgts, int IsWeighted, int numbering)
{
  int i, j, k;
  int lastedge;
  int readew, readvw; 

  if (numbering == 1)
    for (i=0; i<n+1; i++)
      xadj[i]--;

  if (IsWeighted == 2 || IsWeighted == 3)
    readvw = 1;
  else
    readvw = 0;

  if (IsWeighted == 1 || IsWeighted == 3)
    readew = 1;
  else
    readew = 0;

  InitGraph(graph);

  graph->nvtxs = n;
  graph->nedges = xadj[n];

  graph->level = 0;
  graph->tvwgt = 0;
  graph->vtxs = (VertexType **)GKmalloc(sizeof(VertexType *)*graph->nvtxs, "readgraph: graph->vtxs");
  graph->label = imalloc(graph->nvtxs, "readgraph: graph->label");
  __edges = (EdgeType *)GKmalloc(graph->nedges*sizeof(EdgeType), "ConvertGraph: __edges");
  __vtxs = (VertexType *)GKmalloc(graph->nvtxs*sizeof(VertexType), "ConvertGraph: __vtxs");
  lastedge = 0;


  /* Start reading the graph file */
  for (i=0; i<n; i++) {
    graph->label[i] = i;
    graph->vtxs[i] = __vtxs + i;
    graph->vtxs[i]->edges = __edges + lastedge;
    graph->vtxs[i]->cewgt = 0;
    graph->vtxs[i]->ewgtsum = 0;
    graph->vtxs[i]->vwgt = (readvw ? vwgts[i] : 1);
    graph->tvwgt += graph->vtxs[i]->vwgt;
    graph->vtxs[i]->nedges = xadj[i+1] - xadj[i];


    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (numbering == 0)
        __edges[lastedge].edge = adjncy[j];
      else
        __edges[lastedge].edge = adjncy[j]-1;
      __edges[lastedge].ewgt = (readew ? ewgts[j] : 1);
      graph->vtxs[i]->ewgtsum += __edges[lastedge].ewgt;
      lastedge++;
    }
  }

  if (lastedge != graph->nedges)
    errexit("readgraph: Something wrong with the edges from input file %d %d",graph->nedges, lastedge);

  if (numbering == 1)
    for (i=0; i<n+1; i++)
      xadj[i]++;
}


/*************************************************************************
* This function frees whatever memory was left over
**************************************************************************/
void CleanUpRootGraph(void)
{
  free(__edges);
  free(__vtxs);
}

