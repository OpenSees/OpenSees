/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: graph.c
/
/ author        J"urgen Schulze, University of Paderborn
/ created       99sep14
/
/ This file contains functions dealing with the graph object.
/
******************************************************************************

Data type:  struct graph
              int  nvtx;       number of vertices
              int  nedges;     number of edges
              int  type;       vertices can be UNWEIGTHED or WEIGTHED
              int  totvwght;   total vertex weight
              int  *xadj;      xadj[u] points to start of u's adjacency list
              int  *adjncy;    holds the adjacency lists
              int  *vwght;     holds the vertex weights
Comments:
  o no edge weights are stored. In our application weighted graphs re-
    present compressed unweighted graphs and, therefore, ewght[(u,v)] =
    vwght[u] * vwght[v].
Methods in lib/graph.c:
- G = newGraph(int nvtx, int nedges);
    o Initial: we assume that G is unweighted, therefore:
               type = UNWEIGTHED, totvwght = nvtx, and vwght[u] = 1
- void freeGraph(graph_t *G);
- void printGraph(graph_t *G);
- void randomizeGraph(graph_t *G);
- Gsub = setupSubgraph(graph_t *G, int *intvertex, int nvint, int *vtxmap);
    o extracts the subgraph induced by the vertices in array intvertex from G.
      vtxmap maps the vertices in intvertex to the vertices of the subgraph.
- G = setupGraphFromMtx(inputMtx_t *A);
- G = setupGridGraph(int dimX, int dimY, int type);
    o type e {GRID, MESH, TORUS}
- int connectedComponents(graph_t *G);
- cG = compressGraph(graph_t *G, int *vtxmap)
    o cG = NULL, if there are not enough ind. vertices (see COMPRESS_FRACTION)
    o for u in G vtxmap[u] points to representative of u in cG

******************************************************************************/

#include <space.h>


/*****************************************************************************
******************************************************************************/
graph_t*
newGraph(PORD_INT nvtx, PORD_INT nedges)
{ graph_t *G;
  PORD_INT     i;

  mymalloc(G, 1, graph_t);
  mymalloc(G->xadj, (nvtx+1), PORD_INT);
  mymalloc(G->adjncy, nedges, PORD_INT);
  mymalloc(G->vwght, nvtx, PORD_INT);

  G->nvtx = nvtx;
  G->nedges = nedges;
  G->type = UNWEIGHTED;
  G->totvwght = nvtx;
  for (i = 0; i < nvtx; i++)
    G->vwght[i] = 1;

  return(G);
}


/*****************************************************************************
******************************************************************************/
void
freeGraph(graph_t *G)
{ 
  free(G->xadj);
  free(G->adjncy);
  free(G->vwght);
  free(G);
}


/*****************************************************************************
******************************************************************************/
void
printGraph(graph_t *G)
{ PORD_INT count, u, i, istart, istop;

  printf("\n#vertices %d, #edges %d, type %d, totvwght %d\n", G->nvtx,
         G->nedges >> 1, G->type, G->totvwght);
  for (u = 0; u < G->nvtx; u++)
   { count = 0;
     printf("--- adjacency list of vertex %d (weight %d):\n", u, G->vwght[u]);
     istart = G->xadj[u];
     istop = G->xadj[u+1];
     for (i = istart; i < istop; i++)
      { printf("%5d", G->adjncy[i]);
        if ((++count % 16) == 0)
          printf("\n");
      }
     if ((count % 16) != 0)
       printf("\n");
   }
}


/*****************************************************************************
******************************************************************************/
void
randomizeGraph(graph_t *G)
{ PORD_INT *xadj, *adjncy, nvtx, u, v, len, j, i, istart, istop;

  nvtx = G->nvtx;
  xadj = G->xadj;
  adjncy = G->adjncy;

  for (u = 0; u < nvtx; u++)
   { istart = xadj[u];
     istop = xadj[u+1];
     if ((len = istop - istart) > 1)
       for (i = istart; i < istop; i++)
        { j = myrandom(len);
          swap(adjncy[i], adjncy[i+j], v);
          len--;
        }
   }
}


/*****************************************************************************
******************************************************************************/
graph_t*
setupSubgraph(graph_t *G, PORD_INT *intvertex, PORD_INT nvint, PORD_INT *vtxmap)
{ graph_t *Gsub;
  PORD_INT     *xadj, *adjncy, *vwght, *xadjGsub, *adjncyGsub, *vwghtGsub;
  PORD_INT     nvtx, nedgesGsub, totvwght, u, v, i, j, jstart, jstop, ptr;

  nvtx = G->nvtx;
  xadj = G->xadj;
  adjncy = G->adjncy;
  vwght = G->vwght;

  /* -------------------------------------------------------------
     compute number of edges and local indices of vertices in Gsub
     ------------------------------------------------------------- */
  nedgesGsub = 0;
  for (i = 0; i < nvint; i++)
   { u = intvertex[i];
     if ((u < 0) || (u >= nvtx))
      { fprintf(stderr, "\nError in function setupSubgraph\n"
             "  node %d does not belong to graph\n", u);
        quit();
      }
     jstart = xadj[u];
     jstop = xadj[u+1];
     for (j = jstart; j < jstop; j++)
       vtxmap[adjncy[j]] = -1;
     nedgesGsub += (jstop - jstart);
   }
  for (i = 0; i < nvint; i++)
   { u = intvertex[i];
     vtxmap[u] = i;
   }

  Gsub = newGraph(nvint, nedgesGsub);
  xadjGsub = Gsub->xadj;
  adjncyGsub = Gsub->adjncy;
  vwghtGsub = Gsub->vwght;

  /* --------------------------
     build the induced subgraph
     -------------------------- */
  totvwght = 0; ptr = 0;
  for (i = 0; i < nvint; i++)
   { u = intvertex[i];
     xadjGsub[i] = ptr;
     vwghtGsub[i] = vwght[u];
     totvwght += vwght[u];
     jstart = xadj[u];
     jstop = xadj[u+1];
     for (j = jstart; j < jstop; j++)
      { v = adjncy[j];
        if (vtxmap[v] >= 0)
          adjncyGsub[ptr++] = vtxmap[v];
      }
   }
  xadjGsub[nvint] = ptr;
  Gsub->type = G->type;
  Gsub->totvwght = totvwght;
  return(Gsub);
}

/*****************************************************************************
******************************************************************************/
graph_t*
setupGraphFromMtx(inputMtx_t *A)
{ graph_t *G;
  PORD_INT     *xnza, *nzasub, *xadj, *adjncy;
  PORD_INT     neqs, nelem, nvtx, k, h1, h2, j, i, istart, istop;

  neqs = A->neqs;
  nelem = A->nelem;
  xnza = A->xnza;
  nzasub = A->nzasub;

  /* ------------------------------------
     allocate memory for unweighted graph
     ------------------------------------ */
  G = newGraph(neqs, 2*nelem);
  nvtx = G->nvtx;
  xadj = G->xadj;
  adjncy = G->adjncy;

  /* -----------------------------------------
     determine the size of each adjacency list
     ----------------------------------------- */
  for (k = 0; k < neqs; k++)
    xadj[k] = xnza[k+1] - xnza[k];
  for (k = 0; k < nelem; k++)
    xadj[nzasub[k]]++;

  /* -------------------------------------------------------------
     determine for each vertex where its adjacency list will start
     ------------------------------------------------------------- */
  h1 = xadj[0];
  xadj[0] = 0;
  for (k = 1; k <= nvtx; k++)
   { h2 = xadj[k];
     xadj[k] = xadj[k-1] + h1;
     h1 = h2;
   }

  /* ------------------------
     fill the adjacency lists
     ------------------------ */
  for (k = 0; k < neqs; k++)
   { istart = xnza[k];
     istop = xnza[k+1];
     for (i = istart; i < istop; i++)
      { j = nzasub[i];
        adjncy[xadj[k]++] = j;  /* store {k,j} in adjacency list of k */
        adjncy[xadj[j]++] = k;  /* store {j,k} in adjacency list of j */
      }
   }

  /* --------------------------------------------
     restore startpoint of each vertex and return
     -------------------------------------------- */
  for (k = nvtx-1; k > 0; k--)
    xadj[k] = xadj[k-1];
  xadj[0] = 0;
  return(G);
}


/*****************************************************************************
******************************************************************************/
graph_t*
setupGridGraph(PORD_INT dimX, PORD_INT dimY, PORD_INT type)
{ graph_t *G;
  PORD_INT     *xadj, *adjncy, nvtx, nedges, knz, k;

  /* ---------------
     initializations
     --------------- */
  G = NULL;
  knz = 0;
  nvtx = dimX * dimY;

  /* ---------------------------------
     create unweighted grid/mesh graph
     --------------------------------- */
  if ((type == GRID) || (type == MESH))
   { nedges = 8                           /* for edge vertices */
              + 6 * (dimX-2 + dimY-2)     /* for border vertices */
              + 4 * (dimX-2) * (dimY-2);  /* for interior vertices */
     if (type == MESH)
       nedges += 4 * (dimX-1) * (dimY-1); /* diagonals */

     G = newGraph(nvtx, nedges);
     xadj = G->xadj;
     adjncy = G->adjncy;

     for (k = 0; k < nvtx; k++)
      { xadj[k] = knz;
        if ((k+1) % dimX > 0)              /*   / k+1-dimX  (MESH) */
         { adjncy[knz++] = k+1;            /* k - k+1       (GRID) */
           if (type == MESH)               /*   \ k+1+dimX  (MESH) */
            { if (k+1+dimX < nvtx)
                adjncy[knz++] = k+1+dimX;
              if (k+1-dimX >= 0)
                adjncy[knz++] = k+1-dimX;
            }
         }
        if (k % dimX > 0)                  /* k-1-dimX \    (MESH) */
         { adjncy[knz++] = k-1;            /*      k-1 - k  (GRID) */
           if (type == MESH)               /* k-1+dimX /    (MESH) */
            { if (k-1+dimX < nvtx)
                adjncy[knz++] = k-1+dimX;
              if (k-1-dimX >= 0)
                adjncy[knz++] = k-1-dimX;
            }
         }
        if (k+dimX < nvtx)                 /* k-dimX  (GRID) */
          adjncy[knz++] = k+dimX;          /*    |           */
        if (k-dimX >= 0)                   /*    k           */
          adjncy[knz++] = k-dimX;          /*    |           */
      }                                    /* k+dimX  (GRID) */
     xadj[nvtx] = knz;
   }

  /* -----------------------------
     create unweighted torus graph
     ----------------------------- */
  if (type == TORUS)
   { nedges = 4 * dimX * dimY;

     G = newGraph(nvtx, nedges);
     xadj = G->xadj;
     adjncy = G->adjncy;

     for (k = 0; k < nvtx; k++)
      { xadj[k] = knz;
        if (((k+1) % dimX) == 0)           /* k -- k+1 */
          adjncy[knz++] = k+1-dimX;
        else
          adjncy[knz++] = k+1;
        if ((k % dimX) == 0)               /* k-1 -- k */
          adjncy[knz++] = k-1+dimX;
        else
          adjncy[knz++] = k-1;
        adjncy[knz++] = (k+dimX) % nvtx;               /* k-dimX */
        adjncy[knz++] =  (k+dimX*(dimY-1)) % nvtx;     /*    |   */
      }                                                /*    k   */
     xadj[nvtx] = knz;                                 /*    |   */
   }                                                   /* k+dimX */

  return(G);
}


/*****************************************************************************
******************************************************************************/
PORD_INT
connectedComponents(graph_t *G)
{ PORD_INT *xadj, *adjncy, *marker, *queue;
  PORD_INT nvtx, u, v, w, qhead, qtail, comp, i, istart, istop;

  nvtx = G->nvtx;
  xadj = G->xadj;
  adjncy = G->adjncy;
 
  /* ------------------------
     allocate working storage
     ------------------------ */ 
  mymalloc(marker, nvtx, PORD_INT);
  mymalloc(queue, nvtx, PORD_INT);

  /* ---------------
     initializations
     --------------- */
  comp = 0;
  for (u = 0; u < nvtx; u++)
    marker[u] = -1;

  /* --------------------------------------
     get the number of connected components
     -------------------------------------- */
  for (u = 0; u < nvtx; u++)
   if (marker[u] == -1)
    { comp++;
      qhead = 0; qtail = 1;
      queue[0] = u; marker[u] = 0;

      while (qhead != qtail)       /* breadth first search in each comp. */
       { v = queue[qhead++];
         istart = xadj[v];
         istop = xadj[v+1];
         for (i = istart; i < istop; i++)
          { w = adjncy[i];
            if (marker[w] == -1)
             { queue[qtail++] = w;
               marker[w] = 0;
             }
          }
       }
    }

  /* -------------------------------
     free working storage and return
     ------------------------------- */
  free(marker); free(queue);
  return(comp);
}


/*****************************************************************************
private function of compressGraph
******************************************************************************/
static PORD_INT
indNodes(graph_t *G, PORD_INT *vtxmap)
{ PORD_INT *xadj, *adjncy, *deg, *checksum, *tmp;
  PORD_INT nvtx, cnvtx, u, v, i, istart, istop, j, jstart, jstop;

  nvtx = G->nvtx;
  xadj = G->xadj;
  adjncy = G->adjncy;

  /* -------------------------
     set up the working arrays
     ------------------------- */
  mymalloc(deg, nvtx, PORD_INT);
  mymalloc(checksum, nvtx, PORD_INT);
  mymalloc(tmp, nvtx, PORD_INT);

  /* -------------------------------------------------
     compute for each vertex u its degree and checksum
     ------------------------------------------------- */
  for (u = 0; u < nvtx; u++)
   { istart = xadj[u];
     istop = xadj[u+1];
     deg[u] = istop - istart;
     checksum[u] = u;
     tmp[u] = -1;
     vtxmap[u] = u;
     for (i = istart; i < istop; i++)
       checksum[u] += adjncy[i];
   }

  /* -------------------------------------
     search for indistinguishable vertices
     ------------------------------------- */
  cnvtx = nvtx;
  for (u = 0; u < nvtx; u++)
   if (vtxmap[u] == u)
    { tmp[u] = u;
      istart = xadj[u];
      istop = xadj[u+1];
      for (i = istart; i < istop; i++)
        tmp[adjncy[i]] = u;

        /* scan adjacency list of vertex u for indistinguishable vertices */
        for (i = istart; i < istop; i++)
         { v = adjncy[i];
           if ((v > u) && (checksum[v] == checksum[u]) && (deg[v] == deg[u])
              && (vtxmap[v] == v))
            { jstart = xadj[v];
              jstop = xadj[v+1];
              for (j = jstart; j < jstop; j++)
                if (tmp[adjncy[j]] != u) goto FAILURE;

              /* found it!!! map v onto u */
              vtxmap[v] = u;
              cnvtx--;
FAILURE:      ;
            }
         }
    }

  /* ----------------------
     free memory and return
     ---------------------- */
  free(deg); free(checksum); free(tmp);
  return(cnvtx);
}


/*****************************************************************************
******************************************************************************/
graph_t*
compressGraph(graph_t* G, PORD_INT* vtxmap)
{ graph_t *Gc;
  PORD_INT     *xadj, *adjncy, *vwght, *xadjGc, *adjncyGc, *vwghtGc, *perm;
  PORD_INT     nvtx, nvtxGc, nedgesGc, u, v, i, istart, istop;

  nvtx = G->nvtx;
  xadj = G->xadj;
  adjncy = G->adjncy;
  vwght = G->vwght;

  /* --------------------------------------------------------------
     compressed graph small enough? if so, allocate working storage
     -------------------------------------------------------------- */
  /* avoid print statement
   * printf("indNodes(G, vtxmap) = %d",indNodes(G, vtxmap)); */
  if ((nvtxGc = indNodes(G, vtxmap)) > COMPRESS_FRACTION * nvtx)
    return(NULL);
  mymalloc(perm, nvtx, PORD_INT);

  /* -----------------------------------
     count edges of the compressed graph
     ----------------------------------- */
  nedgesGc = 0;
  for (u = 0; u < nvtx; u++)
   if (vtxmap[u] == u)
    { istart = xadj[u];
      istop = xadj[u+1];
      for (i = istart; i < istop; i++)
       { v = adjncy[i];
         if (vtxmap[v] == v) nedgesGc++;
       }
    }

  /* ---------------------------------------------------------
     allocate memory for the compressed graph and construct it
     --------------------------------------------------------- */
  Gc = newGraph(nvtxGc, nedgesGc);
  xadjGc = Gc->xadj;
  adjncyGc = Gc->adjncy;
  vwghtGc = Gc->vwght;

  nvtxGc = nedgesGc = 0;
  for (u = 0; u < nvtx; u++)
   if (vtxmap[u] == u)
    { istart = xadj[u];
      istop = xadj[u+1];
      xadjGc[nvtxGc] = nedgesGc;
      vwghtGc[nvtxGc] = 0;
      perm[u] = nvtxGc++;
      for (i = istart; i < istop; i++)
       { v = adjncy[i];
         if (vtxmap[v] == v) adjncyGc[nedgesGc++] = v;
       }
    }
  xadjGc[nvtxGc] = nedgesGc;

  for (i = 0; i < nedgesGc; i++)
    adjncyGc[i] = perm[adjncyGc[i]];
  for (u = 0; u < nvtx; u++)
   { vtxmap[u] = perm[vtxmap[u]];
     vwghtGc[vtxmap[u]] += vwght[u];
   }
  Gc->type = WEIGHTED;
  Gc->totvwght = G->totvwght;

  /* ----------------------
     free memory and return
     ---------------------- */
  free(perm);
  return(Gc);
}

