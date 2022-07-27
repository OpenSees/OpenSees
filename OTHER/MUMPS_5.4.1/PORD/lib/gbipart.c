/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: gbipart.c
/
/ author        J"urgen Schulze, University of Paderborn
/ created       00dec26
/
/ This file contains functions dealing with bipartite graphs
/
******************************************************************************

Data type:  struct gbipart
              graph_t *G;   pointer to graph object with E c X x Y
              int     nX;   the vertices 0,...,nX-1 belong to X
              int     nY;   the vertices nX,...,nX+nY-1 belong to Y
Comments:
  o Structure used to smooth a separator computed for a subgraph Gbisect.
    The separator is paired with the border vertices in black/white partition,
    thus, resulting in a bipartite graph.
Methods in lib/gbipart.c:
- Gbipart = newBipartiteGraph(int nX, int nY, int nedges);
- void freeBipartiteGraph(gbipart_t *Gbipart);
- void printGbipart(gbipart_t *Gbipart);
- Gbipart = setupBipartiteGraph(graph_t *G, int *bipartvertex, int nX, int nY,
                                int *vtxmap)
    o Gbipart is induced by the vertices in bipartvertex. The first
      nX vertices are the vertices 0...nX-1 and the last nY vertices
      are the vertices nX...nX+nY-1 of Gbipart. Vector vtxmap maps the
      vertices in bipartvertex to the vertices of the bipartite graph.
- void maximumMatching(gbipart_t *Gbipart, int *matching);
- void maximumFlow(gbipart_t *Gbipart, int *flow, int *rc)
    o flow[i] stores the flow over the edge in adjncy[i] of Gbipart. It is
      positive, if the edge is from X to Y, otherwise flow is negative.
    o rc[u] stores the residual capacity of edge (source,u), u e X,
      respectively (u,sink), u e Y. All edges between X and Y have
      infinite capacity, therefore, no rc value must be computed for them.
- void DMviaMatching(gbipart_t *Gbipart, int *matching, int *dmflag,
                     int *dmwght);
    o on return. vector dmflag is filled with the following values:
                    / SI, iff x e X is reachable via exposed node e X
       dmflag[x] = <  SX, iff x e X is reachable via exposed node e Y
                    \ SR, iff x e X - (SI u SX)
                    / BI, iff y e Y is reachable via exposed node e Y
       dmflag[y] = <  BX, iff y e Y is reachable via exposed node e X
                    \ BR, iff y e Y - (BI u BX)
    o on return, vector dmwght is filled with the following values:
         dmwght[SI] - weight of SI    dmwght[BI] - weight of BI
         dmwght[SX] - weight of SX    dmwght[BX] - weight of BX
         dmwght[SR] - weight of SR    dmwght[BR] - weight of BR
- void DMviaFlow(gbipart_t *Gbipart, int *flow, int *rc, int *dmflag,
                 int *dmwght);
    o vectors dmflag and dmwght are filled as described above

******************************************************************************/

#include <space.h>

#define FREE       -1
#define SOURCE     -2
#define SINK       -3


/*****************************************************************************
******************************************************************************/
gbipart_t*
newBipartiteGraph(PORD_INT nX, PORD_INT nY, PORD_INT nedges)
{ gbipart_t *Gbipart;

  mymalloc(Gbipart, 1, gbipart_t);
  Gbipart->G = newGraph(nX+nY, nedges);
  Gbipart->nX = nX;
  Gbipart->nY = nY;

  return(Gbipart);
}


/*****************************************************************************
******************************************************************************/
void
freeBipartiteGraph(gbipart_t *Gbipart)
{
  freeGraph(Gbipart->G);
  free(Gbipart);
}


/*****************************************************************************
******************************************************************************/
void
printGbipart(gbipart_t *Gbipart)
{ graph_t *G;
  PORD_INT     count, u, i, istart, istop;

  G = Gbipart->G;
  printf("\n#vertices %d (nX %d, nY %d), #edges %d, type %d, totvwght %d\n",
         G->nvtx, Gbipart->nX, Gbipart->nY, G->nedges >> 1, G->type,
         G->totvwght);
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
gbipart_t*
setupBipartiteGraph(graph_t *G, PORD_INT *bipartvertex, PORD_INT nX, PORD_INT nY, PORD_INT *vtxmap)
{ gbipart_t *Gbipart;
  PORD_INT       *xadj, *adjncy, *vwght, *xadjGb, *adjncyGb, *vwghtGb;
  PORD_INT       nvtx, nedgesGb, totvwght, u, x, y, i, j, jstart, jstop, ptr;

  nvtx = G->nvtx;
  xadj = G->xadj;
  adjncy = G->adjncy;
  vwght = G->vwght;

  /* ----------------------------------------------------------------
     compute number of edges and local indices of vertices in Gbipart
     ---------------------------------------------------------------- */
  nedgesGb = 0;
  for (i = 0; i < nX+nY; i++)
   { u = bipartvertex[i];
     if ((u < 0) || (u >= nvtx))
      { fprintf(stderr, "\nError in function setupBipartiteGraph\n"
             "  node %d does not belong to graph\n", u);
        quit();
      }
     jstart = xadj[u];
     jstop = xadj[u+1];
     for (j = jstart; j < jstop; j++)
       vtxmap[adjncy[j]] = -1;
     nedgesGb += (jstop - jstart);
   }
  for (i = 0; i < nX+nY; i++)
   { u = bipartvertex[i];
     vtxmap[u] = i;
   }

  Gbipart = newBipartiteGraph(nX, nY, nedgesGb);
  xadjGb = Gbipart->G->xadj;
  adjncyGb = Gbipart->G->adjncy;
  vwghtGb = Gbipart->G->vwght;

  /* ---------------------------------
     build the induced bipartite graph
     --------------------------------- */
  totvwght = 0; ptr = 0;
  for (i = 0; i < nX; i++)
   { x = bipartvertex[i];
     xadjGb[i] = ptr;
     vwghtGb[i] = vwght[x];
     totvwght += vwght[x];
     jstart = xadj[x];
     jstop = xadj[x+1];
     for (j = jstart; j < jstop; j++)
      { y = adjncy[j];
        if (vtxmap[y] >= nX)
          adjncyGb[ptr++] = vtxmap[y];
      }
   }
  for (i = nX; i < nX+nY; i++)
   { y = bipartvertex[i];
     xadjGb[i] = ptr;
     vwghtGb[i] = vwght[y];
     totvwght += vwght[y];
     jstart = xadj[y];
     jstop = xadj[y+1];
     for (j = jstart; j < jstop; j++)
      { x = adjncy[j];
        if ((vtxmap[x] >= 0) && (vtxmap[x] < nX))
          adjncyGb[ptr++] = vtxmap[x];
      }
   }
  xadjGb[nX+nY] = ptr;
  Gbipart->G->type = G->type;
  Gbipart->G->totvwght = totvwght;
  return(Gbipart);
} 


/*****************************************************************************
******************************************************************************/
void
maximumMatching(gbipart_t *Gbipart, PORD_INT *matching)
{ PORD_INT *xadj, *adjncy, *level, *marker, *queue, *stack;
  PORD_INT top, top2, u, x, x2, y, y2, nX, nY, i, istart, istop;
  PORD_INT qhead, qtail, max_level;

  xadj = Gbipart->G->xadj;
  adjncy = Gbipart->G->adjncy;
  nX = Gbipart->nX;
  nY = Gbipart->nY;

  mymalloc(level, (nX+nY), PORD_INT);
  mymalloc(marker, (nX+nY), PORD_INT);
  mymalloc(queue, nX, PORD_INT);
  mymalloc(stack, nY, PORD_INT);

  /* -------------------
     initialize matching
     ------------------- */
  for (u = 0; u < nX+nY; u++)
    matching[u] = FREE;

  /* ---------------------------------------------------
     construct maximal matching in bipartite graph (X,Y)
     --------------------------------------------------- */
  for (x = 0; x < nX; x++)
   { istart = xadj[x];
     istop = xadj[x+1];
     for (i = istart; i < istop; i++)
      { y = adjncy[i];
        if (matching[y] == FREE)
         { matching[x] = y;
           matching[y] = x;
           break;
         }
      }
   }

  /* --------------------------------------------------------------------
     construct maximum matching in bipartite graph (X,Y) (Hopcroft, Karp)
     -------------------------------------------------------------------- */
  while (TRUE)
   { for (u = 0; u < nX+nY; u++)
       level[u] = marker[u] = -1;
     qhead = qtail = 0;                /* fill queue with free X nodes */
     for (x = 0; x < nX; x++)
       if (matching[x] == FREE)
        { queue[qtail++] = x;
          level[x] = 0;
        }

     /* --------------------------------------------------------------
        breadth first search to construct layer network containing all
        vertex disjoint augmenting paths of minimal length
        -------------------------------------------------------------- */
     top = 0;
     max_level = MAX_INT;
     while (qhead != qtail)
      { x = queue[qhead++];                  /* note: queue contains only */
        if (level[x] < max_level)            /*       nodes from X        */
         { istart = xadj[x];
           istop = xadj[x+1];
           for (i = istart; i < istop; i++)
            { y = adjncy[i];
              if (level[y] == -1)
               { level[y] = level[x] + 1;
                 if (matching[y] == FREE)
                  { max_level = level[y];    /* note: stack contains only */
                    stack[top++] = y;        /*       nodes form Y        */
                  }
                 else if (level[y] < max_level)
                  { x2 = matching[y];
                    level[x2] = level[y] + 1;
                    queue[qtail++] = x2;
                  }
               }
            }
         }
      }
     if (top == 0) break;              /* no augmenting path found */ 

     /* ------------------------------------------------------------
        restricted depth first search to construct maximal number of
        vertex disjoint augmenting paths in layer network
        ------------------------------------------------------------ */
     while (top > 0)
      { top2 = top--;
        y = stack[top2-1];             /* get the next exposed node in Y */
        marker[y] = xadj[y];           /* points to next neighbor of y */

        while (top2 > top)
         { y = stack[top2-1];
           i = marker[y]++;
           if (i < xadj[y+1])          /* not all neighbors of y visited */
            { x = adjncy[i];
              if ((marker[x] == -1) && (level[x] == level[y]-1))
               { marker[x] = 0;
                 if (level[x] == 0)        /* augmenting path found */
                   while (top2 > top)      /* pop stack */
                    { y2 = stack[--top2];
                      x2 = matching[y2];   /*    / o == o        */
                      matching[x] = y2;    /*   /                */
                      matching[y2] = x;    /* x -- y2 == x2 -- y */
                      x = x2;              /*   \                */
                    }                      /*    \ o == o        */
                 else
                  { y2 = matching[x];
                    stack[top2++] = y2;
                    marker[y2] = xadj[y2];
                  }
               }
            }
           else top2--;
         }
      }
   }

  /* -------------------------------
     free working storage and return
     ------------------------------- */
  free(level); free(marker);
  free(queue); free(stack);
}


/*****************************************************************************
******************************************************************************/
void
maximumFlow(gbipart_t *Gbipart, PORD_INT *flow, PORD_INT *rc)
{ PORD_INT *xadj, *adjncy, *vwght, *parent, *marker, *queue;
  PORD_INT nedges, u, v, x, y, nX, nY, j, i, istart, istop;
  PORD_INT qhead, qtail, capacity;

  nedges = Gbipart->G->nedges;
  xadj = Gbipart->G->xadj;
  adjncy = Gbipart->G->adjncy;
  vwght = Gbipart->G->vwght;
  nX = Gbipart->nX;
  nY = Gbipart->nY;

  mymalloc(parent, (nX+nY), PORD_INT);
  mymalloc(marker, (nX+nY), PORD_INT);
  mymalloc(queue, (nX+nY), PORD_INT);

  /* -------------------------------------
     initialize flow and residual capacity
     ------------------------------------- */
  for (u = 0; u < nX+nY; u++)
    rc[u] = vwght[u];
  for (i = 0; i < nedges; i++)
    flow[i] = 0;

  /* --------------------------------------------------
     determine an initial flow in the bipartite network
     -------------------------------------------------- */
  for (x = 0; x < nX; x++)
   { istart = xadj[x];
     istop = xadj[x+1];
     for (i = istart; i < istop; i++)
      { y = adjncy[i];
        capacity = min(rc[x], rc[y]);
        if (capacity > 0)
         { rc[x] -= capacity;
           rc[y] -= capacity;
           flow[i] = capacity;
           for (j = xadj[y]; adjncy[j] != x; j++);
           flow[j] = -capacity;
         }
        if (rc[x] == 0) break;
      }
   }

  /* -----------------------------------------------------------
     construct maximum flow in bipartite network (Edmonds, Karp)
     ----------------------------------------------------------- */
  while (TRUE)
   { for (u = 0; u < nX+nY; u++)
       parent[u] = marker[u] = -1;
     qhead = qtail = 0;                /* fill queue with free X nodes */
     for (x = 0; x < nX; x++)
       if (rc[x] > 0)
        { queue[qtail++] = x;
          parent[x] = x;
        }

      /* ---------------------------------------------------------
         breadth first search to find the shortest augmenting path
         --------------------------------------------------------- */
      capacity = 0;
      while (qhead != qtail)
       { u = queue[qhead++];
         istart = xadj[u];
         istop = xadj[u+1];
         for (i = istart; i < istop; i++)
          { v = adjncy[i];
            if ((parent[v] == -1) && ((v >= nX) || (flow[i] < 0)))
            /* v >= nX => u->v is a forward edge having infty capacity     */
            /* otherwise  u<-v is a backward edge and (v,u) must have      */
            /*            positive capacity (i.e. (u,v) has neg. capacity) */
             { parent[v] = u;
               marker[v] = i;
               queue[qtail++] = v;
               if ((v >= nX) && (rc[v] > 0))  /* found it! */
                { u = v;                      /* (v,sink) is below capacity */
                  capacity = rc[u];
                  while (parent[u] != u)      /* get minimal residual capa. */
                   { i = marker[u];
                     u = parent[u];
                     if (u >= nX)
                       capacity = min(capacity, -flow[i]);
                   }
                  capacity = min(capacity, rc[u]);
                  rc[v] -= capacity;          /* augment flow by min. rc */
                  while (parent[v] != v)
                   { i = marker[v];
                     u = parent[v];
                     flow[i] += capacity;
                     for (j = xadj[v]; adjncy[j] != u; j++);
                     flow[j] = -flow[i];
                     v = u;
                   }
                  rc[v] -= capacity;
                  qhead = qtail;              /* escape inner while loop */
                  break;
                }
             }
          }
       }

     if (capacity == 0)
       break;
   }

  free(parent); free(marker);
  free(queue);
}


/*****************************************************************************
******************************************************************************/
void
DMviaMatching(gbipart_t *Gbipart, PORD_INT *matching, PORD_INT *dmflag, PORD_INT *dmwght)
{ PORD_INT *xadj, *adjncy, *vwght, *queue, qhead, qtail;
  PORD_INT u, x, nX, y, nY, i, istart, istop;

  xadj = Gbipart->G->xadj;
  adjncy = Gbipart->G->adjncy;
  vwght = Gbipart->G->vwght;
  nX = Gbipart->nX;
  nY = Gbipart->nY;

  mymalloc(queue, (nX+nY), PORD_INT);

  /* ----------------------------------------------------------------------
     mark all exposed nodes of X with SI and all exposed nodes of Y with BI
     ---------------------------------------------------------------------- */
  qhead = qtail = 0;
  for (x = 0; x < nX; x++)
    if (matching[x] == FREE)
     { queue[qtail++] = x;
       dmflag[x] = SI;
     }
    else dmflag[x] = SR;
  for (y = nX; y < nX+nY; y++)
    if (matching[y] == FREE)
     { queue[qtail++] = y;
       dmflag[y] = BI;
     }
    else dmflag[y] = BR;

  /* ------------------------------------------------------------------
     construct Dulmage-Mendelsohn decomp. starting with SI and BI nodes
     ------------------------------------------------------------------ */
  while (qhead != qtail)
   { u = queue[qhead++];
     istart = xadj[u];
     istop = xadj[u+1];
     switch(dmflag[u])
      { case SI:
          for (i = istart; i < istop; i++)
           { y = adjncy[i];
             if (dmflag[y] == BR)
              { queue[qtail++] = y;
                dmflag[y] = BX;
              }
           }
          break;
        case BX:
          x = matching[u];
          dmflag[x] = SI;
          queue[qtail++] = x;
          break;
        case BI:
          for (i = istart; i < istop; i++)
           { x = adjncy[i];
             if (dmflag[x] == SR)
              { queue[qtail++] = x;
                dmflag[x] = SX;
              }
           }
          break;
        case SX:
          y = matching[u];
          dmflag[y] = BI;
          queue[qtail++] = y;
          break;
      }
   }

  /* ----------------------
     fill the dmwght vector
     ---------------------- */
  dmwght[SI] = dmwght[SX] = dmwght[SR] = 0;
  for (x = 0; x < nX; x++)
    switch(dmflag[x])
     { case SI: dmwght[SI] += vwght[x]; break;
       case SX: dmwght[SX] += vwght[x]; break;
       case SR: dmwght[SR] += vwght[x]; break;
     }
  dmwght[BI] = dmwght[BX] = dmwght[BR] = 0;
  for (y = nX; y < nX+nY; y++)
    switch(dmflag[y])
     { case BI: dmwght[BI] += vwght[y]; break;
       case BX: dmwght[BX] += vwght[y]; break;
       case BR: dmwght[BR] += vwght[y]; break;
     }

  free(queue);
}


/*****************************************************************************
******************************************************************************/
void
DMviaFlow(gbipart_t *Gbipart, PORD_INT *flow, PORD_INT *rc, PORD_INT *dmflag, PORD_INT *dmwght)
{ PORD_INT *xadj, *adjncy, *vwght, *queue, qhead, qtail;
  PORD_INT u, v, x, nX, y, nY, i, istart, istop;

  xadj = Gbipart->G->xadj;
  adjncy = Gbipart->G->adjncy;
  vwght = Gbipart->G->vwght;
  nX = Gbipart->nX;
  nY = Gbipart->nY;

  mymalloc(queue, (nX+nY), PORD_INT);

  /* ----------------------------------------------------------
     mark all nodes reachable from source/sink with SOURCE/SINK
     ---------------------------------------------------------- */
  qhead = qtail = 0;
  for (x = 0; x < nX; x++)
    if (rc[x] > 0)
     { queue[qtail++] = x;
       dmflag[x] = SOURCE;
     }
    else dmflag[x] = FREE;
  for (y = nX; y < nX+nY; y++)
    if (rc[y] > 0)
     { queue[qtail++] = y;
       dmflag[y] = SINK;
     }
    else dmflag[y] = FREE;

  /* --------------------------------------------------------------------
     construct Dulmage-Mendelsohn decomp. starting with SOURCE/SINK nodes
     -------------------------------------------------------------------- */
  while (qhead != qtail)
   { u = queue[qhead++];
     istart = xadj[u];
     istop = xadj[u+1];
     switch(dmflag[u])
      { case SOURCE:
          for (i = istart; i < istop; i++)
           { v = adjncy[i];
             if ((dmflag[v] == FREE) && ((v >= nX) || (flow[i] < 0)))
              { queue[qtail++] = v;
                dmflag[v] = SOURCE;    /* v reachable via forward edge u->v */
              }                        /*         or via backward edge u<-v */
           }
          break;
        case SINK:
          for (i = istart; i < istop; i++)
           { v = adjncy[i];
             if ((dmflag[v] == FREE) && ((v < nX) || (flow[i] > 0)))
              { queue[qtail++] = v;
                dmflag[v] = SINK;    /* u reachable via forward edge v->u */
              }                      /*         or via backward edge v<-u */
           }
          break;
      }
   }

  /* -----------------------------------------------------
     all nodes x in X with dmflag[x] = SOURCE belong to SI
     all nodes x in X with dmflag[x] = SINK   belong to SX
     all nodes x in X with dmflag[x] = FREE   belong to SR
     ----------------------------------------------------- */
  dmwght[SI] = dmwght[SX] = dmwght[SR] = 0;
  for (x = 0; x < nX; x++)
    switch(dmflag[x])
     { case SOURCE: dmflag[x] = SI; dmwght[SI] += vwght[x]; break;
       case SINK:   dmflag[x] = SX; dmwght[SX] += vwght[x]; break;
       default:     dmflag[x] = SR; dmwght[SR] += vwght[x];
     }

  /* -----------------------------------------------------
     all nodes y in Y with dmflag[y] = SOURCE belong to BX
     all nodes y in Y with dmflag[y] = SINK   belong to BI
     all nodes y in Y with dmflag[y] = FREE   belong to BR
     ----------------------------------------------------- */
  dmwght[BI] = dmwght[BX] = dmwght[BR] = 0;
  for (y = nX; y < nX+nY; y++)
    switch(dmflag[y])
     { case SOURCE: dmflag[y] = BX; dmwght[BX] += vwght[y]; break;
       case SINK:   dmflag[y] = BI; dmwght[BI] += vwght[y]; break;
       default:     dmflag[y] = BR; dmwght[BR] += vwght[y];
     }

  free(queue);
}

