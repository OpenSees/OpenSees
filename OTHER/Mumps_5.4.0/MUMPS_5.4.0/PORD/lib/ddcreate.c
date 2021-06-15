/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: ddcreate.c
/
/ author        J"urgen Schulze, University of Paderborn
/ created       00nov28
/
/ This file contains functions dealing with construction/coarsening
/ of a domain decomposition
/
******************************************************************************

Data type:  struct domdec
              graph_t *G;            pointer to graph object
              int     ndom;          number of domains
              int     domwght;       total weight of domains
              int     *vtype;        type of node (see comment below)
              int     *color;        color of node (GRAY, BLACK, or WHITE)
              int     cwght[3];      weights of GRAY, BLACK, WHITE partitions
              int     *map;          maps nodes to next coarser domain decomp.
              struct domdec *prev;   pointer to previous finer domain decomp.
              struct domdec *next;   pointer to next coarser domain decomp.
Comments:
  o Structure holds the domain decompositions constructed by the
    coarsening process; it also holds the colorings of the domain decomp.
    computed by the refinement process
  o vtype[v]: represents the status of a node in the domain decomposition
              0, iff status of v is unknown
              1, iff v is a domain vertex
              2, iff v is a multisector vertex
              3, iff multisec v is eliminated and now forms a domain
              4, iff multisec v is absorbed by another multisec/domain
Methods in lib/ddcreate.c:
- dd = newDomainDecomposition(int nvtx, int nedges);
    o Initial: ndom = domwght = 0,
               cwght[GRAY] = cwght[BLACK] = cwght[WHITE] = 0,
               and prev = next = NULL
- void freeDomainDecomposition(domdec_t *dd);
- void printDomainDecomposition(domdec_t *dd);
- void checkDomainDecomposition(domdec_t *dd);
- void buildInitialDomains(graph_t *G, int *vtxlist, int *vtype, int *rep);
    o determines initial domains according to the order of nodes in vtxlist;
      furthermore, it sets rep[u] = v for all multisecs u that are adjacent
      to only one domain v
    o on start vtype[u] = 0 for all 0 <= u < nvtx, on return
      vtype[u] = 1, iff u belongs to a domain (rep[u]=u => u is seed of domain)
      vtype[u] = 2, iff u belongs to a multisec (rep[u]=u => u is seed)
- void mergeMultisecs(graph_t *G, int *vtype, int *rep);
    o merges all adjacent multisecs that do not share a common domain
    o on return vtype[w] = 4, iff multisec w belongs to multisec cluster
      u = rep[w]
- dd = initialDomainDecomposition(graph_t *G, int *map, int *vtype, int *rep);
    o allocates memory for the initial domain decomposition of G by calling
      newDomainDecomposition and creates the domain decomposition according
      to the vectors vtype and rep; the map vector maps vertices of G onto
      vertices of dd
- dd = constructDomainDecomposition(graph_t *G, int *map);
    o constructs an initial domain decomposition for the graph G by calling
      the functions (a) buildInitialDomains
                    (b) mergeMultisecs
                    (c) initialDomainDecomposition
      vextor map identifies vertices of G in the domain decomposition
- void computePriorities(domdec_t *dd, int *msvtxlist, int *key, int scoretype);
    o computes for each multisec u in msvtxlist its priority key[u] according
      to the node selection strategy scoretype
- void eliminateMultisecs(domdec_t *dd, int *msvtxlist, int *rep);
    o eliminates multisecs according to their order in msvtxlist; furthermore,
      it sets rep[u] = v for all multisecs u that are adjacent to only one
      newly formed domain v
    o on return
      dd->vtype[u] = 1, iff u is a domain (rep[u] = u)
      dd->vtype[u] = 2, iff u is an uneliminated multisec (rep[u] = u)
      dd->vtype[u] = 3, iff u is an eliminated multisec (rep[u] = u)
      dd->vtype[u] = 4, iff multisec u is absorbed by new domain v = rep[u];
- void findIndMultisecs(domdec_t *dd, int *msvtxlist, int *rep);
    o searches all unelim./unabsorbed multisecs in msnvtxlist for
      indistinguishable multisecs; sets dd->vtype[u] = 4 and rep[u] = v, iff
      u, v are indistinguishable and v is the representative of u
- dd2 = coarserDomainDecomposition(domdec_t* dd1, int *rep);
    o allocates memory for the coarser domain decomposition by calling
      newDomainDecomposition and creates the domain decomposition according
      to the vectors dd1->vtype and rep; vector dd1->map identifies the
      vertices of dd1 in dd2
- void shrinkDomainDecomposition(domdec_t *dd, int scoretype);
    o shrinks dd according to a chosen node selection strategy by calling
      the functions (a) computePriorities
                    (b) eliminateMultisecs
                    (c) findIndMultisecs
                    (d) coarserDomainDecomposition
      the coarser domain decomposition is appended to dd via prev/next pointers

******************************************************************************/

#include <space.h>


/*****************************************************************************
******************************************************************************/
domdec_t*
newDomainDecomposition(PORD_INT nvtx, PORD_INT nedges)
{ domdec_t *dd;

  mymalloc(dd, 1, domdec_t);
  mymalloc(dd->vtype, nvtx, PORD_INT);
  mymalloc(dd->color, nvtx, PORD_INT);
  mymalloc(dd->map, nvtx, PORD_INT);

  dd->G = newGraph(nvtx, nedges);
  dd->ndom = dd->domwght = 0;
  dd->cwght[GRAY] = dd->cwght[BLACK] = dd->cwght[WHITE] = 0;
  dd->prev = dd->next = NULL;

  return(dd);
}


/*****************************************************************************
******************************************************************************/
void
freeDomainDecomposition(domdec_t *dd)
{
  freeGraph(dd->G);
  free(dd->vtype);
  free(dd->color);
  free(dd->map);
  free(dd);
}


/*****************************************************************************
******************************************************************************/
void
printDomainDecomposition(domdec_t *dd)
{ graph_t *G;
  PORD_INT   count, u, v, i, istart, istop;

  G = dd->G;
  printf("\n#nodes %d (#domains %d, weight %d), #edges %d, totvwght %d\n",
         G->nvtx, dd->ndom, dd->domwght, G->nedges >> 1, G->totvwght);
  printf("partition weights: S %d, B %d, W %d\n", dd->cwght[GRAY],
         dd->cwght[BLACK], dd->cwght[WHITE]);
  for (u = 0; u < G->nvtx; u++)
   { count = 0;
     printf("--- adjacency list of node %d (vtype %d, color %d, map %d\n",
            u, dd->vtype[u], dd->color[u], dd->map[u]);
     istart = G->xadj[u];
     istop = G->xadj[u+1];
     for (i = istart; i < istop; i++)
      { v = G->adjncy[i];
        printf("%5d (vtype %2d, color %2d)", v, dd->vtype[v], dd->color[v]);
        if ((++count % 3) == 0)
          printf("\n");
      }
     if ((count % 3) != 0)
       printf("\n");
   }
}


/*****************************************************************************
******************************************************************************/
void
checkDomainDecomposition(domdec_t *dd)
{ PORD_INT *xadj, *adjncy, *vwght, *vtype;
  PORD_INT err, nvtx, ndom, domwght, dom, multi, u, v, i, istart, istop;

  nvtx = dd->G->nvtx;
  xadj = dd->G->xadj;
  adjncy = dd->G->adjncy;
  vwght = dd->G->vwght;
  vtype = dd->vtype;

  err = FALSE;
  printf("checking domain decomposition (#nodes %d, #edges %d)\n",
         dd->G->nvtx, dd->G->nedges >> 1);

  ndom = domwght = 0;
  for (u = 0; u < nvtx; u++)
   { /* check node type */
     if ((vtype[u] != 1) && (vtype[u] != 2))
      { printf("ERROR: node %d is neither DOMAIN nor MULTISEC\n", u);
        err = TRUE;
      }
     /* count domains and sum up their weight */
     if (vtype[u] == 1)
      { ndom++;
        domwght += vwght[u];
      } 
     /* check number of neighboring domains and multisecs */
     dom = multi = 0;
     istart = xadj[u];
     istop = xadj[u+1];
     for (i = istart; i < istop; i++)
      { v = adjncy[i];
        if (vtype[v] == 1) dom++;
        if (vtype[v] == 2) multi++;
      }
     if ((vtype[u] == 1) && (dom > 0))
      { printf("ERROR: domain %d is adjacent to other domain\n", u);
        err = TRUE;
      }
     if ((vtype[u] == 2) && (dom < 2))
      { printf("ERROR: less than 2 domains adjacent to multisec node %d\n", u);
        err = TRUE;
      }
     if ((vtype[u] == 2) && (multi > 0))
      { printf("ERROR: multisec %d is adjacent to other multisec nodes\n", u);
        err = TRUE;
      }
   }
  /* check number and weight of domains */
  if ((ndom != dd->ndom) || (domwght != dd->domwght))
   { printf("ERROR: number/size (%d/%d) of domains does not match with those in"
            " domain decomp. (%d/%d)\n", ndom, domwght, dd->ndom, dd->domwght);
     err = TRUE;
   }
  if (err) quit();
}


/*****************************************************************************
******************************************************************************/
void
buildInitialDomains(graph_t *G, PORD_INT *vtxlist, PORD_INT *vtype, PORD_INT *rep)
{ PORD_INT *xadj, *adjncy;
  PORD_INT nvtx, u, v, w, i, j, jstart, jstop;

  xadj = G->xadj;
  adjncy = G->adjncy;
  nvtx = G->nvtx;

  /* --------------------------------------------------------------------
     determine initial domains according to the order of nodes in vtxlist
     -------------------------------------------------------------------- */
  for (i = 0; i < nvtx; i++)
   { u = vtxlist[i];
     if (vtype[u] == 0)
      { vtype[u] = 1;
        jstart = xadj[u];
        jstop = xadj[u+1];
        for (j = jstart; j < jstop; j++)
         { v = adjncy[j];
           vtype[v] = 2;
         }
      }
   }

  /* ------------------------------------------------------------
     eliminate all multisecs that are adjacent to only one domain
     ------------------------------------------------------------ */
  for (i = 0; i < nvtx; i++)
   { u = vtxlist[i];
     if (vtype[u] == 2)
      { v = -1;
        jstart = xadj[u];
        jstop = xadj[u+1];
        for (j = jstart; j < jstop; j++)
         { w = adjncy[j];
           if (vtype[w] == 1)
            { if (v == -1)
                v = rep[w];          /* u adjacent to domain v = rep[w] */
              else if (v != rep[w])
               { v = -1;             /* u adjacent to another domain */
                 break;
               }
            }
         }
        if (v != -1)                 /* u absorbed by domain v */
         { vtype[u] = 1;
           rep[u] = v;
         }
      }
   }
}


/*****************************************************************************
******************************************************************************/
void
mergeMultisecs(graph_t *G, PORD_INT *vtype, PORD_INT *rep)
{ PORD_INT *xadj, *adjncy, *tmp, *queue;
  PORD_INT nvtx, qhead, qtail, flag, keepon, u, v, w, x;
  PORD_INT i, istart, istop, j, jstart, jstop;

  nvtx = G->nvtx;
  xadj = G->xadj;
  adjncy = G->adjncy;

  /* ------------------------
     allocate working storage
     ------------------------ */
  mymalloc(tmp, nvtx, PORD_INT);
  mymalloc(queue, nvtx, PORD_INT);
  for (u = 0; u < nvtx; u++)
    tmp[u] = -1;

  /* -------------------------------------------------------
     merge all adjacent multisecs that do not share a domain
     ------------------------------------------------------- */
  flag = 1;
  for (u = 0; u < nvtx; u++)
    if (vtype[u] == 2)
     { qhead = 0; qtail = 1;
       queue[0] = u;
       vtype[u] = -2;

       /* multisec u is the seed of a new cluster, mark all adj. domains */
       istart = xadj[u];
       istop = xadj[u+1];
       for (i = istart; i < istop; i++)
        { v = adjncy[i];
          if (vtype[v] == 1)
            tmp[rep[v]] = flag;
        }

       /* and now build the cluster */
       while (qhead != qtail)
        { v = queue[qhead++];
          istart = xadj[v];
          istop = xadj[v+1];
          for (i = istart; i < istop; i++)
           { keepon = TRUE;
             w = adjncy[i];
             if (vtype[w] == 2)
              { jstart = xadj[w];
                jstop = xadj[w+1];
                for (j = jstart; j < jstop; j++)
                 { x = adjncy[j];
                   if ((vtype[x] == 1) && (tmp[rep[x]] == flag))
                    { keepon = FALSE;
                      break;
                    }
                 }
                if (keepon)
                /* multisecs v and w have no domain in common; mark */
                /* all domains adjacent to w and put w in cluster u */
                 { for (j = jstart; j < jstop; j++)
                    { x = adjncy[j];
                      if (vtype[x] == 1) tmp[rep[x]] = flag;
                     }
                    queue[qtail++] = w;
                    rep[w] = u;
                    vtype[w] = -2;
                 }
              }
           }
        }

       /* clear tmp vector for next round */
       flag++;
     }

  /* ------------------------------------
     reset vtype and free working storage
     ------------------------------------ */
  for (u = 0; u < nvtx; u++)
    if (vtype[u] == -2)
      vtype[u] = 2;
  free(tmp); free(queue);
}


/*****************************************************************************
******************************************************************************/
domdec_t*
initialDomainDecomposition(graph_t *G, PORD_INT *map, PORD_INT *vtype, PORD_INT *rep)
{ domdec_t *dd;
  PORD_INT      *xadj, *adjncy, *vwght, *xadjdd, *adjncydd, *vwghtdd, *vtypedd;
  PORD_INT      *tmp, *bin, nvtx, nedges, nvtxdd, nedgesdd, ndom, domwght, flag;
  PORD_INT      i, j, jstart, jstop, u, v, w;

  nvtx = G->nvtx;
  nedges = G->nedges;
  xadj = G->xadj;
  adjncy = G->adjncy;
  vwght = G->vwght;

  /* ------------------------
     allocate working storage
     ------------------------ */
  mymalloc(tmp, nvtx, PORD_INT);
  mymalloc(bin, nvtx, PORD_INT);
  for (u = 0; u < nvtx; u++)
   { tmp[u] = -1;
     bin[u] = -1;
   }

  /* -------------------------------------------------------------
     allocate memory for the dd using upper bounds nvtx and nedges
     ------------------------------------------------------------- */
  dd = newDomainDecomposition(nvtx, nedges);
  xadjdd = dd->G->xadj;
  adjncydd = dd->G->adjncy;
  vwghtdd = dd->G->vwght;
  vtypedd = dd->vtype;

  /* -------------------------------------------------------
     put all nodes u belonging to representative v in bin[v]
     ------------------------------------------------------- */
  for (u = 0; u < nvtx; u++)
   { v = rep[u];
     if (u != v)
      { bin[u] = bin[v];
        bin[v] = u;
      }
   }

  /* ----------------------------------------------
     and now build the initial domain decomposition
     ---------------------------------------------- */
  flag = 1;
  nedgesdd = nvtxdd = 0;
  ndom = domwght = 0;
  for (u = 0; u < nvtx; u++)
    if (rep[u] == u)
     { xadjdd[nvtxdd] = nedgesdd;
       vtypedd[nvtxdd] = vtype[u];
       vwghtdd[nvtxdd] = 0;
       tmp[u] = flag;

       /* find all cluster that are adjacent to u in dom. dec. */
       v = u;
       do
        { map[v] = nvtxdd;
          vwghtdd[nvtxdd] += vwght[v];
          jstart = xadj[v];
          jstop = xadj[v+1];
          for (j = jstart; j < jstop; j++)
           { w = adjncy[j];
             if ((vtype[w] != vtype[u]) && (tmp[rep[w]] != flag))
              { tmp[rep[w]] = flag;
                adjncydd[nedgesdd++] = rep[w];
              }
           }
          v = bin[v];
        } while (v != -1);

       if (vtypedd[nvtxdd] == 1)
        { ndom++;
          domwght += vwghtdd[nvtxdd];
        }
       nvtxdd++;
       flag++;
     }

  /* --------------------------------------------
     finalize the new domain decomposition object
     -------------------------------------------- */
  xadjdd[nvtxdd] = nedgesdd;
  dd->G->nvtx = nvtxdd;
  dd->G->nedges = nedgesdd;
  dd->G->type = WEIGHTED;
  dd->G->totvwght = G->totvwght;
  for (i = 0; i < nedgesdd; i++)
    adjncydd[i] = map[adjncydd[i]];
  for (u = 0; u < nvtxdd; u++)
    dd->color[u] = dd->map[u] = -1;
  dd->ndom = ndom;
  dd->domwght = domwght;

  /* -------------------------------
     free working storage and return
     ------------------------------- */
  free(tmp); free(bin);
  return(dd);
}


/*****************************************************************************
******************************************************************************/
domdec_t*
constructDomainDecomposition(graph_t *G, PORD_INT *map)
{ domdec_t *dd;
  PORD_INT      *xadj, *adjncy, *vwght, *vtxlist, *vtype, *key, *rep;
  PORD_INT      nvtx, deg, u, i, istart, istop;

  nvtx = G->nvtx;
  xadj = G->xadj;
  adjncy = G->adjncy;
  vwght = G->vwght;

  /* ---------------------------------------------------------
     sort the vertices in G in ascending order of their degree
     --------------------------------------------------------- */
  mymalloc(vtxlist, nvtx, PORD_INT);
  mymalloc(key, nvtx, PORD_INT);
  for (u = 0; u < nvtx; u++)
   { vtxlist[u] = u;
     istart = xadj[u];
     istop = xadj[u+1];
     switch(G->type)
      { case UNWEIGHTED:
          deg = istop - istart;
          break;
        case WEIGHTED:
          deg = 0;
          for (i = istart; i < istop; i++)
            deg += vwght[adjncy[i]];
          break;
        default:
          fprintf(stderr, "\nError in function constructDomainDecomposition\n"
               "  unrecognized graph type %d\n", G->type);
          quit();
      }
     key[u] = deg;
   }
  distributionCounting(nvtx, vtxlist, key);
  free(key);

   /* -------------------------------------------------------------
      build initial domains and cluster multisecs that do not share
      a common domain
      ------------------------------------------------------------- */
   mymalloc(vtype, nvtx, PORD_INT);
   mymalloc(rep, nvtx, PORD_INT);
   for (u = 0; u < nvtx; u++)
    { vtype[u] = 0;
      rep[u] = u;
    }
   buildInitialDomains(G, vtxlist, vtype, rep);
   mergeMultisecs(G, vtype, rep);
   free(vtxlist);

   /* --------------------------------------------------
      finally, build the domain decomposition and return
      -------------------------------------------------- */
   dd = initialDomainDecomposition(G, map, vtype, rep);
   free(vtype); free(rep);
   return(dd);
}


/*****************************************************************************
******************************************************************************/
void
computePriorities(domdec_t *dd, PORD_INT *msvtxlist, PORD_INT *key, PORD_INT scoretype)
{ PORD_INT *xadj, *adjncy, *vwght, *marker;
  PORD_INT nvtx, nlist, k, weight, deg, u, v, w;
  PORD_INT i, istart, istop, j, jstart, jstop;

  nvtx = dd->G->nvtx;
  xadj = dd->G->xadj;
  adjncy = dd->G->adjncy;
  vwght = dd->G->vwght;
  marker = dd->map;
  nlist = nvtx - dd->ndom;

  switch(scoretype)
   { case QMRDV:  /* maximal relative decrease of variables in quotient graph */
       for (k = 0; k < nlist; k++)
        { u = msvtxlist[k];
          weight = vwght[u];
          istart = xadj[u];
          istop = xadj[u+1];
          for (i = istart; i < istop; i++)
            weight += vwght[adjncy[i]];
          key[u] = weight / vwght[u];
        }
       break;

     case QMD:    /* ----------------------- minimum degree in quotient graph */
       for (k = 0; k < nlist; k++)
        { u = msvtxlist[k];
          marker[u] = -1;
        }
       for (k = 0; k < nlist; k++)
        { u = msvtxlist[k];
          marker[u] = u;
          deg = 0;
          istart = xadj[u];
          istop = xadj[u+1];
          for (i = istart; i < istop; i++)
           { v = adjncy[i];
             jstart = xadj[v];
             jstop = xadj[v+1];
             for (j = jstart; j < jstop; j++)
              { w = adjncy[j];
                if (marker[w] != u)
                 { marker[w] = u;
                   deg += vwght[w];
                 }
              }
           }
          key[u] = deg;
        }
       break;

     case QRAND:  /* ------------------------------------------------- random */
       for (k = 0; k < nlist; k++)
        { u = msvtxlist[k];
          key[u] = myrandom(nvtx);
        }
       break;

     default:
       fprintf(stderr, "\nError in internal function computePriorities\n"
            "  unrecognized node selection strategy %d\n", scoretype);
       quit();
   }
}


/*****************************************************************************
******************************************************************************/
void
eliminateMultisecs(domdec_t *dd, PORD_INT *msvtxlist, PORD_INT *rep)
{ PORD_INT *xadj, *adjncy, *vtype;
  PORD_INT nvtx, nlist, keepon, u, v, w, k, i, istart, istop;

  nvtx = dd->G->nvtx;
  xadj = dd->G->xadj;
  adjncy = dd->G->adjncy;
  vtype = dd->vtype;
  nlist = nvtx - dd->ndom;

  /* -------------------------------------------------------
     eliminate multisecs according to the order in msvtxlist
     ------------------------------------------------------- */
  for (k = 0; k < nlist; k++)
   { u = msvtxlist[k];
     istart = xadj[u];
     istop = xadj[u+1];
     keepon = TRUE;
     for (i = istart; i < istop; i++)
      { v = adjncy[i];
        if (rep[v] != v)        /* domain already absorbed by an eliminated */
         { keepon = FALSE;      /* multisec => multisec u cannot be deleted */
           break;
         }
      }
     if (keepon)
      { vtype[u] = 3;
        for (i = istart; i < istop; i++)
         { v = adjncy[i];
           rep[v] = u;
         }
      }
   }

  /* ------------------------------------------------------------
     eliminate all multisecs that are adjacent to only one domain
     ------------------------------------------------------------ */
  for (k = 0; k < nlist; k++)
   { u = msvtxlist[k];
     if (vtype[u] == 2)
      { v = -1;
        istart = xadj[u];
        istop = xadj[u+1];
        for (i = istart; i < istop; i++)
         { w = adjncy[i];
           if (v == -1)
             v = rep[w];          /* u adjacent to domain v = rep[w] */
           else if (v != rep[w])
            { v = -1;             /* u adjacent to another domain */
              break;
            }
         }
        if (v != -1)              /* u absorbed by domain v */
         { vtype[u] = 4;
           rep[u] = v;
         }
      }
   }
}


/*****************************************************************************
******************************************************************************/
void
findIndMultisecs(domdec_t *dd, PORD_INT *msvtxlist, PORD_INT *rep)
{ PORD_INT *xadj, *adjncy, *vtype, *tmp, *bin, *checksum, *next, *key;
  PORD_INT nvtx, nlist, flag, keepon, deg, chk, ulast, u, v, k, i, istart, istop;

  nvtx = dd->G->nvtx;
  xadj = dd->G->xadj;
  adjncy = dd->G->adjncy;
  vtype = dd->vtype;
  nlist = nvtx - dd->ndom;
  checksum = dd->map;

  /* ------------------------
     allocate working storage
     ------------------------ */
  mymalloc(tmp, nvtx, PORD_INT);
  mymalloc(bin, nvtx, PORD_INT);
  mymalloc(next, nvtx, PORD_INT);
  mymalloc(key, nvtx, PORD_INT);
  for (u = 0; u < nvtx; u++)
   { tmp[u] = -1;
     bin[u] = -1;
   }

  /* -------------------------------------------------------------------
     compute checksums for all unelim./unabsorbed multisecs in msvtxlist
     ------------------------------------------------------------------- */
  flag = 1;
  for (k = 0; k < nlist; k++)
   { u = msvtxlist[k];
     if (vtype[u] == 2)
      { deg = chk = 0;
        istart = xadj[u];
        istop = xadj[u+1];
        for (i = istart; i < istop; i++)
         { v = adjncy[i];
           if (tmp[rep[v]] != flag)
            { tmp[rep[v]] = flag;
              chk += rep[v];
              deg++;
            }
         }
        chk = chk % nvtx;
        checksum[u] = chk;
        key[u] = deg;
        next[u] = bin[chk];
        bin[chk] = u;
        flag++;
      }
   }

  /* ---------------------------------
     merge indistinguishable multisecs
     --------------------------------- */
  for (k = 0; k < nlist; k++)
   { u = msvtxlist[k];
     if (vtype[u] == 2)
      { chk = checksum[u];
        v = bin[chk];          /* examine all multisecs in bin[hash] */
        bin[chk] = -1;         /* do this only once */
        while (v != -1)
         { istart = xadj[v];
           istop = xadj[v+1];
           for (i = istart; i < istop; i++)
             tmp[rep[adjncy[i]]] = flag;
           ulast = v;          /* v is principal and u is a potiential */
           u = next[v];        /* nonprincipal variable */
           while (u != -1)
            { keepon = TRUE;
              if (key[u] != key[v])
                keepon = FALSE;
              if (keepon)
               { istart = xadj[u];
                 istop = xadj[u+1];
                 for (i = istart; i < istop; i++)
                   if (tmp[rep[adjncy[i]]] != flag)
                    { keepon = FALSE;
                      break;
                    }
               }
              if (keepon)           /* found it! mark u as nonprincipal */
               { rep[u] = v;
                 /* printf(" >> mapping %d onto %d\n", u, v); */
                 vtype[u] = 4;
                 u = next[u];
                 next[ulast] = u;  /* remove u from bin */
               }
              else                 /* failed */
               { ulast = u;
                 u = next[u];
               }
            }
           v = next[v];        /* no more variables can be absorbed by v */
           flag++;             /* clear tmp vector for next round */
         }
      }
   }

  /* --------------------
     free working storage
     -------------------- */
  free(tmp); free(bin);
  free(next); free(key);
}


/*****************************************************************************
******************************************************************************/
domdec_t*
coarserDomainDecomposition(domdec_t* dd1, PORD_INT *rep)
{ domdec_t *dd2;
  PORD_INT      *xadjdd1, *adjncydd1, *vwghtdd1, *vtypedd1, *mapdd1;
  PORD_INT      *xadjdd2, *adjncydd2, *vwghtdd2, *vtypedd2;
  PORD_INT      *tmp, *bin, nvtxdd1, nedgesdd1, nvtxdd2, nedgesdd2;
  PORD_INT      ndom, domwght, flag, u, v, w, i, istart, istop;

  nvtxdd1 = dd1->G->nvtx;
  nedgesdd1 = dd1->G->nedges;
  xadjdd1 = dd1->G->xadj;
  adjncydd1 = dd1->G->adjncy;
  vwghtdd1 = dd1->G->vwght;
  vtypedd1 = dd1->vtype;
  mapdd1 = dd1->map;

  /* ------------------------
     allocate working storage
     ------------------------ */
  mymalloc(tmp, nvtxdd1, PORD_INT);
  mymalloc(bin, nvtxdd1, PORD_INT);
  for (u = 0; u < nvtxdd1; u++)
   { tmp[u] = -1;
     bin[u] = -1;
   }

  /* ------------------------------------------------------------
     allocate memory using the upper bounds nvtxdd1 and nedgesdd1
     ------------------------------------------------------------ */
  dd2 = newDomainDecomposition(nvtxdd1, nedgesdd1);
  xadjdd2 = dd2->G->xadj;
  adjncydd2 = dd2->G->adjncy;
  vwghtdd2 = dd2->G->vwght;
  vtypedd2 = dd2->vtype;

  /* -------------------------------------------------------
     put all nodes u belonging to representative v in bin[v]
     ------------------------------------------------------- */
  for (u = 0; u < nvtxdd1; u++)
   { v = rep[u];
     if (u != v)
      { bin[u] = bin[v];
        bin[v] = u;
      }
   }

  /* ----------------------------------------------
     and now build the coarser domain decomposition
     ---------------------------------------------- */
  flag = 1;
  nvtxdd2 = nedgesdd2 = 0;
  ndom = domwght = 0;
  for (u = 0; u < nvtxdd1; u++)
    if (rep[u] == u)
     { xadjdd2[nvtxdd2] = nedgesdd2;
       vwghtdd2[nvtxdd2] = 0;
       vtypedd2[nvtxdd2] = vtypedd1[u];
       if (vtypedd2[nvtxdd2] == 3)
         vtypedd2[nvtxdd2] = 1;
       tmp[u] = flag;

       /* find all cluster that are adjacent to u in dom. dec. */
       v = u;
       do
        { mapdd1[v] = nvtxdd2;
          vwghtdd2[nvtxdd2] += vwghtdd1[v];
          if ((vtypedd1[v] == 1) || (vtypedd1[v] == 2))
           { istart = xadjdd1[v];
             istop = xadjdd1[v+1];
             for (i = istart; i < istop; i++)
              { w = adjncydd1[i];
                if (tmp[rep[w]] != flag)
                 { tmp[rep[w]] = flag;
                   adjncydd2[nedgesdd2++] = rep[w];
                 }
              }
           }
          v = bin[v];
        } while (v != -1);

       if (vtypedd2[nvtxdd2] == 1)
        { ndom++;
          domwght += vwghtdd2[nvtxdd2];
        }
       nvtxdd2++;
       flag++;
     }

  /* --------------------------------------------
     finalize the new domain decomposition object
     -------------------------------------------- */
  xadjdd2[nvtxdd2] = nedgesdd2;
  dd2->G->nvtx = nvtxdd2;
  dd2->G->nedges = nedgesdd2;
  dd2->G->type = WEIGHTED;
  dd2->G->totvwght = dd1->G->totvwght;
  for (i = 0; i < nedgesdd2; i++)
    adjncydd2[i] = mapdd1[adjncydd2[i]];
  for (u = 0; u < nvtxdd2; u++)
    dd2->color[u] = dd2->map[u] = -1;
  dd2->ndom = ndom;
  dd2->domwght = domwght;

  /* --------------------------
     set back node types in dd1
     -------------------------- */
  for (u = 0; u < nvtxdd1; u++)
    if ((vtypedd1[u] == 3) || (vtypedd1[u] == 4))
      vtypedd1[u] = 2;

  /* -------------------------------
     free working storage and return
     ------------------------------- */
  free(tmp); free(bin);
  return(dd2);
}


/*****************************************************************************
******************************************************************************/
void
shrinkDomainDecomposition(domdec_t* dd1, PORD_INT scoretype)
{ domdec_t *dd2;
  PORD_INT      *msvtxlist, *rep, *key;
  PORD_INT      nvtxdd1, nlist, u;

  nvtxdd1 = dd1->G->nvtx;
  mymalloc(msvtxlist, nvtxdd1, PORD_INT);
  mymalloc(rep, nvtxdd1, PORD_INT);
  mymalloc(key, nvtxdd1, PORD_INT);

  /* ---------------
     initializations
     --------------- */
  nlist = 0;
  for (u = 0; u < nvtxdd1; u++)
   { if (dd1->vtype[u] == 2)
       msvtxlist[nlist++] = u;
     rep[u] = u;
   }

  /* -------------------------------------
     compute priorities and sort multisecs
     ------------------------------------- */
  computePriorities(dd1, msvtxlist, key, scoretype);
  distributionCounting(nlist, msvtxlist, key);

  /* ----------------------------------------------------------
     eliminate multisecs and build coarser domain decomposition
     ---------------------------------------------------------- */
  eliminateMultisecs(dd1, msvtxlist, rep);
  findIndMultisecs(dd1, msvtxlist, rep);
  dd2 = coarserDomainDecomposition(dd1, rep);

  /* -----------------------------------
     append coarser domain decomposition
     ----------------------------------- */
  dd1->next = dd2;
  dd2->prev = dd1;

  free(msvtxlist);
  free(rep);
  free(key);
}
