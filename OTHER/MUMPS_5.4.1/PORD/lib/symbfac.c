/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: symbfac.c
/
/ author        J"urgen Schulze, University of Paderborn
/ created       09/15/99
/
/ This file contains code for the symbolical factorization.
/
******************************************************************************

Data type:  struct css
              int   neqs;       number of equations
              int   nind;       number of row subscripts in compressed format
              int   owned;      does the object own vector nzlsub?
              int   *xnzl;      start of column
              int   *nzlsub;    row subscripts
              int   *xnzlsub;   start of column's row subscripts

            struct frontsub
              elimtree_t *PTP;  permuted elimination tree
              int   nind        number of indices
              int   *xnzf;      start of front subscripts
              int   *nzfsub     front subscripts for permuted elimtree PTP

            struct factorMtx
              int   nelem;      number of nonzeros (incl. diagonal entries)
              int   *perm;      permutation vector
              FLOAT *nzl;       vector of nonzeros (incl. diagonal entries)
              css_t *css;       compressed subscript structure of factorMtx
              frontsub_t *frontsub;  front subscripts
Comments:
Methods in lib/symbfac.c:
- css = newCSS(int neqs, int nind, int owned);
- void freeCSS(css_t *css);
- css = setupCSSFromGraph(graph_t *G, int *perm, int *invp);
- css = setupCSSFromFrontSubscripts(frontsub_t *frontsub);
- frontsub = newFrontSubscripts(elimtree_t *PTP);
- void freeFrontSubscripts(frontsub_t *frontsub);
- void printFrontSubscripts(frontsub_t *frontsub);
- frontsub = setupFrontSubscripts(elimtree_t *PTP, inputMtx_t *PAP);
- L = newFactorMtx(int nelem);
- void freeFactorMtx(factorMtx_t *L);
- void printFactorMtx(factorMtx_t *L);
- void initFactorMtx(factorMtx_t *L, inputMtx_t *PAP);
- void initFactorMtxNEW(factorMtx_t *L, inputMtx_t *PAP);

******************************************************************************/

#include <space.h>


/*****************************************************************************
******************************************************************************/
css_t*
newCSS(PORD_INT neqs, PORD_INT nind, PORD_INT owned)
{ css_t *css;

  mymalloc(css, 1, css_t);
  mymalloc(css->xnzl, (neqs+1), PORD_INT);
  mymalloc(css->xnzlsub, neqs, PORD_INT);
  if (owned)
   { mymalloc(css->nzlsub, nind, PORD_INT); }
  else
   { css->nzlsub = NULL; }
  css->neqs = neqs;
  css->nind = nind;
  css->owned = owned;

  return(css);
}


/*****************************************************************************
******************************************************************************/
void
freeCSS(css_t *css)
{
  free(css->xnzl);
  free(css->xnzlsub);
  if (css->owned)
    free(css->nzlsub);
  free(css);
}


/*****************************************************************************
******************************************************************************/
css_t*
setupCSSFromGraph(graph_t *G, PORD_INT *perm, PORD_INT *invp)
{ css_t *css;
  PORD_INT   *marker, *mergelink, *indices, *tmp, *xnzl, *xnzlsub, *nzlsub;
  PORD_INT   neqs, maxmem, u, v, col, mergecol, knz, mrk, beg, end;
  PORD_INT   fast, len, k, p, e, i, istart, istop;

  neqs = G->nvtx;
  maxmem = 2 * neqs;

  /* -------------------------
     set up the working arrays
     ------------------------- */
  mymalloc(marker, neqs, PORD_INT);
  mymalloc(indices, neqs, PORD_INT);
  mymalloc(mergelink, neqs, PORD_INT);
  mymalloc(tmp, neqs, PORD_INT);
  for (k = 0; k < neqs; k++)
    marker[k] = mergelink[k] = -1;

  /* -------------------------------------------------------
     allocate storage for the compressed subscript structure
     ------------------------------------------------------- */
  css = newCSS(neqs, maxmem, TRUE);

  xnzl = css->xnzl;
  nzlsub = css->nzlsub;
  xnzlsub = css->xnzlsub;

  /* ------------------------------------------------------------
     main loop: determine the subdiag. row indices of each column
     ------------------------------------------------------------ */
  xnzl[0] = 0;
  beg = end = 0;
  for (k = 0; k < neqs; k++)
   { indices[0] = k;
     knz = 1;

     if ((mergecol = mergelink[k]) != -1)        /* is k a leaf ??? */
      { mrk = marker[mergecol];
        fast = TRUE;
      }
     else
      { mrk = k;
        fast = FALSE;
      }

     /* --------------------------
        original columns (indices)
        -------------------------- */
     u = invp[k];
     istart = G->xadj[u];
     istop = G->xadj[u+1];
     for (i = istart; i < istop; i++)
      { v = G->adjncy[i];
        if ((col = perm[v]) > k)
         { indices[knz++] = col;
           if (marker[col] != mrk) fast = FALSE;
         }
      }

     /* --------------------------
        external columns (indices)
        -------------------------- */
     if ((fast) && (mergelink[mergecol] == -1))
      { xnzlsub[k] = xnzlsub[mergecol] + 1;
        knz = xnzl[mergecol+1] - xnzl[mergecol] - 1;
      }
     else
      { for (i = 0; i < knz; i++)
          marker[indices[i]] = k;
        while (mergecol != -1)
         { len = xnzl[mergecol+1] - xnzl[mergecol];
           istart = xnzlsub[mergecol];
           istop = istart + len;
           for (i = istart; i < istop; i++)
            { col = nzlsub[i];
              if ((col > k) && (marker[col] != k))
               { indices[knz++] = col;
                 marker[col] = k;
               }
            }
           mergecol = mergelink[mergecol];
         }
        qsortUpInts(knz, indices, tmp);

        /* ---------------------------------------------------
           store indices in nzlsub; resize nzlsub if too small
           --------------------------------------------------- */
        beg = end;
        xnzlsub[k] = beg;
        end = beg + knz;
        if (end > maxmem)
         { maxmem += neqs;
           myrealloc(nzlsub, maxmem, PORD_INT);
         }
        len = 0;
        for (i = beg; i < end; i++)
          nzlsub[i] = indices[len++];
      }

     /* ----------------------------
        append column k to mergelink
        ---------------------------- */
     if (knz > 1)
      { p = xnzlsub[k]+1;
        e = nzlsub[p];
        mergelink[k] = mergelink[e];
        mergelink[e] = k;
      }
     xnzl[k+1] = xnzl[k] + knz;
   }

  /* -----------------------------
     end of main loop: free memory
     ----------------------------- */
  free(marker); free(indices);
  free(tmp); free(mergelink);

  /* ------------------------------------------------------
     finalize the compressed subscript structure and return
     ------------------------------------------------------ */
  css->nind = xnzlsub[neqs-1] + 1;
  myrealloc(nzlsub, css->nind, PORD_INT);
  css->nzlsub = nzlsub;
  return(css);
}


/*****************************************************************************
******************************************************************************/
css_t*
setupCSSFromFrontSubscripts(frontsub_t *frontsub)
{ elimtree_t *PTP;
  css_t      *css;
  PORD_INT        *xnzf, *nzfsub, *ncolfactor, *xnzl, *xnzlsub;
  PORD_INT        nind, nvtx, K, beg, knz, firstcol, col;

  PTP = frontsub->PTP;
  xnzf = frontsub->xnzf;
  nzfsub = frontsub->nzfsub;
  nind = frontsub->nind;

  nvtx = PTP->nvtx;
  ncolfactor = PTP->ncolfactor;

  /* -------------------------------------------------------
     allocate storage for the compressed subscript structure
     ------------------------------------------------------- */
  css = newCSS(nvtx, nind, FALSE);
  css->nzlsub = nzfsub;

  xnzl = css->xnzl;
  xnzlsub = css->xnzlsub;
 
  /* ---------------------------------------
     fill the compressed subscript structure
     --------------------------------------- */ 
  xnzl[0] = 0;
  for (K = firstPostorder(PTP); K != -1; K = nextPostorder(PTP, K))
   { beg = xnzf[K];
     knz = xnzf[K+1] - beg;
     firstcol = nzfsub[beg];
     for (col = firstcol; col < firstcol + ncolfactor[K]; col++)
      { xnzlsub[col] = beg++;
        xnzl[col+1] = xnzl[col] + knz--;
      }
   }

  return(css);
}


/*****************************************************************************
******************************************************************************/
frontsub_t*
newFrontSubscripts(elimtree_t *PTP)
{ frontsub_t *frontsub;
  PORD_INT        nfronts, nind;

  nfronts = PTP->nfronts;
  nind = nFactorIndices(PTP);

  mymalloc(frontsub, 1, frontsub_t);
  mymalloc(frontsub->xnzf, (nfronts+1), PORD_INT);
  mymalloc(frontsub->nzfsub, nind, PORD_INT);

  frontsub->PTP = PTP;
  frontsub->nind = nind;

  return(frontsub);
}


/*****************************************************************************
******************************************************************************/
void
freeFrontSubscripts(frontsub_t *frontsub)
{
  freeElimTree(frontsub->PTP);
  free(frontsub->xnzf);
  free(frontsub->nzfsub);
  free(frontsub);
}


/*****************************************************************************
******************************************************************************/
void
printFrontSubscripts(frontsub_t *frontsub)
{ elimtree_t *PTP;
  PORD_INT        *xnzf, *nzfsub, *ncolfactor, *ncolupdate, *parent;
  PORD_INT        nfronts, root, K, count, i, istart, istop;

  PTP = frontsub->PTP;
  xnzf = frontsub->xnzf;
  nzfsub = frontsub->nzfsub;

  nfronts = PTP->nfronts;
  root = PTP->root;
  ncolfactor = PTP->ncolfactor;
  ncolupdate = PTP->ncolupdate;
  parent = PTP->parent;

  printf("#fronts %d, root %d\n", nfronts, root);
  for (K = firstPostorder(PTP); K != -1; K = nextPostorder(PTP, K))
   { printf("--- front %d, ncolfactor %d, ncolupdate %d, parent %d\n",
            K, ncolfactor[K], ncolupdate[K], parent[K]);
     count = 0;
     istart = xnzf[K];
     istop = xnzf[K+1];
     for (i = istart; i < istop; i++)
      { printf("%5d", nzfsub[i]);
        if ((++count % 16) == 0)
          printf("\n");
      }
     if ((count % 16) != 0)
       printf("\n");
   }
}


/*****************************************************************************
******************************************************************************/
frontsub_t*
setupFrontSubscripts(elimtree_t *PTP, inputMtx_t *PAP)
{ frontsub_t *frontsub;
  PORD_INT        *ncolfactor, *ncolupdate, *firstchild, *silbings, *vtx2front;
  PORD_INT        *xnza, *nzasub, *xnzf, *nzfsub;
  PORD_INT        *marker, *tmp, *first, *indices;
  PORD_INT        nvtx, nfronts, col, firstcol, knz;
  PORD_INT        u, i, istart, istop, K, J;

  nvtx = PTP->nvtx;
  nfronts = PTP->nfronts;
  ncolfactor = PTP->ncolfactor;
  ncolupdate = PTP->ncolupdate;
  firstchild = PTP->firstchild;
  silbings = PTP->silbings;
  vtx2front = PTP->vtx2front;

  xnza = PAP->xnza;
  nzasub = PAP->nzasub;

  /* -------------------------
     set up the working arrays
     ------------------------- */
  mymalloc(marker, nvtx, PORD_INT);
  mymalloc(tmp, nvtx, PORD_INT);
  mymalloc(first, nfronts, PORD_INT);
  for (i = 0; i < nvtx; i++)
    marker[i] = -1;

  /* --------------------------------
     find the first column of a front
     -------------------------------- */
  for (u = nvtx-1; u >= 0; u--)
   { K = vtx2front[u];
     first[K] = u;
   }

  /* -----------------------------------------
     allocate storage for the front subscripts
     ----------------------------------------- */
  frontsub = newFrontSubscripts(PTP);
  xnzf = frontsub->xnzf;
  nzfsub = frontsub->nzfsub;

  knz = 0;
  for (K = 0; K < nfronts; K++)
   { xnzf[K] = knz;
     knz += (ncolfactor[K] + ncolupdate[K]);
   }
  xnzf[K] = knz;

  /* -------------------------------------------
     postorder traversal of the elimination tree
     ------------------------------------------- */
  for (K = firstPostorder(PTP); K != -1; K = nextPostorder(PTP, K))
   { knz = 0;
     indices = nzfsub + xnzf[K];
     firstcol = first[K];

     /* -------------------------------------
        internal columns (indices) of front K
        ------------------------------------- */
     for (col = firstcol; col < firstcol + ncolfactor[K]; col++)
      { indices[knz++] = col;
        marker[col] = K;
      }

     /* -------------------------------------
        external columns (indices) of front K
        ------------------------------------- */
     for (J = firstchild[K]; J != -1; J = silbings[J])
      { istart = xnzf[J];
        istop = xnzf[J+1];
        for (i = istart; i < istop; i++)
         { col = nzfsub[i];
           if ((col > firstcol) && (marker[col] != K))
            { marker[col] = K;
              indices[knz++] = col;
            }
         }
      }

     /* -------------------------------------
        original columns (indices) of front K
        ------------------------------------- */
     for (u = 0; u < ncolfactor[K]; u++)
      { istart = xnza[firstcol + u];
        istop = xnza[firstcol + u + 1];
        for (i = istart; i < istop; i++)
         { col = nzasub[i];
           if ((col > firstcol) && (marker[col] != K))
            { marker[col] = K;
              indices[knz++] = col;
            }
         }
      }

     /* ----------------
        sort the indices
        ---------------- */
     qsortUpInts(knz, indices, tmp);
   }

  /* ----------------------
     free memory and return
     ---------------------- */
  free(marker); free(tmp); free(first);
  return(frontsub);
}


/*****************************************************************************
******************************************************************************/
factorMtx_t*
newFactorMtx(PORD_INT nelem)
{ factorMtx_t *L;

  mymalloc(L, 1, factorMtx_t);
  mymalloc(L->nzl, nelem, FLOAT);

  L->nelem = nelem;
  L->css = NULL;
  L->frontsub = NULL;
  L->perm = NULL;

  return(L);
}


/*****************************************************************************
******************************************************************************/
void
freeFactorMtx(factorMtx_t *L)
{
  freeCSS(L->css);
  freeFrontSubscripts(L->frontsub);
  free(L->nzl);
  free(L->perm);
  free(L);
}


/*****************************************************************************
******************************************************************************/
void
printFactorMtx(factorMtx_t *L)
{ css_t *css;
  FLOAT *nzl;
  PORD_INT   *xnzl, *nzlsub, *xnzlsub;
  PORD_INT   neqs, nelem, nind, k, ksub, i, istart, istop;

  nelem = L->nelem;
  nzl = L->nzl;
  css = L->css;

  neqs = css->neqs;
  nind = css->nind;
  xnzl = css->xnzl;
  nzlsub = css->nzlsub;
  xnzlsub = css->xnzlsub;

  printf("#equations %d, #elements (+diag.) %d, #indices (+diag.) %d\n",
         neqs, nelem, nind);
  for (k = 0; k < neqs; k++)
   { printf("--- column %d\n", k);
     ksub = xnzlsub[k];
     istart = xnzl[k];
     istop = xnzl[k+1];
     for (i = istart; i < istop; i++)
       printf("  row %5d, entry %e\n", nzlsub[ksub++], nzl[i]);
   }
}


/*****************************************************************************
******************************************************************************/
void
initFactorMtx(factorMtx_t *L, inputMtx_t *PAP)
{ elimtree_t *PTP;
  frontsub_t *frontsub;
  css_t      *css;
  PORD_INT        *ncolfactor;
  FLOAT      *nzl, *nza, *diag;
  PORD_INT        *xnzl, *nzlsub, *xnzlsub, *xnza, *nzasub, *xnzf, *nzfsub;
  PORD_INT        nelem, K, k, kstart, h, hstart, dis, i, istart, istop;
  PORD_INT        firstcol, lastcol;
 
  nelem = L->nelem; 
  nzl = L->nzl;
  css = L->css;
  xnzl = css->xnzl;
  nzlsub = css->nzlsub;
  xnzlsub = css->xnzlsub;

  frontsub = L->frontsub;
  PTP = frontsub->PTP;
  ncolfactor = PTP->ncolfactor;
  xnzf = frontsub->xnzf;
  nzfsub = frontsub->nzfsub;

  diag = PAP->diag;
  nza = PAP->nza;
  xnza = PAP->xnza;
  nzasub = PAP->nzasub;

  /* ------------------------------------
     set all numerical values of L to 0.0
     ------------------------------------ */
  for (i = 0; i < nelem; i++)
    nzl[i] = 0.0;

  /* --------------------------------------------
     init. factor matrix with the nonzeros of PAP
     -------------------------------------------- */
  for (K = firstPostorder(PTP); K != -1; K = nextPostorder(PTP, K))
   { firstcol = nzfsub[xnzf[K]];
     lastcol = firstcol + ncolfactor[K];
     for (k = firstcol; k < lastcol; k++)
      { istart = xnza[k];
        istop = xnza[k+1];
        kstart = xnzl[k];
        hstart = xnzlsub[k];
        h = hstart;
        for (i = istart; i < istop; i++)
         { for (; nzlsub[h] != nzasub[i]; h++);
           dis = h - hstart;
           nzl[kstart+dis] = nza[i];
         }
        nzl[kstart] = diag[k];
      }
   }
}


/*****************************************************************************
******************************************************************************/
void
initFactorMtxNEW(factorMtx_t *L, inputMtx_t *PAP)
{ elimtree_t *PTP;
  frontsub_t *frontsub;
  css_t      *css;
  PORD_INT        *ncolfactor;
  FLOAT      *nzl, *nza, *diag, *entriesL;
  PORD_INT        *xnzl, *xnza, *nzasub, *xnzf, *nzfsub;
  PORD_INT        *tmp, neqs, nelem, K, k, len, row, i, istart, istop;
  PORD_INT        firstcol, lastcol;

  nelem = L->nelem;
  nzl = L->nzl;
  css = L->css;
  xnzl = css->xnzl;

  frontsub = L->frontsub;
  PTP = frontsub->PTP;
  ncolfactor = PTP->ncolfactor;
  xnzf = frontsub->xnzf;
  nzfsub = frontsub->nzfsub;

  neqs = PAP->neqs;
  diag = PAP->diag;
  nza = PAP->nza;
  xnza = PAP->xnza;
  nzasub = PAP->nzasub;

  /* ------------------------
     allocate working storage
     ------------------------ */
  mymalloc(tmp, neqs, PORD_INT);

  /* ------------------------------------
     set all numerical values of L to 0.0
     ------------------------------------ */
  for (i = 0; i < nelem; i++)
    nzl[i] = 0.0;

  /* --------------------------------------------
     init. factor matrix with the nonzeros of PAP
     -------------------------------------------- */
  for (K = firstPostorder(PTP); K != -1; K = nextPostorder(PTP, K))
   { len = 0;
     istart = xnzf[K];
     istop = xnzf[K+1];
     for (i = istart; i < istop; i++)
       tmp[nzfsub[i]] = len++;

     firstcol = nzfsub[istart];
     lastcol = firstcol + ncolfactor[K];
     entriesL = nzl + xnzl[firstcol];
     for (k = firstcol; k < lastcol; k++)
      { istart = xnza[k];
        istop = xnza[k+1];
        for (i = istart; i < istop; i++)
         { row = nzasub[i];
           entriesL[tmp[row]] = nza[i];
         }
        entriesL[tmp[k]] = diag[k];
        entriesL += --len;
      }
   }

  /* --------------------
     free working storage
     -------------------- */
  free(tmp);
}

