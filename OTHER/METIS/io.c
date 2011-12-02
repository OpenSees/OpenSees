/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * io.c
 *
 * This file contains routines related to I/O
 *
 * Started 8/28/94
 * George
 *
 * $Id: io.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include "multilevel.h"

extern int IsWeighted;	/* main.c */
extern timer IOTmr;	/* main.c */

/*************************************************************************
* This function reads the spd matrix
**************************************************************************/
void readgraph(CoarseGraphType *graph, char *filename)
{
  char line[MAXLINE+1], *oldstr, *newstr;
  int i, j, k, fmt;
  FILE *fpin;
  EdgeType *edges, aux;
  int lastedge;
  int readew, readvw; 
  int error;


  starttimer(&IOTmr);

  error = 0;
  InitGraph(graph);

  if ((fpin = fopen(filename, "r")) == NULL) 
    errexit("Failed to open file %s", filename);

  do {
    fgets(line, MAXLINE, fpin);
  } while (line[0] == '%');

  if (sscanf(line,"%d %d %d",&(graph->nvtxs), &(graph->nedges), &fmt) == 2)
    fmt = 0;

  readew = (fmt%10 > 0);
  readvw = ((fmt/10)%10 > 0);
  if (fmt >= 100)
    errexit("Cannot read this type of file format!");

  IsWeighted = readvw;

  graph->nedges *=2;

  if (graph->nvtxs > MAXIDX) 
    errexit("\nThe matrix is too big: %d", graph->nvtxs);

  graph->level = 0;
  graph->tvwgt = 0;
  graph->vtxs = (VertexType **)GKmalloc(sizeof(VertexType *)*graph->nvtxs, "readgraph: graph->vtxs");
  graph->label = imalloc(graph->nvtxs, "readgraph: graph->label");
  edges = (EdgeType *)GKmalloc(graph->nedges*sizeof(EdgeType), "readgraph: edges");
  lastedge = 0;
  graph->allvtxs = (VertexType *)GKmalloc(sizeof(VertexType)*graph->nvtxs, "readgraph: allvtxs");

  /* Start reading the graph file */
  for (i=0; i<graph->nvtxs; i++) {
    do {
      fgets(line, MAXLINE, fpin);
    } while (line[0] == '%');
    oldstr = line;
    newstr = NULL;

    if (strlen(line) == MAXLINE) 
      errexit("\nBuffer for fgets not big enough!\n");

    graph->label[i] = i;
    graph->vtxs[i] = graph->allvtxs+i;
    graph->vtxs[i]->edges = edges + lastedge;
    graph->vtxs[i]->nedges = lastedge;
    graph->vtxs[i]->cewgt = 0;
    graph->vtxs[i]->ewgtsum = 0;

    if (readvw) {
      graph->vtxs[i]->vwgt = (int)strtol(oldstr, &newstr, 10);
      oldstr = newstr;
    }
    else
      graph->vtxs[i]->vwgt = 1;

    graph->tvwgt += graph->vtxs[i]->vwgt;
    for (;;) {
      aux.edge = (int)strtol(oldstr, &newstr, 10) - 1;
      oldstr = newstr;

      if (readew) {
        aux.ewgt = (int)strtol(oldstr, &newstr, 10);
        oldstr = newstr;
      }
      else
        aux.ewgt = 1;

      if (aux.edge < 0)
        break;

      edges[lastedge++] = aux;
      graph->vtxs[i]->ewgtsum += aux.ewgt;

      if (aux.edge < i) {
        k = aux.edge;
        for (j=0; j<graph->vtxs[k]->nedges; j++) {
          if (graph->vtxs[k]->edges[j].edge == i) {
            if (graph->vtxs[k]->edges[j].ewgt != aux.ewgt) {
              printf("Edge (%d %d) and (%d %d) do not have the same weight (%d %d)!\n", k+1, i+1, i+1, k+1, graph->vtxs[k]->edges[j].ewgt, aux.ewgt);
              error = 1;
            }
            break;
          }
        }
        if (j == graph->vtxs[k]->nedges) {
          printf("Edge (%d %d) does not exist (even-though edge (%d %d) exists)!\n", k+1, i+1, i+1, k+1);
          error = 1;
        }
      }
    } 
    graph->vtxs[i]->nedges = lastedge-graph->vtxs[i]->nedges;
  }

  fclose(fpin);

  if (error)
    errexit("Input graph file contains errors. Cannot continue...");

  if (lastedge != graph->nedges)
    errexit("readgraph: Something wrong with the edges from input file %d %d",graph->nedges, lastedge);

  stoptimer(&IOTmr);
}


/*************************************************************************
* This function reads the spd matrix
**************************************************************************/
void ReadGraphSKIT(int *xadj, int *adjncy, char *filename)
{
  char line[MAXLINE+1], *oldstr, *newstr;
  int iedges[2*MAXDEGREE];
  int i, j, k;
  int nz;
  FILE *fpin;
  int n, m;			/* The number of vertices and edges */
  int node;


  starttimer(&IOTmr);
  fpin = GKfopen(filename, "r", "ReadGraphSKIT: failed to open file");

  fgets(line, MAXLINE, fpin);

  sscanf(line,"%d %d",&n, &m);

  k = 0;
  /* Start reading the graph file */
  for (i=0; i<n; i++) {
    fgets(line, MAXLINE, fpin);
    if (strlen(line) == MAXLINE) 
      errexit("\nBuffer for fgets not big enough!\n");

    nz = 0;
    oldstr = line;
    iedges[nz] = (int)strtol(oldstr, &newstr, 10) - 1;

    while (oldstr != newstr) {
      nz++; 
      oldstr = newstr;
      iedges[nz] = (int)strtol(oldstr, &newstr, 10) - 1;
      ASSERT(nz < 2*MAXDEGREE-1);
    }

    xadj[i] = k+1;
    for (j=0; j<nz; j++) {
      node = iedges[j];
      if (i != node)
        adjncy[k++] = node+1;
    }
  }
  xadj[n] = k+1;

  GKfclose(fpin);

  stoptimer(&IOTmr);
}



/*************************************************************************
* This function prints a graph
**************************************************************************/
void WriteGraph(CoarseGraphType *graph, char *filename)
{
  char fileout[256];
  FILE *fpout;
  int i, j;

  sprintf(fileout, "%s.%d", filename, graph->nvtxs);
  fpout = GKfopen(fileout, "w", "WriteGraph: failed to open file");

  fprintf(fpout, "%d %d 11\n",graph->nvtxs, graph->nedges/2);
  for (i=0; i<graph->nvtxs; i++) {
    fprintf(fpout, "%d ", graph->vtxs[i]->vwgt);
    for (j=0; j<graph->vtxs[i]->nedges; j++)
      fprintf(fpout,"%d %d ",graph->vtxs[i]->edges[j].edge+1, graph->vtxs[i]->edges[j].ewgt);
    fprintf(fpout, "\n");
  }

  GKfclose(fpout);
}

/*************************************************************************
* This function prints a graph
**************************************************************************/
void printgraph(CoarseGraphType *graph, FILE *fpout)
{
  int i, j;

  fprintf(fpout, "\n%d",graph->nvtxs);
  for (i=0; i<graph->nvtxs; i++) {
    fprintf(fpout, "\n%6d %3d %3d\t", i, graph->vtxs[i]->vwgt, graph->vtxs[i]->nedges);
    for (j=0; j<graph->vtxs[i]->nedges; j++)
      fprintf(fpout,"[%d, %d] ",graph->vtxs[i]->edges[j].edge, graph->vtxs[i]->edges[j].ewgt);
  }

  if (graph->cmap != NULL) {
    fprintf(fpout,"\n\nCoersion Map");
    for (i=0; i<graph->nvtxs; i++) {
      if (i%10 == 0)
        fprintf(fpout,"\n");
      fprintf(fpout,"%6d ",graph->cmap[i]);
    }
    fprintf(fpout,"\n\nMatchings");

    for (i=0; i<graph->nvtxs; i++) {
      if (i%5 == 0)
        fprintf(fpout,"\n");
      fprintf(fpout,"[%6d, %6d] ",i, graph->match[i]);
    }
  }

  fprintf(fpout,"\n");

}



/*************************************************************************
* This function writes out the permutation vector
**************************************************************************/
void WriteOrder(char *fname, int *perm, int n)
{
  FILE *fpout;
  int i;
  char filename[256];

  starttimer(&IOTmr);

  sprintf(filename, "%s.perm", fname);

  fpout = GKfopen(filename, "w", "Problems in opening the permutation file");

  for (i=0; i<n; i++)
    fprintf(fpout,"%d\n",perm[i]);

  GKfclose(fpout);

  stoptimer(&IOTmr);
}



/*************************************************************************
* This function writes out the partition vector
**************************************************************************/
void WritePartition(char *fname, int *part, int n, int nparts)
{
  FILE *fpout;
  int i;
  char filename[256];

  starttimer(&IOTmr);

  sprintf(filename,"%s.part.%d",fname, nparts);

  fpout = GKfopen(filename, "w", "Problems in opening the partition file");

  for (i=0; i<n; i++)
    fprintf(fpout,"%d\n",part[i]);

  fclose(fpout);

  stoptimer(&IOTmr);
}




/*************************************************************************
* This function writes out the graph in .grid format
**************************************************************************/
void WriteCoarseGraph2Grid(char *filename, CoarseGraphType *graph)
{
  int i, k;
  int *rowidx;
  int *colidx;
  char file[80];
  FILE *afp;
  int maxdeg=0;


  colidx = imalloc(graph->nvtxs+1, "WriteCoarseGraph2Matlab: colidx");
  rowidx = imalloc(graph->nedges+graph->nvtxs, "WriteCoarseGraph2Matlab: rowidx");

  colidx[0] = 0;
  for (i=0; i<graph->nvtxs; i++) {
    memcpy(rowidx+colidx[i], graph->vtxs[i]->edges, graph->vtxs[i]->nedges*sizeof(int));
    colidx[i+1] = colidx[i] + graph->vtxs[i]->nedges;

    if (maxdeg < graph->vtxs[i]->nedges)
      maxdeg = graph->vtxs[i]->nedges;
  }


  sprintf(file, "%s.grid",filename);
  afp = fopen(file,"w");
  fprintf(afp, "%10d%10d%10d%10d%10d%10d\n",graph->nvtxs, graph->nedges, 0, maxdeg, 2, 0);

  k = 0;
  for (i=0; i<graph->nvtxs+1; i++) {
    fprintf(afp, "%10d", colidx[i]);
    k++;
    if (k==8) {
      fprintf(afp, "\n");
      k = 0;
    }
  }
  if (k != 0)
    fprintf(afp, "\n");

  k = 0;
  for (i=0; i<graph->nedges; i++) {
    fprintf(afp, "%10d", rowidx[i]);
    k++;
    if (k==8) {
      fprintf(afp, "\n");
      k = 0;
    }
  }
  if (k != 0)
    fprintf(afp, "\n");

  k = 0;
  for (i=0; i<graph->nvtxs; i++) {
    fprintf(afp, "%10d", -1);
    k++;
    if (k==8) {
      fprintf(afp, "\n");
      k = 0;
    }
  }
  if (k != 0)
    fprintf(afp, "\n");

  fclose(afp);


  GKfree(rowidx, colidx, -1);
}


