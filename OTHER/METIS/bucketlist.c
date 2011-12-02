/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * bucketlist.c
 *
 * This file contains functions for manipulating the bucket list
 * representation of the gains associated with each vertex in a graph.
 * These functions are used by the refinement algorithms
 *
 * Started 9/2/94
 * George
 *
 * $Id: bucketlist.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include "multilevel.h"


/*************************************************************************
* External global variables 
**************************************************************************/
extern CtrlType *__Ctrl;		/* mlevelpart.c */



/*************************************************************************
* This function initializes the klpartdef data structures
**************************************************************************/
void initbucket(BucketListType *part, int orgnvtxs, int nvtxs, int maxnodes, int level)
{
  part->nnodes = 0;
  part->maxnodes = maxnodes;
  part->pgainspan = amin(PLUS_GAINSPAN, orgnvtxs/4);
  part->ngainspan = amin(NEG_GAINSPAN, orgnvtxs/8);

  if (orgnvtxs/nvtxs < 10 && level != -1) { 
    part->type = 1;
    if (orgnvtxs == nvtxs) {
      part->pgainspan = MAXDEGREE;
      part->ngainspan = MAXDEGREE;
    }
  }
  else {
    part->type = 2;
  }

  __Ctrl->cgain = 0;	/* Reset the gaincore */

  switch (part->type) {
    case 1:
      part->buckets = (ListNodeType **)
        GKmalloc(sizeof(ListNodeType *)*(part->ngainspan + part->pgainspan + 1), "initbucket: part->buckets");
      memset((void *)part->buckets, (int)NULL, sizeof(ListNodeType *)*(part->ngainspan+part->pgainspan+1));

      part->buckets += part->ngainspan;  /* Advance buckets by the ngainspan proper indexing */

      part->mingain = part->pgainspan;
      part->maxgain = -part->ngainspan;
      break;
    case 2:
      __Ctrl->cbucket = 0;
      part->head = part->tail = NULL;
      break;
  }

}


/*************************************************************************
* This function resets the buckets
**************************************************************************/
void resetbucket(BucketListType *part)
{
  part->nnodes = 0;

  __Ctrl->cgain = 0;	/* Reset the gaincore */

  switch (part->type) {
    case 1:
      part->mingain = part->pgainspan;
      part->maxgain = -part->ngainspan;
      memset((void *)(part->buckets-part->ngainspan), (int)NULL, sizeof(ListNodeType *)*(part->ngainspan + part->pgainspan+1));
      break;
    case 2:
      part->head = part->tail = NULL;
      __Ctrl->cbucket = 0;
      break;
  }

}


/*************************************************************************
* This function frees the buckets
**************************************************************************/
void freebucket(BucketListType *part)
{
  switch (part->type) {
    case 1:
      part->buckets -= part->ngainspan;
      free(part->buckets);
      break;
    case 2:
      break;
  }
}


/*************************************************************************
* This function adds a node of certain gain into a partition
**************************************************************************/
int Add2Part(BucketListType *part, int node, int gain)
{
  ListNodeType *newgain;
  GainBucketType *cptr, *pptr, *nptr;

  switch (part->type) {
    case 1:
      /* printf("Added: %d %d, %d %d\n", node, gain, part->mingain, part->maxgain); */
      if (gain < -part->ngainspan) /* If out of bounds do nothing */
        return 0;	

      if (gain > part->pgainspan) 
        gain = part->pgainspan; /* Set the gain to arbitrarily be at the highest bucket */

      part->nnodes++;
      newgain = __Ctrl->gaincore + __Ctrl->cgain++;
      newgain->pos = node;

      part->buckets[gain] = addthisnode(part->buckets[gain], newgain);
      if (part->maxgain < gain)
        part->maxgain = gain;
      if (part->mingain > gain)
        part->mingain = gain;

      /* printf("Added: %d %d, %d %d\n", node, gain, part->mingain, part->maxgain); */
      break;
    case 2:
      part->nnodes++;
      newgain = __Ctrl->gaincore + __Ctrl->cgain++;
      newgain->pos = node;

      for (cptr = pptr = part->tail; cptr != NULL; pptr = cptr, cptr = cptr->nbucket) 
        if (cptr->gain >= gain)
          break;

      nptr = __Ctrl->gbcore + __Ctrl->cbucket++;
      if (__Ctrl->cbucket == __Ctrl->maxbucket) {
        printf("\n***Warning*** Premature exit from refinement.\n");
        return -1;
      }

      if (pptr == NULL) {  /* First node */
        nptr->gain = gain;
        nptr->gainlist = NULL;
        nptr->gainlist = addthisnode(nptr->gainlist, newgain);
        nptr->nbucket = nptr->pbucket = NULL;
        part->head = part->tail = nptr;
      }
      else if (cptr == NULL) {  /* I got a new max gain */
        nptr->gain = gain;
        nptr->gainlist = NULL;
        nptr->gainlist = addthisnode(nptr->gainlist, newgain);
        nptr->nbucket = NULL;
        nptr->pbucket = pptr;
        pptr->nbucket = nptr;
        part->head = nptr;
      }
      else if (cptr->gain == gain) {  /* Bingo! */
        cptr->gainlist = addthisnode(cptr->gainlist, newgain);

        __Ctrl->cbucket--;
      }
      else { /* Put it somewhere in the middle */
        nptr->gain = gain;
        nptr->gainlist = NULL;
        nptr->gainlist = addthisnode(nptr->gainlist, newgain);

        nptr->nbucket = cptr;
        nptr->pbucket = cptr->pbucket;
        cptr->pbucket = nptr;

        if (pptr == cptr) /* Inserting it in the beginning */
          part->tail = nptr;
        else
          pptr->nbucket = nptr;
      }
      break;
  }

  return 0;
}


/*************************************************************************
* This function deletes a node from a partition and reinserts it with
* an updated gain
**************************************************************************/
int UpdatePart(BucketListType *part, int node, int oldgain, int newgain)
{
  ListNodeType *gainptr;
  GainBucketType *cptr, *pptr, *nptr;

  switch (part->type) {
    case 1:
      if (oldgain < -part->ngainspan) /* Is it realy in? */
        return Add2Part(part, node, newgain);

      if (oldgain > part->pgainspan) /* Is the the very high bucket ? */
        oldgain = part->pgainspan;
      if (newgain > part->pgainspan) /* Is the the very high bucket ? */
        newgain = part->pgainspan;

      part->nnodes--;
      part->buckets[oldgain] = delthisnode(part->buckets[oldgain], node, &gainptr);

      if (part->buckets[oldgain] == NULL) {
        if (oldgain == part->maxgain) {
          if (part->nnodes == 0) {
            part->mingain = part->pgainspan;
            part->maxgain = -part->ngainspan;
          }
          else 
            for (; part->buckets[part->maxgain]==NULL; part->maxgain--);
        }
        else if (oldgain == part->mingain) {
          if (part->nnodes == 0) {
            part->mingain = part->pgainspan;
            part->maxgain = -part->ngainspan;
          }
          else 
            for (; part->buckets[part->mingain]==NULL; part->mingain++);
        }
      }

      if (newgain < -part->ngainspan) /* Shall I put it ? */
        return 0;

      part->nnodes++;
      part->buckets[newgain] = addthisnode(part->buckets[newgain], gainptr);
      if (part->maxgain < newgain)
        part->maxgain = newgain;
      if (part->mingain > newgain)
        part->mingain = newgain;
      break;
    case 2:
      for (cptr = part->tail; cptr != NULL; cptr = cptr->nbucket) 
        if (cptr->gain >= oldgain)
          break;

      if (cptr == NULL || cptr->gain != oldgain || cptr->gainlist == NULL) {
        errexit("UpdatePart failed. GainBucket list is corrupted %d %d %d", node, oldgain, cptr->gain);
      }

      cptr->gainlist = delthisnode(cptr->gainlist, node, &gainptr);

      if (cptr->gainlist == NULL) { 
        if (cptr->nbucket == NULL) {
          part->head = cptr->pbucket;
          if (cptr->pbucket == NULL)
            part->tail = NULL;
          else
            cptr->pbucket->nbucket = NULL;
        }
        else if (cptr->pbucket == NULL) {
          part->tail = cptr->nbucket;
          cptr->nbucket->pbucket = NULL;
        }
        else {
          cptr->nbucket->pbucket = cptr->pbucket;
          cptr->pbucket->nbucket = cptr->nbucket;
        }
      }

      for (cptr = pptr = part->tail; cptr != NULL; pptr = cptr, cptr = cptr->nbucket) 
        if (cptr->gain >= newgain)
          break;

      nptr = __Ctrl->gbcore + __Ctrl->cbucket++;
      if (__Ctrl->cbucket == __Ctrl->maxbucket) {
        printf("\n***Warning*** Premature exit from refinement.\n");
        return -1;
      }

      if (pptr == NULL) {  /* First node */
        nptr->gain = newgain;
        nptr->gainlist = NULL;
        nptr->gainlist = addthisnode(nptr->gainlist, gainptr);
        nptr->nbucket = nptr->pbucket = NULL;
        part->head = part->tail = nptr;
      }
      else if (cptr == NULL) {  /* I got a new max gain */
        nptr->gain = newgain;
        nptr->gainlist = NULL;
        nptr->gainlist = addthisnode(nptr->gainlist, gainptr);
        nptr->nbucket = NULL;
        nptr->pbucket = pptr;
        pptr->nbucket = nptr;
        part->head = nptr;
      }
      else if (cptr->gain == newgain) {  /* Bingo! */
        cptr->gainlist = addthisnode(cptr->gainlist, gainptr);

        __Ctrl->cbucket--;
      }
      else { /* Put it somewhere in the middle */
        nptr->gain = newgain;
        nptr->gainlist = NULL;
        nptr->gainlist = addthisnode(nptr->gainlist, gainptr);

        nptr->nbucket = cptr;
        nptr->pbucket = cptr->pbucket;
        cptr->pbucket = nptr;

        if (pptr == cptr) /* Inserting it in the beginning */
          part->tail = nptr;
        else
          pptr->nbucket = nptr;
      }
      break;
  }

  return 0;
}



/*************************************************************************
* This function returns the vertex with the largest gain from a partition
* and removes the node from the bucket list
**************************************************************************/
int GetMaxGainVtx(BucketListType *part)
{
  int vtx;
  ListNodeType *tptr;
  GainBucketType *cptr;

  if (part->nnodes == 0)
    return -1;

  part->nnodes--;

  switch (part->type) {
    case 1:
      tptr = part->buckets[part->maxgain];
      part->buckets[part->maxgain] = tptr->link;;

      if (tptr->link == NULL) {
        if (part->nnodes == 0) {
          part->mingain = part->pgainspan;
          part->maxgain = -part->ngainspan;
        }
        else 
          for (; part->buckets[part->maxgain]==NULL; part->maxgain--);
      }
      vtx = tptr->pos;
      break;
    case 2:
      cptr = part->head;
      vtx = cptr->gainlist->pos;
      cptr->gainlist = cptr->gainlist->link;

      if (cptr->gainlist == NULL) {
        if (part->nnodes == 0) {
          part->tail = part->head = NULL;
        }
        else {
          part->head = cptr->pbucket;
          cptr->pbucket->nbucket = NULL;
        }
      }
      break;
  }

  return vtx;
}
      

/*************************************************************************
* This function returns the vertex with the largest gain from a partition
**************************************************************************/
int SeeMaxGainVtx(BucketListType *part)
{
  int vtx;

  if (part->nnodes == 0)
    return -1;

  switch (part->type) {
    case 1:
      vtx = part->buckets[part->maxgain]->pos;
      break;
    case 2:
      vtx = part->head->gainlist->pos;
      break;
  }

  return vtx;
}
      



/*************************************************************************
* This functions prints the nodes with positive gains 
**************************************************************************/
void PrintPlusPart(char *title, BucketListType *part)
{
  int i;
  ListNodeType *ptr;
  GainBucketType *cptr;


  switch (part->type) {
    case 1:
      printf("%s min:%d \tmax:%d\n",title, part->mingain, part->maxgain);
      for (i=1; i<=part->maxgain; i++) {
        if (part->buckets[i] != NULL) {
          printf("%d \t",i);
          for (ptr = part->buckets[i]; ptr != NULL; ptr = ptr->link)
            printf("%d ",ptr->pos);
          printf("\n");
        }
      }
      break;
    case 2:
      printf("\n%d\n", part->nnodes);
      for (cptr = part->tail; cptr != NULL; cptr = cptr->nbucket) {
        printf("%d \t",cptr->gain);
        for (ptr = cptr->gainlist; ptr != NULL; ptr = ptr->link)
          printf("%d ",ptr->pos);
        printf("\n");
      }
      printf("----\n");
      for (cptr = part->head; cptr != NULL; cptr = cptr->pbucket) {
        printf("%d \t",cptr->gain);
        for (ptr = cptr->gainlist; ptr != NULL; ptr = ptr->link)
          printf("%d ",ptr->pos);
        printf("\n");
      }
      break;
  }

}


/*************************************************************************
* This functions prints the nodes with positive gains 
**************************************************************************/
void PrintPartGains(BucketListType *part1, BucketListType *part2)
{
  printf("Partitions: [%5d, %5d] [%5d, %5d]\n",
    part1->mingain, part1->maxgain, part2->mingain, part2->maxgain);
}
