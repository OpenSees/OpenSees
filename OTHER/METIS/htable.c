/*
 * Copyright 1995, Regents of the University of Minnesota
 *
 * htable.c
 *
 * This file implements a simple hash table scheme.
 *
 * Started 4/11/95
 * George
 *
 * $Id: htable.c,v 1.1.1.1 2000-09-15 08:23:12 fmk Exp $
 *
 */

#include "multilevel.h"


static PRIMES[30] = {
  51, 	101, 	211, 	307, 	449, 	587, 	797, 	991, 	1301, 	1723, 	2203, 
  2543, 3001, 	4001, 	5003, 	6007, 	7001, 	8009, 	9719, 	12007, 	15277,
  17321, 19181, 20599, 23173, 26227, 29327, 35279, 40289, 50311};


/*************************************************************************
* This function initializes the htable
**************************************************************************/
void CreateHTable(HTableType *htable, int size)
{
  htable->ht = ismalloc(size, HT_EMPTY, "CreateHTable: htable->ht");
  htable->size = size;
  htable->nelem = 0;
}

/*************************************************************************
* This function is given 3 numbers and picks a prime that is higher than
* the max of this numbers
**************************************************************************/
int SelHTSize(int max)
{
  int i;

  for (i=0; i<30; i++)
    if (max <= PRIMES[i])
      return PRIMES[i];

  return max+11;
}


/*************************************************************************
* This function inserts an element into the htable
**************************************************************************/
void AddHTable(HTableType *htable, int val)
{
  int i, k, *ht;

  if (htable->nelem == htable->size)
    IncreaseHTable(htable);

  k = HTVALUE(val, htable->size);
  htable->nelem++;
  ht = htable->ht;

  if (ht[k] == HT_EMPTY) {
    ht[k] = val;
    return;
  }

  for (i=k+1; i<htable->size; i++) {
    if (ht[i] == HT_EMPTY) {
      ht[i] = val;
      return;
    }
  }
  for (i=0; i<k; i++) {
    if (ht[i] == HT_EMPTY) {
      ht[i] = val;
      return;
    }
  }

  errexit("Hashtable is full! %d", htable->size);
}


/*************************************************************************
* This function deletes an element from the htable
**************************************************************************/
void DelHTable(HTableType *htable, int val)
{
  int i, k, *ht;

  k = HTVALUE(val, htable->size);
  htable->nelem--;
  ht = htable->ht;

  if (ht[k] == val) {
    ht[k] = HT_EMPTY;
    return;
  }

  for (i=k+1; i<htable->size; i++) {
    if (ht[i] == val) {
      ht[i] = HT_EMPTY;
      return;
    }
  }
  for (i=0; i<k; i++) {
    if (ht[i] == val) {
      ht[i] = HT_EMPTY;
      return;
    }
  }

  errexit("%d not found in Hashtable! %d", val);
}


/*************************************************************************
* This function inserts an element into the htable
**************************************************************************/
void IncreaseHTable(HTableType *htable)
{
  int ii, i, k, val, *oldht, *newht, newsize;

  newsize = SelHTSize(1.5*htable->size);
  oldht = htable->ht;
  newht = ismalloc(newsize, HT_EMPTY, "IncreaseHTable: newht");

/*
  printf("Rehashing %d to %d\n", htable->size, newsize);
*/

  for (ii=0; ii<htable->size; ii++) {
    val = oldht[ii];
    k = HTVALUE(val, newsize);

    if (newht[k] == HT_EMPTY) {
      newht[k] = val;
      continue;
    }

    for (i=k+1; i<newsize; i++) {
      if (newht[i] == HT_EMPTY) {
        newht[i] = val;
        break;
      }
    }
    if (i == newsize) {
      for (i=0; i<k; i++) {
        if (newht[i] == HT_EMPTY) {
          newht[i] = val;
          break;
        }
      }
      if (i==k)
        errexit("NewHashtable is full! %d %d %d", newsize, val, k);
    }
  }

  free(oldht);
  htable->size = newsize;
  htable->ht = newht;

}
