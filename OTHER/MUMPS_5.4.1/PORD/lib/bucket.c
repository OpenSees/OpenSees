/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: bucket.c
/
/ author        J"urgen Schulze, University of Paderborn
/ created       12/06/00
/
/ This file contains functions dealing with buckets.
/
******************************************************************************

Data type:  struct bucket
              int maxbin;   maximal bin in bucket
              int maxitem;  maximal item that can be stored in bucket
              int offset;   to store items with negative key-value
              int nobj;     number of items in bucket
              int minbin;   leftmost non-empty bin
              int *bin;     there are maxbin+1 bins (bin[0]...bin[maxbin])
              int *next;    next[item] points to next item in bin
              int *last;    last[item] points to previous item in bin
              int *key;     holds key of item (MAX_INT if item not in bucket)
Comments:
  o Any implementation of a bucket should enable insert/remove operations in
    constant time
  o There a two special bins:
    bin[0] contains all items u with key[u] + offset < 0
    bin[maxbin] contains all items u with key[u] + offset > maxbin
Methods in lib/bucket.c:
- bucket = newBucket(int maxbin, int maxitem, int offset);
    o Initial: nobj = 0 and minbin = MAX_INT
- void freeBucket(bucket_t *bucket);
- bucket = setupBucket(int maxbin, int maxitem, int offset);
    o allocates memory for the bucket by calling newBucket and initializes
      the vectors, i.e. bin[i] = -1 for all 0 <= i <= maxbin,
      next[u] = last[u] = -1, and key[u] = MAX_INT for all 0 <= u <= maxitem
- int minBucket(bucket_t *bucket);
    o returns the item whose key-value is minimal; this item is stored in
      bin[minbin]; if minbin = 0 or minbin = maxbin, the whole bin must be
      searched, since the items stored herein may have different keys
    o if nobj = 0, the function returns -1
- void insertBucket(bucket_t *bucket, int k, int item);
    o insert item with key k in bucket; if key[item] != MAX_INT (i.e. item
      already in bucket) or if item > maxitem the program terminates
- void removeBucket(bucket_t *bucket, int item);
    o removes item from bucket; if key[item] == MAX_INT (i.e. item not in
      bucket) the program terminates

******************************************************************************/

#include <space.h>


/******************************************************************************
******************************************************************************/
bucket_t*
newBucket(PORD_INT maxbin, PORD_INT maxitem, PORD_INT offset)
{ bucket_t *bucket;

  mymalloc(bucket, 1, bucket_t);
  mymalloc(bucket->bin, (maxbin+1), PORD_INT);
  mymalloc(bucket->next, (maxitem+1), PORD_INT);
  mymalloc(bucket->last, (maxitem+1), PORD_INT);
  mymalloc(bucket->key, (maxitem+1), PORD_INT);

  bucket->maxbin = maxbin;
  bucket->maxitem = maxitem;
  bucket->offset = offset;
  bucket->nobj = 0;
  bucket->minbin = MAX_INT;

  return(bucket);
}


/******************************************************************************
******************************************************************************/
void
freeBucket(bucket_t *bucket)
{
  free(bucket->bin);
  free(bucket->next);
  free(bucket->last);
  free(bucket->key);
  free(bucket);
}


/******************************************************************************
******************************************************************************/
bucket_t*
setupBucket(PORD_INT maxbin, PORD_INT maxitem, PORD_INT offset)
{ bucket_t *bucket;
  PORD_INT      i, u;

  if (offset < 0)
   { fprintf(stderr, "\nError in function setupBucket\n"
          "  offset must be >= 0\n");
     quit();
   }

  bucket = newBucket(maxbin, maxitem, offset);

  for (i = 0; i <= maxbin; i++)
    bucket->bin[i] = -1;
  for (u = 0; u <= maxitem; u++)
    { bucket->next[u] = bucket->last[u] = -1;
      bucket->key[u] = MAX_INT;
    }

  return(bucket);
}


/******************************************************************************
******************************************************************************/
PORD_INT
minBucket(bucket_t *bucket)
{ PORD_INT *bin, *next, *key, maxbin, minbin, nobj;
  PORD_INT item, bestitem, bestkey;

  maxbin = bucket->maxbin;
  nobj = bucket->nobj;
  minbin = bucket->minbin;
  bin = bucket->bin;
  next = bucket->next;
  key = bucket->key;

  if (nobj > 0)
   { /* ---------------------------------------------
        get the first item from leftmost nonempty bin
        --------------------------------------------- */
     while (bin[minbin] == -1) minbin++;
     bucket->minbin = minbin;
     bestitem = bin[minbin];
     bestkey = minbin;

     /* --------------------------------------------------
        items in bins 0 and maxbin can have different keys
        => search for item with smallest key
        -------------------------------------------------- */
     if ((minbin == 0) || (minbin == maxbin))
      { item = next[bestitem];
        while (item != -1)
         { if (key[item] < bestkey)
            { bestitem = item;
              bestkey = key[item];
            }
           item = next[item];
         }
      }
     /* ---------------------------------
        return the item with smallest key
        --------------------------------- */
     return(bestitem);
   }
  else return(-1);
}


/******************************************************************************
******************************************************************************/
void
insertBucket(bucket_t *bucket, PORD_INT k, PORD_INT item)
{ PORD_INT s, nextitem;

  /* ------------------------------------
     check whether there are any problems
     ------------------------------------ */
  if (abs(k) >= MAX_INT - bucket->offset - 1)
   { fprintf(stderr, "\nError in function insertBucket\n"
          "  key %d too large/small for bucket\n", k);
     quit();
   }
  if (item > bucket->maxitem)
   { fprintf(stderr, "\nError in function insertBucket\n"
          "  item %d too large for bucket (maxitem is %d)\n", item,
          bucket->maxitem);
     quit();
   }
  if (bucket->key[item] != MAX_INT)
   { fprintf(stderr, "\nError in function insertBucket\n"
          "  item %d already in bucket\n", item);
     quit();
   }

  /* -------------------------------------
     determine the bin that holds the item
     ------------------------------------- */
  s = max(0, (k + bucket->offset));
  s = min(s, bucket->maxbin);

  /* --------------------------------------------------------------
     adjust minbin, increase nobj, and mark item as being in bucket
     -------------------------------------------------------------- */
  bucket->minbin = min(bucket->minbin, s);
  bucket->nobj++;
  bucket->key[item] = k;

  /* -----------------------------
     finally, insert item in bin s
     ----------------------------- */
  nextitem = bucket->bin[s];
  if (nextitem != -1)
    bucket->last[nextitem] = item;
  bucket->next[item] = nextitem;
  bucket->last[item] = -1;
  bucket->bin[s] = item;
}


/******************************************************************************
******************************************************************************/
void
removeBucket(bucket_t *bucket, PORD_INT item)
{ PORD_INT s, nextitem, lastitem;

  /* ----------------------------
     check whether item in bucket
     ---------------------------- */
  if (bucket->key[item] == MAX_INT)
   { fprintf(stderr, "\nError in function removeBucket\n"
          "  item %d is not in bucket\n", item);
     quit();
   }

  /* -----------------------
     remove item from bucket
     ----------------------- */
  nextitem = bucket->next[item];
  lastitem = bucket->last[item];
  if (nextitem != -1)
    bucket->last[nextitem] = lastitem;
  if (lastitem != -1)
    bucket->next[lastitem] = nextitem;
  else
   { s = max(0, (bucket->key[item] + bucket->offset));
     s = min(s, bucket->maxbin);
     bucket->bin[s] = nextitem;
   }

  /* --------------------------------------------
     decrease nobj and mark item as being removed
     -------------------------------------------- */
  bucket->nobj--;
  bucket->key[item] = MAX_INT;
}
