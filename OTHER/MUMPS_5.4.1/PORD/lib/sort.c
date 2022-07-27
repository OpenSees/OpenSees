/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: sort.c
/
/ author        J"urgen Schulze, University of Paderborn
/ created       09/15/99
/
/ This file contains some sorting functions. the code is adopted from
/ the book "Algorithms in C" by R. Sedgewick.
/
******************************************************************************/

#include <space.h>
#define THRES 10


/*****************************************************************************
/ insertion sort upwards (INTS, without keys)
******************************************************************************/
void
insertUpInts(PORD_INT n, PORD_INT *array)
{ PORD_INT i, j, v;

  for (i = 1; i < n; i++)
   { v = array[i]; j = i;
     while ((j > 0) && (array[j-1] > v))
      { array[j] = array[j-1];
        j--;
      }
     array[j] = v;
   }
}


/*****************************************************************************
/ insertion sort upwards (INTS, with static INT keys)
******************************************************************************/
void
insertUpIntsWithStaticIntKeys(PORD_INT n, PORD_INT *array, PORD_INT *key)
{ PORD_INT   i, j, ke;
  PORD_INT   e;

  for (i = 1; i < n; i++)
   { e = array[i]; ke = key[e]; j = i;
     while ((j > 0) && (key[array[j-1]] > ke))
      { array[j] = array[j-1];
        j--;
      }
     array[j] = e;
   }
}


/*****************************************************************************
/ insertion sort downwards (INTS, with static INT keys)
******************************************************************************/
void
insertDownIntsWithStaticFloatKeys(PORD_INT n, PORD_INT *array, FLOAT *key)
{ PORD_INT   i, j, e;
  FLOAT ke;

  for (i = 1; i < n; i++)
   { e = array[i]; ke = key[e]; j = i;
     while ((j > 0) && (key[array[j-1]] < ke))
      { array[j] = array[j-1];
        j--;
      }
     array[j] = e;
   }
}


/*****************************************************************************
/ insertion sort upwards (FLOATS, with INT keys)
******************************************************************************/
void
insertUpFloatsWithIntKeys(PORD_INT n, FLOAT *array, PORD_INT *key)
{ PORD_INT   i, j, ke;
  FLOAT e;

  for (i = 1; i < n; i++)
   { e = array[i]; ke = key[i]; j = i;
     while ((j > 0) && (key[j-1] > ke))
      { array[j] = array[j-1];
        key[j] = key[j-1];
        j--;
      }
     array[j] = e;
     key[j] = ke;
   }
}


/*****************************************************************************
/ median-of-three quicksort upwards (INTS, without keys)
******************************************************************************/
void
qsortUpInts(PORD_INT n, PORD_INT *array, PORD_INT *stack)
{ register PORD_INT i, j;
  PORD_INT t, l, m, r, p;

  l = 0; r = n-1; p = 2;
  while (p > 0)
   if ((r-l) > THRES)
    { m = l + ((r-l) >> 1);
      if (array[l] > array[r]) swap(array[l], array[r], t);
      if (array[l] > array[m]) swap(array[l], array[m], t);
      if (array[r] > array[m]) swap(array[m], array[r], t);
      m = array[r]; i = l-1; j = r;
      for (;;)
       { while (array[++i] < m);
         while (array[--j] > m);
         if (i >= j) break;
         swap(array[i], array[j], t);
       }
      swap(array[i], array[r], t);
      if ((i-l) > (r-i))
       { stack[p++] = l;
         stack[p++] = i-1;
         l = i+1;
       }
      else
       { stack[p++] = i+1;
         stack[p++] = r;
         r = i-1;
       }
    }
   else
    { r = stack[--p];
      l = stack[--p];
    }
   if (THRES > 0) insertUpInts(n, array);
}


/*****************************************************************************
/ median-of-three quicksort upwards (FLOATS, with INT keys)
******************************************************************************/
void
qsortUpFloatsWithIntKeys(PORD_INT n, FLOAT *array, PORD_INT *key, PORD_INT *stack)
{ register PORD_INT i, j;
  PORD_INT   t, l, m, r, p;
  FLOAT e;

  l = 0; r = n-1; p = 2;
  while (p > 0)
   if ((r-l) > THRES)
    { m = l + ((r-l) >> 1);
      if (key[l] > key[r])
       { swap(array[l], array[r], e); swap(key[l], key[r], t); }
      if (key[l] > key[m])
       { swap(array[l], array[m], e); swap(key[l], key[m], t); }
      if (key[r] > key[m])
       { swap(array[m], array[r], e); swap(key[m], key[r], t); }
      m = key[r]; i = l-1; j = r;
      for (;;)
       { while (key[++i] < m);
         while (key[--j] > m);
         if (i >= j) break;
         swap(array[i], array[j], e); swap(key[i], key[j], t);
       }
      swap(array[i], array[r], e); swap(key[i], key[r], t);
      if ((i-l) > (r-i))
       { stack[p++] = l;
         stack[p++] = i-1;
         l = i+1;
       }
      else
       { stack[p++] = i+1;
         stack[p++] = r;
         r = i-1;
       }
    }
   else
    { r = stack[--p];
      l = stack[--p];
    }
   if (THRES > 0) insertUpFloatsWithIntKeys(n, array, key);
}


/*****************************************************************************
/ distribution counting (INTS, with static INT keys)
******************************************************************************/
void
distributionCounting(PORD_INT n, PORD_INT *node, PORD_INT *key)
{ register PORD_INT i;
  PORD_INT *tmp, *count, minkey, maxkey, l, u, vk;

  /* determine maximal and minimal key */
  minkey = MAX_INT;
  maxkey = 0;
  for (i = 0; i < n; i++)
   { u = node[i];
     maxkey = max(key[u], maxkey);
     minkey = min(key[u], minkey);
   }
  l = maxkey-minkey;
  /* printf("minkey %d, maxkey %d, range %d\n", minkey, maxkey, l); */
  mymalloc(count, (l+1), PORD_INT);
  mymalloc(tmp, n, PORD_INT);
  for (i = 0; i <= l; i++)
    count[i] = 0;

  /* scale down all key-values */
  for (i = 0; i < n; i++)
   { u = node[i];
     vk = key[u]-minkey;
     key[u] = vk;
     count[vk]++;
   }

  /* now do the sorting */
  for (i = 1; i <= l; i++)
    count[i] += count[i-1];
  for (i = n-1; i >= 0; i--)
   { u = node[i];
     tmp[--count[key[u]]] = u;
   }
  for (i = 0; i < n; i++)
    node[i] = tmp[i];
/*
  for (i = 0; i < n; i++)
   { u = node[i];
     printf("  node %d, key %d\n", u, key[u]);
   }
*/
  free(count);
  free(tmp);
} 

