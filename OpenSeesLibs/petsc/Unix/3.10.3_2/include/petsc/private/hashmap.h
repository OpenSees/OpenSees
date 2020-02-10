#if !defined(_PETSC_HASHMAP_H)
#define _PETSC_HASHMAP_H

#include <petsc/private/hashtable.h>

/*MC
  PETSC_HASH_MAP - Instantiate a PETSc hash table map type

  Synopsis:
  #include <petsc/private/hashmap.h>
  PETSC_HASH_MAP(HMapT, KeyType, ValType, HashFunc, EqualFunc, DefaultValue)

  Input Parameters:
+ HMapT - The hash table map type name suffix
. KeyType - The type of keys
. ValType - The type of values
. HashFunc - Routine or function-like macro computing hash values from keys
. EqualFunc - Routine or function-like macro computing whether two values are equal
- DefaultValue - Default value to use for queries in case of missing keys

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map
.seealso: PetscHMapT, PetscHMapTCreate()
M*/

/*S
  PetscHMapT - Hash table map

  Synopsis:
  typedef khash_t(HMapT) *PetscHMapT;

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map
.seealso:  PETSC_HASH_MAP(), PetscHMapTCreate()
S*/

/*MC
  PetscHMapTCreate - Create a hash table

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTCreate(PetscHMapT *ht)

  Output Parameter:
. ht - The hash table

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map, create
.seealso: PetscHMapTDestroy()
M*/

/*MC
  PetscHMapTDestroy - Destroy a hash table

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTDestroy(PetscHMapT *ht)

  Input Parameter:
. ht - The hash table

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map, destroy
.seealso: PetscHMapTCreate()
M*/

/*MC
  PetscHMapTReset - Reset a hash table

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTReset(PetscHMapT ht)

  Input Parameter:
. ht - The hash table

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map, reset
.seealso: PetscHMapTClear()
M*/

/*MC
  PetscHMapTDuplicate - Duplicate a hash table

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTDuplicate(PetscHMapT ht,PetscHMapT *hd)

  Input Parameter:
. ht - The source hash table

  Output Parameter:
. ht - The duplicated hash table

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map, duplicate
.seealso: PetscHMapTCreate()
M*/

/*MC
  PetscHMapTClear - Clear a hash table

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTClear(PetscHMapT ht)

  Input Parameter:
. ht - The hash table

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map, reset
.seealso: PetscHMapTReset()
M*/

/*MC
  PetscHMapTResize - Set the number of buckets in a hash table

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTResize(PetscHMapT ht,PetscInt nb)

  Input Parameters:
+ ht - The hash table
- nb - The number of buckets

  Level: developer

  Concepts: hash table, map

.seealso: PetscHMapTCreate()
M*/

/*MC
  PetscHMapTGetSize - Get the number of entries in a hash table

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTGetSize(PetscHMapT ht,PetscInt *n)

  Input Parameter:
. ht - The hash table

  Output Parameter:
. n - The number of entries

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map, resize
.seealso: PetscHMapTResize()
M*/

/*MC
  PetscHMapTHas - Query for a key in the hash table

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTHas(PetscHMapT ht,KeyType key,PetscBool *has)

  Input Parameters:
+ ht  - The hash table
- key - The key

  Output Parameter:
. has - Boolean indicating whether key is in the hash table

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map, query
.seealso:  PetscHMapTGet(), PetscHMapTSet(), PetscHMapTFind()
M*/

/*MC
  PetscHMapTGet - Get the value for a key in the hash table

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTGet(PetscHMapT ht,KeyType key,ValType *val)

  Input Parameters:
+ ht  - The hash table
- key - The key

  Output Parameter:
. val - The value

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map, get
.seealso:  PetscHMapTSet(), PetscHMapTIterGet()
M*/

/*MC
  PetscHMapTSet - Set a (key,value) entry in the hash table

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTSet(PetscHMapT ht,KeyType key,ValType val)

  Input Parameters:
+ ht  - The hash table
. key - The key
- val - The value

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map, set
.seealso: PetscHMapTGet(), PetscHMapTIterSet()
M*/

/*MC
  PetscHMapTDel - Remove a key and its value from the hash table

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTDel(PetscHMapT ht,KeyType key)

  Input Parameters:
+ ht  - The hash table
- key - The key

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map, del
.seealso: PetscHMapTHas(), PetscHMapTIterDel()
M*/

/*MC
  PetscHMapTQuerySet - Query and set a (key,value) entry in the hash table

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTQuerySet(PetscHMapT ht,KeyType key,ValType val,PetscBool *missing)

  Input Parameters:
+ ht  - The hash table
. key - The key
- val - The value

  Output Parameter:
. missing - Boolean indicating whether the key was missing

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map, query, set
.seealso: PetscHMapTQueryDel(), PetscHMapTSet()
M*/

/*MC
  PetscHMapTQueryDel - Query and remove a (key,value) entry from the hash table

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTQueryDel(PetscHMapT ht,KeyType key,PetscBool *present)

  Input Parameters:
+ ht  - The hash table
- key - The key

  Output Parameter:
. present - Boolean indicating whether the key was present

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map, query, del
.seealso: PetscHMapTQuerySet(), PetscHMapTDel()
M*/

/*MC
  PetscHMapTFind - Query for key in the hash table

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTFind(PetscHMapT ht,KeyType key,PetscHashIter *iter,PetscBool *found)

  Input Parameters:
+ ht  - The hash table
- key - The key

  Output Parameter:
+ iter - Iterator referencing the value for key
- found - Boolean indicating whether the key was present

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map, iterator, find
.seealso: PetscHMapTIterGet(), PetscHMapTIterDel()
M*/

/*MC
  PetscHMapTPut - Set a key in the hash table

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTPut(PetscHMapT ht,KeyType key,PetscHashIter *iter,PetscBool *missing)

  Input Parameters:
+ ht  - The hash table
- key - The key

  Output Parameter:
+ iter - Iterator referencing the value for key
- missing - Boolean indicating whether the key was missing

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map, iterator, put
.seealso: PetscHMapTIterSet(), PetscHMapTQuerySet(), PetscHMapTSet()
M*/

/*MC
  PetscHMapTIterGet - Get the value referenced by an iterator in the hash table

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTIterGet(PetscHMapT ht,PetscHashIter iter,ValType *val)

  Input Parameters:
+ ht   - The hash table
- iter - The iterator

  Output Parameter:
. val  - The value

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map, iterator, get
.seealso: PetscHMapTFind(), PetscHMapTGet()
M*/

/*MC
  PetscHMapTIterSet - Set the value referenced by an iterator in the hash

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTIterSet(PetscHMapT ht,PetscHashIter iter,ValType val)

  Input Parameters:
+ ht   - The hash table
. iter - The iterator
- val  - The value

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map, iterator, set
.seealso: PetscHMapTPut(), PetscHMapTQuerySet(), PetscHMapTSet()
M*/

/*MC
  PetscHMapTIterDel - Remove the (key,value) referenced by an iterator from the hash table

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTIterDel(PetscHMapT ht,PetscHashIter iter)

  Input Parameters:
+ ht   - The hash table
- iter - The iterator

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map, iterator, del
.seealso: PetscHMapTFind(), PetscHMapTQueryDel(), PetscHMapTDel()
M*/

/*MC
  PetscHMapTGetKeys - Get all keys from a hash table

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTGetKeys(PetscHMapT ht,PetscInt *off,KeyType array[])

  Input Parameters:
+ ht    - The hash table
. off   - Input offset in array (usually zero)
- array - Array where to put hash table keys into

  Output Parameter:
+ off   - Output offset in array (output offset = input offset + hash table size)
- array - Array filled with the hash table keys

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map, array
.seealso: PetscHSetTGetSize(), PetscHMapTGetVals()
M*/

/*MC
  PetscHMapTGetVals - Get all values from a hash table

  Synopsis:
  #include <petsc/private/hashmap.h>
  PetscErrorCode PetscHMapTGetVals(PetscHMapT ht,PetscInt *off,ValType array[])

  Input Parameters:
+ ht    - The hash table
. off   - Input offset in array (usually zero)
- array - Array where to put hash table values into

  Output Parameter:
+ off   - Output offset in array (output offset = input offset + hash table size)
- array - Array filled with the hash table values

  Level: developer

  Concepts: hash table, map

.keywords: hash table, map, array
.seealso: PetscHSetTGetSize(), PetscHMapTGetKeys()
M*/

#define PETSC_HASH_MAP(HashT, KeyType, ValType, HashFunc, EqualFunc, DefaultValue)                   \
                                                                                                     \
KHASH_INIT(HashT, KeyType, ValType, 1, HashFunc, EqualFunc)                                          \
                                                                                                     \
typedef khash_t(HashT) *Petsc##HashT;                                                                \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##Create(Petsc##HashT *ht)                                                \
{                                                                                                    \
  PetscFunctionBegin;                                                                                \
  PetscValidPointer(ht,1);                                                                           \
  *ht = kh_init(HashT);                                                                              \
  PetscHashAssert(*ht!=NULL);                                                                        \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##Destroy(Petsc##HashT *ht)                                               \
{                                                                                                    \
  PetscFunctionBegin;                                                                                \
  PetscValidPointer(ht,1);                                                                           \
  if (!*ht) PetscFunctionReturn(0);                                                                  \
  kh_destroy(HashT,*ht); *ht = NULL;                                                                 \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##Reset(Petsc##HashT ht)                                                  \
{                                                                                                    \
  PetscFunctionBegin;                                                                                \
  PetscValidPointer(ht,1);                                                                           \
  kh_reset(HashT,ht);                                                                                \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##Duplicate(Petsc##HashT ht,Petsc##HashT *hd)                             \
{                                                                                                    \
  int     ret;                                                                                       \
  KeyType key;                                                                                       \
  ValType val;                                                                                       \
  PetscFunctionBegin;                                                                                \
  PetscValidPointer(ht,1);                                                                           \
  PetscValidPointer(hd,2);                                                                           \
  *hd = kh_init(HashT);                                                                              \
  PetscHashAssert(*hd!=NULL);                                                                        \
  ret = kh_resize(HashT,*hd,kh_size(ht));                                                            \
  PetscHashAssert(ret==0);                                                                           \
  kh_foreach(ht,key,val,{ khiter_t i;                                                                \
      i = kh_put(HashT,*hd,key,&ret);                                                                \
      PetscHashAssert(ret>=0);                                                                       \
      kh_val(*hd,i) = val;})                                                                         \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##Clear(Petsc##HashT ht)                                                  \
{                                                                                                    \
  PetscFunctionBegin;                                                                                \
  PetscValidPointer(ht,1);                                                                           \
  kh_clear(HashT,ht);                                                                                \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##Resize(Petsc##HashT ht,PetscInt nb)                                     \
{                                                                                                    \
  int ret;                                                                                           \
  PetscFunctionBegin;                                                                                \
  PetscValidPointer(ht,1);                                                                           \
  ret = kh_resize(HashT,ht,(khint_t)nb);                                                             \
  PetscHashAssert(ret>=0);                                                                           \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##GetSize(Petsc##HashT ht,PetscInt *n)                                    \
{                                                                                                    \
  PetscFunctionBegin;                                                                                \
  PetscValidPointer(ht,1);                                                                           \
  PetscValidIntPointer(n,2);                                                                         \
  *n = (PetscInt)kh_size(ht);                                                                        \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##Has(Petsc##HashT ht,KeyType key,PetscBool *has)                         \
{                                                                                                    \
  khiter_t iter;                                                                                     \
  PetscFunctionBeginHot;                                                                             \
  PetscValidPointer(ht,1);                                                                           \
  PetscValidPointer(has,3);                                                                          \
  iter = kh_get(HashT,ht,key);                                                                       \
  *has = (iter != kh_end(ht)) ? PETSC_TRUE : PETSC_FALSE;                                            \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##Get(Petsc##HashT ht,KeyType key,ValType *val)                           \
{                                                                                                    \
  khiter_t iter;                                                                                     \
  PetscFunctionBeginHot;                                                                             \
  PetscValidPointer(ht,1);                                                                           \
  PetscValidIntPointer(val,3);                                                                       \
  iter = kh_get(HashT,ht,key);                                                                       \
  *val = (iter != kh_end(ht)) ? kh_val(ht,iter) : (DefaultValue);                                    \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##Set(Petsc##HashT ht,KeyType key,ValType val)                            \
{                                                                                                    \
  int      ret;                                                                                      \
  khiter_t iter;                                                                                     \
  PetscFunctionBeginHot;                                                                             \
  PetscValidPointer(ht,1);                                                                           \
  iter = kh_put(HashT,ht,key,&ret);                                                                  \
  PetscHashAssert(ret>=0);                                                                           \
  kh_val(ht,iter) = val;                                                                             \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##Del(Petsc##HashT ht,KeyType key)                                        \
{                                                                                                    \
  khiter_t iter;                                                                                     \
  PetscFunctionBeginHot;                                                                             \
  PetscValidPointer(ht,1);                                                                           \
  iter = kh_get(HashT,ht,key);                                                                       \
  kh_del(HashT,ht,iter);                                                                             \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##QuerySet(Petsc##HashT ht,KeyType key,ValType val,PetscBool *missing)    \
{                                                                                                    \
  int      ret;                                                                                      \
  khiter_t iter;                                                                                     \
  PetscFunctionBeginHot;                                                                             \
  PetscValidPointer(ht,1);                                                                           \
  PetscValidPointer(missing,3);                                                                      \
  iter = kh_put(HashT,ht,key,&ret);                                                                  \
  PetscHashAssert(ret>=0);                                                                           \
  kh_val(ht,iter) = val;                                                                             \
  *missing = ret ? PETSC_TRUE : PETSC_FALSE;                                                         \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##QueryDel(Petsc##HashT ht,KeyType key,PetscBool *present)                \
{                                                                                                    \
  khiter_t iter;                                                                                     \
  PetscFunctionBeginHot;                                                                             \
  PetscValidPointer(ht,1);                                                                           \
  PetscValidPointer(present,3);                                                                      \
  iter = kh_get(HashT,ht,key);                                                                       \
  if (iter != kh_end(ht)) {                                                                          \
    kh_del(HashT,ht,iter);                                                                           \
    *present = PETSC_TRUE;                                                                           \
  } else {                                                                                           \
    *present = PETSC_FALSE;                                                                          \
  }                                                                                                  \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##Find(Petsc##HashT ht,KeyType key,PetscHashIter *iter,PetscBool *found)  \
                                                                                                     \
{                                                                                                    \
  PetscFunctionBeginHot;                                                                             \
  PetscValidPointer(ht,1);                                                                           \
  PetscValidPointer(iter,2);                                                                         \
  PetscValidPointer(found,3);                                                                        \
  *iter = kh_get(HashT,ht,key);                                                                      \
  *found = (*iter != kh_end(ht)) ? PETSC_TRUE : PETSC_FALSE;                                         \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##Put(Petsc##HashT ht,KeyType key,PetscHashIter *iter,PetscBool *missing) \
{                                                                                                    \
  int ret;                                                                                           \
  PetscFunctionBeginHot;                                                                             \
  PetscValidPointer(ht,1);                                                                           \
  PetscValidPointer(iter,2);                                                                         \
  PetscValidPointer(missing,3);                                                                      \
  *iter = kh_put(HashT,ht,key,&ret);                                                                 \
  PetscHashAssert(ret>=0);                                                                           \
  *missing = ret ? PETSC_TRUE : PETSC_FALSE;                                                         \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##IterGet(Petsc##HashT ht,PetscHashIter iter,ValType *val)                \
{                                                                                                    \
  PetscFunctionBeginHot;                                                                             \
  PetscValidPointer(ht,1);                                                                           \
  PetscValidPointer(val,3);                                                                          \
  *val = PetscLikely(iter < kh_end(ht) && kh_exist(ht,iter)) ? kh_val(ht,iter) : (DefaultValue);     \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##IterSet(Petsc##HashT ht,PetscHashIter iter,ValType val)                 \
{                                                                                                    \
  PetscFunctionBeginHot;                                                                             \
  PetscValidPointer(ht,1);                                                                           \
  if (PetscLikely(iter < kh_end(ht) && kh_exist(ht,iter))) kh_val(ht,iter) = val;                    \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##IterDel(Petsc##HashT ht,PetscHashIter iter)                             \
{                                                                                                    \
  PetscFunctionBeginHot;                                                                             \
  PetscValidPointer(ht,1);                                                                           \
  if (PetscLikely(iter < kh_end(ht))) kh_del(HashT,ht,iter);                                         \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##GetKeys(Petsc##HashT ht,PetscInt *off,KeyType array[])                  \
{                                                                                                    \
  KeyType  key;                                                                                      \
  PetscInt pos;                                                                                      \
  PetscFunctionBegin;                                                                                \
  PetscValidPointer(ht,1);                                                                           \
  PetscValidIntPointer(off,2);                                                                       \
  pos = *off;                                                                                        \
  kh_foreach_key(ht,key,array[pos++] = key);                                                         \
  *off = pos;                                                                                        \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##GetVals(Petsc##HashT ht,PetscInt *off,ValType array[])                  \
{                                                                                                    \
  ValType  val;                                                                                      \
  PetscInt pos;                                                                                      \
  PetscFunctionBegin;                                                                                \
  PetscValidPointer(ht,1);                                                                           \
  PetscValidIntPointer(off,2);                                                                       \
  pos = *off;                                                                                        \
  kh_foreach_value(ht,val,array[pos++] = val);                                                       \
  *off = pos;                                                                                        \
  PetscFunctionReturn(0);                                                                            \
}                                                                                                    \

#endif /* _PETSC_HASHMAP_H */
