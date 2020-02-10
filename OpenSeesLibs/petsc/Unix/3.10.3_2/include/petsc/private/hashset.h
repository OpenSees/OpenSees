#if !defined(_PETSC_HASHSET_H)
#define _PETSC_HASHSET_H

#include <petsc/private/hashtable.h>

/*MC
  PETSC_HASH_SET - Instantiate a PETSc hash table set type

  Synopsis:
  #include <petsc/private/hashset.h>
  PETSC_HASH_SET(HSetT, KeyType, HashFunc, EqualFunc)

  Input Parameters:
+ HSetT - The hash table set type name suffix
. KeyType - The type of entries
. HashFunc - Routine or function-like macro computing hash values from entries
- EqualFunc - Routine or function-like macro computing whether two values are equal

  Level: developer

  Concepts: hash table, set

.keywords: hash table, set
.seealso: PetscHSetT, PetscHSetTCreate()
M*/

/*S
  PetscHSetT - Hash table set

  Synopsis:
  typedef khash_t(HSetT) *PetscHSetT;

  Level: developer

  Concepts: hash table, set

.keywords: hash table, set
.seealso:  PETSC_HASH_SET(), PetscHSetTCreate()
S*/

/*MC
  PetscHSetTCreate - Create a hash table

  Synopsis:
  #include <petsc/private/hashset.h>
  PetscErrorCode PetscHSetTCreate(PetscHSetT *ht)

  Output Parameter:
. ht - The hash table

  Level: developer

  Concepts: hash table, set

.keywords: hash table, set, create
.seealso: PetscHSetTDestroy()
M*/

/*MC
  PetscHSetTDestroy - Destroy a hash table

  Synopsis:
  #include <petsc/private/hashset.h>
  PetscErrorCode PetscHSetTDestroy(PetscHSetT *ht)

  Input Parameter:
. ht - The hash table

  Level: developer

  Concepts: hash table, set

.keywords: hash table, set, destroy
.seealso: PetscHSetTCreate()
M*/

/*MC
  PetscHSetTReset - Reset a hash table

  Synopsis:
  #include <petsc/private/hashset.h>
  PetscErrorCode PetscHSetTReset(PetscHSetT ht)

  Input Parameter:
. ht - The hash table

  Level: developer

  Concepts: hash table, set

.keywords: hash table, set, reset
.seealso: PetscHSetTClear()
M*/

/*MC
  PetscHSetTDuplicate - Duplicate a hash table

  Synopsis:
  #include <petsc/private/hashset.h>
  PetscErrorCode PetscHSetTDuplicate(PetscHSetT ht,PetscHSetT *hd)

  Input Parameter:
. ht - The source hash table

  Output Parameter:
. ht - The duplicated hash table

  Level: developer

  Concepts: hash table, set

.keywords: hash table, set, duplicate
.seealso: PetscHSetTCreate()
M*/

/*MC
  PetscHSetTClear - Clear a hash table

  Synopsis:
  #include <petsc/private/hashset.h>
  PetscErrorCode PetscHSetTClear(PetscHSetT ht)

  Input Parameter:
. ht - The hash table

  Level: developer

  Concepts: hash table, set

.keywords: hash table, set, reset
.seealso: PetscHSetTReset()
M*/

/*MC
  PetscHSetTResize - Set the number of buckets in a hash table

  Synopsis:
  #include <petsc/private/hashset.h>
  PetscErrorCode PetscHSetTResize(PetscHSetT ht,PetscInt nb)

  Input Parameters:
+ ht - The hash table
- nb - The number of buckets

  Level: developer

  Concepts: hash table, set

.seealso: PetscHSetTCreate()
M*/

/*MC
  PetscHSetTGetSize - Get the number of entries in a hash table

  Synopsis:
  #include <petsc/private/hashset.h>
  PetscErrorCode PetscHSetTGetSize(PetscHSetT ht,PetscInt *n)

  Input Parameter:
. ht - The hash table

  Output Parameter:
. n - The number of entries

  Level: developer

  Concepts: hash table, set

.keywords: hash table, set, resize
.seealso: PetscHSetTResize()
M*/

/*MC
  PetscHSetTHas - Query for an entry in the hash table

  Synopsis:
  #include <petsc/private/hashset.h>
  PetscErrorCode PetscHSetTHas(PetscHSetT ht,KeyType key,PetscBool *has)

  Input Parameters:
+ ht  - The hash table
- key - The entry

  Output Parameter:
. has - Boolean indicating whether the entry is in the hash table

  Level: developer

  Concepts: hash table, set

.keywords: hash table, set, query
.seealso:  PetscHSetTAdd(), PetscHSetTDel()
M*/

/*MC
  PetscHSetTAdd - Set an entry in the hash table

  Synopsis:
  #include <petsc/private/hashset.h>
  PetscErrorCode PetscHSetTAdd(PetscHSetT ht,KeyType key)

  Input Parameters:
+ ht  - The hash table
- key - The entry

  Level: developer

  Concepts: hash table, set

.keywords: hash table, set, add
.seealso: PetscHSetTDel(), PetscHSetTHas()
M*/

/*MC
  PetscHSetTDel - Remove an entry from the hash table

  Synopsis:
  #include <petsc/private/hashset.h>
  PetscErrorCode PetscHSetTDel(PetscHSetT ht,KeyType key)

  Input Parameters:
+ ht  - The hash table
- key - The entry

  Level: developer

  Concepts: hash table, set

.keywords: hash table, set, del
.seealso: PetscHSetTAdd(), PetscHSetTHas()
M*/

/*MC
  PetscHSetTQueryAdd - Query and add an entry in the hash table

  Synopsis:
  #include <petsc/private/hashset.h>
  PetscErrorCode PetscHSetTQueryAdd(PetscHSetT ht,KeyType key,PetscBool *missing)

  Input Parameters:
+ ht  - The hash table
- key - The entry

  Output Parameter:
. missing - Boolean indicating whether the entry was missing

  Level: developer

  Concepts: hash table, set

.keywords: hash table, set, query, set
.seealso: PetscHSetTQueryDel(), PetscHSetTAdd()
M*/

/*MC
  PetscHSetTQueryDel - Query and remove an entry from the hash table

  Synopsis:
  #include <petsc/private/hashset.h>
  PetscErrorCode PetscHSetTQueryDel(PetscHSetT ht,KeyType key,PetscBool *present)

  Input Parameters:
+ ht  - The hash table
- key - The entry

  Output Parameter:
. present - Boolean indicating whether the entry was present

  Level: developer

  Concepts: hash table, set

.keywords: hash table, set, query, del
.seealso: PetscHSetTQueryAdd(), PetscHSetTDel()
M*/

/*MC
  PetscHSetTGetElems - Get all entries from a hash table

  Synopsis:
  #include <petsc/private/hashset.h>
  PetscErrorCode PetscHSetTGetElems(PetscHSetT ht,PetscInt *off,KeyType array[])

  Input Parameters:
+ ht    - The hash table
. off   - Input offset in array (usually zero)
- array - Array where to put hash table entries into

  Output Parameter:
+ off   - Output offset in array (output offset = input offset + hash table size)
- array - Array filled with the hash table entries

  Level: developer

  Concepts: hash table, set

.keywords: hash table, set, array
.seealso: PetscHSetTGetSize()
M*/

#define PETSC_HASH_SET(HashT, KeyType, HashFunc, EqualFunc)                                          \
                                                                                                     \
KHASH_INIT(HashT, KeyType, char, 0, HashFunc, EqualFunc)                                             \
                                                                                                     \
typedef khash_t(HashT) *Petsc##HashT;                                                                \
                                                                                                     \
PETSC_STATIC_INLINE PETSC_UNUSED                                                                     \
PetscErrorCode Petsc##HashT##Create(Petsc##HashT *ht)                                                \
{                                                                                                    \
  PetscFunctionBegin;                                                                                \
  PetscValidPointer(ht,1);                                                                           \
  *ht = kh_init(HashT);                                                                              \
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
  PetscFunctionBegin;                                                                                \
  PetscValidPointer(ht,1);                                                                           \
  PetscValidPointer(hd,2);                                                                           \
  *hd = kh_init(HashT);                                                                              \
  ret = kh_resize(HashT,*hd,kh_size(ht));                                                            \
  PetscHashAssert(ret==0);                                                                           \
  kh_foreach_key(ht,key,{                                                                            \
      kh_put(HashT,*hd,key,&ret);                                                                    \
      PetscHashAssert(ret>=0);})                                                                     \
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
  PetscHashAssert(ret==0);                                                                           \
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
PetscErrorCode Petsc##HashT##Add(Petsc##HashT ht,KeyType key)                                        \
{                                                                                                    \
  int      ret;                                                                                      \
  khiter_t iter;                                                                                     \
  PetscFunctionBeginHot;                                                                             \
  PetscValidPointer(ht,1);                                                                           \
  iter = kh_put(HashT,ht,key,&ret); (void)iter;                                                      \
  PetscHashAssert(ret>=0);                                                                           \
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
PetscErrorCode Petsc##HashT##QueryAdd(Petsc##HashT ht,KeyType key,PetscBool *missing)                \
{                                                                                                    \
  int      ret;                                                                                      \
  khiter_t iter;                                                                                     \
  PetscFunctionBeginHot;                                                                             \
  PetscValidPointer(ht,1);                                                                           \
  PetscValidPointer(missing,3);                                                                      \
  iter = kh_put(HashT,ht,key,&ret); (void)iter;                                                      \
  PetscHashAssert(ret>=0);                                                                           \
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
PetscErrorCode Petsc##HashT##GetElems(Petsc##HashT ht,PetscInt *off,KeyType array[])                 \
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

#endif /* _PETSC_HASHSET_H */
