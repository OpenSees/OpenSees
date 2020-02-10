
#if !defined(_PETSCFPIMPL_H)
#define _PETSCFPIMPL_H
#include <petscviewertypes.h>
#include <petscsys.h>
/*
    Function pointer table that maps from function pointers to their string representation

    Does not use the PetscFunctionBegin/Return() because these routines are called within those macros
*/
typedef struct _n_PetscFPT* PetscFPT;
struct _n_PetscFPT {
  void     **functionpointer;
  char     **functionname;
  PetscInt count;
  PetscInt tablesize;
};
PETSC_INTERN PetscFPT PetscFPTData;

PETSC_STATIC_INLINE PetscErrorCode  PetscFPTView(PetscViewer viewer)
{
  PetscInt       i;

  if (!PetscFPTData) return(0);
  for (i=0; i<PetscFPTData->tablesize; i++) {
    if (PetscFPTData->functionpointer[i]) {
      printf("%s()\n",PetscFPTData->functionname[i]);
    }
  }
  return(0);
}

PETSC_STATIC_INLINE PetscErrorCode  PetscFPTDestroy(void)
{
  PetscErrorCode ierr;
  PetscFPT       _PetscFPTData = PetscFPTData;

  PetscFPTData = NULL;
  if (!_PetscFPTData) return 0;
  ierr = PetscFree((_PetscFPTData)->functionpointer);CHKERRQ(ierr);
  ierr = PetscFree((_PetscFPTData)->functionname);CHKERRQ(ierr);
  ierr = PetscFree(_PetscFPTData);CHKERRQ(ierr);
  return(0);
}

/*
   PetscFPTCreate  Creates a PETSc look up table from function pointers to strings

   Input Parameters:
.     n - expected number of keys

*/
PETSC_STATIC_INLINE PetscErrorCode  PetscFPTCreate(PetscInt n)
{
  PetscErrorCode ierr;
  PetscInt       i;
  PetscFPT       _PetscFPTData;

  if (n < 0) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"n < 0");
  /* Cannot use PetscNew() here because it is not yet defined in the include file chain */
  ierr          = PetscMalloc(sizeof(struct _n_PetscFPT),&_PetscFPTData);CHKERRQ(ierr);
  _PetscFPTData->tablesize = (3*n)/2 + 17;
  if (_PetscFPTData->tablesize < n) _PetscFPTData->tablesize = PETSC_MAX_INT/4; /* overflow */
  ierr          = PetscMalloc(sizeof(void*)*_PetscFPTData->tablesize,&_PetscFPTData->functionpointer);CHKERRQ(ierr);
  for (i=0; i<_PetscFPTData->tablesize; i++) {
    _PetscFPTData->functionpointer[i] = NULL;
  }
  ierr          = PetscMalloc(sizeof(char**)*_PetscFPTData->tablesize,&_PetscFPTData->functionname);CHKERRQ(ierr);
  _PetscFPTData->count     = 0;
  PetscFPTData = _PetscFPTData;
  return(0);
}

PETSC_STATIC_INLINE unsigned long PetscHashPointer(void *ptr)
{
#define PETSC_FPT_HASH_FACT 79943
  return((PETSC_FPT_HASH_FACT*((size_t)ptr))%PetscFPTData->tablesize);
}

PETSC_STATIC_INLINE PetscErrorCode PetscFPTAdd(void* key,const char* data)
{
  PetscInt       i,hash;

  if (!data) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Null function name");
  if (!PetscFPTData) return(0);
  hash = (PetscInt)PetscHashPointer(key);
  for (i=0; i<PetscFPTData->tablesize; i++) {
    if (PetscFPTData->functionpointer[hash] == key) {
      PetscFPTData->functionname[hash] = (char*) data;
      return(0);
    } else if (!PetscFPTData->functionpointer[hash]) {
      PetscFPTData->count++;
      PetscFPTData->functionpointer[hash] = key;
      PetscFPTData->functionname[hash] = (char*) data;
      return(0);
    }
    hash = (hash == (PetscFPTData->tablesize-1)) ? 0 : hash+1;
  }
  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_PLIB,"Function pointer table is full");
  return(0);
}

/*
    PetscFPTFind - checks if a function pointer is in the table

    If data==0, then no entry exists

*/
PETSC_STATIC_INLINE PetscErrorCode  PetscFPTFind(void* key,char const **data)
{
  PetscInt hash,ii = 0;

  *data = 0;
  if (!PetscFPTData) return(0);
  hash  = PetscHashPointer(key);
  while (ii++ < PetscFPTData->tablesize) {
    if (!PetscFPTData->functionpointer[hash]) break;
    else if (PetscFPTData->functionpointer[hash] == key) {
      *data = PetscFPTData->functionname[hash];
      break;
    }
    hash = (hash == (PetscFPTData->tablesize-1)) ? 0 : hash+1;
  }
  return(0);
}

#endif
