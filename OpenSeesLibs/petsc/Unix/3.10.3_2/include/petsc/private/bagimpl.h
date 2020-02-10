#if !defined(__BAGIMPL)
#define __BAGIMPL
#include <petscbag.h>

#define PETSC_BAG_NAME_LENGTH 64
#define PETSC_BAG_HELP_LENGTH 128

struct _n_PetscBagItem {
  PetscDataType dtype;
  PetscInt      offset;
  PetscInt      msize;
  char          name[PETSC_BAG_NAME_LENGTH],help[PETSC_BAG_HELP_LENGTH];
  char          **list;
  PetscBool     freelist;
  PetscBagItem  next;
};

struct _n_PetscBag {
  MPI_Comm     bagcomm;
  PetscInt     bagsize;
  void         *structlocation;
  PetscInt     count;
  char         bagname[PETSC_BAG_NAME_LENGTH];
  char         baghelp[PETSC_BAG_HELP_LENGTH];
  char         *bagprefix;
  PetscBagItem bagitems;
};


#endif
