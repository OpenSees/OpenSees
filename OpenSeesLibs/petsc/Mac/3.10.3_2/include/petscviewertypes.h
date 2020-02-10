/*
     PetscViewers are objects where other objects can be looked at or stored.
*/

#if !defined(_PETSCVIEWERTYPES_H)
#define _PETSCVIEWERTYPES_H

/*S
     PetscViewer - Abstract PETSc object that helps view (in ASCII, binary, graphically etc)
         other PETSc objects

   Level: beginner

  Concepts: viewing

.seealso:  PetscViewerCreate(), PetscViewerSetType(), PetscViewerType
S*/
typedef struct _p_PetscViewer* PetscViewer;

#endif
