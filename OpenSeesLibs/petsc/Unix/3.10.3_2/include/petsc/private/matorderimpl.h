#ifndef __MATORDERIMPL_H
#define __MATORDERIMPL_H

#include <petscmat.h>
#include <petsc/private/petscimpl.h>

/*
   Defines the interface to the SparsePack routines, translated into C.
*/
PETSC_EXTERN PetscErrorCode SPARSEPACKgen1wd(const PetscInt*,const PetscInt*,const PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode SPARSEPACKgennd(const PetscInt*,const PetscInt*,const PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode SPARSEPACKgenrcm(const PetscInt*,const PetscInt*,const PetscInt*,PetscInt*,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode SPARSEPACKgenqmd(const PetscInt*,const PetscInt*,const PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);

PETSC_EXTERN PetscErrorCode SPARSEPACKqmdrch(const PetscInt*,const PetscInt*,const PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode SPARSEPACKqmdmrg(const PetscInt*,const PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode SPARSEPACKqmdqt(const PetscInt*,const PetscInt*,const PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode SPARSEPACKqmdupd(const PetscInt*,const PetscInt*,const PetscInt*,const PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode SPARSEPACKfnroot(PetscInt*,const PetscInt*,const PetscInt*, PetscInt*, PetscInt*, PetscInt*, PetscInt*);
PETSC_EXTERN PetscErrorCode SPARSEPACKfn1wd(PetscInt*,const PetscInt*,const PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode SPARSEPACKrevrse(const PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode SPARSEPACKrootls(const PetscInt*,const PetscInt*,const PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode SPARSEPACKfndsep(PetscInt*,const PetscInt*,const PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode SPARSEPACKdegree(const PetscInt*,const PetscInt*,const PetscInt*, PetscInt*, PetscInt*, PetscInt*, PetscInt*);
PETSC_EXTERN PetscErrorCode SPARSEPACKrcm(const PetscInt*,const PetscInt*,const PetscInt*,PetscInt*,PetscInt*,PetscInt*,PetscInt*);

#endif
