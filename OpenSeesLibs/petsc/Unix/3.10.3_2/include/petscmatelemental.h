#if !defined(__PETSCMATELEMENTAL_H)
#define __PETSCMATELEMENTAL_H

#include <petscmat.h>

#if defined(PETSC_HAVE_ELEMENTAL) && defined(__cplusplus)
#include <El.hpp>
#if defined(PETSC_USE_COMPLEX)
typedef El::Complex<PetscReal> PetscElemScalar;
#else
typedef PetscScalar PetscElemScalar;
#endif
PETSC_EXTERN PetscErrorCode PetscElementalInitializePackage(void);
PETSC_EXTERN PetscErrorCode PetscElementalFinalizePackage(void);
#endif

#endif /* __PETSCMATELEMENTAL_H */
