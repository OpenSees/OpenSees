#if !defined(__DMFIELDIMPL_H)
#define      __DMFIELDIMPL_H

#include <petscdmfield.h>
#include <petsc/private/petscimpl.h>

PETSC_EXTERN PetscBool      DMFieldRegisterAllCalled;
PETSC_EXTERN PetscErrorCode DMFieldRegisterAll(void);

typedef struct _DMFieldOps *DMFieldOps;
struct _DMFieldOps {
  PetscErrorCode (*create) (DMField);
  PetscErrorCode (*destroy) (DMField);
  PetscErrorCode (*setfromoptions) (PetscOptionItems*,DMField);
  PetscErrorCode (*setup) (DMField);
  PetscErrorCode (*view) (DMField,PetscViewer);
  PetscErrorCode (*evaluate) (DMField,Vec,PetscDataType,void*,void*,void*);
  PetscErrorCode (*evaluateFE) (DMField,IS,PetscQuadrature,PetscDataType,void*,void*,void*);
  PetscErrorCode (*evaluateFV) (DMField,IS,PetscDataType,void*,void*,void*);
  PetscErrorCode (*getDegree) (DMField,IS,PetscInt *,PetscInt *);
  PetscErrorCode (*createDefaultQuadrature) (DMField,IS,PetscQuadrature*);
  PetscErrorCode (*computeFaceData) (DMField,IS,PetscQuadrature,PetscFEGeom *);
};
struct _p_DMField {
  PETSCHEADER(struct _DMFieldOps);
  DM dm;
  DMFieldContinuity continuity;
  PetscInt numComponents;
  void *data;
};

PETSC_INTERN PetscErrorCode DMFieldCreate(DM,PetscInt,DMFieldContinuity,DMField*);
#endif
