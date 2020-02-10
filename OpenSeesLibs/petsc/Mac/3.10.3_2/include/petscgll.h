#ifndef __PETSCGLL_H
#define __PETSCGLL_H
#include <petscsys.h>

/*S
    PetscGLL  - the locations, from [-1,1] and weights of the Gauss-Lobatto-Legendre nodes of a given size

  Level: beginner

  These values are usful in implementing spectral methods based on the Gauss-Lobatto-Legendre nodes

  The array nodes[] contains the vertices of each node
  The array weights[] are the integration weights

  The mass matrix for the element corresponds to the diagonal matrix whose entries are the weights[]

  Developer Notes:
    This may eventually get merged into a more abstract or general object for managing
    integration schemes or discretization schemes.

  References: XXXX

.seealso: PetscGLLCreate(), PetscGLLDestroy(), PetscGLLView()
S*/
typedef struct {
  PetscInt  n;
  PetscReal *nodes;
  PetscReal *weights;
} PetscGLL;

/*E
  PetscGLLCreateType - algorithm used to compute the GLL nodes and weights

  Level: beginner

$  PETSCGLL_VIA_LINEARALGEBRA - compute the nodes via linear algebra
$  PETSCGLL_VIA_NEWTON - compute the nodes by solving a nonlinear equation with Newton's method

.seealso: PetscGLL, PetscGLLCreate()
E*/
typedef enum {PETSCGLL_VIA_LINEARALGEBRA,PETSCGLL_VIA_NEWTON} PetscGLLCreateType;

#endif
PETSC_EXTERN PetscErrorCode PetscGLLCreate(PetscInt,PetscGLLCreateType,PetscGLL*);
PETSC_EXTERN PetscErrorCode PetscGLLDestroy(PetscGLL*);
PETSC_EXTERN PetscErrorCode PetscGLLView(PetscGLL*,PetscViewer);
PETSC_EXTERN PetscErrorCode PetscGLLIntegrate(PetscGLL*,const PetscReal*,PetscReal*);
PETSC_EXTERN PetscErrorCode PetscGLLElementLaplacianCreate(PetscGLL*,PetscReal***);
PETSC_EXTERN PetscErrorCode PetscGLLElementLaplacianDestroy(PetscGLL*,PetscReal***);
PETSC_EXTERN PetscErrorCode PetscGLLElementGradientCreate(PetscGLL*,PetscReal***,PetscReal***);
PETSC_EXTERN PetscErrorCode PetscGLLElementGradientDestroy(PetscGLL*,PetscReal***,PetscReal***);
PETSC_EXTERN PetscErrorCode PetscGLLElementAdvectionCreate(PetscGLL*,PetscReal***);
PETSC_EXTERN PetscErrorCode PetscGLLElementAdvectionDestroy(PetscGLL*,PetscReal***);

PETSC_EXTERN PetscErrorCode PetscGLLElementMassCreate(PetscGLL*,PetscReal***);
PETSC_EXTERN PetscErrorCode PetscGLLElementMassDestroy(PetscGLL*,PetscReal***);
