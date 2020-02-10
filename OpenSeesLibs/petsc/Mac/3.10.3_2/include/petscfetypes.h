#if !defined(_PETSCFETYPES_H)
#define _PETSCFETYPES_H

/*S
  PetscSpace - PETSc object that manages a linear space, e.g. the space of d-dimensional polynomials of given degree

  Level: intermediate

  Concepts: finite element

.seealso: PetscSpaceCreate(), PetscDualSpaceCreate(), PetscSpaceSetType(), PetscSpaceType
S*/
typedef struct _p_PetscSpace *PetscSpace;

/*S
  PetscDualSpace - PETSc object that manages the dual space to a linear space, e.g. the space of evaluation functionals at the vertices of a triangle

  Level: intermediate

  Concepts: finite element

.seealso: PetscDualSpaceCreate(), PetscSpaceCreate(), PetscDualSpaceSetType(), PetscDualSpaceType
S*/
typedef struct _p_PetscDualSpace *PetscDualSpace;

/*S
  PetscFE - PETSc object that manages a finite element space, e.g. the P_1 Lagrange element

  Level: intermediate

  Concepts: finite element

.seealso: PetscFECreate(), PetscSpaceCreate(), PetscDualSpaceCreate(), PetscFESetType(), PetscFEType
S*/
typedef struct _p_PetscFE *PetscFE;

/*MC
  PetscFEJacobianType - indicates which pointwise functions should be used to fill the Jacobian matrix

  Level: intermediate

.seealso: PetscFEIntegrateJacobian()
M*/
typedef enum { PETSCFE_JACOBIAN, PETSCFE_JACOBIAN_PRE, PETSCFE_JACOBIAN_DYN } PetscFEJacobianType;

#endif
