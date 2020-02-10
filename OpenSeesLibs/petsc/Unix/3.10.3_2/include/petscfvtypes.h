#if !defined(_PETSCFVTYPES_H)
#define _PETSCFVTYPES_H

/*S
  PetscLimiter - PETSc object that manages a finite volume slope limiter

  Level: intermediate

  Concepts: finite volume, limiter

.seealso: PetscLimiterCreate(), PetscLimiterSetType(), PetscLimiterType
S*/
typedef struct _p_PetscLimiter *PetscLimiter;

/*S
  PetscFV - PETSc object that manages a finite volume discretization

  Level: intermediate

  Concepts: finite volume

.seealso: PetscFVCreate(), PetscFVSetType(), PetscFVType
S*/
typedef struct _p_PetscFV *PetscFV;

/*S
  PetscFVFaceGeom - Data structure (C struct) for storing information about face geometry for a finite volume method.

  Level: beginner

  Note: The components are
$  PetscReal   normal[3]   - Area-scaled normals
$  PetscReal   centroid[3] - Location of centroid (quadrature point)
$  PetscScalar grad[2][3]  - Face contribution to gradient in left and right cell

  Concepts: finite volume; geometry; unstructured mesh

.seealso: DMPlexComputeGeometryFVM()
S*/
typedef struct {
  PetscReal   normal[3];   /* Area-scaled normals */
  PetscReal   centroid[3]; /* Location of centroid (quadrature point) */
  PetscScalar grad[2][3];  /* Face contribution to gradient in left and right cell */
} PetscFVFaceGeom;

/*S
  PetscFVCellGeom - Data structure (C struct) for storing information about cell geometry for a finite volume method.

  Level: beginner

  Note: The components are
$  PetscReal   centroid[3] - The cell centroid
$  PetscReal   volume      - The cell volume

  Concepts: finite volume; geometry; unstructured mesh

.seealso: DMPlexComputeGeometryFVM()
S*/
typedef struct {
  PetscReal centroid[3];
  PetscReal volume;
} PetscFVCellGeom;

#endif
