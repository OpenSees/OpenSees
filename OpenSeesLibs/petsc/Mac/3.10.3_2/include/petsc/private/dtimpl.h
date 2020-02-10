#if !defined(_DT_H)
#define _DT_H

#include <petscdt.h>

struct _p_PetscQuadrature {
  PETSCHEADER(int);
  PetscInt         dim;       /* The spatial dimension */
  PetscInt         Nc;        /* The number of components */
  PetscInt         order;     /* The order, i.e. the highest degree polynomial that is exactly integrated */
  PetscInt         numPoints; /* The number of quadrature points on an element */
  const PetscReal *points;    /* The quadrature point coordinates */
  const PetscReal *weights;   /* The quadrature weights */
};

#endif
