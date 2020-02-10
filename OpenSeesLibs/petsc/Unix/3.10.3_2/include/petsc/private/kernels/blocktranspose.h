#if !defined(_BLOCKTRANSPOSE_H)
#define _BLOCKTRANSPOSE_H

#include <petscsys.h>
#include <petscblaslapack.h>

#define PetscKernel_A_gets_transpose_A_BODY(a,N)                      \
  int i,j;                                                            \
  for (i=0; i<N; i++) {                                               \
    for (j=i+1; j<N; j++) {                                           \
      MatScalar t = a[i*N+j];                                         \
      a[i*N+j] = a[j*N+i];                                            \
      a[j*N+i] = t;                                                   \
    }                                                                 \
  }                                                                   \
  return 0

PETSC_STATIC_INLINE PetscErrorCode PetscKernel_A_gets_transpose_A_N(MatScalar *a,PetscInt N)
{
  PetscKernel_A_gets_transpose_A_BODY(a,N);
}
#define PetscKernel_A_gets_transpose_A_DECLARE(N)                            \
  PETSC_STATIC_INLINE PetscErrorCode PetscKernel_A_gets_transpose_A_ ## N (MatScalar *a) \
  {                                                                          \
    PetscKernel_A_gets_transpose_A_BODY(a,N);                                \
  }

PetscKernel_A_gets_transpose_A_DECLARE(2)
PetscKernel_A_gets_transpose_A_DECLARE(3)
PetscKernel_A_gets_transpose_A_DECLARE(4)
PetscKernel_A_gets_transpose_A_DECLARE(5)
PetscKernel_A_gets_transpose_A_DECLARE(6)
PetscKernel_A_gets_transpose_A_DECLARE(7)
PetscKernel_A_gets_transpose_A_DECLARE(8)
PetscKernel_A_gets_transpose_A_DECLARE(9)

#endif
