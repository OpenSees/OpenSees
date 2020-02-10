#if !defined(__PETSCVIENNACL_H)
#define __PETSCVIENNACL_H

#include <petscvec.h>

#if defined(PETSC_HAVE_CUDA)
#define VIENNACL_WITH_CUDA
#endif

#if defined(PETSC_HAVE_OPENCL)
#define VIENNACL_WITH_OPENCL
#endif

#if defined(PETSC_HAVE_OPENMP)
#define VIENNACL_WITH_OPENMP
#endif

#include <viennacl/forwards.h>
#include <viennacl/vector_proxy.hpp>
#include <viennacl/vector.hpp>

PETSC_EXTERN PetscErrorCode VecViennaCLGetArrayReadWrite(Vec v, viennacl::vector<PetscScalar> **a);
PETSC_EXTERN PetscErrorCode VecViennaCLRestoreArrayReadWrite(Vec v, viennacl::vector<PetscScalar> **a);

PETSC_EXTERN PetscErrorCode VecViennaCLGetArrayRead(Vec v, const viennacl::vector<PetscScalar> **a);
PETSC_EXTERN PetscErrorCode VecViennaCLRestoreArrayRead(Vec v, const viennacl::vector<PetscScalar> **a);

PETSC_EXTERN PetscErrorCode VecViennaCLGetArrayWrite(Vec v, viennacl::vector<PetscScalar> **a);
PETSC_EXTERN PetscErrorCode VecViennaCLRestoreArrayWrite(Vec v, viennacl::vector<PetscScalar> **a);


#endif
