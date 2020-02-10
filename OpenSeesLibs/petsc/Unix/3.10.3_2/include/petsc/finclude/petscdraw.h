!
!
!  Include file for Fortran use of the PetscDraw package in PETSc
!

#if !defined (__PETSCDRAWDEF_H)
#define __PETSCDRAWDEF_H

#define PetscDraw PetscFortranAddr
#define PetscDrawLG PetscFortranAddr
#define PetscDrawAxis PetscFortranAddr
#define PetscDrawSP PetscFortranAddr
#define PetscDrawHG PetscFortranAddr
#define PetscDrawMesh PetscFortranAddr
#define PetscDrawButton PetscEnum
#define PetscDrawType character*(80)
#define PetscDrawMarkerType PetscEnum
#define PetscDrawBar PetscFortranAddr
!
!  types of draw context
!
#define PETSC_DRAW_X 'x'
#define PETSC_DRAW_NULL 'null'
#define PETSC_DRAW_PS 'ps'
#define PETSC_DRAW_WIN32 'win32'

#endif
