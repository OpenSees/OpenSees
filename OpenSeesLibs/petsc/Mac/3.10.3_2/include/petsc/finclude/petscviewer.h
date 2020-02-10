!
!  Include file for Fortran use of the PetscViewer package in PETSc
!
#if !defined (__PETSCVIEWERDEF_H)
#define __PETSCVIEWERDEF_H


#define PetscViewer type(tPetscViewer)

#define PetscViewerAndFormat PetscFortranAddr

#define PetscViewers PetscFortranAddr
#define PetscFileMode PetscEnum
#define PetscViewerType character*(80)
#define PetscViewerFormat PetscEnum

#define PETSCVIEWERSOCKET 'socket'
#define PETSCVIEWERASCII 'ascii'
#define PETSCVIEWERBINARY 'binary'
#define PETSCVIEWERSTRING 'string'
#define PETSCVIEWERDRAW 'draw'
#define PETSCVIEWERVU 'vu'
#define PETSCVIEWERMATHEMATICA 'mathematica'
#define PETSCVIEWERHDF5 'hdf5'
#define PETSCVIEWERVTK 'vtk'
#define PETSCVIEWERMATLAB 'matlab'
#define PETSCVIEWERSAWS 'saws'
#define PETSCVIEWERGLVIS 'glvis'
#define PETSCVIEWERADIOS 'adios'
#define PETSCVIEWERADIOS2 'adios2'

#endif
