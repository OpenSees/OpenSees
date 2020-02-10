
!
!  Include file for Fortran use of the DM package in PETSc
!
#if !defined (__PETSCDMDEF_H)
#define __PETSCDMDEF_H

#include "petsc/finclude/petscis.h"
#include "petsc/finclude/petscvec.h"
#include "petsc/finclude/petscmat.h"

#define DMType character*(80)
#define DMBoundaryType      PetscEnum
#define DMPointLocationType PetscEnum
#define DMAdaptationType    PetscEnum
#define DMAdaptFlag         PetscEnum
#define PetscUnit           PetscEnum
#define DMAdaptationStrategy PetscEnum

#define DM               type(tDM)

#define DMAdaptor        PetscFortranAddr
#define PetscQuadrature  PetscFortranAddr
#define PetscDS          PetscFortranAddr
#define PetscFE          PetscFortranAddr
#define PetscSpace       PetscFortranAddr
#define PetscDualSpace   PetscFortranAddr
#define PetscFV          PetscFortranAddr
#define PetscLimiter     PetscFortranAddr
#define PetscPartitioner PetscFortranAddr
#define DMField          PetscFortranAddr


#define DMDA        'da'
#define DMCOMPOSITE 'composite'
#define DMSLICED    'sliced'
#define DMSHELL     'shell'
#define DMPLEX      'plex'
#define DMREDUNDANT 'redundant'
#define DMPATCH     'patch'
#define DMMOAB      'moab'
#define DMNETWORK   'network'
#define DMFOREST    'forest'
#define DMP4EST     'p4est'
#define DMP8EST     'p8est'
#define DMSWARM     'swarm'

#endif
