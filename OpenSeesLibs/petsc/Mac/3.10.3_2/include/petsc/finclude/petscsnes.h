!
!  Include file for Fortran use of the SNES package in PETSc
!
#if !defined (__PETSCSNESDEF_H)
#define __PETSCSNESDEF_H

#include "petsc/finclude/petscksp.h"

#define SNES type(tSNES)

#define PetscConvEst type(tPetscConvEst)

#define SNESType character*(80)
#define SNESMSType character*(80)
#define SNESConvergedReason PetscEnum
#define SNESLineSearchReason PetscEnum
#define SNESLineSearchType  character*(80)
#define MatMFFD PetscFortranAddr
#define MatMFFDType PetscFortranAddr
#define SNESLineSearch PetscFortranAddr
#define SNESLineSearchOrder PetscEnum
#define SNESNormSchedule PetscEnum
#define SNESQNType PetscEnum
#define SNESQNRestartType PetscEnum
#define SNESQNCompositionType PetscEnum
#define SNESQNScaleType PetscEnum
#define SNESNCGType PetscEnum
#define SNESNGMRESRestartType PetscEnum
#define SNESNGMRESSelectType PetscEnum

!
!  SNESType
!
#define SNESNEWTONLS     'newtonls'
#define SNESNEWTONTR     'newtontr'
#define SNESPYTHON       'python'
#define SNESNRICHARDSON  'nrichardson'
#define SNESKSPONLY      'ksponly'
#define SNESVINEWTONRSLS 'vinewtonrsls'
#define SNESVINEWTONSSLS 'vinewtonssls'
#define SNESNGMRES       'ngmres'
#define SNESQN           'qn'
#define SNESSHELL        'shell'
#define SNESNCG          'ncg'
#define SNESFAS          'fas'
#define SNESMS           'ms'

!
! SNESLineSearchType
!

#define SNESLINESEARCHBASIC 'basic'
#define SNESLINESEARCHBT    'bt'
#define SNESLINESEARCHL2    'l2'
#define SNESLINESEARCHCP    'cp'
#define SNESLINESEARCHSHELL 'shell'

!
! SNESLineSearchOrder
!

#define SNES_LINESEARCH_ORDER_LINEAR    1
#define SNES_LINESEARCH_ORDER_QUADRATIC 2
#define SNES_LINESEARCH_ORDER_CUBIC     3


!
!  SNESMSType
!
#define SNESMSEULER     'euler'
#define SNESMSM62       'm62'
#define SNESMSJAMESON83 'jameson83'
#define SNESMSVLTP21    'vltp21'
#define SNESMSVLTP31    'vltp31'
#define SNESMSVLTP41    'vltp41'
#define SNESMSVLTP51    'vltp51'
#define SNESMSVLTP61    'vltp61'



#endif
