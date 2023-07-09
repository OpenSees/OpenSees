      module materialAPI
      use iso_c_binding      
      use materialTypes
      implicit none

      integer ::ISW_COMMIT, ISW_INIT, ISW_REVERT_TO_START
      integer ::ISW_REVERT, ISW_FORM_TANG_AND_RESID, ISW_FORM_MASS
      real *8 ::DBL_EPSILON

      PARAMETER (ISW_INIT = 0)
      PARAMETER (ISW_COMMIT = 1)
      PARAMETER (ISW_REVERT = 2)
      PARAMETER (ISW_FORM_TANG_AND_RESID = 3)
      PARAMETER (ISW_FORM_MASS = 4)
      PARAMETER (ISW_REVERT_TO_START = 5)
      PARAMETER (DBL_EPSILON = 1.0e-21)


c      PUBLIC OPS_GetIntInput
c      interface
c         function OPS_GetIntInput(numData, iData)
c         integer       :: OPS_GetIntInput
c         integer       :: numData
c         integer       :: iData(:)
c        end function OPS_GetIntInput
c      end interface

      PUBLIC OPS_GetIntInput
      interface
         function OPS_GetIntInput(numData, iData)
         integer       :: OPS_GetIntInput
         integer       :: numData
         integer       :: iData(*)
        end function OPS_GetIntInput
      end interface


      PUBLIC OPS_GetDoubleInput
      interface
         function OPS_GetDoubleInput(numData, dData)
         integer   :: OPS_GetDoubleInput
         integer   :: numData
         real *8   :: dData(*)
         end function OPS_GetDoubleInput
      end interface

      PUBLIC OPS_AllocateMaterial
      interface
         function OPS_AllocateMaterial(mat) 
         use materialtypes
         integer :: OPS_AllocateMaterial
         type(matObject) :: mat
         end function OPS_AllocateMaterial
      end interface

      PUBLIC OPS_Error
      interface
         function OPS_Error(msg)
         integer :: OPS_Error
         character, dimension(*) :: msg
         end function OPS_Error
      end interface

c     contains
      end module materialAPI
