      module elementAPI
      use iso_c_binding      
      use elementTypes
      implicit none

      integer ::ISW_COMMIT, ISW_INIT, ISW_REVERT_TO_START
      integer ::ISW_REVERT, ISW_FORM_TANG_AND_RESID, ISW_FORM_MASS
      integer ::OPS_UNIAXIAL_MATERIAL_TYPE
      real *8 ::DBL_EPSILON

      PARAMETER (ISW_INIT = 0)
      PARAMETER (ISW_COMMIT = 1)
      PARAMETER (ISW_REVERT = 2)
      PARAMETER (ISW_FORM_TANG_AND_RESID = 3)
      PARAMETER (ISW_FORM_MASS = 4)
      PARAMETER (ISW_REVERT_TO_START = 5)
      PARAMETER (DBL_EPSILON = 1.0e-21)
      PARAMETER (OPS_UNIAXIAL_MATERIAL_TYPE = 1)

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
         use elementtypes
         integer :: OPS_AllocateMaterial
         type(matObject) :: mat
         end function OPS_AllocateMaterial
      end interface

      PUBLIC OPS_GetMaterial
      interface
         function OPS_GetMaterial(matTag, matType)
         use iso_c_binding      
         use elementtypes
         type(c_ptr) :: OPS_GetMaterial
         integer     :: matTag
         integer     :: matType
         end function OPS_GetMaterial
      end interface

c         integer :: matTags(*)
c         type(c_ptr) :: matTags

      PUBLIC OPS_AllocateElement
      interface
         function OPS_AllocateElement(ele, matTags, matType) 
         use elementtypes
         use iso_c_binding      
         integer :: OPS_AllocateElement
         type(eleObject) :: ele
         integer :: matTags(*)
         integer :: matType
         end function OPS_AllocateElement
      end interface


      PUBLIC OPS_GetNodeCrd
      interface
         function OPS_GetNodeCrd(nodeTag, numData, dData)
         integer   :: OPS_GetNodeCrd
         integer   :: nodeTag
         integer   :: numData
         real *8   :: dData(*)
         end function OPS_GetNodeCrd
      end interface

      PUBLIC OPS_GetNodeDisp
      interface
         function OPS_GetNodeDisp(nodeTag, numData, dData)
         integer   :: OPS_GetNodeDisp
         integer   :: nodeTag
         integer   :: numData
         real *8   :: dData(*)
         end function OPS_GetNodeDisp
      end interface

      PUBLIC OPS_GetNodeIncrDisp
      interface
         function OPS_GetNodeIncrDisp(nodeTag, numData, dData)
         integer   :: OPS_GetNodeDisp
         integer   :: nodeTag
         integer   :: numData
         real *8   :: dData(*)
         end function OPS_GetNodeIncrDisp
      end interface

      PUBLIC OPS_GetNodeIncrDeltaDisp
      interface
         function OPS_GetNodeIncrDeltaDisp(nodeTag, numData, dData)
         integer   :: OPS_GetNodeDisp
         integer   :: nodeTag
         integer   :: numData
         real *8   :: dData(*)
         end function OPS_GetNodeIncrDeltaDisp
      end interface

      PUBLIC OPS_InvokeMaterial
      interface
         function OPS_InvokeMaterial(e, mt, md, strn, strs, tang, isw)
         use elementtypes
         integer   :: OPS_InvokeMaterial
         type(eleObject) :: e
         integer :: mt
         type(modelState) :: md
         real *8   :: strn(*)
         real *8   :: strs(*)
         real *8   :: tang(*)
         integer   :: isw
         end function OPS_InvokeMaterial
      end interface

      PUBLIC OPS_InvokeMaterialDirectly
      interface
      function OPS_InvokeMaterialDirectly(mat, md, strn, strs, tang, i)
         use elementtypes
         use iso_c_binding          
         integer   :: OPS_InvokeMaterialDirectly
         type(c_ptr) :: mat
         type(modelState) :: md
         real *8   :: strn(*)
         real *8   :: strs(*)
         real *8   :: tang(*)
         integer   :: i
         end function OPS_InvokeMaterialDirectly
      end interface

      PUBLIC OPS_InvokeMaterialDirectly2
      interface
      function OPS_InvokeMaterialDirectly2(mat, md, strn, strs, tang, i)
         use elementtypes
         use iso_c_binding          
         integer   :: OPS_InvokeMaterialDirectly2
         type(c_ptr), value :: mat
         type(modelState) :: md
         real *8   :: strn(*)
         real *8   :: strs(*)
         real *8   :: tang(*)
         integer   :: i
         end function OPS_InvokeMaterialDirectly2
      end interface

      PUBLIC OPS_Error
      interface
         function OPS_Error(msg)
         integer :: OPS_Error
         character, dimension(*) :: msg
         end function OPS_Error
      end interface

      end module elementAPI
