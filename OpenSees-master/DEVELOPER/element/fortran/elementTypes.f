      module elementTypes
      use iso_c_binding      
      implicit none

      private

      public modelState
      type, bind(c) :: modelState
      real (C_DOUBLE) time
      real (C_DOUBLE) dt
      end type modelState
      
c      type modelState
c      real *8::time;
c      real *8::dt;
c      end type modelState

      public matObject
      type, bind(C) :: matObject
      integer(C_INT) ::tag;
      integer(C_INT) ::nParam;
      integer(C_INT) ::nState
      type (c_ptr) :: theParam;
      type (c_ptr) :: cState;
      type (c_ptr) :: tState;
      type(c_funptr) ::functPtr;
      type (c_ptr) :: matObjPtr;
      end type matObject;

c
c      real (C_DOUBLE) :: theParam;
c      real (C_DOUBLE) :: cState;
c      real (C_DOUBLE) :: tState;

c      public matObject
c      type :: matObject
c      integer::tag
c      integer::nParam
c      integer::nState
c      real *8, pointer ::theParam(:)
c      real *8, pointer ::cState(:)
c      real *8, pointer ::tState(:)
c      type(c_funptr)::functPtr
c      end type matObject

      public eleObject
      type, bind(C) :: eleObject
      integer(C_INT) :: tag;
      integer(C_INT) :: nNode;
      integer(C_INT) :: nDOF;
      integer(C_INT) :: nParam;
      integer(C_INT) :: nState;
      integer(C_INT) :: nMat;
      type (c_ptr) :: node;
      type (c_ptr) :: param;
      type (c_ptr) :: cState;
      type (c_ptr) :: tState;
      type (c_ptr) :: mats;
      type(c_funptr) ::functPtr
      end type eleObject;

      end module
