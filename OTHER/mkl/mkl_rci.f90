!===============================================================================
! Copyright 2005-2018 Intel Corporation All Rights Reserved.
!
! The source code,  information  and material  ("Material") contained  herein is
! owned by Intel Corporation or its  suppliers or licensors,  and  title to such
! Material remains with Intel  Corporation or its  suppliers or  licensors.  The
! Material  contains  proprietary  information  of  Intel or  its suppliers  and
! licensors.  The Material is protected by  worldwide copyright  laws and treaty
! provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
! modified, published,  uploaded, posted, transmitted,  distributed or disclosed
! in any way without Intel's prior express written permission.  No license under
! any patent,  copyright or other  intellectual property rights  in the Material
! is granted to  or  conferred  upon  you,  either   expressly,  by implication,
! inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
! property rights must be express and approved by Intel in writing.
!
! Unless otherwise agreed by Intel in writing,  you may not remove or alter this
! notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
! suppliers or licensors in any way.
!===============================================================================

!  Content:
!    Intel(R) Math Kernel Library (Intel(R) MKL) Fortran 90 interface for preconditioners,
!    RCI ISS and TR solvers routines
!******************************************************************************

MODULE MKL_RCI_TYPE

TYPE , PUBLIC :: HANDLE_TR
    INTEGER N
    INTEGER M
    INTEGER*8 HANDLE
END TYPE HANDLE_TR

END MODULE

MODULE MKL_RCI

      USE MKL_RCI_TYPE
! PARAMETERS
      INTEGER(KIND=4), PARAMETER :: TR_SUCCESS        = 1501
      INTEGER(KIND=4), PARAMETER :: TR_INVALID_OPTION = 1502
      INTEGER(KIND=4), PARAMETER :: TR_OUT_OF_MEMORY  = 1503

! SUBROUTINES

      INTERFACE
       SUBROUTINE DCG(n,x,b,rci_request,ipar,dpar,tmp)
       INTEGER, INTENT(IN)  :: n
       INTEGER, INTENT(OUT) :: rci_request
       INTEGER, INTENT(INOUT), DIMENSION(*) :: ipar
       DOUBLE PRECISION, INTENT(INOUT), DIMENSION(*)   :: x
       DOUBLE PRECISION, INTENT(IN),    DIMENSION(*)   :: b
       DOUBLE PRECISION, INTENT(INOUT),   DIMENSION(*)   :: dpar
       DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n,*) :: tmp
       END SUBROUTINE
      END INTERFACE

      INTERFACE
       SUBROUTINE DCG_INIT(n,x,b,rci_request,ipar,dpar,tmp)
       INTEGER, INTENT(IN)  :: n
       INTEGER, INTENT(OUT) :: rci_request
       INTEGER, INTENT(OUT), DIMENSION(*) :: ipar
       DOUBLE PRECISION, INTENT(IN),  DIMENSION(*)   :: x
       DOUBLE PRECISION, INTENT(IN),  DIMENSION(*)   :: b
       DOUBLE PRECISION, INTENT(OUT), DIMENSION(*)   :: dpar
       DOUBLE PRECISION, INTENT(OUT), DIMENSION(n,*) :: tmp
       END SUBROUTINE
      END INTERFACE

      INTERFACE
       SUBROUTINE DCG_CHECK(n,x,b,rci_request,ipar,dpar,tmp)
       INTEGER, INTENT(IN)  :: n
       INTEGER, INTENT(OUT) :: rci_request
       INTEGER, INTENT(INOUT), DIMENSION(*) :: ipar
       DOUBLE PRECISION, INTENT(IN),  DIMENSION(*)   :: x
       DOUBLE PRECISION, INTENT(IN),  DIMENSION(*)   :: b
       DOUBLE PRECISION, INTENT(OUT), DIMENSION(*)   :: dpar
       DOUBLE PRECISION, INTENT(OUT), DIMENSION(n,*) :: tmp
       END SUBROUTINE
      END INTERFACE

      INTERFACE
       SUBROUTINE DCG_GET(n,x,b,rci_request,ipar,dpar,tmp,itercount)
       INTEGER, INTENT(IN)  :: n
       INTEGER, INTENT(IN)  :: rci_request
       INTEGER, INTENT(OUT) :: itercount
       INTEGER, INTENT(IN), DIMENSION(*) :: ipar
       DOUBLE PRECISION, INTENT(IN),  DIMENSION(*)   :: x
       DOUBLE PRECISION, INTENT(IN),  DIMENSION(*)   :: b
       DOUBLE PRECISION, INTENT(IN),  DIMENSION(*)   :: dpar
       DOUBLE PRECISION, INTENT(IN),  DIMENSION(n,*) :: tmp
       END SUBROUTINE
      END INTERFACE

      INTERFACE
       SUBROUTINE DCSRILU0(n,a,ia,ja,alu,ipar,dpar,ierr)
       INTEGER, INTENT(IN)  :: n
       INTEGER, INTENT(OUT) :: ierr
       INTEGER, INTENT(IN), DIMENSION(*) :: ia, ja, ipar
       DOUBLE PRECISION, INTENT(IN),  DIMENSION(*) :: a, dpar
       DOUBLE PRECISION, INTENT(OUT), DIMENSION(*) :: alu
       END SUBROUTINE
      END INTERFACE

      INTERFACE
       SUBROUTINE DCSRILUT(n,a,ia,ja,alut,ialut,jalut,tol,maxfil,ipar,dpar,ierr)
       INTEGER, INTENT(IN)  :: n, maxfil
       INTEGER, INTENT(OUT) :: ierr
       DOUBLE PRECISION, INTENT(IN) :: tol
       INTEGER, INTENT(IN),  DIMENSION(*) :: ia, ja, ipar
       INTEGER, INTENT(OUT), DIMENSION(*) :: jalut, ialut
       DOUBLE PRECISION, INTENT(OUT), DIMENSION(*) :: alut
       DOUBLE PRECISION, INTENT(IN),  DIMENSION(*) :: a, dpar
       END SUBROUTINE
      END INTERFACE

      INTERFACE
       SUBROUTINE DFGMRES(n,x,b,rci_request,ipar,dpar,tmp)
       INTEGER, INTENT(IN)  :: n
       INTEGER, INTENT(OUT) :: rci_request
       INTEGER, INTENT(INOUT), DIMENSION(*) :: ipar
       DOUBLE PRECISION, INTENT(INOUT),    DIMENSION(*) :: x
       DOUBLE PRECISION, INTENT(INOUT),    DIMENSION(*) :: b
       DOUBLE PRECISION, INTENT(INOUT), DIMENSION(*) :: dpar
       DOUBLE PRECISION, INTENT(INOUT), DIMENSION(*) :: tmp
       END SUBROUTINE
      END INTERFACE

      INTERFACE
       SUBROUTINE DFGMRES_INIT(n,x,b,rci_request,ipar,dpar,tmp)
       INTEGER, INTENT(IN)  :: n
       INTEGER, INTENT(OUT) :: rci_request
       INTEGER, INTENT(OUT), DIMENSION(*) :: ipar
       DOUBLE PRECISION, INTENT(IN),    DIMENSION(*) :: x
       DOUBLE PRECISION, INTENT(IN),    DIMENSION(*) :: b
       DOUBLE PRECISION, INTENT(OUT),   DIMENSION(*) :: dpar
       DOUBLE PRECISION, INTENT(OUT),   DIMENSION(*) :: tmp
       END SUBROUTINE
      END INTERFACE

      INTERFACE
       SUBROUTINE DFGMRES_CHECK(n,x,b,rci_request,ipar,dpar,tmp)
       INTEGER, INTENT(IN)  :: n
       INTEGER, INTENT(OUT) :: rci_request
       INTEGER, INTENT(INOUT), DIMENSION(*) :: ipar
       DOUBLE PRECISION, INTENT(IN) ,    DIMENSION(*) :: x
       DOUBLE PRECISION, INTENT(IN) ,    DIMENSION(*) :: b
       DOUBLE PRECISION, INTENT(INOUT),   DIMENSION(*) :: dpar
       DOUBLE PRECISION, INTENT(OUT)  ,   DIMENSION(*) :: tmp
       END SUBROUTINE
      END INTERFACE

      INTERFACE
       SUBROUTINE DFGMRES_GET(n,x,b,rci_request,ipar,dpar,tmp,itercount)
       INTEGER, INTENT(IN)  :: n
       INTEGER, INTENT(OUT) :: rci_request
       INTEGER, INTENT(OUT) :: itercount
       INTEGER, INTENT(IN), DIMENSION(*) :: ipar
       DOUBLE PRECISION, INTENT(OUT),   DIMENSION(*) :: x
       DOUBLE PRECISION, INTENT(OUT),   DIMENSION(*) :: b
       DOUBLE PRECISION, INTENT(IN),    DIMENSION(*) :: dpar
       DOUBLE PRECISION, INTENT(IN),    DIMENSION(*) :: tmp
       END SUBROUTINE
      END INTERFACE

      INTERFACE
       SUBROUTINE DCGMRHS(n,x,nRhs,b,rci_request,ipar,dpar,tmp)
       INTEGER, INTENT(IN)  :: n
       INTEGER, INTENT(IN)  :: nRhs
       INTEGER, INTENT(OUT) :: rci_request
       INTEGER, INTENT(INOUT), DIMENSION(*) :: ipar
       DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n,*)    :: x
       DOUBLE PRECISION, INTENT(IN),    DIMENSION(nRhs,*) :: b
       DOUBLE PRECISION, INTENT(INOUT), DIMENSION(n,*)    :: tmp
       DOUBLE PRECISION, INTENT(IN),    DIMENSION(*)      :: dpar
       END SUBROUTINE
      END INTERFACE

      INTERFACE
       SUBROUTINE DCGMRHS_INIT(n,x,nRhs,b,method,rci_request,ipar,dpar,tmp)
       INTEGER, INTENT(IN)  :: n
       INTEGER, INTENT(IN)  :: nRhs
       INTEGER, INTENT(IN)  :: method
       INTEGER, INTENT(OUT) :: rci_request
       INTEGER, INTENT(OUT), DIMENSION(*) :: ipar
       DOUBLE PRECISION, INTENT(IN),  DIMENSION(n,*)    :: x
       DOUBLE PRECISION, INTENT(IN),  DIMENSION(nRhs,*) :: b
       DOUBLE PRECISION, INTENT(OUT), DIMENSION(n,*)    :: tmp
       DOUBLE PRECISION, INTENT(OUT), DIMENSION(*)      :: dpar
       END SUBROUTINE
      END INTERFACE

      INTERFACE
       SUBROUTINE DCGMRHS_CHECK(n,x,nRhs,b,rci_request,ipar,dpar,tmp)
       INTEGER, INTENT(IN)  :: n
       INTEGER, INTENT(IN)  :: nRhs
       INTEGER, INTENT(OUT) :: rci_request
       INTEGER, INTENT(INOUT), DIMENSION(*) :: ipar
       DOUBLE PRECISION, INTENT(IN),  DIMENSION(n,*)    :: x
       DOUBLE PRECISION, INTENT(IN),  DIMENSION(nRhs,*) :: b
       DOUBLE PRECISION, INTENT(OUT), DIMENSION(n,*)    :: tmp
       DOUBLE PRECISION, INTENT(INOUT), DIMENSION(*)    :: dpar
       END SUBROUTINE
      END INTERFACE

      INTERFACE
       SUBROUTINE DCGMRHS_GET(n,x,nRhs,b,rci_request,ipar,dpar,tmp,itercount)
       INTEGER, INTENT(IN)  :: n
       INTEGER, INTENT(IN)  :: nRhs
       INTEGER, INTENT(IN)  :: rci_request
       INTEGER, INTENT(OUT) :: itercount
       INTEGER, INTENT(IN), DIMENSION(*) :: ipar
       DOUBLE PRECISION, INTENT(IN),  DIMENSION(n,*)    :: x
       DOUBLE PRECISION, INTENT(IN),  DIMENSION(nRhs,*) :: b
       DOUBLE PRECISION, INTENT(IN),  DIMENSION(n,*)    :: tmp
       DOUBLE PRECISION, INTENT(IN),  DIMENSION(*)      :: dpar
       END SUBROUTINE
      END INTERFACE

! FUNCTIONS

      INTERFACE
        INTEGER FUNCTION DTRNLSP_INIT(handle, n, m, x, eps, iter1, iter2, rs)
        USE MKL_RCI_TYPE
        INTEGER, INTENT(IN)   :: n
        INTEGER, INTENT(IN)   :: m
        INTEGER, INTENT(IN)   :: iter1
        INTEGER, INTENT(IN)   :: iter2
        DOUBLE PRECISION, INTENT(IN), DIMENSION(*)    :: x
        DOUBLE PRECISION, INTENT(IN), DIMENSION(*)    :: eps
        DOUBLE PRECISION, INTENT(IN)    :: rs
        TYPE(HANDLE_TR), INTENT(OUT) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION DTRNLSP_CHECK(handle, n, m, fjac, fvec, eps, info)
        USE MKL_RCI_TYPE
        INTEGER, INTENT(IN)   :: n
        INTEGER, INTENT(IN)   :: m
        INTEGER, INTENT(OUT), DIMENSION(*):: info
        DOUBLE PRECISION, INTENT(IN), DIMENSION(*)    :: fvec
        DOUBLE PRECISION, INTENT(IN), DIMENSION(*)    :: eps
        DOUBLE PRECISION, INTENT(IN), DIMENSION(m, *) :: fjac
        TYPE(HANDLE_TR), INTENT(IN) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION DTRNLSP_SOLVE(handle, fvec, fjac, rci_request)
        USE MKL_RCI_TYPE
        INTEGER, INTENT(OUT)  :: rci_request
        DOUBLE PRECISION, INTENT(INOUT), DIMENSION(*) :: fvec
        TYPE(HANDLE_TR), INTENT(IN) :: handle
        DOUBLE PRECISION, INTENT(IN), DIMENSION(handle%m, *) :: fjac
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION DTRNLSP_GET(handle, iter, st_cr, r1, r2)
        USE MKL_RCI_TYPE
        INTEGER, INTENT(OUT)  :: iter
        INTEGER, INTENT(OUT)  :: st_cr
        DOUBLE PRECISION, INTENT(OUT)     :: r1
        DOUBLE PRECISION, INTENT(OUT)     :: r2
        TYPE(HANDLE_TR), INTENT(IN) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION DTRNLSP_DELETE(handle)
        USE MKL_RCI_TYPE
        TYPE(HANDLE_TR), INTENT(INOUT) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION DTRNLSPBC_INIT(handle, n, m, x, low, up, e, iter1, iter2, rs)
        USE MKL_RCI_TYPE
        INTEGER, INTENT(IN)   :: n
        INTEGER, INTENT(IN)   :: m
        INTEGER, INTENT(IN)   :: iter1
        INTEGER, INTENT(IN)   :: iter2
        DOUBLE PRECISION, INTENT(IN), DIMENSION(*)    :: x
        DOUBLE PRECISION, INTENT(IN), DIMENSION(*)    :: low
        DOUBLE PRECISION, INTENT(IN), DIMENSION(*)    :: up
        DOUBLE PRECISION, INTENT(IN), DIMENSION(*)    :: e
        DOUBLE PRECISION, INTENT(IN) :: rs
        TYPE(HANDLE_TR), INTENT(OUT) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION DTRNLSPBC_CHECK(handle, n, m, fjac, fvec, low, up, eps, info)
        USE MKL_RCI_TYPE
        INTEGER, INTENT(IN)   :: n
        INTEGER, INTENT(IN)   :: m
        DOUBLE PRECISION, INTENT(IN), DIMENSION(*)    :: fvec
        DOUBLE PRECISION, INTENT(IN), DIMENSION(*)    :: low
        DOUBLE PRECISION, INTENT(IN), DIMENSION(*)    :: up
        DOUBLE PRECISION, INTENT(IN), DIMENSION(m, *) :: fjac
        DOUBLE PRECISION, INTENT(IN), DIMENSION(*)    :: eps
        INTEGER, INTENT(OUT), DIMENSION(*):: info
        TYPE(HANDLE_TR), INTENT(IN) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION DTRNLSPBC_SOLVE(handle, fvec, fjac, rci_request)
        USE MKL_RCI_TYPE
        INTEGER, INTENT(OUT)   :: rci_request
        DOUBLE PRECISION, INTENT(INOUT), DIMENSION(*) :: fvec
        TYPE(HANDLE_TR), INTENT(IN) :: handle
        DOUBLE PRECISION, INTENT(IN), DIMENSION(handle%m, *) :: fjac
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION DTRNLSPBC_GET(handle, iter, st_cr, r1, r2)
        USE MKL_RCI_TYPE
        INTEGER, INTENT(OUT)   :: iter
        INTEGER, INTENT(OUT)   :: st_cr
        DOUBLE PRECISION, INTENT(OUT)      :: r1
        DOUBLE PRECISION, INTENT(OUT)      :: r2
        TYPE(HANDLE_TR), INTENT(IN) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION DTRNLSPBC_DELETE(handle)
        USE MKL_RCI_TYPE
        TYPE(HANDLE_TR), INTENT(INOUT) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION DJACOBI_INIT(handle, n, m, x, fjac, eps)
        INTEGER, INTENT(IN)   :: n
        INTEGER, INTENT(IN)   :: m
        DOUBLE PRECISION, INTENT(IN), DIMENSION(*)    :: x
        DOUBLE PRECISION, INTENT(IN), DIMENSION(m, *) :: fjac
        DOUBLE PRECISION, INTENT(IN)      :: eps
        INTEGER(KIND=8), INTENT(OUT) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION DJACOBI_SOLVE(handle, f1, f2, rci_request)
        INTEGER, INTENT(OUT)        :: rci_request
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(*)     :: f1
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(*)     :: f2
        INTEGER(KIND=8), INTENT(IN) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION DJACOBI_DELETE(handle)
        INTEGER(KIND=8), INTENT(INOUT) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION DJACOBI(fcn, n, m, fjac, x, eps)
        INTEGER, INTENT(IN)   :: n
        INTEGER, INTENT(IN)   :: m
        DOUBLE PRECISION, INTENT(IN)      :: eps
        DOUBLE PRECISION, INTENT(IN), DIMENSION(*)     :: x
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(m, *) :: fjac
        INTERFACE
            SUBROUTINE fcn(m, n, x, f)
                INTEGER, INTENT(IN)   :: n
                INTEGER, INTENT(IN)   :: m
                DOUBLE PRECISION, INTENT(IN), DIMENSION(*) :: x
                DOUBLE PRECISION, INTENT(OUT), DIMENSION(*) :: f
            END SUBROUTINE
        END INTERFACE
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION DJACOBIX(fcn, n, m, fjac, x, eps, user_data)
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER, INTENT(IN)   :: n
        INTEGER, INTENT(IN)   :: m
        DOUBLE PRECISION, INTENT(IN)      :: eps
        DOUBLE PRECISION, INTENT(IN), DIMENSION(*)     :: x
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(m, *) :: fjac
        INTEGER(C_INTPTR_T)   :: user_data
        INTERFACE
            SUBROUTINE fcn(m, n, x, f, user_data)
                USE, INTRINSIC :: ISO_C_BINDING
                INTEGER, INTENT(IN) :: n
                INTEGER, INTENT(IN) :: m
                DOUBLE PRECISION, INTENT(IN), DIMENSION(*) :: x
                DOUBLE PRECISION, INTENT(OUT), DIMENSION(*) :: f
                INTEGER(C_INTPTR_T), INTENT(IN) :: user_data
            END SUBROUTINE
        END INTERFACE
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION STRNLSP_INIT(handle, n, m, x, eps, iter1, iter2, rs)
        USE MKL_RCI_TYPE
        INTEGER, INTENT(IN)   :: n
        INTEGER, INTENT(IN)   :: m
        INTEGER, INTENT(IN)   :: iter1
        INTEGER, INTENT(IN)   :: iter2
        REAL, INTENT(IN), DIMENSION(*)    :: x
        REAL, INTENT(IN), DIMENSION(*)    :: eps
        REAL, INTENT(IN)    :: rs
        TYPE(HANDLE_TR), INTENT(OUT) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION STRNLSP_CHECK(handle, n, m, fjac, fvec, eps, info)
        USE MKL_RCI_TYPE
        INTEGER, INTENT(IN)   :: n
        INTEGER, INTENT(IN)   :: m
        INTEGER, INTENT(OUT), DIMENSION(*):: info
        REAL, INTENT(IN), DIMENSION(*)    :: fvec
        REAL, INTENT(IN), DIMENSION(*)    :: eps
        REAL, INTENT(IN), DIMENSION(m, *) :: fjac
        TYPE(HANDLE_TR), INTENT(IN) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION STRNLSP_SOLVE(handle, fvec, fjac, rci_request)
        USE MKL_RCI_TYPE
        INTEGER, INTENT(OUT)  :: rci_request
        REAL, INTENT(INOUT), DIMENSION(*) :: fvec
        TYPE(HANDLE_TR), INTENT(IN) :: handle
        REAL, INTENT(IN), DIMENSION(handle%m, *) :: fjac
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION STRNLSP_GET(handle, iter, st_cr, r1, r2)
        USE MKL_RCI_TYPE
        INTEGER, INTENT(OUT)  :: iter
        INTEGER, INTENT(OUT)  :: st_cr
        REAL, INTENT(OUT)     :: r1
        REAL, INTENT(OUT)     :: r2
        TYPE(HANDLE_TR), INTENT(IN) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION STRNLSP_DELETE(handle)
        USE MKL_RCI_TYPE
        TYPE(HANDLE_TR), INTENT(IN) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION STRNLSPBC_INIT(handle, n, m, x, low, up, e, iter1, iter2, rs)
        USE MKL_RCI_TYPE
        INTEGER, INTENT(IN)   :: n
        INTEGER, INTENT(IN)   :: m
        INTEGER, INTENT(IN)   :: iter1
        INTEGER, INTENT(IN)   :: iter2
        REAL, INTENT(IN), DIMENSION(*)    :: x
        REAL, INTENT(IN), DIMENSION(*)    :: low
        REAL, INTENT(IN), DIMENSION(*)    :: up
        REAL, INTENT(IN), DIMENSION(*)    :: e
        REAL, INTENT(IN)    :: rs
        TYPE(HANDLE_TR), INTENT(OUT) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION STRNLSPBC_CHECK(handle, n, m, fjac, fvec, low, up, eps, info)
        USE MKL_RCI_TYPE
        INTEGER, INTENT(IN)   :: n
        INTEGER, INTENT(IN)   :: m
        REAL, INTENT(IN), DIMENSION(*)    :: fvec
        REAL, INTENT(IN), DIMENSION(*)    :: low
        REAL, INTENT(IN), DIMENSION(*)    :: up
        REAL, INTENT(IN), DIMENSION(m, *) :: fjac
        REAL, INTENT(IN), DIMENSION(*)    :: eps
        INTEGER, INTENT(OUT), DIMENSION(*):: info
        TYPE(HANDLE_TR), INTENT(OUT) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION STRNLSPBC_SOLVE(handle, fvec, fjac, rci_request)
        USE MKL_RCI_TYPE
        INTEGER, INTENT(OUT)   :: rci_request
        REAL, INTENT(INOUT), DIMENSION(*) :: fvec
        TYPE(HANDLE_TR), INTENT(IN) :: handle
        REAL, INTENT(IN), DIMENSION(handle%m, *) :: fjac
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION STRNLSPBC_GET(handle, iter, st_cr, r1, r2)
        USE MKL_RCI_TYPE
        INTEGER, INTENT(OUT)   :: iter
        INTEGER, INTENT(OUT)   :: st_cr
        REAL, INTENT(OUT)      :: r1
        REAL, INTENT(OUT)      :: r2
        TYPE(HANDLE_TR), INTENT(IN) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION STRNLSPBC_DELETE(handle)
        USE MKL_RCI_TYPE
        TYPE(HANDLE_TR), INTENT(IN) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION SJACOBI_INIT(handle, n, m, x, fjac, eps)
        INTEGER, INTENT(IN)   :: n
        INTEGER, INTENT(IN)   :: m
        REAL, INTENT(IN), DIMENSION(*)    :: x
        REAL, INTENT(IN), DIMENSION(m, *) :: fjac
        REAL, INTENT(IN)      :: eps
        INTEGER(KIND=8), INTENT(OUT) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION SJACOBI_SOLVE(handle, f1, f2, rci_request)
        INTEGER, INTENT(OUT)        :: rci_request
        REAL, INTENT(OUT), DIMENSION(*)     :: f1
        REAL, INTENT(OUT), DIMENSION(*)     :: f2
        INTEGER(KIND=8), INTENT(IN) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION SJACOBI_DELETE(handle)
        INTEGER(KIND=8), INTENT(IN) :: handle
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION SJACOBI(fcn, n, m, fjac, x, eps)
        INTEGER, INTENT(IN)   :: n
        INTEGER, INTENT(IN)   :: m
        REAL, INTENT(IN)      :: eps
        REAL, INTENT(IN), DIMENSION(*)     :: x
        REAL, INTENT(OUT), DIMENSION(m, *) :: fjac
        INTERFACE
            SUBROUTINE fcn(m, n, x, f)
                INTEGER, INTENT(IN)   :: n
                INTEGER, INTENT(IN)   :: m
                REAL, INTENT(IN), DIMENSION(*) :: x
                REAL, INTENT(OUT), DIMENSION(*) :: f
            END SUBROUTINE
        END INTERFACE
        END FUNCTION
      END INTERFACE

      INTERFACE
        INTEGER FUNCTION SJACOBIX(fcn, n, m, fjac, x, eps, user_data)
        USE, INTRINSIC :: ISO_C_BINDING
        INTEGER, INTENT(IN)   :: n
        INTEGER, INTENT(IN)   :: m
        REAL, INTENT(IN)      :: eps
        REAL, INTENT(IN), DIMENSION(*)     :: x
        REAL, INTENT(OUT), DIMENSION(m, *) :: fjac
        INTEGER(C_INTPTR_T)   :: user_data
        INTERFACE
            SUBROUTINE fcn(m, n, x, f, user_data)
                USE, INTRINSIC :: ISO_C_BINDING
                INTEGER, INTENT(IN) :: n
                INTEGER, INTENT(IN) :: m
                REAL, INTENT(IN), DIMENSION(*) :: x
                REAL, INTENT(OUT), DIMENSION(*) :: f
                INTEGER(C_INTPTR_T), INTENT(IN) :: user_data
            END SUBROUTINE
        END INTERFACE
        END FUNCTION
      END INTERFACE

END MODULE
