!===============================================================================
! Copyright 2006-2018 Intel Corporation All Rights Reserved.
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
!      Intel(R) Math Kernel Library (Intel(R) MKL) interface for TT routines
!*******************************************************************************

      MODULE  MKL_TT_TYPE

! Parameters definitions for the kind of the Trigonometric Transform
      INTEGER, PARAMETER :: MKL_SINE_TRANSFORM               = 0
      INTEGER, PARAMETER :: MKL_COSINE_TRANSFORM             = 1
      INTEGER, PARAMETER :: MKL_STAGGERED_COSINE_TRANSFORM   = 2
      INTEGER, PARAMETER :: MKL_STAGGERED_SINE_TRANSFORM     = 3
      INTEGER, PARAMETER :: MKL_STAGGERED2_COSINE_TRANSFORM  = 4
      INTEGER, PARAMETER :: MKL_STAGGERED2_SINE_TRANSFORM    = 5

      END MODULE MKL_TT_TYPE

      MODULE  MKL_TRIG_TRANSFORMS

      USE MKL_TT_TYPE
      USE MKL_DFTI

      INTERFACE

        SUBROUTINE D_INIT_TRIG_TRANSFORM(n, tt_type, ipar,dpar, stat)
            USE MKL_DFT_TYPE
            !DEC$ ATTRIBUTES C, ALIAS: '_d_init_trig_transform' :: D_INIT_TRIG_TRANSFORM
            !MS$ATTRIBUTES REFERENCE :: n
            !MS$ATTRIBUTES REFERENCE :: tt_type
            !MS$ATTRIBUTES REFERENCE :: ipar
            !MS$ATTRIBUTES REFERENCE :: dpar
            !MS$ATTRIBUTES REFERENCE :: stat
            INTEGER, INTENT(IN) :: n, tt_type
            INTEGER, INTENT(INOUT) :: ipar(*)
            REAL(8), INTENT(INOUT) :: dpar(*)
            INTEGER, INTENT(OUT) :: stat
        END SUBROUTINE D_INIT_TRIG_TRANSFORM

        SUBROUTINE D_COMMIT_TRIG_TRANSFORM(f, handle, ipar,dpar, stat)
            USE MKL_DFT_TYPE
            !DEC$ ATTRIBUTES C, ALIAS: '_d_commit_trig_transform' :: D_COMMIT_TRIG_TRANSFORM
            !MS$ATTRIBUTES REFERENCE :: f
            !MS$ATTRIBUTES REFERENCE :: handle
            !MS$ATTRIBUTES REFERENCE :: ipar
            !MS$ATTRIBUTES REFERENCE :: dpar
            !MS$ATTRIBUTES REFERENCE :: stat
            REAL(8), INTENT(INOUT) :: f(*)
            TYPE(DFTI_DESCRIPTOR), POINTER :: handle
            INTEGER, INTENT(INOUT) :: ipar(*)
            REAL(8), INTENT(OUT) :: dpar(*)
            INTEGER, INTENT(OUT) :: stat
        END SUBROUTINE D_COMMIT_TRIG_TRANSFORM

        SUBROUTINE D_FORWARD_TRIG_TRANSFORM(f, handle, ipar,dpar, stat)
            USE MKL_DFT_TYPE
            !DEC$ ATTRIBUTES C, ALIAS: '_d_forward_trig_transform' :: D_FORWARD_TRIG_TRANSFORM
            !MS$ATTRIBUTES REFERENCE :: f
            !MS$ATTRIBUTES REFERENCE :: handle
            !MS$ATTRIBUTES REFERENCE :: ipar
            !MS$ATTRIBUTES REFERENCE :: dpar
            !MS$ATTRIBUTES REFERENCE :: stat
            REAL(8), INTENT(INOUT) :: f(*)
            TYPE(DFTI_DESCRIPTOR), POINTER :: handle
            INTEGER, INTENT(INOUT) :: ipar(*)
            REAL(8), INTENT(IN) :: dpar(*)
            INTEGER, INTENT(OUT) :: stat
        END SUBROUTINE D_FORWARD_TRIG_TRANSFORM

        SUBROUTINE D_BACKWARD_TRIG_TRANSFORM(f, handle, ipar,dpar, stat)
            USE MKL_DFT_TYPE
            !DEC$ ATTRIBUTES C, ALIAS: '_d_backward_trig_transform' :: D_BACKWARD_TRIG_TRANSFORM
            !MS$ATTRIBUTES REFERENCE :: f
            !MS$ATTRIBUTES REFERENCE :: handle
            !MS$ATTRIBUTES REFERENCE :: ipar
            !MS$ATTRIBUTES REFERENCE :: dpar
            !MS$ATTRIBUTES REFERENCE :: stat
            REAL(8), INTENT(INOUT) :: f(*)
            TYPE(DFTI_DESCRIPTOR), POINTER :: handle
            INTEGER, INTENT(INOUT) :: ipar(*)
            REAL(8), INTENT(IN) :: dpar(*)
            INTEGER, INTENT(OUT) :: stat
        END SUBROUTINE D_BACKWARD_TRIG_TRANSFORM

        SUBROUTINE S_INIT_TRIG_TRANSFORM(n, tt_type, ipar,spar, stat)
            USE MKL_DFT_TYPE
            !DEC$ ATTRIBUTES C, ALIAS: '_s_init_trig_transform' :: S_INIT_TRIG_TRANSFORM
            !MS$ATTRIBUTES REFERENCE :: n
            !MS$ATTRIBUTES REFERENCE :: tt_type
            !MS$ATTRIBUTES REFERENCE :: ipar
            !MS$ATTRIBUTES REFERENCE :: spar
            !MS$ATTRIBUTES REFERENCE :: stat
            INTEGER, INTENT(IN) :: n, tt_type
            INTEGER, INTENT(INOUT) :: ipar(*)
            REAL(4), INTENT(INOUT) :: spar(*)
            INTEGER, INTENT(OUT) :: stat
        END SUBROUTINE S_INIT_TRIG_TRANSFORM

        SUBROUTINE S_COMMIT_TRIG_TRANSFORM(f, handle, ipar,spar, stat)
            USE MKL_DFT_TYPE
            !DEC$ ATTRIBUTES C, ALIAS: '_s_commit_trig_transform' :: S_COMMIT_TRIG_TRANSFORM
            !MS$ATTRIBUTES REFERENCE :: f
            !MS$ATTRIBUTES REFERENCE :: handle
            !MS$ATTRIBUTES REFERENCE :: ipar
            !MS$ATTRIBUTES REFERENCE :: spar
            !MS$ATTRIBUTES REFERENCE :: stat
            REAL(4), INTENT(INOUT) :: f(*)
            TYPE(DFTI_DESCRIPTOR), POINTER :: handle
            INTEGER, INTENT(INOUT) :: ipar(*)
            REAL(4), INTENT(OUT) :: spar(*)
            INTEGER, INTENT(OUT) :: stat
        END SUBROUTINE S_COMMIT_TRIG_TRANSFORM

        SUBROUTINE S_FORWARD_TRIG_TRANSFORM(f, handle, ipar,spar, stat)
            USE MKL_DFT_TYPE
            !DEC$ ATTRIBUTES C, ALIAS: '_s_forward_trig_transform' :: S_FORWARD_TRIG_TRANSFORM
            !MS$ATTRIBUTES REFERENCE :: f
            !MS$ATTRIBUTES REFERENCE :: handle
            !MS$ATTRIBUTES REFERENCE :: ipar
            !MS$ATTRIBUTES REFERENCE :: spar
            !MS$ATTRIBUTES REFERENCE :: stat
            REAL(4), INTENT(INOUT) :: f(*)
            TYPE(DFTI_DESCRIPTOR), POINTER :: handle
            INTEGER, INTENT(INOUT) :: ipar(*)
            REAL(4), INTENT(IN) :: spar(*)
            INTEGER, INTENT(OUT) :: stat
        END SUBROUTINE S_FORWARD_TRIG_TRANSFORM

        SUBROUTINE S_BACKWARD_TRIG_TRANSFORM(f, handle, ipar,spar, stat)
            USE MKL_DFT_TYPE
            !DEC$ ATTRIBUTES C, ALIAS: '_s_backward_trig_transform' :: S_BACKWARD_TRIG_TRANSFORM
            !MS$ATTRIBUTES REFERENCE :: f
            !MS$ATTRIBUTES REFERENCE :: handle
            !MS$ATTRIBUTES REFERENCE :: ipar
            !MS$ATTRIBUTES REFERENCE :: spar
            !MS$ATTRIBUTES REFERENCE :: stat
            REAL(4), INTENT(INOUT) :: f(*)
            TYPE(DFTI_DESCRIPTOR), POINTER :: handle
            INTEGER, INTENT(INOUT) :: ipar(*)
            REAL(4), INTENT(IN) :: spar(*)
            INTEGER, INTENT(OUT) :: stat
        END SUBROUTINE S_BACKWARD_TRIG_TRANSFORM

        SUBROUTINE FREE_TRIG_TRANSFORM(handle, ipar,stat)
            USE MKL_DFT_TYPE
            !DEC$ ATTRIBUTES C, ALIAS: '_free_trig_transform' :: FREE_TRIG_TRANSFORM
            !MS$ATTRIBUTES REFERENCE :: handle
            !MS$ATTRIBUTES REFERENCE :: ipar
            !MS$ATTRIBUTES REFERENCE :: stat
            INTEGER, INTENT(INOUT) :: ipar(*)
            TYPE(DFTI_DESCRIPTOR), POINTER :: handle
            INTEGER, INTENT(OUT) :: stat
        END SUBROUTINE FREE_TRIG_TRANSFORM

      END INTERFACE

      END MODULE MKL_TRIG_TRANSFORMS
