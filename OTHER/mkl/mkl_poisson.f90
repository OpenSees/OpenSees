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

MODULE MKL_POISSON

USE MKL_DFTI

INTERFACE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!INTERFACES FOR 3D CASE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE D_INIT_HELMHOLTZ_3D(AX,BX,AY,BY,AZ,BZ,NX,NY,NZ,BCTYPE,Q,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER NX, NY, NZ, STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION AX,BX,AY,BY,AZ,BZ,Q
   CHARACTER(6) BCTYPE
   DOUBLE PRECISION DPAR(*)
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE D_COMMIT_HELMHOLTZ_3D(F,BD_AX,BD_BX,BD_AY,BD_BY,BD_AZ,BD_BZ,XHANDLE,YHANDLE,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION DPAR(*)
   DOUBLE PRECISION F(IPAR(11)+1,IPAR(12)+1,*)
   DOUBLE PRECISION BD_AX(IPAR(12)+1,*),BD_BX(IPAR(12)+1,*),BD_AY(IPAR(11)+1,*),BD_BY(IPAR(11)+1,*)
   DOUBLE PRECISION BD_AZ(IPAR(11)+1,*),BD_BZ(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: XHANDLE, YHANDLE
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE D_HELMHOLTZ_3D(F,BD_AX,BD_BX,BD_AY,BD_BY,BD_AZ,BD_BZ,XHANDLE,YHANDLE,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION F(IPAR(11)+1,IPAR(12)+1,*)
   DOUBLE PRECISION BD_AX(IPAR(12)+1,*),BD_BX(IPAR(12)+1,*),BD_AY(IPAR(11)+1,*),BD_BY(IPAR(11)+1,*)
   DOUBLE PRECISION BD_AZ(IPAR(11)+1,*),BD_BZ(IPAR(11)+1,*)
   DOUBLE PRECISION DPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: XHANDLE, YHANDLE
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE FREE_HELMHOLTZ_3D(XHANDLE,YHANDLE,IPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: XHANDLE, YHANDLE
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!INTERFACES FOR 2D CASE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE D_INIT_HELMHOLTZ_2D(AX,BX,AY,BY,NX,NY,BCTYPE,Q,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER NX, NY, STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION AX,BX,AY,BY,Q
   CHARACTER(4) BCTYPE
   DOUBLE PRECISION DPAR(*)
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE D_COMMIT_HELMHOLTZ_2D(F,BD_AX,BD_BX,BD_AY,BD_BY,HANDLE,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION F(IPAR(11)+1,*)
   DOUBLE PRECISION BD_AX(*),BD_BX(*),BD_AY(*),BD_BY(*)
   DOUBLE PRECISION DPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE D_HELMHOLTZ_2D(F,BD_AX,BD_BX,BD_AY,BD_BY,HANDLE,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION F(IPAR(11)+1,*)
   DOUBLE PRECISION BD_AX(*),BD_BX(*),BD_AY(*),BD_BY(*)
   DOUBLE PRECISION DPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE FREE_HELMHOLTZ_2D(HANDLE,IPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE
!**********************************************************************
!***********************SINGLE*****************************************
!**********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!INTERFACES FOR 3D CASE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE S_INIT_HELMHOLTZ_3D(AX,BX,AY,BY,AZ,BZ,NX,NY,NZ,BCTYPE,Q,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER NX, NY, NZ, STAT
   INTEGER IPAR(*)
   REAL AX,BX,AY,BY,AZ,BZ,Q
   CHARACTER(6) BCTYPE
   REAL SPAR(*)
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE S_COMMIT_HELMHOLTZ_3D(F,BD_AX,BD_BX,BD_AY,BD_BY,BD_AZ,BD_BZ,XHANDLE,YHANDLE,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   REAL SPAR(*)
   REAL F(IPAR(11)+1,IPAR(12)+1,*)
   REAL BD_AX(IPAR(12)+1,*),BD_BX(IPAR(12)+1,*),BD_AY(IPAR(11)+1,*),BD_BY(IPAR(11)+1,*)
   REAL BD_AZ(IPAR(11)+1,*),BD_BZ(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: XHANDLE, YHANDLE
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE S_HELMHOLTZ_3D(F,BD_AX,BD_BX,BD_AY,BD_BY,BD_AZ,BD_BZ,XHANDLE,YHANDLE,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   REAL F(IPAR(11)+1,IPAR(12)+1,*)
   REAL BD_AX(IPAR(12)+1,*),BD_BX(IPAR(12)+1,*),BD_AY(IPAR(11)+1,*),BD_BY(IPAR(11)+1,*)
   REAL BD_AZ(IPAR(11)+1,*),BD_BZ(IPAR(11)+1,*)
   REAL SPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: XHANDLE, YHANDLE
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!INTERFACES FOR 2D CASE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


SUBROUTINE S_INIT_HELMHOLTZ_2D(AX,BX,AY,BY,NX,NY,BCTYPE,Q,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER NX, NY, STAT
   INTEGER IPAR(*)
   REAL AX,BX,AY,BY,Q
   CHARACTER(4) BCTYPE
   REAL SPAR(*)
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE S_COMMIT_HELMHOLTZ_2D(F,BD_AX,BD_BX,BD_AY,BD_BY,HANDLE,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   REAL F(IPAR(11)+1,*)
   REAL BD_AX(*),BD_BX(*),BD_AY(*),BD_BY(*)
   REAL SPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE S_HELMHOLTZ_2D(F,BD_AX,BD_BX,BD_AY,BD_BY,HANDLE,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   REAL F(IPAR(11)+1,*)
   REAL BD_AX(*),BD_BX(*),BD_AY(*),BD_BY(*)
   REAL SPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE

SUBROUTINE D_INIT_SPH_P(AX,BX,AY,BY,NX,NY,Q,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER NX, NY, STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION AX,BX,AY,BY,Q
   DOUBLE PRECISION DPAR(*)
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE D_COMMIT_SPH_P(F,HANDLE_S,HANDLE_C,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION DPAR(*)
   DOUBLE PRECISION F(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE_C, HANDLE_S
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE D_SPH_P(F,HANDLE_S,HANDLE_C,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION DPAR(*)
   DOUBLE PRECISION F(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE_C, HANDLE_S
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE FREE_SPH_P(HANDLE_S,HANDLE_C,IPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE_S, HANDLE_C
END SUBROUTINE

!---------------------------------------------------------------------
SUBROUTINE D_INIT_SPH_NP(AX,BX,AY,BY,NX,NY,Q,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER NX, NY, STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION AX,BX,AY,BY,Q
   DOUBLE PRECISION DPAR(*)
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE D_COMMIT_SPH_NP(F,HANDLE,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION DPAR(*)
   DOUBLE PRECISION F(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE D_SPH_NP(F,HANDLE,IPAR,DPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   DOUBLE PRECISION DPAR(*)
   DOUBLE PRECISION F(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE FREE_SPH_NP(HANDLE,IPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE

!======================================================================

SUBROUTINE S_INIT_SPH_P(AX,BX,AY,BY,NX,NY,Q,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER NX, NY, STAT
   INTEGER IPAR(*)
   REAL AX,BX,AY,BY,Q
   REAL SPAR(*)
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE S_COMMIT_SPH_P(F,HANDLE_S,HANDLE_C,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   REAL SPAR(*)
   REAL F(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE_C, HANDLE_S
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE S_SPH_P(F,HANDLE_S,HANDLE_C,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   REAL SPAR(*)
   REAL F(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE_C, HANDLE_S
END SUBROUTINE

!--------------------------------------------------------------------

SUBROUTINE S_INIT_SPH_NP(AX,BX,AY,BY,NX,NY,Q,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER NX, NY, STAT
   INTEGER IPAR(*)
   REAL AX,BX,AY,BY,Q
   REAL SPAR(*)
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE S_COMMIT_SPH_NP(F,HANDLE,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   REAL SPAR(*)
   REAL F(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE

!---------------------------------------------------------------------

SUBROUTINE S_SPH_NP(F,HANDLE,IPAR,SPAR,STAT)
   USE MKL_DFTI

   INTEGER STAT
   INTEGER IPAR(*)
   REAL SPAR(*)
   REAL F(IPAR(11)+1,*)
   TYPE(DFTI_DESCRIPTOR), POINTER :: HANDLE
END SUBROUTINE


END INTERFACE
END MODULE
