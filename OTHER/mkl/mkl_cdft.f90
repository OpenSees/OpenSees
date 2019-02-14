!===============================================================================
! Copyright 2002-2018 Intel Corporation All Rights Reserved.
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
!      Intel(R) Math Kernel Library (Intel(R) MKL) interface for Cluster DFT routines
!*******************************************************************************

! Include to build module MKL_DFTI
INCLUDE 'mkl_dfti.f90'

! Definition of module MKL_CDFT_DM_TYPE. It is used just to define type DFTI_DESCRIPTOR_DM
MODULE MKL_CDFT_DM_TYPE

! Definition of descriptor.
! Structure of this type is not used in Fortran code. The pointer to this type is used only
TYPE DFTI_DESCRIPTOR_DM
    PRIVATE
    INTEGER(4) DESCRIPTOR
END TYPE

END MODULE

! Definition of module MKL_CDFT. It is used to define constants and interfaces of routines
MODULE MKL_CDFT

! Module MKL_CDFT includes definitions from module MKL_DFTI and MKL_CDFT_DM_TYPE
USE MKL_DFTI
USE MKL_CDFT_DM_TYPE

IMPLICIT NONE

! Codes of parameters for DftiGetValueDM / DftiSetValueDM
INTEGER, PARAMETER :: CDFT_LOCAL_SIZE        =1000
INTEGER, PARAMETER :: CDFT_LOCAL_X_START     =1001
INTEGER, PARAMETER :: CDFT_LOCAL_NX          =1002
INTEGER, PARAMETER :: CDFT_MPI_COMM          =1003
INTEGER, PARAMETER :: CDFT_WORKSPACE         =1004
INTEGER, PARAMETER :: CDFT_LOCAL_OUT_X_START =1005
INTEGER, PARAMETER :: CDFT_LOCAL_OUT_NX      =1006

! Codes of errors
INTEGER, PARAMETER :: CDFT_MPI_ERROR     =1000
INTEGER, PARAMETER :: CDFT_SPREAD_ERROR  =1001

! Interfaces of routines
INTERFACE DftiCreateDescriptorDM
    MODULE PROCEDURE DftiCreateDescriptorDM1
    MODULE PROCEDURE DftiCreateDescriptorDMn
    MODULE PROCEDURE DftiCreateDescriptorDM1_s
    MODULE PROCEDURE DftiCreateDescriptorDM1_d
    MODULE PROCEDURE DftiCreateDescriptorDMn_s
    MODULE PROCEDURE DftiCreateDescriptorDMn_d
END INTERFACE
PRIVATE DftiCreateDescriptorDM1
PRIVATE DftiCreateDescriptorDMn
PRIVATE DftiCreateDescriptorDM1_s
PRIVATE DftiCreateDescriptorDM1_d
PRIVATE DftiCreateDescriptorDMn_s
PRIVATE DftiCreateDescriptorDMn_d

INTERFACE DftiGetValueDM
    MODULE PROCEDURE DftiGetValueDMs
    MODULE PROCEDURE DftiGetValueDMd
    MODULE PROCEDURE DftiGetValueDMi
END INTERFACE
PRIVATE DftiGetValueDMs
PRIVATE DftiGetValueDMd
PRIVATE DftiGetValueDMi

INTERFACE DftiSetValueDM
    MODULE PROCEDURE DftiSetValueDMs
    MODULE PROCEDURE DftiSetValueDMd
    MODULE PROCEDURE DftiSetValueDMi
    MODULE PROCEDURE DftiSetValueDMpc
    MODULE PROCEDURE DftiSetValueDMpz
    MODULE PROCEDURE DftiSetValueDMps
    MODULE PROCEDURE DftiSetValueDMpd
    MODULE PROCEDURE DftiSetValueDMpi
END INTERFACE
PRIVATE DftiSetValueDMs
PRIVATE DftiSetValueDMd
PRIVATE DftiSetValueDMi
PRIVATE DftiSetValueDMpc
PRIVATE DftiSetValueDMpz
PRIVATE DftiSetValueDMps
PRIVATE DftiSetValueDMpd
PRIVATE DftiSetValueDMpi

INTERFACE DftiCommitDescriptorDM
    INTEGER FUNCTION DftiCommitDescriptorDM_internal(H)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    END FUNCTION
END INTERFACE
PRIVATE DftiCommitDescriptorDM_internal

INTERFACE DftiComputeForwardDM
    MODULE PROCEDURE DftiComputeForwardDMos
    MODULE PROCEDURE DftiComputeForwardDMod
    MODULE PROCEDURE DftiComputeForwardDMoc
    MODULE PROCEDURE DftiComputeForwardDMoz
    MODULE PROCEDURE DftiComputeForwardDMis
    MODULE PROCEDURE DftiComputeForwardDMid
    MODULE PROCEDURE DftiComputeForwardDMic
    MODULE PROCEDURE DftiComputeForwardDMiz
END INTERFACE
PRIVATE DftiComputeForwardDMos
PRIVATE DftiComputeForwardDMod
PRIVATE DftiComputeForwardDMoc
PRIVATE DftiComputeForwardDMoz
PRIVATE DftiComputeForwardDMis
PRIVATE DftiComputeForwardDMid
PRIVATE DftiComputeForwardDMic
PRIVATE DftiComputeForwardDMiz

INTERFACE DftiComputeBackwardDM
    MODULE PROCEDURE DftiComputeBackwardDMos
    MODULE PROCEDURE DftiComputeBackwardDMod
    MODULE PROCEDURE DftiComputeBackwardDMoc
    MODULE PROCEDURE DftiComputeBackwardDMoz
    MODULE PROCEDURE DftiComputeBackwardDMis
    MODULE PROCEDURE DftiComputeBackwardDMid
    MODULE PROCEDURE DftiComputeBackwardDMic
    MODULE PROCEDURE DftiComputeBackwardDMiz
END INTERFACE
PRIVATE DftiComputeBackwardDMos
PRIVATE DftiComputeBackwardDMod
PRIVATE DftiComputeBackwardDMoc
PRIVATE DftiComputeBackwardDMoz
PRIVATE DftiComputeBackwardDMis
PRIVATE DftiComputeBackwardDMid
PRIVATE DftiComputeBackwardDMic
PRIVATE DftiComputeBackwardDMiz

INTERFACE DftiFreeDescriptorDM
    INTEGER FUNCTION DftiFreeDescriptorDM_internal(H)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
    END FUNCTION
END INTERFACE
PRIVATE DftiFreeDescriptorDM_internal

CONTAINS

!INTERFACE DftiCreateDescriptorDM

! overloading of DftiCreateDescriptorDM for nD DFT
INTEGER FUNCTION DftiCreateDescriptorDMn(C,H,P1,P2,D,L)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        INTEGER(4) C
        INTEGER P1,P2,D,L(*)
        INTENT(IN) :: C,P1,P2,D,L
        INTERFACE
          INTEGER FUNCTION DftiCreateDescriptorDMn_internal(C,H,P1,P2,D,L)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            INTEGER(4) C
            INTEGER P1,P2,D,L(*)
            INTENT(IN) :: C,P1,P2,D,L
          END FUNCTION
        END INTERFACE
        DftiCreateDescriptorDMn = DftiCreateDescriptorDMn_internal(C,H,P1,P2,D,L)
END FUNCTION

! overloading of DftiCreateDescriptorDM for 1D DFT
INTEGER FUNCTION DftiCreateDescriptorDM1(C,H,P1,P2,D,L)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        INTEGER(4) C
        INTEGER P1,P2,D,L
        INTENT(IN) :: C,P1,P2,D,L
        INTERFACE
          INTEGER FUNCTION DftiCreateDescriptorDM1_internal(C,H,P1,P2,D,L)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            INTEGER(4) C
            INTEGER P1,P2,D,L
            INTENT(IN) :: C,P1,P2,D,L
          END FUNCTION
        END INTERFACE
        DftiCreateDescriptorDM1 = DftiCreateDescriptorDM1_internal(C,H,P1,P2,D,L)
END FUNCTION
INTEGER FUNCTION DftiCreateDescriptorDM1_s(C,H,P1R,P2,D,L)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        REAL(4), INTENT(IN) :: P1R
        INTEGER(4) C
        INTEGER P1,P2,D,L
        INTENT(IN) :: C,P2,D,L
        P1 = INT(P1R)
        DftiCreateDescriptorDM1_s = DftiCreateDescriptorDM1(C,H,P1,P2,D,L)
END FUNCTION
INTEGER FUNCTION DftiCreateDescriptorDM1_d(C,H,P1R,P2,D,L)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        REAL(8), INTENT(IN) :: P1R
        INTEGER(4) C
        INTEGER P1,P2,D,L
        INTENT(IN) :: C,P2,D,L
        P1 = INT(P1R)
        DftiCreateDescriptorDM1_d = DftiCreateDescriptorDM1(C,H,P1,P2,D,L)
END FUNCTION
INTEGER FUNCTION DftiCreateDescriptorDMn_s(C,H,P1R,P2,D,L)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        REAL(4), INTENT(IN) :: P1R
        INTEGER(4) C
        INTEGER P1,P2,D,L(*)
        INTENT(IN) :: C,P2,D,L
        P1 = INT(P1R)
        DftiCreateDescriptorDMn_s = DftiCreateDescriptorDMn(C,H,P1,P2,D,L)
END FUNCTION
INTEGER FUNCTION DftiCreateDescriptorDMn_d(C,H,P1R,P2,D,L)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        REAL(8), INTENT(IN) :: P1R
        INTEGER(4) C
        INTEGER P1,P2,D,L(*)
        INTENT(IN) :: C,P2,D,L
        P1 = INT(P1R)
        DftiCreateDescriptorDMn_d = DftiCreateDescriptorDMn(C,H,P1,P2,D,L)
END FUNCTION
!END INTERFACE DftiCreateDescriptorDM
!INTERFACE DftiSetValueDM

! overloading of DftiSetValueDM for SP value
INTEGER FUNCTION DftiSetValueDMs(H,P,V)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        INTEGER P
        REAL(4) V
        INTENT(IN) :: P,V
        INTERFACE
          INTEGER FUNCTION DftiSetValueDMf_internal(H,P,V)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            INTEGER P
            REAL(4) V
            INTENT(IN) :: P,V
          END FUNCTION
        END INTERFACE
        DftiSetValueDMs = DftiSetValueDMf_internal(H,P,V)
END FUNCTION

! overloading of DftiSetValueDM for DP value
INTEGER FUNCTION DftiSetValueDMd(H,P,V)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        INTEGER P
        REAL(8) V
        INTENT(IN) :: P,V
        INTERFACE
          INTEGER FUNCTION DftiSetValueDMd_internal(H,P,V)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            INTEGER P
            REAL(8) V
            INTENT(IN) :: P,V
          END FUNCTION
        END INTERFACE
        DftiSetValueDMd = DftiSetValueDMd_internal(H,P,V)
END FUNCTION

! overloading of DftiSetValueDM for integer value
INTEGER FUNCTION DftiSetValueDMi(H,P,V)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        INTEGER P,V
        INTENT(IN) :: P,V
        INTERFACE
          INTEGER FUNCTION DftiSetValueDMi_internal(H,P,V)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            INTEGER P,V
            INTENT(IN) :: P,V
          END FUNCTION
        END INTERFACE
        DftiSetValueDMi = DftiSetValueDMi_internal(H,P,V)
END FUNCTION

! overloading of DftiSetValueDM for SP complex array
INTEGER FUNCTION DftiSetValueDMpc(H,P,V)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        INTEGER P
        COMPLEX(4) V(*)
        INTENT(IN) :: P,V
        INTERFACE
          INTEGER FUNCTION DftiSetValueDMp_internal(H,P,V)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            INTEGER P
            COMPLEX(4) V(*)
            INTENT(IN) :: P,V
          END FUNCTION
        END INTERFACE
        DftiSetValueDMpc = DftiSetValueDMp_internal(H,P,V)
END FUNCTION

! overloading of DftiSetValueDM for DP complex array
INTEGER FUNCTION DftiSetValueDMpz(H,P,V)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        INTEGER P
        COMPLEX(8) V(*)
        INTENT(IN) :: P,V
        INTERFACE
          INTEGER FUNCTION DftiSetValueDMp_internal(H,P,V)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            INTEGER P
            COMPLEX(8) V(*)
            INTENT(IN) :: P,V
          END FUNCTION
        END INTERFACE
        DftiSetValueDMpz = DftiSetValueDMp_internal(H,P,V)
END FUNCTION

! overloading of DftiSetValueDM for SP array
INTEGER FUNCTION DftiSetValueDMps(H,P,V)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        INTEGER P
        REAL(4) V(*)
        INTENT(IN) :: P,V
        INTERFACE
          INTEGER FUNCTION DftiSetValueDMp_internal(H,P,V)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            INTEGER P
            REAL(4) V(*)
            INTENT(IN) :: P,V
          END FUNCTION
        END INTERFACE
        DftiSetValueDMps = DftiSetValueDMp_internal(H,P,V)
END FUNCTION

! overloading of DftiSetValueDM for DP array
INTEGER FUNCTION DftiSetValueDMpd(H,P,V)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        INTEGER P
        REAL(8) V(*)
        INTENT(IN) :: P,V
        INTERFACE
          INTEGER FUNCTION DftiSetValueDMp_internal(H,P,V)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            INTEGER P
            REAL(8) V(*)
            INTENT(IN) :: P,V
          END FUNCTION
        END INTERFACE
        DftiSetValueDMpd = DftiSetValueDMp_internal(H,P,V)
END FUNCTION

! overloading of DftiSetValueDM for integer array
INTEGER FUNCTION DftiSetValueDMpi(H,P,V)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        INTEGER P,V(*)
        INTENT(IN) :: P,V
        INTERFACE
          INTEGER FUNCTION DftiSetValueDMp_internal(H,P,V)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            INTEGER P,V(*)
            INTENT(IN) :: P,V
          END FUNCTION
        END INTERFACE
        DftiSetValueDMpi = DftiSetValueDMp_internal(H,P,V)
END FUNCTION
!END INTERFACE DftiSetValueDM
!INTERFACE DftiGetValueDM

! overloading of DftiGetValueDM for SP value
INTEGER FUNCTION DftiGetValueDMs(H,P,V)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        INTEGER P
        REAL(4) V
        INTENT(IN)  :: P
        INTENT(OUT) :: V
        INTERFACE
          INTEGER FUNCTION DftiGetValueDM_internal(H,P,V)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            INTEGER P
            REAL(4) V
            INTENT(IN)  :: P
            INTENT(OUT) :: V
          END FUNCTION
        END INTERFACE
        DftiGetValueDMs = DftiGetValueDM_internal(H,P,V)
END FUNCTION

! overloading of DftiGetValueDM for DP value
INTEGER FUNCTION DftiGetValueDMd(H,P,V)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        INTEGER P
        REAL(8) V
        INTENT(IN)  :: P
        INTENT(OUT) :: V
        INTERFACE
          INTEGER FUNCTION DftiGetValueDM_internal(H,P,V)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            INTEGER P
            REAL(8) V
            INTENT(IN)  :: P
            INTENT(OUT) :: V
          END FUNCTION
        END INTERFACE
        DftiGetValueDMd = DftiGetValueDM_internal(H,P,V)
END FUNCTION

! overloading of DftiGetValueDM for integer value
INTEGER FUNCTION DftiGetValueDMi(H,P,V)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        INTEGER P,V
        INTENT(IN)  :: P
        INTENT(OUT) :: V
        INTERFACE
          INTEGER FUNCTION DftiGetValueDM_internal(H,P,V)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            INTEGER P,V
            INTENT(IN)  :: P
            INTENT(OUT) :: V
          END FUNCTION
        END INTERFACE
        DftiGetValueDMi = DftiGetValueDM_internal(H,P,V)
END FUNCTION
!END INTERFACE DftiGetValueDM
!INTERFACE DftiComputeForwardDM

    ! overloading of DftiComputeForwardDM for SP R2C DFT (out-of-place)
    INTEGER FUNCTION DftiComputeForwardDMos(H,IN,OUT)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        REAL(4) IN(*)
        COMPLEX(4) OUT(*)
        INTENT(IN)  :: IN
        INTENT(OUT) :: OUT
        INTERFACE
          INTEGER FUNCTION DftiComputeForwardDMo_internal(H,IN,OUT)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            REAL(4) IN(*)
            COMPLEX(4) OUT(*)
            INTENT(IN)  :: IN
            INTENT(OUT) :: OUT
          END FUNCTION
        END INTERFACE
        DftiComputeForwardDMos = DftiComputeForwardDMo_internal(H,IN,OUT)
    END FUNCTION

    ! overloading of DftiComputeForwardDM for DP R2C DFT (out-of-place)
    INTEGER FUNCTION DftiComputeForwardDMod(H,IN,OUT)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        REAL(8) IN(*)
        COMPLEX(8) OUT(*)
        INTENT(IN)  :: IN
        INTENT(OUT) :: OUT
        INTERFACE
          INTEGER FUNCTION DftiComputeForwardDMo_internal(H,IN,OUT)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            REAL(8) IN(*)
            COMPLEX(8) OUT(*)
            INTENT(IN)  :: IN
            INTENT(OUT) :: OUT
          END FUNCTION
        END INTERFACE
        DftiComputeForwardDMod = DftiComputeForwardDMo_internal(H,IN,OUT)
    END FUNCTION

    ! overloading of DftiComputeForwardDM for SP C2C DFT (out-of-place)
    INTEGER FUNCTION DftiComputeForwardDMoc(H,IN,OUT)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        COMPLEX(4) IN(*),OUT(*)
        INTENT(IN)  :: IN
        INTENT(OUT) :: OUT
        INTERFACE
          INTEGER FUNCTION DftiComputeForwardDMo_internal(H,IN,OUT)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            COMPLEX(4) IN(*),OUT(*)
            INTENT(IN)  :: IN
            INTENT(OUT) :: OUT
          END FUNCTION
        END INTERFACE
        DftiComputeForwardDMoc = DftiComputeForwardDMo_internal(H,IN,OUT)
    END FUNCTION

    ! overloading of DftiComputeForwardDM for DP C2C DFT (out-of-place)
    INTEGER FUNCTION DftiComputeForwardDMoz(H,IN,OUT)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        COMPLEX(8) IN(*),OUT(*)
        INTENT(IN)  :: IN
        INTENT(OUT) :: OUT
        INTERFACE
          INTEGER FUNCTION DftiComputeForwardDMo_internal(H,IN,OUT)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            COMPLEX(8) IN(*),OUT(*)
            INTENT(IN)  :: IN
            INTENT(OUT) :: OUT
          END FUNCTION
        END INTERFACE
        DftiComputeForwardDMoz = DftiComputeForwardDMo_internal(H,IN,OUT)
    END FUNCTION

    ! overloading of DftiComputeForwardDM for SP R2C DFT (inplace)
    INTEGER FUNCTION DftiComputeForwardDMis(H,IN)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        REAL(4) IN(*)
        INTERFACE
          INTEGER FUNCTION DftiComputeForwardDMi_internal(H,IN)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            REAL(4) IN(*)
          END FUNCTION
        END INTERFACE
        DftiComputeForwardDMis = DftiComputeForwardDMi_internal(H,IN)
    END FUNCTION

    ! overloading of DftiComputeForwardDM for DP R2C DFT (inplace)
    INTEGER FUNCTION DftiComputeForwardDMid(H,IN)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        REAL(8) IN(*)
        INTERFACE
          INTEGER FUNCTION DftiComputeForwardDMi_internal(H,IN)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            REAL(8) IN(*)
          END FUNCTION
        END INTERFACE
        DftiComputeForwardDMid = DftiComputeForwardDMi_internal(H,IN)
    END FUNCTION

    ! overloading of DftiComputeForwardDM for SP C2C DFT (inplace)
    INTEGER FUNCTION DftiComputeForwardDMic(H,IN)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        COMPLEX(4) IN(*)
        INTERFACE
          INTEGER FUNCTION DftiComputeForwardDMi_internal(H,IN)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            COMPLEX(4) IN(*)
          END FUNCTION
        END INTERFACE
        DftiComputeForwardDMic = DftiComputeForwardDMi_internal(H,IN)
    END FUNCTION

    ! overloading of DftiComputeForwardDM for DP C2C DFT (inplace)
    INTEGER FUNCTION DftiComputeForwardDMiz(H,IN)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        COMPLEX(8) IN(*)
        INTERFACE
          INTEGER FUNCTION DftiComputeForwardDMi_internal(H,IN)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            COMPLEX(8) IN(*)
          END FUNCTION
        END INTERFACE
        DftiComputeForwardDMiz = DftiComputeForwardDMi_internal(H,IN)
    END FUNCTION
!END INTERFACE DftiComputeForwardDM
!INTERFACE DftiComputeBackwardDM

! overloading of DftiComputeBackwardDM for SP R2C DFT (out-of-place)
INTEGER FUNCTION DftiComputeBackwardDMos(H,IN,OUT)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        COMPLEX(4) IN(*)
        REAL(4) OUT(*)
        INTENT(IN)  :: IN
        INTENT(OUT) :: OUT
        INTERFACE
          INTEGER FUNCTION DftiComputeBackwardDMo_internal(H,IN,OUT)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            COMPLEX(4) IN(*)
            REAL(4) OUT(*)
            INTENT(IN)  :: IN
            INTENT(OUT) :: OUT
          END FUNCTION
        END INTERFACE
        DftiComputeBackwardDMos = DftiComputeBackwardDMo_internal(H,IN,OUT)
END FUNCTION

! overloading of DftiComputeBackwardDM for DP R2C DFT (out-of-place)
INTEGER FUNCTION DftiComputeBackwardDMod(H,IN,OUT)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        COMPLEX(8) IN(*)
        REAL(8) OUT(*)
        INTENT(IN)  :: IN
        INTENT(OUT) :: OUT
        INTERFACE
          INTEGER FUNCTION DftiComputeBackwardDMo_internal(H,IN,OUT)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            COMPLEX(8) IN(*)
            REAL(8) OUT(*)
            INTENT(IN)  :: IN
            INTENT(OUT) :: OUT
          END FUNCTION
        END INTERFACE
        DftiComputeBackwardDMod = DftiComputeBackwardDMo_internal(H,IN,OUT)
END FUNCTION

! overloading of DftiComputeBackwardDM for SP C2C DFT (out-of-place)
INTEGER FUNCTION DftiComputeBackwardDMoc(H,IN,OUT)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        COMPLEX(4) IN(*),OUT(*)
        INTENT(IN)  :: IN
        INTENT(OUT) :: OUT
        INTERFACE
          INTEGER FUNCTION DftiComputeBackwardDMo_internal(H,IN,OUT)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            COMPLEX(4) IN(*),OUT(*)
            INTENT(IN)  :: IN
            INTENT(OUT) :: OUT
          END FUNCTION
        END INTERFACE
        DftiComputeBackwardDMoc = DftiComputeBackwardDMo_internal(H,IN,OUT)
END FUNCTION

! overloading of DftiComputeBackwardDM for DP C2C DFT (out-of-place)
INTEGER FUNCTION DftiComputeBackwardDMoz(H,IN,OUT)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        COMPLEX(8) IN(*),OUT(*)
        INTENT(IN)  :: IN
        INTENT(OUT) :: OUT
        INTERFACE
          INTEGER FUNCTION DftiComputeBackwardDMo_internal(H,IN,OUT)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            COMPLEX(8) IN(*),OUT(*)
            INTENT(IN)  :: IN
            INTENT(OUT) :: OUT
          END FUNCTION
        END INTERFACE
        DftiComputeBackwardDMoz = DftiComputeBackwardDMo_internal(H,IN,OUT)
END FUNCTION

! overloading of DftiComputeBackwardDM for SP R2C DFT (inplace)
INTEGER FUNCTION DftiComputeBackwardDMis(H,IN)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        REAL(4) IN(*)
        INTERFACE
          INTEGER FUNCTION DftiComputeBackwardDMi_internal(H,IN)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            REAL(4) IN(*)
          END FUNCTION
        END INTERFACE
        DftiComputeBackwardDMis = DftiComputeBackwardDMi_internal(H,IN)
END FUNCTION

! overloading of DftiComputeBackwardDM for DP R2C DFT (inplace)
INTEGER FUNCTION DftiComputeBackwardDMid(H,IN)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        REAL(8) IN(*)
        INTERFACE
          INTEGER FUNCTION DftiComputeBackwardDMi_internal(H,IN)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            REAL(8) IN(*)
          END FUNCTION
        END INTERFACE
        DftiComputeBackwardDMid = DftiComputeBackwardDMi_internal(H,IN)
END FUNCTION

! overloading of DftiComputeBackwardDM for SP C2C DFT (inplace)
INTEGER FUNCTION DftiComputeBackwardDMic(H,IN)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        COMPLEX(4) IN(*)
        INTERFACE
          INTEGER FUNCTION DftiComputeBackwardDMi_internal(H,IN)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            COMPLEX(4) IN(*)
          END FUNCTION
        END INTERFACE
        DftiComputeBackwardDMic = DftiComputeBackwardDMi_internal(H,IN)
END FUNCTION

! overloading of DftiComputeBackwardDM for DP C2C DFT (inplace)
INTEGER FUNCTION DftiComputeBackwardDMiz(H,IN)
        USE MKL_CDFT_DM_TYPE
        TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
        COMPLEX(8) IN(*)
        INTERFACE
          INTEGER FUNCTION DftiComputeBackwardDMi_internal(H,IN)
            USE MKL_CDFT_DM_TYPE
            TYPE(DFTI_DESCRIPTOR_DM), POINTER :: H
            COMPLEX(8) IN(*)
          END FUNCTION
        END INTERFACE
        DftiComputeBackwardDMiz = DftiComputeBackwardDMi_internal(H,IN)
END FUNCTION
!END INTERFACE DftiComputeBackwardDM
END MODULE
