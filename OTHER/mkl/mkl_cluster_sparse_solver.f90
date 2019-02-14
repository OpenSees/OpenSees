!===============================================================================
! Copyright 1999-2018 Intel Corporation All Rights Reserved.
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
!      Intel(R) Math Kernel Library (Intel(R) MKL) Fortran 90 interface for
!	   Cluster Sparse Solver
!*******************************************************************************

!DEC$ IF .NOT. DEFINED( __MKL_CLUSTER_SPARSE_SOLVER_F90 )

!DEC$ DEFINE __MKL_CLUSTER_SPARSE_SOLVER_F90

      MODULE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
        TYPE MKL_CLUSTER_SPARSE_SOLVER_HANDLE; INTEGER(KIND=8) DUMMY; END TYPE
      END MODULE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE

      MODULE MKL_CLUSTER_SPARSE_SOLVER
        USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE

!
! Subroutine prototype for CLUSTER_SPARSE_SOLVER
!

      INTERFACE CLUSTER_SPARSE_SOLVER
        SUBROUTINE CLUSTER_SPARSE_SOLVER_D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(IN)    :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          REAL(KIND=8),     INTENT(IN)    :: A(*)
          REAL(KIND=8),     INTENT(INOUT) :: B(*)
          REAL(KIND=8),     INTENT(OUT)   :: X(*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_D

        SUBROUTINE CLUSTER_SPARSE_SOLVER_S(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(IN)    :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          REAL(KIND=4),     INTENT(IN)    :: A(*)
          REAL(KIND=4),     INTENT(INOUT) :: B(*)
          REAL(KIND=4),     INTENT(OUT)   :: X(*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_S

        SUBROUTINE CLUSTER_SPARSE_SOLVER_DC(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(IN)    :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=8),     INTENT(IN)    :: A(*)
          COMPLEX(KIND=8),     INTENT(INOUT) :: B(*)
          COMPLEX(KIND=8),     INTENT(OUT)   :: X(*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_DC

        SUBROUTINE CLUSTER_SPARSE_SOLVER_SC(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(IN)    :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=4),     INTENT(IN)    :: A(*)
          COMPLEX(KIND=4),     INTENT(INOUT) :: B(*)
          COMPLEX(KIND=4),     INTENT(OUT)   :: X(*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_SC

        SUBROUTINE CLUSTER_SPARSE_SOLVER_D_2D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(IN)    :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          REAL(KIND=8),     INTENT(IN)    :: A(*)
          REAL(KIND=8),     INTENT(INOUT) :: B(N,*)
          REAL(KIND=8),     INTENT(OUT)   :: X(N,*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_D_2D

        SUBROUTINE CLUSTER_SPARSE_SOLVER_S_2D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(IN)    :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          REAL(KIND=4),     INTENT(IN)    :: A(*)
          REAL(KIND=4),     INTENT(INOUT) :: B(N,*)
          REAL(KIND=4),     INTENT(OUT)   :: X(N,*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_S_2D

        SUBROUTINE CLUSTER_SPARSE_SOLVER_DC_2D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(IN)    :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=8),     INTENT(IN)    :: A(*)
          COMPLEX(KIND=8),     INTENT(INOUT) :: B(N,*)
          COMPLEX(KIND=8),     INTENT(OUT)   :: X(N,*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_DC_2D

        SUBROUTINE CLUSTER_SPARSE_SOLVER_SC_2D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER,          INTENT(IN)    :: MAXFCT
          INTEGER,          INTENT(IN)    :: MNUM
          INTEGER,          INTENT(IN)    :: MTYPE
          INTEGER,          INTENT(IN)    :: PHASE
          INTEGER,          INTENT(IN)    :: N
          INTEGER,          INTENT(IN)    :: IA(*)
          INTEGER,          INTENT(IN)    :: JA(*)
          INTEGER,          INTENT(IN)    :: PERM(*)
          INTEGER,          INTENT(IN)    :: NRHS
          INTEGER,          INTENT(INOUT) :: IPARM(*)
          INTEGER,          INTENT(IN)    :: MSGLVL
          INTEGER,          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=4),     INTENT(IN)    :: A(*)
          COMPLEX(KIND=4),     INTENT(INOUT) :: B(N,*)
          COMPLEX(KIND=4),     INTENT(OUT)   :: X(N,*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_SC_2D
      END INTERFACE

!
! Subroutine prototype for CLUSTER_SPARSE_SOLVER_64
!

      INTERFACE CLUSTER_SPARSE_SOLVER_64
        SUBROUTINE CLUSTER_SPARSE_SOLVER_D_64(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          REAL(KIND=8),     INTENT(IN)    :: A(*)
          REAL(KIND=8),     INTENT(INOUT) :: B(*)
          REAL(KIND=8),     INTENT(OUT)   :: X(*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_D_64

        SUBROUTINE CLUSTER_SPARSE_SOLVER_S_64(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          REAL(KIND=4),     INTENT(IN)    :: A(*)
          REAL(KIND=4),     INTENT(INOUT) :: B(*)
          REAL(KIND=4),     INTENT(OUT)   :: X(*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_S_64

        SUBROUTINE CLUSTER_SPARSE_SOLVER_DC_64(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=8),     INTENT(IN)    :: A(*)
          COMPLEX(KIND=8),     INTENT(INOUT) :: B(*)
          COMPLEX(KIND=8),     INTENT(OUT)   :: X(*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_DC_64

        SUBROUTINE CLUSTER_SPARSE_SOLVER_SC_64(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=4),     INTENT(IN)    :: A(*)
          COMPLEX(KIND=4),     INTENT(INOUT) :: B(*)
          COMPLEX(KIND=4),     INTENT(OUT)   :: X(*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_SC_64

        SUBROUTINE CLUSTER_SPARSE_SOLVER_D_64_2D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          REAL(KIND=8),     INTENT(IN)    :: A(*)
          REAL(KIND=8),     INTENT(INOUT) :: B(N,*)
          REAL(KIND=8),     INTENT(OUT)   :: X(N,*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_D_64_2D

        SUBROUTINE CLUSTER_SPARSE_SOLVER_S_64_2D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          REAL(KIND=4),     INTENT(IN)    :: A(*)
          REAL(KIND=4),     INTENT(INOUT) :: B(N,*)
          REAL(KIND=4),     INTENT(OUT)   :: X(N,*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_S_64_2D

        SUBROUTINE CLUSTER_SPARSE_SOLVER_DC_64_2D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=8),     INTENT(IN)    :: A(*)
          COMPLEX(KIND=8),     INTENT(INOUT) :: B(N,*)
          COMPLEX(KIND=8),     INTENT(OUT)   :: X(N,*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_DC_64_2D

        SUBROUTINE CLUSTER_SPARSE_SOLVER_SC_64_2D(PT,MAXFCT,MNUM,MTYPE,PHASE,N,A,IA,JA,PERM,NRHS,IPARM,MSGLVL,B,X,COMM,ERROR)
          USE MKL_CLUSTER_SPARSE_SOLVER_PRIVATE
          TYPE(MKL_CLUSTER_SPARSE_SOLVER_HANDLE), INTENT(INOUT) :: PT(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MAXFCT
          INTEGER(KIND=8),          INTENT(IN)    :: MNUM
          INTEGER(KIND=8),          INTENT(IN)    :: MTYPE
          INTEGER(KIND=8),          INTENT(IN)    :: PHASE
          INTEGER(KIND=8),          INTENT(IN)    :: N
          INTEGER(KIND=8),          INTENT(IN)    :: IA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: JA(*)
          INTEGER(KIND=8),          INTENT(IN)    :: PERM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: NRHS
          INTEGER(KIND=8),          INTENT(INOUT) :: IPARM(*)
          INTEGER(KIND=8),          INTENT(IN)    :: MSGLVL
          INTEGER(KIND=8),          INTENT(OUT)   :: ERROR
          COMPLEX(KIND=4),     INTENT(IN)    :: A(*)
          COMPLEX(KIND=4),     INTENT(INOUT) :: B(N,*)
          COMPLEX(KIND=4),     INTENT(OUT)   :: X(N,*)
          INTEGER*4,          INTENT(IN)    :: COMM
        END SUBROUTINE CLUSTER_SPARSE_SOLVER_SC_64_2D
      END INTERFACE

      END MODULE MKL_CLUSTER_SPARSE_SOLVER

!DEC$ ENDIF
