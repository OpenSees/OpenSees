C
C  This file is part of MUMPS 5.4.1, released
C  on Tue Aug  3 09:49:43 UTC 2021
C
C
C  Copyright 1991-2021 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
C  Mumps Technologies, University of Bordeaux.
C
C  This version of MUMPS is provided to you free of charge. It is
C  released under the CeCILL-C license 
C  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
C  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
C
      INTEGER FUNCTION MUMPS_BLOC2_GET_NSLAVESMIN 
     &         ( SLAVEF, K48, K821, K50,
     &         NFRONT, NCB, K375, K119)
      IMPLICIT NONE
      INTEGER,    INTENT (IN) :: SLAVEF, K48, K50, NFRONT, NCB
      INTEGER,    INTENT (IN) :: K375
      INTEGER,    INTENT (IN) :: K119
      INTEGER(8), INTENT (IN) :: K821
      INTEGER NSLAVESMIN, NASS, KMAX
      REAL Wmaster, Wtotal, Wmax
      INTEGER ACC,X
      REAL MUMPS_BLOC2_COUT
      INTEGER MUMPS_REG_GETKMAX    
      EXTERNAL  MUMPS_BLOC2_COUT, MUMPS_REG_GETKMAX
      KMAX = MUMPS_REG_GETKMAX( K821, NCB )
      NASS = NFRONT - NCB
      NSLAVESMIN  = 1
      IF ( K48 .EQ.0 .OR.
     &    (K48.EQ.5 .AND.  (K119.EQ.1. OR.K50.EQ.0))) THEN 
         NSLAVESMIN = max(NCB/max(1,KMAX),1)
      ELSE IF (K48 .EQ. 3 .OR.(K48.EQ.5 .AND.K50.NE.0) ) THEN
        Wmax    = MUMPS_BLOC2_COUT(KMAX,NFRONT,NASS)
        Wtotal  = MUMPS_BLOC2_COUT(NCB,NFRONT,NASS)
        Wmaster = real(NASS)*real(NASS)*real(NASS)/(3.0E0)
        IF ( Wmaster .GT.  Wmax ) THEN 
         NSLAVESMIN = max ( nint ( Wtotal / Wmaster ), 1 )
        ELSE
         NSLAVESMIN = max ( nint ( Wtotal / Wmax ), 1 )
        ENDIF
        IF (K48 .EQ. 5) THEN
          IF (K119.EQ.2) THEN
            NSLAVESMIN = max ( NSLAVESMIN/2, 1 )
          ENDIF
        END IF
      ELSE IF (K48 .EQ. 4 ) THEN
         IF ( K821 > 0_8 ) THEN
           WRITE(*,*) 'Internal Error 1 in MUMPS_BLOC2_GET_NSLAVESMIN'
           CALL MUMPS_ABORT()
         ENDIF
         CALL MUMPS_ABORT_ON_OVERFLOW(K821,
     &           "K821 too large in MUMPS_BLOC2_GET_NSLAVESMIN" )
         KMAX=int(abs(K821))
         IF(K50.EQ.0)THEN
            NSLAVESMIN = max(int(
     &                (int(NCB,8)*int(NCB,8))/int(KMAX,8)
     &                  ),1)
         ELSE
            ACC=0
            NSLAVESMIN=0
            DO WHILE (ACC.NE.NCB)
               X=int((-real(NFRONT-NCB+ACC)
     &              +sqrt(((real(NFRONT-NCB+ACC)*
     &              real(NFRONT-NCB+ACC))+real(4)*
     &              real(KMAX))))/
     &              real(2))
               ACC=ACC+X
               NSLAVESMIN=NSLAVESMIN+1
               IF (((NCB-ACC)*NCB).LT.KMAX)THEN
                  ACC=NCB 
                  NSLAVESMIN=NSLAVESMIN+1
               ENDIF
            ENDDO
         ENDIF
      ENDIF
      NSLAVESMIN = min ( NSLAVESMIN,(SLAVEF-1) )
      MUMPS_BLOC2_GET_NSLAVESMIN = 
     &               min ( NSLAVESMIN, NCB )
      IF (K375 .EQ. 1) THEN
        MUMPS_BLOC2_GET_NSLAVESMIN=1
      ENDIF
      RETURN
      END FUNCTION MUMPS_BLOC2_GET_NSLAVESMIN 
      INTEGER FUNCTION MUMPS_BLOC2_GET_NSLAVESMAX 
     &        ( SLAVEF, K48, K821, K50,
     &          NFRONT, NCB, K375, K119 )
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: SLAVEF, K48, K50,NFRONT, NCB, K375, K119
      INTEGER(8), INTENT(IN) :: K821
      INTEGER NSLAVESMAX, KMAX, KMIN
      INTEGER NSLAVESMIN
      INTEGER MUMPS_REG_GETKMAX,MUMPS_GETKMIN,
     &        MUMPS_BLOC2_GET_NSLAVESMIN,
     &        MUMPS_BLOC2_GET_NS_BLSIZE
      EXTERNAL MUMPS_REG_GETKMAX,MUMPS_GETKMIN,
     &        MUMPS_BLOC2_GET_NSLAVESMIN,
     &        MUMPS_BLOC2_GET_NS_BLSIZE
      IF (K48 .eq. 0 .OR. K48.eq.3.OR.K48.EQ.5) THEN
         KMAX = MUMPS_REG_GETKMAX( K821, NCB )
         KMIN = MUMPS_GETKMIN( K821, K50, KMAX, NCB)
         NSLAVESMAX = MUMPS_BLOC2_GET_NS_BLSIZE(
     &                SLAVEF, K48, K50, KMIN, NFRONT, NCB )
      ELSE
         NSLAVESMAX = SLAVEF-1
      ENDIF
      NSLAVESMIN = MUMPS_BLOC2_GET_NSLAVESMIN(
     &     SLAVEF, K48, K821, K50, NFRONT, NCB, K375, K119 )
      NSLAVESMAX = max ( NSLAVESMAX, NSLAVESMIN )
      MUMPS_BLOC2_GET_NSLAVESMAX = 
     &               min ( NSLAVESMAX, NCB )
      IF (K375 .EQ. 1) THEN
        MUMPS_BLOC2_GET_NSLAVESMAX = SLAVEF-1 
      ENDIF
      RETURN
      END FUNCTION MUMPS_BLOC2_GET_NSLAVESMAX
      SUBROUTINE MUMPS_MAX_SURFCB_NBROWS( WHAT, KEEP,KEEP8,
     &           NCB, NFR, SLAVEF, NBROWMAX, MAXSURFCB8
     &     )
      IMPLICIT NONE
      INTEGER, intent(in) :: WHAT, NCB, NFR, SLAVEF
      INTEGER, intent(in) :: KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER, intent(out) :: NBROWMAX
      INTEGER(8), intent(out) :: MAXSURFCB8
      INTEGER KMAX, KMIN, NSLAVES, SIZEDUMMY, TABDUMMY(1)
      EXTERNAL MUMPS_REG_GETKMAX, MUMPS_GETKMIN,
     &         MUMPS_BLOC2_GET_NSLAVESMIN
      INTEGER MUMPS_REG_GETKMAX, MUMPS_GETKMIN,
     &        MUMPS_BLOC2_GET_NSLAVESMIN
      IF ( WHAT .NE. 1 .and. WHAT .NE. 2 ) THEN
        IF (WHAT .NE. 4 .and. WHAT .NE. 5 .AND.
     &        KEEP(48).NE.5 ) THEN
        WRITE(*,*) "Internal error 1 in MUMPS_MAX_SURFCB_NBROWS"
        CALL MUMPS_ABORT()
        END IF
      ENDIF
      KMAX    = MUMPS_REG_GETKMAX( KEEP8(21), NCB )
      IF (WHAT .EQ.1.OR.WHAT.EQ.2) THEN
        NSLAVES = MUMPS_BLOC2_GET_NSLAVESMIN( SLAVEF, KEEP(48),
     &            KEEP8(21), KEEP(50),
     &            NFR, NCB, KEEP(375), KEEP(119) )
      ELSE
        NSLAVES=SLAVEF
      ENDIF
      IF ( KEEP(48) == 0 .OR. (KEEP(48).EQ.5.AND.KEEP(50).EQ.0)) THEN
        NBROWMAX = NCB / NSLAVES + mod( NCB, NSLAVES )
        IF ( WHAT == 2 .OR. WHAT == 5 )
     &    MAXSURFCB8 = int(NBROWMAX,8) * int(NCB,8)
      ELSE IF (KEEP(48) == 3.OR.(KEEP(48).EQ.5.AND.KEEP(50).NE.0))THEN
        KMIN = MUMPS_GETKMIN( KEEP8(21), KEEP(50), KMAX, NCB )  
        SIZEDUMMY        = 1
        IF (WHAT.GT.3) THEN
           CALL  MUMPS_BLOC2_SET_POSK483(
     &          WHAT-3, NSLAVES, NFR, NCB, 
     &          KMIN, KMAX, SLAVEF,
     &          NBROWMAX, MAXSURFCB8, TABDUMMY, SIZEDUMMY) 
        ELSE
           CALL  MUMPS_BLOC2_SET_POSK483(
     &          WHAT, NSLAVES, NFR, NCB, 
     &          KMIN, KMAX, SLAVEF,
     &          NBROWMAX, MAXSURFCB8, TABDUMMY, SIZEDUMMY) 
        ENDIF
      ELSE IF ( KEEP(48) == 4 ) THEN
         IF (KEEP8(21) > 0_8) THEN
            WRITE(*,*) "Internal error 2 in MUMPS_MAX_SURFCB_NBROWS"
            CALL MUMPS_ABORT()
         END IF
         IF(KEEP(50).EQ.0)THEN
            IF ( abs(KEEP8(21)) * int( SLAVEF - 1,8 ) >
     &                            int( NCB,8) * int(NFR,8) ) THEN
              NBROWMAX = (NCB + SLAVEF -2 ) / ( SLAVEF - 1 )
              IF ( WHAT == 2 ) MAXSURFCB8 = int(NBROWMAX,8) *int(NCB,8)
            ELSE
              NBROWMAX=int(
     &                      (abs(KEEP8(21)) + int(NFR - 1,8))
     &                    /  int(NFR,8)
     &                    )
              IF ( WHAT == 2 ) MAXSURFCB8 = abs(KEEP8(21))
            ENDIF
         ELSE
            NBROWMAX=int((-real(NFR-NCB)
     &              +sqrt((real(NFR-NCB)*
     &              real(NFR-NCB))+real(4)*
     &              real(abs(KEEP8(21)))))/
     &              real(2))
            IF ( WHAT == 2 ) MAXSURFCB8 = abs(KEEP8(21))
         ENDIF
      ELSE
        NBROWMAX = NCB
        IF (WHAT == 2) MAXSURFCB8 = int(NCB,8) * int(NCB,8)
      ENDIF 
      NBROWMAX = min ( max(NBROWMAX, 1), NCB)
      RETURN
      END SUBROUTINE MUMPS_MAX_SURFCB_NBROWS
      INTEGER FUNCTION MUMPS_BLOC2_GET_NS_BLSIZE( SLAVEF, K48, K50,
     &         BLSIZE, NFRONT, NCB)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: SLAVEF, K48, K50, BLSIZE, NFRONT, NCB
      INTEGER NSLAVES, NASS
      REAL Wtotal, Wblsize
      REAL MUMPS_BLOC2_COUT
      EXTERNAL          MUMPS_BLOC2_COUT
      NASS = NFRONT - NCB
      NSLAVES  = SLAVEF-1
      IF ( K48 .EQ.0 .OR. (K48.EQ.5 .AND. K50.EQ.0)) THEN 
         NSLAVES = max(NCB/max(1,BLSIZE),1)
      ELSE IF (K48.EQ.3 .OR. (K48.EQ.5 .AND. K50.NE.0))THEN
        Wblsize = MUMPS_BLOC2_COUT(BLSIZE,NFRONT,NASS)
        Wtotal  = MUMPS_BLOC2_COUT(NCB,NFRONT,NASS)
        NSLAVES = max(nint ( Wtotal / Wblsize ), 1)
      ENDIF
      MUMPS_BLOC2_GET_NS_BLSIZE = 
     &               min ( NSLAVES,(SLAVEF-1) )
      RETURN
      END FUNCTION MUMPS_BLOC2_GET_NS_BLSIZE
      SUBROUTINE  MUMPS_BLOC2_SET_POSK483( 
     &    GETPOSITIONS, NSLAVES, NFRONT, NCB, 
     &    KMIN, KMAX, SLAVEF,
     &    NBROWMAX, MAXSURFCB, TABPOS, SIZETABPOS) 
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: GETPOSITIONS, 
     &    NSLAVES, NFRONT, NCB, 
     &    KMIN, KMAX, SLAVEF, SIZETABPOS 
      INTEGER, INTENT (OUT) :: NBROWMAX
      INTEGER(8), INTENT(OUT) :: MAXSURFCB
      INTEGER, INTENT (OUT) :: TABPOS(SIZETABPOS) 
      REAL W, COSTni
      REAL delta
      INTEGER  SumNi, NCOLim1, I, BLSIZE, NASS
      LOGICAL GETROW, GETSURF, GETPOS, GET_AVGROW, GET_AVGSURF
      REAL MUMPS_BLOC2_COUT
      EXTERNAL          MUMPS_BLOC2_COUT
      GETROW = (GETPOSITIONS.EQ.1)
      GETSURF= (GETPOSITIONS.EQ.2)
      GETPOS = (GETPOSITIONS.EQ.3)
      GET_AVGROW = (GETPOSITIONS.EQ.4)
      GET_AVGSURF = (GETPOSITIONS.EQ.5)
      NBROWMAX  = 0
      MAXSURFCB = 0_8
      IF (GETPOS) THEN 
        TABPOS (1) = 1
        TABPOS (NSLAVES+1)= NCB+1
        TABPOS (SLAVEF+2) = NSLAVES
      ENDIF
      IF (NSLAVES.EQ.1) THEN
       IF ( GETSURF ) THEN
         NBROWMAX  = NCB
         MAXSURFCB = int(NCB,8)*int(NCB,8)
       ELSEIF ( GETROW ) THEN
         NBROWMAX  = NCB
       ENDIF
      ELSE
        NASS    = NFRONT - NCB
        W       = MUMPS_BLOC2_COUT(NCB,NFRONT,NASS)
        SumNi   = 0
        NCOLim1 = NASS
        DO I = 1, NSLAVES-1
          delta   = real(2*NCOLim1-NASS+1)**2 + 
     &                  (real(4)*W)/real(NASS*(NSLAVES-I+1))
          delta   = sqrt(delta)
          delta   = (real(-2*NCOLim1+NASS-1) + delta )/real(2)
          BLSIZE  = max(int(delta), 1)
          IF ( (NFRONT-NCOLim1-BLSIZE) .LE. NSLAVES-I ) THEN
            BLSIZE = 1
          ENDIF
          NCOLim1 = NCOLim1+BLSIZE
          COSTni  = MUMPS_BLOC2_COUT(BLSIZE,NCOLim1,NASS)
          W       = W - COSTni
          IF (GETPOS) TABPOS(I) = SumNi + 1
          IF (GETSURF) THEN
            NBROWMAX  = max ( NBROWMAX, 
     &       BLSIZE )
            MAXSURFCB = max ( MAXSURFCB, 
     &       int(BLSIZE,8)* int(SumNi+BLSIZE,8) )
          ELSEIF ( GETROW ) THEN         
            NBROWMAX  = max ( NBROWMAX, 
     &       BLSIZE )
             RETURN
          ELSEIF (GET_AVGSURF) THEN
            NBROWMAX = NBROWMAX + BLSIZE
            MAXSURFCB = MAXSURFCB + int(BLSIZE,8)*int(SumNi+BLSIZE,8) 
          ELSEIF (GET_AVGROW) THEN
             NBROWMAX = NBROWMAX + BLSIZE
          ENDIF
          SumNi   = SumNi + BLSIZE
        ENDDO
        BLSIZE = NCB - SumNi
        IF (BLSIZE.LE.0) THEN
          write(*,*) ' Error in MUMPS_BLOC2_SET_POSK483: ', 
     &     ' size lastbloc ', BLSIZE
          CALL MUMPS_ABORT()
        ENDIF
        if (NCOLim1+BLSIZE.NE.NFRONT) then
          write(*,*) ' Error in MUMPS_BLOC2_SET_POSK483: ', 
     &     ' NCOLim1, BLSIZE, NFRONT=', 
     &       NCOLim1, BLSIZE, NFRONT
          CALL MUMPS_ABORT()
        endif
        IF (GETPOS) TABPOS(NSLAVES) = SumNi + 1
        IF (GETSURF) THEN
            NBROWMAX  = max ( NBROWMAX, 
     &       BLSIZE )
            MAXSURFCB = max ( MAXSURFCB, 
     &       int(BLSIZE,8)* int(SumNi+BLSIZE,8 ))
        ELSEIF ( GETROW ) THEN         
            NBROWMAX  = max ( NBROWMAX, 
     &       BLSIZE )
        ELSEIF (GET_AVGSURF) THEN
          NBROWMAX = NBROWMAX + BLSIZE
          MAXSURFCB = MAXSURFCB + int(BLSIZE,8)*int(SumNi+BLSIZE,8)
          NBROWMAX=(NBROWMAX+NSLAVES-1)/NSLAVES
          MAXSURFCB=(MAXSURFCB+int(NSLAVES-1,8))/int(NSLAVES,8)
        ELSEIF (GET_AVGROW) THEN
          NBROWMAX = NBROWMAX + BLSIZE
          NBROWMAX=(NBROWMAX+NSLAVES-1)/NSLAVES
        ENDIF
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_BLOC2_SET_POSK483
      SUBROUTINE MUMPS_BLOC2_SETPARTITION(      
     &            KEEP,KEEP8, SLAVEF,
     &            TAB_POS_IN_PERE,
     &            NSLAVES, NFRONT, NCB
     &             )
      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: NCB, NSLAVES, SLAVEF, NFRONT,
     &                         KEEP(500)
      INTEGER(8) KEEP8(150)
      INTEGER TAB_POS_IN_PERE(SLAVEF+2)
      INTEGER :: I, BLSIZE
      INTEGER KMIN, KMAX, NBROWDUMMY,
     &        GETPOSITIONS, SIZECOLTAB
      INTEGER(8) MAXSURFDUMMY8
      INTEGER MUMPS_GETKMIN, MUMPS_REG_GETKMAX 
      EXTERNAL MUMPS_GETKMIN, MUMPS_REG_GETKMAX, 
     &        MUMPS_BLOC2_SET_POSK483
       IF (KEEP(48).EQ.0) THEN
        BLSIZE = NCB / NSLAVES
        TAB_POS_IN_PERE( 1 ) = 1
        DO I = 1, NSLAVES-1
          TAB_POS_IN_PERE( I+1 ) = TAB_POS_IN_PERE(I) +
     &    BLSIZE
        ENDDO
        TAB_POS_IN_PERE(NSLAVES+1) = NCB+1
        TAB_POS_IN_PERE(SLAVEF+2)  = NSLAVES
        RETURN
      ELSE IF (KEEP(48).EQ.3 ) THEN
        KMAX = MUMPS_REG_GETKMAX(KEEP8(21), NCB)
        KMIN = MUMPS_GETKMIN(KEEP8(21), KEEP(50), KMAX, NCB)
        GETPOSITIONS = 3
        SIZECOLTAB       = SLAVEF+2
        CALL  MUMPS_BLOC2_SET_POSK483(
     &    GETPOSITIONS, NSLAVES, NFRONT, NCB, 
     &    KMIN, KMAX, SLAVEF,
     &    NBROWDUMMY, MAXSURFDUMMY8,
     &    TAB_POS_IN_PERE(1), SIZECOLTAB)
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_BLOC2_SETPARTITION
      SUBROUTINE MUMPS_BLOC2_GET_SLAVE_INFO(
     &            KEEP,KEEP8, INODE, STEP, N, SLAVEF,
     &            ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &            ISLAVE, NCB, NSLAVES, SIZE, FIRST_INDEX )
      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: ISLAVE, NCB, NSLAVES, SLAVEF,
     &                         KEEP(500), INODE, N 
      INTEGER(8) KEEP8(150)
      INTEGER, INTENT( IN ) :: STEP(N), 
     &          ISTEP_TO_INIV2(KEEP(71)), 
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER, INTENT( OUT ):: SIZE, FIRST_INDEX
      INTEGER BLSIZE, J
      IF (KEEP(48).EQ.0) THEN
       BLSIZE = NCB / NSLAVES
       IF ( ISLAVE .NE. NSLAVES ) THEN
        SIZE = BLSIZE
       ELSE
        SIZE = BLSIZE + mod( NCB, NSLAVES )
       END IF
       FIRST_INDEX = ( ISLAVE - 1 ) * BLSIZE + 1
      ELSEIF (KEEP(48).EQ.3) THEN
       J = ISTEP_TO_INIV2 ( STEP(INODE) )
       FIRST_INDEX = TAB_POS_IN_PERE (ISLAVE,J)
       SIZE        = TAB_POS_IN_PERE (ISLAVE+1,J) - FIRST_INDEX 
      ELSEIF (KEEP(48).EQ.4) THEN
         J = ISTEP_TO_INIV2 ( STEP(INODE) )
         FIRST_INDEX = TAB_POS_IN_PERE (ISLAVE,J)
         SIZE = TAB_POS_IN_PERE (ISLAVE+1,J) - FIRST_INDEX 
      ELSEIF (KEEP(48).EQ.5) THEN
         J = ISTEP_TO_INIV2 ( STEP(INODE) )
         FIRST_INDEX = TAB_POS_IN_PERE (ISLAVE,J)
         SIZE = TAB_POS_IN_PERE (ISLAVE+1,J) - FIRST_INDEX 
      ELSE
       WRITE(*,*) 'Error in MUMPS_BLOC2 undef strat'  
       CALL MUMPS_ABORT()
      ENDIF
      RETURN
      END SUBROUTINE MUMPS_BLOC2_GET_SLAVE_INFO
      REAL FUNCTION MUMPS_BLOC2_COUT(NROW,NCOL,NASS)
      IMPLICIT NONE
      INTEGER, INTENT (IN) :: NROW,NCOL,NASS
      MUMPS_BLOC2_COUT = real(NASS)*real(NROW)*
     &                 real(2*NCOL - NASS - NROW + 1)
      RETURN
      END FUNCTION MUMPS_BLOC2_COUT
      INTEGER FUNCTION MUMPS_REG_GET_NSLAVES 
     &      (K821, K48, K50, SLAVEF, 
     &      NCB, NFRONT, NSLAVES_less, NMB_OF_CAND, K375, K119)
      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: NCB, NFRONT, NSLAVES_less, 
     &                  K48, K50, SLAVEF, NMB_OF_CAND, K375, K119
      INTEGER(8), INTENT(IN) :: K821
      INTEGER NSLAVES
      INTEGER KMAX, NPIV,
     &        NSLAVES_ref, NSLAVES_max
      REAL WK_MASTER, WK_SLAVE
      INTEGER  MUMPS_REG_GETKMAX, MUMPS_BLOC2_GET_NSLAVESMIN, 
     &         MUMPS_BLOC2_GET_NSLAVESMAX
      REAL  MUMPS_BLOC2_COUT
      EXTERNAL MUMPS_REG_GETKMAX, MUMPS_BLOC2_GET_NSLAVESMIN, 
     &         MUMPS_BLOC2_GET_NSLAVESMAX
      EXTERNAL MUMPS_BLOC2_COUT
      IF (NMB_OF_CAND.LE.0) THEN
      ENDIF
      IF ( (K48.EQ.0).OR. (K48.EQ.3) ) THEN
         KMAX = MUMPS_REG_GETKMAX( K821, NCB )
         NSLAVES_ref = MUMPS_BLOC2_GET_NSLAVESMIN(
     &     SLAVEF, K48, K821, K50, NFRONT, NCB, K375, K119 )
         NSLAVES = NSLAVES_ref
         IF ( NSLAVES_ref.LT.SLAVEF ) THEN
           NSLAVES_max = MUMPS_BLOC2_GET_NSLAVESMAX(
     &       SLAVEF, K48, K821, K50, NFRONT, NCB, K375, K119 )
           IF ( NSLAVES_max .LT. NSLAVES_less ) THEN
            NSLAVES =  NSLAVES_max
           ELSE 
            NSLAVES =  NSLAVES_less
           ENDIF
           NSLAVES = max(NSLAVES_ref,NSLAVES)
         ENDIF
         NSLAVES = min (NSLAVES, NMB_OF_CAND)
         IF ( NSLAVES.GT.NSLAVES_ref) THEN
          NPIV = NFRONT - NCB
          IF ( K50.EQ.0 ) THEN
           WK_SLAVE = real( NPIV ) * real( NCB ) *
     &         ( 2.0E0 * real(NFRONT) - real(NPIV) )
     &         / real(NSLAVES)
           WK_MASTER = 0.66667E0 *
     &                 real(NPIV)*real(NPIV)*real(NPIV)+
     &                 real(NPIV)*real(NPIV)*real(NCB)
          ELSE
           WK_SLAVE = MUMPS_BLOC2_COUT(NCB,NFRONT,NPIV)
     &         / real(NSLAVES)
           WK_MASTER =  real(NPIV)*real(NPIV)*real(NPIV)/3.0E0
          ENDIF
          IF ( (WK_MASTER.GT.WK_SLAVE).AND.
     &     (WK_SLAVE.GT.1.0E0) ) THEN
           NSLAVES = 
     &         int( real(NSLAVES) * (WK_SLAVE/WK_MASTER))
           NSLAVES = max(NSLAVES_ref, NSLAVES)
          ENDIF
         ENDIF
      ELSE
       NSLAVES = NSLAVES_less
      ENDIF
      NSLAVES = min (NSLAVES, NCB)
      NSLAVES = min (NSLAVES, NMB_OF_CAND)
      MUMPS_REG_GET_NSLAVES = NSLAVES
      RETURN
      END FUNCTION MUMPS_REG_GET_NSLAVES
      SUBROUTINE MUMPS_BLOC2_GET_ISLAVE( 
     &   KEEP,KEEP8, INODE, STEP, N, SLAVEF, 
     &   ISTEP_TO_INIV2, TAB_POS_IN_PERE,
     &
     &   NASS, NCB, 
     &   NSLAVES, POSITION, ISLAVE, IPOSSLAVE )
      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: KEEP(500),INODE,N,SLAVEF 
      INTEGER(8) KEEP8(150)
      INTEGER, INTENT( IN ) :: STEP(N), 
     &          ISTEP_TO_INIV2(KEEP(71)), 
     &          TAB_POS_IN_PERE(SLAVEF+2,max(1,KEEP(56)))
      INTEGER, INTENT( IN  ) :: NASS, NCB, 
     &                          NSLAVES, POSITION
      INTEGER, INTENT( OUT ) :: ISLAVE, IPOSSLAVE
      INTEGER BLSIZE, J, ISHIFT
      IF ((NSLAVES.LE.0).OR.(POSITION.LE.NASS)) THEN
       ISLAVE = 0
       IPOSSLAVE = POSITION
       RETURN
      ENDIF
      IF (KEEP(48).NE.0.and.KEEP(48).NE.3.and.KEEP(48).NE.4
     &     .and.KEEP(48).NE.5) THEN
        WRITE(*,*) 'Error in MUMPS_BLOC2_GET_ISLAVE: undef strat'  
        CALL MUMPS_ABORT()
      ENDIF
      IF (KEEP(48).ne.0) THEN
        J = ISTEP_TO_INIV2 ( STEP(INODE) )
        ISHIFT = POSITION - NASS
        DO ISLAVE = NSLAVES,1,-1
          IF ( ISHIFT .GE. TAB_POS_IN_PERE(ISLAVE,J)) THEN
             IPOSSLAVE = ISHIFT - TAB_POS_IN_PERE(ISLAVE,J) + 1
             EXIT
          END IF
        END DO
      ELSE
        BLSIZE = NCB / NSLAVES
        ISLAVE    = min( NSLAVES,
     &               ( POSITION - NASS - 1 ) / BLSIZE + 1 )
         IPOSSLAVE = POSITION - NASS - ( ISLAVE - 1 ) * BLSIZE
       ENDIF
      RETURN
      END SUBROUTINE MUMPS_BLOC2_GET_ISLAVE
      INTEGER FUNCTION MUMPS_GETKMIN( K821, K50, KMAX, NCB )
      IMPLICIT NONE
      INTEGER, INTENT( IN    )  :: KMAX, NCB, K50
      INTEGER(8), INTENT(IN) :: K821
      INTEGER KMIN, MINGRAN
      INTEGER(8) :: KMINSURF
      IF ( ( NCB .LE.0 ).OR. (KMAX.LE.0) ) THEN 
        MUMPS_GETKMIN = 1
        RETURN
      ENDIF
      IF (K50.EQ.0) THEN
       KMINSURF = 60000_8
#if defined(t3e) || defined(sgi)
       MINGRAN = 40 
#else
       MINGRAN = 50 
#endif
      ELSE
       KMINSURF = 30000_8
#if defined(t3e) || defined(sgi)
       MINGRAN = 10 
#else
       MINGRAN = 20 
#endif
      ENDIF
      IF (K821.GT.0_8) THEN
#if defined(t3e) || defined(sgi)
           KMIN = max(MINGRAN,KMAX/10)
#else
           KMIN = max(MINGRAN,KMAX/20)
#endif
      ELSE
           KMINSURF = max( abs(K821)/500_8, KMINSURF )
           KMIN     = max(
     &                     int( KMINSURF / int(max(NCB,1),8) ),
     &                     1
     &                   )
      ENDIF
      KMIN = min(KMIN,KMAX)
      KMIN = max(KMIN,1)  
      MUMPS_GETKMIN = KMIN
      RETURN
      END FUNCTION MUMPS_GETKMIN
      INTEGER FUNCTION MUMPS_REG_GETKMAX( KEEP821, NCB )
      IMPLICIT NONE
      INTEGER,    intent( in    )  :: NCB
      INTEGER(8), intent( in    )  :: KEEP821
      INTEGER KMAX 
      IF ( NCB .LE.0 ) THEN 
        MUMPS_REG_GETKMAX = 1
        RETURN
      ENDIF
      IF ( KEEP821.GT.0_8 ) THEN
       KMAX = int(KEEP821)
      ELSE
       KMAX =  -int(KEEP821/int(NCB,8))
      ENDIF
      KMAX = min (NCB, KMAX)
      MUMPS_REG_GETKMAX = max ( KMAX, 1 )
      RETURN
      END FUNCTION MUMPS_REG_GETKMAX
