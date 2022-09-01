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
      SUBROUTINE DMUMPS_COMPACT_FACTORS(A, LDA, NPIV, NBROW, K50,
     &           SIZEA )
      IMPLICIT NONE
      INTEGER LDA, NPIV, NBROW, K50
      INTEGER(8), INTENT(IN) :: SIZEA
      DOUBLE PRECISION A(SIZEA)
      INTEGER(8) :: IOLD, INEW, J8
      INTEGER I , ILAST
      INTEGER NBROW_L_RECTANGLE_TO_MOVE
      IF ((NPIV.EQ.0).OR.(LDA.EQ.NPIV)) GOTO 500
      IF ( K50.NE.0 ) THEN
        IOLD = int(LDA  + 1,8)
        INEW = int(NPIV + 1,8)
        IF (IOLD .EQ. INEW ) THEN
          INEW = INEW + int(NPIV,8) * int(NPIV - 1,8)
          IOLD = IOLD + int(LDA,8) * int(NPIV - 1,8)
        ELSE
          DO I = 1, NPIV - 1
            IF ( I .LE. NPIV-2 ) THEN
              ILAST = I+1
            ELSE
              ILAST = I
            ENDIF
            DO J8 = 0_8, int(ILAST,8)
              A( INEW + J8 ) = A( IOLD + J8 )
            END DO
            INEW = INEW + int(NPIV,8)
            IOLD = IOLD + int(LDA,8)
          END DO
        ENDIF
        NBROW_L_RECTANGLE_TO_MOVE = NBROW
      ELSE 
        INEW = 1_8 + int(NPIV,8) * int(LDA + 1,8)
        IOLD = 1_8 + int(LDA,8) * int(NPIV +1,8)
        NBROW_L_RECTANGLE_TO_MOVE = NBROW - 1
      ENDIF
      DO I = 1, NBROW_L_RECTANGLE_TO_MOVE
         DO J8 = 0_8, int(NPIV - 1,8)
           A( INEW + J8 ) = A( IOLD + J8 )
         END DO
         INEW = INEW + int(NPIV,8)
         IOLD = IOLD + int(LDA,8)
      ENDDO
 500  RETURN
      END SUBROUTINE DMUMPS_COMPACT_FACTORS
      SUBROUTINE DMUMPS_COMPACT_FACTORS_UNSYM(A, LDA, NPIV, NCONTIG,
     &                                        SIZEA )
      IMPLICIT NONE
      INTEGER,    INTENT(IN)     :: NCONTIG, NPIV, LDA
      INTEGER(8), INTENT(IN)     :: SIZEA
      DOUBLE PRECISION,    INTENT(INOUT)  :: A(SIZEA)
      INTEGER I, J
      INTEGER(8) :: INEW, IOLD
      INEW = int(NPIV+1,8)
      IOLD = int(LDA+1,8)
      DO I = 2, NCONTIG
        DO J = 1, NPIV
          A(INEW)=A(IOLD)
          INEW = INEW + 1_8
          IOLD = IOLD + 1_8
        ENDDO
        IOLD = IOLD + int(LDA - NPIV,8)
      ENDDO
      RETURN
      END SUBROUTINE DMUMPS_COMPACT_FACTORS_UNSYM
      SUBROUTINE DMUMPS_COPY_CB_RIGHT_TO_LEFT( A, LA, LDA, POSELT,
     &           IPTRLU, NPIV,
     &           NBCOL_STACK, NBROW_STACK,
     &           NBROW_SEND, SIZECB, KEEP, PACKED_CB,
     &           LAST_ALLOWED, NBROW_ALREADY_STACKED )
      IMPLICIT NONE
      INTEGER(8), intent (in) :: POSELT, IPTRLU, LA, SIZECB
      LOGICAL, intent (in) :: PACKED_CB
      DOUBLE PRECISION A(LA)
      INTEGER, intent(in):: LDA, NPIV, NBCOL_STACK, NBROW_STACK,
     &                      NBROW_SEND
      INTEGER, intent(inout) :: NBROW_ALREADY_STACKED
      INTEGER(8), intent(in)    :: LAST_ALLOWED
      INTEGER(8) :: APOS, NPOS
      INTEGER NBROW
      INTEGER(8) :: J
      INTEGER I, KEEP(500)
#if defined(ZERO_TRIANGLE)
      DOUBLE PRECISION ZERO
        PARAMETER( ZERO = 0.0D0 )
#endif
      NBROW = NBROW_STACK + NBROW_SEND
      IF (NBROW_STACK .NE. 0 ) THEN
        NPOS = IPTRLU + SIZECB         
        APOS = POSELT + int(NPIV+NBROW,8) 
     &       * int(LDA,8) - 1_8 
        IF ( KEEP(50) .EQ. 0 .OR. .NOT. PACKED_CB ) THEN
          APOS = APOS - int(LDA,8) * int(NBROW_ALREADY_STACKED,8)
          NPOS = NPOS
     &         - int(NBCOL_STACK,8) * int(NBROW_ALREADY_STACKED,8)
        ELSE
          APOS = APOS - int(LDA - 1,8) * int(NBROW_ALREADY_STACKED,8)
          NPOS = NPOS - ( int(NBROW_ALREADY_STACKED,8) *
     &                    int(NBROW_ALREADY_STACKED+1,8) ) / 2_8
        ENDIF
        DO I = NBROW - NBROW_ALREADY_STACKED, NBROW_SEND+1, -1
          IF (KEEP(50).EQ.0) THEN
            IF ( NPOS - int(NBCOL_STACK,8) + 1_8 .LT.
     &                                  LAST_ALLOWED ) THEN
              EXIT
            ENDIF
            DO J= 1_8,int(NBCOL_STACK,8)
              A(NPOS-J+1_8) = A(APOS-J+1_8)
            ENDDO
            NPOS = NPOS - int(NBCOL_STACK,8)
          ELSE
            IF (.NOT. PACKED_CB) THEN
              IF ( NPOS - int(NBCOL_STACK,8) + 1_8 .LT.
     &                                  LAST_ALLOWED ) THEN
                EXIT
              ENDIF
#if defined(ZERO_TRIANGLE)
              DO J = 1_8, int(NBCOL_STACK - I,8)
                A(NPOS - J + 1_8) = ZERO
              END DO
#endif
              NPOS = NPOS + int(- NBCOL_STACK + I,8)
            ENDIF
            IF ( NPOS - int(I,8) + 1_8 .LT. LAST_ALLOWED ) THEN
              EXIT
            ENDIF
            DO J =1_8, int(I,8)
              A(NPOS-J+1_8) = A(APOS-J+1_8)
            ENDDO
            NPOS = NPOS - int(I,8)
          ENDIF
          IF (KEEP(50).EQ.0) THEN
            APOS = APOS - int(LDA,8)
          ELSE
            APOS = APOS - int(LDA + 1,8)
          ENDIF
          NBROW_ALREADY_STACKED = NBROW_ALREADY_STACKED + 1
        ENDDO
      END IF
      RETURN
      END SUBROUTINE DMUMPS_COPY_CB_RIGHT_TO_LEFT
      SUBROUTINE DMUMPS_COPY_CB_LEFT_TO_RIGHT( A, LA, LDA, POSELT,
     &           IPTRLU, NPIV,
     &           NBCOL_STACK, NBROW_STACK,
     &           NBROW_SEND, SIZECB, KEEP, PACKED_CB)
      IMPLICIT NONE
      INTEGER(8), intent (in) :: POSELT, IPTRLU, LA, SIZECB
      LOGICAL, intent (in) :: PACKED_CB
      DOUBLE PRECISION A(LA)
      INTEGER, intent(in):: LDA, NPIV, NBCOL_STACK, NBROW_STACK,
     &                      NBROW_SEND
      INTEGER(8) :: APOS, NPOS, APOS_ini, NPOS_ini
      INTEGER I, KEEP(500)
      INTEGER(8) :: J, LDA8
#if defined(ZERO_TRIANGLE)
      DOUBLE PRECISION ZERO
        PARAMETER( ZERO = 0.0D0 )
#endif
      LDA8 = int(LDA,8)
      NPOS_ini = IPTRLU + 1_8
      APOS_ini = POSELT + int(NPIV+NBROW_SEND,8)* LDA8 + int(NPIV,8)
!$OMP PARALLEL DO PRIVATE(J, NPOS, APOS) IF (NBROW_STACK > KEEP(360))
      DO I = 1, NBROW_STACK
         IF (PACKED_CB) THEN
            NPOS = NPOS_ini + int(I-1,8) * int(I,8)/2_8 +
     &             int(I-1,8) * int(NBROW_SEND,8)
         ELSE
            NPOS = NPOS_ini + int(I-1,8) * int(NBCOL_STACK,8)
        ENDIF
        APOS  =  APOS_ini + int(I-1,8) * LDA8
        IF (KEEP(50).EQ.0) THEN
          DO J = 1_8, int(NBCOL_STACK,8)
            A(NPOS+J-1_8) = A(APOS+J-1_8)
          ENDDO
        ELSE
          DO J  = 1_8, int(I + NBROW_SEND,8)
            A(NPOS+J-1_8)=A(APOS+J-1_8)
          ENDDO
#if defined(ZERO_TRIANGLE)
          IF (.NOT. PACKED_CB) THEN
            A(NPOS+int(I+NBROW_SEND,8):
     &        NPOS+int(NBCOL_STACK-1,8))=ZERO
          ENDIF
#endif
        ENDIF
      ENDDO
!$OMP END PARALLEL DO
      RETURN
      END SUBROUTINE DMUMPS_COPY_CB_LEFT_TO_RIGHT
