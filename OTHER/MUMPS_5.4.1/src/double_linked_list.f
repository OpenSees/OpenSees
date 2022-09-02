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
      MODULE MUMPS_IDLL
      IMPLICIT NONE
      TYPE IDLL_NODE_T
          TYPE ( IDLL_NODE_T ), POINTER :: NEXT, PREV
          INTEGER ELMT
      END TYPE IDLL_NODE_T
      TYPE IDLL_T
          TYPE ( IDLL_NODE_T ), POINTER :: FRONT, BACK
      END TYPE IDLL_T
      CONTAINS
      FUNCTION IDLL_CREATE(DLL)
          INTEGER :: IDLL_CREATE
#if defined(MUMPS_F2003)
          TYPE ( IDLL_T ), POINTER, INTENT ( OUT ) :: DLL
#else
          TYPE ( IDLL_T ), POINTER :: DLL
#endif
          INTEGER IERR
          ALLOCATE ( DLL, STAT=IERR )
          IF ( IERR .NE. 0 ) THEN
              IDLL_CREATE = -2
              RETURN
          END IF
          NULLIFY ( DLL%FRONT )
          NULLIFY ( DLL%BACK )
          IDLL_CREATE = 0
          RETURN
      END FUNCTION IDLL_CREATE
      FUNCTION IDLL_DESTROY(DLL)
          INTEGER :: IDLL_DESTROY
#if defined(MUMPS_F2003)
          TYPE ( IDLL_T ), POINTER, INTENT ( OUT ) :: DLL
#else
          TYPE ( IDLL_T ), POINTER :: DLL
#endif
          TYPE ( IDLL_NODE_T ), POINTER :: AUX
          IF ( .NOT. associated ( DLL ) ) THEN
              IDLL_DESTROY = -1
              RETURN
          END IF
          DO WHILE ( associated ( DLL%FRONT ) )
              AUX => DLL%FRONT
              DLL%FRONT => DLL%FRONT%NEXT
              DEALLOCATE( AUX )
          END DO
          DEALLOCATE( DLL )
          IDLL_DESTROY = 0
      END FUNCTION IDLL_DESTROY
      FUNCTION IDLL_PUSH_FRONT(DLL, ELMT)
          INTEGER :: IDLL_PUSH_FRONT
#if defined(MUMPS_F2003)
          TYPE ( IDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
#else
          TYPE ( IDLL_T ), POINTER :: DLL
#endif
          INTEGER, INTENT ( IN ) :: ELMT
          TYPE ( IDLL_NODE_T ), POINTER :: NODE
          INTEGER IERR
          IF ( .NOT. associated ( DLL ) ) THEN
              IDLL_PUSH_FRONT = -1
              RETURN
          END IF
          ALLOCATE( NODE, STAT=IERR )
          IF ( IERR .NE. 0 ) THEN
              IDLL_PUSH_FRONT = -2
              RETURN
          END IF
          NODE%ELMT = ELMT
          NODE%NEXT => DLL%FRONT
          NULLIFY ( NODE%PREV )
          IF ( associated ( DLL%FRONT ) ) THEN
              DLL%FRONT%PREV => NODE
          END IF
          DLL%FRONT => NODE
          IF ( .NOT. associated ( DLL%BACK ) ) THEN
              DLL%BACK => NODE
          END IF
          IDLL_PUSH_FRONT = 0
      END FUNCTION IDLL_PUSH_FRONT
      FUNCTION IDLL_POP_FRONT(DLL, ELMT)
          INTEGER :: IDLL_POP_FRONT
#if defined(MUMPS_F2003)
          TYPE ( IDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
#else
          TYPE ( IDLL_T ), POINTER :: DLL
#endif
          INTEGER, INTENT ( OUT ) :: ELMT
          TYPE ( IDLL_NODE_T ), POINTER :: AUX
          IF ( .NOT. associated ( DLL ) ) THEN
              IDLL_POP_FRONT = -1
              RETURN
          END IF
          IF ( .NOT. associated ( DLL%FRONT ) ) THEN
              IDLL_POP_FRONT = -3
              RETURN
          END IF
          ELMT = DLL%FRONT%ELMT
          AUX => DLL%FRONT
          DLL%FRONT => DLL%FRONT%NEXT
          IF ( associated ( DLL%FRONT ) ) THEN
              NULLIFY ( DLL%FRONT%PREV )
          END IF
          IF ( associated ( DLL%BACK, AUX ) ) THEN
              NULLIFY ( DLL%BACK )
          END IF
          DEALLOCATE ( AUX )
          IDLL_POP_FRONT = 0
      END FUNCTION IDLL_POP_FRONT
      FUNCTION IDLL_PUSH_BACK(DLL, ELMT)
          INTEGER :: IDLL_PUSH_BACK
#if defined(MUMPS_F2003)
          TYPE ( IDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
#else
          TYPE ( IDLL_T ), POINTER :: DLL
#endif
          INTEGER, INTENT ( IN ) :: ELMT
          TYPE ( IDLL_NODE_T ), POINTER :: NODE
          INTEGER IERR
          IF ( .NOT. associated ( DLL ) ) THEN
              IDLL_PUSH_BACK = -1
              RETURN
          END IF
          ALLOCATE( NODE, STAT=IERR )
          IF ( IERR .NE. 0 ) THEN
              IDLL_PUSH_BACK = -2
              RETURN
          END IF
          NODE%ELMT = ELMT
          NULLIFY ( NODE%NEXT )
          NODE%PREV => DLL%BACK
          IF ( associated ( DLL%BACK ) ) THEN
              DLL%BACK%NEXT => NODE
          END IF
          DLL%BACK => NODE
          IF ( .NOT. associated ( DLL%FRONT ) ) THEN
              DLL%FRONT => NODE
          END IF
          IDLL_PUSH_BACK = 0
      END FUNCTION IDLL_PUSH_BACK
      FUNCTION IDLL_POP_BACK(DLL, ELMT)
          INTEGER :: IDLL_POP_BACK
#if defined(MUMPS_F2003)
          TYPE ( IDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
#else
          TYPE ( IDLL_T ), POINTER :: DLL
#endif
          INTEGER, INTENT ( OUT ) :: ELMT
          TYPE ( IDLL_NODE_T ), POINTER :: AUX
          IF ( .NOT. associated ( DLL ) ) THEN
              IDLL_POP_BACK = -1
              RETURN
          END IF
          IF ( .NOT. associated ( DLL%BACK ) ) THEN
              IDLL_POP_BACK = -3
              RETURN
          END IF
          ELMT = DLL%BACK%ELMT
          AUX => DLL%BACK
          DLL%BACK => DLL%BACK%PREV
          IF ( associated ( DLL%BACK ) ) THEN
              NULLIFY ( DLL%BACK%NEXT )
          END IF
          IF ( associated ( DLL%FRONT, AUX ) ) THEN
              NULLIFY ( DLL%FRONT )
          END IF
          DEALLOCATE ( AUX )
          IDLL_POP_BACK = 0
      END FUNCTION IDLL_POP_BACK
      FUNCTION IDLL_INSERT(DLL, POS, ELMT)
          INTEGER :: IDLL_INSERT
#if defined(MUMPS_F2003)
          TYPE ( IDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
#else
          TYPE ( IDLL_T ), POINTER :: DLL
#endif
          INTEGER, INTENT ( IN ) :: POS, ELMT
          TYPE ( IDLL_NODE_T ), POINTER :: NODE
          TYPE ( IDLL_NODE_T ), POINTER :: NEW_PTR, OLD_PTR
          INTEGER :: IERR, CPT
          IF ( .NOT. associated ( DLL ) ) THEN
              IDLL_INSERT = -1
              RETURN
          END IF
          IF ( POS .LE. 0 ) THEN
              IDLL_INSERT = -4
              RETURN
          END IF
          CPT = 1
          NEW_PTR => DLL%FRONT
          NULLIFY ( OLD_PTR )
          DO WHILE ( ( CPT .LT. POS ) .AND.
     &               ( associated ( NEW_PTR ) ) )
              OLD_PTR => NEW_PTR
              NEW_PTR => NEW_PTR%NEXT
              CPT = CPT + 1
          END DO
          ALLOCATE ( NODE, STAT=IERR )
          IF ( IERR .NE. 0 ) THEN
              IDLL_INSERT = -2
              RETURN
          END IF
          NODE%ELMT = ELMT
          IF ( .NOT. associated ( OLD_PTR ) ) THEN
              IF ( .NOT. associated ( NEW_PTR ) ) THEN
                  NULLIFY ( NODE%PREV )
                  NULLIFY ( NODE%NEXT )
                  DLL%FRONT => NODE
                  DLL%BACK => NODE
              ELSE
                  NULLIFY ( NODE%PREV )
                  NODE%NEXT => NEW_PTR
                  NEW_PTR%PREV => NODE
                  DLL%FRONT => NODE
              END IF
          ELSE
              IF ( .NOT. associated ( NEW_PTR ) ) THEN
                  NODE%PREV => OLD_PTR
                  NULLIFY ( NODE%NEXT )
                  OLD_PTR%NEXT => NODE
                  DLL%BACK => NODE
              ELSE
                  NODE%PREV => OLD_PTR
                  NODE%NEXT => NEW_PTR
                  OLD_PTR%NEXT => NODE
                  NEW_PTR%PREV => NODE
              END IF
          END IF
          IDLL_INSERT = 0
      END FUNCTION IDLL_INSERT
      FUNCTION IDLL_INSERT_BEFORE(DLL, NODE_AFTER, ELMT)
          INTEGER :: IDLL_INSERT_BEFORE
#if defined(MUMPS_F2003)
          TYPE ( IDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
          TYPE ( IDLL_NODE_T ), POINTER, INTENT ( IN ) :: NODE_AFTER
#else
          TYPE ( IDLL_T ), POINTER :: DLL
          TYPE ( IDLL_NODE_T ), POINTER :: NODE_AFTER
#endif
          INTEGER, INTENT ( IN ) :: ELMT
          TYPE ( IDLL_NODE_T ), POINTER :: NODE_BEFORE
          INTEGER :: IERR
          ALLOCATE ( NODE_BEFORE, STAT=IERR )
          IF ( IERR .NE. 0 ) THEN
              IDLL_INSERT_BEFORE = -2
              RETURN
          END IF
          NODE_BEFORE%ELMT = ELMT
          IF ( .NOT. associated ( NODE_AFTER%PREV ) ) THEN
              NODE_AFTER%PREV => NODE_BEFORE
              NODE_BEFORE%NEXT => NODE_AFTER
              NULLIFY ( NODE_BEFORE%PREV )
              DLL%FRONT => NODE_BEFORE
          ELSE
              NODE_BEFORE%NEXT => NODE_AFTER
              NODE_BEFORE%PREV => NODE_AFTER%PREV
              NODE_AFTER%PREV => NODE_BEFORE
              NODE_BEFORE%PREV%NEXT => NODE_BEFORE
          END IF
          IDLL_INSERT_BEFORE = 0
      END FUNCTION IDLL_INSERT_BEFORE
      FUNCTION IDLL_INSERT_AFTER(DLL, NODE_BEFORE, ELMT)
          INTEGER :: IDLL_INSERT_AFTER
#if defined(MUMPS_F2003)
          TYPE ( IDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
          TYPE ( IDLL_NODE_T ), POINTER, INTENT ( IN ) :: NODE_BEFORE
#else
          TYPE ( IDLL_T ), POINTER :: DLL
          TYPE ( IDLL_NODE_T ), POINTER :: NODE_BEFORE
#endif
          INTEGER, INTENT ( IN ) :: ELMT
          TYPE ( IDLL_NODE_T ), POINTER :: NODE_AFTER
          INTEGER :: IERR
          ALLOCATE ( NODE_AFTER, STAT=IERR )
          IF ( IERR .NE. 0 ) THEN
              IDLL_INSERT_AFTER = -2
              RETURN
          END IF
          NODE_AFTER%ELMT = ELMT
          IF ( .NOT. associated ( NODE_BEFORE%NEXT ) ) THEN
              NODE_BEFORE%NEXT => NODE_AFTER
              NODE_AFTER%PREV => NODE_BEFORE
              NULLIFY ( NODE_AFTER%NEXT )
              DLL%BACK => NODE_AFTER
          ELSE
              NODE_AFTER%PREV => NODE_BEFORE
              NODE_AFTER%NEXT => NODE_BEFORE%NEXT
              NODE_BEFORE%NEXT => NODE_AFTER
              NODE_AFTER%NEXT%PREV => NODE_AFTER
          END IF
          IDLL_INSERT_AFTER = 0
      END FUNCTION IDLL_INSERT_AFTER
      FUNCTION IDLL_LOOKUP (DLL, POS, ELMT)
          INTEGER :: IDLL_LOOKUP
#if defined(MUMPS_F2003)
          TYPE ( IDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
#else
          TYPE ( IDLL_T ), POINTER :: DLL
#endif
          INTEGER, INTENT ( IN ) :: POS
          INTEGER, INTENT ( OUT ) :: ELMT
          TYPE ( IDLL_NODE_T ), POINTER :: AUX
          INTEGER :: CPT
          IF ( .NOT. associated ( DLL ) ) THEN
              IDLL_LOOKUP = -1
              RETURN
          END IF
          IF ( POS .LE. 0 ) THEN
              IDLL_LOOKUP = -4
              RETURN
          END IF
          CPT = 1
          AUX => DLL%FRONT
          DO WHILE ( ( CPT .LT. POS ) .AND. ( associated ( AUX ) ) )
              CPT = CPT + 1
              AUX => AUX%NEXT
          END DO
          IF ( .NOT. associated ( AUX ) ) THEN
              IDLL_LOOKUP = -3
              RETURN
          END IF
          ELMT = AUX%ELMT
          IDLL_LOOKUP = 0
      END FUNCTION IDLL_LOOKUP
      FUNCTION IDLL_REMOVE_POS(DLL, POS, ELMT)
          INTEGER :: IDLL_REMOVE_POS
#if defined(MUMPS_F2003)
          TYPE ( IDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
#else
          TYPE ( IDLL_T ), POINTER :: DLL
#endif
          INTEGER, INTENT ( IN ) :: POS
          INTEGER, INTENT ( OUT ) :: ELMT
          TYPE ( IDLL_NODE_T ), POINTER :: AUX
          INTEGER :: CPT
          IF ( .NOT. associated ( DLL ) ) THEN
              IDLL_REMOVE_POS = -1
              RETURN
          END IF
          CPT = 1
          AUX => DLL%FRONT
          DO WHILE ( ( associated ( AUX ) ) .AND.
     &               ( CPT .LT. POS ) )
              CPT = CPT + 1
              AUX => AUX%NEXT
          END DO
          IF ( associated ( AUX ) ) THEN
              IF ( .NOT. associated ( AUX%PREV ) ) THEN
                  IF ( .NOT. associated ( AUX%NEXT ) ) THEN
                      NULLIFY ( DLL%FRONT )
                      NULLIFY ( DLL%BACK )
                  ELSE
                      NULLIFY ( AUX%NEXT%PREV )
                      DLL%FRONT => AUX%NEXT
                  END IF
              ELSE
                  IF ( .NOT. associated ( AUX%NEXT ) ) THEN
                      NULLIFY ( AUX%PREV%NEXT )
                      DLL%BACK => AUX%PREV
                  ELSE
                      AUX%PREV%NEXT => AUX%NEXT
                      AUX%NEXT%PREV => AUX%PREV
                  END IF
              END IF
              ELMT = AUX%ELMT
              DEALLOCATE ( AUX )
          ELSE
              IDLL_REMOVE_POS = -3
              RETURN
          END IF
          IDLL_REMOVE_POS = 0
      END FUNCTION IDLL_REMOVE_POS
      FUNCTION IDLL_REMOVE_ELMT(DLL, ELMT, POS)
          INTEGER :: IDLL_REMOVE_ELMT
#if defined(MUMPS_F2003)
          TYPE ( IDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
#else
          TYPE ( IDLL_T ), POINTER :: DLL
#endif
          INTEGER, INTENT ( IN ) :: ELMT
          INTEGER, INTENT ( OUT ) :: POS
          TYPE ( IDLL_NODE_T ), POINTER :: AUX
          INTEGER :: CPT
          IF ( .NOT. associated ( DLL ) ) THEN
              IDLL_REMOVE_ELMT = -1
              RETURN
          END IF
          CPT = 1
          AUX => DLL%FRONT
          DO WHILE ( ( associated ( AUX ) ) .AND.
     &               ( AUX%ELMT .NE. ELMT ) )
              CPT = CPT + 1
              AUX => AUX%NEXT
          END DO
          IF ( associated ( AUX ) ) THEN
              IF ( .NOT. associated ( AUX%PREV ) ) THEN
                  IF ( .NOT. associated ( AUX%NEXT ) ) THEN
                      NULLIFY ( DLL%FRONT )
                      NULLIFY ( DLL%BACK )
                  ELSE
                      NULLIFY ( AUX%NEXT%PREV )
                      DLL%FRONT => AUX%NEXT
                  END IF
              ELSE
                  IF ( .NOT. associated ( AUX%NEXT ) ) THEN
                      NULLIFY ( AUX%PREV%NEXT )
                      DLL%BACK => AUX%PREV
                  ELSE
                      AUX%PREV%NEXT => AUX%NEXT
                      AUX%NEXT%PREV => AUX%PREV
                  END IF
              END IF
              POS = CPT
              DEALLOCATE ( AUX )
          ELSE
              IDLL_REMOVE_ELMT = -3
              RETURN
          END IF
          IDLL_REMOVE_ELMT = 0
      END FUNCTION IDLL_REMOVE_ELMT
      FUNCTION IDLL_LENGTH(DLL)
          INTEGER :: IDLL_LENGTH
#if defined(MUMPS_F2003)
          TYPE ( IDLL_T ), POINTER, INTENT ( IN ) :: DLL
#else
          TYPE ( IDLL_T ), POINTER :: DLL
#endif
          INTEGER :: LENGTH
          TYPE ( IDLL_NODE_T ), POINTER :: AUX
          LENGTH = 0
          IF ( .NOT. associated ( DLL ) ) THEN
              IDLL_LENGTH = -1
              RETURN
          END IF
          AUX => DLL%FRONT
          DO WHILE ( associated ( AUX ) )
              LENGTH = LENGTH + 1
              AUX => AUX%NEXT
          END DO
          IDLL_LENGTH = LENGTH
      END FUNCTION IDLL_LENGTH
      FUNCTION IDLL_ITERATOR_BEGIN(DLL, PTR)
          INTEGER :: IDLL_ITERATOR_BEGIN
#if defined(MUMPS_F2003)
          TYPE ( IDLL_T ), POINTER, INTENT ( IN ) :: DLL
          TYPE ( IDLL_NODE_T ), POINTER, INTENT ( OUT ) :: PTR
#else
          TYPE ( IDLL_T ), POINTER  :: DLL
          TYPE ( IDLL_NODE_T ), POINTER :: PTR
#endif
          IF ( .NOT. associated ( DLL ) ) THEN
              IDLL_ITERATOR_BEGIN = -1
              RETURN
          END IF
          PTR => DLL%FRONT
          IDLL_ITERATOR_BEGIN = 0
      END FUNCTION IDLL_ITERATOR_BEGIN
      FUNCTION IDLL_ITERATOR_END(DLL, PTR)
          INTEGER :: IDLL_ITERATOR_END
#if defined(MUMPS_F2003)
          TYPE ( IDLL_T ), POINTER, INTENT ( IN ) :: DLL
          TYPE ( IDLL_NODE_T ), POINTER, INTENT ( OUT ) :: PTR
#else
          TYPE ( IDLL_T ), POINTER :: DLL
          TYPE ( IDLL_NODE_T ), POINTER :: PTR
#endif
          IF ( .NOT. associated ( DLL ) ) THEN
              IDLL_ITERATOR_END = -1
              RETURN
          END IF
          PTR => DLL%BACK
          IDLL_ITERATOR_END = 0
      END FUNCTION IDLL_ITERATOR_END
      FUNCTION IDLL_IS_EMPTY(DLL)
          LOGICAL :: IDLL_IS_EMPTY
#if defined(MUMPS_F2003)
          TYPE ( IDLL_T ), POINTER, INTENT ( IN ) :: DLL
#else
          TYPE ( IDLL_T ), POINTER :: DLL
#endif
          IDLL_IS_EMPTY = ( associated ( DLL%FRONT ) )
      END FUNCTION IDLL_IS_EMPTY
      FUNCTION IDLL_2_ARRAY(DLL, ARRAY, LENGTH)
          INTEGER :: IDLL_2_ARRAY
#if defined(MUMPS_F2003)
          TYPE ( IDLL_T ), POINTER, INTENT ( IN ) :: DLL
          INTEGER, POINTER, DIMENSION (:), INTENT ( OUT ) :: ARRAY
#else
          TYPE ( IDLL_T ), POINTER :: DLL
          INTEGER, POINTER, DIMENSION (:) :: ARRAY
#endif
          INTEGER, INTENT ( OUT ) :: LENGTH
          TYPE ( IDLL_NODE_T ), POINTER :: AUX
          INTEGER :: I, IERR
          IF ( .NOT. associated ( DLL ) ) THEN
              IDLL_2_ARRAY = -1
              RETURN
          END IF
          LENGTH = IDLL_LENGTH(DLL)
          ALLOCATE ( ARRAY ( max(1,LENGTH) ), STAT=IERR )
          IF ( IERR .NE. 0 ) THEN
              IDLL_2_ARRAY = -2
              RETURN
          END IF
          I = 1
          AUX => DLL%FRONT
          DO WHILE ( associated ( AUX ) )
              ARRAY ( I ) = AUX%ELMT
              I = I + 1
              AUX => AUX%NEXT
          END DO
          IDLL_2_ARRAY = 0
      END FUNCTION IDLL_2_ARRAY
      END MODULE MUMPS_IDLL
      MODULE MUMPS_DDLL
      IMPLICIT NONE
      TYPE DDLL_NODE_T
          TYPE ( DDLL_NODE_T ), POINTER :: NEXT, PREV
          DOUBLE PRECISION :: ELMT
      END TYPE DDLL_NODE_T
      TYPE DDLL_T
          TYPE ( DDLL_NODE_T ), POINTER :: FRONT, BACK
      END TYPE DDLL_T
      CONTAINS
      FUNCTION DDLL_CREATE(DLL)
          INTEGER :: DDLL_CREATE
#if defined(MUMPS_F2003)
          TYPE ( DDLL_T ), POINTER, INTENT ( OUT ) :: DLL
#else
          TYPE ( DDLL_T ), POINTER :: DLL
#endif
          INTEGER IERR
          ALLOCATE ( DLL, STAT=IERR )
          IF ( IERR .NE. 0 ) THEN
              DDLL_CREATE = -2
              RETURN
          END IF
          NULLIFY ( DLL%FRONT )
          NULLIFY ( DLL%BACK )
          DDLL_CREATE = 0
          RETURN
      END FUNCTION DDLL_CREATE
      FUNCTION DDLL_DESTROY(DLL)
          INTEGER :: DDLL_DESTROY
#if defined(MUMPS_F2003)
          TYPE ( DDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
#else
          TYPE ( DDLL_T ), POINTER :: DLL
#endif
          TYPE ( DDLL_NODE_T ), POINTER :: AUX
          IF ( .NOT. associated ( DLL ) ) THEN
              DDLL_DESTROY = -1
              RETURN
          END IF
          DO WHILE ( associated ( DLL%FRONT ) )
              AUX => DLL%FRONT
              DLL%FRONT => DLL%FRONT%NEXT
              DEALLOCATE( AUX )
          END DO
          DEALLOCATE( DLL )
          DDLL_DESTROY = 0
      END FUNCTION DDLL_DESTROY
      FUNCTION DDLL_PUSH_FRONT(DLL, ELMT)
          INTEGER :: DDLL_PUSH_FRONT
#if defined(MUMPS_F2003)
          TYPE ( DDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
#else
          TYPE ( DDLL_T ), POINTER :: DLL
#endif
          DOUBLE PRECISION, INTENT ( IN ) :: ELMT
          TYPE ( DDLL_NODE_T ), POINTER :: NODE
          INTEGER IERR
          IF ( .NOT. associated ( DLL ) ) THEN
              DDLL_PUSH_FRONT = -1
              RETURN
          END IF
          ALLOCATE( NODE, STAT=IERR )
          IF ( IERR .NE. 0 ) THEN
              DDLL_PUSH_FRONT = -2
              RETURN
          END IF
          NODE%ELMT = ELMT
          NODE%NEXT => DLL%FRONT
          NULLIFY ( NODE%PREV )
          IF ( associated ( DLL%FRONT ) ) THEN
              DLL%FRONT%PREV => NODE
          END IF
          DLL%FRONT => NODE
          IF ( .NOT. associated ( DLL%BACK ) ) THEN
              DLL%BACK => NODE
          END IF
          DDLL_PUSH_FRONT = 0
      END FUNCTION DDLL_PUSH_FRONT
      FUNCTION DDLL_POP_FRONT(DLL, ELMT)
          INTEGER :: DDLL_POP_FRONT
#if defined(MUMPS_F2003)
          TYPE ( DDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
#else
          TYPE ( DDLL_T ), POINTER :: DLL
#endif
          DOUBLE PRECISION, INTENT ( OUT ) :: ELMT
          TYPE ( DDLL_NODE_T ), POINTER :: AUX
          IF ( .NOT. associated ( DLL ) ) THEN
              DDLL_POP_FRONT = -1
              RETURN
          END IF
          IF ( .NOT. associated ( DLL%FRONT ) ) THEN
              DDLL_POP_FRONT = -3
              RETURN
          END IF
          ELMT = DLL%FRONT%ELMT
          AUX => DLL%FRONT
          DLL%FRONT => DLL%FRONT%NEXT
          IF ( associated ( DLL%FRONT ) ) THEN
              NULLIFY ( DLL%FRONT%PREV )
          END IF
          IF ( associated ( DLL%BACK, AUX ) ) THEN
              NULLIFY ( DLL%BACK )
          END IF
          DEALLOCATE ( AUX )
          DDLL_POP_FRONT = 0
      END FUNCTION DDLL_POP_FRONT
      FUNCTION DDLL_PUSH_BACK(DLL, ELMT)
          INTEGER :: DDLL_PUSH_BACK
#if defined(MUMPS_F2003)
          TYPE ( DDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
#else
          TYPE ( DDLL_T ), POINTER :: DLL
#endif
          DOUBLE PRECISION, INTENT ( IN ) :: ELMT
          TYPE ( DDLL_NODE_T ), POINTER :: NODE
          INTEGER IERR
          IF ( .NOT. associated ( DLL ) ) THEN
              DDLL_PUSH_BACK = -1
              RETURN
          END IF
          ALLOCATE( NODE, STAT=IERR )
          IF ( IERR .NE. 0 ) THEN
              DDLL_PUSH_BACK = -2
              RETURN
          END IF
          NODE%ELMT = ELMT
          NULLIFY ( NODE%NEXT )
          NODE%PREV => DLL%BACK
          IF ( associated ( DLL%BACK ) ) THEN
              DLL%BACK%NEXT => NODE
          END IF
          DLL%BACK => NODE
          IF ( .NOT. associated ( DLL%FRONT ) ) THEN
              DLL%FRONT => NODE
          END IF
          DDLL_PUSH_BACK = 0
      END FUNCTION DDLL_PUSH_BACK
      FUNCTION DDLL_POP_BACK(DLL, ELMT)
          INTEGER :: DDLL_POP_BACK
#if defined(MUMPS_F2003)
          TYPE ( DDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
#else
          TYPE ( DDLL_T ), POINTER :: DLL
#endif
          DOUBLE PRECISION, INTENT ( OUT ) :: ELMT
          TYPE ( DDLL_NODE_T ), POINTER :: AUX
          IF ( .NOT. associated ( DLL ) ) THEN
              DDLL_POP_BACK = -1
              RETURN
          END IF
          IF ( .NOT. associated ( DLL%BACK ) ) THEN
              DDLL_POP_BACK = -3
              RETURN
          END IF
          ELMT = DLL%BACK%ELMT
          AUX => DLL%BACK
          DLL%BACK => DLL%BACK%PREV
          IF ( associated ( DLL%BACK ) ) THEN
              NULLIFY ( DLL%BACK%NEXT )
          END IF
          IF ( associated ( DLL%FRONT, AUX ) ) THEN
              NULLIFY ( DLL%FRONT )
          END IF
          DEALLOCATE ( AUX )
          DDLL_POP_BACK = 0
      END FUNCTION DDLL_POP_BACK
      FUNCTION DDLL_INSERT(DLL, POS, ELMT)
          INTEGER :: DDLL_INSERT
#if defined(MUMPS_F2003)
          TYPE ( DDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
#else
          TYPE ( DDLL_T ), POINTER :: DLL
#endif
          INTEGER, INTENT ( IN ) :: POS
          DOUBLE PRECISION , INTENT ( IN ) :: ELMT
          TYPE ( DDLL_NODE_T ), POINTER :: NODE
          TYPE ( DDLL_NODE_T ), POINTER :: NEW_PTR, OLD_PTR
          INTEGER :: IERR, CPT
          IF ( .NOT. associated ( DLL ) ) THEN
              DDLL_INSERT = -1
              RETURN
          END IF
          IF ( POS .LE. 0 ) THEN
              DDLL_INSERT = -4
              RETURN
          END IF
          CPT = 1
          NEW_PTR => DLL%FRONT
          NULLIFY ( OLD_PTR )
          DO WHILE ( ( CPT .LT. POS ) .AND.
     &               ( associated ( NEW_PTR ) ) )
              OLD_PTR => NEW_PTR
              NEW_PTR => NEW_PTR%NEXT
              CPT = CPT + 1
          END DO
          ALLOCATE ( NODE, STAT=IERR )
          IF ( IERR .NE. 0 ) THEN
              DDLL_INSERT = -2
              RETURN
          END IF
          NODE%ELMT = ELMT
          IF ( .NOT. associated ( OLD_PTR ) ) THEN
              IF ( .NOT. associated ( NEW_PTR ) ) THEN
                  NULLIFY ( NODE%PREV )
                  NULLIFY ( NODE%NEXT )
                  DLL%FRONT => NODE
                  DLL%BACK => NODE
              ELSE
                  NULLIFY ( NODE%PREV )
                  NODE%NEXT => NEW_PTR
                  NEW_PTR%PREV => NODE
                  DLL%FRONT => NODE
              END IF
          ELSE
              IF ( .NOT. associated ( NEW_PTR ) ) THEN
                  NODE%PREV => OLD_PTR
                  NULLIFY ( NODE%NEXT )
                  OLD_PTR%NEXT => NODE
                  DLL%BACK => NODE
              ELSE
                  NODE%PREV => OLD_PTR
                  NODE%NEXT => NEW_PTR
                  OLD_PTR%NEXT => NODE
                  NEW_PTR%PREV => NODE
              END IF
          END IF
          DDLL_INSERT = 0
      END FUNCTION DDLL_INSERT
      FUNCTION DDLL_INSERT_BEFORE(DLL, NODE_AFTER, ELMT)
          INTEGER :: DDLL_INSERT_BEFORE
#if defined(MUMPS_F2003)
          TYPE ( DDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
          TYPE ( DDLL_NODE_T ), POINTER, INTENT ( IN ) :: NODE_AFTER
#else
          TYPE ( DDLL_T ), POINTER :: DLL
          TYPE ( DDLL_NODE_T ), POINTER :: NODE_AFTER
#endif
          DOUBLE PRECISION, INTENT ( IN ) :: ELMT
          TYPE ( DDLL_NODE_T ), POINTER :: NODE_BEFORE
          INTEGER :: IERR
          ALLOCATE ( NODE_BEFORE, STAT=IERR )
          IF ( IERR .NE. 0 ) THEN
              DDLL_INSERT_BEFORE = -2
              RETURN
          END IF
          NODE_BEFORE%ELMT = ELMT
          IF ( .NOT. associated ( NODE_AFTER%PREV ) ) THEN
              NODE_AFTER%PREV => NODE_BEFORE
              NODE_BEFORE%NEXT => NODE_AFTER
              NULLIFY ( NODE_BEFORE%PREV )
              DLL%FRONT => NODE_BEFORE
          ELSE
              NODE_BEFORE%NEXT => NODE_AFTER
              NODE_BEFORE%PREV => NODE_AFTER%PREV
              NODE_AFTER%PREV => NODE_BEFORE
              NODE_BEFORE%PREV%NEXT => NODE_BEFORE
          END IF
          DDLL_INSERT_BEFORE = 0
      END FUNCTION DDLL_INSERT_BEFORE
      FUNCTION DDLL_INSERT_AFTER(DLL, NODE_BEFORE, ELMT)
          INTEGER :: DDLL_INSERT_AFTER
#if defined(MUMPS_F2003)
          TYPE ( DDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
          TYPE ( DDLL_NODE_T ), POINTER, INTENT ( IN ) :: NODE_BEFORE
#else
          TYPE ( DDLL_T ), POINTER :: DLL
          TYPE ( DDLL_NODE_T ), POINTER :: NODE_BEFORE
#endif
          DOUBLE PRECISION, INTENT ( IN ) :: ELMT
          TYPE ( DDLL_NODE_T ), POINTER :: NODE_AFTER
          INTEGER :: IERR
          ALLOCATE ( NODE_AFTER, STAT=IERR )
          IF ( IERR .NE. 0 ) THEN
              DDLL_INSERT_AFTER = -2
              RETURN
          END IF
          NODE_AFTER%ELMT = ELMT
          IF ( .NOT. associated ( NODE_BEFORE%NEXT ) ) THEN
              NODE_BEFORE%NEXT => NODE_AFTER
              NODE_AFTER%PREV => NODE_BEFORE
              NULLIFY ( NODE_AFTER%NEXT )
              DLL%BACK => NODE_AFTER
          ELSE
              NODE_AFTER%PREV => NODE_BEFORE
              NODE_AFTER%NEXT => NODE_BEFORE%NEXT
              NODE_BEFORE%NEXT => NODE_AFTER
              NODE_AFTER%NEXT%PREV => NODE_AFTER
          END IF
          DDLL_INSERT_AFTER = 0
      END FUNCTION DDLL_INSERT_AFTER
      FUNCTION DDLL_LOOKUP (DLL, POS, ELMT)
          INTEGER :: DDLL_LOOKUP
#if defined(MUMPS_F2003)
          TYPE ( DDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
#else
          TYPE ( DDLL_T ), POINTER :: DLL
#endif
          INTEGER, INTENT ( IN ) :: POS
          DOUBLE PRECISION, INTENT ( OUT ) :: ELMT
          TYPE ( DDLL_NODE_T ), POINTER :: AUX
          INTEGER :: CPT
          IF ( .NOT. associated ( DLL ) ) THEN
              DDLL_LOOKUP = -1
              RETURN
          END IF
          IF ( POS .LE. 0 ) THEN
              DDLL_LOOKUP = -4
              RETURN
          END IF
          CPT = 1
          AUX => DLL%FRONT
          DO WHILE ( ( CPT .LT. POS ) .AND. ( associated ( AUX ) ) )
              CPT = CPT + 1
              AUX => AUX%NEXT
          END DO
          IF ( .NOT. associated ( AUX ) ) THEN
              DDLL_LOOKUP = -3
              RETURN
          END IF
          ELMT = AUX%ELMT
          DDLL_LOOKUP = 0
      END FUNCTION DDLL_LOOKUP
      FUNCTION DDLL_REMOVE_POS(DLL, POS, ELMT)
          INTEGER :: DDLL_REMOVE_POS
#if defined(MUMPS_F2003)
          TYPE ( DDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
#else
          TYPE ( DDLL_T ), POINTER :: DLL
#endif
          INTEGER, INTENT ( IN ) :: POS
          DOUBLE PRECISION, INTENT ( OUT ) :: ELMT
          TYPE ( DDLL_NODE_T ), POINTER :: AUX
          INTEGER :: CPT
          IF ( .NOT. associated ( DLL ) ) THEN
              DDLL_REMOVE_POS = -1
              RETURN
          END IF
          CPT = 1
          AUX => DLL%FRONT
          DO WHILE ( ( associated ( AUX ) ) .AND.
     &               ( CPT .LT. POS ) )
              CPT = CPT + 1
              AUX => AUX%NEXT
          END DO
          IF ( associated ( AUX ) ) THEN
              IF ( .NOT. associated ( AUX%PREV ) ) THEN
                  IF ( .NOT. associated ( AUX%NEXT ) ) THEN
                      NULLIFY ( DLL%FRONT )
                      NULLIFY ( DLL%BACK )
                  ELSE
                      NULLIFY ( AUX%NEXT%PREV )
                      DLL%FRONT => AUX%NEXT
                  END IF
              ELSE
                  IF ( .NOT. associated ( AUX%NEXT ) ) THEN
                      NULLIFY ( AUX%PREV%NEXT )
                      DLL%BACK => AUX%PREV
                  ELSE
                      AUX%PREV%NEXT => AUX%NEXT
                      AUX%NEXT%PREV => AUX%PREV
                  END IF
              END IF
              ELMT = AUX%ELMT
              DEALLOCATE ( AUX )
          ELSE
              DDLL_REMOVE_POS = -3
              RETURN
          END IF
          DDLL_REMOVE_POS = 0
      END FUNCTION DDLL_REMOVE_POS
      FUNCTION DDLL_REMOVE_ELMT(DLL, ELMT, POS)
          INTEGER :: DDLL_REMOVE_ELMT
#if defined(MUMPS_F2003)
          TYPE ( DDLL_T ), POINTER, INTENT ( INOUT ) :: DLL
#else
          TYPE ( DDLL_T ), POINTER :: DLL
#endif
          DOUBLE PRECISION, INTENT ( IN ) :: ELMT
          INTEGER, INTENT ( OUT ) :: POS
          TYPE ( DDLL_NODE_T ), POINTER :: AUX
          INTEGER :: CPT
          IF ( .NOT. associated ( DLL ) ) THEN
              DDLL_REMOVE_ELMT = -1
              RETURN
          END IF
          CPT = 1
          AUX => DLL%FRONT
          DO WHILE ( ( associated ( AUX ) ) .AND.
     &               ( AUX%ELMT .NE. ELMT ) )
              CPT = CPT + 1
              AUX => AUX%NEXT
          END DO
          IF ( associated ( AUX ) ) THEN
              IF ( .NOT. associated ( AUX%PREV ) ) THEN
                  IF ( .NOT. associated ( AUX%NEXT ) ) THEN
                      NULLIFY ( DLL%FRONT )
                      NULLIFY ( DLL%BACK )
                  ELSE
                      NULLIFY ( AUX%NEXT%PREV )
                      DLL%FRONT => AUX%NEXT
                  END IF
              ELSE
                  IF ( .NOT. associated ( AUX%NEXT ) ) THEN
                      NULLIFY ( AUX%PREV%NEXT )
                      DLL%BACK => AUX%PREV
                  ELSE
                      AUX%PREV%NEXT => AUX%NEXT
                      AUX%NEXT%PREV => AUX%PREV
                  END IF
              END IF
              POS = CPT
              DEALLOCATE ( AUX )
          ELSE
              DDLL_REMOVE_ELMT = -3
              RETURN
          END IF
          DDLL_REMOVE_ELMT = 0
      END FUNCTION DDLL_REMOVE_ELMT
      FUNCTION DDLL_LENGTH(DLL)
          INTEGER :: DDLL_LENGTH
#if defined(MUMPS_F2003)
          TYPE ( DDLL_T ), POINTER, INTENT ( IN ) :: DLL
#else
          TYPE ( DDLL_T ), POINTER :: DLL
#endif
          INTEGER :: LENGTH
          TYPE ( DDLL_NODE_T ), POINTER :: AUX
          IF ( .NOT. associated ( DLL ) ) THEN
              DDLL_LENGTH = -1
              RETURN
          END IF
          LENGTH = 0
          AUX => DLL%FRONT
          DO WHILE ( associated ( AUX ) )
              LENGTH = LENGTH + 1
              AUX => AUX%NEXT
          END DO
          DDLL_LENGTH = LENGTH
      END FUNCTION DDLL_LENGTH
      FUNCTION DDLL_ITERATOR_BEGIN(DLL, PTR)
          INTEGER :: DDLL_ITERATOR_BEGIN
#if defined(MUMPS_F2003)
          TYPE ( DDLL_T ), POINTER, INTENT ( IN ) :: DLL
          TYPE ( DDLL_NODE_T ), POINTER, INTENT ( OUT ) :: PTR
#else
          TYPE ( DDLL_T ), POINTER :: DLL
          TYPE ( DDLL_NODE_T ), POINTER :: PTR
#endif
          IF ( .NOT. associated ( DLL ) ) THEN
              DDLL_ITERATOR_BEGIN = -1
              RETURN
          END IF
          PTR => DLL%FRONT
          DDLL_ITERATOR_BEGIN = 0
      END FUNCTION DDLL_ITERATOR_BEGIN
      FUNCTION DDLL_ITERATOR_END(DLL, PTR)
          INTEGER :: DDLL_ITERATOR_END
#if defined(MUMPS_F2003)
          TYPE ( DDLL_T ), POINTER, INTENT ( IN ) :: DLL
          TYPE ( DDLL_NODE_T ), POINTER, INTENT ( OUT ) :: PTR
#else
          TYPE ( DDLL_T ), POINTER :: DLL
          TYPE ( DDLL_NODE_T ), POINTER :: PTR
#endif
          IF ( .NOT. associated ( DLL ) ) THEN
              DDLL_ITERATOR_END = -1
              RETURN
          END IF
          PTR => DLL%BACK
          DDLL_ITERATOR_END = 0
      END FUNCTION DDLL_ITERATOR_END
      FUNCTION DDLL_IS_EMPTY(DLL)
          LOGICAL :: DDLL_IS_EMPTY
#if defined(MUMPS_F2003)
          TYPE ( DDLL_T ), POINTER, INTENT ( IN ) :: DLL
#else
          TYPE ( DDLL_T ), POINTER :: DLL
#endif
          DDLL_IS_EMPTY = ( associated ( DLL%FRONT ) )
      END FUNCTION DDLL_IS_EMPTY
      FUNCTION DDLL_2_ARRAY(DLL, ARRAY, LENGTH)
          INTEGER :: DDLL_2_ARRAY
#if defined(MUMPS_F2003)
          TYPE ( DDLL_T ), POINTER, INTENT ( IN ) :: DLL
          DOUBLE PRECISION, POINTER, DIMENSION(:), INTENT(OUT) :: ARRAY
#else
          TYPE ( DDLL_T ), POINTER :: DLL
          DOUBLE PRECISION, POINTER, DIMENSION(:) :: ARRAY
#endif
          INTEGER, INTENT ( OUT ) :: LENGTH
          TYPE ( DDLL_NODE_T ), POINTER :: AUX
          INTEGER :: I, IERR
          IF ( .NOT. associated ( DLL ) ) THEN
              DDLL_2_ARRAY = -1
              RETURN
          END IF
          LENGTH = DDLL_LENGTH(DLL)
          ALLOCATE ( ARRAY ( max(1,LENGTH) ), STAT=IERR )
          IF ( IERR .NE. 0 ) THEN
              DDLL_2_ARRAY = -2
              RETURN
          END IF
          I = 1
          AUX => DLL%FRONT
          DO WHILE ( associated ( AUX ) )
              ARRAY ( I ) = AUX%ELMT
              I = I + 1
              AUX => AUX%NEXT
          END DO
          DDLL_2_ARRAY = 0
      END FUNCTION DDLL_2_ARRAY
      END MODULE MUMPS_DDLL
