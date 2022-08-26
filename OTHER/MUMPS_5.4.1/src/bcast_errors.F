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
      SUBROUTINE MUMPS_PROPINFO( ICNTL, INFO, COMM, ID )
      INTEGER ICNTL(60), INFO(80), COMM, ID
      INCLUDE 'mpif.h'
#if defined(WORKAROUNDINTELILP64MPI2INTEGER)
      INTEGER(4) :: IN( 2 ), OUT( 2 )
#else
      INTEGER :: IN( 2 ), OUT( 2 )
#endif
      INTEGER ::  LP, IERR
      LP      = ICNTL( 1 )
      IN( 1 ) = INFO ( 1 )
      IN( 2 ) = ID
      CALL MPI_ALLREDUCE( IN, OUT, 1, MPI_2INTEGER, MPI_MINLOC,
     &                    COMM, IERR)
      IF ( OUT( 1 ) .LT. 0 .and. INFO(1) .GE. 0 ) THEN
        INFO( 1 ) = -001
        INFO( 2 ) = OUT( 2 )
      END IF
      RETURN
      END SUBROUTINE MUMPS_PROPINFO
