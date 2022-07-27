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
      MODULE MUMPS_L0_OMP_M
        LOGICAL, DIMENSION(:), POINTER :: NB_CORE_PER_THREAD_CHANGED
        INTEGER, DIMENSION(:), POINTER :: NB_CORE_PER_THREAD
        INTEGER :: THREAD_ID
        LOGICAL :: IS_ROOT_OF_L0_OMP
!$OMP   THREADPRIVATE ( THREAD_ID , IS_ROOT_OF_L0_OMP )
      END MODULE MUMPS_L0_OMP_M
