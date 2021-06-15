C
C  This file is part of MUMPS 5.4.0, released
C  on Tue Apr 13 15:26:30 UTC 2021
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
       MODULE ZMUMPS_FAC_S_IS_POINTERS_M
C      ----------------------------------
C      This module defines a type used in 
C      ZMUMPS_FAC_DRIVER and ZMUMPS_FAC_B
C      ----------------------------------
       TYPE S_IS_POINTERS_T
          COMPLEX(kind=8), POINTER, DIMENSION(:) :: A
          INTEGER, POINTER, DIMENSION(:) :: IW
       END TYPE S_IS_POINTERS_T
       END MODULE ZMUMPS_FAC_S_IS_POINTERS_M
