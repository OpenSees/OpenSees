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
      MODULE MUMPS_ANA_BLK_M
      TYPE COL_LMATRIX_T
       INTEGER :: NBINCOL   
       INTEGER, POINTER :: IRN(:) => null()  
      END TYPE COL_LMATRIX_T
      TYPE LMATRIX_T
       INTEGER    :: NBCOL
       INTEGER(8) :: NZL
       TYPE(COL_LMATRIX_T), POINTER :: COL(:) => null()
      END TYPE LMATRIX_T
      TYPE COMPACT_GRAPH_T
       INTEGER(8) :: NZG, SIZEADJALLOCATED
       INTEGER    :: NG
       INTEGER(8), POINTER :: IPE(:) => null()  
       INTEGER, POINTER :: ADJ(:) => null()  
      END TYPE COMPACT_GRAPH_T
      END MODULE MUMPS_ANA_BLK_M
