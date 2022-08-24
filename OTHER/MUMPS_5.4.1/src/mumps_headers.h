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
C
C     Common header positions:
C 
C     XXI    ->  size of integer record
C     XXR    ->  size of real record
C     XXS    ->  status of the node
C     XXN    ->  node number
C     XXP    ->  pointer to previous record
C     XXA    ->  active fronts data management
C     XXF    ->  blr data passed from factorization to solve
C     XXLR   ->  Low rank status of a node (0=FR, 
C                                           1=LowRank CB only
C                                           2=LowRank factors/panels only
C                                           3=LowRank CB+factor/panel)
C     XXEBF  ->  End of Blocfacto (0=not yet, 1=finished)  
C     XXD    ->  dynamic data size
C     XXG    ->  GPU information (currently number of pinned rows NFRONT-NBROWS_CPU)
C REMARK: .h file could be replaced by a module with functions to get node status
C          added in the module.
C 
      INTEGER, PARAMETER :: XXI = 0, XXR = 1, XXS = 3, XXN = 4, XXP = 5
      INTEGER, PARAMETER :: XXA = 6, XXF = 7 
      INTEGER, PARAMETER :: XXLR = 8
      INTEGER, PARAMETER :: XXNBPR = 9
      INTEGER, PARAMETER :: XXEBF = 10
      INTEGER, PARAMETER :: XXD = 11
      INTEGER, PARAMETER :: XXG = 13
C 
C     Size of header in incore and out-of-core
C
      INTEGER XSIZE_IC, XSIZE_OOC_SYM, XSIZE_OOC_UNSYM
      INTEGER XSIZE_OOC_NOPANEL ! To store virtual addresses
C     At the moment, all headers are of the same size because
C     no OOC specific information are stored in header.
CM     other OOC specific information directly in the headers.
      PARAMETER (XSIZE_IC=14,XSIZE_OOC_SYM=14,XSIZE_OOC_UNSYM=14,
     &           XSIZE_OOC_NOPANEL=14)
C
C     -------------------------------------------------------
C     Position of header size (formerly XSIZE) in KEEP array.
C     KEEP(IXSZ) is set at the beginning of the factorization
C     to either XSIZE_IC, XSIZE_OOC_SYM or XSIZE_OOC_UNSYM.
C     -------------------------------------------------------
      INTEGER IXSZ
      PARAMETER(IXSZ= 222)    ! KEEP(222) used
      INTEGER, PARAMETER :: S_CB1COMP = 314
      INTEGER S_ACTIVE, S_ALL, S_NOLCBCONTIG,
     &        S_NOLCBNOCONTIG, S_NOLCLEANED,
     &        S_NOLCBNOCONTIG38, S_NOLCBCONTIG38,
     &        S_NOLCLEANED38, 
     &        S_NOLNOCB, S_NOLNOCBCLEANED,
     &        C_FINI
      PARAMETER(S_ACTIVE=400, S_ALL=401, S_NOLCBCONTIG=402,
     &          S_NOLCBNOCONTIG=403, S_NOLCLEANED=404,
     &          S_NOLCBNOCONTIG38=405, S_NOLCBCONTIG38=406,
     &          S_NOLCLEANED38=407, 
     &          S_NOLNOCB=408, S_NOLNOCBCLEANED=409,
     &          C_FINI=1)
      INTEGER, PARAMETER :: S_FREE = 54321
      INTEGER, PARAMETER :: S_NOTFREE = -123
      INTEGER, PARAMETER :: TOP_OF_STACK = -999999
      INTEGER XTRA_SLAVES_SYM, XTRA_SLAVES_UNSYM
      PARAMETER(XTRA_SLAVES_SYM=4, XTRA_SLAVES_UNSYM=2)
         INTEGER S_ROOT2SON_CALLED, S_REC_CONTSTATIC, 
     &  S_ROOTBAND_INIT
         PARAMETER(S_ROOT2SON_CALLED=-341,S_REC_CONTSTATIC=1,
     &             S_ROOTBAND_INIT=0)
