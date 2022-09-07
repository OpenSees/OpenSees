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
      SUBROUTINE DMUMPS_SIMSCALEABS(IRN_loc, JCN_loc, A_loc, NZ_loc,
     &     M, N, NUMPROCS, MYID, COMM,
     &     RPARTVEC, CPARTVEC,
     &     RSNDRCVSZ, CSNDRCVSZ, REGISTRE,
     &     IWRK, IWRKSZ,
     &     INTSZ, RESZ, OP,
     &     ROWSCA, COLSCA, WRKRC, ISZWRKRC,
     &     SYM, NB1, NB2, NB3, EPS,
     &     ONENORMERR,INFNORMERR)
C----------------------------------------------------------------------
C    IF SYM=0 CALLs unsymmetric variant DMUMPS_SIMSCALEABSUNS.
C    IF SYM=2 CALLS symmetric variant where only one of a_ij and a_ji 
C         is stored. DMUMPS_SIMSCALEABSSYM
C---------------------------------------------------------------------
C    For details, see the two subroutines below
C         DMUMPS_SIMSCALEABSUNS and DMUMPS_SIMSCALEABSSYM
C ---------------------------------------------------------------------
C
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER(8) NZ_loc
      INTEGER IWRKSZ, ISZWRKRC
      INTEGER M, N, OP
      INTEGER NUMPROCS, MYID, COMM
      INTEGER INTSZ, RESZ
      INTEGER IRN_loc(NZ_loc)
      INTEGER JCN_loc(NZ_loc)
      DOUBLE PRECISION A_loc(NZ_loc)
      INTEGER RPARTVEC(M)
      INTEGER RSNDRCVSZ(2*NUMPROCS)
      INTEGER CPARTVEC(N) 
      INTEGER CSNDRCVSZ(2*NUMPROCS)
      INTEGER IWRK(IWRKSZ)
      INTEGER REGISTRE(12)
      DOUBLE PRECISION ROWSCA(M)
      DOUBLE PRECISION COLSCA(N)
      DOUBLE PRECISION WRKRC(ISZWRKRC)
      DOUBLE PRECISION ONENORMERR,INFNORMERR
C     LOCALS
C IMPORTANT POINTERS
C     FOR the scaling phase
      INTEGER SYM, NB1, NB2, NB3
      DOUBLE PRECISION EPS
C     EXTERNALS
      EXTERNAL DMUMPS_SIMSCALEABSUNS,DMUMPS_SIMSCALEABSSYM, 
     &     DMUMPS_INITREAL
C     MUST HAVE IT
      INTEGER I
      IF(SYM.EQ.0) THEN
         CALL DMUMPS_SIMSCALEABSUNS(IRN_loc, JCN_loc, A_loc, 
     &        NZ_loc,
     &        M, N, NUMPROCS, MYID, COMM,
     &        RPARTVEC, CPARTVEC,
     &        RSNDRCVSZ, CSNDRCVSZ, REGISTRE,
     &        IWRK, IWRKSZ,
     &        INTSZ, RESZ, OP,
     &        ROWSCA, COLSCA, WRKRC, ISZWRKRC,
     &        NB1, NB2, NB3, EPS,
     &        ONENORMERR, INFNORMERR)  
      ELSE
         CALL DMUMPS_SIMSCALEABSSYM(IRN_loc, JCN_loc, A_loc, 
     &        NZ_loc,
     &        N, NUMPROCS, MYID, COMM,
     &        RPARTVEC, 
     &        RSNDRCVSZ, REGISTRE,
     &        IWRK, IWRKSZ,
     &        INTSZ, RESZ, OP,
     &        ROWSCA, WRKRC, ISZWRKRC,
     &        NB1, NB2, NB3, EPS,
     &        ONENORMERR, INFNORMERR)  
         DO I=1,N
            COLSCA(I) = ROWSCA(I)
         ENDDO
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_SIMSCALEABS
      SUBROUTINE DMUMPS_SIMSCALEABSUNS(IRN_loc, JCN_loc, A_loc, NZ_loc,
     &     M, N, NUMPROCS, MYID, COMM,
     &     RPARTVEC, CPARTVEC,
     &     RSNDRCVSZ, CSNDRCVSZ, REGISTRE,
     &     IWRK, IWRKSZ,
     &     INTSZ, RESZ, OP,
     &     ROWSCA, COLSCA, WRKRC, ISZWRKRC,
     &     NB1, NB2, NB3, EPS,
     &     ONENORMERR, INFNORMERR)    
C----------------------------------------------------------------------
C Input parameters:
C      M, N: size of matrix (in general M=N, but the algorithm 
C            works for rectangular matrices as well (norms other than
C            inf-norm are not possible mathematically in this case).
C      NUMPROCS, MYID, COMM: guess what are those
C      RPARTVEC:  row partvec to be filled when OP=1
C      CPARTVEC:  col partvec to be filled when OP=1
C      RSNDRCVSZ: send recv sizes for row operations. 
C                 to be filled when OP=1
C      CSNDRCVSZ: send recv sizes for col operations. 
C                 to be filled when OP=1
C      REGISTRE:  to store some pointers (size etc)
C      IWRK: working space. when OP=1 IWRKSZ.GE.4*MAXMN
C            when OP=2 INTSZ portion is used. Thus, IWRKSZ>INTSZ 
C            when OP=2
C      IWRKSZ: size
C      INTSZ: to be computed when OP=1, necessary integer space to run 
C             scaling algo when OP=2
C      RESZ:  to be computed when OP=1, necessary real space to run 
C             scaling algo when OP=2
C      OP: 
C          =1 estimation of memory and construction of partvecs
C           writes into RPARTVEC,CPARTVEC,RSNDRCVSZ,CSNDRCVSZ,REGISTRE
C           does not access WRKRC, uses IWRK as workspace
C           computes INTSZ and RESZ.
C          =2 Compute scalings 
C           restores pointers from REGISTRE, 
C           stores communication structure in IWRK (from the start). 
C
C      ROWSCA: space for row scaling factor; has size M
C      COLSCA: space for col scaling factor; has size N
C      WRKRC: real working space. when OP=1, is not accessed. Thus, it
C             can be declared to be of size 1 at OP=1 call.
C      ISZWRKRC: size
C      SYM: is matrix symmetric
C      NB1, NB2, NB3: algo runs 
C                     NB1 iters of inf-norm (default  1/1), 
C                     NB2 iters of 1-norm   (default  3/10),
C                     NB3 iters of inf-norm (default  3/10).
C          in succession.
C      EPS: tolerance for concergence. 
C           IF EPS < 0.R0 then does not test convergence.
C           If convergence occured during the first set of inf-norm
C           iterations, we start performing one-norm iterations.
C           If convergence occured during the one-norm iterations,
C           we start performing the second set of inf-norm iterations.
C           If convergence occured during the second set of inf-norm,
C           we prepare to return.
C     ONENORMERR : error in one norm scaling (associated with the scaling 
C                  arrays of the previous iterations), 
C     INFNORMERR : error in inf norm scaling (associated with the scaling 
C                  arrays of the previous iterations).
C---------------------------------------------------------------------
C On input:
C      OP=1==>Requirements
C             IWRKSZ.GE.4*MAXMN
C             RPARTVEC  of size M
C             CPARTVEC  of size N
C             RSNDRCVSZ of size 2*NUMPROCS
C             CSNDRCVSZ of size 2*NUMPROCS
C             REGISTRE  of size 12
C             
C      OP=2==>Requirements
C             INTSZ .GE. REGISTRE(11)
C             RESZ  .GE. REGISTRE(12) 
C---------------------------------------------------------------------
C On output:
C     ROWSCA and COLSCA 
C            at processor 0 of COMM: complete factors.
C            at other processors   : only the ROWSCA(i) or COLSCA(j) 
C            for which there is a nonzero a_i* or a_*j are useful.
C     ONENORMERR : error in one norm scaling 
C                = -1.0 if iter2=0.
C     INFNORMERR : error in inf norm scaling 
C                = inf norm error at iter3 if iter3 > 0
C                = inf norm error at iter1 if iter1 > 0, iter3=0
C                = -1.0 if iter1=iter3=0
C ---------------------------------------------------------------------
C References:
C     The scaling algorithms are based on those discussed in
C     [1] D. Ruiz, "A scaling algorithm to equilibrate both rows and 
C         columns norms in matrices", Tech. Rep. Rutherford 
C         Appleton Laboratory, Oxon, UK and ENSEEIHT-IRIT, 
C         Toulouse, France, RAL-TR-2001-034 and RT/APO/01/4, 2001.
C     [2] D. Ruiz and B. Ucar, "A symmetry preserving algorithm for
C         matrix scaling", in preparation as of Jan'08.
C
C     The parallelization approach is discussed in
C     [3] P. R. Amestoy, I. S. Duff, D. Ruiz, and B. Ucar,
C         "A parallel matrix scaling algorithm".
C         In proceedings of VECPAR'08-International Meeting-High 
C         Performance Computing for Computational Science, Jan'08.
C     and was supported by ANR-SOLSTICE project (ANR-06-CIS6-010)
C ---------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER(8) NZ_loc
      INTEGER IWRKSZ, INTSZ, ISZWRKRC
      INTEGER M, N, OP
      INTEGER NUMPROCS, MYID, COMM
      INTEGER RESZ
      INTEGER IRN_loc(NZ_loc)
      INTEGER JCN_loc(NZ_loc)
      DOUBLE PRECISION A_loc(NZ_loc)
      INTEGER RPARTVEC(M) 
      INTEGER CPARTVEC(N) 
      INTEGER RSNDRCVSZ(2*NUMPROCS)
      INTEGER CSNDRCVSZ(2*NUMPROCS)
      INTEGER REGISTRE(12)
      INTEGER IWRK(IWRKSZ)
      DOUBLE PRECISION ROWSCA(M)
      DOUBLE PRECISION COLSCA(N)
      DOUBLE PRECISION WRKRC(ISZWRKRC)  
      DOUBLE PRECISION ONENORMERR,INFNORMERR
C     LOCALS
      INTEGER IRSNDRCVNUM, ORSNDRCVNUM
      INTEGER ICSNDRCVNUM, OCSNDRCVNUM
      INTEGER IRSNDRCVVOL, ORSNDRCVVOL
      INTEGER ICSNDRCVVOL, OCSNDRCVVOL
      INTEGER  INUMMYR, INUMMYC
C IMPORTANT POINTERS
      INTEGER IMYRPTR,IMYCPTR 
      INTEGER IRNGHBPRCS, IRSNDRCVIA,IRSNDRCVJA
      INTEGER ORNGHBPRCS, ORSNDRCVIA,ORSNDRCVJA
      INTEGER ICNGHBPRCS, ICSNDRCVIA,ICSNDRCVJA
      INTEGER OCNGHBPRCS, OCSNDRCVIA,OCSNDRCVJA
      INTEGER ISTATUS, REQUESTS, TMPWORK
      INTEGER ITDRPTR, ITDCPTR, ISRRPTR
      INTEGER OSRRPTR, ISRCPTR, OSRCPTR
C     FOR the scaling phase
      INTEGER NB1, NB2, NB3
      DOUBLE PRECISION EPS
C     Iteration vars 
      INTEGER ITER, IR, IC
      INTEGER(8) :: NZIND
      DOUBLE PRECISION ELM
C     COMM TAGS....
      INTEGER TAG_COMM_COL
      PARAMETER(TAG_COMM_COL=100)
      INTEGER TAG_COMM_ROW
      PARAMETER(TAG_COMM_ROW=101)
      INTEGER TAG_ITERS
      PARAMETER(TAG_ITERS=102)
C     FUNCTIONS
      EXTERNAL DMUMPS_CREATEPARTVEC,
     &     DMUMPS_NUMVOLSNDRCV, 
     &     DMUMPS_SETUPCOMMS,
     &     DMUMPS_FINDNUMMYROWCOL, 
     &     DMUMPS_CHKCONVGLO,
     &     DMUMPS_CHK1CONV,
     &     DMUMPS_FILLMYROWCOLINDICES,
     &     DMUMPS_INITREAL,
     &     DMUMPS_INITREALLST,
     &     DMUMPS_DOCOMMINF,
     &     DMUMPS_DOCOMM1N
      INTEGER DMUMPS_CHKCONVGLO 
      INTEGER DMUMPS_CHK1CONV
      DOUBLE PRECISION DMUMPS_ERRSCALOC
      DOUBLE PRECISION DMUMPS_ERRSCA1
      INTRINSIC abs
      DOUBLE PRECISION RONE, RZERO
      PARAMETER(RONE=1.0D0,RZERO=0.0D0)
C     TMP VARS
      INTEGER RESZR, RESZC
      INTEGER INTSZR, INTSZC
      INTEGER MAXMN
      INTEGER I, IERROR
      DOUBLE PRECISION ONEERRROW, ONEERRCOL, ONEERRL, ONEERRG
      DOUBLE PRECISION INFERRROW, INFERRCOL, INFERRL, INFERRG
      INTEGER OORANGEIND
      INFERRG = -RONE
      ONEERRG = -RONE
      OORANGEIND = 0
      MAXMN = M
      IF(MAXMN < N) MAXMN = N
C     Create row partvec and col partvec
      IF(OP == 1) THEN
         IF(NUMPROCS > 1) THEN
C     Check done outside
C     IF(IWRKSZ.LT.4*MAXMN) THEN   ERROR.... 
            CALL DMUMPS_CREATEPARTVEC(MYID, NUMPROCS, COMM,
     &           IRN_loc, JCN_loc, NZ_loc,
     &           RPARTVEC, M, N,
     &           IWRK, IWRKSZ)
            CALL DMUMPS_CREATEPARTVEC(MYID, NUMPROCS, COMM,
     &           JCN_loc, IRN_loc,  NZ_loc,
     &           CPARTVEC, N, M,
     &           IWRK, IWRKSZ)         
C     Compute sndrcv sizes, store them for later use           
            CALL DMUMPS_NUMVOLSNDRCV(MYID, NUMPROCS, M, RPARTVEC,
     &           NZ_loc, IRN_loc, N, JCN_loc,
     &           IRSNDRCVNUM,IRSNDRCVVOL,
     &           ORSNDRCVNUM,ORSNDRCVVOL,
     &           IWRK,IWRKSZ, 
     &           RSNDRCVSZ(1), RSNDRCVSZ(1+NUMPROCS), COMM)
            CALL DMUMPS_NUMVOLSNDRCV(MYID, NUMPROCS, N, CPARTVEC,
     &           NZ_loc, JCN_loc, M, IRN_loc, 
     &           ICSNDRCVNUM,ICSNDRCVVOL,
     &           OCSNDRCVNUM,OCSNDRCVVOL,
     &           IWRK,IWRKSZ, 
     &           CSNDRCVSZ(1), CSNDRCVSZ(1+NUMPROCS), COMM)         
            CALL DMUMPS_FINDNUMMYROWCOL(MYID, NUMPROCS, COMM,
     &           IRN_loc, JCN_loc, NZ_loc,
     &           RPARTVEC, CPARTVEC, M, N,
     &           INUMMYR,
     &           INUMMYC,     
     &           IWRK, IWRKSZ)
            INTSZR =  IRSNDRCVNUM + ORSNDRCVNUM + 
     &           IRSNDRCVVOL + ORSNDRCVVOL +
     &           2*(NUMPROCS+1) + INUMMYR
            INTSZC = ICSNDRCVNUM + OCSNDRCVNUM + 
     &           ICSNDRCVVOL + OCSNDRCVVOL +
     &           2*(NUMPROCS+1) + INUMMYC
            INTSZ = INTSZR + INTSZC + MAXMN + 
     &           (MPI_STATUS_SIZE +1) *  NUMPROCS
         ELSE
C     NUMPROCS IS 1
             IRSNDRCVNUM = 0
             ORSNDRCVNUM = 0
             IRSNDRCVVOL = 0
             ORSNDRCVVOL = 0
             INUMMYR = 0
             ICSNDRCVNUM  = 0
             OCSNDRCVNUM  = 0
             ICSNDRCVVOL = 0
             OCSNDRCVVOL  = 0
             INUMMYC = 0
             INTSZ = 0
          ENDIF
C     CALCULATE NECESSARY DOUBLE PRECISION SPACE
          RESZR = M + IRSNDRCVVOL + ORSNDRCVVOL
          RESZC = N + ICSNDRCVVOL + OCSNDRCVVOL
          RESZ = RESZR  + RESZC  
C     CALCULATE NECESSARY INT SPACE
C     The last maxmn is tmpwork for setup comm and fillmyrowcol
          REGISTRE(1) = IRSNDRCVNUM 
          REGISTRE(2) = ORSNDRCVNUM
          REGISTRE(3) = IRSNDRCVVOL 
          REGISTRE(4) = ORSNDRCVVOL
          REGISTRE(5) = ICSNDRCVNUM 
          REGISTRE(6) = OCSNDRCVNUM
          REGISTRE(7) = ICSNDRCVVOL
          REGISTRE(8) = OCSNDRCVVOL
          REGISTRE(9) = INUMMYR
          REGISTRE(10) = INUMMYC
          REGISTRE(11) = INTSZ
          REGISTRE(12) = RESZ
       ELSE
C     else of op=1. That is op=2 now.
C     restore the numbers
          IRSNDRCVNUM = REGISTRE(1) 
          ORSNDRCVNUM = REGISTRE(2) 
          IRSNDRCVVOL = REGISTRE(3)
          ORSNDRCVVOL = REGISTRE(4) 
          ICSNDRCVNUM = REGISTRE(5) 
          OCSNDRCVNUM = REGISTRE(6) 
          ICSNDRCVVOL = REGISTRE(7) 
          OCSNDRCVVOL = REGISTRE(8) 
          INUMMYR = REGISTRE(9) 
          INUMMYC = REGISTRE(10)
          IF(NUMPROCS > 1) THEN
C     Check done outsize
C             IF(INTSZ < REGISTRE(11)) THEN ERROR
C             IF(RESZ < REGISTRE(12)) THEN ERROR
C     Fill up myrows and my colsX
             CALL DMUMPS_FILLMYROWCOLINDICES(MYID, NUMPROCS,COMM,    
     &            IRN_loc, JCN_loc, NZ_loc,
     &            RPARTVEC, CPARTVEC, M, N,
     &            IWRK(1), INUMMYR,
     &            IWRK(1+INUMMYR), INUMMYC,     
     &            IWRK(1+INUMMYR+INUMMYC), IWRKSZ-INUMMYR-INUMMYC )
             IMYRPTR = 1
             IMYCPTR = IMYRPTR + INUMMYR
C     Set up comm and run.
C     set pointers in iwrk (4 parts)
C     
C     ROWS    [---------------------------------------------]
             IRNGHBPRCS = IMYCPTR+ INUMMYC
             IRSNDRCVIA = IRNGHBPRCS+IRSNDRCVNUM
             IRSNDRCVJA = IRSNDRCVIA + NUMPROCS+1
             ORNGHBPRCS = IRSNDRCVJA + IRSNDRCVVOL
             ORSNDRCVIA = ORNGHBPRCS + ORSNDRCVNUM
             ORSNDRCVJA = ORSNDRCVIA + NUMPROCS + 1
C     COLS    [---------------------------------------------]
             ICNGHBPRCS = ORSNDRCVJA + ORSNDRCVVOL
             ICSNDRCVIA = ICNGHBPRCS + ICSNDRCVNUM
             ICSNDRCVJA = ICSNDRCVIA + NUMPROCS+1
             OCNGHBPRCS = ICSNDRCVJA + ICSNDRCVVOL
             OCSNDRCVIA = OCNGHBPRCS + OCSNDRCVNUM
             OCSNDRCVJA = OCSNDRCVIA +  NUMPROCS + 1
C     
C     MPI     [-----------------]
             REQUESTS = OCSNDRCVJA + OCSNDRCVVOL
             ISTATUS = REQUESTS + NUMPROCS
C     
C     TMPWRK  [-----------------]
             TMPWORK = ISTATUS + MPI_STATUS_SIZE *  NUMPROCS
             CALL DMUMPS_SETUPCOMMS(MYID, NUMPROCS, M, RPARTVEC,
     &            NZ_loc, IRN_loc,N, JCN_loc,
     &            IRSNDRCVNUM, IRSNDRCVVOL, 
     &            IWRK(IRNGHBPRCS),IWRK(IRSNDRCVIA),IWRK(IRSNDRCVJA),
     &            ORSNDRCVNUM, ORSNDRCVVOL, 
     &            IWRK(ORNGHBPRCS),IWRK(ORSNDRCVIA),IWRK(ORSNDRCVJA),
     &            RSNDRCVSZ(1), RSNDRCVSZ(1+NUMPROCS),
     &            IWRK(TMPWORK), 
     &            IWRK(ISTATUS), IWRK(REQUESTS),
     &            TAG_COMM_ROW, COMM)
             CALL DMUMPS_SETUPCOMMS(MYID, NUMPROCS, N, CPARTVEC,
     &            NZ_loc, JCN_loc, M, IRN_loc,
     &            ICSNDRCVNUM, ICSNDRCVVOL, 
     &            IWRK(ICNGHBPRCS),
     &            IWRK(ICSNDRCVIA),
     &            IWRK(ICSNDRCVJA),
     &            OCSNDRCVNUM, OCSNDRCVVOL, 
     &            IWRK(OCNGHBPRCS),IWRK(OCSNDRCVIA),IWRK(OCSNDRCVJA),
     &            CSNDRCVSZ(1), CSNDRCVSZ(1+NUMPROCS),
     &            IWRK(TMPWORK), 
     &            IWRK(ISTATUS),  IWRK(REQUESTS),
     &            TAG_COMM_COL, COMM)
             CALL DMUMPS_INITREAL(ROWSCA, M, RZERO)
             CALL DMUMPS_INITREAL(COLSCA, N, RZERO)
             CALL DMUMPS_INITREALLST(ROWSCA, M, 
     &            IWRK(IMYRPTR),INUMMYR, RONE)
             CALL DMUMPS_INITREALLST(COLSCA, N, 
     &            IWRK(IMYCPTR),INUMMYC, RONE)   
          ELSE
             CALL DMUMPS_INITREAL(ROWSCA, M, RONE)
             CALL DMUMPS_INITREAL(COLSCA, N, RONE)            
          ENDIF
          ITDRPTR = 1
          ITDCPTR = ITDRPTR + M
C     
          ISRRPTR = ITDCPTR + N
          OSRRPTR = ISRRPTR + IRSNDRCVVOL
C     
          ISRCPTR = OSRRPTR + ORSNDRCVVOL
          OSRCPTR = ISRCPTR + ICSNDRCVVOL
C     To avoid bound check errors...
          IF(NUMPROCS == 1)THEN
             OSRCPTR = OSRCPTR - 1
             ISRCPTR = ISRCPTR - 1
             OSRRPTR = OSRRPTR - 1
             ISRRPTR = ISRRPTR - 1
          ELSE
             IF(IRSNDRCVVOL == 0) ISRRPTR = ISRRPTR - 1
             IF(ORSNDRCVVOL == 0) OSRRPTR = OSRRPTR - 1
             IF(ICSNDRCVVOL == 0) ISRCPTR = ISRCPTR - 1
             IF(OCSNDRCVVOL == 0) OSRCPTR = OSRCPTR - 1
          ENDIF
          ITER = 1
          DO WHILE (ITER.LE.NB1+NB2+NB3)
C     CLEAR temporary Dr and Dc
             IF(NUMPROCS > 1) THEN
                CALL DMUMPS_ZEROOUT(WRKRC(ITDRPTR),M,
     &               IWRK(IMYRPTR),INUMMYR)
                CALL DMUMPS_ZEROOUT(WRKRC(ITDCPTR),N,
     &               IWRK(IMYCPTR),INUMMYC)
             ELSE
                CALL DMUMPS_INITREAL(WRKRC(ITDRPTR),M, RZERO)
                CALL DMUMPS_INITREAL(WRKRC(ITDCPTR),N, RZERO)
             ENDIF
             IF((ITER.LE.NB1).OR.(ITER > NB1+NB2)) THEN
C     INF-NORM ITERATION
                IF((ITER.EQ.1).OR.(OORANGEIND.EQ.1)) THEN
                   DO NZIND=1_8,NZ_loc
                      IR = IRN_loc(NZIND)
                      IC = JCN_loc(NZIND)
                      IF((IR.GE.1).AND.(IR.LE.M).AND.
     &                     (IC.GE.1).AND.(IC.LE.N)) THEN
                         ELM = abs(A_loc(NZIND))*ROWSCA(IR)*COLSCA(IC)
                         IF(WRKRC(ITDRPTR-1+IR)<ELM) THEN
                            WRKRC(ITDRPTR-1+IR)= ELM
                         ENDIF
                         IF(WRKRC(ITDCPTR-1+IC)<ELM) THEN
                            WRKRC(ITDCPTR-1+IC)= ELM
                         ENDIF
                      ELSE
                         OORANGEIND = 1
                      ENDIF
                   ENDDO
                ELSEIF(OORANGEIND.EQ.0) THEN
                   DO NZIND=1,NZ_loc
                      IR = IRN_loc(NZIND)
                      IC = JCN_loc(NZIND)
                      ELM = abs(A_loc(NZIND))*ROWSCA(IR)*COLSCA(IC)
                      IF(WRKRC(ITDRPTR-1+IR)<ELM) THEN
                         WRKRC(ITDRPTR-1+IR)= ELM
                      ENDIF
                      IF(WRKRC(ITDCPTR-1+IC)<ELM) THEN
                         WRKRC(ITDCPTR-1+IC)= ELM
                      ENDIF
                   ENDDO
                ENDIF
                IF(NUMPROCS > 1) THEN
                   CALL DMUMPS_DOCOMMINF(MYID, NUMPROCS,
     &                  WRKRC(ITDCPTR), N, TAG_ITERS+ITER, 
     &                  ICSNDRCVNUM,IWRK(ICNGHBPRCS),
     &                  ICSNDRCVVOL,IWRK(ICSNDRCVIA), IWRK(ICSNDRCVJA), 
     &                  WRKRC(ISRCPTR),
     &                  OCSNDRCVNUM,IWRK(OCNGHBPRCS),
     &                  OCSNDRCVVOL,IWRK(OCSNDRCVIA), IWRK(OCSNDRCVJA),
     &                  WRKRC( OSRCPTR),
     &                  IWRK(ISTATUS),IWRK(REQUESTS),
     &                  COMM)
C     
                  CALL DMUMPS_DOCOMMINF(MYID, NUMPROCS,
     &                  WRKRC(ITDRPTR), M, TAG_ITERS+2+ITER, 
     &                  IRSNDRCVNUM,IWRK(IRNGHBPRCS),
     &                  IRSNDRCVVOL,IWRK(IRSNDRCVIA), IWRK(IRSNDRCVJA), 
     &                  WRKRC(ISRRPTR),
     &                  ORSNDRCVNUM,IWRK(ORNGHBPRCS),
     &                  ORSNDRCVVOL,IWRK(ORSNDRCVIA), IWRK(ORSNDRCVJA),
     &                  WRKRC( OSRRPTR),
     &                  IWRK(ISTATUS),IWRK(REQUESTS),
     &                  COMM)
                  IF((EPS .GT. RZERO) .OR. 
     &                 (ITER.EQ.NB1).OR.
     &                 ((ITER.EQ.NB1+NB2+NB3).AND.
     &                 (NB1+NB3.GT.0))) THEN
                     INFERRROW = DMUMPS_ERRSCALOC(ROWSCA, 
     &                    WRKRC(ITDRPTR), M,
     &                    IWRK(IMYRPTR),INUMMYR)
C     find error for the cols
                     INFERRCOL = DMUMPS_ERRSCALOC(COLSCA,  
     &                    WRKRC(ITDCPTR), N,
     &                    IWRK(IMYCPTR),INUMMYC)
C     get max of those two errors
                     INFERRL = INFERRCOL
                     IF(INFERRROW > INFERRL ) THEN
                        INFERRL = INFERRROW                   
                     ENDIF
C     
                     CALL MPI_ALLREDUCE(INFERRL, INFERRG, 
     &                    1, MPI_DOUBLE_PRECISION,
     &                    MPI_MAX, COMM, IERROR)   
                     IF(INFERRG.LE.EPS) THEN
                        CALL DMUMPS_UPDATESCALE(COLSCA,  
     &                       WRKRC(ITDCPTR),N,
     &                       IWRK(IMYCPTR),INUMMYC)
                        CALL DMUMPS_UPDATESCALE(ROWSCA,  
     &                       WRKRC(ITDRPTR),M,
     &                       IWRK(IMYRPTR),INUMMYR)         
                        IF(ITER .LE. NB1) THEN
                           ITER = NB1+1
                           CYCLE
                        ELSE
                           EXIT
                        ENDIF
                     ENDIF
                  ENDIF                  
               ELSE
C     SINGLE PROCESSOR CASE: INF-NORM ERROR COMPUTATION
                  IF((EPS .GT. RZERO) .OR. 
     &                 (ITER.EQ.NB1).OR.
     &                 ((ITER.EQ.NB1+NB2+NB3).AND.
     &                 (NB1+NB3.GT.0))) THEN
                     INFERRROW = DMUMPS_ERRSCA1(ROWSCA, 
     &                    WRKRC(ITDRPTR), M)
C     find error for the cols
                     INFERRCOL = DMUMPS_ERRSCA1(COLSCA,  
     &                    WRKRC(ITDCPTR), N)
C     get max of those two errors
                     INFERRL = INFERRCOL
                     IF(INFERRROW > INFERRL) THEN
                        INFERRL = INFERRROW                    
                     ENDIF                     
                     INFERRG = INFERRL
                     IF(INFERRG.LE.EPS) THEN
                        CALL DMUMPS_UPSCALE1(COLSCA,  WRKRC(ITDCPTR), N)
                        CALL DMUMPS_UPSCALE1(ROWSCA,  WRKRC(ITDRPTR), M)
                        IF(ITER .LE. NB1) THEN
                           ITER = NB1+1
                           CYCLE
                        ELSE
                           EXIT
                        ENDIF
                     ENDIF 
                  ENDIF
               ENDIF
            ELSE
C     WE HAVE ITER.GT.NB1 AND ITER.LE.NB1+NB2. 
C     ONE-NORM ITERATION
               IF((ITER .EQ.1).OR.(OORANGEIND.EQ.1))THEN
                  DO NZIND=1_8,NZ_loc
                     IR = IRN_loc(NZIND)
                     IC = JCN_loc(NZIND)
                     IF((IR.GE.1).AND.(IR.LE.M).AND.
     &                    (IC.GE.1).AND.(IC.LE.N)) THEN
                        ELM = abs(A_loc(NZIND))*ROWSCA(IR)*COLSCA(IC)
                        WRKRC(ITDRPTR-1+IR) = WRKRC(ITDRPTR-1+IR) + ELM
                        WRKRC(ITDCPTR-1+IC) = WRKRC(ITDCPTR-1+IC) + ELM
                     ELSE
                        OORANGEIND = 1
                     ENDIF
                  ENDDO
               ELSEIF(OORANGEIND.EQ.0) THEN
                  DO NZIND=1,NZ_loc
                     IR = IRN_loc(NZIND)
                     IC = JCN_loc(NZIND)
                     ELM = abs(A_loc(NZIND))*ROWSCA(IR)*COLSCA(IC)
                     WRKRC(ITDRPTR-1+IR) = WRKRC(ITDRPTR-1+IR) + ELM
                     WRKRC(ITDCPTR-1+IC) = WRKRC(ITDCPTR-1+IC) + ELM
                  ENDDO
               ENDIF
               IF(NUMPROCS > 1) THEN                 
                  CALL DMUMPS_DOCOMM1N(MYID, NUMPROCS,
     &                 WRKRC(ITDCPTR), N, TAG_ITERS+ITER, 
     &                 ICSNDRCVNUM, IWRK(ICNGHBPRCS),
     &                 ICSNDRCVVOL, IWRK(ICSNDRCVIA), IWRK(ICSNDRCVJA), 
     &                 WRKRC(ISRCPTR),
     &                 OCSNDRCVNUM, IWRK(OCNGHBPRCS),
     &                 OCSNDRCVVOL, IWRK(OCSNDRCVIA), IWRK(OCSNDRCVJA),
     &                 WRKRC( OSRCPTR),
     &                 IWRK(ISTATUS), IWRK(REQUESTS),
     &                 COMM)
C     
                  CALL DMUMPS_DOCOMM1N(MYID, NUMPROCS,
     &                 WRKRC(ITDRPTR), M, TAG_ITERS+2+ITER, 
     &                 IRSNDRCVNUM, IWRK(IRNGHBPRCS),
     &                 IRSNDRCVVOL, IWRK(IRSNDRCVIA), IWRK(IRSNDRCVJA), 
     &                 WRKRC(ISRRPTR),
     &                 ORSNDRCVNUM, IWRK(ORNGHBPRCS),
     &                 ORSNDRCVVOL, IWRK(ORSNDRCVIA), IWRK(ORSNDRCVJA),
     &                 WRKRC( OSRRPTR),
     &                 IWRK(ISTATUS), IWRK(REQUESTS),
     &                 COMM)
                  IF((EPS .GT. RZERO) .OR. 
     &                 ((ITER.EQ.NB1+NB2).AND.
     &                 (NB2.GT.0))) THEN
                     ONEERRROW = DMUMPS_ERRSCALOC(ROWSCA, 
     &                    WRKRC(ITDRPTR), M,
     &                    IWRK(IMYRPTR),INUMMYR)
C     find error for the cols
                     ONEERRCOL = DMUMPS_ERRSCALOC(COLSCA,  
     &                    WRKRC(ITDCPTR), N,
     &                    IWRK(IMYCPTR),INUMMYC)
C     get max of those two errors
                     ONEERRL = ONEERRCOL
                     IF(ONEERRROW > ONEERRL ) THEN
                        ONEERRL = ONEERRROW                   
                     ENDIF
C     
                     CALL MPI_ALLREDUCE(ONEERRL, ONEERRG, 
     &                    1, MPI_DOUBLE_PRECISION,
     &                    MPI_MAX, COMM, IERROR)   
                     IF(ONEERRG.LE.EPS) THEN
                        CALL DMUMPS_UPDATESCALE(COLSCA,
     &                       WRKRC(ITDCPTR),N,
     &                       IWRK(IMYCPTR),INUMMYC)
                        CALL DMUMPS_UPDATESCALE(ROWSCA,
     &                       WRKRC(ITDRPTR),M,
     &                       IWRK(IMYRPTR),INUMMYR)          
                        ITER = NB1+NB2+1
                        CYCLE
                     ENDIF
                  ENDIF                            
               ELSE
C     SINGLE-PROCESSOR CASE: ONE-NORM ERROR COMPUTATION
                  IF((EPS .GT. RZERO) .OR. 
     &                 ((ITER.EQ.NB1+NB2).AND.
     &                 (NB2.GT.0))) THEN
                     ONEERRROW = DMUMPS_ERRSCA1(ROWSCA, 
     &                    WRKRC(ITDRPTR), M)
C     find error for the cols
                     ONEERRCOL = DMUMPS_ERRSCA1(COLSCA,  
     &                    WRKRC(ITDCPTR), N)
C     get max of those two errors
                     ONEERRL = ONEERRCOL
                     IF(ONEERRROW > ONEERRL) THEN
                        ONEERRL = ONEERRROW                    
                     ENDIF                     
                     ONEERRG = ONEERRL
                     IF(ONEERRG.LE.EPS) THEN
                        CALL DMUMPS_UPSCALE1(COLSCA,  WRKRC(ITDCPTR), N)
                        CALL DMUMPS_UPSCALE1(ROWSCA,  WRKRC(ITDRPTR), M)
                        ITER = NB1+NB2+1                        
                        CYCLE
                     ENDIF
                  ENDIF                  
               ENDIF 
            ENDIF
            IF(NUMPROCS > 1) THEN               
               CALL DMUMPS_UPDATESCALE(COLSCA,  WRKRC(ITDCPTR), N,
     &              IWRK(IMYCPTR),INUMMYC)
               CALL DMUMPS_UPDATESCALE(ROWSCA,  WRKRC(ITDRPTR), M,
     &              IWRK(IMYRPTR),INUMMYR)               
C   
            ELSE
C     SINGLE PROCESSOR CASE: Conv check and update of sca arrays
               CALL DMUMPS_UPSCALE1(COLSCA,  WRKRC(ITDCPTR), N)
               CALL DMUMPS_UPSCALE1(ROWSCA,  WRKRC(ITDRPTR), M)
            ENDIF
            ITER = ITER + 1
         ENDDO
         ONENORMERR = ONEERRG 
         INFNORMERR = INFERRG 
         IF(NUMPROCS > 1) THEN
            CALL MPI_REDUCE(ROWSCA, WRKRC(1), M, MPI_DOUBLE_PRECISION,
     &           MPI_MAX, 0, 
     &           COMM, IERROR)
            IF(MYID.EQ.0) THEN
               DO I=1, M
                  ROWSCA(I) = WRKRC(I)
               ENDDO
            ENDIF
C     Scaling factors are printed
C     WRITE (6,*) MYID, 'ROWSCA=',ROWSCA
C     WRITE (6,*) MYID, 'COLSCA=',COLSCA
C     CALL FLUSH(6)
c     REduce the whole scaling factors to processor 0 of COMM
            CALL MPI_REDUCE(COLSCA, WRKRC(1+M), N, MPI_DOUBLE_PRECISION,
     &           MPI_MAX, 0, 
     &           COMM, IERROR)
            If(MYID.EQ.0) THEN
               DO I=1, N
                  COLSCA(I) = WRKRC(I+M)
               ENDDO
            ENDIF         
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_SIMSCALEABSUNS
C
C 
C     SEPARATOR: Another function begins
C
C 
      SUBROUTINE DMUMPS_SIMSCALEABSSYM(IRN_loc, JCN_loc, A_loc, NZ_loc,
     &     N, NUMPROCS, MYID, COMM,
     &     PARTVEC, 
     &     RSNDRCVSZ, 
     &     REGISTRE,
     &     IWRK, IWRKSZ,
     &     INTSZ, RESZ, OP,
     &     SCA, WRKRC, ISZWRKRC,
     &     NB1, NB2, NB3, EPS,
     &     ONENORMERR, INFNORMERR)    
C----------------------------------------------------------------------
C Input parameters:
C     N: size of matrix (sym matrix, square).
C     NUMPROCS, MYID, COMM: guess what are those
C     PARTVEC:  row/col partvec to be filled when OP=1
C     RSNDRCVSZ:send recv sizes for row/col operations. 
C               to be filled when OP=1
C     REGISTRE: to store some pointers (size etc). Its size is 12,
C               but we do not use all in this routine.
C     IWRK: working space. when OP=1 IWRKSZ.GE.2*MAXMN
C           when OP=2 INTSZ portion is used. Donc, IWRKSZ>INTSZ 
C           when OP=2
C      IWRKSZ: size
C      INTSZ: to be computed when OP=1, necessary integer space to run 
C             scaling algo when OP=2
C      RESZ:  to be computed when OP=1, necessary real space to run 
C             scaling algo when OP=2
C      OP: 
C          =1 estimation of memory and construction of partvecs
C           writes into PARTVEC,RSNDRCVSZ,REGISTRE
C           does not access WRKRC, uses IWRK as workspace
C           computes INTSZ and RESZ.
C          =2 Compute scalings 
C           restores pointers from REGISTRE, 
C           stores communication structure in IWRK (from the start). 
C
C      SCA: space for row/col scaling factor; has size M
C      WRKRC: real working space. when OP=1, is not accessed. Donc, it
C             can be declared to be of size 1 at OP=1 call.
C      ISZWRKRC: size
C      SYM: is matrix symmetric
C      NB1, NB2, NB3: algo runs 
C                     NB1 iters of inf-norm (default  1/1), 
C                     NB2 iters of 1-norm   (default  3/10),
C                     NB3 iters of inf-norm (default  3/10).
C          in succession.
C      EPS: tolerance for concergence. 
C           IF EPS < 0.R0 then does not test convergence.
C           See comments for the uns case above.      
C     ONENORMERR : error in one norm scaling (see comments for the 
C                  uns case above), 
C     INFNORMERR : error in inf norm scaling (see comments for the 
C                  uns case above).
C---------------------------------------------------------------------
C On input:
C      OP=1==>Requirements
C             IWRKSZ.GE.2*MAXMN   XXXX compare with uns variant.
C             PARTVEC  of size N
C             SNDRCVSZ of size 2*NUMPROCS
C             REGISTRE  of size 12
C             
C      OP=2==>Requirements
C             INTSZ .GE. REGISTRE(11)
C             RESZ  .GE. REGISTRE(12) 
C---------------------------------------------------------------------
C On output:
C     SCA
C            at processor 0 of COMM: complete factors.
C            at other processors   : only the SCA(i) and SCA(j)
C            for which there is a nonzero a_ij.
C     ONENORMERR : error in one norm scaling 
C                = -1.0 if iter2=0.
C     INFNORMERR : error in inf norm scaling 
C                = inf norm error at iter3 if iter3 > 0
C                = inf norm error at iter1 if iter1 > 0, iter3=0
C                = -1.0 if iter1=iter3=0
C ---------------------------------------------------------------------
C NOTE: some variables are named in such a way that they correspond
C       to the row variables in unsym case. They are used for both 
C       row and col communications.
C ---------------------------------------------------------------------
C References:
C     The scaling algorithms are based on those discussed in
C     [1] D. Ruiz, "A scaling algorithm to equilibrate both rows and 
C         columns norms in matrices", Tech. Rep. Rutherford 
C         Appleton Laboratory, Oxon, UK and ENSEEIHT-IRIT, 
C         Toulouse, France, RAL-TR-2001-034 and RT/APO/01/4, 2001.
C     [2] D. Ruiz and B. Ucar, "A symmetry preserving algorithm for
C         matrix scaling", in preparation as of Jan'08.
C
C     The parallelization approach is based on discussion in
C     [3] P. R. Amestoy, I. S. Duff, D. Ruiz, and B. Ucar, "A parallel
C         matrix scaling algorithm", accepted for publication, 
C         In proceedings of VECPAR'08-International Meeting-High 
C         Performance Computing for Computational Science, Jan'08.
C     and was supported by ANR-SOLSTICE project (ANR-06-CIS6-010)
C ---------------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'mpif.h'
      INTEGER(8) :: NZ_loc 
      INTEGER N, IWRKSZ, OP
      INTEGER NUMPROCS, MYID, COMM
      INTEGER INTSZ, RESZ
      INTEGER IRN_loc(NZ_loc)
      INTEGER JCN_loc(NZ_loc)
      DOUBLE PRECISION A_loc(NZ_loc)
      INTEGER PARTVEC(N), RSNDRCVSZ(2*NUMPROCS)
      INTEGER IWRK(IWRKSZ)
      INTEGER REGISTRE(12)
      DOUBLE PRECISION SCA(N)
      INTEGER ISZWRKRC
      DOUBLE PRECISION WRKRC(ISZWRKRC)
C     LOCALS
      INTEGER IRSNDRCVNUM, ORSNDRCVNUM
      INTEGER IRSNDRCVVOL, ORSNDRCVVOL
      INTEGER  INUMMYR
C IMPORTANT POINTERS
      INTEGER IMYRPTR,IMYCPTR 
      INTEGER IRNGHBPRCS, IRSNDRCVIA,IRSNDRCVJA
      INTEGER ORNGHBPRCS, ORSNDRCVIA,ORSNDRCVJA
      INTEGER ISTATUS, REQUESTS, TMPWORK
      INTEGER ITDRPTR, ISRRPTR, OSRRPTR
      DOUBLE PRECISION ONENORMERR,INFNORMERR
C     FOR the scaling phase  
      INTEGER NB1, NB2, NB3
      DOUBLE PRECISION EPS
C     Iteration vars 
      INTEGER ITER, IR, IC
      INTEGER(8) :: NZIND
      DOUBLE PRECISION ELM
C     COMM TAGS....
      INTEGER TAG_COMM_ROW
      PARAMETER(TAG_COMM_ROW=101)
      INTEGER TAG_ITERS
      PARAMETER(TAG_ITERS=102)
C     FUNCTIONS
      EXTERNAL DMUMPS_CREATEPARTVECSYM,
     &     DMUMPS_NUMVOLSNDRCVSYM, 
     &     DMUMPS_SETUPCOMMSSYM,
     &     DMUMPS_FINDNUMMYROWCOLSYM, 
     &     DMUMPS_CHKCONVGLOSYM,
     &     DMUMPS_CHK1CONV,
     &     DMUMPS_FILLMYROWCOLINDICESSYM,
     &     DMUMPS_DOCOMMINF,
     &     DMUMPS_DOCOMM1N,
     &     DMUMPS_INITREAL,
     &     DMUMPS_INITREALLST
      INTEGER DMUMPS_CHKCONVGLOSYM 
      INTEGER DMUMPS_CHK1CONV
      DOUBLE PRECISION DMUMPS_ERRSCALOC
      DOUBLE PRECISION DMUMPS_ERRSCA1
      INTRINSIC abs
      DOUBLE PRECISION RONE, RZERO
      PARAMETER(RONE=1.0D0,RZERO=0.0D0)
C     TMP VARS
      INTEGER INTSZR
      INTEGER MAXMN
      INTEGER I, IERROR
      DOUBLE PRECISION ONEERRL, ONEERRG
      DOUBLE PRECISION INFERRL, INFERRG
      INTEGER OORANGEIND
      OORANGEIND = 0
      INFERRG = -RONE
      ONEERRG = -RONE
      MAXMN = N
      IF(OP == 1) THEN
         IF(NUMPROCS > 1) THEN
C     Check done outside
C     IF(IWRKSZ.LT.2*MAXMN) THEN   ERROR.... 
            CALL DMUMPS_CREATEPARTVECSYM(MYID, NUMPROCS, COMM,
     &           IRN_loc, JCN_loc, NZ_loc,
     &           PARTVEC, N,
     &           IWRK, IWRKSZ)
C     
            CALL DMUMPS_NUMVOLSNDRCVSYM(MYID, NUMPROCS, N, PARTVEC,
     &           NZ_loc, IRN_loc, JCN_loc, IRSNDRCVNUM,IRSNDRCVVOL,
     &           ORSNDRCVNUM, ORSNDRCVVOL,
     &           IWRK,IWRKSZ, 
     &           RSNDRCVSZ(1), RSNDRCVSZ(1+NUMPROCS), COMM)
C     
            CALL DMUMPS_FINDNUMMYROWCOLSYM(MYID, NUMPROCS, COMM,
     &           IRN_loc, JCN_loc, NZ_loc,
     &           PARTVEC, N,
     &           INUMMYR,
     &           IWRK, IWRKSZ)
C     
            INTSZR =  IRSNDRCVNUM + ORSNDRCVNUM + 
     &           IRSNDRCVVOL + ORSNDRCVVOL +
     &           2*(NUMPROCS+1) + INUMMYR
            INTSZ = INTSZR + N + 
     &           (MPI_STATUS_SIZE +1) *  NUMPROCS
         ELSE
C     NUMPROCS IS 1
            IRSNDRCVNUM = 0
            ORSNDRCVNUM = 0
            IRSNDRCVVOL = 0 
            ORSNDRCVVOL = 0
            INUMMYR = 0
            INTSZ = 0
         ENDIF
C     CALCULATE NECESSARY DOUBLE PRECISION SPACE
         RESZ = N + IRSNDRCVVOL + ORSNDRCVVOL
         REGISTRE(1) = IRSNDRCVNUM 
         REGISTRE(2) = ORSNDRCVNUM
         REGISTRE(3) = IRSNDRCVVOL 
         REGISTRE(4) = ORSNDRCVVOL
         REGISTRE(9) = INUMMYR
         REGISTRE(11) = INTSZ
         REGISTRE(12) = RESZ
      ELSE
C     else of op=1. That is op=2 now.
C     restore the numbers
         IRSNDRCVNUM = REGISTRE(1) 
         ORSNDRCVNUM = REGISTRE(2) 
         IRSNDRCVVOL = REGISTRE(3)
         ORSNDRCVVOL = REGISTRE(4) 
         INUMMYR = REGISTRE(9) 
          IF(NUMPROCS > 1) THEN
C     Check done outsize
C             IF(INTSZ < REGISTRE(11)) THEN ERROR
C             IF(RESZ < REGISTRE(12)) THEN ERROR
C     Fill up myrows and my colsX
             CALL DMUMPS_FILLMYROWCOLINDICESSYM(MYID, NUMPROCS,COMM,    
     &            IRN_loc, JCN_loc, NZ_loc,
     &            PARTVEC, N,
     &            IWRK(1), INUMMYR,
     &            IWRK(1+INUMMYR), IWRKSZ-INUMMYR)
             IMYRPTR = 1
             IMYCPTR = IMYRPTR + INUMMYR
C     Set up comm and run.
C     set pointers in iwrk (3 parts)
C     
C     ROWS    [---------------------------------------------]
             IRNGHBPRCS = IMYCPTR 
             IRSNDRCVIA = IRNGHBPRCS+IRSNDRCVNUM
             IRSNDRCVJA = IRSNDRCVIA + NUMPROCS+1
             ORNGHBPRCS = IRSNDRCVJA + IRSNDRCVVOL
             ORSNDRCVIA = ORNGHBPRCS + ORSNDRCVNUM
             ORSNDRCVJA = ORSNDRCVIA + NUMPROCS + 1
C     MPI     [-----------------]
             REQUESTS = ORSNDRCVJA + ORSNDRCVVOL 
             ISTATUS = REQUESTS + NUMPROCS
C     TMPWRK  [-----------------]
             TMPWORK = ISTATUS + MPI_STATUS_SIZE *  NUMPROCS
             CALL DMUMPS_SETUPCOMMSSYM(MYID, NUMPROCS, N, PARTVEC,
     &            NZ_loc, IRN_loc, JCN_loc,
     &            IRSNDRCVNUM, IRSNDRCVVOL, 
     &            IWRK(IRNGHBPRCS),IWRK(IRSNDRCVIA),IWRK(IRSNDRCVJA),
     &            ORSNDRCVNUM, ORSNDRCVVOL, 
     &            IWRK(ORNGHBPRCS),IWRK(ORSNDRCVIA),IWRK(ORSNDRCVJA),
     &            RSNDRCVSZ(1), RSNDRCVSZ(1+NUMPROCS),
     &            IWRK(TMPWORK), 
     &            IWRK(ISTATUS), IWRK(REQUESTS),
     &            TAG_COMM_ROW, COMM)
             CALL DMUMPS_INITREAL(SCA, N, RZERO)
             CALL DMUMPS_INITREALLST(SCA, N, 
     &            IWRK(IMYRPTR),INUMMYR, RONE)
          ELSE
             CALL DMUMPS_INITREAL(SCA, N, RONE)
          ENDIF
          ITDRPTR = 1
          ISRRPTR = ITDRPTR + N
          OSRRPTR = ISRRPTR + IRSNDRCVVOL
C     
C     To avoid bound check errors...
          IF(NUMPROCS == 1)THEN
             OSRRPTR = OSRRPTR - 1
             ISRRPTR = ISRRPTR - 1
          ELSE
             IF(IRSNDRCVVOL == 0) ISRRPTR = ISRRPTR - 1
             IF(ORSNDRCVVOL == 0) OSRRPTR = OSRRPTR - 1
          ENDIF
C     computation starts
          ITER = 1
          DO WHILE(ITER.LE.NB1+NB2+NB3)
C     CLEAR temporary Dr and Dc
             IF(NUMPROCS > 1) THEN
                CALL DMUMPS_ZEROOUT(WRKRC(ITDRPTR),N,
     &               IWRK(IMYRPTR),INUMMYR)
             ELSE
                CALL DMUMPS_INITREAL(WRKRC(ITDRPTR),N, RZERO)
             ENDIF
             IF((ITER.LE.NB1).OR.(ITER > NB1+NB2)) THEN
C     INF-NORM ITERATION
                IF((ITER.EQ.1).OR.(OORANGEIND.EQ.1)) THEN
                   DO NZIND=1_8,NZ_loc
                      IR = IRN_loc(NZIND)
                      IC = JCN_loc(NZIND)
                      IF((IR.GE.1).AND.(IR.LE.N).AND.
     &                     (IC.GE.1).AND.(IC.LE.N)) THEN
                         ELM = abs(A_loc(NZIND))*SCA(IR)*SCA(IC)
                         IF(WRKRC(ITDRPTR-1+IR)<ELM) THEN
                            WRKRC(ITDRPTR-1+IR)= ELM
                         ENDIF
                         IF(WRKRC(ITDRPTR-1+IC)<ELM) THEN
                            WRKRC(ITDRPTR-1+IC)= ELM
                         ENDIF
                      ELSE
                         OORANGEIND = 1
                      ENDIF
                   ENDDO
                ELSEIF(OORANGEIND.EQ.0) THEN
                   DO NZIND=1_8,NZ_loc
                      IR = IRN_loc(NZIND)
                      IC = JCN_loc(NZIND)
                      ELM = abs(A_loc(NZIND))*SCA(IR)*SCA(IC)
                      IF(WRKRC(ITDRPTR-1+IR)<ELM) THEN
                         WRKRC(ITDRPTR-1+IR)= ELM
                      ENDIF
                      IF(WRKRC(ITDRPTR-1+IC)<ELM) THEN
                         WRKRC(ITDRPTR-1+IC)= ELM
                      ENDIF
                   ENDDO
                ENDIF                      
                IF(NUMPROCS > 1) THEN
                  CALL DMUMPS_DOCOMMINF(MYID, NUMPROCS,
     &                  WRKRC(ITDRPTR), N, TAG_ITERS+2+ITER, 
     &                  IRSNDRCVNUM,IWRK(IRNGHBPRCS),
     &                  IRSNDRCVVOL,IWRK(IRSNDRCVIA), IWRK(IRSNDRCVJA), 
     &                  WRKRC(ISRRPTR),
     &                  ORSNDRCVNUM,IWRK(ORNGHBPRCS),
     &                  ORSNDRCVVOL,IWRK(ORSNDRCVIA), IWRK(ORSNDRCVJA),
     &                  WRKRC( OSRRPTR),
     &                  IWRK(ISTATUS),IWRK(REQUESTS),
     &                  COMM)
                  IF((EPS .GT. RZERO) .OR. 
     &                 (ITER.EQ.NB1).OR.
     &                 ((ITER.EQ.NB1+NB2+NB3).AND.
     &                 (NB1+NB3.GT.0))) THEN
                     INFERRL = DMUMPS_ERRSCALOC(SCA,  
     &                    WRKRC(ITDRPTR), N,
     &                    IWRK(IMYRPTR),INUMMYR)                  
                     CALL MPI_ALLREDUCE(INFERRL, INFERRG, 
     &                    1, MPI_DOUBLE_PRECISION,
     &                    MPI_MAX, COMM, IERROR)   
                     IF(INFERRG.LE.EPS) THEN
                        CALL DMUMPS_UPDATESCALE(SCA,  WRKRC(ITDRPTR), N,
     &                       IWRK(IMYRPTR),INUMMYR)
                        IF(ITER .LE. NB1) THEN
                           ITER = NB1+1
                           CYCLE
                        ELSE
                           EXIT
                        ENDIF
                     ENDIF
                  ENDIF
               ELSE
C     SINGLE PROCESSOR CASE: INF-NORM ERROR COMPUTATION
                  IF((EPS .GT. RZERO) .OR. 
     &                 (ITER.EQ.NB1).OR.
     &                 ((ITER.EQ.NB1+NB2+NB3).AND.
     &                 (NB1+NB3.GT.0))) THEN
                     INFERRL = DMUMPS_ERRSCA1(SCA, 
     &                    WRKRC(ITDRPTR), N)
                     INFERRG = INFERRL
                     IF(INFERRG.LE.EPS) THEN
                        CALL DMUMPS_UPSCALE1(SCA,  WRKRC(ITDRPTR), N)
                        IF(ITER .LE. NB1) THEN
                           ITER = NB1+1
                           CYCLE
                        ELSE
                           EXIT
                        ENDIF
                     ENDIF 
                  ENDIF
               ENDIF
            ELSE
C     WE HAVE ITER.GT.NB1 AND ITER.LE.NB1+NB2. 
C     ONE-NORM ITERATION
               IF((ITER.EQ.1).OR.(OORANGEIND.EQ.1))THEN
                  DO NZIND=1_8,NZ_loc
                     IR = IRN_loc(NZIND)
                     IC = JCN_loc(NZIND)
                     IF((IR.GE.1).AND.(IR.LE.N).AND.
     &                    (IC.GE.1).AND.(IC.LE.N)) THEN
                        ELM = abs(A_loc(NZIND))*SCA(IR)*SCA(IC)
                        WRKRC(ITDRPTR-1+IR) = WRKRC(ITDRPTR-1+IR) + ELM
                        IF(IR.NE.IC) THEN
                           WRKRC(ITDRPTR-1+IC) = 
     &                          WRKRC(ITDRPTR-1+IC) + ELM
                        ENDIF
                     ELSE
                        OORANGEIND = 1
                     ENDIF
                  ENDDO
               ELSEIF(OORANGEIND.EQ.0)THEN
                  DO NZIND=1_8,NZ_loc
                     IR = IRN_loc(NZIND)
                     IC = JCN_loc(NZIND)
                     ELM = abs(A_loc(NZIND))*SCA(IR)*SCA(IC)
                     WRKRC(ITDRPTR-1+IR) = WRKRC(ITDRPTR-1+IR) + ELM
                     IF(IR.NE.IC) THEN
                        WRKRC(ITDRPTR-1+IC) = WRKRC(ITDRPTR-1+IC) + ELM
                     ENDIF
                  ENDDO
               ENDIF
               IF(NUMPROCS > 1) THEN
                  CALL DMUMPS_DOCOMM1N(MYID, NUMPROCS,
     &                 WRKRC(ITDRPTR), N, TAG_ITERS+2+ITER, 
     &                 IRSNDRCVNUM, IWRK(IRNGHBPRCS),
     &                 IRSNDRCVVOL, IWRK(IRSNDRCVIA), IWRK(IRSNDRCVJA), 
     &                 WRKRC(ISRRPTR),
     &                 ORSNDRCVNUM, IWRK(ORNGHBPRCS),
     &                 ORSNDRCVVOL, IWRK(ORSNDRCVIA), IWRK(ORSNDRCVJA),
     &                 WRKRC( OSRRPTR),
     &                 IWRK(ISTATUS), IWRK(REQUESTS),
     &                 COMM)
                  IF((EPS .GT. RZERO) .OR. 
     &                 ((ITER.EQ.NB1+NB2).AND.
     &                 (NB2.GT.0))) THEN
                     ONEERRL = DMUMPS_ERRSCALOC(SCA,  
     &                    WRKRC(ITDRPTR), N,
     &                    IWRK(IMYRPTR),INUMMYR) 
C     mpi allreduce.
                     CALL MPI_ALLREDUCE(ONEERRL, ONEERRG, 
     &                    1, MPI_DOUBLE_PRECISION,
     &                    MPI_MAX, COMM, IERROR)
                     IF(ONEERRG.LE.EPS) THEN
                        CALL DMUMPS_UPDATESCALE(SCA,  WRKRC(ITDRPTR), N,
     &                       IWRK(IMYRPTR),INUMMYR)
                        ITER = NB1+NB2+1
                        CYCLE
                     ENDIF
                  ENDIF
               ELSE
C     SINGLE-PROCESSOR CASE: ONE-NORM ERROR COMPUTATION
                  IF((EPS .GT. RZERO) .OR. 
     &                 ((ITER.EQ.NB1+NB2).AND.
     &                 (NB2.GT.0))) THEN
                     ONEERRL = DMUMPS_ERRSCA1(SCA, 
     &                    WRKRC(ITDRPTR), N)
                     ONEERRG = ONEERRL
                     IF(ONEERRG.LE.EPS) THEN
                        CALL DMUMPS_UPSCALE1(SCA,  WRKRC(ITDRPTR), N)
                        ITER = NB1+NB2+1
                        CYCLE
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
            IF(NUMPROCS > 1) THEN
               CALL DMUMPS_UPDATESCALE(SCA,  WRKRC(ITDRPTR), N,
     &              IWRK(IMYRPTR),INUMMYR)
            ELSE
               CALL DMUMPS_UPSCALE1(SCA,  WRKRC(ITDRPTR), N)
            ENDIF     
            ITER = ITER + 1
         ENDDO
         ONENORMERR = ONEERRG 
         INFNORMERR = INFERRG 
         IF(NUMPROCS > 1) THEN
            CALL MPI_REDUCE(SCA, WRKRC(1), N, MPI_DOUBLE_PRECISION,
     &           MPI_MAX, 0, 
     &           COMM, IERROR)
            IF(MYID.EQ.0) THEN
               DO I=1, N
                  SCA(I) = WRKRC(I)
               ENDDO
            ENDIF
         ENDIF
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_SIMSCALEABSSYM
