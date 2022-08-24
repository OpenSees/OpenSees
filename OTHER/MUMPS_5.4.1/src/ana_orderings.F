C =========================================================
C 
C  This file includes various modifications of an original 
C  routine MUMPS_ANA_H. The main reference for the approach 
C  used in this routine is
C   Patrick Amestoy, Timothy A. Davis, and Iain S. Duff,
C    "An approximate minimum degree ordering algorithm,"
C    SIAM J. Matrix Analysis  vol 17, pages=886--905 (1996)
C    MUMPS_ANA_H is based on the original AMD code:
C
C    AMD, Copyright (c), 1996-2016, Timothy A. Davis,
C    Patrick R. Amestoy, and Iain S. Duff.  All Rights Reserved.
C    Used in MUMPS under the BSD 3-clause license.
C
C All other routines are modifications of this original routine
C done by MUMPS developers over the years (1996-2020) and are 
C used in MUMPS under the BSD 3-clause license.
C
C BSD 3-clause licence:
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions 
C are met:
C  * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C  * Redistributions in binary form must reproduce the above 
C    copyright notice, this list of conditions and the following 
C    disclaimer in the documentation and/or other materials provided 
C    with the distribution.
C  * Neither the name of the University of California, Berkeley nor 
C    the names of its contributors may be used to endorse or promote 
C    products derived from this software without specific prior 
C    written permission.
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND 
C    CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, 
C    INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
C    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
C    DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR 
C    CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
C    NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
C    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
C    HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
C    CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
C    OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
C    EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C
C   MUMPS_AMD_ELT is a modification 
C   designed to handle amalgamated and compressed
C   graphs and was developed in 1999 by Patrick Amestoy 
C   in the context of the PARASOL project (1997-1999).
C 
C   MUMPS_HAMD is a modification 
C   designed to take into account a halo in the graph. 
C   The graph is composed is partitioned in two types of nodes
C   the so called internal nodes and the so called halo nodes.
C   Halo nodes cannot be selected the both the initial degrees 
C   and updated degrees of internal node should be taken 
C   into account.  
C   This routine also referred to as HALOAMD in MUMPS comments
C   is used for both Schur functionality and in the coupling with 
C   partitioners such as SCOTCH.
C   This code was developed for MUMPS platform 
C   by Patrick Amestoy between 1997 and 1999.
C
C   MUMPS_HAMF4 is a major modification of MUMPS_HAMD 
C   since metric used to select pivots in not anymore the 
C   degree but an approximation of the fill-in.
C   In this approximation 
C   all cliques of elements adjacent to the variable are deducted.
C   Written by Patrick Amestoy between 1999 and 2000.
C   It is also used by F. Pellegrini in SCOTCH since 2000.
C
C   MUMPS_QAMD: modified version of reference AMD routine MUMPS_ANA_H 
C   designed to automatically detect and exploit dense or quasi dense
C   rows in the reduced matrix at any step of the minimum degree.
C   Written in 1997 by Patrick Amestoy.
C   References:
C    P.R. AMESTOY, Recent progress in parallel multifrontal solvers
C      for unsymmetric sparse matrices,
C      Proceedings of the 15th World Congress on Scientific Computation,
C      Modelling and Applied Mathematics, IMACS, Berlin (1997).
C    P.R. AMESTOY (1999), Methodes directes paralleles de
C      resolution des systemes creux de grande taille.
C      Rapport de these d'habilitation de l'INPT.
C
C   MUMPS_CST_AMF: modified version of MUMPS_HAMF4 routine 
C   implementing constraint minimum fill-in based ordering.
C   Written by Stephane Pralet for MUMPS platform 
C   during his post-doctorate at INPT-IRIT (Oct. 2004- Oct. 2005)
C
C  ----------------------------------------
C  To suppress aggressive absorption in ...
C      MUMPS_ANA_H   : Historical AMD
C        define NOAGG1
C      MUMPS_AMD_ELT  : (work on compressed graphs)
C        define NOAGG2
C      MUMPS_HAMD     : AMD with Halo and used for Schur
C        define NOAGG3
C      MUMPS_HAMF4 : Halo AMF version
C        define NOAGG4
C      MUMPS_QAMD     : Quasi dense
C        define NOAGG5
C      MUMPS_SYMQAMD  : Symbolic facto based on quasi dense
C        In the case of MUMPS_SYMQAMD, the aggressive absorption
C        is controlled by a parameter, AGG6.
C
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C MUMPS_ANA_H:  Approximate Minimum Degree AMD approach.
C
C Description of MUMPS_ANA_H
C   Given a representation of the nonzero pattern of a symmetric matrix,
C   A, (excluding the diagonal) perform an approximate minimum
C   degree ordering to compute a pivot order
C   such that fill-in in the Cholesky factors A = LL^T is kept low. 
C   Aggressive absorption might be used to
C   tighten the bound on the degree.  This can result a
C   significant improvement in the quality of the ordering for
C   some matrices.
C
C     References and definitions:  
C     [1] Timothy A. Davis and Iain Duff, "An unsymmetric-pattern
C              multifrontal method for sparse LU factorization",
C              SIAM J. Matrix Analysis and Applications, 
C              volume=18, pages=140-158 (1997)
C     [2] Patrick R. Amestoy, Timothy A. Davis, and Iain S. Duff,
C              "An approximate minimum degree ordering algorithm,"
C              SIAM J. Matrix Analysis  vol 17, pages=886--905 (1996)
C     [3] Alan George and Joseph Liu, "The evolution of the
C              minimum degree ordering algorithm," SIAM Review, vol.
C              31, no. 1, pp. 1-19, March 1989.  We list below the
C              features mentioned in that paper that this code
C              includes:
C       mass elimination:
C               Yes.  supervariable detection for mass elimination.
C       indistinguishable nodes:
C               Yes (we call these "supervariables").  
C               We modified the approach used by Duff and Reid to 
C               detect them (the previous hash was the true degree,
C               which we no longer keep track of).  A supervariable is
C               a set of rows with identical nonzero pattern.  All
C               variables in a supervariable are eliminated together.
C               Each supervariable has as its numerical name that of
C               one of its variables (its principal variable).
C       quotient graph representation:
C               Yes.  We use the term "element" for the cliques formed
C               during elimination.  
C               The algorithm can operate in place, but it will work
C               more efficiently if given some "elbow room."
C       element absorption:
C               Yes. Similar to Duff,Reid and  George,Liu approaches
C       external degree:
C               Yes. Similar to Duff, Reid and  George, Liu approaches
C       incomplete degree update and multiple elimination:
C               No implemented.  Our method of
C               degree update within MUMPS_ANA_H is element-based, not
C               variable-based.  It is thus not well-suited for use
C               with incomplete degree update or multiple elimination.
C
C-----------------------------------------------------------------------
      SUBROUTINE MUMPS_ANA_H(TOTEL, COMPUTE_PERM,
     &                   N, IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
     &                   LAST, NCMPA, DEGREE, HEAD, NEXT, W, PARENT)
C 
C    Restrictive integer 64 bit variant :
C    it is assumed that IW array size can exceed 32-bit integer
C
C     Input not modified
      INTEGER, INTENT(IN)     :: TOTEL, N
      INTEGER(8), INTENT(IN)  :: IWLEN
      LOGICAL, INTENT(IN)     :: COMPUTE_PERM
C     Input undefined on output 
      INTEGER, INTENT(INOUT)  :: LEN(N), IW(IWLEN)
C 
C     Output only 
      INTEGER, INTENT(OUT)   :: NCMPA
      INTEGER, INTENT(OUT)   :: ELEN(N), LAST(N), PARENT(N)
C 
C     Input/output
      INTEGER(8), INTENT(INOUT) :: PFREE
      INTEGER(8), INTENT(INOUT) :: PE(N)
C     NV also meaningful as input to encode compressed graphs
      INTEGER, INTENT(INOUT) :: NV(N)
C 
C     Internal Workspace only
      INTEGER         :: NEXT(N), DEGREE(N), HEAD(TOTEL), W(N)
C ---------------------
C Interface Description
C ---------------------
C INPUT ARGUMENTS (unaltered):
C-----------------------------
C n     : The matrix order.
C         number of supervariables if compress/blocked format
C         Restriction:  n .ge. 1
C totel : Number of variables to eliminate
C         In case of blocked format:
C         each variable i is a supervariable of size nv(i)
C         totel is computed as the sum(nv(i)) for i \in [1:n]
C         the algorithm stops when totel variables are
C         eliminated.
C compute_perm : indicates if permutations should be computed 
C         on output in last/elen 
C iwlen:        The length of iw (1..iwlen).  On input, the matrix is
C       stored in iw (1..pfree-1).  However, iw (1..iwlen) should be
C       slightly larger than what is required to hold the matrix, at
C       least iwlen .ge. pfree + n is recommended.  Otherwise,
C       excessive compressions will take place.
C       *** We do not recommend running this algorithm with ***
C       ***      iwlen .lt. pfree + n.                      ***
C       *** Better performance will be obtained if          ***
C       ***      iwlen .ge. pfree + n                       ***
C       *** or better yet                                   ***
C       ***      iwlen .gt. 1.2 * pfree                     ***
C       *** (where pfree is its value on input).            ***
C       The algorithm will not run at all if iwlen .lt. pfree-1.
C
C       Restriction: iwlen .ge. pfree-1
C-----------------------------------------------------------------------
C INPUT/OUPUT ARGUMENTS:
C-----------------------------------------------------------------------
C pe:   On input, pe (i) is the index in iw of the start of row i, or
C       zero if row i has no off-diagonal non-zeros.
C
C       During execution, it is used for both supervariables and
C       elements:
C
C       * Principal supervariable i:  index into iw of the
C               description of supervariable i.  A supervariable
C               represents one or more rows of the matrix
C               with identical nonzero pattern.
C       * Non-principal supervariable i:  if i has been absorbed
C               into another supervariable j, then pe (i) = -j.
C               That is, j has the same pattern as i.
C               Note that j might later be absorbed into another
C               supervariable j2, in which case pe (i) is still -j,
C               and pe (j) = -j2.
C       * Unabsorbed element e:  the index into iw of the description
C               of element e, if e has not yet been absorbed by a
C               subsequent element.  Element e is created when
C               the supervariable of the same name is selected as
C               the pivot.
C       * Absorbed element e:  if element e is absorbed into element
C               e2, then pe (e) = -e2.  This occurs when the pattern of
C               e (that is, Le) is found to be a subset of the pattern
C               of e2 (that is, Le2).  If element e is "null" (it has
C               no nonzeros outside its pivot block), then pe (e) = 0.
C
C       On output, pe holds the assembly tree/forest, which implicitly
C       represents a pivot order with identical fill-in as the actual
C       order (via a depth-first search of the tree).
C
C       On output:  (PE is copied on output into PARENT array)
C       If nv (i) .gt. 0, then i represents a node in the assembly tree,
C       and the parent of i is -pe (i), or zero if i is a root.
C       If nv (i) = 0, then (i,-pe (i)) represents an edge in a
C       subtree, the root of which is a node in the assembly tree.
C pfree:On input, the matrix is stored in iw (1..pfree-1) and
C       the rest of the array iw is free.
C       During execution, additional data is placed in iw, and pfree
C       is modified so that components  of iw from pfree are free.
C       On output, pfree is set equal to the size of iw that
C       would have been needed for no compressions to occur.  If
C       ncmpa is zero, then pfree (on output) is less than or equal to
C       iwlen, and the space iw (pfree+1 ... iwlen) was not used.
C       Otherwise, pfree (on output) is greater than iwlen, and all the
C       memory in iw was used.
C nv:   On input, encoding of compressed graph:
C       if nv(1) = -1 then graph is not compressed otherwise
C       nv(I) holds the weight of node I. 
C       During execution, abs (nv (i)) is equal to the number of rows
C       that are represented by the principal supervariable i.  If i is
C       a nonprincipal variable, then nv (i) = 0.  
C       nv (i) .lt. 0 signifies that i is a
C       principal variable in the pattern Lme of the current pivot
C       element me.  
C       On output, nv (e) holds the true degree of element
C       e at the time it was created (including the diagonal part).
C-----------------------------------------------------------------------
C INPUT/MODIFIED (undefined on output):
C-----------------------------------------------------------------------
C len:  On input, len (i) holds the number of entries in row i of the
C       matrix, excluding the diagonal.  The contents of len (1..n)
C       are undefined on output.
C iw:   On input, iw (1..pfree-1) holds the description of each row i
C       in the matrix.  The matrix must be symmetric, and both upper
C       and lower triangular parts must be present.  The diagonal must
C       not be present.  
C       Row i is held as follows:
C               len (i):  the length of the row i data structure
C               iw (pe (i) ... pe (i) + len (i) - 1):
C                       the list of column indices for nonzeros
C                       in row i (simple supervariables), excluding
C                       the diagonal.  All supervariables start with
C                       one row/column each (supervariable i is just
C                       row i).
C               if len (i) is zero on input, then pe (i) is ignored
C               on input.
C
C               Note that the rows need not be in any particular order,
C               and there may be empty space between the rows.
C
C       During execution, the supervariable i experiences fill-in.
C       This is represented by placing in i a list of the elements
C       that cause fill-in in supervariable i:
C
C               len (i):  the length of supervariable i
C               iw (pe (i) ... pe (i) + elen (i) - 1):
C                       the list of elements that contain i.  This list
C                       is kept short by removing absorbed elements.
C               iw (pe (i) + elen (i) ... pe (i) + len (i) - 1):
C                       the list of supervariables in i.  This list
C                       is kept short by removing nonprincipal
C                       variables, and any entry j that is also
C                       contained in at least one of the elements
C                       (j in Le) in the list for i (e in row i).
C
C       When supervariable i is selected as pivot, we create an
C       element e of the same name (e=i):
C
C               len (e):  the length of element e
C               iw (pe (e) ... pe (e) + len (e) - 1):
C                       the list of supervariables in element e.
C
C       An element represents the fill-in that occurs when supervariable
C       i is selected as pivot (which represents the selection of row i
C       and all non-principal variables whose principal variable is i).
C       We use the term Le to denote the set of all supervariables
C       in element e.  Absorbed supervariables and elements are pruned
C       from these lists when computationally convenient.
C
C       CAUTION:  THE INPUT MATRIX IS OVERWRITTEN DURING COMPUTATION.
C       The contents of iw are undefined on output.
C-----------------------------------------------------------------------
C OUTPUT (need not be set on input):
C-----------------------------------------------------------------------
C elen: 
C       See the description of iw above.  At the start of execution,
C       elen (i) is set to zero.  During execution, elen (i) is the
C       number of elements in the list for supervariable i.  When e
C       becomes an element, elen (e) = -nel is set, where nel is the
C       current step of factorization.  elen (i) = 0 is done when i
C       becomes nonprincipal.
C
C       For variables, elen (i) .ge. 0 holds until just before the
C       permutation vectors are computed.  For elements,
C       elen (e) .lt. 0 holds.
C
C       On output elen (1..n) holds the inverse permutation (the same
C       as the 'INVP' argument in Sparspak).  That is, if k = elen (i),
C       then row i is the kth pivot row.  Row i of A appears as the
C       (elen(i))-th row in the permuted matrix, PAP^T.
C last: 
C       In a degree list, last (i) is the supervariable preceding i,
C       or zero if i is the head of the list.  In a hash bucket,
C       last (i) is the hash key for i.  last (head (hash)) is also
C       used as the head of a hash bucket if head (hash) contains a
C       degree list (see head, below).
C
C       On output, last (1..n) holds the permutation (the same as the
C       'PERM' argument in Sparspak).  That is, if i = last (k), then
C       row i is the kth pivot row.  Row last (k) of A is the k-th row
C       in the permuted matrix, PAP^T.
C ncmpa:        The number of times iw was compressed.  If this is
C       excessive, then the execution took longer than what could have
C       been.  To reduce ncmpa, try increasing iwlen to be 10% or 20%
C       larger than the value of pfree on input (or at least
C       iwlen .ge. pfree + n).  The fastest performance will be
C       obtained when ncmpa is returned as zero.  If iwlen is set to
C       the value returned by pfree on *output*, then no compressions
C       will occur.
C-----------------------------------------------------------------------
C LOCAL (not input or output - used only during execution):
C-----------------------------------------------------------------------
C degree:       If i is a supervariable, then degree (i) holds the
C       current approximation of the external degree of row i (an upper
C       bound).  The external degree is the number of nonzeros in row i,
C       minus abs (nv (i)) (the diagonal part).  The bound is equal to
C       the external degree if elen (i) is less than or equal to two.
C
C       We also use the term "external degree" for elements e to refer
C       to |Le \ Lme|.  If e is an element, then degree (e) holds |Le|,
C       which is the degree of the off-diagonal part of the element e
C       (not including the diagonal part).
C head: head is used for degree lists.  head (deg) is the first
C       supervariable in a degree list (all supervariables i in a
C       degree list deg have the same approximate degree, namely,
C       deg = degree (i)).  If the list deg is empty then
C       head (deg) = 0.
C
C       During supervariable detection head (hash) also serves as a
C       pointer to a hash bucket.
C       If head (hash) .gt. 0, there is a degree list of degree hash.
C               The hash bucket head pointer is last (head (hash)).
C       If head (hash) = 0, then the degree list and hash bucket are
C               both empty.
C       If head (hash) .lt. 0, then the degree list is empty, and
C               -head (hash) is the head of the hash bucket.
C       After supervariable detection is complete, all hash buckets
C       are empty, and the (last (head (hash)) = 0) condition is
C       restored for the non-empty degree lists.
C next: next (i) is the supervariable following i in a link list, or
C       zero if i is the last in the list.  Used for two kinds of
C       lists:  degree lists and hash buckets (a supervariable can be
C       in only one kind of list at a time).
C w:    The flag array w determines the status of elements and
C       variables, and the external degree of elements.
C
C       for elements:
C          if w (e) = 0, then the element e is absorbed
C          if w (e) .ge. wflg, then w (e) - wflg is the size of
C               the set |Le \ Lme|, in terms of nonzeros (the
C               sum of abs (nv (i)) for each principal variable i that
C               is both in the pattern of element e and NOT in the
C               pattern of the current pivot element, me).
C          if wflg .gt. w (e) .gt. 0, then e is not absorbed and has
C               not yet been seen in the scan of the element lists in
C               the computation of |Le\Lme| in loop 150 below.
C
C       for variables:
C          during supervariable detection, if w (j) .ne. wflg then j is
C          not in the pattern of variable i
C
C       The w array is initialized by setting w (i) = 1 for all i,
C       and by setting wflg = 2.  It is reinitialized if wflg becomes
C       too large (to ensure that wflg+n does not cause integer
C       overflow).
C-----------------------------------------------------------------------
C LOCAL INTEGERS:
C-----------------------------------------------------------------------
      INTEGER   :: DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     &        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     &        LENJ, LN,  ME,  MINDEG, NEL, 
     &        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X
      INTEGER KNT1_UPDATED, KNT2_UPDATED
      INTEGER(8) :: MAXMEM, MEM, NEWMEM
      INTEGER    :: MAXINT_N
      INTEGER(8) :: HASH, HMOD
C deg:        the degree of a variable or element
C degme:      size, |Lme|, of the current element, me (= degree (me))
C dext:       external degree, |Le \ Lme|, of some element e
C dmax:       largest |Le| seen so far
C e:          an element
C elenme:     the length, elen (me), of element list of pivotal var.
C eln:        the length, elen (...), of an element list
C hash:       the computed value of the hash function
C hmod:       the hash function is computed modulo hmod = max (1,n-1)
C i:          a supervariable
C ilast:      the entry in a link list preceding i
C inext:      the entry in a link list following i
C j:          a supervariable
C jlast:      the entry in a link list preceding j
C jnext:      the entry in a link list, or path, following j
C k:          the pivot order of an element or variable
C knt1:       loop counter used during element construction
C knt2:       loop counter used during element construction
C knt3:       loop counter used during compression
C lenj:       len (j)
C ln:         length of a supervariable list
C maxint_n    large integer to test risk of overflow on wflg
C maxmem:     amount of memory needed for no compressions
C me:         current supervariable being eliminated, and the
C                     current element created by eliminating that
C                     supervariable
C mem:        memory in use assuming no compressions have occurred
C mindeg:     current minimum degree
C nel:        number of pivots selected so far
C newmem:     amount of new memory needed for current pivot element
C nleft:      n - nel, the number of nonpivotal rows/columns remaining
C nvi:        the number of variables in a supervariable i (= nv (i))
C nvj:        the number of variables in a supervariable j (= nv (j))
C nvpiv:      number of pivots in current element
C slenme:     number of variables in variable list of pivotal variable
C we:         w (e)
C wflg:       used for flagging the w array.  See description of iw.
C wnvi:       wflg - nv (i)
C x:          either a supervariable or an element
C-----------------------------------------------------------------------
C LOCAL POINTERS:
C-----------------------------------------------------------------------
      INTEGER(8) P, P1, P2, P3, PDST, PEND, PJ, PME, 
     &           PME1, PME2, PN, PSRC
C             Any parameter (pe (...) or pfree) or local variable
C             starting with "p" (for Pointer) is an index into iw,
C             and all indices into iw use variables starting with
C             "p."  The only exception to this rule is the iwlen
C             input argument.
C p:          pointer into lots of things
C p1:         pe (i) for some variable i (start of element list)
C p2:         pe (i) + elen (i) -  1 for some var. i (end of el. list)
C p3:         index of first supervariable in clean list
C pdst:       destination pointer, for compression
C pend:       end of memory to compress
C pj:         pointer into an element or variable
C pme:        pointer into the current element (pme1...pme2)
C pme1:       the current element, me, is stored in iw (pme1...pme2)
C pme2:       the end of the current element
C pn:         pointer into a "clean" variable, also used to compress
C psrc:       source pointer, for compression
      LOGICAL COMPRESS
C-----------------------------------------------------------------------
C  FUNCTIONS CALLED:
C-----------------------------------------------------------------------
      INTRINSIC max, min, mod
C=======================================================================
C  INITIALIZATIONS
C=======================================================================
      WFLG = 2
      MAXINT_N=huge(WFLG)-N
      MINDEG = 1
      NCMPA = 0
      NEL = 0
      HMOD = int(max (1, N-1),kind=8)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      DO I = 1, N
        LAST (I) = 0
        HEAD (I) = 0
        W (I) = 1
        ELEN (I) = 0
      ENDDO
      DO I = 1, TOTEL
        HEAD(I) = 0
      ENDDO
      IF(NV(1) .LT. 0) THEN
         COMPRESS = .FALSE.
      ELSE
         COMPRESS = .TRUE.
      ENDIF
      IF (COMPRESS) THEN
         DO I=1,N
            DEGREE(I) = 0
            DO P= PE(I) , PE(I)+int(LEN(I)-1,8)
               DEGREE(I) = DEGREE(I) + NV(IW(P))
            ENDDO
         ENDDO
      ELSE
         DO I=1,N
            NV(I) = 1
            DEGREE (I) = LEN (I)
         ENDDO
      ENDIF
C
C     ----------------------------------------------------------------
C     initialize degree lists and eliminate rows with no off-diag. nz.
C     ----------------------------------------------------------------
      DO 20 I = 1, N
        DEG = DEGREE (I)
        IF (DEG .GT. 0) THEN
C         ----------------------------------------------------------
C         place i in the degree list corresponding to its degree
C         ----------------------------------------------------------
          INEXT = HEAD (DEG)
          IF (INEXT .NE. 0) LAST (INEXT) = I
          NEXT (I) = INEXT
          HEAD (DEG) = I
        ELSE
C         ----------------------------------------------------------
C         we have a variable that can be eliminated at once because
C         there is no off-diagonal non-zero in its row.
C         ----------------------------------------------------------
            NEL = NEL + NV(I)
          ELEN (I) = -NEL
          PE (I) = 0
          W (I) = 0
        ENDIF
   20 CONTINUE
C =====================================================================
C  WHILE (selecting pivots) DO
C =====================================================================
   30 IF (NEL .LT. TOTEL) THEN
C =====================================================================
C  GET PIVOT OF MINIMUM DEGREE
C ======================================================================
C       -------------------------------------------------------------
C       find next supervariable for elimination
C       -------------------------------------------------------------
        DO 40 DEG = MINDEG, TOTEL
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
C       -------------------------------------------------------------
C       remove chosen variable from link list
C       -------------------------------------------------------------
        INEXT = NEXT (ME)
        IF (INEXT .NE. 0) LAST (INEXT) = 0
        HEAD (DEG) = INEXT
C       -------------------------------------------------------------
C       me represents the elimination of pivots nel+1 to nel+nv(me).
C       place me itself as the first in this set.  It will be moved
C       to the nel+nv(me) position when the permutation vectors are
C       computed.
C       -------------------------------------------------------------
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
C=======================================================================
C  CONSTRUCT NEW ELEMENT
C=======================================================================
C       -------------------------------------------------------------
C       At this point, me is the pivotal supervariable.  It will be
C       converted into the current element.  Scan list of the
C       pivotal supervariable, me, setting tree pointers and
C       constructing new list of supervariables for the new element,
C       me.  p is a pointer to the current position in the old list.
C       -------------------------------------------------------------
C       flag the variable "me" as being in Lme by negating nv (me)
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
C         ----------------------------------------------------------
C         construct the new element in place
C         ----------------------------------------------------------
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
C             ----------------------------------------------------
C             i is a principal variable not yet placed in Lme.
C             store i in new list
C             ----------------------------------------------------
              DEGME = DEGME + NVI
C             flag i as being in Lme by negating nv (i)
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
C             ----------------------------------------------------
C             remove variable i from degree list.
C             ----------------------------------------------------
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
C               i is at the head of the degree list
                HEAD (DEGREE (I)) = INEXT
              ENDIF
            ENDIF
   60     CONTINUE
C         this element takes no new memory in iw:
          NEWMEM = 0
        ELSE
C         ----------------------------------------------------------
C         construct the new element in empty space, iw (pfree ...)
C         ----------------------------------------------------------
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          KNT1_UPDATED = 0
          DO 120 KNT1 = 1, ELENME + 1
            KNT1_UPDATED = KNT1_UPDATED +1
            IF (KNT1 .GT. ELENME) THEN
C             search the supervariables in me.
              E = ME
              PJ = P
              LN = SLENME
            ELSE
C             search the elements in me.
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
C           -------------------------------------------------------
C           search for different supervariables and add them to the
C           new list, compressing when necessary. this loop is
C           executed once for each element in the list and once for
C           all the supervariables in the list.
C           -------------------------------------------------------
            KNT2_UPDATED = 0
            DO 110 KNT2 = 1, LN
              KNT2_UPDATED = KNT2_UPDATED+1
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
C               -------------------------------------------------
C               compress iw, if necessary
C               -------------------------------------------------
                IF (PFREE .GT. IWLEN) THEN
C                 prepare for compressing iw by adjusting
C                 pointers and lengths so that the lists being
C                 searched in the inner and outer loops contain
C                 only the remaining entries.
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1_UPDATED
C                 Reset KNT1_UPDATED in case of recompress 
C                 at same iteration of the loop 120
                  KNT1_UPDATED = 0
C                 Check if anything left in supervariable ME
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2_UPDATED
C                 Reset KNT2_UPDATED in case of recompress 
C                 at same iteration of the loop 110
                  KNT2_UPDATED = 0
C                 Check if anything left in element E
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
C                 store first item in pe
C                 set first entry to -item
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = int(IW (PN), 8)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
C                 psrc/pdst point to source/destination
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
C                 while loop:
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
C                   search for next negative entry
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = int(PE (J))
                      PE (J) = PDST
                      PDST = PDST + 1
C                     copy from source to destination
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
C                 move the new partially-constructed element
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
C               -------------------------------------------------
C               i is a principal variable not yet placed in Lme
C               store i in new list
C               -------------------------------------------------
                DEGME = DEGME + NVI
C               flag i as being in Lme by negating nv (i)
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
C               -------------------------------------------------
C               remove variable i from degree link list
C               -------------------------------------------------
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
C                 i is at the head of the degree list
                  HEAD (DEGREE (I)) = INEXT
                ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
C             set tree pointer and flag to indicate element e is
C             absorbed into new element me (the parent of e is me)
              PE (E) = int(-ME,8)
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
C         this element takes newmem new memory in iw (possibly zero)
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = max (MAXMEM, MEM)
        ENDIF
C       -------------------------------------------------------------
C       me has now been converted into an element in iw (pme1..pme2)
C       -------------------------------------------------------------
C       degme holds the external degree of new element
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = int(PME2 - PME1 + 1)
C       -------------------------------------------------------------
C       make sure that wflg is not too large.  With the current
C       value of wflg, wflg+n must not cause integer overflow
C       -------------------------------------------------------------
        IF (WFLG .GT. MAXINT_N) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
C=======================================================================
C  COMPUTE (w (e) - wflg) = |Le\Lme| FOR ALL ELEMENTS
C=======================================================================
C       -------------------------------------------------------------
C       Scan 1:  compute the external degrees of previous elements
C       with respect to the current element.  That is:
C            (w (e) - wflg) = |Le \ Lme|
C       for each element e that appears in any supervariable in Lme.
C       The notation Le refers to the pattern (list of
C       supervariables) of a previous element e, where e is not yet
C       absorbed, stored in iw (pe (e) + 1 ... pe (e) + iw (pe (e))).
C       The notation Lme refers to the pattern of the current element
C       (stored in iw (pme1..pme2)).   If (w (e) - wflg) becomes
C       zero, then the element e will be absorbed in scan 2.
C       -------------------------------------------------------------
        DO 150 PME = PME1, PME2
          I = IW (PME)
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
C           note that nv (i) has been negated to denote i in Lme:
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + ELN - 1
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
C               unabsorbed element e has been seen in this loop
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
C               e is an unabsorbed element
C               this is the first we have seen e in all of Scan 1
                WE = DEGREE (E) + WNVI
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
C=======================================================================
C  DEGREE UPDATE AND ELEMENT ABSORPTION
C=======================================================================
C       -------------------------------------------------------------
C       Scan 2:  for each i in Lme, sum up the degree of Lme 
C       (which is degme), 
C       plus the sum of the external degrees of each Le
C       for the elements e appearing within i, plus the
C       supervariables in i.  Place i in hash list.
C       -------------------------------------------------------------
        DO 180 PME = PME1, PME2
          I = IW (PME)
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0_8
          DEG = 0
C         ----------------------------------------------------------
C         scan the element list associated with supervariable i
C         ----------------------------------------------------------
          DO 160 P = P1, P2
            E = IW (P)
C           dext = | Le \ Lme |
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
            ELSE IF (DEXT .EQ. 0) THEN
#if  defined (NOAGG1)
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
#else
C             aggressive absorption: e is not adjacent to me, but
C             the |Le \ Lme| is 0, so absorb it into me
              PE (E) = int(-ME,8)
              W (E) = 0
#endif
            ENDIF
  160     CONTINUE
C         count the number of elements in i (including me):
          ELEN (I) = int(PN - P1 + 1)
C         ----------------------------------------------------------
C         scan the supervariables in the list associated with i
C         ----------------------------------------------------------
          P3 = PN
          DO 170 P = P2 + 1, P1 + int(LEN (I) - 1,8)
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
C             j is unabsorbed, and not in Lme.
C             add to degree and add to new list
              DEG = DEG + NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + int(J,kind=8)
            ENDIF
  170     CONTINUE
C         ----------------------------------------------------------
C         update the degree and check for mass elimination
C         ----------------------------------------------------------
#if  defined (NOAGG1)
          IF (DEG.EQ.0.AND.(ELEN(I).GT.1)) THEN
C         When DEG is zero we need to 
C         absorb in ME all elements adjacent to I 
                   P1 = PE (I)
C                  exclude ME --> -2 
                   P2 = P1 + int(ELEN (I),8) - 2_8
                   DO P =P1,P2  
                     E      = IW(P)
                     PE (E) = int(-ME,8)
                     W (E)  = 0
                   ENDDO
          ENDIF
C              .... Ready for mass elimination
#endif
          IF (DEG .EQ. 0) THEN
C           -------------------------------------------------------
C           mass elimination
C           -------------------------------------------------------
C           There is nothing left of this node except for an
C           edge to the current pivot element.  elen (i) is 1,
C           and there are no variables adjacent to node i.
C           Absorb i into the current pivot element, me.
            PE (I) = int(-ME,8)
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
C           -------------------------------------------------------
C           update the upper-bound degree of i
C           -------------------------------------------------------
C           the following degree does not yet include the size
C           of the current element, which is added later:
            DEGREE (I) = min (DEGREE (I), DEG)
C           -------------------------------------------------------
C           add me to the list for i
C           -------------------------------------------------------
C           move first supervariable to end of list
            IW (PN) = IW (P3)
C           move first element to end of element part of list
            IW (P3) = IW (P1)
C           add new element to front of list.
            IW (P1) = ME
C           store the new length of the list in len (i)
            LEN (I) = int(PN - P1 + 1)
C           -------------------------------------------------------
C           place in hash bucket.  Save hash key of i in last (i).
C           -------------------------------------------------------
            HASH = mod (HASH, HMOD) + 1_8
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
C             the degree list is empty, hash head is -j
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
C             degree list is not empty
C             use last (head (hash)) as hash head
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = int(HASH,kind=kind(LAST))
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
C       -------------------------------------------------------------
C       Clear the counter array, w (...), by incrementing wflg.
C       -------------------------------------------------------------
        DMAX = max (DMAX, DEGME)
        WFLG = WFLG + DMAX
C       make sure that wflg+n does not cause integer overflow
        IF (WFLG .GT. MAXINT_N) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
C       at this point, w (1..n) .lt. wflg holds
C=======================================================================
C  SUPERVARIABLE DETECTION
C=======================================================================
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF (NV (I) .LT. 0) THEN
C           i is a principal variable in Lme
C           -------------------------------------------------------
C           examine all hash buckets with 2 or more variables.  We
C           do this by examing all unique hash keys for super-
C           variables in the pattern Lme of the current element, me
C           -------------------------------------------------------
            HASH = int(LAST (I),kind=8)
C           let i = head of hash bucket, and empty the hash bucket
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
C             degree list is empty
              I = -J
              HEAD (HASH) = 0
            ELSE
C             degree list is not empty, restore last () of head
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
C           while loop:
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
C             ----------------------------------------------------
C             this bucket has one or more variables following i.
C             scan all of them to see if i can absorb any entries
C             that follow i in hash bucket.  Scatter i into w.
C             ----------------------------------------------------
              LN = LEN (I)
              ELN = ELEN (I)
C             do not flag the first element in the list (me)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
C             ----------------------------------------------------
C             scan every other entry j following i in bucket
C             ----------------------------------------------------
              JLAST = I
              J = NEXT (I)
C             while loop:
  220         CONTINUE
              IF (J .NE. 0) THEN
C               -------------------------------------------------
C               check if j and i have identical nonzero pattern
C               -------------------------------------------------
C               jump if i and j do not have same size data structure
                IF (LEN (J) .NE. LN) GO TO 240
C               jump if i and j do not have same number adj elts
                IF (ELEN (J) .NE. ELN) GO TO 240
C               do not flag the first element in the list (me)
                DO 230 P = PE (J) + 1, PE (J) + LN - 1
C                 jump if an entry (iw(p)) is in j but not in i
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
C               -------------------------------------------------
C               found it!  j can be absorbed into i
C               -------------------------------------------------
                PE (J) = int(-I,8)
C               both nv (i) and nv (j) are negated since they
C               are in Lme, and the absolute values of each
C               are the number of variables in i and j:
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
C               delete j from hash bucket
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
C               -------------------------------------------------
  240           CONTINUE
C               j cannot be absorbed into i
C               -------------------------------------------------
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
C             ----------------------------------------------------
C             no more variables can be absorbed into i
C             go to next i in bucket and clear flag array
C             ----------------------------------------------------
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
C=======================================================================
C  RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVAR. FROM ELEMENT
C=======================================================================
        P = PME1
        NLEFT = TOTEL - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
C           i is a principal variable in Lme
C           restore nv (i) to signify that i is principal
            NV (I) = NVI
C           -------------------------------------------------------
C           compute the external degree (add size of current elem)
C           -------------------------------------------------------
            DEG = min (DEGREE (I) + DEGME - NVI, NLEFT - NVI)
C           -------------------------------------------------------
C           place the supervariable at the head of the degree list
C           -------------------------------------------------------
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (DEG) = I
C           -------------------------------------------------------
C           save the new degree, and find the minimum degree
C           -------------------------------------------------------
            MINDEG = min (MINDEG, DEG)
            DEGREE (I) = DEG
C           -------------------------------------------------------
C           place the supervariable in the element pattern
C           -------------------------------------------------------
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
C=======================================================================
C  FINALIZE THE NEW ELEMENT
C=======================================================================
        NV (ME) = NVPIV + DEGME
C       nv (me) is now the degree of pivot (including diagonal part)
C       save the length of the list for the new element me
        LEN (ME) = int(P - PME1)
        IF (LEN (ME) .EQ. 0) THEN
C         there is nothing left of the current pivot element
          PE (ME) = 0_8
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
C         element was not constructed in place: deallocate part
C         of it (final size is less than or equal to newmem,
C         since newly nonprincipal variables have been removed).
          PFREE = P
          MEM = MEM - NEWMEM + LEN (ME)
        ENDIF
C=======================================================================
C       END WHILE (selecting pivots)
      GO TO 30
      ENDIF
C=======================================================================
C=======================================================================
C     COMPUTE THE PERMUTATION VECTORS and update TREE
C=======================================================================
C     ----------------------------------------------------------------
C     The time taken by the following code is O(n).  At this
C     point, elen (e) = -k has been done for all elements e,
C     and elen (i) = 0 has been done for all nonprincipal
C     variables i.  At this point, there are no principal
C     supervariables left, and all elements are absorbed.
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
C     compute the ordering of unordered nonprincipal variables
C     ----------------------------------------------------------------
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
C         ----------------------------------------------------------
C         i is an un-ordered row.  Traverse the tree from i until
C         reaching an element, e.  The element, e, was the
C         principal supervariable of i and all nodes in the path
C         from i to when e was selected as pivot.
C         ----------------------------------------------------------
          J = int(-PE (I))
C         while (j is a variable) do:
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = int(-PE (J))
              GO TO 270
            ENDIF
            E = J
C           ----------------------------------------------------------
C           get the current pivot ordering of e
C           ----------------------------------------------------------
            K = -ELEN (E)
C           ----------------------------------------------------------
C           traverse the path again from i to e, and compress the
C           path (all nodes point to e).  Path compression allows
C           this code to compute in O(n) time.  Order the unordered
C           nodes in the path, and place the element e at the end.
C           ----------------------------------------------------------
            J = I
C           while (j is a variable) do:
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = int(-PE (J))
              PE (J) = int(-E,8)
              IF (ELEN (J) .EQ. 0) THEN
C               j is an unordered row
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
C         leave elen (e) negative, so we know it is an element
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
C
      IF (COMPUTE_PERM) THEN
C     ----------------------------------------------------------------
C     reset the inverse permutation (elen (1..n)) to be positive,
C     and compute the permutation (last (1..n)).
C     ----------------------------------------------------------------
      IF(COMPRESS) THEN
        LAST(1:N) = 0
        HEAD(1:TOTEL-N)=0  
        DO I = 1, N
          K = abs (ELEN (I))
          IF ( K <= N ) THEN
            LAST (K) = I
          ELSE
            HEAD(K-N)=I
          ENDIF
        ENDDO
        I = 1
        DO K = 1, N
          IF(LAST (K) .NE. 0) THEN
            LAST(I) = LAST(K)
            ELEN(LAST(K)) = I
            I = I + 1
          ENDIF
        ENDDO
        DO K = N+1, TOTEL
          IF (HEAD(K-N) .NE. 0) THEN
            LAST(I)=HEAD(K-N)
            ELEN(HEAD(K-N)) = I
            I = I + 1
          ENDIF
        END DO
      ELSE
        DO 300 I = 1, N
          K = abs (ELEN (I))
          LAST (K) = I
          ELEN (I) = K
  300   CONTINUE
       ENDIF
C=======================================================================
C      END OF COMPUTING PERMUTATIONS
C=======================================================================
       ENDIF
C=======================================================================
C  RETURN THE MEMORY USAGE IN IW
C=======================================================================
C     If maxmem is less than or equal to iwlen, then no compressions
C     occurred, and iw (maxmem+1 ... iwlen) was unused.  Otherwise
C     compressions did occur, and iwlen would have had to have been
C     greater than or equal to maxmem for no compressions to occur.
C     Return the value of maxmem in the pfree argument.
      PFREE = MAXMEM
C===============================
C     Save IPE in PARENT array
      DO I=1,N
       PARENT(I) = int(PE(I))
      ENDDO
C===============================
      RETURN
      END SUBROUTINE MUMPS_ANA_H
C-----------------------------------------------------------------------
C MUMPS_AMD_ELT: modified version of reference AMD routine MUMPS_ANA_H
C capable of processing already amalgamated or compressed graph. 
C Used within MUMPS process for the elemental input format of matrices
C Input data is in this context modified to be a graph of supervariables.
C   
C Modifications of the interface : 
C ------------------------------
C  INPUT:
C  -----
C     1/ LEN(I) < 0   <=> i is a secondary variable whose principal
C                         variable is -LEN(I)
C     2/ For all secondary variables the adj list MUST not be provided. 
C        THAT is:
C        -------
C           if pe(isecondary) = 0 then 
C                 adjacency list of isecondary is not provided
C           else
C             pe(isecondary) >0 
C             len(isecondary) must be equal to len(iprincipal_associated)
C             then the corresponding space wil not be used and 
C             will be freed by amd if necessary.
C           endif
C REMARK:
C ------
C 1/ N must be still set to the order of the matrix 
C    (not of the amalgamated gragh)
C 2/ For each supervariable S only supervariables adjacent to S are provided
C    len(S) is then the number of such supervariables
C    NV(S) is however updated during the initialisation phase to represent 
C    the size of the supervariable 
C    ( increment nv(s) for each i / len(i) =-s )
C 3/ If (len(i) >=0  for all i ) then we get the classical AMD code
C ------------------
      SUBROUTINE MUMPS_AMD_ELT(N,IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
     &                   LAST, NCMPA, DEGREE, HEAD, NEXT, W, PARENT)
C 
C    Restrictive integer 64 bit variant :
C    it is assumed that IW array size can exceed 32-bit integer
C 
C     Input not modified
      INTEGER, INTENT(IN)    :: N
      INTEGER(8), INTENT(IN) :: IWLEN
C     Input undefined on output 
      INTEGER, INTENT(INOUT)  :: LEN(N), IW(IWLEN)
C 
C     Output only 
      INTEGER, INTENT(OUT)   :: NCMPA
      INTEGER, INTENT(OUT)   :: NV(N), ELEN(N), LAST(N), PARENT(N)
C 
C     Input/output
      INTEGER(8), INTENT(INOUT) :: PFREE
      INTEGER(8), INTENT(INOUT) :: PE(N)
C
C     Internal Workspace only
      INTEGER NEXT(N), DEGREE(N), HEAD(N), W(N)
C
C Description:
C   Given a representation of the nonzero pattern of a symmetric matrix,
C   A, (excluding the diagonal) perform an approximate minimum
C   degree ordering to compute a pivot order
C   such that fill-in in the Cholesky factors A = LL^T is kept low. 
C ---------------------
C Interface Description
C ---------------------
C INPUT ARGUMENTS (unaltered):
C-----------------------------
C n:    The matrix order.
C
C       Restriction:  n .ge. 1
C iwlen:        The length of iw (1..iwlen).  On input, the matrix is
C       stored in iw (1..pfree-1).  However, iw (1..iwlen) should be
C       slightly larger than what is required to hold the matrix, at
C       least iwlen .ge. pfree + n is recommended.  Otherwise,
C       excessive compressions will take place.
C       *** We do not recommend running this algorithm with ***
C       ***      iwlen .lt. pfree + n.                      ***
C       *** Better performance will be obtained if          ***
C       ***      iwlen .ge. pfree + n                       ***
C       *** or better yet                                   ***
C       ***      iwlen .gt. 1.2 * pfree                     ***
C       *** (where pfree is its value on input).            ***
C       The algorithm will not run at all if iwlen .lt. pfree-1.
C
C       Restriction: iwlen .ge. pfree-1
C-----------------------------------------------------------------------
C INPUT/OUPUT ARGUMENTS:
C-----------------------------------------------------------------------
C pe:   On input, pe (i) is the index in iw of the start of row i, or
C       zero if row i has no off-diagonal non-zeros.
C
C       During execution, it is used for both supervariables and
C       elements:
C
C       * Principal supervariable i:  index into iw of the
C               description of supervariable i.  A supervariable
C               represents one or more rows of the matrix
C               with identical nonzero pattern.
C       * Non-principal supervariable i:  if i has been absorbed
C               into another supervariable j, then pe (i) = -j.
C               That is, j has the same pattern as i.
C               Note that j might later be absorbed into another
C               supervariable j2, in which case pe (i) is still -j,
C               and pe (j) = -j2.
C       * Unabsorbed element e:  the index into iw of the description
C               of element e, if e has not yet been absorbed by a
C               subsequent element.  Element e is created when
C               the supervariable of the same name is selected as
C               the pivot.
C       * Absorbed element e:  if element e is absorbed into element
C               e2, then pe (e) = -e2.  This occurs when the pattern of
C               e (that is, Le) is found to be a subset of the pattern
C               of e2 (that is, Le2).  If element e is "null" (it has
C               no nonzeros outside its pivot block), then pe (e) = 0.
C
C       On output, pe holds the assembly tree/forest, which implicitly
C       represents a pivot order with identical fill-in as the actual
C       order (via a depth-first search of the tree).
C
C       On output:
C       If nv (i) .gt. 0, then i represents a node in the assembly tree,
C       and the parent of i is -pe (i), or zero if i is a root.
C       If nv (i) = 0, then (i,-pe (i)) represents an edge in a
C       subtree, the root of which is a node in the assembly tree.
C 
C       On output:  (PE is copied on output into PARENT array)
C pfree:        On input, the matrix is stored in iw (1..pfree-1) and
C       the rest of the array iw is free.
C       During execution, additional data is placed in iw, and pfree
C       is modified so that components  of iw from pfree are free.
C       On output, pfree is set equal to the size of iw that
C       would have been needed for no compressions to occur.  If
C       ncmpa is zero, then pfree (on output) is less than or equal to
C       iwlen, and the space iw (pfree+1 ... iwlen) was not used.
C       Otherwise, pfree (on output) is greater than iwlen, and all the
C       memory in iw was used.
C-----------------------------------------------------------------------
C INPUT/MODIFIED (undefined on output):
C-----------------------------------------------------------------------
C len:  On input, len (i) holds the number of entries in row i of the
C       matrix, excluding the diagonal.  The contents of len (1..n)
C       are undefined on output.
C iw:   On input, iw (1..pfree-1) holds the description of each row i
C       in the matrix.  The matrix must be symmetric, and both upper
C       and lower triangular parts must be present.  The diagonal must
C       not be present.  Row i is held as follows:
C
C               len (i):  the length of the row i data structure
C               iw (pe (i) ... pe (i) + len (i) - 1):
C                       the list of column indices for nonzeros
C                       in row i (simple supervariables), excluding
C                       the diagonal.  All supervariables start with
C                       one row/column each (supervariable i is just
C                       row i).
C               if len (i) is zero on input, then pe (i) is ignored
C               on input.
C
C               Note that the rows need not be in any particular order,
C               and there may be empty space between the rows.
C
C       During execution, the supervariable i experiences fill-in.
C       This is represented by placing in i a list of the elements
C       that cause fill-in in supervariable i:
C
C               len (i):  the length of supervariable i
C               iw (pe (i) ... pe (i) + elen (i) - 1):
C                       the list of elements that contain i.  This list
C                       is kept short by removing absorbed elements.
C               iw (pe (i) + elen (i) ... pe (i) + len (i) - 1):
C                       the list of supervariables in i.  This list
C                       is kept short by removing nonprincipal
C                       variables, and any entry j that is also
C                       contained in at least one of the elements
C                       (j in Le) in the list for i (e in row i).
C
C       When supervariable i is selected as pivot, we create an
C       element e of the same name (e=i):
C
C               len (e):  the length of element e
C               iw (pe (e) ... pe (e) + len (e) - 1):
C                       the list of supervariables in element e.
C
C       An element represents the fill-in that occurs when supervariable
C       i is selected as pivot (which represents the selection of row i
C       and all non-principal variables whose principal variable is i).
C       We use the term Le to denote the set of all supervariables
C       in element e.  Absorbed supervariables and elements are pruned
C       from these lists when computationally convenient.
C
C       CAUTION:  THE INPUT MATRIX IS OVERWRITTEN DURING COMPUTATION.
C       The contents of iw are undefined on output.
C-----------------------------------------------------------------------
C OUTPUT (need not be set on input):
C-----------------------------------------------------------------------
C nv:   During execution, abs (nv (i)) is equal to the number of rows
C       that are represented by the principal supervariable i.  If i is
C       a nonprincipal variable, then nv (i) = 0.  Initially,
C       nv (i) = 1 for all i.  nv (i) .lt. 0 signifies that i is a
C       principal variable in the pattern Lme of the current pivot
C       element me.  On output, nv (e) holds the true degree of element
C       e at the time it was created (including the diagonal part).
C elen: See the description of iw above.  At the start of execution,
C       elen (i) is set to zero.  During execution, elen (i) is the
C       number of elements in the list for supervariable i.  When e
C       becomes an element, elen (e) = -nel is set, where nel is the
C       current step of factorization.  elen (i) = 0 is done when i
C       becomes nonprincipal.
C
C       For variables, elen (i) .ge. 0 holds until just before the
C       permutation vectors are computed.  For elements,
C       elen (e) .lt. 0 holds.
C
C       On output elen (1..n) holds the inverse permutation (the same
C       as the 'INVP' argument in Sparspak).  That is, if k = elen (i),
C       then row i is the kth pivot row.  Row i of A appears as the
C       (elen(i))-th row in the permuted matrix, PAP^T.
C last: In a degree list, last (i) is the supervariable preceding i,
C       or zero if i is the head of the list.  In a hash bucket,
C       last (i) is the hash key for i.  last (head (hash)) is also
C       used as the head of a hash bucket if head (hash) contains a
C       degree list (see head, below).
C
C       On output, last (1..n) holds the permutation (the same as the
C       'PERM' argument in Sparspak).  That is, if i = last (k), then
C       row i is the kth pivot row.  Row last (k) of A is the k-th row
C       in the permuted matrix, PAP^T.
C ncmpa:        The number of times iw was compressed.  If this is
C       excessive, then the execution took longer than what could have
C       been.  To reduce ncmpa, try increasing iwlen to be 10% or 20%
C       larger than the value of pfree on input (or at least
C       iwlen .ge. pfree + n).  The fastest performance will be
C       obtained when ncmpa is returned as zero.  If iwlen is set to
C       the value returned by pfree on *output*, then no compressions
C       will occur.
C-----------------------------------------------------------------------
C LOCAL (not input or output - used only during execution):
C-----------------------------------------------------------------------
C degree:       If i is a supervariable, then degree (i) holds the
C       current approximation of the external degree of row i (an upper
C       bound).  The external degree is the number of nonzeros in row i,
C       minus abs (nv (i)) (the diagonal part).  The bound is equal to
C       the external degree if elen (i) is less than or equal to two.
C
C       We also use the term "external degree" for elements e to refer
C       to |Le \ Lme|.  If e is an element, then degree (e) holds |Le|,
C       which is the degree of the off-diagonal part of the element e
C       (not including the diagonal part).
C head: head is used for degree lists.  head (deg) is the first
C       supervariable in a degree list (all supervariables i in a
C       degree list deg have the same approximate degree, namely,
C       deg = degree (i)).  If the list deg is empty then
C       head (deg) = 0.
C
C       During supervariable detection head (hash) also serves as a
C       pointer to a hash bucket.
C       If head (hash) .gt. 0, there is a degree list of degree hash.
C               The hash bucket head pointer is last (head (hash)).
C       If head (hash) = 0, then the degree list and hash bucket are
C               both empty.
C       If head (hash) .lt. 0, then the degree list is empty, and
C               -head (hash) is the head of the hash bucket.
C       After supervariable detection is complete, all hash buckets
C       are empty, and the (last (head (hash)) = 0) condition is
C       restored for the non-empty degree lists.
C next: next (i) is the supervariable following i in a link list, or
C       zero if i is the last in the list.  Used for two kinds of
C       lists:  degree lists and hash buckets (a supervariable can be
C       in only one kind of list at a time).
C w:    The flag array w determines the status of elements and
C       variables, and the external degree of elements.
C
C       for elements:
C          if w (e) = 0, then the element e is absorbed
C          if w (e) .ge. wflg, then w (e) - wflg is the size of
C               the set |Le \ Lme|, in terms of nonzeros (the
C               sum of abs (nv (i)) for each principal variable i that
C               is both in the pattern of element e and NOT in the
C               pattern of the current pivot element, me).
C          if wflg .gt. w (e) .gt. 0, then e is not absorbed and has
C               not yet been seen in the scan of the element lists in
C               the computation of |Le\Lme| in loop 150 below.
C
C       for variables:
C          during supervariable detection, if w (j) .ne. wflg then j is
C          not in the pattern of variable i
C
C       The w array is initialized by setting w (i) = 1 for all i,
C       and by setting wflg = 2.  It is reinitialized if wflg becomes
C       too large (to ensure that wflg+n does not cause integer
C       overflow).
C-----------------------------------------------------------------------
C LOCAL INTEGERS:
C-----------------------------------------------------------------------
      INTEGER :: DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     &        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     &        LENJ, LN, ME, MINDEG, NEL, 
     &        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X, 
     &        NPRINC
      INTEGER KNT1_UPDATED, KNT2_UPDATED
      INTEGER(8) :: MAXMEM, MEM, NEWMEM
      INTEGER    :: MAXINT_N
      INTEGER(8) :: HASH, HMOD
C deg:        the degree of a variable or element
C degme:      size, |Lme|, of the current element, me (= degree (me))
C dext:       external degree, |Le \ Lme|, of some element e
C dmax:       largest |Le| seen so far
C e:          an element
C elenme:     the length, elen (me), of element list of pivotal var.
C eln:        the length, elen (...), of an element list
C hash:       the computed value of the hash function
C hmod:       the hash function is computed modulo hmod = max (1,n-1)
C i:          a supervariable
C ilast:      the entry in a link list preceding i
C inext:      the entry in a link list following i
C j:          a supervariable
C jlast:      the entry in a link list preceding j
C jnext:      the entry in a link list, or path, following j
C k:          the pivot order of an element or variable
C knt1:       loop counter used during element construction
C knt2:       loop counter used during element construction
C knt3:       loop counter used during compression
C lenj:       len (j)
C ln:         length of a supervariable list
C maxint_n    large integer to test risk of overflow on wflg
C maxmem:     amount of memory needed for no compressions
C me:         current supervariable being eliminated, and the
C                     current element created by eliminating that
C                     supervariable
C mem:        memory in use assuming no compressions have occurred
C mindeg:     current minimum degree
C nel:        number of pivots selected so far
C newmem:     amount of new memory needed for current pivot element
C nleft:      n - nel, the number of nonpivotal rows/columns remaining
C nvi:        the number of variables in a supervariable i (= nv (i))
C nvj:        the number of variables in a supervariable j (= nv (j))
C nvpiv:      number of pivots in current element
C slenme:     number of variables in variable list of pivotal variable
C we:         w (e)
C wflg:       used for flagging the w array.  See description of iw.
C wnvi:       wflg - nv (i)
C x:          either a supervariable or an element
C nprinc :    number of principal variables = number of varialbles
C             of the compressed graph.
C             (if the graph is not compressed then nprinc = n)
C-----------------------------------------------------------------------
C LOCAL POINTERS:
C-----------------------------------------------------------------------
      INTEGER(8) :: P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, 
     &              PN, PSRC 
C             Any parameter (pe (...) or pfree) or local variable
C             starting with "p" (for Pointer) is an index into iw,
C             and all indices into iw use variables starting with
C             "p."  The only exception to this rule is the iwlen
C             input argument.
C p:          pointer into lots of things
C p1:         pe (i) for some variable i (start of element list)
C p2:         pe (i) + elen (i) -  1 for some var. i (end of el. list)
C p3:         index of first supervariable in clean list
C pdst:       destination pointer, for compression
C pend:       end of memory to compress
C pj:         pointer into an element or variable
C pme:        pointer into the current element (pme1...pme2)
C pme1:       the current element, me, is stored in iw (pme1...pme2)
C pme2:       the end of the current element
C pn:         pointer into a "clean" variable, also used to compress
C psrc:       source pointer, for compression
C-----------------------------------------------------------------------
C  FUNCTIONS CALLED:
C-----------------------------------------------------------------------
      INTRINSIC max, min, mod
C=======================================================================
C  INITIALIZATIONS
C=======================================================================
      WFLG = 2
      MAXINT_N=huge(WFLG)-N
      MINDEG = 1
      NCMPA = 0
      NEL = 0
      HMOD = int(max (1, N-1),kind=8)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      NPRINC = 0
      DO I = 1, N
        LAST (I) = 0
        HEAD (I) = 0
        NV (I) = 1
        W (I) = 1
        ELEN (I) = 0
      ENDDO
      DO I=1, N
        IF (LEN (I).GE.0) THEN
           DEGREE (I) = LEN (I)
           NPRINC = NPRINC + 1
        ELSE
C          i is a secondary variable belonging 
C          to supervariable j=-len (i)
           J        = -LEN (I)
C          used only to skip secondary variables in loop 20
           DEGREE (I) = - 1
           IF ( PE(I) .NE. 0_8 ) THEN
C            adjacency list of secondary variable was 
C            provided by the user, 
C            the space will be compressed if necessary
             LEN (I) = LEN(J)
           ELSE
             LEN (I) = 0
           ENDIF
           PE (I)   = int(-J,8)
           NV (J)   = NV (J) + NV (I)
           NV (I)   = 0
           ELEN (I) = 0
        ENDIF
      ENDDO
C     ----------------------------------------------------------------
C     initialize degree lists and eliminate rows with no off-diag. nz.
C     ----------------------------------------------------------------
      DO 20 I = 1, N
        DEG = DEGREE (I)
C       degree(i) < 0 corresponds to secondary variables
C       that need be skipped.
        IF (DEG .GT. 0) THEN
C         ----------------------------------------------------------
C         place i in the degree list corresponding to its degree
C         ----------------------------------------------------------
          INEXT = HEAD (DEG)
          IF (INEXT .NE. 0) LAST (INEXT) = I
          NEXT (I) = INEXT
          HEAD (DEG) = I
        ELSE IF ( DEG.EQ. 0) THEN
C         ----------------------------------------------------------
C         we have a variable that can be eliminated at once because
C         there is no off-diagonal non-zero in its row.
C         ----------------------------------------------------------
C
C         We have a graph of supervariable and thus need to update 
C         singleton that might already be supervariables with nv(i)
C         When a supervariable is eliminated its 
C        principal variable must be set to the current step
C         (NEL+1) which must be stored (negated) in ELEN  
C         ONLY THEN (current step) NEL should be incremented.
C         This will be exploited when computing the global ordering
C         of all (secondary and principal) variables at the end of the AMD routine.
          ELEN (I) = - (NEL + 1)
          NEL = NEL + NV(I)
          PE (I) = 0_8
          W (I) = 0
        ENDIF
   20 CONTINUE
C=======================================================================
C  WHILE (selecting pivots) DO
C=======================================================================
C
C  Note that we do want to loop until NEL = N since
C  we update NEL with the size of the eliminated supervariable
C  
   30 IF (NEL .LT. N) THEN
C=======================================================================
C  GET PIVOT OF MINIMUM DEGREE
C=======================================================================
C       -------------------------------------------------------------
C       find next supervariable for elimination
C       -------------------------------------------------------------
        DO 40 DEG = MINDEG, N
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
C       -------------------------------------------------------------
C       remove chosen variable from link list
C       -------------------------------------------------------------
        INEXT = NEXT (ME)
        IF (INEXT .NE. 0) LAST (INEXT) = 0
        HEAD (DEG) = INEXT
C       -------------------------------------------------------------
C       me represents the elimination of pivots nel+1 to nel+nv(me).
C       place me itself as the first in this set.  It will be moved
C       to the nel+nv(me) position when the permutation vectors are
C       computed.
C       -------------------------------------------------------------
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
C=======================================================================
C  CONSTRUCT NEW ELEMENT
C=======================================================================
C       -------------------------------------------------------------
C       At this point, me is the pivotal supervariable.  It will be
C       converted into the current element.  Scan list of the
C       pivotal supervariable, me, setting tree pointers and
C       constructing new list of supervariables for the new element,
C       me.  p is a pointer to the current position in the old list.
C       -------------------------------------------------------------
C       flag the variable "me" as being in Lme by negating nv (me)
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
C         ----------------------------------------------------------
C         construct the new element in place
C         ----------------------------------------------------------
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + int(LEN (ME) - 1,8)
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
C             ----------------------------------------------------
C             i is a principal variable not yet placed in Lme.
C             store i in new list
C             ----------------------------------------------------
              DEGME = DEGME + NVI
C             flag i as being in Lme by negating nv (i)
              NV (I) = -NVI
              PME2 = PME2 + 1_8
              IW (PME2) = I
C             ----------------------------------------------------
C             remove variable i from degree list.
C             ----------------------------------------------------
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
C               i is at the head of the degree list
                HEAD (DEGREE (I)) = INEXT
              ENDIF
            ENDIF
   60     CONTINUE
C         this element takes no new memory in iw:
          NEWMEM = 0
        ELSE
C         ----------------------------------------------------------
C         construct the new element in empty space, iw (pfree ...)
C         ----------------------------------------------------------
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          KNT1_UPDATED = 0
          DO 120 KNT1 = 1, ELENME + 1
            KNT1_UPDATED = KNT1_UPDATED +1
            IF (KNT1 .GT. ELENME) THEN
C             search the supervariables in me.
              E = ME
              PJ = P
              LN = SLENME
            ELSE
C             search the elements in me.
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
C           -------------------------------------------------------
C           search for different supervariables and add them to the
C           new list, compressing when necessary. this loop is
C           executed once for each element in the list and once for
C           all the supervariables in the list.
C           -------------------------------------------------------
            KNT2_UPDATED = 0
            DO 110 KNT2 = 1, LN
              KNT2_UPDATED = KNT2_UPDATED+1
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
C               -------------------------------------------------
C               compress iw, if necessary
C               -------------------------------------------------
                IF (PFREE .GT. IWLEN) THEN
C                 prepare for compressing iw by adjusting
C                 pointers and lengths so that the lists being
C                 searched in the inner and outer loops contain
C                 only the remaining entries.
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1_UPDATED
C                 Reset KNT1_UPDATED in case of recompress 
C                 at same iteration of the loop 120
                  KNT1_UPDATED = 0
C                 Check if anything left in supervariable ME
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0_8
                  PE (E) = PJ
                  LEN (E) = LN - KNT2_UPDATED
C                 Reset KNT2_UPDATED in case of recompress 
C                 at same iteration of the loop 110
                  KNT2_UPDATED = 0
C                 Check if anything left in element E
                  IF (LEN (E) .EQ. 0) PE (E) = 0_8
                  NCMPA = NCMPA + 1
C                 store first item in pe
C                 set first entry to -item
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0_8) THEN
                      PE (J) = int(IW (PN),8)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
C                 psrc/pdst point to source/destination
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
C                 while loop:
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
C                   search for next negative entry
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = int(PE (J))
                      PE (J) = PDST
                      PDST = PDST + 1_8
C                     copy from source to destination
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + int(LENJ - 1,8)
                      PSRC = PSRC + int(LENJ - 1,8)
                    ENDIF
                    GO TO 80
                  ENDIF
C                 move the new partially-constructed element
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
C               -------------------------------------------------
C               i is a principal variable not yet placed in Lme
C               store i in new list
C               -------------------------------------------------
                DEGME = DEGME + NVI
C               flag i as being in Lme by negating nv (i)
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
C               -------------------------------------------------
C               remove variable i from degree link list
C               -------------------------------------------------
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
C                 i is at the head of the degree list
                  HEAD (DEGREE (I)) = INEXT
                ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
C             set tree pointer and flag to indicate element e is
C             absorbed into new element me (the parent of e is me)
              PE (E) = int(-ME,8)
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
C         this element takes newmem new memory in iw (possibly zero)
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = max (MAXMEM, MEM)
        ENDIF
C       -------------------------------------------------------------
C       me has now been converted into an element in iw (pme1..pme2)
C       -------------------------------------------------------------
C       degme holds the external degree of new element
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = int(PME2 - PME1 + 1)
C       -------------------------------------------------------------
C       make sure that wflg is not too large.  With the current
C       value of wflg, wflg+n must not cause integer overflow
C       -------------------------------------------------------------
        IF (WFLG .GT. MAXINT_N) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
C=======================================================================
C  COMPUTE (w (e) - wflg) = |Le\Lme| FOR ALL ELEMENTS
C=======================================================================
C       -------------------------------------------------------------
C       Scan 1:  compute the external degrees of previous elements
C       with respect to the current element.  That is:
C            (w (e) - wflg) = |Le \ Lme|
C       for each element e that appears in any supervariable in Lme.
C       The notation Le refers to the pattern (list of
C       supervariables) of a previous element e, where e is not yet
C       absorbed, stored in iw (pe (e) + 1 ... pe (e) + iw (pe (e))).
C       The notation Lme refers to the pattern of the current element
C       (stored in iw (pme1..pme2)).   If (w (e) - wflg) becomes
C       zero, then the element e will be absorbed in scan 2.
C       -------------------------------------------------------------
        DO 150 PME = PME1, PME2
          I = IW (PME)
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
C           note that nv (i) has been negated to denote i in Lme:
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + int(ELN - 1,8)
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
C               unabsorbed element e has been seen in this loop
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
C               e is an unabsorbed element
C               this is the first we have seen e in all of Scan 1
                WE = DEGREE (E) + WNVI
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
C=======================================================================
C  DEGREE UPDATE AND ELEMENT ABSORPTION
C=======================================================================
C       -------------------------------------------------------------
C       Scan 2:  for each i in Lme, sum up the degree of Lme (which
C       is degme), plus the sum of the external degrees of each Le
C       for the elements e appearing within i, plus the
C       supervariables in i.  Place i in hash list.
C       -------------------------------------------------------------
        DO 180 PME = PME1, PME2
          I = IW (PME)
          P1 = PE (I)
          P2 = P1 + int(ELEN (I) - 1,8)
          PN = P1
          HASH = 0_8
          DEG = 0
C         ----------------------------------------------------------
C         scan the element list associated with supervariable i
C         ----------------------------------------------------------
          DO 160 P = P1, P2
            E = IW (P)
C           dext = | Le \ Lme |
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
            ELSE IF (DEXT .EQ. 0) THEN
#if defined (NOAGG2)
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
#else
C             aggressive absorption: e is not adjacent to me, but
C             the |Le \ Lme| is 0, so absorb it into me
              PE (E) = int(-ME,8)
              W (E) = 0
#endif
            ENDIF
  160     CONTINUE
C         count the number of elements in i (including me):
          ELEN (I) = int(PN - P1 + 1_8)
C         ----------------------------------------------------------
C         scan the supervariables in the list associated with i
C         ----------------------------------------------------------
          P3 = PN
          DO 170 P = P2 + 1, P1 + int(LEN (I) - 1,8)
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
C             j is unabsorbed, and not in Lme.
C             add to degree and add to new list
              DEG = DEG + NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + int(J,kind=8)
            ENDIF
  170     CONTINUE
C         ----------------------------------------------------------
C         update the degree and check for mass elimination
C         ----------------------------------------------------------
#if defined (NOAGG2)
          IF (DEG.EQ.0.AND.(ELEN(I).GT.1)) THEN
C         When DEG is zero we need to 
C         absorb in ME all elements adjacent to I 
                   P1 = PE (I)
C                  exclude ME --> -2 
                   P2 = P1 + int(ELEN (I),8) - 2_8
                   DO P =P1,P2  
                     E      = IW(P)
                     PE (E) = int(-ME,8)
                     W (E)  = 0
                   ENDDO
          ENDIF
C              .... Ready for mass elimination
#endif
          IF (DEG .EQ. 0) THEN
C           -------------------------------------------------------
C           mass elimination
C           -------------------------------------------------------
C           There is nothing left of this node except for an
C           edge to the current pivot element.  elen (i) is 1,
C           and there are no variables adjacent to node i.
C           Absorb i into the current pivot element, me.
            PE (I) = int(-ME,8)
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
C           -------------------------------------------------------
C           update the upper-bound degree of i
C           -------------------------------------------------------
C           the following degree does not yet include the size
C           of the current element, which is added later:
            DEGREE (I) = min (DEGREE (I), DEG)
C           -------------------------------------------------------
C           add me to the list for i
C           -------------------------------------------------------
C           move first supervariable to end of list
            IW (PN) = IW (P3)
C           move first element to end of element part of list
            IW (P3) = IW (P1)
C           add new element to front of list.
            IW (P1) = ME
C           store the new length of the list in len (i)
            LEN (I) = int(PN - P1 + 1_8)
C           -------------------------------------------------------
C           place in hash bucket.  Save hash key of i in last (i).
C           -------------------------------------------------------
            HASH = mod (HASH, HMOD) + 1_8
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
C             the degree list is empty, hash head is -j
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
C             degree list is not empty
C             use last (head (hash)) as hash head
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = int(HASH,kind=kind(LAST))
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
C       -------------------------------------------------------------
C       Clear the counter array, w (...), by incrementing wflg.
C       -------------------------------------------------------------
        DMAX = max (DMAX, DEGME)
        WFLG = WFLG + DMAX
C       make sure that wflg+n does not cause integer overflow
        IF (WFLG .GT. MAXINT_N) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
C       at this point, w (1..n) .lt. wflg holds
C=======================================================================
C  SUPERVARIABLE DETECTION
C=======================================================================
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF (NV (I) .LT. 0) THEN
C           i is a principal variable in Lme
C           -------------------------------------------------------
C           examine all hash buckets with 2 or more variables.  We
C           do this by examing all unique hash keys for super-
C           variables in the pattern Lme of the current element, me
C           -------------------------------------------------------
            HASH = int(LAST (I),kind=8)
C           let i = head of hash bucket, and empty the hash bucket
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
C             degree list is empty
              I = -J
              HEAD (HASH) = 0
            ELSE
C             degree list is not empty, restore last () of head
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
C           while loop:
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
C             ----------------------------------------------------
C             this bucket has one or more variables following i.
C             scan all of them to see if i can absorb any entries
C             that follow i in hash bucket.  Scatter i into w.
C             ----------------------------------------------------
              LN = LEN (I)
              ELN = ELEN (I)
C             do not flag the first element in the list (me)
              DO 210 P = PE (I) + 1, PE (I) + int(LN - 1,8)
                W (IW (P)) = WFLG
  210         CONTINUE
C             ----------------------------------------------------
C             scan every other entry j following i in bucket
C             ----------------------------------------------------
              JLAST = I
              J = NEXT (I)
C             while loop:
  220         CONTINUE
              IF (J .NE. 0) THEN
C               -------------------------------------------------
C               check if j and i have identical nonzero pattern
C               -------------------------------------------------
C               jump if i and j do not have same size data structure
                IF (LEN (J) .NE. LN) GO TO 240
C               jump if i and j do not have same number adj elts
                IF (ELEN (J) .NE. ELN) GO TO 240
C               do not flag the first element in the list (me)
                DO 230 P = PE (J) + 1, PE (J) + int(LN - 1,8)
C                 jump if an entry (iw(p)) is in j but not in i
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
C               -------------------------------------------------
C               found it!  j can be absorbed into i
C               -------------------------------------------------
                PE (J) = int(-I,8)
C               both nv (i) and nv (j) are negated since they
C               are in Lme, and the absolute values of each
C               are the number of variables in i and j:
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
C               delete j from hash bucket
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
C               -------------------------------------------------
  240           CONTINUE
C               j cannot be absorbed into i
C               -------------------------------------------------
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
C             ----------------------------------------------------
C             no more variables can be absorbed into i
C             go to next i in bucket and clear flag array
C             ----------------------------------------------------
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
C=======================================================================
C  RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVAR. FROM ELEMENT
C=======================================================================
        P = PME1
        NLEFT = N - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
C           i is a principal variable in Lme
C           restore nv (i) to signify that i is principal
            NV (I) = NVI
C           -------------------------------------------------------
C           compute the external degree (add size of current elem)
C           -------------------------------------------------------
            DEG = min (DEGREE (I) + DEGME - NVI, NLEFT - NVI)
C           -------------------------------------------------------
C           place the supervariable at the head of the degree list
C           -------------------------------------------------------
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (DEG) = I
C           -------------------------------------------------------
C           save the new degree, and find the minimum degree
C           -------------------------------------------------------
            MINDEG = min (MINDEG, DEG)
            DEGREE (I) = DEG
C           -------------------------------------------------------
C           place the supervariable in the element pattern
C           -------------------------------------------------------
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
C=======================================================================
C  FINALIZE THE NEW ELEMENT
C=======================================================================
        NV (ME) = NVPIV + DEGME
C       nv (me) is now the degree of pivot (including diagonal part)
C       save the length of the list for the new element me
        LEN (ME) = int(P - PME1)
        IF (LEN (ME) .EQ. 0) THEN
C         there is nothing left of the current pivot element
          PE (ME) = 0_8
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
C         element was not constructed in place: deallocate part
C         of it (final size is less than or equal to newmem,
C         since newly nonprincipal variables have been removed).
          PFREE = P
          MEM = MEM - NEWMEM + int(LEN (ME),8)
        ENDIF
C=======================================================================
C       END WHILE (selecting pivots)
      GO TO 30
      ENDIF
C=======================================================================
C=======================================================================
C  COMPUTE THE PERMUTATION VECTORS
C=======================================================================
C     ----------------------------------------------------------------
C     The time taken by the following code is O(n).  At this
C     point, elen (e) = -k has been done for all elements e,
C     and elen (i) = 0 has been done for all nonprincipal
C     variables i.  At this point, there are no principal
C     supervariables left, and all elements are absorbed.
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
C     compute the ordering of unordered nonprincipal variables
C     ----------------------------------------------------------------
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
C         ----------------------------------------------------------
C         i is an un-ordered row.  Traverse the tree from i until
C         reaching an element, e.  The element, e, was the
C         principal supervariable of i and all nodes in the path
C         from i to when e was selected as pivot.
C         ----------------------------------------------------------
          J = int(-PE (I))
C         while (j is a variable) do:
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = int(-PE (J))
              GO TO 270
            ENDIF
            E = J
C           ----------------------------------------------------------
C           get the current pivot ordering of e
C           ----------------------------------------------------------
            K = -ELEN (E)
C           ----------------------------------------------------------
C           traverse the path again from i to e, and compress the
C           path (all nodes point to e).  Path compression allows
C           this code to compute in O(n) time.  Order the unordered
C           nodes in the path, and place the element e at the end.
C           ----------------------------------------------------------
            J = I
C           while (j is a variable) do:
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = int(-PE (J))
              PE (J) = int(-E,8)
              IF (ELEN (J) .EQ. 0) THEN
C               j is an unordered row
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
C         leave elen (e) negative, so we know it is an element
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
C     ----------------------------------------------------------------
C     reset the inverse permutation (elen (1..n)) to be positive,
C     and compute the permutation (last (1..n)).
C     ----------------------------------------------------------------
      DO 300 I = 1, N
        K = abs (ELEN (I))
        LAST (K) = I
        ELEN (I) = K
  300 CONTINUE
C=======================================================================
C  RETURN THE MEMORY USAGE IN IW
C=======================================================================
C     If maxmem is less than or equal to iwlen, then no compressions
C     occurred, and iw (maxmem+1 ... iwlen) was unused.  Otherwise
C     compressions did occur, and iwlen would have had to have been
C     greater than or equal to maxmem for no compressions to occur.
C     Return the value of maxmem in the pfree argument.
      PFREE = MAXMEM
C===============================
C     Save PE in PARENT array
      DO I=1,N
       PARENT(I) = int(PE(I))
      ENDDO
C===============================
      RETURN
      END SUBROUTINE MUMPS_AMD_ELT
C ----------------------------------------------------------------------
C Description of MUMPS_HAMD:
C   MUMPS_HAMD is a modification of AMD reference code (MUMPS_ANA_H) 
C   designed to take into account a halo in the graph. 
C   The graph is composed is partitioned in two types of nodes
C   the so called internal nodes and the so called halo nodes.
C   Halo nodes cannot be selected the both the inital degrees 
C   and updated degrees of internal node should be taken 
C   into account.  
C   This routine also referred to as HALOAMD in MUMPS comments
C   is used for both Schur functionality and in the coupling with 
C   partitioners such as SCOTCH.
C
C   Restrictive integer 64 bit variant :
C   it is assumed that IW array size can exceed 32-bit integer
C
C      
      SUBROUTINE MUMPS_HAMD(N, IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
     &                   LAST, NCMPA, DEGREE, HEAD, NEXT, W, PARENT,
     &                   LISTVAR_SCHUR, SIZE_SCHUR)
C
C Parameters
C    Input not modified
      INTEGER, intent(in) :: SIZE_SCHUR
      INTEGER, intent(in) :: LISTVAR_SCHUR(SIZE_SCHUR)
      INTEGER, INTENT(IN)    :: N
      INTEGER(8), INTENT(IN) :: IWLEN
C     Input undefined on output 
      INTEGER, INTENT(INOUT)  :: LEN(N), IW(IWLEN)
C 
C     Output only 
      INTEGER, INTENT(OUT)   :: NCMPA
      INTEGER, INTENT(OUT)   :: NV(N), ELEN(N), LAST(N), PARENT(N)
C 
C     Input/output
      INTEGER(8), INTENT(INOUT) :: PFREE
      INTEGER(8), INTENT(INOUT) :: PE(N)
C 
C     Internal Workspace only
      INTEGER     :: NEXT(N), DEGREE(N), HEAD(N), W(N)
C
C ---------------------
C Interface Description
C ---------------------
C    HAMD (short for HALOAMD)
C    The initial version (so called HALOAMD_V1, developped in September 1997)
C    is designed to experiment the numerical (fill-in) impact 
C    of taking into account the halo. This code should be able 
C    to experiment no-halo, partial halo, complete halo.
C    DATE: September 17th 1997
C
C    HALOAMD is designed to process a gragh composed of two types
C            of nodes, V0 and V1, extracted from a larger gragh. 
C            V0^V1 = {}, 
C       
C            We used Min. degree heuristic to order only 
C            nodes in V0, but the adjacency to nodes
C            in V1 is taken into account during ordering.
C            Nodes in V1 are odered at last.
C            Adjacency between nodes of V1 need not be provided,
C            however |len(i)| must always corresponds to the number of 
C            edges effectively provided in the adjacency list of i.
C          On input :
C          ********
C            Nodes INODE in V1 are flagged with len(INODE) = -degree 
C            modif version HALO V3 (August 1998): 
C                                 if len(i) =0 and i \in V1 then 
C                                 len(i) must be set on input to -N-1
C          ERROR return (negative values in ncmpa)
C          ************
C            negative value in ncmpa indicates an error detected 
C               by HALOAMD.
C
C            The graph provided MUST follow the rule:
C             if (i,j) is an edge in the gragh then 
C             j must be in the adjacency list of i AND 
C             i must be in the adjacency list of j.
C    REMARKS:
C        1/  Providing edges between nodes of V1 should not 
C            affect the final ordering, only the amount of edges 
C            of the halo should effectively affect the solution.
C            This code should work in the following cases:
C              1/ halo not provided
C              2/ halo partially provided
C              3/ complete halo
C              4/ complete halo+interconnection between nodes of V1.
C
C              1/ should run and provide identical results (w.r.t to current 
C               implementation of AMD in SCOTCH).
C             3/ and 4 should provide identical results.
C
C        2/ All modifications of the AMD initial code are indicated
C           with begin HALO .. end HALO
C
C            
C   Ordering of nodes in V0 is based on 
C   Approximate Minimum Degree ordering algorithm, 
C   with aggressive absorption:
C   Given a representation of the nonzero pattern of a symmetric matrix,
C   A, (excluding the diagonal) perform an approximate minimum
C   degree ordering to compute a pivot order
C   such that fill-in in the Cholesky factors A = LL^T is kept low. 
C
C ------------------------------
C Modification history:
C ---------------------
C Date:  September, 1997 (V1)
C        April, 1998 (V2)
C        August, 1998 (V3)
C        Octobre, 1998 (V4)
C        December, 1998 (V5)
C        January, 1999 (V6)
C   HALOAMD_V6: 
C   ----------
C    1/ ERROR 2 detection followed by stop statement suppressed 
C  . 2/ pb 1  identified in V5 was not correctly solved
C
C   HALOAMD_V5:
C   ----------
C    1/ Pb with matrix psmigr 1, because upper bound 
C      degree DEG  >N was considered as a node in V1 
C
C   HALOAMD_V4: 
C   ---------- 
C    Only UnsymetrizedMultifrontal interface 
C    (ok for both scotch and UnsymetricMultifrontal) is 
C    included in this file
C
C   HALOAMD_V3: 
C   ----------
C    Problem in version 2 : variables of V1 with len(i) =0 
C                        are not well processed. 
C                        See modification of the 
C                        input to characterize those variables.
C 
C    Problem detected by Jacko Koster while experimenting with 
C    version 2 of haloAMD in the context of multiple front method :
C    "if for an interface variable i, row i in the matrix has only a
C    nonzero entry on the diagonal, we first remove this entry and len(i)  
C    is set to zero on input to HALOAMD. However, this means that HALOAMD 
C    will treat variable i as an interior variable (in V0) instead as an 
C    interface variable (in V1). (It is indeed a bit strange to have such 
C    interface variables but we encountered some in our debugging 
C    experiments with some random partitionings.)
C
C    Solution : 
C     IF on input i \in V1 and len(i) =0 (that is adjlist(i)={}) THEN
C      len(i) must be set on input to -N-1. 
C     ENDIF
C     therefore all variables i / len(i) < 0 an only those are in V1
C     variable with len(i) = -N-1 are then processed differently at 
C     the beginning of the code
C
C   HALOAMD_V2:
C   ----------
C    The end of the tree (including links to block of flagged indices
C                       is built) . The list of flagged indices is 
C                       considered as a dense amalgamated node.
C
C    Comments on the OUTPUT:
C    ----------------------
C    Let V= V0 U V1 the nodes of the initial graph (|V|=n). 
C    The assembly tree corresponds to the tree 
C    of the supernodes (or supervariables). Each node of the 
C    assembly tree is then composed of one principal variable 
C    and a list of secondary variables. The list of 
C    variable of a node (principal + secondary variables) then 
C    describes the structure of the diagonal bloc of the 
C    supernode. 
C    The elimination tree denotes the tree of all the variables(=node) and 
C    is therefore of order n.
C
C    The arrays NV(N) and PE(N) give a description of the 
C    assembly tree. 
C     1/ Description of array nv(N) (on OUPUT)
C        nv(i)=0 i is a secondary variable 
C        N+1> nv(i) >0 i is a principal variable, nv(i) holds the 
C                the number of elements in column i of L (true degree of i)
C     2/ Description of array PE(N) (on OUPUT)
C       pe(i) = -(father of variable/node i) in the elimination tree:
C       If nv (i) .gt. 0, then i represents a node in the assembly tree,
C       and the parent of i is -pe (i), or zero if i is a root.
C       If nv (i) = 0, then (i,-pe (i)) represents an edge in a
C       subtree, the root of which is a node in the assembly tree.
C     3/ Example:
C       Let If be a root node father of Is in the assembly tree. 
C       If is the principal 
C       variable of the node If and let If1, If2, If3 be the 
C       secondary variables of node If.
C       Is is the principal 
C       variable of the node Is and let Is1, Is2 be the secondary variables
C       of node Is.
C      
C       THEN: 
C        NV(If1)=NV(If2)=NV(If3) = 0  (secondary variables)
C        NV(Is1)=NV(Is2) = 0  (secondary variables)
C        NV(If) > 0  ( principal variable)
C        NV(Is) > 0  ( principal variable)
C        PE(If)  = 0 (root node)
C        PE(Is)  = -If (If is the father of Is in the assembly tree)
C        PE(If1)=PE(If2)=PE(If3)= -If  ( If is the principal variable)
C        PE(Is1)=PE(Is2)= -Is  ( Is is the principal variable)
C-----------------------------------------------------------------------
C INPUT ARGUMENTS (unaltered):
C-----------------------------------------------------------------------
C n:    The matrix order.
C
C       Restriction:  n .ge. 1
C iwlen:        The length of iw (1..iwlen).  On input, the matrix is
C       stored in iw (1..pfree-1).  However, iw (1..iwlen) should be
C       slightly larger than what is required to hold the matrix, at
C       least iwlen .ge. pfree + n is recommended.  Otherwise,
C       excessive compressions will take place.
C       *** We do not recommend running this algorithm with ***
C       ***      iwlen .lt. pfree + n.                      ***
C       *** Better performance will be obtained if          ***
C       ***      iwlen .ge. pfree + n                       ***
C       *** or better yet                                   ***
C       ***      iwlen .gt. 1.2 * pfree                     ***
C       *** (where pfree is its value on input).            ***
C       The algorithm will not run at all if iwlen .lt. pfree-1.
C
C       Restriction: iwlen .ge. pfree-1
C-----------------------------------------------------------------------
C INPUT/OUPUT ARGUMENTS:
C-----------------------------------------------------------------------
C pe:   On input, pe (i) is the index in iw of the start of row i, or
C       zero if row i has no off-diagonal non-zeros.
C
C       During execution, it is used for both supervariables and
C       elements:
C
C       * Principal supervariable i:  index into iw of the
C               description of supervariable i.  A supervariable
C               represents one or more rows of the matrix
C               with identical nonzero pattern.
C       * Non-principal supervariable i:  if i has been absorbed
C               into another supervariable j, then pe (i) = -j.
C               That is, j has the same pattern as i.
C               Note that j might later be absorbed into another
C               supervariable j2, in which case pe (i) is still -j,
C               and pe (j) = -j2.
C       * Unabsorbed element e:  the index into iw of the description
C               of element e, if e has not yet been absorbed by a
C               subsequent element.  Element e is created when
C               the supervariable of the same name is selected as
C               the pivot.
C       * Absorbed element e:  if element e is absorbed into element
C               e2, then pe (e) = -e2.  This occurs when the pattern of
C               e (that is, Le) is found to be a subset of the pattern
C               of e2 (that is, Le2).  If element e is "null" (it has
C               no nonzeros outside its pivot block), then pe (e) = 0.
C
C       On output, pe holds the assembly tree/forest, which implicitly
C       represents a pivot order with identical fill-in as the actual
C       order (via a depth-first search of the tree).
C
C       On output:
C       If nv (i) .gt. 0, then i represents a node in the assembly tree,
C       and the parent of i is -pe (i), or zero if i is a root.
C       If nv (i) = 0, then (i,-pe (i)) represents an edge in a
C       subtree, the root of which is a node in the assembly tree.
C       On output:  (PE is copied on output into PARENT array)
C
C pfree:        On input, the matrix is stored in iw (1..pfree-1) and
C       the rest of the array iw is free.
C       During execution, additional data is placed in iw, and pfree
C       is modified so that components  of iw from pfree are free.
C       On output, pfree is set equal to the size of iw that
C       would have been needed for no compressions to occur.  If
C       ncmpa is zero, then pfree (on output) is less than or equal to
C       iwlen, and the space iw (pfree+1 ... iwlen) was not used.
C       Otherwise, pfree (on output) is greater than iwlen, and all the
C       memory in iw was used.
C-----------------------------------------------------------------------
C INPUT/MODIFIED (undefined on output):
C-----------------------------------------------------------------------
C len:  On input, len (i) 
C           positive or null (>=0) : i \in V0 and 
C                     len(i) holds the number of entries in row i of the
C                     matrix, excluding the diagonal.  
C           negative (<0) : i \in V1, and 
C                     -len(i) hold the number of entries in row i of the
C                     matrix, excluding the diagonal.
C                     len(i) = - | Adj(i) | if i \in V1                    
C                              or -N -1 if  | Adj(i) | = 0 and i \in V1 
C       The contents of len (1..n)
C       are undefined on output.
C iw:   On input, iw (1..pfree-1) holds the description of each row i
C       in the matrix.  The matrix must be symmetric, and both upper
C       and lower triangular parts must be present.  The diagonal must
C       not be present.  Row i is held as follows:
C
C               len (i):  the length of the row i data structure
C               iw (pe (i) ... pe (i) + len (i) - 1):
C                       the list of column indices for nonzeros
C                       in row i (simple supervariables), excluding
C                       the diagonal.  All supervariables start with
C                       one row/column each (supervariable i is just
C                       row i).
C               if len (i) is zero on input, then pe (i) is ignored
C               on input.
C
C               Note that the rows need not be in any particular order,
C               and there may be empty space between the rows.
C
C       During execution, the supervariable i experiences fill-in.
C       This is represented by placing in i a list of the elements
C       that cause fill-in in supervariable i:
C
C               len (i):  the length of supervariable i
C               iw (pe (i) ... pe (i) + elen (i) - 1):
C                       the list of elements that contain i.  This list
C                       is kept short by removing absorbed elements.
C               iw (pe (i) + elen (i) ... pe (i) + len (i) - 1):
C                       the list of supervariables in i.  This list
C                       is kept short by removing nonprincipal
C                       variables, and any entry j that is also
C                       contained in at least one of the elements
C                       (j in Le) in the list for i (e in row i).
C
C       When supervariable i is selected as pivot, we create an
C       element e of the same name (e=i):
C
C               len (e):  the length of element e
C               iw (pe (e) ... pe (e) + len (e) - 1):
C                       the list of supervariables in element e.
C
C       An element represents the fill-in that occurs when supervariable
C       i is selected as pivot (which represents the selection of row i
C       and all non-principal variables whose principal variable is i).
C       We use the term Le to denote the set of all supervariables
C       in element e.  Absorbed supervariables and elements are pruned
C       from these lists when computationally convenient.
C
C       CAUTION:  THE INPUT MATRIX IS OVERWRITTEN DURING COMPUTATION.
C       The contents of iw are undefined on output.
C-----------------------------------------------------------------------
C OUTPUT (need not be set on input):
C-----------------------------------------------------------------------
C nv:   During execution, abs (nv (i)) is equal to the number of rows
C       that are represented by the principal supervariable i.  If i is
C       a nonprincipal variable, then nv (i) = 0.  Initially,
C       nv (i) = 1 for all i.  nv (i) .lt. 0 signifies that i is a
C       principal variable in the pattern Lme of the current pivot
C       element me.  On output, nv (e) holds the true degree of element
C       e at the time it was created (including the diagonal part).
C begin HALO
C       On output, nv(I) can be used to find node in set V1.
C       nv(I) = N+1 characterizes nodes in V1. 
C end HALO
C elen: See the description of iw above.  At the start of execution,
C       elen (i) is set to zero.  During execution, elen (i) is the
C       number of elements in the list for supervariable i.  When e
C       becomes an element, elen (e) = -nel is set, where nel is the
C       current step of factorization.  elen (i) = 0 is done when i
C       becomes nonprincipal.
C
C       For variables, elen (i) .ge. 0 holds until just before the
C       permutation vectors are computed.  For elements,
C       elen (e) .lt. 0 holds.
C
C       On output elen (1..n) holds the inverse permutation (the same
C       as the 'INVP' argument in Sparspak).  That is, if k = elen (i),
C       then row i is the kth pivot row.  Row i of A appears as the
C       (elen(i))-th row in the permuted matrix, PAP^T.
C last: In a degree list, last (i) is the supervariable preceding i,
C       or zero if i is the head of the list.  In a hash bucket,
C       last (i) is the hash key for i.  last (head (hash)) is also
C       used as the head of a hash bucket if head (hash) contains a
C       degree list (see head, below).
C
C       On output, last (1..n) holds the permutation (the same as the
C       'PERM' argument in Sparspak).  That is, if i = last (k), then
C       row i is the kth pivot row.  Row last (k) of A is the k-th row
C       in the permuted matrix, PAP^T.
C ncmpa:        The number of times iw was compressed.  If this is
C       excessive, then the execution took longer than what could have
C       been.  To reduce ncmpa, try increasing iwlen to be 10% or 20%
C       larger than the value of pfree on input (or at least
C       iwlen .ge. pfree + n).  The fastest performance will be
C       obtained when ncmpa is returned as zero.  If iwlen is set to
C       the value returned by pfree on *output*, then no compressions
C       will occur.
C begin HALO
C        on output ncmpa <0 --> error detected during HALO_AMD:
C           error 1: ncmpa = -N , ordering was stopped.
C end HALO
C
C-----------------------------------------------------------------------
C LOCAL (not input or output - used only during execution):
C-----------------------------------------------------------------------
C degree:       If i is a supervariable, then degree (i) holds the
C       current approximation of the external degree of row i (an upper
C       bound).  The external degree is the number of nonzeros in row i,
C       minus abs (nv (i)) (the diagonal part).  The bound is equal to
C       the external degree if elen (i) is less than or equal to two.
C       We also use the term "external degree" for elements e to refer
C       to |Le \ Lme|.  If e is an element, then degree (e) holds |Le|,
C       which is the degree of the off-diagonal part of the element e
C       (not including the diagonal part).
C begin HALO
C       degree(I) = n+1 indicates that i belongs to V1
C end HALO
C
C head: head is used for degree lists.  head (deg) is the first
C       supervariable in a degree list (all supervariables i in a
C       degree list deg have the same approximate degree, namely,
C       deg = degree (i)).  If the list deg is empty then
C       head (deg) = 0.
C
C       During supervariable detection head (hash) also serves as a
C       pointer to a hash bucket.
C       If head (hash) .gt. 0, there is a degree list of degree hash.
C               The hash bucket head pointer is last (head (hash)).
C       If head (hash) = 0, then the degree list and hash bucket are
C               both empty.
C       If head (hash) .lt. 0, then the degree list is empty, and
C               -head (hash) is the head of the hash bucket.
C       After supervariable detection is complete, all hash buckets
C       are empty, and the (last (head (hash)) = 0) condition is
C       restored for the non-empty degree lists.
C next: next (i) is the supervariable following i in a link list, or
C       zero if i is the last in the list.  Used for two kinds of
C       lists:  degree lists and hash buckets (a supervariable can be
C       in only one kind of list at a time).
C w:    The flag array w determines the status of elements and
C       variables, and the external degree of elements.
C
C       for elements:
C          if w (e) = 0, then the element e is absorbed
C          if w (e) .ge. wflg, then w (e) - wflg is the size of
C               the set |Le \ Lme|, in terms of nonzeros (the
C               sum of abs (nv (i)) for each principal variable i that
C               is both in the pattern of element e and NOT in the
C               pattern of the current pivot element, me).
C          if wflg .gt. w (e) .gt. 0, then e is not absorbed and has
C               not yet been seen in the scan of the element lists in
C               the computation of |Le\Lme| in loop 150 below.
C
C       for variables:
C          during supervariable detection, if w (j) .ne. wflg then j is
C          not in the pattern of variable i
C
C       The w array is initialized by setting w (i) = 1 for all i,
C       and by setting wflg = 2.  It is reinitialized if wflg becomes
C       too large (to ensure that wflg+n does not cause integer
C       overflow).
C-----------------------------------------------------------------------
C LOCAL INTEGERS:
C-----------------------------------------------------------------------
      INTEGER    :: DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     &        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     &        LENJ, LN, ME, MINDEG, NEL, 
     &        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X,
     &        NBFLAG, NREAL, LASTD, NELME
      INTEGER KNT1_UPDATED, KNT2_UPDATED
      INTEGER(8) ::  MAXMEM, MEM, NEWMEM
      INTEGER    :: MAXINT_N
      INTEGER(8) :: HASH, HMOD
C deg:        the degree of a variable or element
C degme:      size, |Lme|, of the current element, me (= degree (me))
C dext:       external degree, |Le \ Lme|, of some element e
C dmax:       largest |Le| seen so far
C e:          an element
C elenme:     the length, elen (me), of element list of pivotal var.
C eln:        the length, elen (...), of an element list
C hash:       the computed value of the hash function
C hmod:       the hash function is computed modulo hmod = max (1,n-1)
C i:          a supervariable
C ilast:      the entry in a link list preceding i
C inext:      the entry in a link list following i
C j:          a supervariable
C jlast:      the entry in a link list preceding j
C jnext:      the entry in a link list, or path, following j
C k:          the pivot order of an element or variable
C knt1:       loop counter used during element construction
C knt2:       loop counter used during element construction
C knt3:       loop counter used during compression
C lenj:       len (j)
C ln:         length of a supervariable list
C maxint_n:   large integer to test risk of overflow on wflg
C maxmem:     amount of memory needed for no compressions
C me:         current supervariable being eliminated, and the
C                     current element created by eliminating that
C                     supervariable
C mem:        memory in use assuming no compressions have occurred
C mindeg:     current minimum degree
C nel:        number of pivots selected so far
C newmem:     amount of new memory needed for current pivot element
C nleft:      n - nel, the number of nonpivotal rows/columns remaining
C nvi:        the number of variables in a supervariable i (= nv (i))
C nvj:        the number of variables in a supervariable j (= nv (j))
C nvpiv:      number of pivots in current element
C slenme:     number of variables in variable list of pivotal variable
C we:         w (e)
C wflg:       used for flagging the w array.  See description of iw.
C wnvi:       wflg - nv (i)
C x:          either a supervariable or an element
C begin HALO
C nbflag:     number of flagged entries in the initial gragh.
C nreal :     number of entries on which ordering must be perfomed
C             (nreel = N- nbflag)
C nelme number of pivots selected when reaching the root
C lastd index of the last row in the list of dense rows
C end HALO
C-----------------------------------------------------------------------
C LOCAL POINTERS:
C-----------------------------------------------------------------------
      INTEGER(8) :: P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, 
     &              PN, PSRC
C             Any parameter (pe (...) or pfree) or local variable
C             starting with "p" (for Pointer) is an index into iw,
C             and all indices into iw use variables starting with
C             "p."  The only exception to this rule is the iwlen
C             input argument.
C p:          pointer into lots of things
C p1:         pe (i) for some variable i (start of element list)
C p2:         pe (i) + elen (i) -  1 for some var. i (end of el. list)
C p3:         index of first supervariable in clean list
C pdst:       destination pointer, for compression
C pend:       end of memory to compress
C pj:         pointer into an element or variable
C pme:        pointer into the current element (pme1...pme2)
C pme1:       the current element, me, is stored in iw (pme1...pme2)
C pme2:       the end of the current element
C pn:         pointer into a "clean" variable, also used to compress
C psrc:       source pointer, for compression
C-----------------------------------------------------------------------
C  FUNCTIONS CALLED:
C-----------------------------------------------------------------------
      INTRINSIC max, min, mod
C=======================================================================
C  INITIALIZATIONS
C=======================================================================
      WFLG = 2
      MAXINT_N=huge(WFLG)-N
      MINDEG = 1
      NCMPA = 0
      NEL = 0
      HMOD = int(max (1, N-1),kind=8)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
C begin HALO 
      NBFLAG = 0
      LASTD  = 0
C end HALO 
      DO 10 I = 1, N
        LAST (I) = 0
        HEAD (I) = 0
        NV (I) = 1
        W (I) = 1
        ELEN (I) = 0
        DEGREE(I) = LEN(I)
   10 CONTINUE
C
C begin HALO-SCHUR
      NBFLAG = SIZE_SCHUR
C
      DO K=1,SIZE_SCHUR
C
       I = LISTVAR_SCHUR(K)
       DEGREE(I) = N+1
       IF ((LEN(I) .EQ.0).OR.(LEN(I).EQ.-N-1)) THEN
C      Both ways of characterizing i \in Schur with Adj(I) = 0
C        Because of compress, we force skipping this
C        entry which is anyway empty
         PE (I)     = 0_8
         LEN(I)     = 0
       ENDIF
C      insert I at the end of degree list of n
C                  (safe: because max external degree is N-1)
       DEG = N
       IF (LASTD.EQ.0) THEN
C              degree list is empty
               LASTD     = I
               HEAD(DEG) = I
               NEXT(I)   = 0
               LAST(I)   = 0
       ELSE
               NEXT(LASTD) = I
               LAST(I)     = LASTD
               LASTD       = I
               NEXT(I)     = 0
       ENDIF
C
      ENDDO
C     number of entries to be ordered.
      NREAL = N - NBFLAG
C end HALO-SCHUR
C     ----------------------------------------------------------------
C     initialize degree lists and eliminate rows with no off-diag. nz.
C     ----------------------------------------------------------------
      DO 20 I = 1, N
        DEG = DEGREE (I)
C begin HALO-SCHUR
        IF (DEG.EQ.N+1)  GOTO 20
C end HALO-SCHUR
C
        IF (DEG .GT. 0) THEN
C         ----------------------------------------------------------
C         place i in the degree list corresponding to its degree
C         ----------------------------------------------------------
          INEXT = HEAD (DEG)
          IF (INEXT .NE. 0) LAST (INEXT) = I
          NEXT (I) = INEXT
          HEAD (DEG) = I
        ELSE
C         ----------------------------------------------------------
C         we have a variable that can be eliminated at once because
C         there is no off-diagonal non-zero in its row.
C         ----------------------------------------------------------
          NEL = NEL + NV(I)
          ELEN (I) = -NEL
          PE (I) = 0_8
          W (I) = 0
        ENDIF
   20 CONTINUE
C=======================================================================
C  WHILE (selecting pivots) DO
C=======================================================================
C begin HALO V5
      NLEFT = N-NEL
C end HALO V5
C begin HALO
C AMD test:   30 IF (NEL .LT. N) THEN
   30 IF (NEL .LT. NREAL) THEN
C end HALO
C=======================================================================
C  GET PIVOT OF MINIMUM DEGREE
C=======================================================================
C       -------------------------------------------------------------
C       find next supervariable for elimination
C       -------------------------------------------------------------
        DO 40 DEG = MINDEG, N
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
C begin HALO
        IF (ME.LE.0) THEN
          write (*,*) ' ERROR 1 in HALO_AMD '
C         return to calling program with error return
          NCMPA = -N
          GOTO 500
        ENDIF
C end HALO
C         -------------------------------------------------------------
C         remove chosen variable from link list
C         -------------------------------------------------------------
          INEXT = NEXT (ME)
          IF (INEXT .NE. 0) LAST (INEXT) = 0
          HEAD (DEG) = INEXT
C       -------------------------------------------------------------
C       me represents the elimination of pivots nel+1 to nel+nv(me).
C       place me itself as the first in this set.  It will be moved
C       to the nel+nv(me) position when the permutation vectors are
C       computed.
C       -------------------------------------------------------------
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
C=======================================================================
C  CONSTRUCT NEW ELEMENT
C=======================================================================
C       -------------------------------------------------------------
C       At this point, me is the pivotal supervariable.  It will be
C       converted into the current element.  Scan list of the
C       pivotal supervariable, me, setting tree pointers and
C       constructing new list of supervariables for the new element,
C       me.  p is a pointer to the current position in the old list.
C       -------------------------------------------------------------
C       flag the variable "me" as being in Lme by negating nv (me)
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
C         ----------------------------------------------------------
C         construct the new element in place
C         ----------------------------------------------------------
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
C             ----------------------------------------------------
C             i is a principal variable not yet placed in Lme.
C             store i in new list
C             ----------------------------------------------------
              DEGME = DEGME + NVI
C             flag i as being in Lme by negating nv (i)
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
C begin HALO
              IF (DEGREE(I).LE.N) THEN
C end HALO
C             ----------------------------------------------------
C             remove variable i from degree list. (only if i \in V0)
C             ----------------------------------------------------
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
C               i is at the head of the degree list
                HEAD (DEGREE (I)) = INEXT
              ENDIF
C begin HALO
              ENDIF
C end HALO
            ENDIF
   60     CONTINUE
C         this element takes no new memory in iw:
          NEWMEM = 0
        ELSE
C         ----------------------------------------------------------
C         construct the new element in empty space, iw (pfree ...)
C         ----------------------------------------------------------
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          KNT1_UPDATED = 0
          DO 120 KNT1 = 1, ELENME + 1
            KNT1_UPDATED = KNT1_UPDATED +1
            IF (KNT1 .GT. ELENME) THEN
C             search the supervariables in me.
              E = ME
              PJ = P
              LN = SLENME
            ELSE
C             search the elements in me.
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
C           -------------------------------------------------------
C           search for different supervariables and add them to the
C           new list, compressing when necessary. this loop is
C           executed once for each element in the list and once for
C           all the supervariables in the list.
C           -------------------------------------------------------
            KNT2_UPDATED = 0
            DO 110 KNT2 = 1, LN
              KNT2_UPDATED = KNT2_UPDATED+1
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
C               -------------------------------------------------
C               compress iw, if necessary
C               -------------------------------------------------
                IF (PFREE .GT. IWLEN) THEN
C                 prepare for compressing iw by adjusting
C                 pointers and lengths so that the lists being
C                 searched in the inner and outer loops contain
C                 only the remaining entries.
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1_UPDATED
C                 Reset KNT1_UPDATED in case of recompress 
C                 at same iteration of the loop 120
                  KNT1_UPDATED = 0
C                 Check if anything left in supervariable ME
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0_8
                  PE (E) = PJ
                  LEN (E) = LN - KNT2_UPDATED
C                 Reset KNT2_UPDATED in case of recompress 
C                 at same iteration of the loop 110
                  KNT2_UPDATED = 0
C                 Check if anything left in element E
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
C                 store first item in pe
C                 set first entry to -item
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = int(IW (PN),8)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
C                 psrc/pdst point to source/destination
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
C                 while loop:
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
C                   search for next negative entry
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = int(PE (J))
                      PE (J) = PDST
                      PDST = PDST + 1
C                     copy from source to destination
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
C                 move the new partially-constructed element
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
C               -------------------------------------------------
C               i is a principal variable not yet placed in Lme
C               store i in new list
C               -------------------------------------------------
                DEGME = DEGME + NVI
C               flag i as being in Lme by negating nv (i)
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
C begin HALO
              IF (DEGREE(I).LE.N) THEN
C end HALO
C               -------------------------------------------------
C               remove variable i from degree link list 
C                            (only if i in V0)
C               -------------------------------------------------
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
C                 i is at the head of the degree list
                  HEAD (DEGREE (I)) = INEXT
                ENDIF
C begin HALO
              ENDIF
C end HALO
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
C             set tree pointer and flag to indicate element e is
C             absorbed into new element me (the parent of e is me)
              PE (E) = int(-ME,8)
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
C         this element takes newmem new memory in iw (possibly zero)
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = max (MAXMEM, MEM)
        ENDIF
C       -------------------------------------------------------------
C       me has now been converted into an element in iw (pme1..pme2)
C       -------------------------------------------------------------
C       degme holds the external degree of new element
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = int(PME2 - PME1 + 1_8)
C       -------------------------------------------------------------
C       make sure that wflg is not too large.  With the current
C       value of wflg, wflg+n must not cause integer overflow
C       -------------------------------------------------------------
        IF (WFLG .GT. MAXINT_N) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
C=======================================================================
C  COMPUTE (w (e) - wflg) = |Le\Lme| FOR ALL ELEMENTS
C=======================================================================
C       -------------------------------------------------------------
C       Scan 1:  compute the external degrees of previous elements
C       with respect to the current element.  That is:
C            (w (e) - wflg) = |Le \ Lme|
C       for each element e that appears in any supervariable in Lme.
C       The notation Le refers to the pattern (list of
C       supervariables) of a previous element e, where e is not yet
C       absorbed, stored in iw (pe (e) + 1 ... pe (e) + iw (pe (e))).
C       The notation Lme refers to the pattern of the current element
C       (stored in iw (pme1..pme2)).   If (w (e) - wflg) becomes
C       zero, then the element e will be absorbed in scan 2.
C       -------------------------------------------------------------
        DO 150 PME = PME1, PME2
          I = IW (PME)
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
C           note that nv (i) has been negated to denote i in Lme:
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + int(ELN - 1,8)
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
C               unabsorbed element e has been seen in this loop
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
C               e is an unabsorbed element
C               this is the first we have seen e in all of Scan 1
                WE = DEGREE (E) + WNVI
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
C=======================================================================
C  DEGREE UPDATE AND ELEMENT ABSORPTION
C=======================================================================
C       -------------------------------------------------------------
C       Scan 2:  for each i in Lme, sum up the degree of Lme (which
C       is degme), plus the sum of the external degrees of each Le
C       for the elements e appearing within i, plus the
C       supervariables in i.  Place i in hash list.
C       -------------------------------------------------------------
        DO 180 PME = PME1, PME2
          I = IW (PME)
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0_8
          DEG = 0
C         ----------------------------------------------------------
C         scan the element list associated with supervariable i
C         ----------------------------------------------------------
          DO 160 P = P1, P2
            E = IW (P)
C           dext = | Le \ Lme |
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
            ELSE IF (DEXT .EQ. 0) THEN
#if defined (NOAGG3)
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + E
#else
C             aggressive absorption: e is not adjacent to me, but
C             the |Le \ Lme| is 0, so absorb it into me
              PE (E) = int(-ME,8)
              W (E) = 0
#endif
            ENDIF
  160     CONTINUE
C         count the number of elements in i (including me):
          ELEN (I) = int(PN - P1 + 1_8)
C         ----------------------------------------------------------
C         scan the supervariables in the list associated with i
C         ----------------------------------------------------------
          P3 = PN
          DO 170 P = P2 + 1, P1 + int(LEN (I) - 1,8)
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
C             j is unabsorbed, and not in Lme.
C             add to degree and add to new list
              DEG = DEG + NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + int(J,kind=8)
            ENDIF
  170     CONTINUE
C begin HALO
          IF (DEGREE(I).EQ.N+1) DEG = N+1
C end HALO
C         ----------------------------------------------------------
C         update the degree and check for mass elimination
C         ----------------------------------------------------------
#if defined (NOAGG3)
          IF (ELEN(I).EQ.1 .AND. P3.EQ.PN) THEN
#else
          IF (DEG .EQ. 0) THEN
#endif
C           -------------------------------------------------------
C           mass elimination
C           -------------------------------------------------------
C           There is nothing left of this node except for an
C           edge to the current pivot element.  elen (i) is 1,
C           and there are no variables adjacent to node i.
C           Absorb i into the current pivot element, me.
            PE (I) = int(-ME,8)
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
C           -------------------------------------------------------
C           update the upper-bound degree of i
C           -------------------------------------------------------
C           the following degree does not yet include the size
C           of the current element, which is added later:
C begin HALO V6
            IF (DEGREE(I).NE.N+1) THEN
C                I does not belong to halo
                 DEG        = min (DEG, NLEFT)
                 DEGREE (I) = min (DEGREE (I), DEG)
            ENDIF
C end HALO V6
C           -------------------------------------------------------
C           add me to the list for i
C           -------------------------------------------------------
C           move first supervariable to end of list
            IW (PN) = IW (P3)
C           move first element to end of element part of list
            IW (P3) = IW (P1)
C           add new element to front of list.
            IW (P1) = ME
C           store the new length of the list in len (i)
            LEN (I) = int(PN - P1 + 1)
C begin HALO
            IF (DEG.LE.N) THEN
C end HALO
C           -------------------------------------------------------
C           place in hash bucket.  Save hash key of i in last (i).
C           -------------------------------------------------------
            HASH = mod (HASH, HMOD) + 1_8
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
C             the degree list is empty, hash head is -j
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
C             degree list is not empty
C             use last (head (hash)) as hash head
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = int(HASH, kind=kind(LAST))
C begin HALO
            ENDIF
C end HALO
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
C       -------------------------------------------------------------
C       Clear the counter array, w (...), by incrementing wflg.
C       -------------------------------------------------------------
        DMAX = max (DMAX, DEGME)
        WFLG = WFLG + DMAX
C       make sure that wflg+n does not cause integer overflow
        IF (WFLG .GT. MAXINT_N) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
C       at this point, w (1..n) .lt. wflg holds
C=======================================================================
C  SUPERVARIABLE DETECTION
C=======================================================================
        DO 250 PME = PME1, PME2
          I = IW (PME)
C begin HALO
C old AMD          IF (NV (I) .LT. 0) THEN
          IF ( (NV (I) .LT. 0) .AND. (DEGREE(I) .LE. N) ) THEN
C end HALO
C           i is a principal variable in Lme
C           -------------------------------------------------------
C           examine all hash buckets with 2 or more variables.  We
C           do this by examing all unique hash keys for super-
C           variables in the pattern Lme of the current element, me
C           -------------------------------------------------------
            HASH = int(LAST (I),kind=8)
C           let i = head of hash bucket, and empty the hash bucket
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
C             degree list is empty
              I = -J
              HEAD (HASH) = 0
            ELSE
C             degree list is not empty, restore last () of head
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
C           while loop:
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
C             ----------------------------------------------------
C             this bucket has one or more variables following i.
C             scan all of them to see if i can absorb any entries
C             that follow i in hash bucket.  Scatter i into w.
C             ----------------------------------------------------
              LN = LEN (I)
              ELN = ELEN (I)
C             do not flag the first element in the list (me)
              DO 210 P = PE (I) + 1, PE (I) + LN - 1
                W (IW (P)) = WFLG
  210         CONTINUE
C             ----------------------------------------------------
C             scan every other entry j following i in bucket
C             ----------------------------------------------------
              JLAST = I
              J = NEXT (I)
C             while loop:
  220         CONTINUE
              IF (J .NE. 0) THEN
C               -------------------------------------------------
C               check if j and i have identical nonzero pattern
C               -------------------------------------------------
C               jump if i and j do not have same size data structure
                IF (LEN (J) .NE. LN) GO TO 240
C               jump if i and j do not have same number adj elts
                IF (ELEN (J) .NE. ELN) GO TO 240
C               do not flag the first element in the list (me)
                DO 230 P = PE (J) + 1, PE (J) + LN - 1
C                 jump if an entry (iw(p)) is in j but not in i
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
C               -------------------------------------------------
C               found it!  j can be absorbed into i
C               -------------------------------------------------
                PE (J) = int(-I,8)
C               both nv (i) and nv (j) are negated since they
C               are in Lme, and the absolute values of each
C               are the number of variables in i and j:
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
C               delete j from hash bucket
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
C               -------------------------------------------------
  240           CONTINUE
C               j cannot be absorbed into i
C               -------------------------------------------------
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
C             ----------------------------------------------------
C             no more variables can be absorbed into i
C             go to next i in bucket and clear flag array
C             ----------------------------------------------------
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
C=======================================================================
C  RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVAR. FROM ELEMENT
C=======================================================================
        P = PME1
        NLEFT = N - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
C           i is a principal variable in Lme
C           restore nv (i) to signify that i is principal
            NV (I) = NVI
C begin HALO
            IF (DEGREE(I).LE.N) THEN
C end HALO
C           -------------------------------------------------------
C           compute the external degree (add size of current elem)
C           -------------------------------------------------------
            DEG = min (DEGREE (I) + DEGME - NVI, NLEFT - NVI)
C           -------------------------------------------------------
C           place the supervariable at the head of the degree list
C           -------------------------------------------------------
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (DEG) = I
C           -------------------------------------------------------
C           save the new degree, and find the minimum degree
C           -------------------------------------------------------
            MINDEG = min (MINDEG, DEG)
            DEGREE (I) = DEG
C begin HALO
              ENDIF
C end HALO
C           -------------------------------------------------------
C           place the supervariable in the element pattern
C           -------------------------------------------------------
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
C=======================================================================
C  FINALIZE THE NEW ELEMENT
C=======================================================================
        NV (ME) = NVPIV + DEGME
C       nv (me) is now the degree of pivot (including diagonal part)
C       save the length of the list for the new element me
        LEN (ME) = int(P - PME1)
        IF (LEN (ME) .EQ. 0) THEN
C         there is nothing left of the current pivot element
          PE (ME) = 0_8
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
C         element was not constructed in place: deallocate part
C         of it (final size is less than or equal to newmem,
C         since newly nonprincipal variables have been removed).
          PFREE = P
          MEM = MEM - NEWMEM + LEN (ME)
        ENDIF
C=======================================================================
C       END WHILE (selecting pivots)
      GO TO 30
      ENDIF
C=======================================================================
C begin HALO V2
      IF (NEL.LT.N) THEN 
C
C     All possible pivots (not flagged have been eliminated).
C     We amalgamate all flagged variables at the root and 
C     we finish the elimination tree.
C          1/ Go through all
C          non absorbed elements (root of the subgraph)
C          and absorb in ME
C          2/ perform mass elimination of all dense rows
           DO DEG = MINDEG, N
             ME = HEAD (DEG)
             IF (ME .GT. 0) GO TO 51
           ENDDO
   51      MINDEG = DEG
C
           IF (ME.NE.LISTVAR_SCHUR(1)) THEN
             write(6,*) ' ERROR 2 in MUMPS_HAMD '
             write(6,*) ' wrong principal var for Schur !!'
             NCMPA = -N - 2
             CALL MUMPS_ABORT()
           ENDIF
C
           NELME    = -(NEL+1)
           DO X=1,N
            IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
C            X is an unabsorbed element
             PE(X) = int(-ME,8)
C            W(X) = 0 could be suppressed ?? check it
            ELSEIF (DEGREE(X).EQ.N+1) THEN
C            X is a dense row, absorb it in ME (mass elimination)
             NEL   = NEL + NV(X)
             PE(X) = int(-ME,8)
             ELEN(X) = 0
C            Correct value of NV is (secondary variable)
             NV(X) = 0
            ENDIF
           ENDDO
C          ME is the root node
           ELEN(ME) = NELME
C          Correct value of NV is (principal variable)
           NV(ME)   = N-NREAL
           PE(ME)   = 0
C end HALO V2
C
C begin HALO
           IF (NEL.NE.N) THEN
         write(*,*) ' ERROR 2 in MUMPS_HAMD NEL, N=', NEL,N
         NCMPA = -N - 1
        ENDIF
      ENDIF
C end HALO
C=======================================================================
C  COMPUTE THE PERMUTATION VECTORS
C=======================================================================
C     ----------------------------------------------------------------
C     The time taken by the following code is O(n).  At this
C     point, elen (e) = -k has been done for all elements e,
C     and elen (i) = 0 has been done for all nonprincipal
C     variables i.  At this point, there are no principal
C     supervariables left, and all elements are absorbed.
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
C     compute the ordering of unordered nonprincipal variables
C     ----------------------------------------------------------------
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
C         ----------------------------------------------------------
C         i is an un-ordered row.  Traverse the tree from i until
C         reaching an element, e.  The element, e, was the
C         principal supervariable of i and all nodes in the path
C         from i to when e was selected as pivot.
C         ----------------------------------------------------------
          J = int(-PE (I))
C         while (j is a variable) do:
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = int(-PE (J))
              GO TO 270
            ENDIF
            E = J
C           ----------------------------------------------------------
C           get the current pivot ordering of e
C           ----------------------------------------------------------
            K = -ELEN (E)
C           ----------------------------------------------------------
C           traverse the path again from i to e, and compress the
C           path (all nodes point to e).  Path compression allows
C           this code to compute in O(n) time.  Order the unordered
C           nodes in the path, and place the element e at the end.
C           ----------------------------------------------------------
            J = I
C           while (j is a variable) do:
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = int(-PE (J))
              PE (J) = int(-E,8)
              IF (ELEN (J) .EQ. 0) THEN
C               j is an unordered row
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
C         leave elen (e) negative, so we know it is an element
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
C     ----------------------------------------------------------------
C     reset the inverse permutation (elen (1..n)) to be positive,
C     and compute the permutation (last (1..n)).
C     ----------------------------------------------------------------
      DO 300 I = 1, N
        K = abs (ELEN (I))
        LAST (K) = I
        ELEN (I) = K
  300 CONTINUE
C=======================================================================
C  RETURN THE MEMORY USAGE IN IW
C=======================================================================
C     If maxmem is less than or equal to iwlen, then no compressions
C     occurred, and iw (maxmem+1 ... iwlen) was unused.  Otherwise
C     compressions did occur, and iwlen would have had to have been
C     greater than or equal to maxmem for no compressions to occur.
C     Return the value of maxmem in the pfree argument.
 500  PFREE = MAXMEM
C===============================
C     Save IPE in PARENT array
      DO I=1,N
       PARENT(I) = int(PE(I))
      ENDDO
C===============================
      RETURN
      END SUBROUTINE MUMPS_HAMD
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C Description of MUMPS_HAMF4: 
C MUMPS_HAMF4 is a modified version of halo AMD routine MUMPS_HAMD
C implementing an approximate minimum fill-in heuritic.
C Version provided to F. Pellegrini on Nov 2000 to be used in SCOTCH.
C Approximation of level4 of the minimum fill heuristic
C
C   Restrictive integer 64 bit variant :
C   it is assumed that IW array size can exceed 32-bit integer
C
      SUBROUTINE MUMPS_HAMF4
     &           (NORIG, N, COMPUTE_PERM, NBBUCK, 
     &                   IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
     &                   LAST, NCMPA, DEGREE, WF, NEXT, W, HEAD
     &                   , PARENT
     &                   )
      IMPLICIT NONE
C
C Parameters
C    Input not modified
C    N : number of nodes in the complete graph including halo
C    NORIG :
C       if compressed graph (nv(1).ne-1) then
C         NORIG is the sum(nv(i)) for i \in [1:N]
C       else NORIG = N 
      INTEGER, INTENT(IN)    :: NORIG, N, NBBUCK
      LOGICAL, INTENT(IN)     :: COMPUTE_PERM
      INTEGER(8), INTENT(IN) :: IWLEN
C     Input undefined on output 
      INTEGER, INTENT(INOUT)  :: LEN(N), IW(IWLEN)
C     NV also meaningful as input to encode compressed graphs
      INTEGER, INTENT(INOUT)  :: NV(N)
C 
C     Output only 
      INTEGER, INTENT(OUT)   :: NCMPA
      INTEGER, INTENT(OUT)   :: ELEN(N), LAST(N)
      INTEGER, INTENT(OUT)   :: PARENT(N)
C 
C     Input/output
      INTEGER(8), INTENT(INOUT) :: PFREE
      INTEGER(8), INTENT(INOUT) :: PE(N)
C 
C     Internal Workspace only
C       Min fill approximation one extra array of size NBBUCK+2 
C       is also needed
      INTEGER     :: NEXT(N), DEGREE(N), W(N)
      INTEGER     :: HEAD(0:NBBUCK+1), WF(N)
C
C  Comments on the OUTPUT:
C  ----------------------
C  Let V= V0 U V1 the nodes of the initial graph (|V|=n). 
C  The assembly tree corresponds to the tree 
C    of the supernodes (or supervariables). Each node of the 
C    assembly tree is then composed of one principal variable 
C    and a list of secondary variables. The list of 
C    variable of a node (principal + secondary variables) then 
C    describes the structure of the diagonal bloc of the 
C    supernode. 
C  The elimination tree denotes the tree of all the variables(=node) and 
C    is therefore of order n.
C
C  The arrays NV(N) and PE(N) give a description of the 
C  assembly tree. 
C  Note that on output 
C   INTEGER(8) PE array is copied on output into
C   INTEGER PARENT array
C  
C   1/ Description of array nv(N) (on OUTPUT)
C    nv(i)=0 i is a secondary variable 
C    nv(i) >0 i is a principal variable, nv(i) holds the 
C                  the number of elements in column i of L (true degree of i)
C    With compressed graph (nv(1).ne.-1 on input), 
C    nv(i) can be greater than N since degree can be as large as NORIG
C
C   2/ Description of array PE(N) (on OUTPUT)
C       Note that on 
C       pe(i) = -(father of variable/node i) in the elimination tree:
C       If nv (i) .gt. 0, then i represents a node in the assembly tree,
C       and the parent of i is -pe (i), or zero if i is a root.
C       If nv (i) = 0, then (i,-pe (i)) represents an edge in a
C       subtree, the root of which is a node in the assembly tree.
C   
C   3/ Example:
C      Let If be a root node father of Is in the assembly tree. 
C      If is the principal 
C      variable of the node If and let If1, If2, If3 be the secondary variables
C      of node If.
C      Is is the principal 
C      variable of the node Is and let Is1, Is2 be the secondary variables
C      of node Is.
C      
C      THEN: 
C        NV(If1)=NV(If2)=NV(If3) = 0  (secondary variables)
C        NV(Is1)=NV(Is2) = 0  (secondary variables)
C        NV(If) > 0  ( principal variable)
C        NV(Is) > 0  ( principal variable)
C        PE(If)  = 0 (root node)
C        PE(Is)  = -If (If is the father of Is in the assembly tree)
C        PE(If1)=PE(If2)=PE(If3)= -If  ( If is the principal variable)
C        PE(Is1)=PE(Is2)= -Is  ( Is is the principal variable)
C      
C
C
C HALOAMD_V1: (September 1997)
C **********
C Initial version designed to experiment the numerical (fill-in) impact 
C of taking into account the halo. This code should be able 
C to experiment no-halo, partial halo, complete halo.
C DATE: September 17th 1997
C
C HALOAMD is designed to process a gragh composed of two types
C            of nodes, V0 and V1, extracted from a larger gragh. 
C            V0^V1 = {}, 
C       
C            We used Min. degree heuristic to order only 
C            nodes in V0, but the adjacency to nodes
C            in V1 is taken into account during ordering.
C            Nodes in V1 are odered at last.
C            Adjacency between nodes of V1 need not be provided,
C            however |len(i)| must always corresponds to the number of 
C            edges effectively provided in the adjacency list of i.
C          On input :
c          ********
C            Nodes INODE in V1 are flagged with len(INODE) = -degree 
C                           if len(i) =0 and i \in V1 then 
C                           len(i) must be set on input to -NORIG-1
C          ERROR return (negative values in ncmpa)
C          ************
C            negative value in ncmpa indicates an error detected 
C               by HALOAMD.
C
C            The graph provided MUST follow the rule:
C             if (i,j) is an edge in the gragh then 
C             j must be in the adjacency list of i AND 
C             i must be in the adjacency list of j.
C    REMARKS
C    -------
C        
C        1/  Providing edges between nodes of V1 should not 
C            affect the final ordering, only the amount of edges 
C            of the halo should effectively affect the solution.
C            This code should work in the following cases:
C              1/ halo not provided
C              2/ halo partially provided
C              3/ complete halo
C              4/ complete halo+interconnection between nodes of V1.
C
C              1/ should run and provide identical results (w.r.t to current 
C               implementation of AMD in SCOTCH).
C             3/ and 4 should provide identical results.
C
C        2/ All modifications of the AMD initial code are indicated
C           with begin HALO .. end HALO
C
C            
C   Given a representation of the nonzero pattern of a symmetric matrix,
C       A, (excluding the diagonal) perform an approximate minimum
C       fill-in heuristic. Aggresive absorption is
C       used to tighten the bound on the degree.  This can result an
C       significant improvement in the quality of the ordering for
C       some matrices.
C-----------------------------------------------------------------------
C INPUT ARGUMENTS (unaltered):
C-----------------------------------------------------------------------
C n:    The matrix order.
C       Restriction:  n .ge. 1
C compute_perm : indicates if permutations should be computed 
C         on output in last/elen 
C iwlen:        The length of iw (1..iwlen).  On input, the matrix is
C       stored in iw (1..pfree-1).  However, iw (1..iwlen) should be
C       slightly larger than what is required to hold the matrix, at
C       least iwlen .ge. pfree + n is recommended.  Otherwise,
C       excessive compressions will take place.
C       *** We do not recommend running this algorithm with ***
C       ***      iwlen .lt. pfree + n.                      ***
C       *** Better performance will be obtained if          ***
C       ***      iwlen .ge. pfree + n                       ***
C       *** or better yet                                   ***
C       ***      iwlen .gt. 1.2 * pfree                     ***
C       *** (where pfree is its value on input).            ***
C       The algorithm will not run at all if iwlen .lt. pfree-1.
C
C       Restriction: iwlen .ge. pfree-1
C-----------------------------------------------------------------------
C INPUT/OUPUT ARGUMENTS:
C-----------------------------------------------------------------------
C pe:   On input, pe (i) is the index in iw of the start of row i, or
C       zero if row i has no off-diagonal non-zeros.
C
C       During execution, it is used for both supervariables and
C       elements:
C
C       * Principal supervariable i:  index into iw of the
C               description of supervariable i.  A supervariable
C               represents one or more rows of the matrix
C               with identical nonzero pattern.
C       * Non-principal supervariable i:  if i has been absorbed
C               into another supervariable j, then pe (i) = -j.
C               That is, j has the same pattern as i.
C               Note that j might later be absorbed into another
C               supervariable j2, in which case pe (i) is still -j,
C               and pe (j) = -j2.
C       * Unabsorbed element e:  the index into iw of the description
C               of element e, if e has not yet been absorbed by a
C               subsequent element.  Element e is created when
C               the supervariable of the same name is selected as
C               the pivot.
C       * Absorbed element e:  if element e is absorbed into element
C               e2, then pe (e) = -e2.  This occurs when the pattern of
C               e (that is, Le) is found to be a subset of the pattern
C               of e2 (that is, Le2).  If element e is "null" (it has
C               no nonzeros outside its pivot block), then pe (e) = 0.
C
C       On output, pe holds the assembly tree/forest, which implicitly
C       represents a pivot order with identical fill-in as the actual
C       order (via a depth-first search of the tree).
C
C       On output:
C       If nv (i) .gt. 0, then i represents a node in the assembly tree,
C       and the parent of i is -pe (i), or zero if i is a root.
C       If nv (i) = 0, then (i,-pe (i)) represents an edge in a
C       subtree, the root of which is a node in the assembly tree.
C       On output:  (PE is copied on output into PARENT array)
C
C pfree:        On input, the matrix is stored in iw (1..pfree-1) and
C       the rest of the array iw is free.
C       During execution, additional data is placed in iw, and pfree
C       is modified so that components  of iw from pfree are free.
C       On output, pfree is set equal to the size of iw that
C       would have been needed for no compressions to occur.  If
C       ncmpa is zero, then pfree (on output) is less than or equal to
C       iwlen, and the space iw (pfree+1 ... iwlen) was not used.
C       Otherwise, pfree (on output) is greater than iwlen, and all the
C       memory in iw was used.
C
C nv:   On input, encoding of compressed graph:
C        if NV(1) = -1 then graph is not compressed otherwise
C        NV(I) holds the weight of node I. 
C       During execution, abs (nv (i)) is equal to the number of rows
C       that are represented by the principal supervariable i.  If i is
C       a nonprincipal variable, then nv (i) = 0.  Initially,
C       nv (i) = 1 for all i.  nv (i) .lt. 0 signifies that i is a
C       principal variable in the pattern Lme of the current pivot
C       element me.  On output, nv (e) holds the true degree of element
C       e at the time it was created (including the diagonal part).
C begin HALO
C       On output, nv(I) can be used to find node in set V1.
C       Not true anymore : ( nv(I) = N+1 characterizes nodes in V1
C                 instead nodes in V1 are considered as a dense root node )
C end HALO
C-----------------------------------------------------------------------
C INPUT/MODIFIED (undefined on output):
C-----------------------------------------------------------------------
C len:  On input, len (i) 
C           positive or null (>=0) : i \in V0 and 
C                     len(i) holds the number of entries in row i of the
C                     matrix, excluding the diagonal.  
C           negative (<0) : i \in V1, and 
C                     -len(i) hold the number of entries in row i of the
C                     matrix, excluding the diagonal.
C       The contents of len (1..n)
C       are undefined on output.
C iw:   On input, iw (1..pfree-1) holds the description of each row i
C       in the matrix.  The matrix must be symmetric, and both upper
C       and lower triangular parts must be present.  The diagonal must
C       not be present.  Row i is held as follows:
C
C               len (i):  the length of the row i data structure
C               iw (pe (i) ... pe (i) + len (i) - 1):
C                       the list of column indices for nonzeros
C                       in row i (simple supervariables), excluding
C                       the diagonal.  All supervariables start with
C                       one row/column each (supervariable i is just
C                       row i).
C               if len (i) is zero on input, then pe (i) is ignored
C               on input.
C
C               Note that the rows need not be in any particular order,
C               and there may be empty space between the rows.
C
C       During execution, the supervariable i experiences fill-in.
C       This is represented by placing in i a list of the elements
C       that cause fill-in in supervariable i:
C
C               len (i):  the length of supervariable i
C               iw (pe (i) ... pe (i) + elen (i) - 1):
C                       the list of elements that contain i.  This list
C                       is kept short by removing absorbed elements.
C               iw (pe (i) + elen (i) ... pe (i) + len (i) - 1):
C                       the list of supervariables in i.  This list
C                       is kept short by removing nonprincipal
C                       variables, and any entry j that is also
C                       contained in at least one of the elements
C                       (j in Le) in the list for i (e in row i).
C
C       When supervariable i is selected as pivot, we create an
C       element e of the same name (e=i):
C
C               len (e):  the length of element e
C               iw (pe (e) ... pe (e) + len (e) - 1):
C                       the list of supervariables in element e.
C
C       An element represents the fill-in that occurs when supervariable
C       i is selected as pivot (which represents the selection of row i
C       and all non-principal variables whose principal variable is i).
C       We use the term Le to denote the set of all supervariables
C       in element e.  Absorbed supervariables and elements are pruned
C       from these lists when computationally convenient.
C
C       CAUTION:  THE INPUT MATRIX IS OVERWRITTEN DURING COMPUTATION.
C       The contents of iw are undefined on output.
C
C-----------------------------------------------------------------------
C OUTPUT (need not be set on input):
C-----------------------------------------------------------------------
C elen: See the description of iw above.  At the start of execution,
C       elen (i) is set to zero.  During execution, elen (i) is the
C       number of elements in the list for supervariable i.  When e
C       becomes an element, elen (e) = -nel is set, where nel is the
C       current step of factorization.  elen (i) = 0 is done when i
C       becomes nonprincipal.
C
C       For variables, elen (i) .ge. 0 holds 
C       until just before the permutation vectors are computed.  
C       For elements, elen (e) .lt. 0 holds.
C
C       On output elen (1..n) holds the inverse permutation (the same
C       as the 'INVP' argument in Sparspak).  That is, if k = elen (i),
C       then row i is the kth pivot row.  Row i of A appears as the
C       (elen(i))-th row in the permuted matrix, PAP^T.
C last: In a degree list, last (i) is the supervariable preceding i,
C       or zero if i is the head of the list.  In a hash bucket,
C       last (i) is the hash key for i.  last (head (hash)) is also
C       used as the head of a hash bucket if head (hash) contains a
C       degree list (see head, below).
C
C       On output, last (1..n) holds the permutation (the same as the
C       'PERM' argument in Sparspak).  That is, if i = last (k), then
C       row i is the kth pivot row.  Row last (k) of A is the k-th row
C       in the permuted matrix, PAP^T.
C ncmpa:        The number of times iw was compressed.  If this is
C       excessive, then the execution took longer than what could have
C       been.  To reduce ncmpa, try increasing iwlen to be 10% or 20%
C       larger than the value of pfree on input (or at least
C       iwlen .ge. pfree + n).  The fastest performance will be
C       obtained when ncmpa is returned as zero.  If iwlen is set to
C       the value returned by pfree on *output*, then no compressions
C       will occur.
C begin HALO
C        on output ncmpa <0 --> error detected during HALO_AMD:
C           error 1: ncmpa = -N , ordering was stopped.
C end HALO
C
C-----------------------------------------------------------------------
C LOCAL (not input or output - used only during execution):
C-----------------------------------------------------------------------
C degree:       If i is a supervariable, then degree (i) holds the
C       current approximation of the external degree of row i (an upper
C       bound).  The external degree is the number of nonzeros in row i,
C       minus abs (nv (i)) (the diagonal part).  The bound is equal to
C       the external degree if elen (i) is less than or equal to two.
C       We also use the term "external degree" for elements e to refer
C       to |Le \ Lme|.  If e is an element, then degree (e) holds |Le|,
C       which is the degree of the off-diagonal part of the element e
C       (not including the diagonal part).
C begin HALO
C       while processing variables degree(I) =   -NBBUCK-1 (=N2) 
C                                  indicates that i belongs to V1
C end HALO
C
C head: head is used for degree lists.  head (deg) is the first
C       supervariable in a degree list (all supervariables i in a
C       degree list deg have the same approximate degree, namely,
C       deg = degree (i)).  If the list deg is empty then
C       head (deg) = 0.
C
C       During supervariable detection head (hash) also serves as a
C       pointer to a hash bucket.
C       If head (hash) .gt. 0, there is a degree list of degree hash.
C               The hash bucket head pointer is last (head (hash)).
C       If head (hash) = 0, then the degree list and hash bucket are
C               both empty.
C       If head (hash) .lt. 0, then the degree list is empty, and
C               -head (hash) is the head of the hash bucket.
C       After supervariable detection is complete, all hash buckets
C       are empty, and the (last (head (hash)) = 0) condition is
C       restored for the non-empty degree lists.
C next: next (i) is the supervariable following i in a link list, or
C       zero if i is the last in the list.  Used for two kinds of
C       lists:  degree lists and hash buckets (a supervariable can be
C       in only one kind of list at a time).
C w:    The flag array w determines the status of elements and
C       variables, and the external degree of elements.
C
C       for elements:
C          if w (e) = 0, then the element e is absorbed
C          if w (e) .ge. wflg, then w (e) - wflg is the size of
C               the set |Le \ Lme|, in terms of nonzeros (the
C               sum of abs (nv (i)) for each principal variable i that
C               is both in the pattern of element e and NOT in the
C               pattern of the current pivot element, me).
C          if wflg .gt. w (e) .gt. 0, then e is not absorbed and has
C               not yet been seen in the scan of the element lists in
C               the computation of |Le\Lme| in loop 150 below.
C
C       for variables:
C          during supervariable detection, if w (j) .ne. wflg then j is
C          not in the pattern of variable i
C
C       The w array is initialized by setting w (i) = 1 for all i,
C       and by setting wflg = 2.  It is reinitialized if wflg becomes
C       too large (to ensure that wflg+n does not cause integer
C       overflow).
C
C wf : integer array  used to store the already filled area of 
C      the variables adajcent to current pivot. 
C      wf is then used to update the score of variable i.
C
C-----------------------------------------------------------------------
C LOCAL INTEGERS:
C-----------------------------------------------------------------------
      INTEGER :: DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     &        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     &        LENJ, LN, ME, MINDEG, NEL,
     &        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X,
     &        NBFLAG, LASTD, NELME, WF3, WF4, N2, PAS
      INTEGER :: NLEFT_V1
       INTEGER KNT1_UPDATED, KNT2_UPDATED
       INTEGER(8) :: MAXMEM, MEM, NEWMEM
       INTEGER    :: MAXINT_N
       INTEGER(8) :: HASH, HMOD
       DOUBLE PRECISION RMF, RMF1 
       DOUBLE PRECISION dummy
       INTEGER idummy
C deg:        the degree of a variable or element
C degme:      size, |Lme|, of the current element, me (= degree (me))
C dext:       external degree, |Le \ Lme|, of some element e
C dmax:       largest |Le| seen so far
C e:          an element
C elenme:     the length, elen (me), of element list of pivotal var.
C eln:        the length, elen (...), of an element list
C hash:       the computed value of the hash function
C hmod:       the hash function is computed modulo hmod = max (1,n-1)
C i:          a supervariable
C ilast:      the entry in a link list preceding i
C inext:      the entry in a link list following i
C j:          a supervariable
C jlast:      the entry in a link list preceding j
C jnext:      the entry in a link list, or path, following j
C k:          the pivot order of an element or variable
C knt1:       loop counter used during element construction
C knt2:       loop counter used during element construction
C knt3:       loop counter used during compression
C lenj:       len (j)
C ln:         length of a supervariable list
C maxint_n:   large integer to test risk of overflow on wflg
C maxmem:     amount of memory needed for no compressions
C me:         current supervariable being eliminated, and the
C                     current element created by eliminating that
C                     supervariable
C mem:        memory in use assuming no compressions have occurred
C mindeg:     current minimum degree
C nel:        number of pivots selected so far
C newmem:     amount of new memory needed for current pivot element
C nleft:      n - nel, the number of nonpivotal rows/columns remaining
C nvi:        the number of variables in a supervariable i (= nv (i))
C nvj:        the number of variables in a supervariable j (= nv (j))
C nvpiv:      number of pivots in current element
C slenme:     number of variables in variable list of pivotal variable
C we:         w (e)
C wflg:       used for flagging the w array.  See description of iw.
C wnvi:       wflg - nv (i)
C x:          either a supervariable or an element
C wf3:  off diagoanl block area
C wf4:  diagonal block area
C mf : Minimum fill
C begin HALO
C nbflag:     number of flagged entries in the initial gragh.
C nreal :     number of entries on which ordering must be perfomed
C             (nreel = N- nbflag)
C nelme number of pivots selected when reaching the root
C lastd index of the last row in the list of dense rows
C end HALO
C-----------------------------------------------------------------------
C LOCAL POINTERS:
C-----------------------------------------------------------------------
      INTEGER(8) :: P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, 
     &              PN, PSRC
C             Any parameter (pe (...) or pfree) or local variable
C             starting with "p" (for Pointer) is an index into iw,
C             and all indices into iw use variables starting with
C             "p."  The only exception to this rule is the iwlen
C             input argument.
C p:          pointer into lots of things
C p1:         pe (i) for some variable i (start of element list)
C p2:         pe (i) + elen (i) -  1 for some var. i (end of el. list)
C p3:         index of first supervariable in clean list
C pdst:       destination pointer, for compression
C pend:       end of memory to compress
C pj:         pointer into an element or variable
C pme:        pointer into the current element (pme1...pme2)
C pme1:       the current element, me, is stored in iw (pme1...pme2)
C pme2:       the end of the current element
C pn:         pointer into a "clean" variable, also used to compress
C psrc:       source pointer, for compression
C-----------------------------------------------------------------------
C  FUNCTIONS CALLED:
C-----------------------------------------------------------------------
      INTRINSIC max, min, mod, huge
      INTEGER TOTEL
      LOGICAL COMPRESS
C=======================================================================
C  INITIALIZATIONS
C=======================================================================
C     HEAD (0:NBBUCK+1)
C
C idummy holds the largest integer - 1
C dummy  = dble (idummy)
      idummy = huge(idummy) - 1
      dummy = dble(idummy)
C     variable with degree equal to N2 are in halo
C     bucket NBBUCK+1 used for HALO variables
      N2 = -NBBUCK-1
C Distance betweeen elements of the N, ..., NBBUCK entries of HEAD
C
      PAS = max((N/8), 1)
      WFLG = 2
      MAXINT_N=huge(WFLG)-N
      NCMPA = 0
      NEL = 0
      HMOD = int(max (1, NBBUCK-1),kind=8)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      MINDEG = 0
      NLEFT_V1 = 0
C
      NBFLAG = 0
      LASTD  = 0
      HEAD(0:NBBUCK+1) = 0
      DO 10 I = 1, N
        LAST(I) = 0
C        NV(I) = 1
        W(I) = 1
        ELEN (I) = 0
   10 CONTINUE
      IF(NV(1) .LT. 0) THEN
         COMPRESS = .FALSE.
      ELSE
         COMPRESS = .TRUE.
      ENDIF
      IF(COMPRESS) THEN
         TOTEL = 0
         DO I=1,N
            IF (LEN(I).LT.0) THEN 
               DEGREE (I) = N2
               NBFLAG     = NBFLAG +1
               NLEFT_V1   = NLEFT_V1 + NV(I)
               IF (LEN(I).EQ.-NORIG-1) THEN
C     variable in V1 with empty adj list 
                  LEN (I)    = 0
C     Because of compress, we force skipping this
C     entry which is anyway empty
                  PE (I)     = 0_8
               ELSE
                  LEN (I)    = - LEN(I)
               ENDIF
C       end HALO V3
            ELSE
               TOTEL = TOTEL + NV(I)
               DEGREE(I) = 0
               DO P= PE(I) , PE(I)+int(LEN(I)-1,8)
                  DEGREE(I) = DEGREE(I) + NV(IW(P))
               ENDDO
C     DEGREE (I) = LEN (I)
            ENDIF
         ENDDO
      ELSE
         DO I=1,N
            NV(I) = 1
            IF (LEN(I).LT.0) THEN 
               DEGREE (I) = N2
               NBFLAG     = NBFLAG +1
               NLEFT_V1   = NLEFT_V1 + NV(I)
               IF (LEN(I).EQ.-N-1) THEN
                  LEN (I)    = 0
C     Because of compress, we force skipping this
C     entry which is anyway empty
                  PE (I)     = 0_8
               ELSE
                  LEN (I)    = - LEN(I)
               ENDIF
C     end HALO V3
            ELSE
               DEGREE (I) = LEN (I)
            ENDIF
         ENDDO
         TOTEL = N - NBFLAG
      ENDIF
C
C
C     ----------------------------------------------------------------
C     initialize degree lists and eliminate rows with no off-diag. nz.
C     ----------------------------------------------------------------
      DO 20 I = 1, N
        DEG = DEGREE (I)
        IF (DEG.EQ.N2) THEN
C            DEG = N2 (flagged variables are stored 
C                  in the degree list of NBBUCK + 1
C                  (safe: because max 
C                         max value of degree is NBBUCK)
C
             DEG = NBBUCK + 1
             IF (LASTD.EQ.0) THEN
C              degree list is empty
               LASTD     = I
               HEAD(DEG) = I
               NEXT(I)   = 0
               LAST(I)   = 0
             ELSE
               NEXT(LASTD) = I
               LAST(I)     = LASTD
               LASTD       = I
               NEXT(I)     = 0
             ENDIF
         GOTO 20
        ENDIF
C
C
        IF (DEG .GT. 0) THEN
          WF(I) = DEG
C   version 1
           IF (DEG.GT.NORIG) THEN
            DEG = min(((DEG-NORIG)/PAS) + NORIG, NBBUCK)
           ENDIF
C           Note that if deg=0 then 
C           No fill-in will occur, 
C           but one variable is adjacent to I
C          ----------------------------------------------------------
C          place i in the degree list corresponding to its degree
C          ----------------------------------------------------------
           INEXT = HEAD (DEG)
           IF (INEXT .NE. 0) LAST (INEXT) = I
           NEXT (I) = INEXT
           HEAD (DEG) = I
        ELSE
C         ----------------------------------------------------------
C         we have a variable that can be eliminated at once because
C         there is no off-diagonal non-zero in its row.
C         ----------------------------------------------------------
          NEL = NEL + NV(I)
          ELEN (I) = -NEL
          PE (I) = 0_8
          W (I) = 0
        ENDIF
C=======================================================================
C
   20 CONTINUE
C=======================================================================
C  WHILE (selecting pivots) DO
C=======================================================================
      NLEFT = TOTEL-NEL + NLEFT_V1
C=======================================================================
C =====================================================================
   30 IF (NEL .LT. TOTEL) THEN
C =====================================================================
C  GET PIVOT OF MINIMUM DEGREE
C=======================================================================
C       -------------------------------------------------------------
C       find next supervariable for elimination
C       -------------------------------------------------------------
        DO 40 DEG = MINDEG, NBBUCK
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        IF (ME.LE.0) THEN
          NCMPA = -N
          CALL MUMPS_ABORT()
        ENDIF
       IF (DEG.GT.NORIG) THEN
C        -------------------------------
C        Linear search to find variable 
C        with best score in the list
C        -------------------------------
C        While end of list list not reached
C         NEXT(J) = 0
         J = NEXT(ME)
         K = WF(ME)
   55    CONTINUE
         IF (J.GT.0) THEN
          IF (WF(J).LT.K) THEN
           ME = J
           K  = WF(ME)
          ENDIF
          J= NEXT(J)
          GOTO 55
         ENDIF
         ILAST = LAST(ME)
         INEXT = NEXT(ME)
         IF (INEXT .NE. 0) LAST (INEXT) = ILAST
         IF (ILAST .NE. 0) THEN
           NEXT (ILAST) = INEXT
         ELSE
C          me is at the head of the degree list
           HEAD (DEG) = INEXT
         ENDIF
C
        ELSE
C         -------------------------------------------------------------
C         remove chosen variable from link list
C         -------------------------------------------------------------
          INEXT = NEXT (ME)
          IF (INEXT .NE. 0) LAST (INEXT) = 0
          HEAD (DEG) = INEXT
        ENDIF
C       -------------------------------------------------------------
C       me represents the elimination of pivots nel+1 to nel+nv(me).
C       place me itself as the first in this set.  It will be moved
C       to the nel+nv(me) position when the permutation vectors are
C       computed.
C       -------------------------------------------------------------
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
C=======================================================================
C  CONSTRUCT NEW ELEMENT
C=======================================================================
C       -------------------------------------------------------------
C       At this point, me is the pivotal supervariable.  It will be
C       converted into the current element.  Scan list of the
C       pivotal supervariable, me, setting tree pointers and
C       constructing new list of supervariables for the new element,
C       me.  p is a pointer to the current position in the old list.
C       -------------------------------------------------------------
C       flag the variable "me" as being in Lme by negating nv (me)
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
C         ----------------------------------------------------------
C         construct the new element in place
C         ----------------------------------------------------------
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + LEN (ME) - 1
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
C             ----------------------------------------------------
C             i is a principal variable not yet placed in Lme.
C             store i in new list
C             ----------------------------------------------------
              DEGME = DEGME + NVI
C             flag i as being in Lme by negating nv (i)
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
              IF (DEGREE(I).NE.N2) THEN
C             ----------------------------------------------------
C             remove variable i from degree list. (only if i \in V0)
C             ----------------------------------------------------
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
C               i is at the head of the degree list
                IF (WF(I).GT.NORIG) THEN
                 DEG = min(((WF(I)-NORIG)/PAS) + NORIG, NBBUCK)
                ELSE
                 DEG = WF(I)
                ENDIF
                HEAD (DEG) = INEXT
              ENDIF
              ENDIF
            ENDIF
   60     CONTINUE
C         this element takes no new memory in iw:
          NEWMEM = 0
        ELSE
C         ----------------------------------------------------------
C         construct the new element in empty space, iw (pfree ...)
C         ----------------------------------------------------------
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          KNT1_UPDATED = 0
          DO 120 KNT1 = 1, ELENME + 1
            KNT1_UPDATED = KNT1_UPDATED +1
            IF (KNT1 .GT. ELENME) THEN
C             search the supervariables in me.
              E = ME
              PJ = P
              LN = SLENME
            ELSE
C             search the elements in me.
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
C           -------------------------------------------------------
C           search for different supervariables and add them to the
C           new list, compressing when necessary. this loop is
C           executed once for each element in the list and once for
C           all the supervariables in the list.
C           -------------------------------------------------------
            KNT2_UPDATED = 0
            DO 110 KNT2 = 1, LN
              KNT2_UPDATED = KNT2_UPDATED+1
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
C               -------------------------------------------------
C               compress iw, if necessary
C               -------------------------------------------------
                IF (PFREE .GT. IWLEN) THEN
C                 prepare for compressing iw by adjusting
C                 pointers and lengths so that the lists being
C                 searched in the inner and outer loops contain
C                 only the remaining entries.
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1_UPDATED
C                 Reset KNT1_UPDATED in case of recompress 
C                 at same iteration of the loop 120
                  KNT1_UPDATED = 0
C                 Check if anything left in supervariable ME
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0_8
                  PE (E) = PJ
                  LEN (E) = LN - KNT2_UPDATED
C                 Reset KNT2_UPDATED in case of recompress 
C                 at same iteration of the loop 110
                  KNT2_UPDATED = 0
C                 Check if anything left in element E
                  IF (LEN (E) .EQ. 0) PE (E) = 0_8
                  NCMPA = NCMPA + 1
C                 store first item in pe
C                 set first entry to -item
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = int(IW (PN),8)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
C                 psrc/pdst point to source/destination
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
C                 while loop:
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
C                   search for next negative entry
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = int(PE (J))
                      PE (J) = PDST
                      PDST = PDST + 1_8
C                     copy from source to destination
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
C                 move the new partially-constructed element
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
C               -------------------------------------------------
C               i is a principal variable not yet placed in Lme
C               store i in new list
C               -------------------------------------------------
                DEGME = DEGME + NVI
C               flag i as being in Lme by negating nv (i)
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
              IF (DEGREE(I).NE.N2) THEN
C               -------------------------------------------------
C               remove variable i from degree link list 
C                            (only if i in V0)
C               -------------------------------------------------
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  IF (WF(I).GT.NORIG) THEN
                   DEG = min(((WF(I)-NORIG)/PAS) + NORIG , NBBUCK)
                  ELSE
                   DEG = WF(I)
                  ENDIF
C                 i is at the head of the degree list
                  HEAD (DEG) = INEXT
                ENDIF
              ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
C             set tree pointer and flag to indicate element e is
C             absorbed into new element me (the parent of e is me)
              PE (E) = int(-ME,8)
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
C         this element takes newmem new memory in iw (possibly zero)
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = max (MAXMEM, MEM)
        ENDIF
C       -------------------------------------------------------------
C       me has now been converted into an element in iw (pme1..pme2)
C       -------------------------------------------------------------
C       degme holds the external degree of new element
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = int(PME2 - PME1 + 1_8)
C       -------------------------------------------------------------
C       make sure that wflg is not too large.  With the current
C       value of wflg, wflg+n must not cause integer overflow
C       -------------------------------------------------------------
        IF (WFLG .GT. MAXINT_N) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
C=======================================================================
C  COMPUTE (w (e) - wflg) = |Le\Lme| FOR ALL ELEMENTS
C=======================================================================
C       -------------------------------------------------------------
C       Scan 1:  compute the external degrees of previous elements
C       with respect to the current element.  That is:
C            (w (e) - wflg) = |Le \ Lme|
C       for each element e that appears in any supervariable in Lme.
C       The notation Le refers to the pattern (list of
C       supervariables) of a previous element e, where e is not yet
C       absorbed, stored in iw (pe (e) + 1 ... pe (e) + iw (pe (e))).
C       The notation Lme refers to the pattern of the current element
C       (stored in iw (pme1..pme2)).   If (w (e) - wflg) becomes
C       zero, then the element e will be absorbed in scan 2.
C       -------------------------------------------------------------
        DO 150 PME = PME1, PME2
          I = IW (PME)
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
C           note that nv (i) has been negated to denote i in Lme:
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + int(ELN - 1,8)
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
C               unabsorbed element e has been seen in this loop
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
C               e is an unabsorbed element
C               this is the first we have seen e in all of Scan 1
                WE = DEGREE (E) + WNVI
                WF(E) = 0
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
C=======================================================================
C  DEGREE UPDATE AND ELEMENT ABSORPTION
C=======================================================================
C       -------------------------------------------------------------
C       Scan 2:  for each i in Lme, sum up the degree of Lme (which
C       is degme), plus the sum of the external degrees of each Le
C       for the elements e appearing within i, plus the
C       supervariables in i.  Place i in hash list.
C       -------------------------------------------------------------
        DO 180 PME = PME1, PME2
          I = IW (PME)
          P1 = PE (I)
          P2 = P1 + ELEN (I) - 1
          PN = P1
          HASH = 0_8
          DEG  = 0
          WF3  = 0
          WF4  = 0
          NVI  = -NV(I)
C         ----------------------------------------------------------
C         scan the element list associated with supervariable i
C         ----------------------------------------------------------
          DO 160 P = P1, P2
            E = IW (P)
C           dext = | Le \ Lme |
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              IF ( WF(E) .EQ. 0 ) THEN
C              First time we meet e : compute wf(e) 
C              which holds the surface associated to element e 
C              it will later be deducted from fill-in 
C              area of all variables adjacent to e
               WF(E) = DEXT * ( (2 * DEGREE(E))  -  DEXT - 1)
              ENDIF
              WF4 = WF4 + WF(E)
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E, kind=8)
            ELSE IF (DEXT .EQ. 0) THEN
#if defined (NOAGG4)
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
#else
C             aggressive absorption: e is not adjacent to me, but
C             the |Le \ Lme| is 0, so absorb it into me
              PE (E) = int(-ME,8)
              W (E) = 0
#endif
            ENDIF
  160     CONTINUE
C         count the number of elements in i (including me):
          ELEN (I) = int(PN - P1 + 1_8)
C         ----------------------------------------------------------
C         scan the supervariables in the list associated with i
C         ----------------------------------------------------------
          P3 = PN
          DO 170 P = P2 + 1_8, P1 + int(LEN (I) - 1,8)
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
C             j is unabsorbed, and not in Lme.
C             add to degree and add to new list
              DEG = DEG + NVJ
              WF3 = WF3 + NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + int(J,kind=8)
            ENDIF
  170     CONTINUE
C
          IF (DEGREE(I).EQ.N2) DEG = N2
C         ----------------------------------------------------------
C         update the degree and check for mass elimination
C         ----------------------------------------------------------
#if defined (NOAGG4)
          IF (ELEN(I).EQ.1 .AND. P3.EQ.PN) THEN
#else
          IF (DEG .EQ. 0) THEN
#endif
C           -------------------------------------------------------
C           mass elimination
C           -------------------------------------------------------
C           There is nothing left of this node except for an
C           edge to the current pivot element.  elen (i) is 1,
C           and there are no variables adjacent to node i.
C           Absorb i into the current pivot element, me.
            PE (I) = int(-ME,8)
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
C           -------------------------------------------------------
C           update the upper-bound degree of i
C           -------------------------------------------------------
C           the following degree does not yet include the size
C           of the current element, which is added later:
            IF (DEGREE(I).NE.N2) THEN
C                I does not belong to halo
                 IF ( DEGREE (I).LT.DEG ) THEN
C                  Our appox degree is loose.
C                  we keep old value. Note that in 
C                  this case we cannot substract WF(I)
C                  for min-fill score.
                   WF4 = 0
                   WF3 = 0
                 ELSE
                   DEGREE(I)  = DEG
                 ENDIF
            ENDIF
C
C           compute WF(I) taking into account size of block 3.0
            WF(I)      = WF4 + 2*NVI*WF3
C           -------------------------------------------------------
C           add me to the list for i
C           -------------------------------------------------------
C           move first supervariable to end of list
            IW (PN) = IW (P3)
C           move first element to end of element part of list
            IW (P3) = IW (P1)
C           add new element to front of list.
            IW (P1) = ME
C           store the new length of the list in len (i)
            LEN (I) = int(PN - P1 + 1)
            IF (DEG.NE.N2) THEN
C           -------------------------------------------------------
C           place in hash bucket.  Save hash key of i in last (i).
C           -------------------------------------------------------
            HASH = mod (HASH, HMOD) + 1_8
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
C             the degree list is empty, hash head is -j
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
C             degree list is not empty
C             use last (head (hash)) as hash head
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = int(HASH,kind=kind(LAST))
            ENDIF
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
C       -------------------------------------------------------------
C       Clear the counter array, w (...), by incrementing wflg.
C       -------------------------------------------------------------
        DMAX = max (DMAX, DEGME)
        WFLG = WFLG + DMAX
C       make sure that wflg+n does not cause integer overflow
        IF (WFLG .GT. MAXINT_N) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
C       at this point, w (1..n) .lt. wflg holds
C=======================================================================
C  SUPERVARIABLE DETECTION
C=======================================================================
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF ( (NV (I) .LT. 0) .AND. (DEGREE(I).NE.N2) ) THEN
C           i is a principal variable in Lme
C           -------------------------------------------------------
C           examine all hash buckets with 2 or more variables.  We
C           do this by examing all unique hash keys for super-
C           variables in the pattern Lme of the current element, me
C           -------------------------------------------------------
            HASH = int(LAST (I),kind=8)
C           let i = head of hash bucket, and empty the hash bucket
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
C             degree list is empty
              I = -J
              HEAD (HASH) = 0
            ELSE
C             degree list is not empty, restore last () of head
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
C           while loop:
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
C             ----------------------------------------------------
C             this bucket has one or more variables following i.
C             scan all of them to see if i can absorb any entries
C             that follow i in hash bucket.  Scatter i into w.
C             ----------------------------------------------------
              LN = LEN (I)
              ELN = ELEN (I)
C             do not flag the first element in the list (me)
              DO 210 P = PE (I) + 1_8, PE (I) + int(LN - 1,8)
                W (IW (P)) = WFLG
  210         CONTINUE
C             ----------------------------------------------------
C             scan every other entry j following i in bucket
C             ----------------------------------------------------
              JLAST = I
              J = NEXT (I)
C             while loop:
  220         CONTINUE
              IF (J .NE. 0) THEN
C               -------------------------------------------------
C               check if j and i have identical nonzero pattern
C               -------------------------------------------------
C               jump if i and j do not have same size data structure
                IF (LEN (J) .NE. LN) GO TO 240
C               jump if i and j do not have same number adj elts
                IF (ELEN (J) .NE. ELN) GO TO 240
C               do not flag the first element in the list (me)
                DO 230 P = PE (J) + 1_8, PE (J) + int(LN - 1,8)
C                 jump if an entry (iw(p)) is in j but not in i
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
C               -------------------------------------------------
C               found it!  j can be absorbed into i
C               -------------------------------------------------
                PE (J) = int(-I,8)
                WF(I)  = max(WF(I),WF(J))
C               both nv (i) and nv (j) are negated since they
C               are in Lme, and the absolute values of each
C               are the number of variables in i and j:
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
C               delete j from hash bucket
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
C               -------------------------------------------------
  240           CONTINUE
C               j cannot be absorbed into i
C               -------------------------------------------------
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
C             ----------------------------------------------------
C             no more variables can be absorbed into i
C             go to next i in bucket and clear flag array
C             ----------------------------------------------------
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
C=======================================================================
C  RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVAR. FROM ELEMENT
C=======================================================================
        P = PME1
        NLEFT = TOTEL - NEL + NLEFT_V1
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
C           i is a principal variable in Lme
C           restore nv (i) to signify that i is principal
            NV (I) = NVI
            IF (DEGREE(I).NE.N2) THEN
C           -------------------------------------------------------
C           compute the external degree (add size of current elem)
C           -------------------------------------------------------
C--------------------------
C--------------------------
            IF (DEGREE (I) + DEGME .GT. NLEFT ) THEN
C
              DEG = DEGREE(I)
              RMF1  = dble(DEG)*dble( (DEG-1) + 2*DEGME )
     &              - dble(WF(I))
              DEGREE(I) = NLEFT - NVI
              DEG       = DEGREE(I) 
              RMF = dble(DEG)*dble(DEG-1) 
     &         -  dble(DEGME-NVI)*dble(DEGME-NVI-1)
              RMF = min(RMF, RMF1)
            ELSE 
              DEG = DEGREE(I)
              DEGREE(I) = DEGREE (I) + DEGME - NVI
C             All previous cliques taken into account (AMF4)
              RMF  = dble(DEG)*dble( (DEG-1) + 2*DEGME ) 
     &              - dble(WF(I))
            ENDIF
C
            RMF =  RMF / dble(NVI+1)
C
            IF (RMF.LT.dummy) THEN
             WF(I) = int ( anint( RMF ))
            ELSEIF (RMF / dble(N) .LT. dummy) THEN 
             WF(I) = int ( anint( RMF/dble(N) ))
            ELSE
             WF(I) = idummy
            ENDIF
            WF(I) = max(1,WF(I))
            DEG = WF(I)
            IF (DEG.GT.NORIG) THEN
              DEG = min(((DEG-NORIG)/PAS) + NORIG, NBBUCK)
            ENDIF
            INEXT = HEAD (DEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (DEG) = I
C           -------------------------------------------------------
C           save the new degree, and find the minimum degree
C           -------------------------------------------------------
            MINDEG = min (MINDEG, DEG)
C begin HALO
              ENDIF
C end HALO
C           -------------------------------------------------------
C           place the supervariable in the element pattern
C           -------------------------------------------------------
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
C=======================================================================
C  FINALIZE THE NEW ELEMENT
C=======================================================================
        NV (ME) = NVPIV + DEGME
C       fill_est = fill_est + nvpiv * (nvpiv + 2 * degme)
C       nv (me) is now the degree of pivot (including diagonal part)
C       save the length of the list for the new element me
        LEN (ME) = int(P - PME1)
        IF (LEN (ME) .EQ. 0) THEN
C         there is nothing left of the current pivot element
          PE (ME) = 0_8
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
C         element was not constructed in place: deallocate part
C         of it (final size is less than or equal to newmem,
C         since newly nonprincipal variables have been removed).
          PFREE = P
          MEM = MEM - NEWMEM + int(LEN (ME),8)
        ENDIF
C=======================================================================
C       END WHILE (selecting pivots)
      GO TO 30
      ENDIF
C=======================================================================
C begin HALO V2
      IF (NEL.LT.NORIG) THEN 
C
C     All possible pivots (not flagged have been eliminated).
C     We amalgamate all flagged variables at the root and 
C     we finish the elimination tree.
C          1/ Go through all
C          non absorbed elements (root of the subgraph)
C          and absorb in ME
C          2/ perform mass elimination of all dense rows
           DO DEG = MINDEG, NBBUCK+1
             ME = HEAD (DEG)
             IF (ME .GT. 0) GO TO 51
           ENDDO
   51      MINDEG = DEG
           NELME    = -(NEL+1)
           DO X=1,N
            IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
C            X is an unabsorbed element
             PE(X) = int(-ME,8)
C            W(X) = 0 could be suppressed ?? check it
            ELSEIF (DEGREE(X).EQ.N2) THEN
C            X is a dense row, absorb it in ME (mass elimination)
             NEL   = NEL + NV(X)
             PE(X) = int(-ME,8)
             ELEN(X) = 0
C            Correct value of NV is (secondary variable)
             NV(X) = 0
            ENDIF
           ENDDO
C          ME is the root node
           ELEN(ME) = NELME
C          Correct value of NV is (principal variable)
           NV(ME)   = NBFLAG
           PE(ME)   = 0_8
        IF (NEL.NE.NORIG) THEN
         NCMPA = -NORIG - 1
         GOTO 500
        ENDIF
      ENDIF
C end HALO
C=======================================================================
C  COMPUTE THE PERMUTATION VECTORS and update TREE
C=======================================================================
C     ----------------------------------------------------------------
C     The time taken by the following code is O(n).  At this
C     point, elen (e) = -k has been done for all elements e,
C     and elen (i) = 0 has been done for all nonprincipal
C     variables i.  At this point, there are no principal
C     supervariables left, and all elements are absorbed.
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
C     compute the ordering of unordered nonprincipal variables
C     ----------------------------------------------------------------
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
C         ----------------------------------------------------------
C         i is an un-ordered row.  Traverse the tree from i until
C         reaching an element, e.  The element, e, was the
C         principal supervariable of i and all nodes in the path
C         from i to when e was selected as pivot.
C         ----------------------------------------------------------
          J = int(-PE (I))
C         while (j is a variable) do:
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = int(-PE (J))
              GO TO 270
            ENDIF
            E = J
C           ----------------------------------------------------------
C           get the current pivot ordering of e
C           ----------------------------------------------------------
            K = -ELEN (E)
C           ----------------------------------------------------------
C           traverse the path again from i to e, and compress the
C           path (all nodes point to e).  Path compression allows
C           this code to compute in O(n) time.  Order the unordered
C           nodes in the path, and place the element e at the end.
C           ----------------------------------------------------------
            J = I
C           while (j is a variable) do:
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = int(-PE (J))
              PE (J) = int(-E,8)
              IF (ELEN (J) .EQ. 0) THEN
C               j is an unordered row
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
C         leave elen (e) negative, so we know it is an element
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      IF (COMPUTE_PERM) THEN
C     ----------------------------------------------------------------
C     reset the inverse permutation (elen (1..n)) to be positive,
C     and compute the pivot order (last (1..n)).
C     ----------------------------------------------------------------
C begin COMPRESS
      IF(COMPRESS) THEN
C       N is the size of the compressed graph.
C       If the graph was compressed on input then
C       indices in ELEN are in [1,TOTEL]
C       We build the inverse of ELEN in LAST (similar to
C       the pivot order but has zeros in it) and then compress
C       it. Since LAST is assumed to be of size N at the
C       interface level, we need another array to store
C       the inverse of ELEN for entries greater than N
C       We use DEGREE.
        LAST(1:N) = 0
        HEAD(1:TOTEL-N)=0
        DO I = 1, N
          K = abs (ELEN (I))
          IF ( K <= N ) THEN
            LAST (K) = I
          ELSE
            HEAD(K-N)=I
          ENDIF
        ENDDO
        I = 1
        DO K = 1, N
          IF(LAST (K) .NE. 0) THEN
            LAST(I) = LAST(K)
            ELEN(LAST(K)) = I
            I = I + 1
          ENDIF
        ENDDO
        DO K = N+1, TOTEL
          IF (HEAD(K-N) .NE. 0) THEN
            LAST(I)=HEAD(K-N)
            ELEN(HEAD(K-N)) = I
            I = I + 1
          ENDIF
        END DO
      ELSE
        DO 300 I = 1, N
           K = abs (ELEN (I))
           LAST (K) = I
           ELEN (I) = K
300     CONTINUE
      ENDIF
C end COMPRESS
      ENDIF
C=======================================================================
C  RETURN THE MEMORY USAGE IN IW
C=======================================================================
C     If maxmem is less than or equal to iwlen, then no compressions
C     occurred, and iw (maxmem+1 ... iwlen) was unused.  Otherwise
C     compressions did occur, and iwlen would have had to have been
C     greater than or equal to maxmem for no compressions to occur.
C     Return the value of maxmem in the pfree argument.
 500  PFREE = MAXMEM
C===============================
C     Save IPE in PARENT array
      DO I=1,N
       PARENT(I) = int(PE(I))
      ENDDO
C===============================
      RETURN
      END SUBROUTINE MUMPS_HAMF4
C
C-----------------------------------------------------------------------
C MUMPS_QAMD: modified version of reference AMD routine MUMPS_ANA_H 
C designed to automatically detect and exploit dense or quasi dense
C rows in the reduced matrix at any step of the minimum degree.
C
C References:
C    P.R. AMESTOY, Recent progress in parallel multifrontal solvers
C      for unsymmetric sparse matrices,
C      Proceedings of the 15th World Congress on Scientific Computation,
C      Modelling and Applied Mathematics, IMACS, Berlin (1997).
C    P.R. AMESTOY (1999), Methodes directes paralleles de
C      resolution des systemes creux de grande taille.
C      Rapport de these d'habilitation de l'INPT.
C
C Date 1997
C ---------
C
      SUBROUTINE MUMPS_QAMD
     &                (TOTEL, COMPUTE_PERM, IVersion, THRESH, NDENSE, 
     &                 N, IWLEN, PE, PFREE, LEN, IW, NV, 
     &                 ELEN, LAST, NCMPA, DEGREE, HEAD, NEXT, W,
     &                 PARENT) 
C    Input not modified
      INTEGER, INTENT(IN)    :: TOTEL, N
      LOGICAL, INTENT(IN)    :: COMPUTE_PERM
      INTEGER, INTENT(IN)    :: IVersion, THRESH
      INTEGER(8), INTENT(IN) :: IWLEN
      INTEGER, INTENT(INOUT)  :: LEN(N), IW(IWLEN)
      INTEGER, INTENT(OUT)   :: NCMPA
      INTEGER, INTENT(OUT)   :: ELEN(N), PARENT(N)
      INTEGER, INTENT(OUT)   :: LAST(N)
      INTEGER(8), INTENT(INOUT) :: PFREE
      INTEGER(8), INTENT(INOUT) :: PE(N)
C     NV also meaningful as input to encode compressed graphs
      INTEGER, INTENT(INOUT)  :: NV(N)
      INTEGER, INTENT(OUT) :: NEXT(N), DEGREE(N), HEAD(TOTEL), W(N)
      INTEGER, INTENT(OUT) :: NDENSE(N)
C The input integer parameter THRESH defines the quasi density:
C THRESH : input parameter (not modified) 
C  THRESH is used to compute THRESM 
C   <=0 or N Only exactly dense rows in the reduced matrix are selected.
C   >1 and <=N THRESH correspond to the munimum density requirement.
C
C       IVersion =
C                   1 : No dense row detection during elimination
C                       Suppressing dense row selection after 1st
C                       and final restrart (Using initial degree of 
C                         quasi dense
C                         rows when restarting and suppress
C                         dense row selection)
C                   else  : All functionalities enabled
C Additionnal parameters/variables due to dense row manipulation:
C PARAMETERS:
C ----------
C          
C Local variables:
C ---------------
      INTEGER THRESM, MINDEN, MAXDEN, NDME
      INTEGER NBD,NBED, NBDM, LASTD, NELME
C      INTEGER DEG1
      LOGICAL IDENSE
      DOUBLE PRECISION RELDEN
C 
C THRESM : Local Integer holding a 
C          potentially modified value of THRESH.
C          When quasi dense rows are reintegrated in the 
C          graph to be processed then THRESM is modified.
C   Note that if one sets THRESM to negative value then
C       <0 Classical AMD algorithm (no dense row detection)
C RELDEN : holds average density to set THRESM automatically
C MINDEN: min degree of quasi-dense rows when restarting
C MAXDEN: max degree of quasi-dense rows when restarting
C NDME  : number of dense row adjacent to me
C NELME number of pivots selected when reching the root
C LASTD index of the last row in the list of dense rows
C NBD is the total number of dense rows selected 
C NBED is the total number of exactly dense rows detected. 
C NBDM is the maximum number of dense rows selected 
C IDENSE is used to indicate that the supervariable I is a dense or
C        quasi-dense row.
C-----------------------------------------------------------------------
C   Given a representation of the nonzero pattern of a symmetric matrix,
C   A, (excluding the diagonal) perform an approximate minimum
C   degree ordering to compute a pivot order
C   such that fill-in in the Cholesky factors A = LL^T is kept low. 
C   Aggressive absorption might be used to
C   tighten the bound on the degree.  This can result a
C   significant improvement in the quality of the ordering for
C   some matrices.
C-----------------------------------------------------------------------
C INPUT ARGUMENTS (unaltered):
C-----------------------------------------------------------------------
C n     : The matrix order.
C         number of supervariables if compress/blocked format
C         Restriction:  n .ge. 1
C totel : Number of variables to eliminate
C         In case of blocked format:
C         each variable i is a supervariable of size nv(i)
C         totel is computed as the sum(nv(i)) for i \in [1:n]
C         the algorithm stops when totel variables are
C         eliminated.
C compute_perm : indicates if permutations should be computed 
C         on output in last/elen 
C iwlen:        The length of iw (1..iwlen).  On input, the matrix is
C       stored in iw (1..pfree-1).  However, iw (1..iwlen) should be
C       slightly larger than what is required to hold the matrix, at
C       least iwlen .ge. pfree + n is recommended.  Otherwise,
C       excessive compressions will take place.
C       *** We do not recommend running this algorithm with ***
C       ***      iwlen .lt. pfree + n.                      ***
C       *** Better performance will be obtained if          ***
C       ***      iwlen .ge. pfree + n                       ***
C       *** or better yet                                   ***
C       ***      iwlen .gt. 1.2 * pfree                     ***
C       *** (where pfree is its value on input).            ***
C       The algorithm will not run at all if iwlen .lt. pfree-1.
C
C       Restriction: iwlen .ge. pfree-1
C-----------------------------------------------------------------------
C INPUT/OUPUT ARGUMENTS:
C-----------------------------------------------------------------------
C pe:   On input, pe (i) is the index in iw of the start of row i, or
C       zero if row i has no off-diagonal non-zeros.
C
C       During execution, it is used for both supervariables and
C       elements:
C
C       * Principal supervariable i:  index into iw of the
C               description of supervariable i.  A supervariable
C               represents one or more rows of the matrix
C               with identical nonzero pattern.
C       * Non-principal supervariable i:  if i has been absorbed
C               into another supervariable j, then pe (i) = -j.
C               That is, j has the same pattern as i.
C               Note that j might later be absorbed into another
C               supervariable j2, in which case pe (i) is still -j,
C               and pe (j) = -j2.
C       * Unabsorbed element e:  the index into iw of the description
C               of element e, if e has not yet been absorbed by a
C               subsequent element.  Element e is created when
C               the supervariable of the same name is selected as
C               the pivot.
C       * Absorbed element e:  if element e is absorbed into element
C               e2, then pe (e) = -e2.  This occurs when the pattern of
C               e (that is, Le) is found to be a subset of the pattern
C               of e2 (that is, Le2).  If element e is "null" (it has
C               no nonzeros outside its pivot block), then pe (e) = 0.
C
C       On output, pe holds the assembly tree/forest, which implicitly
C       represents a pivot order with identical fill-in as the actual
C       order (via a depth-first search of the tree).
C
C       On output:
C       If nv (i) .gt. 0, then i represents a node in the assembly tree,
C       and the parent of i is -pe (i), or zero if i is a root.
C       If nv (i) = 0, then (i,-pe (i)) represents an edge in a
C       subtree, the root of which is a node in the assembly tree.
C
C       On output:  (PE is copied on output into PARENT array)
C
C pfree:        On input, the matrix is stored in iw (1..pfree-1) and
C       the rest of the array iw is free.
C       During execution, additional data is placed in iw, and pfree
C       is modified so that components  of iw from pfree are free.
C       On output, pfree is set equal to the size of iw that
C       would have been needed for no compressions to occur.  If
C       ncmpa is zero, then pfree (on output) is less than or equal to
C       iwlen, and the space iw (pfree+1 ... iwlen) was not used.
C       Otherwise, pfree (on output) is greater than iwlen, and all the
C       memory in iw was used.
C
C nv:   On input, encoding of compressed graph:
C       if nv(1) = -1 then graph is not compressed otherwise
C       nv(I) holds the weight of node I. 
C       During execution, abs (nv (i)) is equal to the number of rows
C       that are represented by the principal supervariable i.  If i is
C       a nonprincipal variable, then nv (i) = 0.  
C       nv (i) .lt. 0 signifies that i is a
C       principal variable in the pattern Lme of the current pivot
C       element me.  
C       On output, nv (e) holds the true degree of element
C       e at the time it was created (including the diagonal part).
C begin HALO
C       On output, nv(I) can be used to find node in set V1.
C       Not true anymore : ( nv(I) = N+1 characterizes nodes in V1. 
C                 instead nodes in V1 are considered as a dense root node )
C end HALO
C-----------------------------------------------------------------------
C INPUT/MODIFIED (undefined on output):
C-----------------------------------------------------------------------
C len:  On input, len (i) holds the number of entries in row i of the
C       matrix, excluding the diagonal.  The contents of len (1..n)
C       are undefined on output.
C iw:   On input, iw (1..pfree-1) holds the description of each row i
C       in the matrix.  The matrix must be symmetric, and both upper
C       and lower triangular parts must be present.  The diagonal must
C       not be present.  Row i is held as follows:
C
C               len (i):  the length of the row i data structure
C               iw (pe (i) ... pe (i) + len (i) - 1):
C                       the list of column indices for nonzeros
C                       in row i (simple supervariables), excluding
C                       the diagonal.  All supervariables start with
C                       one row/column each (supervariable i is just
C                       row i).
C               if len (i) is zero on input, then pe (i) is ignored
C               on input.
C
C               Note that the rows need not be in any particular order,
C               and there may be empty space between the rows.
C
C       During execution, the supervariable i experiences fill-in.
C       This is represented by placing in i a list of the elements
C       that cause fill-in in supervariable i:
C
C               len (i):  the length of supervariable i
C               iw (pe (i) ... pe (i) + elen (i) - 1):
C                       the list of elements that contain i.  This list
C                       is kept short by removing absorbed elements.
C               iw (pe (i) + elen (i) ... pe (i) + len (i) - 1):
C                       the list of supervariables in i.  This list
C                       is kept short by removing nonprincipal
C                       variables, and any entry j that is also
C                       contained in at least one of the elements
C                       (j in Le) in the list for i (e in row i).
C
C       When supervariable i is selected as pivot, we create an
C       element e of the same name (e=i):
C
C               len (e):  the length of element e
C               iw (pe (e) ... pe (e) + len (e) - 1):
C                       the list of supervariables in element e.
C
C       An element represents the fill-in that occurs when supervariable
C       i is selected as pivot (which represents the selection of row i
C       and all non-principal variables whose principal variable is i).
C       We use the term Le to denote the set of all supervariables
C       in element e.  Absorbed supervariables and elements are pruned
C       from these lists when computationally convenient.
C
C       CAUTION:  THE INPUT MATRIX IS OVERWRITTEN DURING COMPUTATION.
C       The contents of iw are undefined on output.
C-----------------------------------------------------------------------
C OUTPUT (need not be set on input):
C-----------------------------------------------------------------------
C elen: See the description of iw above.  At the start of execution,
C       elen (i) is set to zero.  During execution, elen (i) is the
C       number of elements in the list for supervariable i.  When e
C       becomes an element, elen (e) = -nel is set, where nel is the
C       current step of factorization.  elen (i) = 0 is done when i
C       becomes nonprincipal.
C
C       For variables, elen (i) .ge. 0 holds until just before the
C       permutation vectors are computed.  For elements,
C       elen (e) .lt. 0 holds.
C
C       On output elen (1..n) holds the inverse permutation (the same
C       as the 'INVP' argument in Sparspak).  That is, if k = elen (i),
C       then row i is the kth pivot row.  Row i of A appears as the
C       (elen(i))-th row in the permuted matrix, PAP^T.
C last: In a degree list, last (i) is the supervariable preceding i,
C       or zero if i is the head of the list.  In a hash bucket,
C       last (i) is the hash key for i.  last (head (hash)) is also
C       used as the head of a hash bucket if head (hash) contains a
C       degree list (see head, below).
C
C       On output, last (1..n) holds the permutation (the same as the
C       'PERM' argument in Sparspak).  That is, if i = last (k), then
C       row i is the kth pivot row.  Row last (k) of A is the k-th row
C       in the permuted matrix, PAP^T.
C ncmpa:        The number of times iw was compressed.  If this is
C       excessive, then the execution took longer than what could have
C       been.  To reduce ncmpa, try increasing iwlen to be 10% or 20%
C       larger than the value of pfree on input (or at least
C       iwlen .ge. pfree + n).  The fastest performance will be
C       obtained when ncmpa is returned as zero.  If iwlen is set to
C       the value returned by pfree on *output*, then no compressions
C       will occur.
C-----------------------------------------------------------------------
C LOCAL (not input or output - used only during execution):
C-----------------------------------------------------------------------
C degree:       If i is a supervariable, then degree (i) holds the
C       current approximation of the external degree of row i (an upper
C       bound).  The external degree is the number of nonzeros in row i,
C       minus abs (nv (i)) (the diagonal part).  The bound is equal to
C       the external degree if elen (i) is less than or equal to two.
C
C       We also use the term "external degree" for elements e to refer
C       to |Le \ Lme|.  If e is an element, then degree (e) holds |Le|,
C       which is the degree of the off-diagonal part of the element e
C       (not including the diagonal part).
Cdense
C degree (I) =N+1 if I is an exactly dense row in reduced matrix.
C            =N+1+LAST_approximate_external_deg of I   
C                      if I is a quasi dense row in reduced matrix.
C All dense or quasi dense rows are stored in the list pointed 
C       by head(n). Quasi-dense rows (degree(I)=n) are stored first, 
C       and are followed by exactly dense rows in the reduced matrix.
C       LASTD holds the last row in this list of dense rows or is zero
C       if the list is empty.
Cdense
C head: head is used for degree lists.  head (deg) is the first
C       supervariable in a degree list (all supervariables i in a
C       degree list deg have the same approximate degree, namely,
C       deg = degree (i)).  If the list deg is empty then
C       head (deg) = 0.
C
C       During supervariable detection head (hash) also serves as a
C       pointer to a hash bucket.
C       If head (hash) .gt. 0, there is a degree list of degree hash.
C               The hash bucket head pointer is last (head (hash)).
C       If head (hash) = 0, then the degree list and hash bucket are
C               both empty.
C       If head (hash) .lt. 0, then the degree list is empty, and
C               -head (hash) is the head of the hash bucket.
C       After supervariable detection is complete, all hash buckets
C       are empty, and the (last (head (hash)) = 0) condition is
C       restored for the non-empty degree lists.
C next: next (i) is the supervariable following i in a link list, or
C       zero if i is the last in the list.  Used for two kinds of
C       lists:  degree lists and hash buckets (a supervariable can be
C       in only one kind of list at a time).
C w:    The flag array w determines the status of elements and
C       variables, and the external degree of elements.
C
C       for elements:
C          if w (e) = 0, then the element e is absorbed
C          if w (e) .ge. wflg, then w (e) - wflg is the size of
C               the set |Le \ Lme|, in terms of nonzeros (the
C               sum of abs (nv (i)) for each principal variable i that
C               is both in the pattern of element e and NOT in the
C               pattern of the current pivot element, me).
C          if wflg .gt. w (e) .gt. 0, then e is not absorbed and has
C               not yet been seen in the scan of the element lists in
C               the computation of |Le\Lme| in loop 150 below.
C
C       for variables:
C          during supervariable detection, if w (j) .ne. wflg then j is
C          not in the pattern of variable i
C
C       The w array is initialized by setting w (i) = 1 for all i,
C       and by setting wflg = 2.  It is reinitialized if wflg becomes
C       too large (to ensure that wflg+n does not cause integer
C       overflow).
C-----------------------------------------------------------------------
C LOCAL INTEGERS:
C-----------------------------------------------------------------------
      INTEGER :: DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     &        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     &        LENJ, LN, ME, MINDEG, NEL, 
     &        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X
      INTEGER KNT1_UPDATED, KNT2_UPDATED
      INTEGER(8) MAXMEM, MEM, NEWMEM
      INTEGER   :: MAXINT_N
      INTEGER(8):: HASH, HMOD
C deg:        the degree of a variable or element
C degme:      size, |Lme|, of the current element, me (= degree (me))
C dext:       external degree, |Le \ Lme|, of some element e
C dmax:       largest |Le| seen so far
C e:          an element
C elenme:     the length, elen (me), of element list of pivotal var.
C eln:        the length, elen (...), of an element list
C hash:       the computed value of the hash function
C hmod:       the hash function is computed modulo hmod = max (1,n-1)
C i:          a supervariable
C ilast:      the entry in a link list preceding i
C inext:      the entry in a link list following i
C j:          a supervariable
C jlast:      the entry in a link list preceding j
C jnext:      the entry in a link list, or path, following j
C k:          the pivot order of an element or variable
C knt1:       loop counter used during element construction
C knt2:       loop counter used during element construction
C knt3:       loop counter used during compression
C lenj:       len (j)
C ln:         length of a supervariable list
C maxint_n:   large integer to test risk of overflow on wflg
C maxmem:     amount of memory needed for no compressions
C me:         current supervariable being eliminated, and the
C                     current element created by eliminating that
C                     supervariable
C mem:        memory in use assuming no compressions have occurred
C mindeg:     current minimum degree
C nel:        number of pivots selected so far
C newmem:     amount of new memory needed for current pivot element
C nleft:      n - nel, the number of nonpivotal rows/columns remaining
C nvi:        the number of variables in a supervariable i (= nv (i))
C nvj:        the number of variables in a supervariable j (= nv (j))
C nvpiv:      number of pivots in current element
C slenme:     number of variables in variable list of pivotal variable
C we:         w (e)
C wflg:       used for flagging the w array.  See description of iw.
C wnvi:       wflg - nv (i)
C x:          either a supervariable or an element
C-----------------------------------------------------------------------
C LOCAL POINTERS:
C-----------------------------------------------------------------------
      INTEGER(8) P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, 
     &           PN, PSRC, PLN, PELN
C             Any parameter (pe (...) or pfree) or local variable
C             starting with "p" (for Pointer) is an index into iw,
C             and all indices into iw use variables starting with
C             "p."  The only exception to this rule is the iwlen
C             input argument.
C p:          pointer into lots of things
C p1:         pe (i) for some variable i (start of element list)
C p2:         pe (i) + elen (i) -  1 for some var. i (end of el. list)
C p3:         index of first supervariable in clean list
C pdst:       destination pointer, for compression
C pend:       end of memory to compress
C pj:         pointer into an element or variable
C pme:        pointer into the current element (pme1...pme2)
C pme1:       the current element, me, is stored in iw (pme1...pme2)
C pme2:       the end of the current element
C pn:         pointer into a "clean" variable, also used to compress
C psrc:       source pointer, for compression
      LOGICAL COMPRESS
C-----------------------------------------------------------------------
C  FUNCTIONS CALLED:
C-----------------------------------------------------------------------
      INTRINSIC max, min, mod
C=======================================================================
C  INITIALIZATIONS
C=======================================================================
C     ------------------------------------------------------
C     Experiments with automatic setting of parameter THRESH.
C     ------------------------------------------------------
      IF (THRESH.GT.0) THEN 
         THRESM  = min(N,THRESH)
         DO I=1,N
             THRESM = max(THRESM, LEN(I))
          ENDDO
           RELDEN = dble(PFREE-1)/dble(N)
C      RELDEN holds the average density, THRESM the maximum density
         THRESM =  int(RELDEN)*10 + (THRESM-int(RELDEN))/10 + 1
C     ------------------------------------------------------
C     end automatic setting of THRESM
C     ------------------------------------------------------
      ELSE
C        only exactly dense row will be selected
         THRESM = TOTEL
      ENDIF
      IF (THRESM.GE.0) THEN
       IF ((THRESM.GT.TOTEL).OR.(THRESM.LT.2)) THEN 
C      exactly dense rows only
          THRESM = TOTEL
       ENDIF
      ENDIF
      LASTD = 0
      NBD   = 0
      NBED  = 0
      NBDM  = 0
      WFLG = 2
      MAXINT_N=huge(WFLG)-N
      MINDEG = 1
      NCMPA = 0
      NEL = 0
      HMOD = int(max (1, N-1),kind=8)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      DO I = 1, N
        NDENSE(I)= 0
        W (I) = 1
        ELEN (I) = 0
        LAST(I) = 0
      ENDDO
      DO I = 1, TOTEL
        HEAD(I) = 0
      ENDDO
      IF(NV(1) .LT. 0) THEN
         COMPRESS = .FALSE.
      ELSE
         COMPRESS = .TRUE.
      ENDIF
      IF (COMPRESS) THEN
         DO I=1,N
            DEGREE(I) = 0
            DO P= PE(I) , PE(I)+int(LEN(I)-1,8)
               DEGREE(I) = DEGREE(I) + NV(IW(P))
            ENDDO
         ENDDO
      ELSE
         DO I=1,N
            NV(I) = 1
            DEGREE (I) = LEN (I)
         ENDDO
      ENDIF
C     ----------------------------------------------------------------
C     initialize degree lists and eliminate rows with no off-diag. nz.
C     ----------------------------------------------------------------
C         NEXT = 0
      DO 20 I = 1, N
         DEG = DEGREE (I)
         IF (DEG .GT. 0) THEN
C         ----------------------------------------------------------
C         place i in the degree list corresponding to its degree
C         or in the dense row list if i is dense or quasi dense.
C         ----------------------------------------------------------
C         test for row density
            IF ( (THRESM.GE.0) .AND.
     &           (DEG+NV(I).GE.THRESM) ) THEN
C           I will be inserted in the degree list of N
               NBD = NBD+NV(I)
               IF (DEG+NV(I).NE.TOTEL-NEL) THEN
                  DEGREE(I) = DEGREE(I)+TOTEL+1
C            insert I at the beginning of degree list of n
                  DEG = TOTEL
                  INEXT = HEAD (DEG)
                  IF (INEXT .NE. 0) LAST (INEXT) = I
                  NEXT (I) = INEXT
                  HEAD (DEG) = I 
                  LAST(I)  = 0
                  IF (LASTD.EQ.0) LASTD=I
               ELSE
                  NBED = NBED+NV(I)
                  DEGREE(I) = TOTEL+1
C            insert I at the end of degree list of n
                  DEG = TOTEL
                  IF (LASTD.EQ.0) THEN
C              degree list is empty
                     LASTD     = I 
                     HEAD(DEG) = I
                     NEXT(I)   = 0 
                     LAST(I)   = 0
                  ELSE
                     NEXT(LASTD) = I
                     LAST(I)     = LASTD
                     LASTD       = I
                     NEXT(I)     = 0
                  ENDIF
               ENDIF
            ELSE
C           place i in the degree list corresponding to its degree
               INEXT = HEAD (DEG)
               IF (INEXT .NE. 0) LAST (INEXT) = I
               NEXT (I) = INEXT
               HEAD (DEG) = I
            ENDIF
         ELSE
C         ----------------------------------------------------------
C         we have a variable that can be eliminated at once because
C         there is no off-diagonal non-zero in its row.
C         ----------------------------------------------------------
            NEL = NEL + NV(I)
C          NEL = NEL + 1
            ELEN (I) = -NEL
            PE (I) = 0_8
            W (I) = 0
         ENDIF
 20   CONTINUE
C         We suppress dense row selection if none of them was found in A 
C         in the 1st pass
          IF (NBD.EQ.0) THRESM = TOTEL
C
C=======================================================================
C  WHILE (selecting pivots) DO
C=======================================================================
 30       IF (NEL .LT. TOTEL) THEN
C=======================================================================
C  GET PIVOT OF MINIMUM DEGREE
C=======================================================================
C       -------------------------------------------------------------
C       find next supervariable for elimination
C       -------------------------------------------------------------
        DO 40 DEG = MINDEG, TOTEL
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
        IF (DEG.LT.TOTEL)  THEN
C       -------------------------------------------------------------
C       remove chosen variable from link list
C       -------------------------------------------------------------
          INEXT = NEXT (ME)
          IF (INEXT .NE. 0) LAST (INEXT) = 0
          HEAD (DEG) = INEXT
        ELSE
          NBDM = max(NBDM,NBD)
          IF (DEGREE(ME).GT.TOTEL+1) THEN
            MINDEN = NBD
            MAXDEN = 0
            IF (WFLG .GT. MAXINT_N) THEN
             DO  52 X = 1, N
              IF (W (X) .NE. 0) W (X) = 1
  52         CONTINUE
             WFLG = 2
            ENDIF
            WFLG = WFLG + 1
  51        CONTINUE
C           ---------------------------------------------------------
C           remove chosen variable from link list
C           ---------------------------------------------------------
            INEXT = NEXT (ME)
            IF (INEXT .NE. 0) THEN 
               LAST (INEXT) = 0
            ELSE
               LASTD = 0
            ENDIF
C           ----------------------------------------------------------
c           build adjacency list of ME in quotient gragh
C           and calculate its external degree  in ndense(me)
C           ----------------------------------------------------------
            NDENSE(ME) = 0
            W(ME)      = WFLG
            P1 = PE(ME)
            P2 = P1 + int(LEN(ME) -1,8)
C           PLN-1 holds the pointer in IW to the last elet/var in adj list
C              of ME.  LEN(ME) will then be set to PLN-P1
C           PELN-1 hold the pointer in IW to the last elet in adj list
C              of ME.  ELEN(ME) will then be set to PELN-P1
C           element adjacent to ME
            PLN       = P1
            PELN      = P1
            DO 55 P=P1,P2
              E= IW(P)
              IF (W(E).EQ.WFLG) GOTO 55
              W(E) = WFLG
              IF (PE(E).LT.0_8) THEN
C              E is a nonprincipal variable or absorbed element
                X = E
  53            X = int(-PE(X))
                IF (W(X) .EQ.WFLG) GOTO 55
                W(X) = WFLG
                IF ( PE(X) .LT. 0_8 ) GOTO 53
                E = X
              ENDIF
C             -------------------------------------------
C             E is an unabsorbed element or a "dense" row
C                 (NOT already flagged)
C             -------------------------------------------
              IF (ELEN(E).LT.0) then
C              E is a new element in adj(ME)
               NDENSE(E) = NDENSE(E) - NV(ME)
               IW(PLN) = IW(PELN)
               IW(PELN) = E
               PLN  = PLN+1_8
               PELN = PELN + 1_8
C              update ndense of ME with all unflagged dense
C              rows in E
               PME1 = PE(E)
               DO 54 PME = PME1, PME1+int(LEN(E)-1,8)
                X = IW(PME)
                IF ((ELEN(X).GE.0).AND.(W(X).NE.WFLG)) THEN
C                X is a dense row 
                 NDENSE(ME) = NDENSE(ME) + NV(X)
                 W(X) = WFLG
                ENDIF
 54            CONTINUE
              ELSE
C              E is a dense row 
               NDENSE(ME) = NDENSE(ME) + NV(E)
               IW(PLN)=E
               PLN = PLN+1_8
              ENDIF
  55        CONTINUE
C           ----------------------------------------------
C           DEGREE(ME)-(N+1) holds last external degree computed
C           when Me was detected as dense
C           NDENSE(ME) is the exact external degree of ME
C           ----------------------------------------------
            WFLG     = WFLG + 1
            LEN(ME)  = int(PLN-P1)
            ELEN(ME) = int(PELN-P1)
            NDME = NDENSE(ME)+NV(ME)
            MINDEN = min (MINDEN, NDME)
            MAXDEN = max (MAXDEN, NDME)
C            If we want to select ME as exactly dense (NDME.EQ.NBD)
C            of quasi dense NDME.GE.THRESMupdated then 
C            ndense(of elements adjacent to ME) sould be updated
            IF (NDENSE(ME).EQ.0) NDENSE(ME) =1
            IF (IVersion.EQ.1) THEN
C              ------------------------------------------------
C              place ME in the degree list of  DEGREE(ME)-(N+1)
C              NDENSE is not used in this case (simulate of 
C                      preprocessing )
C              ------------------------------------------------
              DEG = max (DEGREE(ME)-(TOTEL+1), 1)
            ELSE
C              -----------------------------------------
C              place ME in the degree list of NDENSE(ME)
C              -----------------------------------------
              DEG = NDENSE(ME)
            ENDIF
            DEGREE(ME) = DEG
            MINDEG = min(DEG,MINDEG)
            JNEXT = HEAD(DEG)
            IF (JNEXT.NE. 0) LAST (JNEXT) = ME
            NEXT(ME) = JNEXT
            HEAD(DEG) = ME
C           ------------------------------
C           process next quasi dense row
C           ------------------------------
            ME    = INEXT
            IF (ME.NE.0) THEN
              IF (DEGREE(ME).GT.(TOTEL+1) ) GOTO 51
            ENDIF
            HEAD (TOTEL) = ME
C           ---------------------------------------
C           update dense row selection strategy
C           -------------------------------------
C   
            IF (IVersion .EQ.1 ) THEN
             THRESM = TOTEL
            ELSE
             THRESM=max(THRESM*2,MINDEN+(MAXDEN-MINDEN)/2)
C            THRESM = max(THRESM*2, MINDEN*2)
             THRESM = min(THRESM,NBD)
             IF (THRESM.GE.NBD) THRESM=TOTEL
            ENDIF
            NBD    = NBED
C
            GOTO 30
          ENDIF
C         -------------------------------------------------------------
C         -------------------------------------------------------------
          IF (DEGREE(ME).EQ.TOTEL+1) THEN
C         we have only  exactly "dense" rows that we
C         amalgamate at the root node
             IF (NBD.NE.NBED) THEN
                write(6,*) ' Internal ERROR quasi dense rows remains'
                CALL MUMPS_ABORT()
             ENDIF
C          1/ Go through all
C          non absorbed elements (root of the subgraph) 
C          and absorb in ME
C          2/ perform mass elimination of all dense rows
C          RMK: we could compute sum(NVPIV(d)) to check if = NBD
           NELME    = -(NEL+1)
           DO 59 X=1,N
            IF ((PE(X).GT.0_8) .AND. (ELEN(X).LT.0)) THEN
C            X is an unabsorbed element
             PE(X) = int(-ME,8)
C            W(X) = 0 could be suppressed ?? check it
            ELSEIF (DEGREE(X).EQ.TOTEL+1) THEN
C            X is a dense row, absorb it in ME (mass elimination)
             NEL   = NEL + NV(X)
             PE(X) = int(-ME,8)
             ELEN(X) = 0
             NV(X) = 0
            ENDIF
   59      CONTINUE
C          ME is the root node
           ELEN(ME) = NELME
           NV(ME)   = NBD
           PE(ME)   = 0_8
           IF (NEL.NE.TOTEL) THEN
            write(6,*) 'Internal ERROR 2 detected in QAMD'
            write(6,*) ' NEL not equal to N: N, NEL =',N,NEL
            CALL MUMPS_ABORT()
           ENDIF
           GOTO 265
          ENDIF
        ENDIF
C       -------------------------------------------------------------
C       me represents the elimination of pivots nel+1 to nel+nv(me).
C       place me itself as the first in this set.  It will be moved
C       to the nel+nv(me) position when the permutation vectors are
C       computed.
C       -------------------------------------------------------------
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
        NDENSE(ME) = 0
C=======================================================================
C  CONSTRUCT NEW ELEMENT
C=======================================================================
C       -------------------------------------------------------------
C       At this point, me is the pivotal supervariable.  It will be
C       converted into the current element.  Scan list of the
C       pivotal supervariable, me, setting tree pointers and
C       constructing new list of supervariables for the new element,
C       me.  p is a pointer to the current position in the old list.
C       -------------------------------------------------------------
C       flag the variable "me" as being in Lme by negating nv (me)
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
C         ----------------------------------------------------------
C         construct the new element in place
C         ----------------------------------------------------------
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + int(LEN (ME) - 1,8)
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
C             ----------------------------------------------------
C             i is a principal variable not yet placed in Lme.
C             store i in new list
C             ----------------------------------------------------
              DEGME = DEGME + NVI
C             flag i as being in Lme by negating nv (i)
              NV (I) = -NVI
              PME2 = PME2 + 1_8
              IW (PME2) = I
C             ----------------------------------------------------
C             remove variable i from degree list.
C             ----------------------------------------------------
C             only done for non "dense" rows
              IF (DEGREE(I).LE.TOTEL) THEN
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                 NEXT (ILAST) = INEXT
              ELSE
C               i is at the head of the degree list
                 HEAD (DEGREE (I)) = INEXT
              ENDIF
              ELSE
               NDENSE(ME) = NDENSE(ME) + NVI
              ENDIF
            ENDIF
   60     CONTINUE
C         this element takes no new memory in iw:
          NEWMEM = 0
        ELSE
C         ----------------------------------------------------------
C         construct the new element in empty space, iw (pfree ...)
C         ----------------------------------------------------------
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          KNT1_UPDATED = 0
          DO 120 KNT1 = 1, ELENME + 1
            KNT1_UPDATED = KNT1_UPDATED +1
            IF (KNT1 .GT. ELENME) THEN
C             search the supervariables in me.
              E = ME
              PJ = P
              LN = SLENME
            ELSE
C             search the elements in me.
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
C           -------------------------------------------------------
C           search for different supervariables and add them to the
C           new list, compressing when necessary. this loop is
C           executed once for each element in the list and once for
C           all the supervariables in the list.
C           -------------------------------------------------------
            KNT2_UPDATED = 0
            DO 110 KNT2 = 1, LN
              KNT2_UPDATED = KNT2_UPDATED+1
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
C               -------------------------------------------------
C               compress iw, if necessary
C               -------------------------------------------------
                IF (PFREE .GT. IWLEN) THEN
C                 prepare for compressing iw by adjusting
C                 pointers and lengths so that the lists being
C                 searched in the inner and outer loops contain
C                 only the remaining entries.
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1_UPDATED
C                 Reset KNT1_UPDATED in case of recompress 
C                 at same iteration of the loop 120
                  KNT1_UPDATED = 0
C                 Check if anything left in supervariable ME
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0_8
                  PE (E) = PJ
                  LEN (E) = LN - KNT2_UPDATED
C                 Reset KNT2_UPDATED in case of recompress 
C                 at same iteration of the loop 110
                  KNT2_UPDATED = 0
C                 Check if anything left in element E
                  IF (LEN (E) .EQ. 0) PE (E) = 0_8
                  NCMPA = NCMPA + 1
C                 store first item in pe
C                 set first entry to -item
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = int(IW (PN),8)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
C                 psrc/pdst point to source/destination
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
C                 while loop:
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
C                   search for next negative entry
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = int(PE (J))
                      PE (J) = PDST
                      PDST = PDST + 1_8
C                     copy from source to destination
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
C                 move the new partially-constructed element
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
C               -------------------------------------------------
C               i is a principal variable not yet placed in Lme
C               store i in new list
C               -------------------------------------------------
                DEGME = DEGME + NVI
C               flag i as being in Lme by negating nv (i)
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
C               -------------------------------------------------
C               remove variable i from degree link list
C               -------------------------------------------------
C             only done for non "dense" rows
                IF (DEGREE(I).LE.TOTEL) THEN
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                   NEXT (ILAST) = INEXT
                ELSE
C                 i is at the head of the degree list
                   HEAD (DEGREE (I)) = INEXT
                ENDIF
                ELSE
                 NDENSE(ME) = NDENSE(ME) + NVI
                ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
C             set tree pointer and flag to indicate element e is
C             absorbed into new element me (the parent of e is me)
              PE (E) = int(-ME,8)
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1_8
C         this element takes newmem new memory in iw (possibly zero)
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = max (MAXMEM, MEM)
        ENDIF
C       -------------------------------------------------------------
C       me has now been converted into an element in iw (pme1..pme2)
C       -------------------------------------------------------------
C       degme holds the external degree of new element
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = int(PME2 - PME1 + 1_8)
C       -------------------------------------------------------------
C       make sure that wflg is not too large.  With the current
C       value of wflg, wflg+n must not cause integer overflow
C       -------------------------------------------------------------
        IF (WFLG .GT. MAXINT_N) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
C=======================================================================
C  COMPUTE (w (e) - wflg) = |Le\Lme| FOR ALL ELEMENTS
C=======================================================================
C       -------------------------------------------------------------
C       Scan 1:  compute the external degrees of previous elements
C       with respect to the current element.  That is:
C            (w (e) - wflg) = |Le \ Lme|
C       for each element e that appears in any supervariable in Lme.
C       The notation Le refers to the pattern (list of
C       supervariables) of a previous element e, where e is not yet
C       absorbed, stored in iw (pe (e) + 1 ... pe (e) + iw (pe (e))).
C       The notation Lme refers to the pattern of the current element
C       (stored in iw (pme1..pme2)).   If (w (e) - wflg) becomes
C       zero, then the element e will be absorbed in scan 2.
C       -------------------------------------------------------------
        DO 150 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).GT.TOTEL) GOTO 150
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
C           note that nv (i) has been negated to denote i in Lme:
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + int(ELN - 1,8)
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
C               unabsorbed element e has been seen in this loop
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
C               e is an unabsorbed element
C               this is the first we have seen e in all of Scan 1
                WE = DEGREE (E) + WNVI - NDENSE(E)
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
C=======================================================================
C  DEGREE UPDATE AND ELEMENT ABSORPTION
C=======================================================================
C       -------------------------------------------------------------
C       Scan 2:  for each i in Lme, sum up the degree of Lme (which
C       is degme), plus the sum of the external degrees of each Le
C       for the elements e appearing within i, plus the
C       supervariables in i.  Place i in hash list.
C       -------------------------------------------------------------
        DO 180 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).GT.TOTEL) GOTO 180
          P1 = PE (I)
          P2 = P1 + int(ELEN (I) - 1,8)
          PN = P1
          HASH = 0_8
          DEG = 0
C         ----------------------------------------------------------
C         scan the element list associated with supervariable i
C         ----------------------------------------------------------
          DO 160 P = P1, P2
            E = IW (P)
C           dext = | Le \ Lme |
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
#if defined (NOAGG5)
C        ------------------------------
C        suppress aggressive absorption
C        ------------------------------
            ELSE IF (DEXT .EQ. 0) THEN
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
#else
C
C        ------------------------------
C        try aggressive absorption 
C         when possible
C
            ELSE IF ((DEXT .EQ. 0) .AND.
     &                (NDENSE(ME).EQ.NBD)) THEN
C             aggressive absorption: e is not adjacent to me, but
C             |Le(G') \ Lme(G')| is 0 and all dense rows
C             are in me, so absorb it into me
                PE (E) = int(-ME,8)
                W (E)  = 0
            ELSE IF (DEXT.EQ.0) THEN
                  IW(PN) = E
                  PN     = PN+1
                  HASH   = HASH + int(E,kind=8)
#endif
            ENDIF
  160     CONTINUE
C         count the number of elements in i (including me):
          ELEN (I) = int(PN - P1 + 1)
C         ----------------------------------------------------------
C         scan the supervariables in the list associated with i
C         ----------------------------------------------------------
          P3 = PN
          DO 170 P = P2 + 1, P1 + int(LEN (I) - 1,8)
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
C             j is unabsorbed, and not in Lme.
C             add to degree and add to new list
C             add degree only of non-dense rows.
              IF (DEGREE(J).LE.TOTEL) DEG=DEG+NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + int(J,kind=8)
            ENDIF
  170     CONTINUE
C         ----------------------------------------------------------
C         update the degree and check for mass elimination
C         ----------------------------------------------------------
#if defined (NOAGG5)
          IF (DEG.EQ.0.AND.(NDENSE(ME).EQ.NBD).AND.(ELEN(I).GT.1)) THEN
C         When mass elimination will be performed then
C         absorb in ME all element adjacent to I
                   P1 = PE (I)
C                  exclude ME --> -2
                   P2 = P1 + int(ELEN (I),8) - 2_8
                   DO P =P1,P2  
                     E      = IW(P)
                     PE (E) = int(-ME,8)
                     W (E)  = 0
                   ENDDO
          ENDIF
C              .... Ready for mass elimination
#endif
          IF ((DEG .EQ. 0).AND.(NDENSE(ME).EQ.NBD)) THEN
C           -------------------------------------------------------
C           mass elimination
C           -------------------------------------------------------
C           There is nothing left of this node except for an
C           edge to the current pivot element.  elen (i) is 1,
C           and there are no variables adjacent to node i.
C           Absorb i into the current pivot element, me.
            PE (I) = int(-ME,8)
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
C           -------------------------------------------------------
C           update the upper-bound degree of i
C           -------------------------------------------------------
C           the following degree does not yet include the size
C           of the current element, which is added later:
            DEGREE(I) = min (DEG+NBD-NDENSE(ME), 
     &                       DEGREE(I))
C           -------------------------------------------------------
C           add me to the list for i
C           -------------------------------------------------------
C           move first supervariable to end of list
            IW (PN) = IW (P3)
C           move first element to end of element part of list
            IW (P3) = IW (P1)
C           add new element to front of list.
            IW (P1) = ME
C           store the new length of the list in len (i)
            LEN (I) = int(PN - P1 + 1)
C           -------------------------------------------------------
C           place in hash bucket.  Save hash key of i in last (i).
C           -------------------------------------------------------
            HASH = mod (HASH, HMOD) + 1_8
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
C             the degree list is empty, hash head is -j
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
C             degree list is not empty
C             use last (head (hash)) as hash head
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = int(HASH,kind=kind(LAST))
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
C       -------------------------------------------------------------
C       Clear the counter array, w (...), by incrementing wflg.
C       -------------------------------------------------------------
        DMAX = max (DMAX, DEGME)
        WFLG = WFLG + DMAX
C       make sure that wflg+n does not cause integer overflow
        IF (WFLG .GT. MAXINT_N) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
C       at this point, w (1..n) .lt. wflg holds
C=======================================================================
C  SUPERVARIABLE DETECTION
C=======================================================================
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF ( (NV(I).LT.0) .AND. (DEGREE(I).LE.TOTEL) ) THEN
C           only done for nondense rows
C           i is a principal variable in Lme
C           -------------------------------------------------------
C           examine all hash buckets with 2 or more variables.  We
C           do this by examing all unique hash keys for super-
C           variables in the pattern Lme of the current element, me
C           -------------------------------------------------------
            HASH = int(LAST (I),kind=8)
C           let i = head of hash bucket, and empty the hash bucket
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
C             degree list is empty
              I = -J
              HEAD (HASH) = 0
            ELSE
C             degree list is not empty, restore last () of head
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
C           while loop:
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
C             ----------------------------------------------------
C             this bucket has one or more variables following i.
C             scan all of them to see if i can absorb any entries
C             that follow i in hash bucket.  Scatter i into w.
C             ----------------------------------------------------
              LN = LEN (I)
              ELN = ELEN (I)
C             do not flag the first element in the list (me)
              DO 210 P = PE (I) + 1, PE (I) + int(LN - 1,8)
                W (IW (P)) = WFLG
  210         CONTINUE
C             ----------------------------------------------------
C             scan every other entry j following i in bucket
C             ----------------------------------------------------
              JLAST = I
              J = NEXT (I)
C             while loop:
  220         CONTINUE
              IF (J .NE. 0) THEN
C               -------------------------------------------------
C               check if j and i have identical nonzero pattern
C               -------------------------------------------------
C               jump if i and j do not have same size data structure
                IF (LEN (J) .NE. LN) GO TO 240
C               jump if i and j do not have same number adj elts
                IF (ELEN (J) .NE. ELN) GO TO 240
C               do not flag the first element in the list (me)
                DO 230 P = PE (J) + 1, PE (J) + int(LN - 1,8)
C                 jump if an entry (iw(p)) is in j but not in i
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
C               -------------------------------------------------
C               found it!  j can be absorbed into i
C               -------------------------------------------------
                PE (J) = int(-I,8)
C               both nv (i) and nv (j) are negated since they
C               are in Lme, and the absolute values of each
C               are the number of variables in i and j:
                NV (I) = NV (I) + NV (J)
                NV (J) = 0
                ELEN (J) = 0
C               delete j from hash bucket
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
C               -------------------------------------------------
  240           CONTINUE
C               j cannot be absorbed into i
C               -------------------------------------------------
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
C             ----------------------------------------------------
C             no more variables can be absorbed into i
C             go to next i in bucket and clear flag array
C             ----------------------------------------------------
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
C=======================================================================
C  RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVAR. FROM ELEMENT
C=======================================================================
        P = PME1
        NLEFT = TOTEL - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
C           i is a principal variable in Lme
C           restore nv (i) to signify that i is principal
            NV (I) = NVI
            IF (DEGREE(I).LE.TOTEL) THEN
C           -------------------------------------------------------
C           compute the external degree (add size of current elem)
C           -------------------------------------------------------
            DEG = min (DEGREE (I)+ DEGME - NVI, NLEFT - NVI)
            DEGREE (I) = DEG
            IDENSE = .FALSE.
C           
       IF ( (IVersion .NE. 1).AND. (THRESM.GE.0)) THEN
C           -------------------
C           Dense row detection
C           -------------------
C           DEGME is exact external degree of pivot ME |Le\Ve|, 
C           DEG is is approx external degree of I
C           Relaxed dense row selection based on:
C            1/ We want to avoid selecting dense  rows that are
C               almost completely represented by adj(ME)
C            1/ its density in reduced matrix and 
          IF (DEG+NVI .GE. THRESM) THEN
             IF (THRESM.EQ.TOTEL) THEN
C             We must be sure that I is exactly dense in reduced matrix
                IF ((ELEN(I).LE.2) .AND. ((DEG+NVI).EQ.NLEFT) ) THEN
C              DEG approximation is exact and I is dense 
                   DEGREE(I) = TOTEL+1
                   IDENSE = .TRUE.
                ENDIF
             ELSE
C             relaxed dense row detection
                IDENSE = .TRUE.
                IF ((ELEN(I).LE.2).AND.((DEG+NVI).EQ.NLEFT) ) THEN
                   DEGREE(I) = TOTEL+1
                ELSE
                   DEGREE(I) = TOTEL+1+DEGREE(I)
                ENDIF
             ENDIF
          ENDIF
          IF (IDENSE) THEN
C            update NDENSE of all elements in the list of element
C            adjacent to I (including ME).
             P1 = PE(I)
             P2 = P1 + int(ELEN(I) - 1,8)
             IF (P2.GE.P1) THEN
                DO 264 PJ=P1,P2
                   E= IW(PJ)
                   NDENSE (E) = NDENSE(E) + NVI
 264            CONTINUE
             ENDIF
C            insert I in the list of dense rows
             NBD = NBD+NVI
             DEG = TOTEL
             IF (DEGREE(I).EQ.TOTEL+1) THEN
c              insert I at the end of the list
                NBED = NBED +NVI
                IF (LASTD.EQ.0) THEN
C                degree list is empty
                   LASTD     = I
                   HEAD(DEG) = I
                   NEXT(I)   = 0
                   LAST(I)   = 0
                ELSE
                   NEXT(LASTD) = I
                   LAST(I)     = LASTD
                   LASTD       = I
                   NEXT(I)     = 0
                ENDIF
             ELSE
C              insert I at the beginning of the list
                INEXT = HEAD(DEG)
                IF (INEXT .NE. 0) LAST (INEXT) = I
                NEXT (I) = INEXT
                HEAD (DEG) = I
                LAST(I)    = 0
                IF (LASTD.EQ.0) LASTD=I
             ENDIF
C            end of IDENSE=true
          ENDIF
C           end of THRESM>0
       ENDIF
C             
       IF (.NOT.IDENSE) THEN
C           -------------------------------------------------------
C           place the supervariable at the head of the degree list
C           -------------------------------------------------------
          INEXT = HEAD (DEG)
          IF (INEXT .NE. 0) LAST (INEXT) = I
          NEXT (I) = INEXT
          LAST (I) = 0
          HEAD (DEG) = I
       ENDIF
C           -------------------------------------------------------
C           save the new degree, and find the minimum degree
C           -------------------------------------------------------
       MINDEG = min (MINDEG, DEG)
            ENDIF
C           -------------------------------------------------------
C           place the supervariable in the element pattern
C           -------------------------------------------------------
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
C=======================================================================
C  FINALIZE THE NEW ELEMENT
C=======================================================================
        NV (ME) = NVPIV + DEGME
C       nv (me) is now the degree of pivot (including diagonal part)
C       save the length of the list for the new element me
        LEN (ME) = int(P - PME1)
        IF (LEN (ME) .EQ. 0) THEN
C         there is nothing left of the current pivot element
          PE (ME) = 0_8
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
C         element was not constructed in place: deallocate part
C         of it (final size is less than or equal to newmem,
C         since newly nonprincipal variables have been removed).
          PFREE = P
          MEM = MEM - NEWMEM + int(LEN (ME),8)
        ENDIF
C=======================================================================
C       END WHILE (selecting pivots)
      GO TO 30
      ENDIF
C=======================================================================
  265 CONTINUE
C=======================================================================
C  COMPUTE THE PERMUTATION VECTORS and update TREE
C=======================================================================
C     ----------------------------------------------------------------
C     The time taken by the following code is O(n).  At this
C     point, elen (e) = -k has been done for all elements e,
C     and elen (i) = 0 has been done for all nonprincipal
C     variables i.  At this point, there are no principal
C     supervariables left, and all elements are absorbed.
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
C     compute the ordering of unordered nonprincipal variables
C     ----------------------------------------------------------------
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
C         ----------------------------------------------------------
C         i is an un-ordered row.  Traverse the tree from i until
C         reaching an element, e.  The element, e, was the
C         principal supervariable of i and all nodes in the path
C         from i to when e was selected as pivot.
C         ----------------------------------------------------------
          J = int(-PE (I))
C         while (j is a variable) do:
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = int(-PE (J))
              GO TO 270
            ENDIF
            E = J
C           ----------------------------------------------------------
C           get the current pivot ordering of e
C           ----------------------------------------------------------
            K = -ELEN (E)
C           ----------------------------------------------------------
C           traverse the path again from i to e, and compress the
C           path (all nodes point to e).  Path compression allows
C           this code to compute in O(n) time.  Order the unordered
C           nodes in the path, and place the element e at the end.
C           ----------------------------------------------------------
            J = I
C           while (j is a variable) do:
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = int(-PE (J))
              PE (J) = int(-E,8)
              IF (ELEN (J) .EQ. 0) THEN
C               j is an unordered row
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
C         leave elen (e) negative, so we know it is an element
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
      IF (COMPUTE_PERM) THEN
C     ----------------------------------------------------------------
C     reset the inverse permutation (elen (1..n)) to be positive,
C     and compute the permutation (last (1..n)).
C     ----------------------------------------------------------------
      IF(COMPRESS) THEN
        LAST(1:N) = 0
        HEAD(1:TOTEL-N)=0  
        DO I = 1, N
          K = abs (ELEN (I))
          IF ( K <= N ) THEN
            LAST (K) = I
          ELSE
            HEAD(K-N)=I
          ENDIF
        ENDDO
        I = 1
        DO K = 1, N
          IF(LAST (K) .NE. 0) THEN
            LAST(I) = LAST(K)
            ELEN(LAST(K)) = I
            I = I + 1
          ENDIF
        ENDDO
        DO K = N+1, TOTEL
          IF (HEAD(K-N) .NE. 0) THEN
            LAST(I)=HEAD(K-N)
            ELEN(HEAD(K-N)) = I
            I = I + 1
          ENDIF
        END DO
      ELSE
         DO 300 I = 1, N
            K = abs (ELEN (I))
            LAST (K) = I
            ELEN (I) = K
 300     CONTINUE
      ENDIF
C=======================================================================
C      END OF COMPUTING PERMUTATIONS
C=======================================================================
       ENDIF
C=======================================================================
C  RETURN THE MEMORY USAGE IN IW
C=======================================================================
C     If maxmem is less than or equal to iwlen, then no compressions
C     occurred, and iw (maxmem+1 ... iwlen) was unused.  Otherwise
C     compressions did occur, and iwlen would have had to have been
C     greater than or equal to maxmem for no compressions to occur.
C     Return the value of maxmem in the pfree argument.
      PFREE = MAXMEM
C===============================
C     Save PE in PARENT array
      DO I=1,N
       PARENT(I) = int(PE(I))
      ENDDO
C===============================
      RETURN
      END SUBROUTINE MUMPS_QAMD
C-----------------------------------------------------------------------
C MUMPS_CST_AMF: modified version of MUMPS_HAMF4 routine 
C implementing constraint minimum fill-in based 
C ordering.
C Written by Stephane Pralet iduring his post-doctorate at INPT-IRIT 
C (Oct. 2004- Oct. 2005)
C
C   Restrictive integer 64 bit variant :
C   it is assumed that IW array size can exceed 32-bit integer
C
      SUBROUTINE MUMPS_CST_AMF (N, NBBUCK, 
     &     IWLEN, PE, PFREE, LEN, IW, NV, ELEN,
     &     LAST, NCMPA, DEGREE, WF, NEXT, W, HEAD,
     &     CONSTRAINT,THESON, PARENT)
      IMPLICIT NONE
C
C Parameters
C    Input not modified
      INTEGER, INTENT(IN)    :: N, NBBUCK
      INTEGER(8), INTENT(IN) :: IWLEN
C     Input undefined on output 
      INTEGER, INTENT(INOUT)  :: LEN(N), IW(IWLEN)
C     NV meaningful as input to encode compressed graphs
      INTEGER, INTENT(INOUT)  :: NV(N)
C 
C     Output only 
      INTEGER, INTENT(OUT)   :: NCMPA
      INTEGER, INTENT(OUT)   :: ELEN(N), LAST(N), PARENT(N)
C 
C     Input/output
      INTEGER(8), INTENT(INOUT) :: PFREE
      INTEGER(8), INTENT(INOUT) :: PE(N)
C 
C     Internal Workspace only
C       Min fill approximation one extra array of size NBBUCK+2 
C       is also needed
      INTEGER     :: NEXT(N), DEGREE(N), W(N)
      INTEGER     :: HEAD(0:NBBUCK+1), WF(N)
C
C  Comments on the OUTPUT:
C  ----------------------
C  Let V= V0 U V1 the nodes of the initial graph (|V|=n). 
C  The assembly tree corresponds to the tree 
C    of the supernodes (or supervariables). Each node of the 
C    assembly tree is then composed of one principal variable 
C    and a list of secondary variables. The list of 
C    variable of a node (principal + secondary variables) then 
C    describes the structure of the diagonal bloc of the 
C    supernode. 
C  The elimination tree denotes the tree of all the variables(=node) and 
C    is therefore of order n.
C
C  The arrays NV(N) and PE(N) give a description of the 
C  assembly tree. 
C  
C   1/ Description of array nv(N) (on OUPUT)
C    nv(i)=0 i is a secondary variable 
C    N+1> nv(i) >0 i is a principal variable, nv(i) holds the 
C                  the number of elements in column i of L (true degree of i)
C
C   2/ Description of array PE(N) (on OUPUT)
C       pe(i) = -(father of variable/node i) in the elimination tree:
C       If nv (i) .gt. 0, then i represents a node in the assembly tree,
C       and the parent of i is -pe (i), or zero if i is a root.
C       If nv (i) = 0, then (i,-pe (i)) represents an edge in a
C       subtree, the root of which is a node in the assembly tree.
C   
C   3/ Example:
C      Let If be a root node father of Is in the assembly tree. 
C      If is the principal 
C      variable of the node If and let If1, If2, If3 be the secondary variables
C      of node If.
C      Is is the principal 
C      variable of the node Is and let Is1, Is2 be the secondary variables
C      of node Is.
C      
C      THEN: 
C        NV(If1)=NV(If2)=NV(If3) = 0  (secondary variables)
C        NV(Is1)=NV(Is2) = 0  (secondary variables)
C        NV(If) > 0  ( principal variable)
C        NV(Is) > 0  ( principal variable)
C        PE(If)  = 0 (root node)
C        PE(Is)  = -If (If is the father of Is in the assembly tree)
C        PE(If1)=PE(If2)=PE(If3)= -If  ( If is the principal variable)
C        PE(Is1)=PE(Is2)= -Is  ( Is is the principal variable)
C      
C
C
C HALOAMD_V1: (September 1997)
C **********
C Initial version designed to experiment the numerical (fill-in) impact 
C of taking into account the halo. This code should be able 
C to experiment no-halo, partial halo, complete halo.
C DATE: September 17th 1997
C
C HALOAMD is designed to process a gragh composed of two types
C            of nodes, V0 and V1, extracted from a larger gragh. 
C            V0^V1 = {}, 
C       
C            We used Min. degree heuristic to order only 
C            nodes in V0, but the adjacency to nodes
C            in V1 is taken into account during ordering.
C            Nodes in V1 are odered at last.
C            Adjacency between nodes of V1 need not be provided,
C            however |len(i)| must always corresponds to the number of 
C            edges effectively provided in the adjacency list of i.
C          On input :
c          ********
C            Nodes INODE in V1 are flagged with len(INODE) = -degree 
C            modif version HALO V3 (August 1998): 
C                                 if len(i) =0 and i \in V1 then 
C                                 len(i) must be set on input to -N-1
C          ERROR return (negative values in ncmpa)
C          ************
C            negative value in ncmpa indicates an error detected 
C               by HALOAMD.
C
C            The graph provided MUST follow the rule:
C             if (i,j) is an edge in the gragh then 
C             j must be in the adjacency list of i AND 
C             i must be in the adjacency list of j.
C    REMARKS
C    -------
C        
C        1/  Providing edges between nodes of V1 should not 
C            affect the final ordering, only the amount of edges 
C            of the halo should effectively affect the solution.
C            This code should work in the following cases:
C              1/ halo not provided
C              2/ halo partially provided
C              3/ complete halo
C              4/ complete halo+interconnection between nodes of V1.
C
C              1/ should run and provide identical results (w.r.t to current 
C               implementation of AMD in SCOTCH).
C             3/ and 4 should provide identical results.
C
C        2/ All modifications of the AMD initial code are indicated
C           with begin HALO .. end HALO
C
C            
C   Ordering of nodes in V0 is based on  approximate minimum
C       fill-in heuristic.
C   
C-----------------------------------------------------------------------
C begin CONSTRAINT
C CONSTRAINT(I) >= 0 : I can be selected
C                < 0 : I cannot be selected
C                > 0 : I release CONSTRAINT(I) 
C THESON(I) = 0 : I is a leaf in the supervariable representation
C THESON(I) > I : THESON(I) belongs to the same supervariable as I
C  Parameters:
      INTEGER, INTENT(INOUT) :: CONSTRAINT(N)
      INTEGER, INTENT(out) :: THESON(N)
      INTEGER PREV,TOTO
C end CONSTRAINT
C-----------------------------------------------------------------------
C INPUT ARGUMENTS (unaltered):
C-----------------------------------------------------------------------
C n:    The matrix order.
C
C       Restriction:  n .ge. 1
C iwlen:        The length of iw (1..iwlen).  On input, the matrix is
C       stored in iw (1..pfree-1).  However, iw (1..iwlen) should be
C       slightly larger than what is required to hold the matrix, at
C       least iwlen .ge. pfree + n is recommended.  Otherwise,
C       excessive compressions will take place.
C       *** We do not recommend running this algorithm with ***
C       ***      iwlen .lt. pfree + n.                      ***
C       *** Better performance will be obtained if          ***
C       ***      iwlen .ge. pfree + n                       ***
C       *** or better yet                                   ***
C       ***      iwlen .gt. 1.2 * pfree                     ***
C       *** (where pfree is its value on input).            ***
C       The algorithm will not run at all if iwlen .lt. pfree-1.
C
C       Restriction: iwlen .ge. pfree-1
C-----------------------------------------------------------------------
C INPUT/OUPUT ARGUMENTS:
C-----------------------------------------------------------------------
C pe:   On input, pe (i) is the index in iw of the start of row i, or
C       zero if row i has no off-diagonal non-zeros.
C
C       During execution, it is used for both supervariables and
C       elements:
C
C       * Principal supervariable i:  index into iw of the
C               description of supervariable i.  A supervariable
C               represents one or more rows of the matrix
C               with identical nonzero pattern.
C       * Non-principal supervariable i:  if i has been absorbed
C               into another supervariable j, then pe (i) = -j.
C               That is, j has the same pattern as i.
C               Note that j might later be absorbed into another
C               supervariable j2, in which case pe (i) is still -j,
C               and pe (j) = -j2.
C       * Unabsorbed element e:  the index into iw of the description
C               of element e, if e has not yet been absorbed by a
C               subsequent element.  Element e is created when
C               the supervariable of the same name is selected as
C               the pivot.
C       * Absorbed element e:  if element e is absorbed into element
C               e2, then pe (e) = -e2.  This occurs when the pattern of
C               e (that is, Le) is found to be a subset of the pattern
C               of e2 (that is, Le2).  If element e is "null" (it has
C               no nonzeros outside its pivot block), then pe (e) = 0.
C
C       On output, pe holds the assembly tree/forest, which implicitly
C       represents a pivot order with identical fill-in as the actual
C       order (via a depth-first search of the tree).
C
C       On output:
C       If nv (i) .gt. 0, then i represents a node in the assembly tree,
C       and the parent of i is -pe (i), or zero if i is a root.
C       If nv (i) = 0, then (i,-pe (i)) represents an edge in a
C       subtree, the root of which is a node in the assembly tree.
C       On output:  (PE is copied on output into PARENT array)
C
C pfree:        On input, the matrix is stored in iw (1..pfree-1) and
C       the rest of the array iw is free.
C       During execution, additional data is placed in iw, and pfree
C       is modified so that components  of iw from pfree are free.
C       On output, pfree is set equal to the size of iw that
C       would have been needed for no compressions to occur.  If
C       ncmpa is zero, then pfree (on output) is less than or equal to
C       iwlen, and the space iw (pfree+1 ... iwlen) was not used.
C       Otherwise, pfree (on output) is greater than iwlen, and all the
C       memory in iw was used.
C
C nv: On input, encoding of compressed graph:
C        if NV(1) = -1 then graph is not compressed otherwise
C        NV(I) holds the weight of node I. 
C       During execution, abs (nv (i)) is equal to the number of rows
C       that are represented by the principal supervariable i.  If i is
C       a nonprincipal variable, then nv (i) = 0.  Initially,
C       nv (i) = 1 for all i.  nv (i) .lt. 0 signifies that i is a
C       principal variable in the pattern Lme of the current pivot
C       element me.  On output, nv (e) holds the true degree of element
C       e at the time it was created (including the diagonal part).
C begin HALO
C       On output, nv(I) can be used to find node in set V1.
C       Not true anymore : ( nv(I) = N+1 characterizes nodes in V1. 
C                 instead nodes in V1 are considered as a dense root node )
C end HALO
C-----------------------------------------------------------------------
C INPUT/MODIFIED (undefined on output):
C-----------------------------------------------------------------------
C len:  On input, len (i) 
C           positive or null (>=0) : i \in V0 and 
C                     len(i) holds the number of entries in row i of the
C                     matrix, excluding the diagonal.  
C           negative (<0) : i \in V1, and 
C                     -len(i) hold the number of entries in row i of the
C                     matrix, excluding the diagonal.
C                     len(i) = - | Adj(i) | if i \in V1                    
C                              or -N -1 if  | Adj(i) | = 0 and i \in V1 
C       The contents of len (1..n)
C       are undefined on output.
C iw:   On input, iw (1..pfree-1) holds the description of each row i
C       in the matrix.  The matrix must be symmetric, and both upper
C       and lower triangular parts must be present.  The diagonal must
C       not be present.  Row i is held as follows:
C
C               len (i):  the length of the row i data structure
C               iw (pe (i) ... pe (i) + len (i) - 1):
C                       the list of column indices for nonzeros
C                       in row i (simple supervariables), excluding
C                       the diagonal.  All supervariables start with
C                       one row/column each (supervariable i is just
C                       row i).
C               if len (i) is zero on input, then pe (i) is ignored
C               on input.
C
C               Note that the rows need not be in any particular order,
C               and there may be empty space between the rows.
C
C       During execution, the supervariable i experiences fill-in.
C       This is represented by placing in i a list of the elements
C       that cause fill-in in supervariable i:
C
C               len (i):  the length of supervariable i
C               iw (pe (i) ... pe (i) + elen (i) - 1):
C                       the list of elements that contain i.  This list
C                       is kept short by removing absorbed elements.
C               iw (pe (i) + elen (i) ... pe (i) + len (i) - 1):
C                       the list of supervariables in i.  This list
C                       is kept short by removing nonprincipal
C                       variables, and any entry j that is also
C                       contained in at least one of the elements
C                       (j in Le) in the list for i (e in row i).
C
C       When supervariable i is selected as pivot, we create an
C       element e of the same name (e=i):
C
C               len (e):  the length of element e
C               iw (pe (e) ... pe (e) + len (e) - 1):
C                       the list of supervariables in element e.
C
C       An element represents the fill-in that occurs when supervariable
C       i is selected as pivot (which represents the selection of row i
C       and all non-principal variables whose principal variable is i).
C       We use the term Le to denote the set of all supervariables
C       in element e.  Absorbed supervariables and elements are pruned
C       from these lists when computationally convenient.
C
C       CAUTION:  THE INPUT MATRIX IS OVERWRITTEN DURING COMPUTATION.
C       The contents of iw are undefined on output.
C-----------------------------------------------------------------------
C OUTPUT (need not be set on input):
C-----------------------------------------------------------------------
C elen: See the description of iw above.  At the start of execution,
C       elen (i) is set to zero.  During execution, elen (i) is the
C       number of elements in the list for supervariable i.  When e
C       becomes an element, elen (e) = -nel is set, where nel is the
C       current step of factorization.  elen (i) = 0 is done when i
C       becomes nonprincipal.
C
C       For variables, elen (i) .ge. 0 holds until just before the
C       permutation vectors are computed.  For elements,
C       elen (e) .lt. 0 holds.
C
C       On output elen (1..n) holds the inverse permutation (the same
C       as the 'INVP' argument in Sparspak).  That is, if k = elen (i),
C       then row i is the kth pivot row.  Row i of A appears as the
C       (elen(i))-th row in the permuted matrix, PAP^T.
C last: In a degree list, last (i) is the supervariable preceding i,
C       or zero if i is the head of the list.  In a hash bucket,
C       last (i) is the hash key for i.  last (head (hash)) is also
C       used as the head of a hash bucket if head (hash) contains a
C       degree list (see head, below).
C
C       On output, last (1..n) holds the permutation (the same as the
C       'PERM' argument in Sparspak).  That is, if i = last (k), then
C       row i is the kth pivot row.  Row last (k) of A is the k-th row
C       in the permuted matrix, PAP^T.
C ncmpa:        The number of times iw was compressed.  If this is
C       excessive, then the execution took longer than what could have
C       been.  To reduce ncmpa, try increasing iwlen to be 10% or 20%
C       larger than the value of pfree on input (or at least
C       iwlen .ge. pfree + n).  The fastest performance will be
C       obtained when ncmpa is returned as zero.  If iwlen is set to
C       the value returned by pfree on *output*, then no compressions
C       will occur.
C begin HALO
C        on output ncmpa <0 --> error detected during HALO_AMD:
C           error 1: ncmpa = -N , ordering was stopped.
C end HALO
C
C-----------------------------------------------------------------------
C LOCAL (not input or output - used only during execution):
C-----------------------------------------------------------------------
C degree:       If i is a supervariable, then degree (i) holds the
C       current approximation of the external degree of row i (an upper
C       bound).  The external degree is the number of nonzeros in row i,
C       minus abs (nv (i)) (the diagonal part).  The bound is equal to
C       the external degree if elen (i) is less than or equal to two.
C       We also use the term "external degree" for elements e to refer
C       to |Le \ Lme|.  If e is an element, then degree (e) holds |Le|,
C       which is the degree of the off-diagonal part of the element e
C       (not including the diagonal part).
C begin HALO
C       degree(I) = n+1 indicates that i belongs to V1
C end HALO
C
C head: head is used for degree lists.  head (deg) is the first
C       supervariable in a degree list (all supervariables i in a
C       degree list deg have the same approximate degree, namely,
C       deg = degree (i)).  If the list deg is empty then
C       head (deg) = 0.
C
C       During supervariable detection head (hash) also serves as a
C       pointer to a hash bucket.
C       If head (hash) .gt. 0, there is a degree list of degree hash.
C               The hash bucket head pointer is last (head (hash)).
C       If head (hash) = 0, then the degree list and hash bucket are
C               both empty.
C       If head (hash) .lt. 0, then the degree list is empty, and
C               -head (hash) is the head of the hash bucket.
C       After supervariable detection is complete, all hash buckets
C       are empty, and the (last (head (hash)) = 0) condition is
C       restored for the non-empty degree lists.
C next: next (i) is the supervariable following i in a link list, or
C       zero if i is the last in the list.  Used for two kinds of
C       lists:  degree lists and hash buckets (a supervariable can be
C       in only one kind of list at a time).
C w:    The flag array w determines the status of elements and
C       variables, and the external degree of elements.
C
C       for elements:
C          if w (e) = 0, then the element e is absorbed
C          if w (e) .ge. wflg, then w (e) - wflg is the size of
C               the set |Le \ Lme|, in terms of nonzeros (the
C               sum of abs (nv (i)) for each principal variable i that
C               is both in the pattern of element e and NOT in the
C               pattern of the current pivot element, me).
C          if wflg .gt. w (e) .gt. 0, then e is not absorbed and has
C               not yet been seen in the scan of the element lists in
C               the computation of |Le\Lme| in loop 150 below.
C
C       for variables:
C          during supervariable detection, if w (j) .ne. wflg then j is
C          not in the pattern of variable i
C
C       The w array is initialized by setting w (i) = 1 for all i,
C       and by setting wflg = 2.  It is reinitialized if wflg becomes
C       too large (to ensure that wflg+n does not cause integer
C       overflow).
C
C wf : integer array  used to store the already filled area of 
C      the variables adajcent to current pivot. 
C      wf is then used to update the score of variable i.
C
C-----------------------------------------------------------------------
C LOCAL INTEGERS:
C-----------------------------------------------------------------------
      INTEGER :: DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     &        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     &        LENJ, LN, ME, MINDEG, NEL, 
     &        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X,
     &        NBFLAG, NREAL, LASTD, NELME, WF3, WF4, N2, PAS
       INTEGER KNT1_UPDATED, KNT2_UPDATED
       INTEGER(8) :: MAXMEM, MEM, NEWMEM
       INTEGER   :: MAXINT_N
       INTEGER(8)::  HASH, HMOD
       DOUBLE PRECISION ::   RMF, RMF1 
       DOUBLE PRECISION ::   dummy
       INTEGER :: idummy
C deg:        the degree of a variable or element
C degme:      size, |Lme|, of the current element, me (= degree (me))
C dext:       external degree, |Le \ Lme|, of some element e
C dmax:       largest |Le| seen so far
C e:          an element
C elenme:     the length, elen (me), of element list of pivotal var.
C eln:        the length, elen (...), of an element list
C hash:       the computed value of the hash function
C hmod:       the hash function is computed modulo hmod = max (1,n-1)
C i:          a supervariable
C ilast:      the entry in a link list preceding i
C inext:      the entry in a link list following i
C j:          a supervariable
C jlast:      the entry in a link list preceding j
C jnext:      the entry in a link list, or path, following j
C k:          the pivot order of an element or variable
C knt1:       loop counter used during element construction
C knt2:       loop counter used during element construction
C knt3:       loop counter used during compression
C lenj:       len (j)
C ln:         length of a supervariable list
C maxint_n:   large integer to test risk of overflow on wflg
C maxmem:     amount of memory needed for no compressions
C me:         current supervariable being eliminated, and the
C                     current element created by eliminating that
C                     supervariable
C mem:        memory in use assuming no compressions have occurred
C mindeg:     current minimum degree
C nel:        number of pivots selected so far
C newmem:     amount of new memory needed for current pivot element
C nleft:      n - nel, the number of nonpivotal rows/columns remaining
C nvi:        the number of variables in a supervariable i (= nv (i))
C nvj:        the number of variables in a supervariable j (= nv (j))
C nvpiv:      number of pivots in current element
C slenme:     number of variables in variable list of pivotal variable
C we:         w (e)
C wflg:       used for flagging the w array.  See description of iw.
C wnvi:       wflg - nv (i)
C x:          either a supervariable or an element
C wf3:  off diagonal block area
C wf4:  diagonal block area
C mf : Minimum fill
C begin HALO
C nbflag:     number of flagged entries in the initial gragh.
C nreal :     number of entries on which ordering must be perfomed
C             (nreel = N- nbflag)
C nelme number of pivots selected when reaching the root
C lastd index of the last row in the list of dense rows
C end HALO
C-----------------------------------------------------------------------
C LOCAL POINTERS:
C-----------------------------------------------------------------------
      INTEGER(8) P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, 
     &           PN, PSRC
C             Any parameter (pe (...) or pfree) or local variable
C             starting with "p" (for Pointer) is an index into iw,
C             and all indices into iw use variables starting with
C             "p."  The only exception to this rule is the iwlen
C             input argument.
C p:          pointer into lots of things
C p1:         pe (i) for some variable i (start of element list)
C p2:         pe (i) + elen (i) -  1 for some var. i (end of el. list)
C p3:         index of first supervariable in clean list
C pdst:       destination pointer, for compression
C pend:       end of memory to compress
C pj:         pointer into an element or variable
C pme:        pointer into the current element (pme1...pme2)
C pme1:       the current element, me, is stored in iw (pme1...pme2)
C pme2:       the end of the current element
C pn:         pointer into a "clean" variable, also used to compress
C psrc:       source pointer, for compression
C-----------------------------------------------------------------------
C  FUNCTIONS CALLED:
C-----------------------------------------------------------------------
      INTRINSIC max, min, mod, huge
      INTEGER TOTEL
C=======================================================================
C  INITIALIZATIONS
C=======================================================================
C     HEAD (0:NBBUCK+1)
C begin HALO
C
C idummy holds the largest integer - 1
C dummy  = dble (idummy)
      idummy = huge(idummy) - 1
      dummy = dble(idummy)
C     variable with degree equal to N2 are in halo
C     bucket NBBUCK+1 used for HALO variables
      N2 = -NBBUCK-1
C end HALO
C Distance betweeen elements of the N, ..., NBBUCK entries of HEAD
C
C update done on 20 Feb 2002 (PAS>= 1)
      PAS = max((N/8), 1)
      WFLG = 2
      MAXINT_N=huge(WFLG)-N
      NCMPA = 0
      NEL = 0
      HMOD = int(max (1, NBBUCK-1),kind=8)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      MINDEG = 0
C
      NBFLAG = 0
      LASTD  = 0
      HEAD(0:NBBUCK+1) = 0
      DO 10 I = 1, N
         THESON(I) = 0
         LAST (I) = 0
C        NV (I) = 1
         W (I) = 1
         ELEN (I) = 0
   10 CONTINUE
      TOTEL = 0
      DO I=1,N
         IF (LEN(I).LT.0) THEN 
            DEGREE (I) = N2
            NBFLAG     = NBFLAG +1
            IF (LEN(I).EQ.-N-1) THEN
C     variable in V1 with empty adj list 
               LEN (I)    = 0
C     Because of compress, we force skipping this
C     entry which is anyway empty
               PE (I)     = 0_8
            ELSE
               LEN (I)    = - LEN(I)
            ENDIF
C       end HALO V3
         ELSE
            TOTEL = TOTEL + NV(I)
            DEGREE(I) = 0
            DO P= PE(I) , PE(I)+int(LEN(I)-1,8)
               DEGREE(I) = DEGREE(I) + NV(IW(P))
            ENDDO
         ENDIF
      ENDDO
C
C
C     number of entries to be ordered.
      NREAL = N - NBFLAG
C     ----------------------------------------------------------------
C     initialize degree lists and eliminate rows with no off-diag. nz.
C     ----------------------------------------------------------------
      DO 20 I = 1, N
        DEG = DEGREE (I)
        IF (DEG.EQ.N2) THEN
C            DEG = N2 (flagged variables are stored 
C                  in the degree list of NBBUCK + 1
C                  (safe: because max 
C                         max value of degree is NBBUCK)
C
             DEG = NBBUCK + 1
             IF (LASTD.EQ.0) THEN
C              degree list is empty
               LASTD     = I
               HEAD(DEG) = I
               NEXT(I)   = 0
               LAST(I)   = 0
             ELSE
               NEXT(LASTD) = I
               LAST(I)     = LASTD
               LASTD       = I
               NEXT(I)     = 0
             ENDIF
         GOTO 20
        ENDIF
C
C
        IF (DEG .GT. 0) THEN
           WF(I) = DEG
           IF (DEG.GT.N) THEN
            DEG = min(((DEG-N)/PAS) + N , NBBUCK)
           ENDIF
C           Note that if deg=0 then 
C           No fill-in will occur, 
C           but one variable is adjacent to I
C          ----------------------------------------------------------
C          place i in the degree list corresponding to its degree
C          ----------------------------------------------------------
           INEXT = HEAD (DEG)
           IF (INEXT .NE. 0) LAST (INEXT) = I
           NEXT (I) = INEXT
           HEAD (DEG) = I
        ELSE
C         ----------------------------------------------------------
C         we have a variable that can be eliminated at once because
C         there is no off-diagonal non-zero in its row.
C         ----------------------------------------------------------
           NEL = NEL + NV(I)
           ELEN (I) = -NEL
           PE (I) = 0_8
           W (I) = 0
        ENDIF
C=======================================================================
C
   20 CONTINUE
C=======================================================================
C  WHILE (selecting pivots) DO
C=======================================================================
      NLEFT = TOTEL-NEL
C=======================================================================
C =====================================================================
 30   IF (NEL .LT. TOTEL) THEN
C =====================================================================
C  GET PIVOT OF MINIMUM DEGREE
C=======================================================================
C       -------------------------------------------------------------
C       find next supervariable for elimination
C       -------------------------------------------------------------
         DO 40 DEG = MINDEG, NBBUCK
            ME = HEAD (DEG)
            IF (ME .GT. 0) GO TO 50
 40      CONTINUE
 50      MINDEG = DEG
         IF (ME.LE.0) THEN
            NCMPA = -N
            CALL MUMPS_ABORT()
         ENDIF
         IF (DEG.GT.N) THEN
C        -------------------------------
C        Linear search to find variable 
C        with best score in the list
C        -------------------------------
C        While end of list list not reached
C         NEXT(J) = 0
            J = NEXT(ME)
            K = WF(ME)
C     if ME is not available
            IF(CONSTRAINT(ME) .LT. 0) THEN
               K = -1
            ENDIF
 55         CONTINUE
            IF (J.GT.0) THEN
C     j is  available
               IF(CONSTRAINT(J) .GE. 0) THEN
                  IF (WF(J).LT.K .OR. K .LT. 0) THEN
                     ME = J
                     K  = WF(ME)
                  ENDIF
               ENDIF
               J= NEXT(J)
               GOTO 55
            ENDIF
            ILAST = LAST(ME)
            INEXT = NEXT(ME)
            IF (INEXT .NE. 0) LAST (INEXT) = ILAST
            IF (ILAST .NE. 0) THEN
               NEXT (ILAST) = INEXT
            ELSE
C          me is at the head of the degree list
               HEAD (DEG) = INEXT
            ENDIF
C     
         ELSE
C     select ME which verify the constraint           
C     if it is directly ok
            IF(CONSTRAINT(ME) .GE. 0) GOTO 59
 56         CONTINUE
C     if ME has a successor exaine it
            IF(NEXT(ME) .NE. 0) THEN
               ME = NEXT(ME)
               IF(CONSTRAINT(ME) .GE. 0) THEN
                  GOTO 59
               ELSE
                  GOTO 56
               ENDIF
            ELSE
C     ME has no successor -> increase deg till finding a valid ME
C     57: increase deg till a non empty list is found
 57            DEG = DEG+1
               ME = HEAD(DEG) 
C     no empty found 
               IF(ME .GT. 0) THEN
C     good piv found
                  IF(CONSTRAINT(ME) .GE. 0) THEN
                     GOTO 59
                  ELSE
C     else loop on next
                     GOTO 56
                  ENDIF
               ELSE
C     increase degree
                  GOTO 57
               ENDIF
            ENDIF   
 59         PREV = LAST (ME)
            INEXT = NEXT (ME)
            IF(PREV .NE. 0) THEN
               NEXT(PREV) = INEXT
            ELSE
               HEAD (DEG) = INEXT
            ENDIF
C     remove ME from the x2 linked lists
            IF (INEXT .NE. 0) LAST (INEXT) = PREV
         ENDIF
C         -------------------------------------------------------------
C         remove chosen variable from link list
C         -------------------------------------------------------------
         TOTO = ME
 5910    IF(TOTO .NE. 0) THEN
            J = CONSTRAINT(TOTO)
            IF(J .GT. 0) THEN
               CONSTRAINT(J) = 0
            ENDIF
            TOTO = THESON(TOTO)
            GOTO 5910
         ENDIF
C       -------------------------------------------------------------
C       me represents the elimination of pivots nel+1 to nel+nv(me).
C       place me itself as the first in this set.  It will be moved
C       to the nel+nv(me) position when the permutation vectors are
C       computed.
C       -------------------------------------------------------------
            ELENME = ELEN (ME)
            ELEN (ME) = - (NEL + 1)
            NVPIV = NV (ME)
            NEL = NEL + NVPIV
C=======================================================================
C  CONSTRUCT NEW ELEMENT
C=======================================================================
C       -------------------------------------------------------------
C       At this point, me is the pivotal supervariable.  It will be
C       converted into the current element.  Scan list of the
C       pivotal supervariable, me, setting tree pointers and
C       constructing new list of supervariables for the new element,
C       me.  p is a pointer to the current position in the old list.
C       -------------------------------------------------------------
C       flag the variable "me" as being in Lme by negating nv (me)
            NV (ME) = -NVPIV
            DEGME = 0
            IF (ELENME .EQ. 0) THEN
C         ----------------------------------------------------------
C         construct the new element in place
C         ----------------------------------------------------------
               PME1 = PE (ME)
               PME2 = PME1 - 1
               DO 60 P = PME1, PME1 + LEN (ME) - 1
                  I = IW (P)
                  NVI = NV (I)
                  IF (NVI .GT. 0) THEN
C             ----------------------------------------------------
C             i is a principal variable not yet placed in Lme.
C             store i in new list
C             ----------------------------------------------------
                     DEGME = DEGME + NVI
C             flag i as being in Lme by negating nv (i)
                     NV (I) = -NVI
                     PME2 = PME2 + 1
                     IW (PME2) = I
                     IF (DEGREE(I).NE.N2) THEN
C             ----------------------------------------------------
C             remove variable i from degree list. (only if i \in V0)
C             ----------------------------------------------------
                        ILAST = LAST (I)
                        INEXT = NEXT (I)
                        IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                        IF (ILAST .NE. 0) THEN
                           NEXT (ILAST) = INEXT
                        ELSE
C               i is at the head of the degree list
                           IF (WF(I).GT.N) THEN
                              DEG = min(((WF(I)-N)/PAS) + N , NBBUCK)
                           ELSE
                              DEG = WF(I)
                           ENDIF
                           HEAD (DEG) = INEXT
                        ENDIF
                     ENDIF
                  ENDIF
 60            CONTINUE
C         this element takes no new memory in iw:
               NEWMEM = 0
            ELSE
C         ----------------------------------------------------------
C         construct the new element in empty space, iw (pfree ...)
C         ----------------------------------------------------------
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          KNT1_UPDATED = 0
          DO 120 KNT1 = 1, ELENME + 1
            KNT1_UPDATED = KNT1_UPDATED +1
            IF (KNT1 .GT. ELENME) THEN
C             search the supervariables in me.
              E = ME
              PJ = P
              LN = SLENME
            ELSE
C             search the elements in me.
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
C           -------------------------------------------------------
C           search for different supervariables and add them to the
C           new list, compressing when necessary. this loop is
C           executed once for each element in the list and once for
C           all the supervariables in the list.
C           -------------------------------------------------------
            KNT2_UPDATED = 0
            DO 110 KNT2 = 1, LN
              KNT2_UPDATED = KNT2_UPDATED+1
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
C               -------------------------------------------------
C               compress iw, if necessary
C               -------------------------------------------------
                IF (PFREE .GT. IWLEN) THEN
C                 prepare for compressing iw by adjusting
C                 pointers and lengths so that the lists being
C                 searched in the inner and outer loops contain
C                 only the remaining entries.
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1_UPDATED
C                 Reset KNT1_UPDATED in case of recompress 
C                 at same iteration of the loop 120
                  KNT1_UPDATED = 0
C                 Check if anything left in supervariable ME
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0_8
                  PE (E) = PJ
                  LEN (E) = LN - KNT2_UPDATED
C                 Reset KNT2_UPDATED in case of recompress 
C                 at same iteration of the loop 110
                  KNT2_UPDATED = 0
C                 Check if anything left in element E
                  IF (LEN (E) .EQ. 0) PE (E) = 0_8
                  NCMPA = NCMPA + 1
C                 store first item in pe
C                 set first entry to -item
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = int(IW (PN),8)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
C                 psrc/pdst point to source/destination
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
C                 while loop:
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
C                   search for next negative entry
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = int(PE (J))
                      PE (J) = PDST
                      PDST = PDST + 1_8
C                     copy from source to destination
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + int(LENJ - 1,8)
                      PSRC = PSRC + int(LENJ - 1,8)
                    ENDIF
                    GO TO 80
                  ENDIF
C                 move the new partially-constructed element
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
C               -------------------------------------------------
C               i is a principal variable not yet placed in Lme
C               store i in new list
C               -------------------------------------------------
                DEGME = DEGME + NVI
C               flag i as being in Lme by negating nv (i)
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
              IF (DEGREE(I).NE.N2) THEN
C               -------------------------------------------------
C               remove variable i from degree link list 
C                            (only if i in V0)
C               -------------------------------------------------
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
                  IF (WF(I).GT.N) THEN
                   DEG = min(((WF(I)-N)/PAS) + N , NBBUCK)
                  ELSE
                   DEG = WF(I)
                  ENDIF
C                 i is at the head of the degree list
                  HEAD (DEG) = INEXT
                ENDIF
              ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
C             set tree pointer and flag to indicate element e is
C             absorbed into new element me (the parent of e is me)
              PE (E) = int(-ME,8)
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
C         this element takes newmem new memory in iw (possibly zero)
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = max (MAXMEM, MEM)
        ENDIF
C       -------------------------------------------------------------
C       me has now been converted into an element in iw (pme1..pme2)
C       -------------------------------------------------------------
C       degme holds the external degree of new element
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = int(PME2 - PME1 + 1_8)
C       -------------------------------------------------------------
C       make sure that wflg is not too large.  With the current
C       value of wflg, wflg+n must not cause integer overflow
C       -------------------------------------------------------------
        IF (WFLG .GT. MAXINT_N) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
C=======================================================================
C  COMPUTE (w (e) - wflg) = |Le\Lme| FOR ALL ELEMENTS
C=======================================================================
C       -------------------------------------------------------------
C       Scan 1:  compute the external degrees of previous elements
C       with respect to the current element.  That is:
C            (w (e) - wflg) = |Le \ Lme|
C       for each element e that appears in any supervariable in Lme.
C       The notation Le refers to the pattern (list of
C       supervariables) of a previous element e, where e is not yet
C       absorbed, stored in iw (pe (e) + 1 ... pe (e) + iw (pe (e))).
C       The notation Lme refers to the pattern of the current element
C       (stored in iw (pme1..pme2)).   If (w (e) - wflg) becomes
C       zero, then the element e will be absorbed in scan 2.
C       -------------------------------------------------------------
        DO 150 PME = PME1, PME2
          I = IW (PME)
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
C           note that nv (i) has been negated to denote i in Lme:
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + int(ELN - 1,8)
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
C               unabsorbed element e has been seen in this loop
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
C               e is an unabsorbed element
C               this is the first we have seen e in all of Scan 1
                WE = DEGREE (E) + WNVI
                WF(E) = 0
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
C=======================================================================
C  DEGREE UPDATE AND ELEMENT ABSORPTION
C=======================================================================
C       -------------------------------------------------------------
C       Scan 2:  for each i in Lme, sum up the degree of Lme (which
C       is degme), plus the sum of the external degrees of each Le
C       for the elements e appearing within i, plus the
C       supervariables in i.  Place i in hash list.
C       -------------------------------------------------------------
        DO 180 PME = PME1, PME2
          I = IW (PME)
          P1 = PE (I)
          P2 = P1 + int(ELEN (I) - 1,8)
          PN = P1
          HASH = 0_8
          DEG  = 0
          WF3  = 0
          WF4  = 0
          NVI  = -NV(I)
C         ----------------------------------------------------------
C         scan the element list associated with supervariable i
C         ----------------------------------------------------------
          DO 160 P = P1, P2
            E = IW (P)
C           dext = | Le \ Lme |
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              IF ( WF(E) .EQ. 0 ) THEN
C              First time we meet e : compute wf(e) 
C              which holds the surface associated to element e 
C              it will later be deducted from fill-in 
C              area of all variables adjacent to e
               WF(E) = DEXT * ( (2 * DEGREE(E))  -  DEXT - 1)
              ENDIF
              WF4 = WF4 + WF(E)
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
            ELSE IF (DEXT .EQ. 0) THEN
#if defined (NOAGG4)
              IW (PN) = E
              PN = PN + 1
              HASH = HASH + int(E,kind=8)
#else
C             aggressive absorption: e is not adjacent to me, but
C             the |Le \ Lme| is 0, so absorb it into me
              PE (E) = int(-ME,8)
              W (E) = 0
#endif
            ENDIF
  160     CONTINUE
C         count the number of elements in i (including me):
          ELEN (I) = int(PN - P1 + 1_8)
C         ----------------------------------------------------------
C         scan the supervariables in the list associated with i
C         ----------------------------------------------------------
          P3 = PN
          DO 170 P = P2 + 1, P1 + int(LEN (I) - 1,8)
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
C             j is unabsorbed, and not in Lme.
C             add to degree and add to new list
              DEG = DEG + NVJ
              WF3 = WF3 + NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + int(J,kind=8)
            ENDIF
  170     CONTINUE
C
          IF (DEGREE(I).EQ.N2) DEG = N2
C         ----------------------------------------------------------
C         update the degree and check for mass elimination
C         ----------------------------------------------------------
#if defined (NOAGG4)
          IF (ELEN(I).EQ.1 .AND. P3.EQ.PN) THEN
#else
          IF (DEG .EQ. 0) THEN
#endif
C           -------------------------------------------------------
C           mass elimination
C           -------------------------------------------------------
C           There is nothing left of this node except for an
C           edge to the current pivot element.  elen (i) is 1,
C           and there are no variables adjacent to node i.
C           Absorb i into the current pivot element, me.
             TOTO = I
 5911        IF(TOTO .NE. 0) THEN
                J = CONSTRAINT(TOTO)
                IF(J .GT. 0) THEN
                   CONSTRAINT(J) = 0
                ENDIF
                TOTO = THESON(TOTO)
                GOTO 5911
             ENDIF
            PE (I) = int(-ME,8)
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
C           -------------------------------------------------------
C           update the upper-bound degree of i
C           -------------------------------------------------------
C           the following degree does not yet include the size
C           of the current element, which is added later:
C AMD            DEGREE (I) = min (DEGREE (I), DEG)
            IF (DEGREE(I).NE.N2) THEN
C                I does not belong to halo
C                dk = min (d(k-1)+degme, deg+degme)
                 IF ( DEGREE (I).LT.DEG ) THEN
C                  Our appox degree is loose.
C                  we keep old value. Note that in 
C                  this case we cannot substract WF(I)
C                  for min-fill score.
                   WF4 = 0
                   WF3 = 0
                 ELSE
                   DEGREE(I)  = DEG
                 ENDIF
            ENDIF
C
C           compute WF(I) taking into account size of block 3.0
            WF(I)      = WF4 + 2*NVI*WF3
C           -------------------------------------------------------
C           add me to the list for i
C           -------------------------------------------------------
C           move first supervariable to end of list
            IW (PN) = IW (P3)
C           move first element to end of element part of list
            IW (P3) = IW (P1)
C           add new element to front of list.
            IW (P1) = ME
C           store the new length of the list in len (i)
            LEN (I) = int(PN - P1 + 1_8)
            IF (DEG.NE.N2) THEN
C           -------------------------------------------------------
C           place in hash bucket.  Save hash key of i in last (i).
C           -------------------------------------------------------
            HASH = mod (HASH, HMOD) + 1_8
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
C             the degree list is empty, hash head is -j
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
C             degree list is not empty
C             use last (head (hash)) as hash head
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = int(HASH,kind=kind(LAST))
            ENDIF
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
C       -------------------------------------------------------------
C       Clear the counter array, w (...), by incrementing wflg.
C       -------------------------------------------------------------
        DMAX = max (DMAX, DEGME)
        WFLG = WFLG + DMAX
C       make sure that wflg+n does not cause integer overflow
        IF (WFLG .GT. MAXINT_N) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
C       at this point, w (1..n) .lt. wflg holds
C=======================================================================
C  SUPERVARIABLE DETECTION
C=======================================================================
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF ( (NV (I) .LT. 0) .AND. (DEGREE(I).NE.N2) ) THEN
C           i is a principal variable in Lme
C           -------------------------------------------------------
C           examine all hash buckets with 2 or more variables.  We
C           do this by examing all unique hash keys for super-
C           variables in the pattern Lme of the current element, me
C           -------------------------------------------------------
            HASH = int(LAST (I),kind=8)
C           let i = head of hash bucket, and empty the hash bucket
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
C             degree list is empty
              I = -J
              HEAD (HASH) = 0
            ELSE
C             degree list is not empty, restore last () of head
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
C           while loop:
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
C             ----------------------------------------------------
C             this bucket has one or more variables following i.
C             scan all of them to see if i can absorb any entries
C             that follow i in hash bucket.  Scatter i into w.
C             ----------------------------------------------------
              LN = LEN (I)
              ELN = ELEN (I)
C             do not flag the first element in the list (me)
              DO 210 P = PE (I) + 1, PE (I) + int(LN - 1,8)
                W (IW (P)) = WFLG
  210         CONTINUE
C             ----------------------------------------------------
C             scan every other entry j following i in bucket
C             ----------------------------------------------------
              JLAST = I
              J = NEXT (I)
C             while loop:
  220         CONTINUE
              IF (J .NE. 0) THEN
                 IF(CONSTRAINT(J) .LT. 0
     &                .AND. CONSTRAINT(I) .LT. 0) THEN
                    GOTO 240
                 ENDIF
                 IF(CONSTRAINT(I) .GE. 0) THEN
                    IF(CONSTRAINT(J) .LT. 0) THEN
                       TOTO = I
 221                   IF(TOTO .NE. 0) THEN
                          IF(CONSTRAINT(TOTO) .EQ. J) THEN
                             GOTO 225
                          ENDIF
                          TOTO =THESON(TOTO)
                          GOTO 221
                       ENDIF
                    ELSE
                       GOTO 225
                    ENDIF
                 ELSE
C     if I is locked see if it is freed thanks to J
                    IF(CONSTRAINT(J) .GE. 0) THEN
                       TOTO = J
 222                   IF(TOTO .NE. 0) THEN
                          IF(CONSTRAINT(TOTO) .EQ. I) THEN
                             GOTO 225
                          ENDIF
                          TOTO =THESON(TOTO)
                          GOTO 222
                       ENDIF
                    ENDIF
                 ENDIF
                 GOTO 240
 225             CONTINUE
C               -------------------------------------------------
C               check if j and i have identical nonzero pattern
C               -------------------------------------------------
C               jump if i and j do not have same size data structure
                 IF (LEN (J) .NE. LN) GO TO 240
C               jump if i and j do not have same number adj elts
                 IF (ELEN (J) .NE. ELN) GO TO 240
C               do not flag the first element in the list (me)
                 DO 230 P = PE (J) + 1_8, PE (J) + int(LN - 1,8)
C                 jump if an entry (iw(p)) is in j but not in i
                    IF (W (IW (P)) .NE. WFLG) GO TO 240
 230             CONTINUE
C               -------------------------------------------------
C               found it!  j can be absorbed into i
C               -------------------------------------------------
C     update the supervariable composition
                 TOTO = I
 231             IF(THESON(TOTO) .NE. 0) THEN
                    TOTO = THESON(TOTO)
                    GOTO 231
                 ENDIF
                 THESON(TOTO) = J
                 IF(CONSTRAINT(I) .LT. 0) THEN
                    CONSTRAINT(I) = 0
                 ENDIF
                 PE (J) = int(-I,8)
                 WF(I)  = max(WF(I),WF(J))
C               both nv (i) and nv (j) are negated since they
C               are in Lme, and the absolute values of each
C               are the number of variables in i and j:
                 NV (I) = NV (I) + NV (J)
                 NV (J) = 0
                 ELEN (J) = 0
C               delete j from hash bucket
                 J = NEXT (J)
                 NEXT (JLAST) = J
                 GO TO 220
C               -------------------------------------------------
 240             CONTINUE
C               j cannot be absorbed into i
C               -------------------------------------------------
                 JLAST = J
                 J = NEXT (J)
                 GO TO 220
              ENDIF
C             ----------------------------------------------------
C             no more variables can be absorbed into i
C             go to next i in bucket and clear flag array
C             ----------------------------------------------------
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
           ENDIF
          ENDIF
 250   CONTINUE
C=======================================================================
C  RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVAR. FROM ELEMENT
C=======================================================================
        P = PME1
        NLEFT = TOTEL - NEL
        DO 260 PME = PME1, PME2
           I = IW (PME)
           NVI = -NV (I)
           IF (NVI .GT. 0) THEN
C           i is a principal variable in Lme
C           restore nv (i) to signify that i is principal
              NV (I) = NVI
              IF (DEGREE(I).NE.N2) THEN
C           -------------------------------------------------------
C           compute the external degree (add size of current elem)
C           -------------------------------------------------------
                 DEG = min (DEGREE (I) + DEGME - NVI, NLEFT - NVI)
                 IF (DEGREE (I) + DEGME .GT. NLEFT ) THEN
C
                  DEG = DEGREE(I)
                  RMF1  = dble(DEG)*dble( (DEG-1) + 2*DEGME )
     &                 - dble(WF(I))
                  DEGREE(I) = NLEFT - NVI
                  DEG       = DEGREE(I) 
                  RMF = dble(DEG)*dble(DEG-1) 
     &                 -  dble(DEGME-NVI)*dble(DEGME-NVI-1)
                  RMF = min(RMF, RMF1)
               ELSE 
                  DEG = DEGREE(I)
                  DEGREE(I) = DEGREE (I) + DEGME - NVI
                  RMF  = dble(DEG)*dble( (DEG-1) + 2*DEGME ) 
     &                 - dble(WF(I))
               ENDIF
               RMF =  RMF / dble(NVI+1)
C
               IF (RMF.LT.dummy) THEN
                  WF(I) = int ( anint( RMF ))
               ELSEIF (RMF / dble(N) .LT. dummy) THEN 
                  WF(I) = int ( anint( RMF/dble(N) ))
               ELSE
                  WF(I) = idummy
               ENDIF
               WF(I) = max(1,WF(I))
C           -------------------------------------------------------
C           place the supervariable at the head of the degree list
C           -------------------------------------------------------
               DEG = WF(I)
               IF (DEG.GT.N) THEN
                  DEG = min(((DEG-N)/PAS) + N , NBBUCK)
               ENDIF
               INEXT = HEAD (DEG)
               IF (INEXT .NE. 0) LAST (INEXT) = I
               NEXT (I) = INEXT
               LAST (I) = 0
               HEAD (DEG) = I
C           -------------------------------------------------------
C           save the new degree, and find the minimum degree
C           -------------------------------------------------------
               MINDEG = min (MINDEG, DEG)
            ENDIF
C           -------------------------------------------------------
C           place the supervariable in the element pattern
C           -------------------------------------------------------
            IW (P) = I
            P = P + 1
         ENDIF
 260  CONTINUE
C=======================================================================
C  FINALIZE THE NEW ELEMENT
C=======================================================================
      NV (ME) = NVPIV + DEGME
C       fill_est = fill_est + nvpiv * (nvpiv + 2 * degme)
C       nv (me) is now the degree of pivot (including diagonal part)
C       save the length of the list for the new element me
      LEN (ME) = int(P - PME1)
      IF (LEN (ME) .EQ. 0) THEN
C         there is nothing left of the current pivot element
         PE (ME) = 0_8
         W (ME) = 0
      ENDIF
      IF (NEWMEM .NE. 0) THEN
C         element was not constructed in place: deallocate part
C         of it (final size is less than or equal to newmem,
C         since newly nonprincipal variables have been removed).
         PFREE = P
         MEM = MEM - NEWMEM + int(LEN (ME),8)
      ENDIF
C=======================================================================
C       END WHILE (selecting pivots)
      GO TO 30
      ENDIF
C=======================================================================
C begin HALO V2
      IF (NBFLAG.GT.0) THEN
C
C     All possible pivots (not flagged have been eliminated).
C     We amalgamate all flagged variables at the root and 
C     we finish the elimination tree.
C          1/ Go through all
C          non absorbed elements (root of the subgraph)
C          and absorb in ME
C          2/ perform mass elimination of all dense rows
           DO DEG = MINDEG, NBBUCK+1
             ME = HEAD (DEG)
             IF (ME .GT. 0) GO TO 51
           ENDDO
   51      MINDEG = DEG
           NELME    = -(NEL+1)
           DO X=1,N
            IF ((PE(X).GT.0_8) .AND. (ELEN(X).LT.0)) THEN
C            X is an unabsorbed element
             PE(X) = int(-ME,8)
C            W(X) = 0 could be suppressed ?? check it
            ELSEIF (DEGREE(X).EQ.N2) THEN
C            X is a dense row, absorb it in ME (mass elimination)
             NEL   = NEL + NV(X)
             PE(X) = int(-ME,8)
             ELEN(X) = 0
C            Correct value of NV is (secondary variable)
             NV(X) = 0
            ENDIF
           ENDDO
C          ME is the root node
           ELEN(ME) = NELME
C          Correct value of NV is (principal variable)
           NV(ME)   = N-NREAL
           PE(ME)   = 0_8
C
      ENDIF
C end HALO
C=======================================================================
C  COMPUTE THE PERMUTATION VECTORS
C=======================================================================
C     ----------------------------------------------------------------
C     The time taken by the following code is O(n).  At this
C     point, elen (e) = -k has been done for all elements e,
C     and elen (i) = 0 has been done for all nonprincipal
C     variables i.  At this point, there are no principal
C     supervariables left, and all elements are absorbed.
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
C     compute the ordering of unordered nonprincipal variables
C     ----------------------------------------------------------------
      DO 290 I = 1, N
         IF (ELEN (I) .EQ. 0) THEN
C         ----------------------------------------------------------
C         i is an un-ordered row.  Traverse the tree from i until
C         reaching an element, e.  The element, e, was the
C         principal supervariable of i and all nodes in the path
C         from i to when e was selected as pivot.
C         ----------------------------------------------------------
            J = int(-PE (I))
C         while (j is a variable) do:
 270        CONTINUE
            IF (ELEN (J) .GE. 0) THEN
               J = int(-PE (J))
               GO TO 270
            ENDIF
            E = J
C           ----------------------------------------------------------
C           get the current pivot ordering of e
C           ----------------------------------------------------------
            K = -ELEN (E)
C           ----------------------------------------------------------
C           traverse the path again from i to e, and compress the
C           path (all nodes point to e).  Path compression allows
C           this code to compute in O(n) time.  Order the unordered
C           nodes in the path, and place the element e at the end.
C           ----------------------------------------------------------
            J = I
C           while (j is a variable) do:
 280        CONTINUE
            IF (ELEN (J) .GE. 0) THEN
               JNEXT = int(-PE (J))
               PE (J) = int(-E,8)
               IF (ELEN (J) .EQ. 0) THEN
C               j is an unordered row
                  ELEN (J) = K
                  K = K + 1
               ENDIF
               J = JNEXT
               GO TO 280
            ENDIF
C         leave elen (e) negative, so we know it is an element
            ELEN (E) = -K
         ENDIF
 290  CONTINUE
C     ----------------------------------------------------------------
C     reset the inverse permutation (elen (1..n)) to be positive,
C     and compute the permutation (last (1..n)).
C     ----------------------------------------------------------------
      IF(.TRUE.) THEN
C       N is the size of the compressed graph.
C       If the graph was compressed on input then
C       indices in ELEN are in [1,TOTEL]
C       We build the inverse of ELEN in LAST (similar to
C       the pivot order but has zeros in it) and then compress
C       it. Since LAST is assumed to be of size N at the
C       interface level, we need another array to store
C       the inverse of ELEN for entries greater than N
C       We use DEGREE.
        LAST(1:N) = 0
        DEGREE(1:TOTEL-N)=0
        DO I = 1, N
          K = abs (ELEN (I))
          IF ( K <= N ) THEN
            LAST (K) = I
          ELSE
            DEGREE(K-N)=I
          ENDIF
        ENDDO
        I = 1
        DO K = 1, N
          IF(LAST (K) .NE. 0) THEN
            LAST(I) = LAST(K)
            ELEN(LAST(K)) = I
            I = I + 1
          ENDIF
        ENDDO
        DO K = N+1, TOTEL
          IF (DEGREE(K-N) .NE. 0) THEN
            LAST(I)=DEGREE(K-N)
            ELEN(DEGREE(K-N)) = I
            I = I + 1
          ENDIF
        END DO
      ELSE
        DO 300 I = 1, N
           K = abs (ELEN (I))
           LAST (K) = I
           ELEN (I) = K
300     CONTINUE
      ENDIF
C=======================================================================
C  RETURN THE MEMORY USAGE IN IW
C=======================================================================
C     If maxmem is less than or equal to iwlen, then no compressions
C     occurred, and iw (maxmem+1 ... iwlen) was unused.  Otherwise
C     compressions did occur, and iwlen would have had to have been
C     greater than or equal to maxmem for no compressions to occur.
C     Return the value of maxmem in the pfree argument.
      PFREE = MAXMEM
C===============================
C     Save PE in PARENT array
      DO I=1,N
       PARENT(I) = int(PE(I))
      ENDDO
C===============================
      RETURN
      END SUBROUTINE MUMPS_CST_AMF
C-----------------------------------------------------------------------
C MUMPS_SYMQAMD: modified version of MUMPS_QAMD code to 
C designed to compute a symbolic factorization given 
C an input ordering (provided in PERM array) and possibly
C a schur area.
C ---------
      SUBROUTINE MUMPS_SYMQAMD
     &                ( THRESH, NDENSE, 
     &                 N, TOTEL, IWLEN, PE, PFREE, LEN, IW, NV, 
     &                 ELEN, LAST, NCMPA, DEGREE, HEAD, NEXT, W, 
     &                 PERM, LISTVAR_SCHUR, SIZE_SCHUR, 
     &                 AGG6, PARENT ) 
      IMPLICIT NONE
C     Input not modified
      INTEGER, INTENT(IN)    :: N, TOTEL, SIZE_SCHUR
      LOGICAL, INTENT(IN)    :: AGG6
      INTEGER, INTENT(IN)    :: THRESH
      INTEGER(8), INTENT(IN) :: IWLEN
      INTEGER, INTENT(IN)    :: LISTVAR_SCHUR(max(1,SIZE_SCHUR))
C     Input undefined on output 
      INTEGER, INTENT(INOUT)  :: LEN(N), IW(IWLEN)
C 
C     Output only 
      INTEGER, INTENT(OUT)   :: NCMPA
      INTEGER, INTENT(OUT)   :: ELEN(N), LAST(TOTEL), PARENT(N)
C 
C     Input/output
      INTEGER, INTENT(INOUT)    :: NV(N)
      INTEGER(8), INTENT(INOUT) :: PFREE
      INTEGER(8), INTENT(INOUT) :: PE(N)
      INTEGER, INTENT(INOUT)    :: PERM(N)
C 
C     Internal Workspace only
      INTEGER, INTENT(OUT) :: NDENSE(N), DEGREE(N), 
     &                        HEAD(TOTEL), NEXT(N), W(N)
C
C  =======================
C  INTERFACE DOCUMENTATION
C  SPECIFIC TO SYMQAMD.
C  =======================
C  (more details are sometimes
C   available in the
C   PREVIOUS DOCUMENTATION
C   section)
C
C N (in): the size of the matrix 
C         number of supervariables if blocked format
C TOTEL (in) : Number of variables to eliminate
C
C IWLEN (in): the length of the workspace IW
C
C PFREE (inout): says that IW(1:PFREE-1) contains the graph on input, see
C             below. (on output see meaning bellow)
C IW (inout):
C    On input, IW(1:PFREE-1) contains the orginal graph
C    On output it has been corrupted because IW(1:IWLEN) has been
C    used as workspace.
C
C LEN(inout):  On input, 
C       LEN (i) holds the number of entries in row i of the
C       matrix, excluding the diagonal.  The contents of LEN(1..N)
C       are undefined on output.
C
C PE(inout): On input PE(i) contains the pointers in IW to (the column
C       indices of) row i of the matrix.
C       On output it contains the tree:
C       - if I is a principal variable (NV(I) >0) then -pe(I) is the principal
C         variable of the father, or 0 if I is a root node.
C       - if I is a secondary variable (NV(I)=0) then -pe(I) is the principal
C         variable of the node it belongs to.
C
C       On output:  (PE is copied on output into PARENT array)
C
C  
C NV(inout): 
C          On input:  encoding of a blocked matrix 
C            if NV(1).NE.-1 the NV(I) holds the weight of node I. 
C          During execution, 
C            abs (nv (i)) is equal to the number of rows
C            that are represented by the principal supervariable i.  
C            If i is a nonprincipal variable, then nv (i) = 0.  
C            nv (i) .lt. 0 signifies that i is a
C            principal variable in the pattern Lme of the current pivot
C            element me. 
C          On output: 
C          - if i is a principal variable, NV(i) is the size of the front
C          in the multifrontal terminology.
C          - if i is a secondary variable, NV(i)=0
C
C PERM  (inout) : MUST BE SET TO HOLD THE POSITION OF VARIABLE I IN THE
C     PERMUTED ORDER.
C     PERM(I) = J means that I is the Jth pivot.
C     PERM IS NOT ALTERED IF SIZE_SCHUR = 0.
C     IF SIZE_SCHUR > 0 and variable I is part of the Schur,
C     then PERM(I) must be greater than N - SIZE_SCHUR.
C     In that case, PERM(I) is altered: it is set to N+1 internally !
C
C SIZE_SCHUR (in) :   > 0 means that the last SIZE_SCHUR variable 
C                in the order (such that PERM(I) >  N-SIZE_SCHUR) 
C                are part of the schur decompositon
C                and should remain ordered last and amalgamated
C                at the root of the elimination tree.
C
C LISTVAR_SCHUR(1:SIZE_SCHUR) (in): should be set on entry to the list of
C                variables (original indices) in the Schur complement
C
C THRESH (in): is used to set the local variable THRESM, corresponding
C      to the internal restarting feature.
C      <= 0 Recommended value. Automatic setting will be done.
C          Note that this does not correspond to the historical
C          documentation further below.
C       =  N Only exactly dense rows in the reduced matrix are selected.
C       >  1 and <= N THRESH correspond to the minimum density requirement.
C
C      At the moment if SIZE_SCHUR > 0 restarting functionality is disabled,
C      which means that performance is not optimal. It should work again with
C      a small modification but this has to be tested when it is re-enabled.
C
C ELEN (out) needs not be set on entry.
C      It contains the inverse
C      permutation on output. Not sure what it contains for the Schur
C      variables.
C      (it should be ok for the Schur too).
C
C LAST  used internally as working space; 
C       On output, last (1..n) holds the permutation,  i = last (k), then
C       row i is the kth pivot row.
C       Not used on output and 
C       Computation has been suppressed 
C       since in the context of blocked matrix format
C       one cannot so easily compute last out of elen 
C       (see end of MUMPS_QAMD in case of COMRPESS, 
C       because elen(i) \in [1:TOTEL] and not \in [1:N])
C
C AGG6 (in): controls if aggressive absorption should be authorized.
C
C  -------------------------------------------
C  ARGUMENTS USED INTERNALLY AS WORKARRAYS
C  Maybe some things are significant on output
C  but not in the normal cases of usage.
C  -------------------------------------------
C
C  NDENSE, LAST, NEXT, HEAD, DEGREE, W
C
C  ------
C  OUTPUT
C  ------
C
C  NCMPA (out): number of compressions.
C
C
C  ======================
C  PREVIOUS DOCUMENTATION
C  ======================
C
C NDENSE of an element is the number of dense rows in the element.
C-----------------------------------------------------------------------
C It is a modified version of MUMPS_QAMD
C designed to automatically detect and exploit dense or quasi dense 
C rows in the reduced matrix at any step of the minimum degree. 
C The input integer parameter THRESH defines the quasi density:
C THRESH : input parameter (not modified) 
C  THRESH is used to compute THRESM 
C   <=0 or N Only exactly dense rows in the reduced matrix are selected.
C   >1 and <=N THRESH correspond to the munimum density requirement.
C      Version 0: All dense and quasi dense rows are amalgamated at the 
C                 root node.
C      Version 1: Restart AMD with all quasi dense rows, and 
C                 increase density requirement.
C-----------------------------------------------------------------------
C Additionnal parameters/variables due to dense row manipulation:
C          
C Local variables:
C ---------------
      INTEGER THRESM, NDME, PERMeqN
      INTEGER NBD,NBED, NBDM, LASTD, NELME
      LOGICAL IDENSE
C THRESM : Local Integer holding a 
C          potentially modified value of THRESH.
C          When quasi dense rows are reintegrated in the 
C          graph to be processed then THRESM is modified.
C   Note that if one sets THRESM to negative value then
C       <0 Classical AMD algorithm (no dense row detection)
C NDME  : number of dense row adjacent to me
C NELME number of pivots selected when reching the root
C LASTD index of the last row in the list of dense rows
C NBD is the total number of dense rows selected 
C NBED is the total number of exactly dense rows detected. 
C NBDM is the maximum number of dense rows selected 
C IDENSE is used to indicate that the supervariable I is a dense or
C        quasi-dense row.
C-----------------------------------------------------------------------
C INPUT ARGUMENTS (unaltered):
C-----------------------------------------------------------------------
C n:    The matrix order.
C
C       Restriction:  n .ge. 1
C iwlen:        The length of iw (1..iwlen).  On input, the matrix is
C       stored in iw (1..pfree-1).  However, iw (1..iwlen) should be
C       slightly larger than what is required to hold the matrix, at
C       least iwlen .ge. pfree + n is recommended.  Otherwise,
C       excessive compressions will take place.
C       *** We do not recommend running this algorithm with ***
C       ***      iwlen .lt. pfree + n.                      ***
C       *** Better performance will be obtained if          ***
C       ***      iwlen .ge. pfree + n                       ***
C       *** or better yet                                   ***
C       ***      iwlen .gt. 1.2 * pfree                     ***
C       *** (where pfree is its value on input).            ***
C       The algorithm will not run at all if iwlen .lt. pfree-1.
C
C       Restriction: iwlen .ge. pfree-1
C-----------------------------------------------------------------------
C INPUT/OUPUT ARGUMENTS:
C-----------------------------------------------------------------------
C pe:   On input, pe (i) is the index in iw of the start of row i, or
C       zero if row i has no off-diagonal non-zeros.
C
C       During execution, it is used for both supervariables and
C       elements:
C
C       * Principal supervariable i:  index into iw of the
C               description of supervariable i.  A supervariable
C               represents one or more rows of the matrix
C               with identical nonzero pattern.
C       * Non-principal supervariable i:  if i has been absorbed
C               into another supervariable j, then pe (i) = -j.
C               That is, j has the same pattern as i.
C               Note that j might later be absorbed into another
C               supervariable j2, in which case pe (i) is still -j,
C               and pe (j) = -j2.
C       * Unabsorbed element e:  the index into iw of the description
C               of element e, if e has not yet been absorbed by a
C               subsequent element.  Element e is created when
C               the supervariable of the same name is selected as
C               the pivot.
C       * Absorbed element e:  if element e is absorbed into element
C               e2, then pe (e) = -e2.  This occurs when the pattern of
C               e (that is, Le) is found to be a subset of the pattern
C               of e2 (that is, Le2).  If element e is "null" (it has
C               no nonzeros outside its pivot block), then pe (e) = 0.
C
C       On output, pe holds the assembly tree/forest, which implicitly
C       represents a pivot order with identical fill-in as the actual
C       order (via a depth-first search of the tree).
C
C       On output:
C       If nv (i) .gt. 0, then i represents a node in the assembly tree,
C       and the parent of i is -pe (i), or zero if i is a root.
C       If nv (i) = 0, then (i,-pe (i)) represents an edge in a
C       subtree, the root of which is a node in the assembly tree.
C pfree:        On input, the matrix is stored in iw (1..pfree-1) and
C       the rest of the array iw is free.
C       During execution, additional data is placed in iw, and pfree
C       is modified so that components  of iw from pfree are free.
C       On output, pfree is set equal to the size of iw that
C       would have been needed for no compressions to occur.  If
C       ncmpa is zero, then pfree (on output) is less than or equal to
C       iwlen, and the space iw (pfree+1 ... iwlen) was not used.
C       Otherwise, pfree (on output) is greater than iwlen, and all the
C       memory in iw was used.
C-----------------------------------------------------------------------
C INPUT/MODIFIED (undefined on output):
C-----------------------------------------------------------------------
C len:  On input, len (i) holds the number of entries in row i of the
C       matrix, excluding the diagonal.  The contents of len (1..n)
C       are undefined on output.
C iw:   On input, iw (1..pfree-1) holds the description of each row i
C       in the matrix.  The matrix must be symmetric, and both upper
C       and lower triangular parts must be present.  The diagonal must
C       not be present.  Row i is held as follows:
C
C               len (i):  the length of the row i data structure
C               iw (pe (i) ... pe (i) + len (i) - 1):
C                       the list of column indices for nonzeros
C                       in row i (simple supervariables), excluding
C                       the diagonal.  All supervariables start with
C                       one row/column each (supervariable i is just
C                       row i).
C               if len (i) is zero on input, then pe (i) is ignored
C               on input.
C
C               Note that the rows need not be in any particular order,
C               and there may be empty space between the rows.
C
C       During execution, the supervariable i experiences fill-in.
C       This is represented by placing in i a list of the elements
C       that cause fill-in in supervariable i:
C
C               len (i):  the length of supervariable i
C               iw (pe (i) ... pe (i) + elen (i) - 1):
C                       the list of elements that contain i.  This list
C                       is kept short by removing absorbed elements.
C               iw (pe (i) + elen (i) ... pe (i) + len (i) - 1):
C                       the list of supervariables in i.  This list
C                       is kept short by removing nonprincipal
C                       variables, and any entry j that is also
C                       contained in at least one of the elements
C                       (j in Le) in the list for i (e in row i).
C
C       When supervariable i is selected as pivot, we create an
C       element e of the same name (e=i):
C
C               len (e):  the length of element e
C               iw (pe (e) ... pe (e) + len (e) - 1):
C                       the list of supervariables in element e.
C
C       An element represents the fill-in that occurs when supervariable
C       i is selected as pivot (which represents the selection of row i
C       and all non-principal variables whose principal variable is i).
C       We use the term Le to denote the set of all supervariables
C       in element e.  Absorbed supervariables and elements are pruned
C       from these lists when computationally convenient.
C
C       CAUTION:  THE INPUT MATRIX IS OVERWRITTEN DURING COMPUTATION.
C       The contents of iw are undefined on output.
C-----------------------------------------------------------------------
C OUTPUT (need not be set on input):
C-----------------------------------------------------------------------
C nv:   During execution, abs (nv (i)) is equal to the number of rows
C       that are represented by the principal supervariable i.  If i is
C       a nonprincipal variable, then nv (i) = 0.  Initially,
C       nv (i) = 1 for all i.  nv (i) .lt. 0 signifies that i is a
C       principal variable in the pattern Lme of the current pivot
C       element me.  On output, nv (e) holds the true degree of element
C       e at the time it was created (including the diagonal part).
C elen: See the description of iw above.  At the start of execution,
C       elen (i) is set to zero.  During execution, elen (i) is the
C       number of elements in the list for supervariable i.  When e
C       becomes an element, elen (e) = -nel is set, where nel is the
C       current step of factorization.  elen (i) = 0 is done when i
C       becomes nonprincipal.
C
C       For variables, elen (i) .ge. 0 holds until just before the
C       permutation vectors are computed.  For elements,
C       elen (e) .lt. 0 holds.
C
C       On output elen (1..n) holds the inverse permutation (the same
C       as the 'INVP' argument in Sparspak).  That is, if k = elen (i),
C       then row i is the kth pivot row.  Row i of A appears as the
C       (elen(i))-th row in the permuted matrix, PAP^T.
C last: In a degree list, last (i) is the supervariable preceding i,
C       or zero if i is the head of the list.  In a hash bucket,
C       last (i) is the hash key for i.  last (head (hash)) is also
C       used as the head of a hash bucket if head (hash) contains a
C       degree list (see head, below).
C
C       On output, last (1..n) holds the permutation (the same as the
C       'PERM' argument in Sparspak).  That is, if i = last (k), then
C       row i is the kth pivot row.  Row last (k) of A is the k-th row
C       in the permuted matrix, PAP^T.
C ncmpa:        The number of times iw was compressed.  If this is
C       excessive, then the execution took longer than what could have
C       been.  To reduce ncmpa, try increasing iwlen to be 10% or 20%
C       larger than the value of pfree on input (or at least
C       iwlen .ge. pfree + n).  The fastest performance will be
C       obtained when ncmpa is returned as zero.  If iwlen is set to
C       the value returned by pfree on *output*, then no compressions
C       will occur.
C-----------------------------------------------------------------------
C LOCAL (not input or output - used only during execution):
C-----------------------------------------------------------------------
C degree:       If i is a supervariable, then degree (i) holds the
C       current approximation of the external degree of row i (an upper
C       bound).  The external degree is the number of nonzeros in row i,
C       minus abs (nv (i)) (the diagonal part).  The bound is equal to
C       the external degree if elen (i) is less than or equal to two.
C
C       We also use the term "external degree" for elements e to refer
C       to |Le \ Lme|.  If e is an element, then degree (e) holds |Le|,
C       which is the degree of the off-diagonal part of the element e
C       (not including the diagonal part).
C degree (I) =N+1 if I is an exactly dense row in reduced matrix.
C            =N+1+LAST_approximate_external_deg of I   
C                      if I is a quasi dense row in reduced matrix.
C All dense or quasi dense rows are stored in the list pointed 
C       by head(n). Quasi-dense rows (degree(I)=n) are stored first, 
C       and are followed by exactly dense rows in the reduced matrix.
C       LASTD holds the last row in this list of dense rows or is zero
C       if the list is empty.
C head: head is used for degree lists.  head (deg) is the first
C       supervariable in a degree list (all supervariables i in a
C       degree list deg have the same approximate degree, namely,
C       deg = degree (i)).  If the list deg is empty then
C       head (deg) = 0.
C
C       During supervariable detection head (hash) also serves as a
C       pointer to a hash bucket.
C       If head (hash) .gt. 0, there is a degree list of degree hash.
C               The hash bucket head pointer is last (head (hash)).
C       If head (hash) = 0, then the degree list and hash bucket are
C               both empty.
C       If head (hash) .lt. 0, then the degree list is empty, and
C               -head (hash) is the head of the hash bucket.
C       After supervariable detection is complete, all hash buckets
C       are empty, and the (last (head (hash)) = 0) condition is
C       restored for the non-empty degree lists.
C next: next (i) is the supervariable following i in a link list, or
C       zero if i is the last in the list.  Used for two kinds of
C       lists:  degree lists and hash buckets (a supervariable can be
C       in only one kind of list at a time).
C w:    The flag array w determines the status of elements and
C       variables, and the external degree of elements.
C
C       for elements:
C          if w (e) = 0, then the element e is absorbed
C          if w (e) .ge. wflg, then w (e) - wflg is the size of
C               the set |Le \ Lme|, in terms of nonzeros (the
C               sum of abs (nv (i)) for each principal variable i that
C               is both in the pattern of element e and NOT in the
C               pattern of the current pivot element, me).
C          if wflg .gt. w (e) .gt. 0, then e is not absorbed and has
C               not yet been seen in the scan of the element lists in
C               the computation of |Le\Lme| in loop 150 below.
C
C       for variables:
C          during supervariable detection, if w (j) .ne. wflg then j is
C          not in the pattern of variable i
C
C       The w array is initialized by setting w (i) = 1 for all i,
C       and by setting wflg = 2.  It is reinitialized if wflg becomes
C       too large (to ensure that wflg+n does not cause integer
C       overflow).
C-----------------------------------------------------------------------
C LOCAL INTEGERS:
C-----------------------------------------------------------------------
C     THRESM is used to
C        accelerate symolic factorization 
C         THRESM is dynamically updated to 
C                allow more quasi-dense row selection
C     ThresPrev holds last starting value 
C                   at the beginning of one iteration
C     ThresMin  holds minimum value of THRESH 
      INTEGER :: FDEG, ThresMin, ThresPrev, IBEGSchur, NbSchur, 
     &        ThresMinINIT
      INTEGER :: DEGMAX,THD, THDperm, THD_AGG
      DOUBLE PRECISION :: RELDEN
      LOGICAL :: AGG6_loc, DenseRows
      LOGICAL :: SchurON
      INTEGER :: DEG, DEGME, DEXT, DMAX, E, ELENME, ELN, I,
     &        ILAST, INEXT, J, JLAST, JNEXT, K, KNT1, KNT2, KNT3,
     &        LENJ, LN, ME, MINDEG, NEL, 
     &        NLEFT, NVI, NVJ, NVPIV, SLENME, WE, WFLG, WNVI, X
      INTEGER KNT1_UPDATED, KNT2_UPDATED
      INTEGER :: SIZE_SCHUR_LOC
      INTEGER(8) MAXMEM, MEM, NEWMEM
      INTEGER :: MAXINT_N
      INTEGER(8) :: HASH, HMOD 
      LOGICAL :: COMPRESS
C deg:        the degree of a variable or element
C degme:      size, |Lme|, of the current element, me (= degree (me))
C dext:       external degree, |Le \ Lme|, of some element e
C dmax:       largest |Le| seen so far
C e:          an element
C elenme:     the length, elen (me), of element list of pivotal var.
C eln:        the length, elen (...), of an element list
C hash:       the computed value of the hash function
C hmod:       the hash function is computed modulo hmod = max (1,n-1)
C i:          a supervariable
C ilast:      the entry in a link list preceding i
C inext:      the entry in a link list following i
C j:          a supervariable
C jlast:      the entry in a link list preceding j
C jnext:      the entry in a link list, or path, following j
C k:          the pivot order of an element or variable
C knt1:       loop counter used during element construction
C knt2:       loop counter used during element construction
C knt3:       loop counter used during compression
C lenj:       len (j)
C ln:         length of a supervariable list
C maxint_n:   large integer to test risk of overflow on wflg
C maxmem:     amount of memory needed for no compressions
C me:         current supervariable being eliminated, and the
C                     current element created by eliminating that
C                     supervariable
C mem:        memory in use assuming no compressions have occurred
C mindeg:     current minimum degree
C nel:        number of pivots selected so far
C newmem:     amount of new memory needed for current pivot element
C nleft:      n - nel, the number of nonpivotal rows/columns remaining
C nvi:        the number of variables in a supervariable i (= nv (i))
C nvj:        the number of variables in a supervariable j (= nv (j))
C nvpiv:      number of pivots in current element
C slenme:     number of variables in variable list of pivotal variable
C we:         w (e)
C wflg:       used for flagging the w array.  See description of iw.
C wnvi:       wflg - nv (i)
C x:          either a supervariable or an element
C-----------------------------------------------------------------------
C LOCAL POINTERS:
C-----------------------------------------------------------------------
      INTEGER(8) P, P1, P2, P3, PDST, PEND, PJ, PME, PME1, PME2, 
     &           PN, PSRC, PLN, PELN
C             Any parameter (pe (...) or pfree) or local variable
C             starting with "p" (for Pointer) is an index into iw,
C             and all indices into iw use variables starting with
C             "p."  The only exception to this rule is the iwlen
C             input argument.
C p:          pointer into lots of things
C p1:         pe (i) for some variable i (start of element list)
C p2:         pe (i) + elen (i) -  1 for some var. i (end of el. list)
C p3:         index of first supervariable in clean list
C pdst:       destination pointer, for compression
C pend:       end of memory to compress
C pj:         pointer into an element or variable
C pme:        pointer into the current element (pme1...pme2)
C pme1:       the current element, me, is stored in iw (pme1...pme2)
C pme2:       the end of the current element
C pn:         pointer into a "clean" variable, also used to compress
C psrc:       source pointer, for compression
C-----------------------------------------------------------------------
C  FUNCTIONS CALLED:
C-----------------------------------------------------------------------
      INTRINSIC max, min, mod, maxval
C=======================================================================
C  INITIALIZATIONS
C=======================================================================
      IF (N.EQ.1) THEN
           ELEN(1) = 1
           LAST(1) = 1
           PE(1) = 0_8
           IF (NV(1).LT.0) NV(1) = 1
           NCMPA = 0
           PARENT(1) = 0
           RETURN
      ENDIF
      AGG6_loc = AGG6
      DenseRows = .FALSE.
C
C       We can now assume that N>1
C
CSymbolic  Intialize degrees with the order given by PERM
C
      SIZE_SCHUR_LOC = SIZE_SCHUR
      SIZE_SCHUR_LOC = min(N,SIZE_SCHUR_LOC)
      SIZE_SCHUR_LOC = max(0,SIZE_SCHUR_LOC)
      SchurON   = (SIZE_SCHUR_LOC > 0)
      IBEGSchur = N-SIZE_SCHUR_LOC+1
      THRESM    = THRESH  ! local value of THRESH
      IF (THRESM.GT.N) THRESM = N
      IF (THRESM.LT.0) THRESM = 0
C     Variables in the schur are considered as exactly dense
C     (Schur variables are ordered last, we check it here)
      IF ( SchurON )  THEN 
           DO I= 1, N
             IF ( PERM(I) .GE. IBEGSchur) THEN 
                 PERM(I) = N + 1
C               Because of compress, we force skipping this
C               entry which is anyway empty
                IF (LEN(I) .EQ.0) THEN
                  PE(I) = 0_8
                ENDIF
             ENDIF
           ENDDO
      ENDIF
C          
      IF (SchurON) THEN
C
C         Only restriction is n>= THRESM > 0 
C
C         only exactly dense row will be selected
C         It should also work ok combined to 
C         quasi dense row selection. 
C           (To be Tested it seperately)
             THRESM    = N
             ThresMin  = N
             ThresPrev = N
      ELSE
             THRESM    = max(int(31*N/32),THRESM)
             THRESM    = max(THRESM,1)
C
             DEGMAX= maxval(LEN)
             RELDEN=dble(PFREE-1)/dble(N)
             THD = int(RELDEN)*10 + (DEGMAX-int(RELDEN))/10 + 1
             IF (THD.LT.DEGMAX) THEN
              DenseRows = .TRUE.
              THDperm = N  
              DO I = 1,N
               IF (LEN(I) .GT. THD) THEN
                THDperm =  min(THDperm,PERM(I))
               ENDIF
              ENDDO
              THRESM  = min(THRESM, THDperm)
             ENDIF
C   Compute ThresMin and initialise  ThresPrev
             ThresMin  = max( 3*THRESM / 4, 1)
             ThresPrev = THRESM
C
      ENDIF  ! test on SchurON
C
      ThresMinINIT = ThresMin/4
      THD_AGG = max(128, min(TOTEL/2048, 1024))
      IF (THRESM.GT.0) THEN
       IF ((THRESM.GT.N).OR.(THRESM.LT.2)) THEN 
C      exactly dense rows only
          THRESM = N
       ENDIF
      ENDIF
      LASTD = 0
      NBD   = 0
      NBED  = 0
      NBDM  = 0
      WFLG = 2
      MAXINT_N=huge(WFLG)-TOTEL
      MINDEG = 1
      NCMPA = 0
      NEL = 0
      HMOD = int(max (1, N-1),kind=8)
      DMAX = 0
      MEM = PFREE - 1
      MAXMEM = MEM
      DO I = 1, N
        NDENSE(I)= 0
        W (I) = 1
        ELEN (I) = 0
C        NV (I) = 1
C        DEGREE (I) = LEN (I)
      ENDDO
      DO I=1, N
        LAST (I) = 0
        HEAD (I) = 0
      ENDDO
C     initialize degree
      IF(NV(1) .LT. 0) THEN
         COMPRESS = .FALSE.
      ELSE
         COMPRESS = .TRUE.
      ENDIF
      IF (COMPRESS) THEN
         DO I=1,N
            DEGREE(I) = 0
            DO P= PE(I) , PE(I)+int(LEN(I)-1,8)
               DEGREE(I) = DEGREE(I) + NV(IW(P))
            ENDDO
         ENDDO
      ELSE
         DO I=1,N
            NV(I) = 1
            DEGREE (I) = LEN (I)
         ENDDO
      ENDIF
C     ----------------------------------------------------------------
C     initialize degree lists and eliminate rows with no off-diag. nz.
C     ----------------------------------------------------------------
      DO 20 I = 1, N
        DEG = DEGREE (I)
        IF (PERM(I).EQ.N) THEN
C          save that I is last in the order
           PERMeqN = I
           PERM(I) = N-1
        ENDIF
        FDEG = PERM(I)
        IF ( (DEG .GT. 0).OR.(PERM(I).EQ.N+1) ) THEN
C         ----------------------------------------------------------
C         place i in the degree list corresponding to its degree
C         or in the dense row list if i is dense or quasi dense.
C         ----------------------------------------------------------
C         test for row density
          IF ( (THRESM.GT.0) .AND.
     &         (FDEG .GT.THRESM) ) THEN
C           I will be inserted in the degree list of N
            NBD = NBD+NV(I)
            IF (FDEG.NE.N+1) THEN
C
             DEGREE(I) = DEGREE(I)+TOTEL+2
C            insert I at the beginning of degree list of n
             DEG = N
             INEXT = HEAD (DEG)
             IF (INEXT .NE. 0) LAST (INEXT) = I
             NEXT (I) = INEXT
             HEAD (DEG) = I 
             LAST(I)  = 0
             IF (LASTD.EQ.0) LASTD=I
            ELSE
C            Only Schur variables are concerned here
C            Property: LISTVAR_SCHUR (1) will 
C            be first in the list of schur variables
             NBED = NBED+NV(I)
             DEGREE(I) = TOTEL+1
C            insert I at the end of degree list of n
             DEG = N
             IF (LASTD.EQ.0) THEN
C              degree list is empty
               LASTD     = I 
               HEAD(DEG) = I
               NEXT(I)   = 0 
               LAST(I)   = 0
             ELSE
               NEXT(LASTD) = I
               LAST(I)     = LASTD
               LASTD       = I
               NEXT(I)     = 0
             ENDIF
            ENDIF
          ELSE
C           place i in the degree list corresponding to its degree
            INEXT = HEAD (FDEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            HEAD (FDEG) = I
          ENDIF
        ELSE
C         ----------------------------------------------------------
C         we have a variable that can be eliminated at once because
C         there is no off-diagonal non-zero in its row.
C         ----------------------------------------------------------
          NEL = NEL + NV(I)
          ELEN (I) = -NEL
          PE (I) = 0_8
          W (I) = 0
        ENDIF
   20 CONTINUE
C         We suppress dense row selection if none of them was found in A 
C         in the 1st pass
          IF ((NBD.EQ.0).AND.(THRESM.GT.0)) THRESM = N
C
C=======================================================================
C  WHILE (selecting pivots) DO
C=======================================================================
   30 IF (NEL .LT. TOTEL) THEN
C=======================================================================
C  GET PIVOT OF MINIMUM DEGREE
C=======================================================================
C       -------------------------------------------------------------
C       find next supervariable for elimination
C       -------------------------------------------------------------
        DO 40 DEG = MINDEG, N
          ME = HEAD (DEG)
          IF (ME .GT. 0) GO TO 50
   40   CONTINUE
   50   MINDEG = DEG
C       -------------------------------------------------------------
C       We want to respect the ordering provided by the user
C       Therefefore if (DEG > THRESM .and. NBD.ge.0) then 
C       A quasi-dense variable might have a perm value 
C       smaller than ME. 
C       We thus in this case force restarting.
C       -------------------------------------------------------------
        IF ( (DEG.NE.N) .AND.
     &    (DEG.GT.THRESM+1) .AND. (NBD.GT.0) ) THEN
           MINDEG = N
           GOTO 30
        ENDIF
        IF (DEGREE(ME).LE.TOTEL)  THEN
C       -------------------------------------------------------------
C       remove chosen variable from link list
C       -------------------------------------------------------------
          INEXT = NEXT (ME)
          IF (INEXT .NE. 0) LAST (INEXT) = 0
          HEAD (DEG) = INEXT
        ELSE
C
C         Because of restarting forced even if 
C         variable (not yest quasi dense) but of 
C         value of perm larger thatn thresm still
C         to be eliminated we have to reset MINDEB to 1
          MINDEG = 1
          NBDM = max(NBDM,NBD)
          IF (DEGREE(ME).GT.TOTEL+1) THEN
            IF (WFLG .GT. MAXINT_N) THEN
             DO  52 X = 1, N
              IF (W (X) .NE. 0) W (X) = 1
  52         CONTINUE
             WFLG = 2
            ENDIF
            WFLG = WFLG + 1
  51        CONTINUE
C           ---------------------------------------------------------
C           remove chosen variable from link list
C           ---------------------------------------------------------
            INEXT = NEXT (ME)
            IF (INEXT .NE. 0) THEN 
               LAST (INEXT) = 0
            ELSE
               LASTD = 0
            ENDIF
C           ----------------------------------------------------------
c           build adjacency list of ME in quotient gragh
C           and calculate its external degree  in ndense(me)
C           ----------------------------------------------------------
            NDENSE(ME) = 0
            W(ME)      = WFLG
            P1 = PE(ME)
            P2 = P1 + int(LEN(ME) -1,8)
C           PLN-1 holds the pointer in IW to the last elet/var in adj list
C              of ME.  LEN(ME) will then be set to PLN-P1
C           PELN-1 hold the pointer in IW to the last elet in in adj list
C              of ME.  ELEN(ME) will then be set to PELN-P1
C           element adjacent to ME
            PLN       = P1
            PELN      = P1
            DO 55 P=P1,P2
              E= IW(P)
              IF (W(E).EQ.WFLG) GOTO 55
              W(E) = WFLG
              IF (PE(E).LT.0_8) THEN
C              E is a nonprincipal variable or absorbed element
                X = E
  53            X = int(-PE(X))
                IF (W(X) .EQ.WFLG) GOTO 55
                W(X) = WFLG
                IF ( PE(X) .LT. 0_8 ) GOTO 53
                E = X
              ENDIF
C             -------------------------------------------
C             E is an unabsorbed element or a "dense" row
C                 (NOT already flagged)
C             -------------------------------------------
              IF (ELEN(E).LT.0) THEN
C              E is a new element in adj(ME)
               NDENSE(E) = NDENSE(E) - NV(ME)
               IW(PLN) = IW(PELN)  
               IW(PELN) = E
               PLN  = PLN+1_8
               PELN = PELN + 1_8
C              update ndense of ME with all unflagged dense
C              rows in E
               PME1 = PE(E)
               DO 54 PME = PME1, PME1+int(LEN(E)-1,8)
                X = IW(PME)
                IF ((ELEN(X).GE.0).AND.(W(X).NE.WFLG)) THEN
C                X is a dense row 
                 NDENSE(ME) = NDENSE(ME) + NV(X)
                 W(X) = WFLG
                ENDIF
 54            CONTINUE
              ELSE
C              E is a dense row 
               NDENSE(ME) = NDENSE(ME) + NV(E)
               IW(PLN)=E
               PLN = PLN+1_8
              ENDIF
  55        CONTINUE
C           ----------------------------------------------
C           DEGREE(ME)-(TOTEL+2) holds last external degree computed
C           when Me was detected as dense
C           NDENSE(ME) is the exact external degree of ME
C           ----------------------------------------------
            WFLG     = WFLG + 1
            LEN(ME)  = int(PLN-P1)
            ELEN(ME) = int(PELN-P1)
            NDME = NDENSE(ME)+NV(ME)
            IF (NDENSE(ME).EQ.0) NDENSE(ME) =1
C           ---------------------------------------------------------
C           place ME in the degree list of NDENSE(ME), update DEGREE
C           ---------------------------------------------------------
            DEGREE(ME) = NDENSE(ME)
            DEG = PERM(ME)
            MINDEG = min(DEG,MINDEG)
            JNEXT = HEAD(DEG)
            IF (JNEXT.NE. 0) LAST (JNEXT) = ME
            NEXT(ME) = JNEXT
            HEAD(DEG) = ME
C           ------------------------------
C           process next quasi dense row
C           ------------------------------
            ME    = INEXT
            IF (ME.NE.0) THEN
              IF (DEGREE(ME).GT.(TOTEL+1) ) GOTO 51
            ENDIF
            HEAD (N) = ME
C           ---------------------------------------
C           update dense row selection strategy
C           -------------------------------------
            IF (THRESM.LT.N) THEN
             ThresMin  = max(THRESM+ThresMin,ThresPrev+ThresMin/2+1)
             ThresMin  = min(ThresMin, N)
             ThresPrev = ThresPrev+(N-ThresPrev)/2+ThresMinINIT
             THRESM    = max(
     &         THRESM + int(sqrt(dble(ThresMin)))+ ThresMinINIT ,
     &         ThresPrev)
             THRESM    = min(THRESM,N) 
             ThresMin  = min(THRESM, ThresMin)
             ThresPrev = THRESM
            ENDIF
            NBD    = NBED
C           get back to Min degree elimination loop
C
            GOTO 30
          ENDIF
C         -------------------------------------------------------------
C         -------------------------------------------------------------
          IF (DEGREE(ME).EQ.TOTEL+1) THEN
C         we have only  exactly "dense" rows that we
C         amalgamate at the root node
             IF (NBD.NE.NBED) THEN
          write(6,*) ' ERROR in MUMPS_SYMQAMD quasi dense rows remains'
          CALL MUMPS_ABORT()
           ENDIF
           NbSchur = 0   ! Only for checking
           NELME    = -(NEL+1)
           DO 59 X=1,N
            IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
              PE(X) = int(-LISTVAR_SCHUR(1),8)
            ELSE IF ((PE(X).GT.0) .AND. (ELEN(X).LT.0)) THEN
C            X is an unabsorbed element
C            -- Force sons to be linked to first node in Schur
             PE(X) = int(-LISTVAR_SCHUR(1),8)
C            W(X) = 0 could be suppressed ?? check it
            ELSEIF (DEGREE(X).EQ.TOTEL+1) THEN
C            X is a dense row, absorb it in ME (mass elimination)
             NEL   = NEL + NV(X)
             PE(X) = int(-ME,8)
             ELEN(X) = 0
             NV(X) = 0
             NbSchur = NbSchur+ 1
            ENDIF
   59      CONTINUE
           IF (NbSchur.NE.SIZE_SCHUR_LOC) then
             write(6,*) ' Internal error 2 in QAMD :',
     &         ' Schur size expected:',SIZE_SCHUR_LOC, 'Real:', NbSchur
             CALL MUMPS_ABORT()
           ENDIF
C          ME is the root node
           ELEN(ME) = NELME
           NV(ME)   = NBD
           PE(ME)   = 0_8
           IF (NEL.NE.N) THEN
            write(6,*) 'Internal ERROR 2 detected in QAMD'
            write(6,*) ' NEL not equal to N: N, NEL =',N,NEL
            CALL MUMPS_ABORT()
           ENDIF
           IF (ME.NE. LISTVAR_SCHUR(1)) THEN
C          -- Set all node in Schur list to point to LISTVAR_SCHUR(1)
             DO I=1, SIZE_SCHUR_LOC
               PE(LISTVAR_SCHUR(I)) = int(-LISTVAR_SCHUR(1),8)
             ENDDO
             PE(LISTVAR_SCHUR(1)) = 0_8
             NV( LISTVAR_SCHUR(1))= NV(ME)
             NV(ME)               = 0
             ELEN( LISTVAR_SCHUR(1)) = ELEN(ME)
             ELEN(ME)             = 0
           ENDIF
           GOTO 265
          ENDIF
        ENDIF
C       -------------------------------------------------------------
C       me represents the elimination of pivots nel+1 to nel+nv(me).
C       place me itself as the first in this set.  It will be moved
C       to the nel+nv(me) position when the permutation vectors are
C       computed.
C       -------------------------------------------------------------
        ELENME = ELEN (ME)
        ELEN (ME) = - (NEL + 1)
        NVPIV = NV (ME)
        NEL = NEL + NVPIV
        NDENSE(ME) = 0
C=======================================================================
C  CONSTRUCT NEW ELEMENT
C=======================================================================
C       -------------------------------------------------------------
C       At this point, me is the pivotal supervariable.  It will be
C       converted into the current element.  Scan list of the
C       pivotal supervariable, me, setting tree pointers and
C       constructing new list of supervariables for the new element,
C       me.  p is a pointer to the current position in the old list.
C       -------------------------------------------------------------
C       flag the variable "me" as being in Lme by negating nv (me)
        NV (ME) = -NVPIV
        DEGME = 0
        IF (ELENME .EQ. 0) THEN
C         ----------------------------------------------------------
C         construct the new element in place
C         ----------------------------------------------------------
          PME1 = PE (ME)
          PME2 = PME1 - 1
          DO 60 P = PME1, PME1 + int(LEN (ME) - 1,8)
            I = IW (P)
            NVI = NV (I)
            IF (NVI .GT. 0) THEN
C             ----------------------------------------------------
C             i is a principal variable not yet placed in Lme.
C             store i in new list
C             ----------------------------------------------------
              DEGME = DEGME + NVI
C             flag i as being in Lme by negating nv (i)
              NV (I) = -NVI
              PME2 = PME2 + 1
              IW (PME2) = I
C             ----------------------------------------------------
C             remove variable i from degree list.
C             ----------------------------------------------------
C             only done for non "dense" rows
              IF (DEGREE(I).LE.TOTEL) THEN
              ILAST = LAST (I)
              INEXT = NEXT (I)
              IF (INEXT .NE. 0) LAST (INEXT) = ILAST
              IF (ILAST .NE. 0) THEN
                NEXT (ILAST) = INEXT
              ELSE
C               i is at the head of the degree list
                HEAD (PERM(I)) = INEXT
              ENDIF
              ELSE
               NDENSE(ME) = NDENSE(ME) + NVI
              ENDIF
            ENDIF
   60     CONTINUE
C         this element takes no new memory in iw:
          NEWMEM = 0
        ELSE
C         ----------------------------------------------------------
C         construct the new element in empty space, iw (pfree ...)
C         ----------------------------------------------------------
          P = PE (ME)
          PME1 = PFREE
          SLENME = LEN (ME) - ELENME
          KNT1_UPDATED = 0
          DO 120 KNT1 = 1, ELENME + 1
            KNT1_UPDATED = KNT1_UPDATED +1
            IF (KNT1 .GT. ELENME) THEN
C             search the supervariables in me.
              E = ME
              PJ = P
              LN = SLENME
            ELSE
C             search the elements in me.
              E = IW (P)
              P = P + 1
              PJ = PE (E)
              LN = LEN (E)
            ENDIF
C           -------------------------------------------------------
C           search for different supervariables and add them to the
C           new list, compressing when necessary. this loop is
C           executed once for each element in the list and once for
C           all the supervariables in the list.
C           -------------------------------------------------------
            KNT2_UPDATED = 0
            DO 110 KNT2 = 1, LN
              KNT2_UPDATED = KNT2_UPDATED+1
              I = IW (PJ)
              PJ = PJ + 1
              NVI = NV (I)
              IF (NVI .GT. 0) THEN
C               -------------------------------------------------
C               compress iw, if necessary
C               -------------------------------------------------
                IF (PFREE .GT. IWLEN) THEN
C                 prepare for compressing iw by adjusting
C                 pointers and lengths so that the lists being
C                 searched in the inner and outer loops contain
C                 only the remaining entries.
                  PE (ME) = P
                  LEN (ME) = LEN (ME) - KNT1_UPDATED
C                 Reset KNT1_UPDATED in case of recompress 
C                 at same iteration of the loop 120
                  KNT1_UPDATED = 0
C                 Check if anything left in supervariable ME
                  IF (LEN (ME) .EQ. 0) PE (ME) = 0
                  PE (E) = PJ
                  LEN (E) = LN - KNT2_UPDATED
C                 Reset KNT2_UPDATED in case of recompress 
C                 at same iteration of the loop 110
                  KNT2_UPDATED = 0
C                 Check if anything left in element E
                  IF (LEN (E) .EQ. 0) PE (E) = 0
                  NCMPA = NCMPA + 1
C                 store first item in pe
C                 set first entry to -item
                  DO 70 J = 1, N
                    PN = PE (J)
                    IF (PN .GT. 0) THEN
                      PE (J) = int(IW (PN),8)
                      IW (PN) = -J
                    ENDIF
   70             CONTINUE
C                 psrc/pdst point to source/destination
                  PDST = 1
                  PSRC = 1
                  PEND = PME1 - 1
C                 while loop:
   80             CONTINUE
                  IF (PSRC .LE. PEND) THEN
C                   search for next negative entry
                    J = -IW (PSRC)
                    PSRC = PSRC + 1
                    IF (J .GT. 0) THEN
                      IW (PDST) = int(PE (J))
                      PE (J) = PDST
                      PDST = PDST + 1
C                     copy from source to destination
                      LENJ = LEN (J)
                      DO 90 KNT3 = 0, LENJ - 2
                        IW (PDST + KNT3) = IW (PSRC + KNT3)
   90                 CONTINUE
                      PDST = PDST + LENJ - 1
                      PSRC = PSRC + LENJ - 1
                    ENDIF
                    GO TO 80
                  ENDIF
C                 move the new partially-constructed element
                  P1 = PDST
                  DO 100 PSRC = PME1, PFREE - 1
                    IW (PDST) = IW (PSRC)
                    PDST = PDST + 1
  100             CONTINUE
                  PME1 = P1
                  PFREE = PDST
                  PJ = PE (E)
                  P = PE (ME)
                ENDIF
C               -------------------------------------------------
C               i is a principal variable not yet placed in Lme
C               store i in new list
C               -------------------------------------------------
                DEGME = DEGME + NVI
C               flag i as being in Lme by negating nv (i)
                NV (I) = -NVI
                IW (PFREE) = I
                PFREE = PFREE + 1
C               -------------------------------------------------
C               remove variable i from degree link list
C               -------------------------------------------------
C             only done for non "dense" rows
                IF (DEGREE(I).LE.TOTEL) THEN
                ILAST = LAST (I)
                INEXT = NEXT (I)
                IF (INEXT .NE. 0) LAST (INEXT) = ILAST
                IF (ILAST .NE. 0) THEN
                  NEXT (ILAST) = INEXT
                ELSE
C                 i is at the head of the degree list
                  HEAD (PERM(I)) = INEXT
                ENDIF
                ELSE
                 NDENSE(ME) = NDENSE(ME) + NVI
                ENDIF
              ENDIF
  110       CONTINUE
            IF (E .NE. ME) THEN
C             set tree pointer and flag to indicate element e is
C             absorbed into new element me (the parent of e is me)
              PE (E) = int(-ME,8)
              W (E) = 0
            ENDIF
  120     CONTINUE
          PME2 = PFREE - 1
C         this element takes newmem new memory in iw (possibly zero)
          NEWMEM = PFREE - PME1
          MEM = MEM + NEWMEM
          MAXMEM = max (MAXMEM, MEM)
        ENDIF
C       -------------------------------------------------------------
C       me has now been converted into an element in iw (pme1..pme2)
C       -------------------------------------------------------------
C       degme holds the external degree of new element
        DEGREE (ME) = DEGME
        PE (ME) = PME1
        LEN (ME) = int(PME2 - PME1 + 1_8)
C       -------------------------------------------------------------
C       make sure that wflg is not too large.  With the current
C       value of wflg, wflg+n must not cause integer overflow
C       -------------------------------------------------------------
        IF (WFLG .GT. MAXINT_N) THEN
          DO 130 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  130     CONTINUE
          WFLG = 2
        ENDIF
C=======================================================================
C  COMPUTE (w (e) - wflg) = |Le\Lme| FOR ALL ELEMENTS
Cdense
C   COMPUTE (w(e) - wflg) = |Le(G')\Lme(G')| FOR ALL ELEMENTS
C   where G' is the subgraph of G excluding ''dense" rows)
Cdense
C=======================================================================
C       -------------------------------------------------------------
C       Scan 1:  compute the external degrees of previous elements
C       with respect to the current element.  That is:
C            (w (e) - wflg) = |Le \ Lme|
C       for each element e that appears in any supervariable in Lme.
C       The notation Le refers to the pattern (list of
C       supervariables) of a previous element e, where e is not yet
C       absorbed, stored in iw (pe (e) + 1 ... pe (e) + iw (pe (e))).
C       The notation Lme refers to the pattern of the current element
C       (stored in iw (pme1..pme2)).   If (w (e) - wflg) becomes
C       zero, then the element e will be absorbed in scan 2.
C       aggressive absorption is possible only if NDENSE(ME) = NBD
C       which is true when only exactly dense rows have been selected.
C       -------------------------------------------------------------
        DO 150 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).GT.TOTEL) GOTO 150
          ELN = ELEN (I)
          IF (ELN .GT. 0) THEN
C           note that nv (i) has been negated to denote i in Lme:
            NVI = -NV (I)
            WNVI = WFLG - NVI
            DO 140 P = PE (I), PE (I) + int(ELN - 1,8)
              E = IW (P)
              WE = W (E)
              IF (WE .GE. WFLG) THEN
C               unabsorbed element e has been seen in this loop
                WE = WE - NVI
              ELSE IF (WE .NE. 0) THEN
C               e is an unabsorbed element
C               this is the first we have seen e in all of Scan 1
                WE = DEGREE (E) + WNVI - NDENSE(E)
Cn dense
              ENDIF
              W (E) = WE
  140       CONTINUE
          ENDIF
  150   CONTINUE
C=======================================================================
C  DEGREE UPDATE AND ELEMENT ABSORPTION
C=======================================================================
C       -------------------------------------------------------------
C       Scan 2:  for each i in Lme, sum up the degree of Lme (which
C       is degme), plus the sum of the external degrees of each Le
C       for the elements e appearing within i, plus the
C       supervariables in i.  Place i in hash list.
C       -------------------------------------------------------------
        AGG6_loc = (AGG6 .OR. (DEGREE(ME) .LT. THD_AGG))
        DO 180 PME = PME1, PME2
          I = IW (PME)
          IF (DEGREE(I).GT.TOTEL) GOTO 180
          P1 = PE (I)
          P2 = P1 + int(ELEN (I) - 1,8)
          PN = P1
          HASH = 0_8
          DEG = 0
C         ----------------------------------------------------------
C         scan the element list associated with supervariable i
C         ----------------------------------------------------------
          DO 160 P = P1, P2
            E = IW (P)
C           dext = | Le \ Lme |
            DEXT = W (E) - WFLG
            IF (DEXT .GT. 0) THEN
              DEG = DEG + DEXT
              IW (PN) = E
              PN = PN + 1_8
              HASH = HASH + int(E,kind=8)
C        ------------------------------
C        suppress aggressive absorption
C        ------------------------------
            ELSE IF (.NOT. AGG6_loc .AND. DEXT .EQ. 0) THEN
              IW (PN) = E
              PN = PN + 1_8
              HASH = HASH + int(E,kind=8)
C
C        ------------------------------
C        try aggressive absorption
C         when possible
            ELSE IF (AGG6_loc .AND. (DEXT .EQ. 0) .AND.
     &            ((NDENSE(ME).EQ.NBD).OR.(NDENSE(E).EQ.0))) THEN
C             aggressive absorption: e is not adjacent to me, but
C             |Le(G') \ Lme(G')| is 0 and all dense rows
C             are in me, so absorb it into me
                PE (E) = int(-ME,8)
                W (E)  = 0
             ELSE IF (AGG6_loc .AND. DEXT.EQ.0) THEN
                  IW(PN) = E
                  PN     = PN+1
                  HASH   = HASH + int(E,kind=8)
            ENDIF
  160     CONTINUE
C         count the number of elements in i (including me):
          ELEN (I) = int(PN - P1 + 1)
C         ----------------------------------------------------------
C         scan the supervariables in the list associated with i
C         ----------------------------------------------------------
          P3 = PN
          DO 170 P = P2 + 1, P1 + int(LEN (I) - 1,8)
            J = IW (P)
            NVJ = NV (J)
            IF (NVJ .GT. 0) THEN
C             j is unabsorbed, and not in Lme.
C             add to degree and add to new list
C             add degree only of non-dense rows.
              IF (DEGREE(J).LE.TOTEL) DEG=DEG+NVJ
              IW (PN) = J
              PN = PN + 1
              HASH = HASH + int(J,kind=8)
            ENDIF
  170     CONTINUE
C         ----------------------------------------------------------
C         update the degree and check for mass elimination
C         ----------------------------------------------------------
          IF (((ELEN(I).EQ.1).AND.(P3.EQ.PN))
     &     .OR.
     &         (AGG6_loc.AND.(DEG .EQ. 0).AND.(NDENSE(ME).EQ.NBD))
     &       )
     &    THEN
C           -------------------------------------------------------
C           mass elimination
C           -------------------------------------------------------
C           There is nothing left of this node except for an
C           edge to the current pivot element.  elen (i) is 1,
C           and there are no variables adjacent to node i.
C           Absorb i into the current pivot element, me.
            PE (I) = int(-ME,8)
            NVI = -NV (I)
            DEGME = DEGME - NVI
            NVPIV = NVPIV + NVI
            NEL = NEL + NVI
            NV (I) = 0
            ELEN (I) = 0
          ELSE
C           -------------------------------------------------------
C           update the upper-bound degree of i
C           -------------------------------------------------------
C           the following degree does not yet include the size
C           of the current element, which is added later:
            DEGREE(I) = min (DEG+NBD-NDENSE(ME), 
     &                       DEGREE(I))
C           -------------------------------------------------------
C           add me to the list for i
C           -------------------------------------------------------
C           move first supervariable to end of list
            IW (PN) = IW (P3)
C           move first element to end of element part of list
            IW (P3) = IW (P1)
C           add new element to front of list.
            IW (P1) = ME
C           store the new length of the list in len (i)
            LEN (I) = int(PN - P1 + 1)
C           -------------------------------------------------------
C           place in hash bucket.  Save hash key of i in last (i).
C           -------------------------------------------------------
            HASH = mod (HASH, HMOD) + 1_8
            J = HEAD (HASH)
            IF (J .LE. 0) THEN
C             the degree list is empty, hash head is -j
              NEXT (I) = -J
              HEAD (HASH) = -I
            ELSE
C             degree list is not empty
C             use last (head (hash)) as hash head
              NEXT (I) = LAST (J)
              LAST (J) = I
            ENDIF
            LAST (I) = int(HASH,kind=kind(LAST))
          ENDIF
  180   CONTINUE
        DEGREE (ME) = DEGME
C       -------------------------------------------------------------
C       Clear the counter array, w (...), by incrementing wflg.
C       -------------------------------------------------------------
        DMAX = max (DMAX, DEGME)
        WFLG = WFLG + DMAX
C       make sure that wflg+n does not cause integer overflow
        IF (WFLG .GT. MAXINT_N) THEN
          DO 190 X = 1, N
            IF (W (X) .NE. 0) W (X) = 1
  190     CONTINUE
          WFLG = 2
        ENDIF
C       at this point, w (1..n) .lt. wflg holds
C=======================================================================
C  SUPERVARIABLE DETECTION
C=======================================================================
        DO 250 PME = PME1, PME2
          I = IW (PME)
          IF ( (NV(I).LT.0) .AND. (DEGREE(I).LE.TOTEL) ) THEN
C           only done for nondense rows
C           i is a principal variable in Lme
C           -------------------------------------------------------
C           examine all hash buckets with 2 or more variables.  We
C           do this by examing all unique hash keys for super-
C           variables in the pattern Lme of the current element, me
C           -------------------------------------------------------
            HASH = int(LAST (I),kind=8)
C           let i = head of hash bucket, and empty the hash bucket
            J = HEAD (HASH)
            IF (J .EQ. 0) GO TO 250
            IF (J .LT. 0) THEN
C             degree list is empty
              I = -J
              HEAD (HASH) = 0
            ELSE
C             degree list is not empty, restore last () of head
              I = LAST (J)
              LAST (J) = 0
            ENDIF
            IF (I .EQ. 0) GO TO 250
C           while loop:
  200       CONTINUE
            IF (NEXT (I) .NE. 0) THEN
             X = I 
C             ----------------------------------------------------
C             this bucket has one or more variables following i.
C             scan all of them to see if i can absorb any entries
C             that follow i in hash bucket.  Scatter i into w.
C             ----------------------------------------------------
              LN = LEN (I)
              ELN = ELEN (I)
C             do not flag the first element in the list (me)
              DO 210 P = PE (I) + 1, PE (I) + int(LN - 1,8)
                W (IW (P)) = WFLG
  210         CONTINUE
C             ----------------------------------------------------
C             scan every other entry j following i in bucket
C             ----------------------------------------------------
              JLAST = I
              J = NEXT (I)
C             while loop:
  220         CONTINUE
              IF (J .NE. 0) THEN
C               -------------------------------------------------
C               check if j and i have identical nonzero pattern
C               -------------------------------------------------
C               jump if i and j do not have same size data structure
                IF (LEN (J) .NE. LN) GO TO 240
C               jump if i and j do not have same number adj elts
                IF (ELEN (J) .NE. ELN) GO TO 240
C               do not flag the first element in the list (me)
                DO 230 P = PE (J) + 1, PE (J) + int(LN - 1,8)
C                 jump if an entry (iw(p)) is in j but not in i
                  IF (W (IW (P)) .NE. WFLG) GO TO 240
  230           CONTINUE
C               -------------------------------------------------
C               found it!  j can be absorbed into i
C               -------------------------------------------------
                IF (PERM(J).GT.PERM(X)) THEN
                ! J is absorbed by X
                  PE (J) = int(-X,8)
                  NV (X) = NV (X) + NV (J)
                  NV (J) = 0
                  ELEN (J) = 0
                ELSE
                ! X is absorbed by J
                  PE (X) = int(-J,8)
                  NV (J) = NV (X) + NV (J)
                  NV (X) = 0
                  ELEN (X) = 0
                  X = J
                ENDIF
C               both nv (i) and nv (j) are negated since they
C               are in Lme, and the absolute values of each
C               are the number of variables in i and j:
C               delete j from hash bucket
                J = NEXT (J)
                NEXT (JLAST) = J
                GO TO 220
C               -------------------------------------------------
  240           CONTINUE
C               j cannot be absorbed into i
C               -------------------------------------------------
                JLAST = J
                J = NEXT (J)
              GO TO 220
              ENDIF
C             ----------------------------------------------------
C             no more variables can be absorbed into i
C             go to next i in bucket and clear flag array
C             ----------------------------------------------------
              WFLG = WFLG + 1
              I = NEXT (I)
              IF (I .NE. 0) GO TO 200
            ENDIF
          ENDIF
  250   CONTINUE
C=======================================================================
C  RESTORE DEGREE LISTS AND REMOVE NONPRINCIPAL SUPERVAR. FROM ELEMENT
C=======================================================================
C       ------------------------------
C       Update thresm for having more
C       quasi dense rows to select
C       ------------------------------
        IF ( .NOT.DenseRows.AND.(THRESM .GT. 0).AND.(THRESM.LT.N) ) 
     &       THEN 
          THRESM = max(ThresMin, THRESM-NVPIV)
        ENDIF
        P = PME1
        NLEFT = TOTEL - NEL
        DO 260 PME = PME1, PME2
          I = IW (PME)
          NVI = -NV (I)
          IF (NVI .GT. 0) THEN
C           i is a principal variable in Lme
C           restore nv (i) to signify that i is principal
            NV (I) = NVI
            IF (DEGREE(I).LE.TOTEL) THEN
C           -------------------------------------------------------
C           compute the external degree (add size of current elem)
C           -------------------------------------------------------
            DEG = min (DEGREE (I)+ DEGME - NVI, NLEFT - NVI)
            DEGREE (I) = DEG
            IDENSE = .FALSE.
C           
C           -------------------
C           Dense row detection
C           -------------------
            IF (THRESM.GT.0) THEN
             IF (PERM(I) .GT. THRESM) THEN
C             relaxed dense row detection
               IDENSE = .TRUE.
C
               DEGREE(I) = DEGREE(I)+TOTEL+2
             ENDIF
             IF (IDENSE) THEN
C            update NDENSE of all elements in the list of element
C            adjacent to I (including ME).
               P1 = PE(I)
               P2 = P1 + int(ELEN(I) - 1,8)
               IF (P2.GE.P1) THEN
               DO 264 PJ=P1,P2
                 E= IW(PJ)
                 NDENSE (E) = NDENSE(E) + NVI
 264           CONTINUE
               ENDIF
C            insert I in the list of dense rows
               NBD = NBD+NVI
               FDEG = N
               DEG = N
C              insert I at the beginning of the list
               INEXT = HEAD(DEG)
               IF (INEXT .NE. 0) LAST (INEXT) = I
               NEXT (I) = INEXT
               HEAD (DEG) = I
               LAST(I)    = 0
               IF (LASTD.EQ.0) LASTD=I
C            end of IDENSE=true
             ENDIF
C           end of THRESM>0
            ENDIF
C             
            IF (.NOT.IDENSE) THEN
            FDEG = PERM(I)
C           -------------------------------------------------------
C           place the supervariable at the head of the degree list
C           -------------------------------------------------------
            INEXT = HEAD (FDEG)
            IF (INEXT .NE. 0) LAST (INEXT) = I
            NEXT (I) = INEXT
            LAST (I) = 0
            HEAD (FDEG) = I
            ENDIF
C           -------------------------------------------------------
C           save the new degree, and find the minimum degree
C           -------------------------------------------------------
            MINDEG = min (MINDEG, FDEG)
            ENDIF
C           -------------------------------------------------------
C           place the supervariable in the element pattern
C           -------------------------------------------------------
            IW (P) = I
            P = P + 1
          ENDIF
  260   CONTINUE
C=======================================================================
C  FINALIZE THE NEW ELEMENT
C=======================================================================
        NV (ME) = NVPIV + DEGME
C       nv (me) is now the degree of pivot (including diagonal part)
C       save the length of the list for the new element me
        LEN (ME) = int(P - PME1)
        IF (LEN (ME) .EQ. 0) THEN
C         there is nothing left of the current pivot element
          PE (ME) = 0_8
          W (ME) = 0
        ENDIF
        IF (NEWMEM .NE. 0) THEN
C         element was not constructed in place: deallocate part
C         of it (final size is less than or equal to newmem,
C         since newly nonprincipal variables have been removed).
          PFREE = P
          MEM = MEM - NEWMEM + int(LEN (ME),8)
        ENDIF
C=======================================================================
C       END WHILE (selecting pivots)
      GO TO 30
      ENDIF
C=======================================================================
  265 CONTINUE
C=======================================================================
C  COMPUTE THE PERMUTATION VECTORS
C=======================================================================
C     ----------------------------------------------------------------
C     The time taken by the following code is O(n).  At this
C     point, elen (e) = -k has been done for all elements e,
C     and elen (i) = 0 has been done for all nonprincipal
C     variables i.  At this point, there are no principal
C     supervariables left, and all elements are absorbed.
C     ----------------------------------------------------------------
C     ----------------------------------------------------------------
C     compute the ordering of unordered nonprincipal variables
C     ----------------------------------------------------------------
      DO 290 I = 1, N
        IF (ELEN (I) .EQ. 0) THEN
C         ----------------------------------------------------------
C         i is an un-ordered row.  Traverse the tree from i until
C         reaching an element, e.  The element, e, was the
C         principal supervariable of i and all nodes in the path
C         from i to when e was selected as pivot.
C         ----------------------------------------------------------
          J = int(-PE (I))
C         while (j is a variable) do:
  270     CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              J = int(-PE (J))
              GO TO 270
            ENDIF
            E = J
C           ----------------------------------------------------------
C           get the current pivot ordering of e
C           ----------------------------------------------------------
            K = -ELEN (E)
C           ----------------------------------------------------------
C           traverse the path again from i to e, and compress the
C           path (all nodes point to e).  Path compression allows
C           this code to compute in O(n) time.  Order the unordered
C           nodes in the path, and place the element e at the end.
C           ----------------------------------------------------------
            J = I
C           while (j is a variable) do:
  280       CONTINUE
            IF (ELEN (J) .GE. 0) THEN
              JNEXT = int(-PE (J))
              PE (J) = int(-E,8)
              IF (ELEN (J) .EQ. 0) THEN
C               j is an unordered row
                ELEN (J) = K
                K = K + 1
              ENDIF
              J = JNEXT
            GO TO 280
            ENDIF
C         leave elen (e) negative, so we know it is an element
          ELEN (E) = -K
        ENDIF
  290 CONTINUE
C     ----------------------------------------------------------------
C     reset the inverse permutation (elen (1..n)) to be positive,
C     and compute the permutation (last (1..n)).
C     ----------------------------------------------------------------
      DO 300 I = 1, N
        K = abs (ELEN (I))
C        LAST (K) = I
C        LAST (K) = I
        ELEN (I) = K
  300 CONTINUE
      IF (.NOT.SchurON) THEN
C       -----------------------------
C       restore PERM(I)=N  for PERMeqN
C       -----------------------------
        PERM(PERMeqN) = N
      ENDIF
C=======================================================================
C  RETURN THE MEMORY USAGE IN IW
C=======================================================================
C     If maxmem is less than or equal to iwlen, then no compressions
C     occurred, and iw (maxmem+1 ... iwlen) was unused.  Otherwise
C     compressions did occur, and iwlen would have had to have been
C     greater than or equal to maxmem for no compressions to occur.
C     Return the value of maxmem in the pfree argument.
      PFREE = MAXMEM
C===============================
C     Save PE in PARENT array
      DO I=1,N
       PARENT(I) = int(PE(I))
      ENDDO
C===============================
      RETURN
      END SUBROUTINE MUMPS_SYMQAMD
