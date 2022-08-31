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
C  ===========================
C  FORTRAN 90 Driver for DMUMPS
C       (MPI based code) 
C  ===========================
C
      SUBROUTINE DMUMPS( id )
      USE DMUMPS_OOC
      USE MUMPS_MEMORY_MOD
      USE DMUMPS_STRUC_DEF
      USE DMUMPS_STATIC_PTR_M ! For Schur pointer
      USE DMUMPS_SAVE_RESTORE
C      
!$    USE OMP_LIB
C
      IMPLICIT NONE
C
C  =======
C  Purpose
C  =======
C
C  TO SOLVE a SPARSE SYSTEM OF LINEAR EQUATIONS.
C  GIVEN AN UNSYMMETRIC, SYMMETRIC, OR SYMMETRIC POSITIVE DEFINITE 
C  SPARSE MATRIX A AND AN N-VECTOR B, THIS SUBROUTINE SOLVES THE 
C  SYSTEM A x = b or ATRANSPOSE x = b. 
C
C  List of main functionalities provided by the package:
C  ----------------------------------------------------
C        -Unsymmetric solver with partial pivoting (LU factorization)
C        -Symmetric positive definite solver (LDLT factorization)
C        -General symmetric solver with pivoting
C        -Either elemental or assembled matrix input
C        -Analysis/Factorization/Solve callable separately
C        -Deficient matrices (symmetric or unsymmetric) 
C          -Rank revealing 
C          -Null space basis computation
C          -Solution 
C        -Return the Schur complement matrix while 
C          also providing solution of interior problem
C        -Distributed input matrix and analysis phase
C        -Sequential or parallel MPI version (any number of processors)
C        -Error analysis and iterative refinement
C        -Out-of-Core factorization and solution
C        -Solution phase:
C          -Multiple Right-Hand-sides (RHS)
C          -Sparse RHS
C          -Distributed RHS
C          -Computation of selected entries of the inverse of 
C           original matrix.
C        - Block Low-Rank (BLR) approximation based factorization
C
C Method
C ------
C  The method used is a parallel direct method
C  based on a sparse multifrontal variant
C  of Gaussian elimination with partial numerical pivoting. 
C  An initial ordering for the pivotal sequence
C  is chosen using the pattern of the matrix A + A^T and is
C  later modified for reasons of numerical stability.  Thus this code
C  performs best on matrices whose pattern is symmetric, or nearly so.
C  For symmetric sparse matrices or for very unsymmetric and
C  very sparse matrices, other software might be more appropriate.
C
C
C References :
C -----------
C
C  P. Amestoy, J.-Y. L'Excellent, G. Moreau, On exploiting sparsity of 
C    multiple right-hand sides in sparse direct solvers,
C    SIAM Journal on Scientific Computing, volume 41, number 2,
C    pages A269-A291 (2019)
C
C  G. Moreau, PhD Thesis, ENS-Lyon, University of Lyon,  
C    On the solution phase of direct methods for sparse linear systems
C    with multiple sparse right-hand sides, December 10th, 2018 
C
C  P. Amestoy, A. Buttari, J.-Y. L'Excellent and T. Mary, 
C    Performance and scalability of the block low-rank multifrontal 
C    factorization on multicore architectures, 
C    ACM Transactions on Mathematical Software (2018)
C
C  T. Mary, PhD Thesis, University of Toulouse,
C    Block Low-Rank multifrontal solvers: complexity, performance, and 
C    scalability, November 2017.
C
C  S. de la Kethulle de Ryhove, P. Jaysaval and D.V. Shantsev, 
C    P. R. Amestoy, J.-Y. L'Excellent and T. Mary, 
C    Large-scale 3D EM modeling with a Block Low-Rank MUMPS solver, 
C    Geophysical Journal International, volume 209, number 3, 
C    pages 1558-1571 (2017) .
C
C  P. Amestoy, A. Buttari, J.-Y. L'Excellent and T. Mary, 
C    On the complexity of the Block Low-Rank multifrontal factorization, 
C    SIAM Journal on Scientific Computing, volume 39, 
C    number 4, pages A1710-A1740 (2017).
C
C  P. Amestoy, R. Brossier, A. Buttari, J.-Y. L'Excellent, T. Mary,
C   L. Metivier, A. Miniussi, and S. Operto.
C   Fast 3D frequency-domain full waveform inversion with a parallel
C   Block Low-Rank multifrontal direct solver: application to OBC data 
C   from the North Sea, Geophysics, 81(6):R363--R383, (2016).
C  
C  P. Amestoy, C. Ashcraft, O. Boiteau, A. Buttari, J.-Y. L'Excellent, 
C   and C. Weisbecker.
C   Improving multifrontal methods by means of block low-rank representations.
C   SIAM Journal on Scientific Computing, 37(3):A1451--A1474 (2015).
C
C  W. M. Sid-Lakhdar, PhD Thesis from Universite de Lyon prepared at ENS Lyon,
C   Scaling the solution of large sparse linear systems using multifrontal
C   methods on hybrid shared-distributed memory architectures (2014).
C
C  P. Amestoy, J.-Y. L'Excellent, W. Sid-Lakhdar,
C   Characterizing asynchronous broadcast trees for multifrontal factorizations,
C   Workshop on Combinatorial Scientific Computing,
C   Lyon, France, July 21-23 (2014).
C
C  P. Amestoy, J.-Y. L'Excellent, F.-H. Rouet, W. Sid-Lakhdar,
C   Modeling 1D distributed-memory dense kernels for an asynchronous
C   multifrontal sparse solver, High-Performance Computing for Computational
C   Science, VECPAR 2014, Eugene, Oregon, USA, June 30 - July 3 (2014).
C
C  J.-Y. L'Excellent and W. M. Sid-Lakhdar,
C   Introduction of shared-memory parallelism in a distributed-memroy
C   multifrontal solver, Parallel Computing (40):3-4, pages 34-46 (2014).
C  
C  C. Weisbecker, PhD Thesis supported by EDF, INPT-IRIT,
C   Improving multifrontal solvers by means of algebraic block low-rank
C   representations (2013).
C
C  E. Agullo, P. Amestoy, A. Buttari, A. Guermouche,  G. Joslin, J.-Y.
C  L'Excellent, X. S. Li, A. Napov, F.-H. Rouet, M. Sid-Lakhdar, S. Wang, C.
C  Weisbecker, I. Yamazaki,
C   Recent Advances in Sparse Direct Solvers, 22nd Conference on Structural
C   Mechanics in Reactor Technology, San Francisco (2013).
C
C  P. Amestoy, A. Buttari, G. Joslin, J.-Y. L'Excellent, W. Sid-Lakhdar, C.
C   Weisbecker, M. Forzan, C. Pozza, R. Perrin, V. Pellissier,
C   Shared memory parallelism and low-rank approximation techniques applied
C   applied to direct solvers in FEM simulation in <i>IEEE Transactions on
C   Magnetics</i>, IEEE, Special issue, Compumag 2013 (2013).
C
C  L. Boucher, P. Amestoy, A, Buttari, F.-H. Rouet and M. Chauvin,
C   INTEGRAL/SPI data segmentation to retrieve sources intensity variations,
C   Astronomy & Astrophysics, Article 52, 20 pages,
C   http://dx.doi.org/10.1051/0004-6361/201219605 (2013).
C
C  F.-H. Rouet, PhD thesis from INPT, Toulouse, France,
C   Memory and Performance issues in parallel multifrontal factorization and
C   triangular solutions with sparse right-hand sides (2014).
C
C  J.-Y. L'Excellent, Habilitation thesis from ENS Lyon,
C   Multifrontal methods: Parallelism, Memory Usage and Numerical
C   Aspects (2012).
C
C  P. Amestoy, I.S. Duff, J.-Y. L'Excellent, Y. Robert, F.H. Rouet
C   and B. Ucar, On computing inverse entries of a sparse matrix in 
C   an out-of-core environment,
C   SIAM J. on Scientific Computing Vol. 34 N. 4, p. 1975-1999 (2012).
C
C  Amestoy, Buttari, Duff, Guermouche, L'Excellent, and Ucar
C   The Multifrontal Method, Encyclopedia of Parallel Computing,
C   editor David Padua, Springer (2011).
C
C  Amestoy, Buttari, Duff, Guermouche, L'Excellent, and Ucar
C   MUMPS, Encyclopedia of Parallel Computing,
C   editor David Padua, Springer (2011).
C
C  Agullo, Guermouche and L'Excellent, Reducing the {I/O} Volume in 
C   Sparse Out-of-core Multifrontal Methods}, SIAM SISC, Vol 31, Nb. 6, 
C   4774-4794 (2010).
C
C  Amestoy, Duff, Guermouche, Slavova, Analysis of the Solution Phase of a
C   Parallel Multifrontal Approach, Parallel Computing, Vol. 36, 3--15 (2010). 
C
C  Tzvetomila Slavova, PhD from INPT prepared at CERFACS,
C   Parallel triangular solution in the out-of-core multifrontal approach 
C   for solving large sparse linear systems, available as CERFACS 
C   Report TH/PA/09/59 (2009).
C
C  Agullo, Guermouche and L'Excellent, A Parallel Out-of-core Multifrontal 
C   Method: Storage of Factors on Disk and Analysis of Models for an 
C   Out-of-core Active Memory, Parallel Computing, Special Issue on Parallel 
C   Matrix Algorithms, Vol. 34, Nb 6-8, 296--317 (2008).
C
C  Emmanuel Agullo, PhD Thesis from LIP-Ecole Normale Superieure de Lyon,
C   On the Out-of-core Factorization of Large Sparse Matrices (Nov 2008).
C
C  Amestoy, Duff, Ruiz, and Ucar, "A parallel
C   matrix scaling algorithm".
C   In proceedings of VECPAR'08-International Meeting-High 
C   Performance Computing for Computational Science, (Jan 2008).
C
C  Guermouche and L'Excellent, Constructing Memory-minimizing Schedules 
C   for Multifrontal Methods, ACM TOMS, Vol. 32, Nb. 1, 17--32 (2006).
C
C  Amestoy, Guermouche, L'Excellent, and Pralet, 
C   Hybrid scheduling for the parallel solution
C   of linear systems.  Vol 32 (2), pp 136-156 (2006).
C
C  Stephane Pralet, PhD from INPT prepared at CERFACS,
C   Constrained orderings and scheduling for parallel sparse linear algebra,
C   available as CERFACS technical report, TH/PA/04/105, (Sept 2004).
C
C  Abdou Guermouche, PhD Thesis from LIP-Ecole Normale Superieure de Lyon,
C   Etude et optimisation du comportement memoire dans les methodes paralleles 
C   de factorisation de matrices creuses (2004).
C
C  Guermouche, L'Excellent and Utard, Impact of Reordering on the Memory of a
C   Multifrontal Solver, Parallel Computing, Vol. 29, Nb. 9, 1191--1218 (2003).
C
C  Amestoy, Duff, L'Excellent and Xiaoye S. Li, Impact of the Implementation 
C   of MPI Point-to-Point Communications on the Performance of Two General 
C   Sparse Solvers, Parallel Computing, Vol. 29, Nb 7, 833--847 (2003).
C
C  Amestoy, Duff, L'Excellent and Xiaoye S. Li, Analysis and Comparison of 
C   Two General Sparse Solvers for Distributed Memory Computers, ACM TOMS,
C   Vol. 27, Nb 4, 388--421 (2001).
C
C  Amestoy, Duff, Koster and  L'Excellent (2001), 
C   A fully asynchronous multifrontal solver using distributed dynamic
C   scheduling, SIAM Journal of Matrix Analysis and Applications,
C   Vol 23, No 1, pp 15-41 (2001).
C
C  Amestoy, Duff and  L'Excellent (2000),
C   Multifrontal parallel distributed symmetric and unsymmetric solvers,
C   Comput. Methods in Appl. Mech. Eng., 184,  501-520 (2000)
C
C  Amestoy, Duff and L'Excellent (1998),
C   Parallelisation de la factorisation LU de matrices
C   creuses non-symmetriques pour des architectures a memoire distribuee,
C   Calculateurs Paralleles Reseaux et systemes repartis, 
C   Vol 10(5), 509-520 (1998).
C
C  PARASOL Deliverable D2.1d (final report), 
C   DMUMPS Version 3.1, A MUltifrontal Massively Parallel Solver,
C   PARASOL project, EU ESPRIT IV LTR project 20160, (June 1999).
C
C  Jacko Koster, PhD from INPT prepared at CERFACS, On the parallel solution 
C   and the reordering of unsymmetric sparse linear systems (1997).
C
C  Vincent Espirat, Master's thesis from INPT(ENSEEIHT)-IRIT, Developpement 
C   d'une approche multifrontale pour machines a memoire distribuee et 
C   reseau heterogene de stations de travail (1996).
C
C  Patrick Amestoy, PhD from INPT prepared at CERFACS, Factorization of large 
C  sparse matrices based on a multifrontal approach in a multiprocessor 
C  environment, Available as CERFACS report TH/PA/91/2 (1991).
C
C============================================
C Argument lists and calling sequences
C============================================
C
C There is only one entry:
*
*  A Fortran 90 driver subroutine DMUMPS has been designed as a user
*   friendly interface to the multifrontal code. 
*   This driver, in addition to providing the 
*   normal functionality of a sparse solver, incorporates some
*   pre- and post-processing.
*   This driver enables the user to preprocess the matrix to obtain a 
*   maximum
*   transversal so that the permuted matrix has a zero-free diagonal,
*   to perform prescaling
*   of the original matrix (a choice of scaling strategies is provided),
*   to use iterative refinement to improve the solution,
*   and finally to perform error analysis.
* 
* The driver routine DMUMPS offers similar functionalities to other 
* sparse direct solvers, depending on the value of one of 
* its parameters (JOB).  The main ones are:
*
* (i)  JOB = -1 
C    initializes an instance of the package. This must be
C    called before any other call to the package concerning that instance.
C    It sets default values for other
C    components of DMUMPS_STRUC, which may then be altered before
C    subsequent calls to DMUMPS.
C    Note that three components of the structure must always be set by the
C    user (on all processors) before a call with JOB=-1. These are
C        id%COMM,
C        id%SYM, and
C        id%PAR.
C    CNTL, ICNTL can then be modified (see documentation) by the user.
C
* A value of JOB = -1 cannot be combined with other values for JOB
*
* (ii) JOB = 1 accepts the pattern of matrix A and chooses pivots
* from the diagonal using a selection criterion to
* preserve sparsity.  It uses the pattern of A + A^T 
* but ignores numerical values. It subsequently constructs subsidiary
* information for the actual factorization by a call with JOB_=_2.  
* An option exists for the user to
* input the pivot sequence, in which case only the necessary
* information for a JOB = 2 entry will be generated.  We call the JOB=1
* entry, the analysis phase.
C The following components of the structure define the centralized matrix 
C pattern and must be set by the user (on the host only) 
C before a call with JOB=1:
C   --- id%N, id%NZ (32-bit int) or id%NNZ (64-bit int),
C       id%IRN, and id%JCN
C       if the user wishes to input the structure of the
C       matrix in assembled format (ICNTL(5)=0, and ICNTL(18) $\neq$ 3),
C   --- id%ELTPTR, and id%ELTVAR
C       if the user wishes to input the matrix in elemental
C       format (ICNTL(5)=1).
C A distributed matrix format is also available (see documentation)
C
* (iii) JOB = 2 factorizes a matrix A using the information
* from a previous call with JOB = 1. The actual pivot sequence
* used may differ slightly from that of this earlier call if A is not
* diagonally dominant.
*
* (iv) JOB = 3 uses the factors generated by a JOB = 2 call to solve
* a system of equations A X = B or A^T X =B, where X and B are matrices
* that can be either dense or sparse.
* The sparsity of B is exploited to limit the number of operations 
* performed during solution. When only part of the solution is
* also needed (such as when computing selected entries of A^1) then
* further reduction of the number of operations is performed.
* This is particularly beneficial in the context of an 
* out-of-core factorization.
*
* (v) JOB = -2 frees all internal data allocated by the package.
*
* A call with JOB=3 must be preceded by a call with JOB=2,
* which in turn must be preceded by a call with JOB=1, which
* in turn must be preceded by a call with JOB=-1. Since the
* information passed from one call to the next is not
* corrupted by the second, several calls with JOB=2 for matrices
* with the same sparsity pattern but different values may follow
* a single call with JOB=1, and similarly several calls with JOB=3 
* can be used for different right-hand sides.
* Values 4, 5, 6 for the parameter JOB can invoke combinations
* of the three basic operations corresponding to JOB=1, 2 or 3.
* 
C
*********
C     --------------------------------------
C     Explicit interface needed for routines
C     using a target argument if they appear
C     in the same compilation unit.
C     --------------------------------------
      INTERFACE
      SUBROUTINE DMUMPS_CHECK_DENSE_RHS
     &(idRHS, idINFO, idN, idNRHS, idLRHS)
      DOUBLE PRECISION, DIMENSION(:), POINTER :: idRHS
      INTEGER, intent(in)    :: idN, idNRHS, idLRHS
      INTEGER, intent(inout) :: idINFO(:)
      END SUBROUTINE DMUMPS_CHECK_DENSE_RHS
      SUBROUTINE DMUMPS_ANA_DRIVER( id )
      USE DMUMPS_STRUC_DEF
      TYPE (DMUMPS_STRUC), TARGET :: id
      END SUBROUTINE DMUMPS_ANA_DRIVER
      SUBROUTINE DMUMPS_FAC_DRIVER( id )
      USE DMUMPS_STRUC_DEF
      TYPE (DMUMPS_STRUC), TARGET :: id
      END SUBROUTINE DMUMPS_FAC_DRIVER
      SUBROUTINE DMUMPS_SOLVE_DRIVER( id )
      USE DMUMPS_STRUC_DEF
      TYPE (DMUMPS_STRUC), TARGET :: id
      END SUBROUTINE DMUMPS_SOLVE_DRIVER
      SUBROUTINE DMUMPS_PRINT_ICNTL(id, LP)
      USE DMUMPS_STRUC_DEF
      TYPE (DMUMPS_STRUC), TARGET, INTENT(IN) :: id
      INTEGER  :: LP
      END SUBROUTINE DMUMPS_PRINT_ICNTL
      END INTERFACE
*  MPI
*  ===
      INCLUDE 'mpif.h'
      INTEGER MASTER
      PARAMETER ( MASTER = 0 )
      INTEGER IERR
*
*  ==========
*  Parameters
*  ==========
      TYPE (DMUMPS_STRUC) :: id
C
C  Main components of the structure are:
C  ------------------------------------
C
C   (see documentation for a complete description)
C
C  JOB is an INTEGER variable which must be set by the user to
C    characterize the factorization step.  Possible values of JOB
C    are given below
C
C     1   Analysis: Ordering and symbolic factorization steps.
C     2   Scaling and Numerical Factorization
C     3   Solve and Error analysis
C     4   Analysis followed by numerical factorization
C     5   Numerical factorization followed by Solving step
C     6   Analysis, Numerical factorization and Solve
C
C  N is an INTEGER variable which must be set by the user to the
C    order n of the matrix A.  It is not altered by the
C     subroutine.  
C
C  NZ / NNZ are INTEGER / INTEGER(8) variables which must be set by the user
C    to the number of entries being input, in case of centralized assembled
C    entry.  It is not altered by the subroutine. Only used if
C    ICNTL(5).eq.0 and ICNTL(18) .ne. 3 (assembled matrix entry,
C    or, at least, centralized matrix graph during analysis).
C
C    Restriction: NZ > 0 or NNZ > 0.
C    If NNZ is different from 0, NNZ is used. Otherwise, NZ is used.
C
C  NELT is an INTEGER variable which must be set by the user to the
C    number of elements being input.  It is not altered by the
C    subroutine. Only used if ICNTL(5).eq.1 (elemental matrix entry).
C    Restriction: NELT > 0.
C
C  IRN and JCN  are INTEGER  arrays of length [N]NZ.
C    IRN(k) and JCN(k), k=1..[N]NZ must be set on entry to hold 
C    the row and column indices respectively.
C    They are not altered by the subroutine except when ICNTL(6) = 1.
C    (in which case only the column indices are modified).
C    The arrays are only used if ICNTL(5).eq.0 (assembled entry)
C    or out-of-range.
C
C  ELTPTR is an INTEGER array of length NELT+1. 
C  ELTVAR is an INTEGER array of length ELTPTR(NELT+1)-1.
C    ELTPTR(I) points in ELTVAR to the first variable in the list of
C    variables that correspond to element I. ELTPTR(NELT+1) points
C    to the first unused location in ELTVAR.
C    The positions ELTVAR(I) .. ELTPTR(I+1)-1 contain the variables
C    for element I. No free space is allowed between variable lists.
C    ELTPTR/ELTVAR are not altered by the subroutine.
C    The arrays are only used if ICNTL(5).ne.0 (element entry).
C
C  A is a DOUBLE PRECISION array of length [N]NZ. 
C     The user must set A(k) to the value 
C     of the entry in row IRN(k) and column JCN(k) of the matrix.
C     It is not altered by the subroutine.
C     (Note that the matrix can also be provided in a distributed 
C      assembled input format)
C
C  RHS is a DOUBLE PRECISION array of length N that is only accessed when
C    JOB = 3, 5, or 6. On entry, RHS(i)
C     must hold the i th component of the right-hand side of the
C     equations being solved.
C     On exit, RHS(i) will hold the i th component of the
C     solution vector.  For other values of JOB, RHS is not accessed and
C     can be declared to have size one.
C     RHS should only be available on the host processor. If
C     it is associated on other processors, an error is raised.
C     (Note that the right-hand sides can also be provided in a 
C      sparse format).
C
C COLSCA, ROWSCA are DOUBLE PRECISION
C     arrays of length N that are used to hold
C     the values used to scale the columns and the rows
C     of the original matrix, respectively. 
C     These arrays need to be set by the user
C     only if ICNTL(8) is set to -1. If ICNTL(8)=0,
C     COLSCA and ROWSCA are not accessed and 
C     so can be declared to have size one.
C     For any other values of ICNTL(8),
C     the scaling arrays are computed before
C     numerical factorization.  The factors of the scaled matrix
C     diag(ROWSCA(i)) <A diag(COLSCA(i)) are computed.
C 
C  The workspace is automatically allocated by the package.
C  At the beginning of the numerical phase. If the user wants to increase
C   the allocated workspace (typically, numerical pivoting that leads to extra
C   storage, or previous call to MUMPS that failed because of 
C   a lack of allocated memory), 
C   we describe in the following how the user can modify the size 
C   of the workspace:
C    1/ The memory relaxation parameter
C       ICNTL(14) is designed to control the increase, with respect to the 
C       estimations performed during analysis, in the size of the workspace 
C       allocated during the numerical phase.
C    2/ The user can also provide 
C       a unique parameter,  ICNTL(23),  holding the maximum size of the total 
C       workspace (in Megabytes) that the package is allowed to use internally.
C       In this case we try as much as possible to follow the indication given
C       by the relaxation parameter (ICNTL(14)).
C
C   If ICNTL(23) is greater than 0 
C   then MUMPS automatically computes the size of the internal working arrays
C   such that the storage for all MUMPS internal data is equal to ICNTL(23).
C   The relaxation ICNTL(14) is first applied to
C   the internal integer working array and communication buffer sizes;
C   the remaining available space is given to the real/complex 
C   internal working arrays.
C   A lower bound of ICNTL(23) (if ICNTL(14) has not
C   been modified since the analysis) is given by INFOG(26).
C   
C   If ICNTL(23) is left to its default value 0 
C   then each processor will allocate workspace based on
C   the estimates computed during the analysis (INFO(17)
C   if ICNTL(14) has not been modified since analysis,
C   or larger if ICNTL(14) was increased). 
C   Note that these estimates are accurate in the sequential
C   version of {\tt MUMPS}, but that they can be inaccurate
C   in the parallel case. Therefore, in parallel, we recommend
C   to use ICNTL(23) and provide a value significantly larger
C   than INFOG(26).
C --------------------------------------------------------------------------------
C    
C CNTL is a DOUBLE PRECISION array of length 15
C  that contains control parameters and must be set by the user. Default
C  values for the components may be set by a call to DMUMPS(JOB=-1)
C  Details of the control parameters are given in DMUMPSID.
C
C ICNTL is an INTEGER array of length 60
C  that contains control parameters and must be set by the user. Default
C  values for the components may be set by a call to DMUMPS(JOB=-1)
C  Details of the control parameters are given in DMUMPSID.
C
C INFO is an INTEGER array of length 80 that need not be set by the
C  user.  On return from DMUMPS, a value of zero for INFO(1)
C  indicates that the subroutine has performed successfully. 
C  Details of the control parameters are given in DMUMPSID.
C 
C RINFO is a DOUBLE PRECISION  array of length 40 that need not be set by the
C  user.  This array supplies information on the execution of DMUMPS.
C  Details of the control parameters are given in DMUMPSID.
C
C
*
*
*   ====================
*    .. Error Return ..
*   ====================
*
C MUMPS uses the following mechanism to process errors that
C may occur during the parallel execution of the code. 
C If, during a call to MUMPS, an error occurs on a processor, 
C this processor informs all the other processors before they
C return from the call.
C In parts of the code where messages are sent asynchronously 
C (for example the factorization and solve phases), 
C the processor on which the error occurs sends a message 
C to the other processors with a specific error tag. 
C On the other hand, if the error occurs in a subroutine that
C does not use asynchronous communication, the processor propagates 
C the error to the other processors.
C On successful completion, a call to MUMPS will exit with the 
C parameter id%INFOG(1) set to zero.
C A negative value for id%INFOG(1) indicates that an 
C error has been detected on one of the processors.
C For example, if processor s returns with
C INFO(1)= -8 and INFO(2)=1000, then processor s ran out of integer 
C workspace during the factorization and the size of the workspace 
C should be increased by 1000 at least. 
C The other processors are informed about this error and return with
C INFO(1)=-1 (i.e., an error occurred on another processor) and 
C INFO(2)=s (i.e., the error occurred on processor s).
C If several processors raised an error, those processors do not overwrite 
C INFO(1), i.e., only processors that did not produce an error will set 
C INFO(1) to -1 and INFO(2) to the rank of the processor having the most 
C negative error code.
C
C The behaviour is slightly different for the global information
C parameters INFOG(1) and INFOG(2):
C in the previous example, all processors would return with
C INFOG(1)=-8 and INFOG(2)=1000.
C
C The possible error codes returned in INFO(1) (and INFOG(1))
C are fully described in the documentation.
C
C A positive value of INFO(1) is associated with a warning message 
C which  will be output on unit ICNTL(2) (see documentation).
C
C
C      .. Local variables ..
C
      INTEGER JOBMIN, JOBMAX, OLDJOB
!$    INTEGER NOMP, NOMPMIN, NOMPMAX
      INTEGER I, J, MP, LP, MPG, KEEP235SAVE, KEEP242SAVE,
     &        KEEP243SAVE, KEEP495SAVE, KEEP497SAVE
      INTEGER(8) :: I8
      LOGICAL LANA, LFACTO, LSOLVE, PROK, LPOK, FLAG, PROKG
      LOGICAL NOERRORBEFOREPERM
      LOGICAL UNS_PERM_DONE,I_AM_SLAVE
C     Saved communicator (pb of interference)
      INTEGER COMM_SAVE
C     Local copy of JOB
      INTEGER JOB
      CHARACTER(LEN=20) :: FROM_C_INTERFACE_STRING
      INTEGER, PARAMETER :: ICNTL18DIST_MIN = 1
      INTEGER, PARAMETER :: ICNTL18DIST_MAX = 3
      INTEGER, DIMENSION(:), ALLOCATABLE :: UNS_PERM_INV
C     TIMINGS
      DOUBLE PRECISION TIMEG, TIMETOTAL
      INTEGER(8) :: FILE_SIZE,STRUC_SIZE
      INTEGER:: ICNTL16_LOC
!$    INTEGER:: PREVIOUS_OMP_THREADS_NUM
      IF (id%JOB .EQ. -200) THEN
C       -----------------
C       QUICK TERMINATION:
C       -----------------
C       In case of very serious error that cannot be recovered
C       (example MPI crash, signal to kill a process, etc.),
C       a user error handler may want to call MUMPS with
C       JOB=-200 to just suppress the existing OOC files (if any)
C       before terminating the application and killing all
C       processes. No checks are performed in this case,
C       data is not freed and MPI should not be called
C       (since MPI might no longer work)
        CALL DMUMPS_OOC_CLEAN_FILES(id,IERR)
        RETURN
      ENDIF
      NOERRORBEFOREPERM = .FALSE.
      UNS_PERM_DONE = .FALSE.
      LANA  = .FALSE.
      LFACTO = .FALSE.
      LSOLVE = .FALSE.
      JOB  = id%JOB
C
C     Initialize error return codes to 0.
      id%INFO(1) = 0
      id%INFO(2) = 0
C     -----------------------------------
C     Check that MPI has been initialized
C     -----------------------------------
      CALL MPI_INITIALIZED( FLAG, IERR )
      IF ( .NOT. FLAG ) THEN
        id%INFO(1) = -23
        id%INFO(2) =   0
        WRITE(6,990)
 990  FORMAT(' Unrecoverable Error in DMUMPS initialization: ',
     &       ' MPI is not running.')
        RETURN               
      END IF
C     ---------------------------
C     Duplicate user communicator
C     to avoid communications not
C     related to DMUMPS
C     ---------------------------
       COMM_SAVE = id%COMM
       CALL MPI_COMM_DUP( COMM_SAVE, id%COMM, IERR )
C
C     Default setting for printing
      LP = 6
      MP = 0
      MPG = 0
      LPOK  = .TRUE.
      PROK  = .FALSE.
      PROKG = .FALSE.
      ICNTL16_LOC = 0
C     -------------------------
C     Check if value of JOB is
C     the same on all processes
C     -------------------------
      CALL MPI_ALLREDUCE(JOB,JOBMIN,1,MPI_INTEGER,MPI_MAX,
     &                   id%COMM,IERR)
      CALL MPI_ALLREDUCE(JOB,JOBMAX,1,MPI_INTEGER,MPI_MIN,
     &                   id%COMM,IERR)
      IF ( JOBMIN .NE. JOBMAX ) THEN
        id%INFO(1) = -3 
        id%INFO(2) = JOB
        GOTO 499
      END IF
C   
C     Check value of JOB and previous value of JOB
C
      IF ((JOB.LT.-3.OR.JOB.EQ.0.OR.JOB.GT.8)
     &    .AND. JOB.NE.9
     &  ) THEN
C       Out of range value
        id%INFO(1) = -3 
        id%INFO(2) = JOB
        GOTO 499
      END IF
      IF (JOB.NE.-1) THEN
C      Check the previous value of JOB
C      One should be able to test for old job value
C      Warning: non initialized value
       OLDJOB = id%KEEP( 40 ) + 456789
       IF (OLDJOB.NE.-1.AND.OLDJOB.NE.-2.AND.
     &    OLDJOB.NE.1.AND.OLDJOB.NE.2.AND.
     &    OLDJOB.NE.3) THEN
        id%INFO(1) = -3 
        id%INFO(2) = JOB
        GOTO 499
       END IF
      END IF
      IF(JOB.NE.-1) THEN
C     Always allow JOB=-1
         IF((JOB.GT.-2).AND.(id%KEEP(140).EQ.1)) then
C     If restore failed (keep(140)=1), the only allowed jobs (except -1)
C     are -2 or -3
            id%INFO(1) = -3 
            id%INFO(2) = JOB
            GOTO 499
         END IF
      ENDIF
C     ----------------------------------
C     Initialize, LANA, LFACTO, LSOLVE
C     LANA indicates if analysis must be performed
C     LFACTO indicates if factorization must be performed
C     LSOLVE indicates if solution must be performed
C     ----------------------------------
      IF ((JOB.EQ.1).OR.(JOB.EQ.4).OR.
     &    (JOB.EQ.6))               LANA  = .TRUE.
      IF ((JOB.EQ.2).OR.(JOB.EQ.4).OR.
     &    (JOB.EQ.5).OR.(JOB.EQ.6)) LFACTO = .TRUE.
      IF ((JOB.EQ.3).OR.(JOB.EQ.5).OR.
     &    (JOB.EQ.6))               LSOLVE = .TRUE.
      IF ( LANA .OR. LFACTO .OR. LSOLVE) THEN
C       Set some specific experimental KEEP entries,
C       Value defined on MASTER is the reference
        CALL MPI_BCAST( id%KEEP(370), 2, MPI_INTEGER, MASTER, id%COMM,
     &         IERR )
        IF (LANA) THEN
           CALL MPI_BCAST( id%KEEP(66), 1, MPI_INTEGER, MASTER, id%COMM,
     &         IERR )
          IF (id%KEEP(370) .EQ. 1) THEN     
            id%KEEP(77)=0
            id%KEEP(78)=0
            id%KEEP(375)=1           
          ENDIF
        ENDIF
        IF (id%KEEP(371) .EQ. 1) THEN 
          IF (LANA) THEN
            IF (id%KEEP(50) .EQ. 0 .AND. id%NSLAVES .GE. 32) THEN
              id%KEEP(376) = 1
            ENDIF
          ENDIF
          IF (LFACTO) THEN
            id%KEEP(351) = 1 
            id%KEEP(206) = 2 
          ENDIF
          IF (LSOLVE) THEN
            id%KEEP(350) = 2 
          ENDIF
        ENDIF
        IF (LANA) THEN
          IF (id%KEEP(66).NE.0) THEN
            id%KEEP(77)=0
            id%KEEP(78)=0
            id%KEEP(375)=1           
            IF ((id%KEEP(50).EQ.0) .AND. (id%NSLAVES.GT.1)) THEN
              id%KEEP(376) = 1
            ENDIF
            IF (id%KEEP(66).EQ.2) THEN
             id%KEEP(83) = 0
             id%KEEP(91) = 0
            ENDIF
          ENDIF
        ENDIF
      ENDIF ! LANA .OR. LFACTO .OR. LSOLVE
      IF (LANA) THEN
C       Decode NZ / NNZ into KEEP8(28) now
C       since we may want to print it.
C       NZ and NNZ only accessed at analysis
C       phase
        CALL MUMPS_GET_NNZ_INTERNAL( id%NNZ, id%NZ, id%KEEP8(28) )
      ENDIF
C     Also decode NNZ_loc in the same way
      IF (LANA .OR. LFACTO) THEN
C       NZ_loc/NNZ_loc may be accessed both analysis
C       and factorization:
        CALL MUMPS_GET_NNZ_INTERNAL( id%NNZ_loc, id%NZ_loc,
     &                               id%KEEP8(29))
      ENDIF
C
      IF (JOB.EQ.-2.OR.JOB.EQ.1.OR.JOB.EQ.2.OR.JOB.EQ.3.OR.
     &    JOB.EQ.4.OR.JOB.EQ.5.OR.JOB.EQ.6
     &     .OR.JOB.EQ.7.OR.JOB.EQ.8.OR.JOB.EQ.-3
     *     .OR. JOB.EQ.9 ! LATER .OR.JOB.EQ.10
     &     ) THEN
C       Value of JOB is correct and differnet from -1
C       ICNTL should have been initialized and can be used
        LP      = id%ICNTL(1)
        MP      = id%ICNTL(2)
        MPG     = id%ICNTL(3)
        LPOK    = ((LP.GT.0).AND.(id%ICNTL(4).GE.1))
        PROK    = ((MP.GT.0).AND.(id%ICNTL(4).GE.2))
        PROKG   = ( MPG .GT. 0 .and. id%MYID .eq. MASTER )
        PROKG   = (PROKG.AND.(id%ICNTL(4).GE.2))
        IF (id%KEEP(500).EQ.1) THEN
          FROM_C_INTERFACE_STRING=" from C interface"
        ELSE
          FROM_C_INTERFACE_STRING=" "
        ENDIF
C       -----------------------------
C       Set the number of threads if required.
C       -----------------------------
        ICNTL16_LOC = id%ICNTL(16)
        CALL MPI_BCAST( ICNTL16_LOC, 1, MPI_INTEGER, MASTER, id%COMM,
     &          IERR )
!$      IF (ICNTL16_LOC .GT. 0) THEN
!$        PREVIOUS_OMP_THREADS_NUM = omp_get_max_threads()
#if defined(WORKAROUNDINTELILP64OPENMPLIMITATION)
!$        CALL omp_set_num_threads(int(ICNTL16_LOC,4))
#else
!$        CALL omp_set_num_threads(ICNTL16_LOC)
#endif
!$      ENDIF
C       -----------------------------
C       Check if number of threads is
C       the same on all processes
C       -----------------------------
!$      NOMP = OMP_GET_MAX_THREADS()
!$      CALL MPI_ALLREDUCE(NOMP,NOMPMIN,1,MPI_INTEGER,MPI_MIN,
!$   &                     id%COMM,IERR)
!$      CALL MPI_ALLREDUCE(NOMP,NOMPMAX,1,MPI_INTEGER,MPI_MAX,
!$   &                     id%COMM,IERR)
        id%KEEP(249) = 1
!$      id%KEEP(249) = OMP_GET_MAX_THREADS()
        IF (PROKG) THEN
C          Print basic information on MUMPS call
           IF (JOB .EQ. -2
     &       .OR. JOB .EQ. 8
     &       ) THEN
C            N, NELT, NNZ not meaningful
C            or not defined yet (JOB=8).
             WRITE(MPG,'(/A,A,A,A,I4)') 
     &               'Entering DMUMPS ',
     &               trim(adjustl(id%VERSION_NUMBER)),
     &               trim(FROM_C_INTERFACE_STRING),
     &               ' with JOB =', JOB
           ELSE IF (id%ICNTL(5) .NE. 1) THEN
C            Assembled format
             IF (id%ICNTL(18) .EQ. 0
     &            ) THEN
                 WRITE(MPG,'(/A,A,A,A,I4,I12,I15)') 
     &                 'Entering DMUMPS ',
     &                 trim(adjustl(id%VERSION_NUMBER)),
     &                 trim(FROM_C_INTERFACE_STRING),
     &                 ' with JOB, N, NNZ =', JOB,id%N,id%KEEP8(28)
             ELSE
                 WRITE(MPG,'(/A,A,A,A,I4,I12)') 
     &                 'Entering DMUMPS ',
     &                 trim(adjustl(id%VERSION_NUMBER)),
     &                 trim(FROM_C_INTERFACE_STRING),
     &                 ' with JOB, N =', JOB,id%N
             ENDIF
           ELSE
C            Elemental format
             WRITE(MPG,'(/A,A,A,A,I4,I12,I15)') 
     &                'Entering DMUMPS ',
     &                trim(adjustl(id%VERSION_NUMBER)),
     &                trim(FROM_C_INTERFACE_STRING),
     &                ' driver with JOB, N, NELT =', JOB,id%N,id%NELT
           ENDIF
C          MPI and OpenMP information
!$         IF (.TRUE.) THEN
!$           WRITE(MPG, '(A,I6,A,I6)') '      executing #MPI = ',
!$   &                 id%NPROCS, ' and #OMP = ', NOMP
!$           IF ( NOMPMIN .NE. NOMPMAX ) THEN
!$             WRITE(MPG, '(A,I4,A,I4,A)')
!$   &  '      WARNING detected: different number of threads (max ',
!$   &         NOMPMAX, ', min ', NOMPMIN, ')'
!$           END IF
!$         ELSE
             WRITE(MPG, '(A,I6,A)')    '      executing #MPI = ',
     &                 id%NPROCS, ', without OMP'
!$         ENDIF
          IF (JOB.GE.1 .AND. JOB.LE.6) THEN
            WRITE(MPG, '(A)')
          ENDIF
        ENDIF
      END IF
C
C----------------------------------------------------------------
C
C     JOB = -1 : START INITIALIZATION PHASE
C                (NEW INSTANCE)
C
C     JOB = -2 : TERMINATE AN INSTANCE
C----------------------------------------------------------------
C
      IF ( JOB .EQ. -1 ) THEN
C
C       ------------------------------------------
C       Check that we have called (JOB=-2), ie
C       that the previous JOB is not 1 2 or 3,
C       before calling the initialization routine.
C       --------------------------------------------
        id%INFO(1)=0
        id%INFO(2)=0
        OLDJOB = id%KEEP( 40 ) + 456789
        IF ( OLDJOB .EQ. 1 .OR.
     &       OLDJOB .EQ. 2 .OR.
     &       OLDJOB .EQ. 3  ) THEN
          IF ( id%N > 0 ) THEN
           id%INFO(1)=-3
           id%INFO(2)=JOB
          ENDIF
        ENDIF
C       Initialize id%MYID now because it is
C       required by MUMPS_PROPINFO. id%MYID
C       used to be initialized inside DMUMPS_INI_DRIVER,
C       leading to an uninitialized access here.
        CALL MPI_COMM_RANK(id%COMM, id%MYID, IERR)
        CALL MUMPS_PROPINFO( id%ICNTL(1),
     &                       id%INFO(1),
     &                       id%COMM, id%MYID )
        IF ( id%INFO(1) .LT. 0 ) THEN
C
C         If there was an error, then initialization
C         was already called and we can rely on the null
C         or non null value of the pointers related to OOC
C         stuff.
C         We use DMUMPS_CLEAN_OOC_DATA that should work even
C         on the master. Note that KEEP(201) was also
C         initialized in a previous call to Mumps.
C
C         If DMUMPS_END_DRIVER or DMUMPS_FAC_DRIVER is called after
C         this error, then DMUMPS_CLEAN_OOC_DATA will be called
C         a second time, though.
C
           IF (id%KEEP(201).GT.0) THEN
             CALL DMUMPS_CLEAN_OOC_DATA(id, IERR)
           ENDIF
           GOTO 499
        ENDIF
C       ----------------------------------------
C       Initialization DMUMPS_INI_DRIVER 
C       ----------------------------------------
C       - Default values for ICNTL, KEEP,KEEP8, CNTL
C       - Attach emission buffer for buffered Send
C       - Nullify pointers in the structure
C       - Get rank and size of the communicator
C       ----------------------------------------
        CALL DMUMPS_INI_DRIVER( id )
        IF ( id%INFO(1) .LT. 0 ) GOTO 499       
        GOTO 500
      END IF
      IF ( JOB .EQ. -2 ) THEN
C       -------------------------------------
C       Deallocation of the instance id
C       -------------------------------------
        id%KEEP(40)= -2 - 456789
        CALL DMUMPS_END_DRIVER( id )
        GOTO 500
      END IF
C
C     TIMINGS: for JOBS different from -1 and -2,
C     we measure TIMETOTAL:
C
      IF (id%MYID.EQ.MASTER) THEN
        id%DKEEP(70)=0.0D0
        CALL MUMPS_SECDEB(TIMETOTAL)
      ENDIF
C
C----------------------------------------------------------------
C
C     JOB = 7 : SAVE THE INSTANCE
C
C     JOB = 8 : RESTORE THE INSTANCE
C----------------------------------------------------------------
C
      IF ( JOB .EQ. 7 .OR. JOB .EQ. 8 ) THEN
         IF( JOB.EQ.8 .AND. OLDJOB.NE.-1) THEN
            id%INFO(1) = -3 
            id%INFO(2) = JOB
            GOTO 499
         END IF
         IF (id%MYID.EQ.MASTER) THEN
C       -----------------------------
C       Check incompatibility between
C       par (=0) and nprocs (=1)
C       -----------------------------
            IF ( (id%KEEP(46).EQ.0).AND.(id%NPROCS.LE.1) ) 
     &           THEN
               id%INFO(1) = -21
               id%INFO(2) = id%NPROCS
            ENDIF
         ENDIF
         CALL MUMPS_PROPINFO( id%ICNTL(1),
     &        id%INFO(1),
     &        id%COMM, id%MYID )
         IF ( id%INFO(1) .LT. 0 ) GOTO 499
         IF ( JOB .EQ. 7 ) THEN
            IF (id%MYID.EQ.MASTER) THEN
             CALL MUMPS_SECDEB(TIMEG)
            ENDIF
            CALL DMUMPS_SAVE( id )
            IF (id%MYID.EQ.MASTER) THEN
               CALL MUMPS_SECFIN(TIMEG)
               IF (PROKG) THEN
                  WRITE( MPG,'(/A,F12.4)')
     &                 ' Elapsed time in save structure driver= ', TIMEG
               END IF
            ENDIF
         ELSE
            IF (id%MYID.EQ.MASTER) THEN
             CALL MUMPS_SECDEB(TIMEG)
            ENDIF
            CALL DMUMPS_RESTORE( id )
            IF (id%MYID.EQ.MASTER) THEN
               CALL MUMPS_SECFIN(TIMEG)
               IF (PROKG) THEN
                  WRITE( MPG,'(/A,F12.4)')
     &                 ' Elapsed time in restore structure driver= '
     &                 , TIMEG
               ENDIF
            END IF
         ENDIF
         IF ( id%INFO(1) .LT. 0 ) GOTO 499
         GOTO 500
      ENDIF
C
C----------------------------------------------------------------
C
C     JOB = -3 : REMOVE SAVED INSTANCE
C
C----------------------------------------------------------------
C
      IF (JOB .EQ. -3) THEN
        CALL DMUMPS_REMOVE_SAVED(id)
        IF ( id%INFO(1) .LT. 0 ) GOTO 499
        GOTO 500
      ENDIF
      IF (JOB.EQ.9) THEN
C       Check that factorization was performed
        IF ( OLDJOB .LT. 2 ) THEN
          id%INFO(1)=-3
          id%INFO(2)=JOB
        ELSE
          CALL DMUMPS_SOL_INIT_IRHS_loc(id)
        ENDIF
        IF ( id%INFO(1) .LT. 0 ) GOTO 499
        GOTO 500
      ENDIF
C
C----------------------------------------------------------------
C
C     MAIN DRIVER
C     OTHER VALUES OF JOB : 1 to 6
C
C----------------------------------------------------------------
      CALL MUMPS_MEMORY_SET_DATA_SIZES()
      IF (id%MYID.EQ.MASTER) THEN
C       -----------------------------
C       Check incompatibility between
C       par (=0) and nprocs (=1)
C       -----------------------------
         IF ( (id%KEEP(46).EQ.0).AND.(id%NPROCS.LE.1) ) 
     &        THEN
            id%INFO(1) = -21
            id%INFO(2) = id%NPROCS
         ENDIF
      END IF
C
C     Propagate possible error to all nodes
      CALL MUMPS_PROPINFO( id%ICNTL(1),
     &                    id%INFO(1),
     &                    id%COMM, id%MYID )
      IF ( id%INFO(1) .LT. 0 ) GOTO 499
C
C     Print ICNTL and KEEP
C
      IF (PROK) CALL DMUMPS_PRINT_ICNTL(id, MP)
C-----------------------------------------------------------------------
C
C           CHECK SEQUENCE
C
C-----------------------------------------------------------------------
      IF ( LANA ) THEN
        IF ( PROKG .AND. OLDJOB .EQ. -1 ) THEN
C         Print compilation options at first call to analysis
          CALL  MUMPS_PRINT_IF_DEFINED(MPG)
        ENDIF
C
C       User wants to perform analysis. Previous value of
C       JOB must be -1, 1, 2 or 3.
C
        IF ( OLDJOB .EQ. 0 .OR. OLDJOB .GT. 3 .OR. OLDJOB .LT. -1 ) THEN
          id%INFO(1) = -3
          id%INFO(2) = JOB
          GOTO 499
        END IF
        IF ( OLDJOB .GE. 2 ) THEN
C         -----------------------------------------
C         Previous step was factorization or solve.
C         As analysis is now performed, deallocate
C         at least some big arrays from facto.
C         -----------------------------------------
          IF (associated(id%IS)) THEN
            DEALLOCATE  (id%IS)
            NULLIFY     (id%IS)
          END IF
          IF (associated(id%S)) THEN
            DEALLOCATE  (id%S)
            NULLIFY     (id%S)
          END IF
        END IF   
      END IF
      IF ( LFACTO ) THEN
C        ------------------------------------
C        User wants to perform factorization.
C        Analysis must have been performed.
C        ------------------------------------
         IF ( OLDJOB .LT. 1 .and. .NOT. LANA ) THEN
            id%INFO(1) = -3
            id%INFO(2) = JOB
            GOTO 499
         END IF
      END IF
      IF ( LSOLVE ) THEN
C        -------------------------------
C        User wants to perform solve.
C        Facto must have been performed.
C        -------------------------------
         IF ( OLDJOB .LT. 2 .AND. .NOT. LFACTO ) THEN
            id%INFO(1) = -3
            id%INFO(2) = JOB
            GOTO 499
         END IF
      END IF
C     ------------------------------------------
C     Permute JCN on entry to JOB if no analysis
C     to be performed and IRN/JCN are needed.
C     (facto: arrowheads + solve: iterative
C      refinement and error analysis)
C     ------------------------------------------
#if ! defined (LARGEMATRICES)
      NOERRORBEFOREPERM =.TRUE.
      UNS_PERM_DONE=.FALSE.
      IF (id%MYID .eq. MASTER .AND. id%KEEP(23) .NE. 0) THEN
        IF ( id%JOB .EQ. 2 .OR. id%JOB .EQ. 5 .OR.
     &       (id%JOB .EQ. 3 .AND. (id%ICNTL(10) .NE.0 .OR.
     &        id%ICNTL(11).NE. 0))) THEN
          UNS_PERM_DONE = .TRUE.
          ALLOCATE(UNS_PERM_INV(id%N),stat=IERR)
          IF (IERR .GT. 0) THEN
C             --------------------------------
C             Exit with an error.
C             We are not able to permute
C             JCN correctly after a MAX-TRANS
C             permutation resulting from a
C             previous call to DMUMPS.
C             --------------------------------
              id%INFO(1)=-13
              id%INFO(2)=id%N
              IF (LPOK) WRITE(LP,99993)
              GOTO 510
          ENDIF
          DO I = 1, id%N
            UNS_PERM_INV(id%UNS_PERM(I))=I
          END DO
          DO I8 = 1_8, id%KEEP8(28)
            J = id%JCN(I8)
C           -- skip out-of range (that are ignored in ANA_O)
            IF (J.LE.0.OR.J.GT.id%N) CYCLE
            id%JCN(I8)=UNS_PERM_INV(J)
          END DO
          DEALLOCATE(UNS_PERM_INV)
        END IF
      END IF
#endif
C
C       Propagate possible error
        CALL MUMPS_PROPINFO( id%ICNTL(1),
     &                    id%INFO(1),
     &                    id%COMM, id%MYID )
        IF ( id%INFO( 1 ) .LT. 0 ) GO TO 499
*
*********
* MaxTrans-Analysis-Distri, Scale-Arrowhead-factorize, and
* Solve-IR-Error_Analysis (depending on the value of JOB)
*********
*
C
      IF ( LANA ) THEN
C-----------------------------------------------------
C-
C-       ANALYSIS : Max-Trans, Analysis, Distribution
C-
C-----------------------------------------------------
C
C        Few checks + allocations
C
C        IS : will be allocated on the slaves later
C        PROCNODE : on the master only,
C             because slave does not know N yet.
C             Will be allocated in analysis for the slave.
C
C        For assembled entry: 
C        IRN, JCN : check that they have been allocated by the
C             user on the master, and if their size is adequate
C
C        For element entry:
C        ELTPTR, ELTVAR : check that they have been allocated by the
C             user on the master, and if their size is adequate
C       ----------------------------
C       Reset KEEP(40) to -1 for the
C       case where an error occurs
C       ----------------------------
        id%KEEP(40)=-1 -456789
C
        IF (id%MYID.EQ.MASTER) THEN
C         Check N, [N]NZ, NELT
          IF ((id%N.LE.0).OR.((id%N+id%N+id%N)/3.NE.id%N)) THEN
             id%INFO(1) = -16
             id%INFO(2) = id%N
             GOTO 100
          END IF
          IF (id%ICNTL(5).NE.1) THEN
C           Assembled input
            IF (id%ICNTL(18) .LT. 1 .OR. id%ICNTL(18) .GT. 3) THEN
C             Centralized input
              IF (id%KEEP8(28) .LE. 0_8) THEN
                id%INFO(1) = -2
                CALL MUMPS_SET_IERROR(id%KEEP8(28), id%INFO(2))
                GOTO 100
              ENDIF
            ENDIF
          ELSE
C           Element entry: check NELT on the master
            IF (id%NELT .LE. 0) THEN
              id%INFO(1) = -24
              id%INFO(2) = id%NELT
              GOTO 100
            ENDIF
          ENDIF
C     -- initialize values of respectively
C     icntl(6), (7) and (12) to not done/chosen
          id%INFOG(7) = -9999
          id%INFOG(23) = 0
          id%INFOG(24) = 1
C         ---------------------------------------
C         Element entry: allocate ELTPROC(1:NELT)
C         ---------------------------------------
          IF ( id%ICNTL(5) .EQ. 1 ) THEN ! Elemental matrix
            IF ( associated( id%ELTPROC ) )
     &           DEALLOCATE( id%ELTPROC )
            ALLOCATE( id%ELTPROC(id%NELT), stat=IERR )
            IF (IERR.gt.0) THEN
              id%INFO(1) = -7
              id%INFO(2) = id%NELT
              IF ( LPOK ) WRITE(LP,'(A)')
     &          'Problem in allocating work array ELTPROC'
              GOTO 100
            END IF
          END IF
C         ---------------------------------------------------
C         Assembled centralized entry: check input parameters
C         IRN/JCN
C         Element entry: check input parameters ELTPTR/ELTVAR
C         ---------------------------------------------------
          IF ( id%ICNTL(5) .NE. 1 ) THEN ! Assembled matrix
            id%KEEP8(30)=0_8
            IF ( id%ICNTL(18) .LT. ICNTL18DIST_MIN
     &           .OR. id%ICNTL(18) .GT. ICNTL18DIST_MAX ) THEN
              IF ( .not. associated( id%IRN ) ) THEN
                id%INFO(1) = -22
                id%INFO(2) = 1
#if defined(MUMPS_F2003)
              ELSE IF ( size( id%IRN, KIND=8 ) < id%KEEP8(28) ) THEN
#else
C             size with kind=8 output not available before f2002. One can
C             still check that if NZ can be stored in a 32-bit integer,
C             the 32-bit size(id%IRN) is large enough
              ELSE IF ( id%KEEP8(28) .LE. int(huge(id%NZ),8) .AND.
     &                  size(id%IRN) < int(id%KEEP8(28)) ) THEN
#endif
                id%INFO(1) = -22
                id%INFO(2) = 1
              ELSE IF ( .not. associated( id%JCN ) ) THEN
                id%INFO(1) = -22
                id%INFO(2) = 2
#if defined(MUMPS_F2003)
              ELSE IF ( size( id%JCN, KIND=8 ) < id%KEEP8(28) ) THEN
#else
C             Same as for IRN above
              ELSE IF ( id%KEEP8(28) .LE. int(huge(id%NZ),8) .AND.
     &                  size(id%JCN) < int(id%KEEP8(28)) ) THEN
#endif
                id%INFO(1) = -22
                id%INFO(2) = 2
              END IF
            END IF
            IF ( id%INFO( 1 ) .eq. -22 ) THEN
              IF ( LPOK ) WRITE(LP,'(A)')
     &           'Error in analysis: IRN/JCN badly allocated.'
            END IF
          ELSE
            IF ( .not. associated( id%ELTPTR ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 1
            ELSE IF ( size( id%ELTPTR ) < id%NELT+1 ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 1
            ELSE IF ( .not. associated( id%ELTVAR ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 2
            ELSE 
              id%LELTVAR = id%ELTPTR( id%NELT+1 ) - 1
              IF ( size( id%ELTVAR ) < id%LELTVAR ) THEN 
                id%INFO(1) = -22
                id%INFO(2) = 2
              ELSE
C               If no error, we compute KEEP8(30) (formerly NA_ELT),
C               required for DMUMPS_MAX_MEM already in analysis, and
C               then later during facto to check the size of A_ELT
                id%KEEP8(30) = 0_8
                IF ( id%KEEP(50) .EQ. 0 ) THEN
C                 Unsymmetric elements (but symmetric structure)
                  DO I = 1,id%NELT
                    J = id%ELTPTR(I+1) - id%ELTPTR(I)
                    id%KEEP8(30) = id%KEEP8(30) + int(J,8) * int(J,8)
                  ENDDO
                ELSE
C                 Symmetric elements
                  DO I = 1,id%NELT
                    J = id%ELTPTR(I+1) - id%ELTPTR(I)
                    id%KEEP8(30) = id%KEEP8(30) +
     &                             (int(J,8) *int(J+1,8))/2_8
                  ENDDO
                ENDIF
              ENDIF
            END IF
            IF ( id%INFO( 1 ) .eq. -22 ) THEN
              IF ( LPOK ) WRITE(LP,'(A)')
     &           'Error in analysis: ELTPTR/ELTVAR badly allocated.'
            END IF
          ENDIF
 100      CONTINUE
        END IF
C
C       Propagate possible error
        CALL MUMPS_PROPINFO( id%ICNTL(1),
     &                    id%INFO(1),
     &                    id%COMM, id%MYID )
        IF ( id%INFO( 1 ) .LT. 0 ) GO TO 499
C       -----------------------------------------
C       Call analysis procedure DMUMPS_ANA_DRIVER
C       -----------------------------------------
        IF (id%MYID .eq. MASTER) THEN
          id%DKEEP(71)=0.0D0
          CALL MUMPS_SECDEB(TIMEG)
        END IF 
C       -------------------------------------------------
C       Set scaling option for analysis in KEEP(52)
C       (ICNTL(8) only defined on host at analysis phase)
C       -------------------------------------------------
        IF (id%MYID.EQ.MASTER) THEN
C{
         id%KEEP(52) = id%ICNTL(8)
C        Out-of-range values => automatic choice
         IF ( id%KEEP(52) .GT. 8 .OR. id%KEEP(52).LT.-2)
     &       id%KEEP(52) = 77
         IF ( id%KEEP(52) .EQ. 2 .OR. id%KEEP(52).EQ.5 
     &       .OR. id%KEEP(52) .EQ. 6 )
     &       id%KEEP(52) = 77
         IF ((id%KEEP(52).EQ.77).AND.(id%KEEP(50).EQ.1)) THEN
          ! for SPD matrices default is no scaling
           id%KEEP(52) = 0
         ENDIF
         IF ( id%KEEP(52).EQ.77 .OR. id%KEEP(52).LE.-2) THEN 
C          -- suppress scaling computed during analysis 
C          -- if centralized matrix is not associated
          IF (.not.associated(id%A)) id%KEEP(52) = 0
         ENDIF
C        deactivate analysis scaling if scaling given
         IF(id%KEEP(52) .EQ. -1) id%KEEP(52) = 0
C
C        deactivate analysis scaling if 
C        permutation to zero-free diagonal not requested
         IF (id%ICNTL(6).EQ.0) id%KEEP(52) = 0
C        deactivate analysis scaling for SPD matrices
         IF (id%KEEP(50).EQ.1) id%KEEP(52) = 0
C
         IF (id%KEEP(52).EQ.-2) THEN
C        deallocate scalings in case of ordering allocated/computed 
C        during analysis. This is needed because in case of 
C        KEEP(52)=-2 then one cannot be sure that 
C        scaling will be effectivly computed during analysis
C        Thus to test if scaling was effectively allocated/computed
C        during analysis after DMUMPS_ANA_DRIVER one must
C        be sure that scaling arrays are nullified.
           IF ( associated(id%COLSCA)) THEN
                DEALLOCATE( id%COLSCA )
                NULLIFY(id%COLSCA)
           ENDIF
           IF ( associated(id%ROWSCA)) THEN
                DEALLOCATE( id%ROWSCA )
                NULLIFY(id%ROWSCA)
           ENDIF
         ENDIF
C     
C}
        ENDIF
C
C       ANALYSIS PHASE:
        CALL DMUMPS_ANA_DRIVER( id )
C
C Check and save scaling option in INFOG(33)
        IF (id%MYID .eq. MASTER) THEN
C{
          IF (id%KEEP(52).EQ.0) id%INFOG(33)=id%ICNTL(8)
          IF (id%KEEP(52).EQ.-2) THEN
C         Scaling should have been computed during
           IF (.not.associated(id%COLSCA).OR.
     &         .not.associated(id%ROWSCA)
     &        ) THEN
C            scaling was not computed reset KEEP(52)
C            the user can then decide during factorization
C            to activate scaling
             id%KEEP(52) =0
             id%INFOG(33)=0
             IF ( MPG .GT. 0 ) THEN
                WRITE(MPG,'(A)') 
     &          ' Warning;  scaling was not computed during analysis'
             ENDIF
             IF ( associated(id%COLSCA)) THEN
                  DEALLOCATE( id%COLSCA )
                  NULLIFY(id%COLSCA)
             ENDIF
             IF ( associated(id%ROWSCA)) THEN
                  DEALLOCATE( id%ROWSCA )
                  NULLIFY(id%ROWSCA)
             ENDIF
           ENDIF
          ENDIF
          IF (id%KEEP(52) .NE. 0) THEN
            id%INFOG(33)=id%KEEP(52)
          ENDIF
C}
        ENDIF
C       return value of ICNTL(12) effectively used
C       that was saved on the master in KEEP(95)
        IF (id%MYID .eq. MASTER) id%INFOG(24)=id%KEEP(95)
C       TIMINGS:
        IF (id%MYID .eq. MASTER) THEN
          CALL MUMPS_SECFIN(TIMEG)
          id%DKEEP(71) = TIMEG
        ENDIF
        IF (PROKG) THEN
          WRITE( MPG,'(/A,F12.4)')
     &         ' Elapsed time in analysis driver= ', TIMEG
        END IF 
C       -----------------------
C     Return in case of error
C     -----------------------
        IF ( id%INFO( 1 ) .LT. 0 ) GO TO 499
        id%KEEP(40) = 1 -456789
      END IF
C
C-------------------------------------------------------
C-
C      
C      BEGIN FACTORIZATION PHASE
C
C-
C-------------------------------------------------------
      IF ( LFACTO ) THEN
         IF (id%MYID .eq. MASTER) THEN
            id%DKEEP(91)=0.0D0
            CALL MUMPS_SECDEB(TIMEG)
         END IF 
C        ----------------------
C        Reset KEEP(40) to 1 in
C        case of error in facto
C        ----------------------
         id%KEEP(40) = 1 - 456789
C
C-------------------------------------------------------
C-
C-      CHECKS, SCALING, ARROWHEAD + FACTORIZATION PHASE
C-
C-------------------------------------------------------
C
        IF ( id%MYID .EQ. MASTER ) THEN
C         -------------------------
C         Check if Schur complement
C         is allocated.
C         -------------------------
          IF (id%KEEP(60).EQ.1) THEN
             IF ( associated( id%SCHUR_CINTERFACE)) THEN
C              Called from C interface...
C              Since id%SCHUR_CINTERFACE is of size 1,
C              instruction below which causes bound check
C              errors should be avoided. We cheat by first
C              setting a static pointer with a routine with
C              implicit interface, and then copying this pointer
C              into id%SCHUR.
               CALL DMUMPS_SET_TMP_PTR(id%SCHUR_CINTERFACE(1),
     &         int(id%SIZE_SCHUR,8)*int(id%SIZE_SCHUR,8))
               CALL DMUMPS_GET_TMP_PTR(id%SCHUR)
               NULLIFY(id%SCHUR_CINTERFACE)
             ENDIF
             IF ( .NOT. associated (id%SCHUR)) THEN
              IF (LP.GT.0) 
     &        write(LP,'(A)') 
     &                      ' SCHUR not associated'
              id%INFO(1)=-22
              id%INFO(2)=9
             ELSE IF ( size(id%SCHUR) .LT.
     &                id%SIZE_SCHUR * id%SIZE_SCHUR ) THEN
                IF (LP.GT.0) 
     &          write(LP,'(A)') 
     &                ' SCHUR allocated but too small' 
                id%INFO(1)=-22
                id%INFO(2)=9
             END IF
          END IF
C     ------------------------------------------------------------
C     Assembled entry: check input parameterd IRN,JCN,A
C     Element entry: check input parameters ELTPTR,ELTVAR,A_ELT
C     ------------------------------------------------------------
          IF ( id%KEEP(54) .EQ. 0 ) THEN
             IF ( id%KEEP(55).eq.0 ) THEN
C     Assembled entry
                IF ( .not. associated( id%IRN ) ) THEN
                   id%INFO(1) = -22
                   id%INFO(2) = 1
#if defined(MUMPS_F2003)
                ELSE IF ( size( id%IRN, KIND=8 ) < id%KEEP8(28) ) THEN
#else
C     size with kind=8 output not available. One can still
C     check that if NZ can be stored in a 32-bit integer,
C     the 32-bit size(id%IRN) (which we then assume not
C     to overflow...) is large enough
                ELSE IF ( id%KEEP8(28) .LE. int(huge(id%NZ),8) .AND.
     &                  size(id%IRN) < int(id%KEEP8(28)) ) THEN
#endif
                   id%INFO(1) = -22
                   id%INFO(2) = 1
                ELSE IF ( .not. associated( id%JCN ) ) THEN
                   id%INFO(1) = -22
                   id%INFO(2) = 2
#if defined(MUMPS_F2003)
                ELSE IF ( size( id%JCN, KIND=8 ) < id%KEEP8(28) ) THEN
#else
C     Same as for IRN above
                ELSE IF ( id%KEEP8(28) .LE. int(huge(id%NZ),8) .AND.
     &                  size(id%JCN) < int(id%KEEP8(28)) ) THEN
#endif
                   id%INFO(1) = -22
                   id%INFO(2) = 2
                ELSEIF ( .not. associated( id%A ) ) THEN
                   id%INFO( 1 ) = -22
                   id%INFO( 2 ) = 4
#if defined(MUMPS_F2003)
                ELSE IF ( size( id%A, KIND=8 ) < id%KEEP8(28) ) THEN
#else
C     Same as for IRN/JCN above
                ELSE IF ( id%KEEP8(28) .LE. int(huge(id%NZ),8) .AND.
     &                  size( id%A ) < int(id%KEEP8(28)) ) THEN
#endif
                   id%INFO( 1 ) = -22
                   id%INFO( 2 ) = 4
                END IF
             ELSE
C     Element entry                
                IF ( .not. associated( id%ELTPTR ) ) THEN
                   id%INFO(1) = -22
                   id%INFO(2) = 1
                ELSE IF ( size( id%ELTPTR ) < id%NELT+1 ) THEN
                   id%INFO(1) = -22
                   id%INFO(2) = 1
                ELSE IF ( .not. associated( id%ELTVAR ) ) THEN
                   id%INFO(1) = -22
                   id%INFO(2) = 2
                ELSEIF ( size( id%ELTVAR ) < id%LELTVAR ) THEN 
                   id%INFO(1) = -22
                   id%INFO(2) = 2
                ELSEIF ( .not. associated( id%A_ELT ) ) THEN
                   id%INFO( 1 ) = -22
                   id%INFO( 2 ) = 4
                ELSE 
#if defined(MUMPS_F2003)
                   IF ( size( id%A_ELT, KIND=8 ) < id%KEEP8(30) ) THEN
#else
                   IF ( id%KEEP8(30) < int(huge(id%NZ),8) .AND.
     &                     size( id%A_ELT ) < int(id%KEEP8(30)) ) THEN
#endif
                      id%INFO( 1 ) = -22
                      id%INFO( 2 ) = 4
                   ENDIF
                END IF
             ENDIF
          ENDIF
C         ----------------------
C         Get the value of PERLU
C         ----------------------
          CALL MUMPS_GET_PERLU(id%KEEP(12),id%ICNTL(14),
     &         id%KEEP(50),id%KEEP(54),id%ICNTL(6),id%ICNTL(8))
C
C         ----------------------
C         Get null space options
C         Note that nullspace is forbidden in case of Schur complement
C         ----------------------
          CALL DMUMPS_GET_NS_OPTIONS_FACTO(id%N,id%KEEP(1),
     &                                     id%ICNTL(1),MPG)
C         ========================================
C         Decode and set scaling options for facto
C         ========================================
          IF (.NOT. ((id%KEEP(52).EQ.-2).AND.(id%ICNTL(8).EQ.77)) ) 
     &    THEN
C           if scaling was computed during analysis and automatic
C           choice of scaling then we do not recompute scaling
            id%KEEP(52)=id%ICNTL(8)
          ENDIF
          IF ( id%KEEP(52) .GT. 8 .OR. id%KEEP(52).LT.-2)
     &    id%KEEP(52) = 77
          IF ( id%KEEP(52) .EQ. 2 .OR. id%KEEP(52).EQ.5 
     &        .OR. id%KEEP(52) .EQ. 6 )
     &        id%KEEP(52) = 77
          IF (id%KEEP(52).EQ.77) THEN
            IF (id%KEEP(50).EQ.1) THEN
              ! for SPD matrices the default is "no scaling"
              id%KEEP(52) = 0
            ELSE
              ! SYM .ne. 1  the default is cheap SIMSCA
              id%KEEP(52) = 7 
            ENDIF
          ENDIF
          IF (id%KEEP(23) .NE. 0 .AND. id%ICNTL(8) .EQ. -1) THEN
             IF ( MPG .GT. 0 ) THEN
                WRITE(MPG,'(A)') ' ** WARNING : SCALING'
                WRITE(MPG,'(A)') 
     &               ' ** column permutation applied:'
                WRITE(MPG,'(A)') 
     &               ' ** column scaling has to be permuted'
             ENDIF 
          ENDIF
C         ------------------------
C         If Schur has been asked
C         for, scaling is disabled
C         ------------------------
          IF ( id%KEEP(60) .ne. 0 .and. id%KEEP(52) .ne. 0 ) THEN
            id%KEEP(52) = 0
            IF ( MPG .GT. 0 .AND. id%ICNTL(8) .NE. 0 ) THEN
              WRITE(MPG,'(A)') ' ** Warning: Scaling not applied.'
              WRITE(MPG,'(A)') ' ** (incompatibility with Schur)'
            END IF
          END IF
C         -------------------------------
C         If matrix is distributed on
C         entry, only options 7 and 8
C         of scaling are allowed.
C         -------------------------------
          IF (id%KEEP(54) .NE. 0 .AND. 
     &        id%KEEP(52).NE.7 .AND. id%KEEP(52).NE.8 .AND.
     &        id%KEEP(52) .NE. 0 ) THEN
             id%KEEP(52) = 0
             IF ( MPG .GT. 0 .and. id%ICNTL(8) .ne. 0 ) THEN
               WRITE(MPG,'(A)')
     &         ' ** Warning: This scaling option not available'
               WRITE(MPG,'(A)') ' ** for distributed matrix entry'
             END IF
          END IF
C         ------------------------------------
C         If matrix is symmetric, only scaling
C         options -1 (given scaling), 1
C         (diagonal scaling), 7 and 8 (SIMSCALING)
C         are allowed.
C         ------------------------------------
          IF ( id%KEEP(50) .NE. 0 ) THEN
             IF ( id%KEEP(52).ne.  1 .and.
     &            id%KEEP(52).ne. -1 .and.
     &            id%KEEP(52).ne.  0 .and.
     &            id%KEEP(52).ne.  7 .and.
     &            id%KEEP(52).ne.  8 .and.
     &            id%KEEP(52).ne. -2 .and.
     &            id%KEEP(52).ne. 77) THEN
              IF ( MPG .GT. 0 ) THEN
                WRITE(MPG,'(A)')
     &  ' ** Warning: Scaling option n.a. for symmetric matrix'
              END IF
              id%KEEP(52) = 0
            END IF
          END IF
C         ----------------------------------
C         If matrix is elemental on entry, 
C         automatic scaling is now forbidden
C         ----------------------------------
          IF (id%KEEP(55) .NE. 0 .AND. 
     &        ( id%KEEP(52) .gt. 0 ) ) THEN
            id%KEEP(52) = 0
            IF ( MPG .GT. 0 ) THEN
              WRITE(MPG,'(A)') ' ** Warning: Scaling not applied.'
              WRITE(MPG,'(A)')
     &        ' ** (only user scaling av. for elt. entry)'
            END IF
          END IF
C         --------------------------------------
C         Check input parameters ROWSCA / COLSCA
C         --------------------------------------
          IF ( id%KEEP(52) .eq. -1 ) THEN
            IF ( .not. associated( id%ROWSCA ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 5
            ELSE IF ( size( id%ROWSCA ) < id%N ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 5
            ELSE IF ( .not. associated( id%COLSCA ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 6
            ELSE IF ( size( id%COLSCA ) < id%N ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 6
            END IF
          END IF
C
C  Allocate -- if required,
C  ROWSCA and COLSCA on the master
C
C  Allocation of scaling arrays.
C  IF (KEEP(52)==-2 then scaling should have been allocated
C  and computed during analysis
C
C  If ICNTL(8) == -1, ROWSCA and COLSCA must have been associated and
C  filled by the user. If ICNTL(8) is >0 and <= 8, the scaling is
C  computed at the beginning of DMUMPS_FAC_DRIVER and is allocated now.
C
          IF (id%KEEP(52).GT.0 .AND.
     &        id%KEEP(52) .LE.8) THEN
            IF ( associated(id%COLSCA))
     &             DEALLOCATE( id%COLSCA )
            IF ( associated(id%ROWSCA))
     &             DEALLOCATE( id%ROWSCA )
            ALLOCATE( id%COLSCA(id%N), stat=IERR)
            IF (IERR .GT.0) THEN
               id%INFO(1)=-13
               id%INFO(2)=id%N
            ENDIF
            ALLOCATE( id%ROWSCA(id%N), stat=IERR)
            IF (IERR .GT.0) THEN
               id%INFO(1)=-13
               id%INFO(2)=id%N
            ENDIF
          END IF
C
C         Allocate scaling arrays of size 1 if
C         they are not used to avoid problems
C         when passing them in arguments
C
          IF (.NOT. associated(id%COLSCA)) THEN
            ALLOCATE( id%COLSCA(1), stat=IERR)
          END IF
          IF (IERR .GT.0) THEN
             id%INFO(1)=-13
             id%INFO(2)=1
          ENDIF
          IF (.NOT. associated(id%ROWSCA))
     &    ALLOCATE( id%ROWSCA(1), stat=IERR)
          IF (IERR .GT.0) THEN
            id%INFO(1)=-13
            id%INFO(2)=1
            IF ( LPOK ) WRITE(LP,'(A)')
     &         'Problems in allocations before facto'
            GOTO 200
          END IF
          IF (id%KEEP(252) .EQ. 1) THEN
             CALL DMUMPS_CHECK_DENSE_RHS
     &       (id%RHS,id%INFO,id%N,id%NRHS,id%LRHS)
C            Sets KEEP(221) and do some checks
             CALL DMUMPS_SET_K221(id)
             CALL DMUMPS_CHECK_REDRHS(id)
          ENDIF
 200      CONTINUE
        END IF        ! End of IF (MYID .eq. MASTER)
C       KEEP(221) was set in DMUMPS_SET_K221 but not broadcast
        CALL MPI_BCAST( id%KEEP(221), 1, MPI_INTEGER, MASTER, id%COMM,
     &       IERR )
C
C       Check distributed matrices on all processors.
        I_AM_SLAVE = ( id%MYID .ne. MASTER  .OR.
     &     ( id%MYID .eq. MASTER .AND.
     &     id%KEEP(46) .eq. 1 ) )
        IF (I_AM_SLAVE .AND.
     &      id%KEEP(54).NE.0 .AND. id%KEEP8(29).GT.0_8) THEN
           IF ( .not. associated( id%IRN_loc ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 16
#if defined(MUMPS_F2003)
           ELSE IF ( size( id%IRN_loc, KIND=8 ) < id%KEEP8(29) ) THEN
#else
C     size with kind=8 output not available. One can still
C     check that if NZ_loc can be stored in a 32-bit integer,
C     the 32-bit size(id%IRN_loc) (which we then assume not
C     to overflow...) is large enough
           ELSE IF ( id%KEEP8(29) .LE. int(huge(id%NZ_loc),8) .AND.
     &             size(id%IRN_loc) < int(id%KEEP8(29)) ) THEN
#endif
              id%INFO(1) = -22
              id%INFO(2) = 16
           ELSE IF ( .not. associated( id%JCN_loc ) ) THEN
              id%INFO(1) = -22
              id%INFO(2) = 16
#if defined(MUMPS_F2003)
           ELSE IF ( size( id%JCN_loc, KIND=8 ) < id%KEEP8(29) ) THEN
#else
C     Same as for IRN_loc above
           ELSE IF ( id%KEEP8(29) .LE. int(huge(id%NZ_loc),8) .AND.
     &             size(id%JCN_loc) < int(id%KEEP8(29)) ) THEN
#endif
              id%INFO(1) = -22
              id%INFO(2) = 16
           ELSEIF ( .not. associated( id%A_loc ) ) THEN
              id%INFO( 1 ) = -22
              id%INFO( 2 ) = 16
#if defined(MUMPS_F2003)
           ELSE IF ( size( id%A_loc, KIND=8 ) < id%KEEP8(29) ) THEN
#else
C     Same as for IRN_loc/JCN_loc above
           ELSE IF ( id%KEEP8(29) .LE. int(huge(id%NZ_loc),8) .AND.
     &             size( id%A_loc ) < int(id%KEEP8(29)) ) THEN
#endif
              id%INFO( 1 ) = -22
              id%INFO( 2 ) = 16
           END IF
        ENDIF
C
C  Check Schur complement on all processors.
C  DMUMPS_PROPINFO will be called right after those checks.
C
        IF (id%KEEP(60).EQ.2.OR.id%KEEP(60).EQ.3) THEN
          IF ( id%root%yes ) THEN
            IF ( associated( id%SCHUR_CINTERFACE )) THEN
C             Called from C interface...
C             The next instruction may cause
C             bound check errors at runtime
C             id%SCHUR=>id%SCHUR_CINTERFACE
C    &          (1:id%SCHUR_LLD*(id%root%SCHUR_NLOC-1)+
C    &          id%root%SCHUR_MLOC)
C             Instead, we set a temporary
C             pointer and then retrieve it
              CALL DMUMPS_SET_TMP_PTR(id%SCHUR_CINTERFACE(1),
     &         int(id%SCHUR_LLD,8)*int(id%root%SCHUR_NLOC-1,8)+
     &         int(id%root%SCHUR_MLOC,8))
              CALL DMUMPS_GET_TMP_PTR(id%SCHUR)
              NULLIFY(id%SCHUR_CINTERFACE)
            ENDIF
C           Check that SCHUR_LLD is large enough
            IF (id%SCHUR_LLD < id%root%SCHUR_MLOC) THEN
              IF (LP.GT.0) write(LP,*) 
     &          ' SCHUR leading dimension SCHUR_LLD ', 
     &          id%SCHUR_LLD, 'too small with respect to', 
     &          id%root%SCHUR_MLOC
              id%INFO(1)=-30
              id%INFO(2)=id%SCHUR_LLD
            ELSE IF ( .NOT. associated (id%SCHUR)) THEN
              IF (LP.GT.0) write(LP,'(A)') 
     &                      ' SCHUR not associated'
              id%INFO(1)=-22
              id%INFO(2)=9
            ELSE IF (size(id%SCHUR) <
     &          id%SCHUR_LLD*(id%root%SCHUR_NLOC-1)+
     &          id%root%SCHUR_MLOC) THEN
              IF (LP.GT.0) THEN 
                write(LP,'(A)') 
     &                      ' SCHUR allocated but too small'
                write(LP,*) id%MYID, ' : Size Schur=', 
     &          size(id%SCHUR), 
     &          ' SCHUR_LLD= ', id%SCHUR_LLD, 
     &          ' SCHUR_MLOC=', id%root%SCHUR_NLOC, 
     &          ' SCHUR_NLOC=', id%root%SCHUR_NLOC
              ENDIF
              id%INFO(1)=-22
              id%INFO(2)= 9
            ELSE
C              We initialize the pointer that
C              we will use within DMUMPS here.
               id%root%SCHUR_LLD=id%SCHUR_LLD
               IF (id%root%SCHUR_NLOC==0) THEN
                  ALLOCATE(id%root%SCHUR_POINTER(1), stat=IERR)
                  IF (IERR .GT.0) THEN
                     id%INFO(1)=-13
                     id%INFO(2)=1
                     IF ( LPOK ) THEN
                        WRITE(LP,'(A)')
     &                       'Problems in allocations before facto'
                     ENDIF
                  END IF
               ELSE
                id%root%SCHUR_POINTER=>id%SCHUR
               ENDIF
            ENDIF
          ENDIF
        ENDIF
C       -------------------------
C       Propagate possible errors
C       -------------------------
        CALL MUMPS_PROPINFO( id%ICNTL(1),
     &                      id%INFO(1),
     &                      id%COMM, id%MYID )
        IF ( id%INFO(1) .LT. 0 ) GO TO 499
C       -----------------------------------------------
C       Call factorization procedure DMUMPS_FAC_DRIVER
C       -----------------------------------------------
        CALL DMUMPS_FAC_DRIVER(id)
C       Save scaling in INFOG(33)
        IF (id%MYID .eq. MASTER) id%INFOG(33)=id%KEEP(52)
C
C       In the case of Schur, free or not associated
C       id%root%SCHUR_POINTER now rather than in end_driver.F
C       (Case of repeated factorizations).
        IF (id%KEEP(60).EQ.2.OR.id%KEEP(60).EQ.3) THEN
           IF (id%root%yes) THEN
              IF (id%root%SCHUR_NLOC==0) THEN
                 DEALLOCATE(id%root%SCHUR_POINTER)
                 NULLIFY(id%root%SCHUR_POINTER)
              ELSE
                 NULLIFY(id%root%SCHUR_POINTER)
              ENDIF
           ENDIF
        ENDIF
C     root%RG2L_ROW and root%RG2L_COL
C     are not used outside of the facto
        IF (associated(id%root%RG2L_ROW))THEN
           DEALLOCATE(id%root%RG2L_ROW)
           NULLIFY(id%root%RG2L_ROW)
        ENDIF
        IF (associated(id%root%RG2L_COL))THEN
           DEALLOCATE(id%root%RG2L_COL)
           NULLIFY(id%root%RG2L_COL)
        ENDIF
        IF (id%MYID .eq. MASTER) THEN
           CALL MUMPS_SECFIN(TIMEG)
           id%DKEEP(91) = TIMEG
        ENDIF
        IF (PROKG) THEN
            WRITE( MPG,'(/A,F12.4)')
     &         ' Elapsed time in factorization driver= ', TIMEG
        END IF 
C
C       Check for errors after FACTO
C       (it was propagated inside)
        IF(id%INFO(1).LT.0) THEN
C     Free id%S if facto failed     
           if (associated(id%S))  then 
              DEALLOCATE(id%S)
              NULLIFY(id%S)
           endif
           GO TO 499
        ENDIF
C
C       Update last successful step
C
        id%KEEP(40) = 2 - 456789
      END IF
C-------------------------------------------------------
C-
C      
C      BEGIN SOLVE PHASE
C
C-
C-------------------------------------------------------     
      IF (LSOLVE) THEN
        IF (id%MYID .eq. MASTER) THEN
           id%DKEEP(111)=0.0D0
           CALL MUMPS_SECDEB(TIMEG)
        END IF 
C       ---------------------
C       Reset KEEP(40) to 2.
C       (last successful step
C       was facto)
C       ---------------------
        id%KEEP(40) = 2 -456789
C       ------------------------------------------
C       Call solution procedure DMUMPS_SOLVE_DRIVER
C       ------------------------------------------
        IF (id%MYID .eq. MASTER) THEN
           KEEP235SAVE = id%KEEP(235)
           KEEP242SAVE = id%KEEP(242)
           KEEP243SAVE = id%KEEP(243)
           KEEP495SAVE = id%KEEP(495)
           KEEP497SAVE = id%KEEP(497)
           ! if no permutation of RHS asked then suppress request
           ! to interleave the RHS
           ! to interleave the RHS on ordering given then 
           ! using option to set permutation to identity should be 
           ! used (note though that 
           ! they # with A-1/sparseRHS and Null Space)
           IF (id%KEEP(242).EQ.0) id%KEEP(243)=0
C     --------------------------------------
C     Check input parameters ROWSCA / COLSCA
C     Only if KEEP(52).NE.0 because
C     only 0 means that no colsca/rowsca are needed
C     --------------------------------------
           IF ( id%KEEP(52) .ne. 0) THEN
              IF ( .not. associated( id%ROWSCA ) ) THEN
                 id%INFO(1) = -22
                 id%INFO(2) = 5
              ELSE IF ( size( id%ROWSCA ) < id%N ) THEN
                 id%INFO(1) = -22
                 id%INFO(2) = 5
              ELSE IF ( .not. associated( id%COLSCA ) ) THEN
                 id%INFO(1) = -22
                 id%INFO(2) = 6
              ELSE IF ( size( id%COLSCA ) < id%N ) THEN
                 id%INFO(1) = -22
                 id%INFO(2) = 6
              END IF
           ENDIF
        ENDIF
C     -------------------------
C     Propagate possible errors
C     -------------------------
        CALL MUMPS_PROPINFO( id%ICNTL(1),
     &       id%INFO(1),
     &       id%COMM, id%MYID )
        IF ( id%INFO(1) .LT. 0 ) GO TO 499
        CALL DMUMPS_SOLVE_DRIVER(id)
        IF (id%MYID .eq. MASTER) THEN
           CALL MUMPS_SECFIN(TIMEG)
            id%DKEEP(111) = TIMEG
        ENDIF
        IF (PROKG) THEN
            WRITE( MPG,'(/A,F12.4)')
     &         ' Elapsed time in solve driver= ', TIMEG
        END IF 
        IF (id%MYID .eq. MASTER) THEN
           id%KEEP(235) = KEEP235SAVE
           id%KEEP(242) = KEEP242SAVE
           id%KEEP(243) = KEEP243SAVE
           id%KEEP(495) = KEEP495SAVE
           id%KEEP(497) = KEEP497SAVE
        ENDIF
        IF (id%INFO(1).LT.0) GOTO 499
C       ---------------------------
C       Update last successful step
C       ---------------------------
        id%KEEP(40) = 3 -456789
      ENDIF
C
C  What was actually done is saved in KEEP(40)
C
      IF (PROK) CALL DMUMPS_PRINT_ICNTL(id, MP)
      GOTO 500
*
*=================
* ERROR section
*=================
  499 CONTINUE
*     Print error message if PROK
      IF (LPOK) WRITE (LP,99995) id%INFO(1)
      IF (LPOK) WRITE (LP,99994) id%INFO(2)
*
500   CONTINUE
#if ! defined(LARGEMATRICES)
C     ---------------------------------
C     Permute JCN on output to DMUMPS if
C     KEEP(23) is different from 0.
C     ---------------------------------
      IF (id%MYID .eq. MASTER .AND. id%KEEP(23) .NE. 0
     &    .AND. NOERRORBEFOREPERM) THEN
C       -------------------------------
C       IF JOB=3 and PERM was not
C       done (no iterative refinement/
C       error analysis), then we do not
C       permute JCN back.
C       -------------------------------
        IF (id%JOB .NE. 3 .OR. UNS_PERM_DONE) THEN
         IF (.not.associated(id%UNS_PERM)) THEN
C         I may happen 
C           (for ex in case of error -7 during analysis:
C           UNS_PERM can be not associated, 
C           KEEP(23) was set to to automatic choice(=7) and
C           an error of memory allocation occurs during analysis
C           before having decided value of KEEP(23))
C         UNS_PERM not associated and KEEP(23).NE.0
C         Permuting JCN back does not make sense and KEEP(23) 
C         should be reset to zero
          id%KEEP(23) = 0
         ELSE
           DO I8 = 1_8, id%KEEP8(28)
            J=id%JCN(I8)
C           -- skip out-of range (that are ignored in ANA_O)
            IF (J.LE.0.OR.J.GT.id%N) CYCLE
            id%JCN(I8)=id%UNS_PERM(J)
           END DO
          ENDIF
        END IF
      END IF
#endif
 510  CONTINUE
C     ------------------------------------
C     Set INFOG(1:2): same value on all
C     processors + broadcast other entries
C     ------------------------------------
      CALL DMUMPS_SET_INFOG(id%INFO(1), id%INFOG(1), id%COMM, id%MYID)
C
C     --------------------------------
C     Broadcast RINFOG entries to make
C     them available on all procs.
C     --------------------------------
      CALL MPI_BCAST( id%RINFOG(1), 40, MPI_DOUBLE_PRECISION, MASTER,
     &                    id%COMM, IERR )
      IF (id%INFOG(1).GE.0 .AND. JOB.NE.-1  
     &     .AND. JOB.NE.-2 ) THEN
         IF (id%MYID .eq. MASTER) THEN
            CALL MUMPS_SECFIN(TIMETOTAL)
            id%DKEEP(70) = TIMETOTAL
         ENDIF
      ENDIF
*=======================
* Compute space for save
*=======================
      IF (id%INFOG(1).GE.0) THEN
         CALL DMUMPS_COMPUTE_MEMORY_SAVE(id,FILE_SIZE,STRUC_SIZE)
         id%KEEP8(55)=FILE_SIZE
         call MPI_ALLREDUCE(id%KEEP8(55),id%KEEP8(57),1,
     &        MPI_INTEGER8, MPI_SUM,id%COMM,IERR)
         id%KEEP8(56)=STRUC_SIZE
         call MPI_ALLREDUCE(id%KEEP8(56),id%KEEP8(58),1,
     &        MPI_INTEGER8, MPI_SUM,id%COMM,IERR)
         id%RINFO(7)=dble(id%KEEP8(55))/1D6
         id%RINFO(8)=dble(id%KEEP8(56))/1D6
         id%RINFOG(17)=dble(id%KEEP8(57))/1D6
         id%RINFOG(18)=dble(id%KEEP8(58))/1D6
      ENDIF
!$    IF (ICNTL16_LOC .GT. 0) THEN
#if defined(WORKAROUNDINTELILP64OPENMPLIMITATION)
!$      CALL omp_set_num_threads(int(PREVIOUS_OMP_THREADS_NUM,4))
#else
!$      CALL omp_set_num_threads(PREVIOUS_OMP_THREADS_NUM)
#endif
!$      ICNTL16_LOC = 0
!$    ENDIF
*===============
* ERRORG section
*===============
      IF (id%MYID.EQ.MASTER.and.MPG.GT.0.and.
     & id%INFOG(1).lt.0) THEN
        WRITE(MPG,'(A,I16)') ' On return from DMUMPS, INFOG(1)=',
     &      id%INFOG(1)
        WRITE(MPG,'(A,I16)') ' On return from DMUMPS, INFOG(2)=',
     &      id%INFOG(2)
      END IF
C     -------------------------
C     Restore user communicator
C     -------------------------
       CALL MPI_COMM_FREE( id%COMM, IERR )
       id%COMM = COMM_SAVE
      RETURN
*
99995 FORMAT (' ** ERROR RETURN ** FROM DMUMPS INFO(1)=', I5)
99994 FORMAT (' ** INFO(2)=', I16)
99993 FORMAT (' ** Allocation error: could not permute JCN.')
      END SUBROUTINE DMUMPS
*
      SUBROUTINE DMUMPS_SET_INFOG( INFO, INFOG, COMM, MYID )
      IMPLICIT NONE
      INCLUDE 'mpif.h'
C
C  Purpose:
C  =======
C
C  If one proc has INFO(1).lt.0 and INFO(1) .ne. -1,
C  puts INFO(1:2) of this proc on all procs in INFOG
C
C  Arguments:
C  =========
C
      INTEGER, PARAMETER :: SIZE_INFOG = 80
      INTEGER :: INFO(80)
      INTEGER :: INFOG(SIZE_INFOG)  ! INFOG(80)
      INTEGER :: COMM, MYID
C
C  Local variables
C  ===============
C
#if defined(WORKAROUNDINTELILP64MPI2INTEGER)
      INTEGER(4) :: TMP1(2),TMP(2)
#else
      INTEGER :: TMP1(2),TMP(2)
#endif
      INTEGER ROOT, IERR
      INTEGER MASTER
      PARAMETER (MASTER=0)
C
C
      IF ( INFO(1) .ge. 0 ) THEN
C
C       This can only happen if the phase was successful
C       on all procs. If one proc failed, then all other
C       procs would have INFO(1)=-1.
C
        INFOG(1) = INFO(1)
        INFOG(2) = INFO(2)
      ELSE
C       ---------------------
C       Find who has smallest
C       error code INFO(1)
C       ---------------------
        INFOG(1) = INFO(1)
C        INFOG(2) = MYID
        TMP1(1) = INFO(1)
        TMP1(2) = MYID
        CALL MPI_ALLREDUCE(TMP1,TMP,1,MPI_2INTEGER,
     &                     MPI_MINLOC,COMM,IERR )
        INFOG(2) = INFO(2)
        ROOT = TMP(2)
        CALL MPI_BCAST( INFOG(1), 1, MPI_INTEGER, ROOT, COMM, IERR )
        CALL MPI_BCAST( INFOG(2), 1, MPI_INTEGER, ROOT, COMM, IERR )
      END IF
C
C    Make INFOG available on all procs:
C
      CALL MPI_BCAST(INFOG(3), SIZE_INFOG-2, MPI_INTEGER,
     &               MASTER, COMM, IERR )
      RETURN
      END SUBROUTINE DMUMPS_SET_INFOG
C--------------------------------------------------------------------
      SUBROUTINE DMUMPS_PRINT_ICNTL (id, LP)
      USE DMUMPS_STRUC_DEF
*
* Purpose: 
*   Print main control parameters CNTL and ICNTL 
*
*  ==========
*  Parameters
*  ==========
      TYPE (DMUMPS_STRUC), TARGET, INTENT(IN) :: id
      INTEGER  :: LP
** Local Variables
      INTEGER, POINTER :: JOB 
      INTEGER,DIMENSION(:),POINTER::ICNTL
      DOUBLE PRECISION,   DIMENSION(:),POINTER::CNTL
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      IF (LP.LE.0) RETURN
      JOB=>id%JOB
      ICNTL=>id%ICNTL
      CNTL=>id%CNTL
      IF (id%MYID.EQ.MASTER) THEN
         SELECT CASE (JOB)
         CASE(1);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           IF (id%SYM.EQ.2) THEN
            WRITE (LP,991) ICNTL(5),ICNTL(6),ICNTL(7),ICNTL(12),
     &          ICNTL(13),
     &          ICNTL(18),ICNTL(19),ICNTL(22)
           ELSE
            WRITE (LP,891) ICNTL(5),ICNTL(6),ICNTL(7),
     &          ICNTL(13),
     &          ICNTL(18),ICNTL(19),ICNTL(22)
           ENDIF
           IF ((ICNTL(6).EQ.5).OR.(ICNTL(6).EQ.6).OR.
     &          (ICNTL(12).NE.1) )  THEN
              WRITE (LP,992) ICNTL(8)
           ENDIF   
           IF (id%ICNTL(19).NE.0)
     &      WRITE(LP,998) id%SIZE_SCHUR
           WRITE (LP,993) ICNTL(14)
         CASE(2);
           WRITE (LP,980) 
           WRITE (LP,981) CNTL(1), CNTL(3), CNTL(4), CNTL(5), CNTL(7)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,992) ICNTL(8)
           WRITE (LP,993) ICNTL(14) 
           WRITE (LP,923) ICNTL(24), ICNTL(31), ICNTL(32), ICNTL(33),
     &                    ICNTL(35), ICNTL(36)
         CASE(3);
           WRITE (LP,980)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,995)
     &     ICNTL(9),ICNTL(10),ICNTL(11),ICNTL(20),ICNTL(21)
         CASE(4);
           WRITE (LP,980) 
           WRITE (LP,981) CNTL(1), CNTL(3), CNTL(4), CNTL(5), CNTL(7)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,992) ICNTL(8)
           IF (id%ICNTL(19).NE.0)
     &      WRITE(LP,998) id%SIZE_SCHUR
           WRITE (LP,993) ICNTL(14)
           WRITE (LP,923) ICNTL(24), ICNTL(31), ICNTL(32), ICNTL(33),
     &                    ICNTL(35), ICNTL(36)
         CASE(5);
           WRITE (LP,980) 
           WRITE (LP,981) CNTL(1), CNTL(3), CNTL(4), CNTL(5), CNTL(7)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           IF (id%SYM.EQ.2) THEN
            WRITE (LP,991) ICNTL(5),ICNTL(6),ICNTL(7),ICNTL(12),
     &          ICNTL(13),
     &          ICNTL(18),ICNTL(19),ICNTL(22)
           ELSE
            WRITE (LP,891) ICNTL(5),ICNTL(6),ICNTL(7),
     &          ICNTL(13),
     &          ICNTL(18),ICNTL(19),ICNTL(22)
           ENDIF
           WRITE (LP,992) ICNTL(8)
           WRITE (LP,993) ICNTL(14)
           WRITE (LP,995)
     &     ICNTL(9),ICNTL(10),ICNTL(11),ICNTL(20),ICNTL(21)
           WRITE (LP,923) ICNTL(24), ICNTL(31), ICNTL(32), ICNTL(33),
     &                    ICNTL(35), ICNTL(36)
         CASE(6);
           WRITE (LP,980)
           WRITE (LP,981) CNTL(1), CNTL(3), CNTL(4), CNTL(5), CNTL(7)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           IF (id%SYM.EQ.2) THEN
            WRITE (LP,991) ICNTL(5),ICNTL(6),ICNTL(7),ICNTL(12),
     &          ICNTL(13),
     &          ICNTL(18),ICNTL(19),ICNTL(22)
           ELSE
            WRITE (LP,891) ICNTL(5),ICNTL(6),ICNTL(7),
     &          ICNTL(13),
     &          ICNTL(18),ICNTL(19),ICNTL(22)
           ENDIF
           IF (id%ICNTL(19).NE.0)
     &      WRITE(LP,998) id%SIZE_SCHUR
           WRITE (LP,992) ICNTL(8)
           WRITE (LP,995)
     &     ICNTL(9),ICNTL(10),ICNTL(11),ICNTL(20),ICNTL(21)
           WRITE (LP,993) ICNTL(14)
           WRITE (LP,923) ICNTL(24), ICNTL(31), ICNTL(32), ICNTL(33),
     &                    ICNTL(35), ICNTL(36)
        END SELECT
      ENDIF
 980  FORMAT (/'***********CONTROL PARAMETERS (ICNTL)**************'/)
 981  FORMAT (
     &     ' CNTL(1)   Threshold for numerical pivoting        =',D16.4/
     &     ' CNTL(3)   Null pivot detection threshold          =',D16.4/
     &     ' CNTL(4)   Threshold for static pivoting           =',D16.4/
     &     ' CNTL(5)   Fixation for null pivots                =',D16.4/
     &     ' CNTL(7)   Dropping threshold for BLR compression  =',D16.4)
 990  FORMAT (
     &     'ICNTL(1)   Output stream for error messages        =',I10/
     &     'ICNTL(2)   Output stream for diagnostic messages   =',I10/
     &     'ICNTL(3)   Output stream for global information    =',I10/
     &     'ICNTL(4)   Level of printing                       =',I10)
 991  FORMAT (
     &     'ICNTL(5)   Matrix format  ( keep(55) )             =',I10/
     &     'ICNTL(6)   Maximum transversal  ( keep(23) )       =',I10/
     &     'ICNTL(7)   Ordering                                =',I10/
     &     'ICNTL(12)  LDLT ordering strat ( keep(95) )        =',I10/
     &     'ICNTL(13)  Parallel root (0=on, 1=off)             =',I10/
     &     'ICNTL(18)  Distributed matrix  ( keep(54) )        =',I10/
     &     'ICNTL(19)  Schur option ( keep(60) 0=off,else=on ) =',I10/
     &     'ICNTL(22)  Out-off-core option (0=off, >0=on)      =',I10)
 891  FORMAT (
     &     'ICNTL(5)   Matrix format  ( keep(55) )             =',I10/
     &     'ICNTL(6)   Maximum transversal  ( keep(23) )       =',I10/
     &     'ICNTL(7)   Ordering                                =',I10/
     &     'ICNTL(13)  Parallel root (0=on, 1=off)             =',I10/
     &     'ICNTL(18)  Distributed matrix  ( keep(54) )        =',I10/
     &     'ICNTL(19)  Schur option ( keep(60) 0=off,else=on ) =',I10/
     &     'ICNTL(22)  Out-off-core option (0=off, >0=on)      =',I10)
 992  FORMAT (
     &     'ICNTL(8)   Scaling strategy                        =',I10)
 923  FORMAT (
     &     'ICNTL(24)  Null pivot detection (0=off)            =',I10/
     &     'ICNTL(31)  Discard factors (0=off, else=on)        =',I10/
     &     'ICNTL(32)  Forward elimination during facto (0=off)=',I10/
     &     'ICNTL(33)  Compute determinant (0=off)             =',I10/
     &     'ICNTL(35)  Block Low Rank (BLR, 0=off >0=on)       =',I10/
     &     'ICNTL(36)  BLR variant                             =',I10)
 993  FORMAT (
     &     'ICNTL(14)  Percent of memory increase              =',I10)
 995  FORMAT (
     &     'ICNTL(9)   Solve A x=b (1) or A''x = b (else)       =',I10/
     &     'ICNTL(10)  Max steps iterative refinement          =',I10/
     &     'ICNTL(11)  Error analysis (1=all,2=some,else=off)  =',I10/
     &     'ICNTL(20)  Den.(0)/sparse(1,2,3)/dist.(10,11) RHS  =',I10/
     &     'ICNTL(21)  Gathered (0) or distributed(1) solution =',I10)
 998  FORMAT (
     &     '      Size of SCHUR matrix (SIZE_SHUR)             =',I10)
      END SUBROUTINE DMUMPS_PRINT_ICNTL
C--------------------------------------------------------------------
      SUBROUTINE DMUMPS_PRINT_KEEP(id, LP)
      USE DMUMPS_STRUC_DEF
*
*  ==========
*  Parameters
*  ==========
      TYPE (DMUMPS_STRUC), TARGET, INTENT(IN) :: id
      INTEGER ::LP
** Local Variables
      INTEGER, POINTER :: JOB 
      INTEGER,DIMENSION(:),POINTER::ICNTL, KEEP
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      IF (LP.LE.0) RETURN
      JOB=>id%JOB
      ICNTL=>id%ICNTL
      KEEP=>id%KEEP
      IF (id%MYID.EQ.MASTER) THEN
         SELECT CASE (JOB)
         CASE(1);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) KEEP(55),KEEP(23),ICNTL(7),KEEP(95),
     &          ICNTL(13),KEEP(54),KEEP(60),ICNTL(22)
           IF ((KEEP(23).EQ.5).OR.(KEEP(23).EQ.6))THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           WRITE (LP,993) KEEP(12)
         CASE(2);
           WRITE (LP,980)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           IF (KEEP(23).EQ.0)THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           WRITE (LP,993) KEEP(12)
         CASE(3);
           WRITE (LP,980)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4) 
           WRITE (LP,995)
     &     ICNTL(9),ICNTL(10),ICNTL(11),ICNTL(20),ICNTL(21)
         CASE(4);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           IF (KEEP(23).NE.0)THEN
              WRITE (LP,992) KEEP(52)
           ENDIF  
           WRITE (LP,991) KEEP(55),KEEP(23),ICNTL(7),KEEP(95),
     &          ICNTL(13),KEEP(54),KEEP(60),ICNTL(22)
           WRITE (LP,995)
     &     ICNTL(9),ICNTL(10),ICNTL(11),ICNTL(20),ICNTL(21)
           WRITE (LP,993) KEEP(12)
         CASE(5);
           WRITE (LP,980) 
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) KEEP(55),KEEP(23),ICNTL(7),KEEP(95),
     &          ICNTL(13),KEEP(54),KEEP(60),ICNTL(22)
           IF ((KEEP(23).EQ.5).OR.(KEEP(23).EQ.6)
     &       .OR. (KEEP(23).EQ.7)) THEN
              WRITE (LP,992) KEEP(52)
           ENDIF              
           IF (KEEP(23).EQ.0)THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           WRITE (LP,993) KEEP(12)
         CASE(6);
           WRITE (LP,980)
           WRITE (LP,990) ICNTL(1),ICNTL(2),ICNTL(3),ICNTL(4)
           WRITE (LP,991) KEEP(55),KEEP(23),ICNTL(7),KEEP(95),
     &          ICNTL(13),KEEP(54),KEEP(60),ICNTL(22)
           IF ((KEEP(23).EQ.5).OR.(KEEP(23).EQ.6)
     &       .OR. (KEEP(23).EQ.7)) THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           IF (KEEP(23).EQ.0)THEN
              WRITE (LP,992) KEEP(52)
           ENDIF   
           WRITE (LP,995)
     &     ICNTL(9),ICNTL(10),ICNTL(11),KEEP(248),ICNTL(21)
           WRITE (LP,993) KEEP(12)
        END SELECT
      ENDIF
 980  FORMAT (/'******INTERNAL VALUE OF PARAMETERS (ICNTL/KEEP)****'/)
 990  FORMAT (
     &     'ICNTL(1)   Output stream for error messages        =',I10/
     &     'ICNTL(2)   Output stream for diagnostic messages   =',I10/
     &     'ICNTL(3)   Output stream for global information    =',I10/
     &     'ICNTL(4)   Level of printing                       =',I10)
 991  FORMAT (
     &     'ICNTL(5)   Matrix format  ( keep(55) )             =',I10/
     &     'ICNTL(6)   Maximum transversal  ( keep(23) )       =',I10/
     &     'ICNTL(7)   Ordering                                =',I10/
     &     'ICNTL(12)  LDLT ordering strat ( keep(95) )        =',I10/
     &     'ICNTL(13)  Parallel root (0=on, 1=off)             =',I10/
     &     'ICNTL(18)  Distributed matrix  ( keep(54) )        =',I10/
     &     'ICNTL(19)  Schur option ( keep(60) 0=off,else=on ) =',I10/
     &     'ICNTL(22)  Out-off-core option (0=Off, >0=ON)      =',I10)
 992  FORMAT (
     &     'ICNTL(8)   Scaling strategy ( keep(52) )           =',I10)
 993  FORMAT (
     &     'ICNTL(14)  Percent of memory increase ( keep(12) ) =',I10)
 995  FORMAT (
     &     'ICNTL(9)   Solve A x=b (1) or A''x = b (else)       =',I10/
     &     'ICNTL(10)  Max steps iterative refinement          =',I10/
     &     'ICNTL(11)  Error analysis ( 0= off, else=on)       =',I10/
     &     'ICNTL(20)  Den.(0)/sparse(1,2,3)/dist.(10,11) RHS  =',I10/
     &     'ICNTL(21)  Gathered (0) or distributed(1) solution =',I10)
      END SUBROUTINE DMUMPS_PRINT_KEEP
      SUBROUTINE DMUMPS_CHECK_DENSE_RHS
     &       (idRHS, idINFO, idN, idNRHS, idLRHS)
      IMPLICIT NONE
C
C  Purpose:
C  =======
C
C     Check that the dense RHS is associated and of
C     correct size. Called on master only, when dense
C     RHS is supposed to be allocated. This can be used
C     either at the beginning of the solve phase or
C     at the beginning of the factorization phase
C     if forward solve is done during factorization
C     (see ICNTL(32)) ; idINFO(1), idINFO(2) may be
C     modified.
C
C
C  Arguments:
C  =========
C
C     id* : see corresponding components of the main
C     MUMPS structure.
C
      DOUBLE PRECISION, DIMENSION(:), POINTER :: idRHS
      INTEGER, intent(in)    :: idN, idNRHS, idLRHS
      INTEGER, intent(inout) :: idINFO(:)
      IF ( .not. associated( idRHS ) ) THEN
              idINFO( 1 ) = -22
              idINFO( 2 ) = 7
      ELSE IF (idNRHS.EQ.1) THEN
               IF ( size( idRHS ) < idN ) THEN
                  idINFO( 1 ) = -22
                  idINFO( 2 ) = 7
               ENDIF
      ELSE IF (idLRHS < idN) 
     &            THEN
                  idINFO( 1 ) = -26
                  idINFO( 2 ) = idLRHS
      ELSE IF 
#if defined(MUMPS_F2003)
     &      (size(idRHS,kind=8) <
     &       int(idNRHS,8)*int(idLRHS,8)-idLRHS+idN)
#else
C           size with kind=8 not available. One can still
C           perform the check if minimal size small enough.
     &      (int(idNRHS,8)*int(idLRHS,8)-idLRHS+idN
     &                                           .LE. int(huge(idN),8)
     &      .and.
     &      size(idRHS) < int(int(idNRHS,8)*int(idLRHS,8)-idLRHS+idN))
#endif
     &            THEN
                  idINFO( 1 ) = -22
                  idINFO( 2 ) = 7
      END IF
      RETURN
      END SUBROUTINE DMUMPS_CHECK_DENSE_RHS
C
      SUBROUTINE DMUMPS_SET_K221(id)
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
C
C     Purpose:
C     =======
C
C     Sets KEEP(221) on master.
C     Constraint: must be called before DMUMPS_CHECK_REDRHS.
C     Can be called at factorization or solve phase
C
      TYPE (DMUMPS_STRUC) :: id
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      IF (id%MYID.EQ.MASTER) THEN
        id%KEEP(221)=id%ICNTL(26)
        IF (id%KEEP(221).ne.0 .and. id%KEEP(221) .NE.1
     &      .AND.id%KEEP(221).ne.2) id%KEEP(221)=0
      ENDIF
      RETURN
      END SUBROUTINE DMUMPS_SET_K221
C
      SUBROUTINE DMUMPS_CHECK_REDRHS(id)
      USE DMUMPS_STRUC_DEF
      IMPLICIT NONE
C
C  Purpose:
C  =======
C
C  * Decode API related to REDRHS and check REDRHS
C  * Can be called at factorization or solve phase
C  * Constraints:
C    - Must be called after solve phase.
C    - KEEP(60) must have been set (ok to check
C    since KEEP(60) was set during analysis phase)
C  * Remark that during solve phase, ICNTL(26)=1 is
C    forbidden in case of fwd in facto.
C
      TYPE (DMUMPS_STRUC) :: id
      INTEGER MASTER
      PARAMETER( MASTER = 0 )
      IF (id%MYID .EQ. MASTER) THEN
          IF ( id%KEEP(221) == 1 .or. id%KEEP(221) == 2 ) THEN
            IF (id%KEEP(221) == 2 .and. id%JOB == 2) THEN
              id%INFO(1)=-35
              id%INFO(2)=id%KEEP(221)
              GOTO 333
            ENDIF
            IF (id%KEEP(221) == 1 .and. id%KEEP(252) == 1
     &          .and. id%JOB == 3) THEN
              id%INFO(1)=-35
              id%INFO(2)=id%KEEP(221)
            ENDIF
            IF ( id%KEEP(60).eq. 0 .or. id%SIZE_SCHUR.EQ.0 ) THEN
              id%INFO(1)=-33
              id%INFO(2)=id%KEEP(221)
              GOTO 333
            ENDIF
            IF ( .NOT. associated( id%REDRHS)) THEN
              id%INFO(1)=-22
              id%INFO(2)=15
              GOTO 333
            ELSE IF (id%NRHS.EQ.1) THEN
              IF (size(id%REDRHS) < id%SIZE_SCHUR ) THEN
                id%INFO(1)=-22
                id%INFO(2)=15
                GOTO 333
              ENDIF
            ELSE IF (id%LREDRHS < id%SIZE_SCHUR) THEN
              id%INFO(1)=-34
              id%INFO(2)=id%LREDRHS
              GOTO 333
            ELSE IF
     &      (size(id%REDRHS)<
     &         id%NRHS*id%LREDRHS-id%LREDRHS+id%SIZE_SCHUR)
     &      THEN
              id%INFO(1)=-22
              id%INFO(2)=15
              GOTO 333
            ENDIF
          ENDIF
      ENDIF
 333  CONTINUE
C     Error is not propagated. It should be propagated outside.
C     The reason to propagate it outside is that there can be
C     one call to PROPINFO instead of several ones.
      RETURN
      END SUBROUTINE DMUMPS_CHECK_REDRHS
