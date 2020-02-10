/*
     Include file for the matrix component of PETSc
*/
#ifndef __PETSCMAT_H
#define __PETSCMAT_H
#include <petscvec.h>

/*S
     Mat - Abstract PETSc matrix object used to manage all linear operators in PETSc, even those without
           an explicit sparse representation (such as matrix-free operators)

   Level: beginner

  Concepts: matrix; linear operator

.seealso:  MatCreate(), MatType, MatSetType(), MatDestroy()
S*/
typedef struct _p_Mat*           Mat;

/*J
    MatType - String with the name of a PETSc matrix type

   Level: beginner

.seealso: MatSetType(), Mat, MatSolverType, MatRegister()
J*/
typedef const char* MatType;
#define MATSAME            "same"
#define MATMAIJ            "maij"
#define MATSEQMAIJ         "seqmaij"
#define MATMPIMAIJ         "mpimaij"
#define MATIS              "is"
#define MATAIJ             "aij"
#define MATSEQAIJ          "seqaij"
#define MATMPIAIJ          "mpiaij"
#define MATAIJCRL          "aijcrl"
#define MATSEQAIJCRL       "seqaijcrl"
#define MATMPIAIJCRL       "mpiaijcrl"
#define MATAIJCUSPARSE     "aijcusparse"
#define MATSEQAIJCUSPARSE  "seqaijcusparse"
#define MATMPIAIJCUSPARSE  "mpiaijcusparse"
#define MATAIJVIENNACL     "aijviennacl"
#define MATSEQAIJVIENNACL  "seqaijviennacl"
#define MATMPIAIJVIENNACL  "mpiaijviennacl"
#define MATAIJPERM         "aijperm"
#define MATSEQAIJPERM      "seqaijperm"
#define MATMPIAIJPERM      "mpiaijperm"
#define MATAIJSELL         "aijsell"
#define MATSEQAIJSELL      "seqaijsell"
#define MATMPIAIJSELL      "mpiaijsell"
#define MATAIJMKL          "aijmkl"
#define MATSEQAIJMKL       "seqaijmkl"
#define MATMPIAIJMKL       "mpiaijmkl"
#define MATBAIJMKL         "baijmkl"
#define MATSEQBAIJMKL      "seqbaijmkl"
#define MATMPIBAIJMKL      "mpibaijmkl"
#define MATSHELL           "shell"
#define MATDENSE           "dense"
#define MATSEQDENSE        "seqdense"
#define MATMPIDENSE        "mpidense"
#define MATELEMENTAL       "elemental"
#define MATBAIJ            "baij"
#define MATSEQBAIJ         "seqbaij"
#define MATMPIBAIJ         "mpibaij"
#define MATMPIADJ          "mpiadj"
#define MATSBAIJ           "sbaij"
#define MATSEQSBAIJ        "seqsbaij"
#define MATMPISBAIJ        "mpisbaij"
#define MATDAAD            "daad"
#define MATMFFD            "mffd"
#define MATNORMAL          "normal"
#define MATNORMALHERMITIAN "normalh"
#define MATLRC             "lrc"
#define MATSCATTER         "scatter"
#define MATBLOCKMAT        "blockmat"
#define MATCOMPOSITE       "composite"
#define MATFFT             "fft"
#define MATFFTW            "fftw"
#define MATSEQCUFFT        "seqcufft"
#define MATTRANSPOSEMAT    "transpose"
#define MATSCHURCOMPLEMENT "schurcomplement"
#define MATPYTHON          "python"
#define MATHYPRE           "hypre"
#define MATHYPRESTRUCT     "hyprestruct"
#define MATHYPRESSTRUCT    "hypresstruct"
#define MATSUBMATRIX       "submatrix"
#define MATLOCALREF        "localref"
#define MATNEST            "nest"
#define MATPREALLOCATOR    "preallocator"
#define MATSELL            "sell"
#define MATSEQSELL         "seqsell"
#define MATMPISELL         "mpisell"
#define MATDUMMY           "dummy"
#define MATLMVM            "lmvm"
#define MATLMVMDFP         "lmvmdfp"
#define MATLMVMBFGS        "lmvmbfgs"
#define MATLMVMSR1         "lmvmsr1"
#define MATLMVMBRDN        "lmvmbrdn"
#define MATLMVMBADBRDN     "lmvmbadbrdn"
#define MATLMVMSYMBRDN     "lmvmsymbrdn"
#define MATLMVMSYMBADBRDN  "lmvmsymbadbrdn"
#define MATLMVMDIAGBRDN    "lmvmdiagbrdn"

/*J
    MatSolverType - String with the name of a PETSc matrix solver type.

    For example: "petsc" indicates what PETSc provides, "superlu_dist" the parallel SuperLU_DIST package etc

   Level: beginner

.seealso: MatGetFactor(), Mat, MatSetType(), MatType
J*/
typedef const char* MatSolverType;
#define MATSOLVERSUPERLU          "superlu"
#define MATSOLVERSUPERLU_DIST     "superlu_dist"
#define MATSOLVERSTRUMPACK        "strumpack"
#define MATSOLVERUMFPACK          "umfpack"
#define MATSOLVERCHOLMOD          "cholmod"
#define MATSOLVERKLU              "klu"
#define MATSOLVERSPARSEELEMENTAL  "sparseelemental"
#define MATSOLVERELEMENTAL        "elemental"
#define MATSOLVERESSL             "essl"
#define MATSOLVERLUSOL            "lusol"
#define MATSOLVERMUMPS            "mumps"
#define MATSOLVERMKL_PARDISO      "mkl_pardiso"
#define MATSOLVERMKL_CPARDISO     "mkl_cpardiso"
#define MATSOLVERPASTIX           "pastix"
#define MATSOLVERMATLAB           "matlab"
#define MATSOLVERPETSC            "petsc"
#define MATSOLVERBAS              "bas"
#define MATSOLVERCUSPARSE         "cusparse"

/*E
    MatFactorType - indicates what type of factorization is requested

    Level: beginner

   Any additions/changes here MUST also be made in include/petsc/finclude/petscmat.h

.seealso: MatSolverType, MatGetFactor()
E*/
typedef enum {MAT_FACTOR_NONE, MAT_FACTOR_LU, MAT_FACTOR_CHOLESKY, MAT_FACTOR_ILU, MAT_FACTOR_ICC,MAT_FACTOR_ILUDT} MatFactorType;
PETSC_EXTERN const char *const MatFactorTypes[];

PETSC_EXTERN PetscErrorCode MatGetFactor(Mat,MatSolverType,MatFactorType,Mat*);
PETSC_EXTERN PetscErrorCode MatGetFactorAvailable(Mat,MatSolverType,MatFactorType,PetscBool *);
PETSC_EXTERN PetscErrorCode MatFactorGetSolverType(Mat,MatSolverType*);
PETSC_EXTERN PetscErrorCode MatGetFactorType(Mat,MatFactorType*);
PETSC_EXTERN PetscErrorCode MatSolverTypeRegister(MatSolverType,MatType,MatFactorType,PetscErrorCode(*)(Mat,MatFactorType,Mat*));
PETSC_EXTERN PetscErrorCode MatSolverTypeGet(MatSolverType,MatType,MatFactorType,PetscBool*,PetscBool*,PetscErrorCode (**)(Mat,MatFactorType,Mat*));
typedef MatSolverType MatSolverPackage PETSC_DEPRECATED("Use MatSolverType");
PETSC_DEPRECATED("Use MatSolverTypeRegister()") PETSC_STATIC_INLINE PetscErrorCode MatSolverPackageRegister(MatSolverType stype,MatType mtype,MatFactorType ftype,PetscErrorCode(*f)(Mat,MatFactorType,Mat*))
{ return MatSolverTypeRegister(stype,mtype,ftype,f); }
PETSC_DEPRECATED("Use MatSolverTypeGet()") PETSC_STATIC_INLINE PetscErrorCode MatSolverPackageGet(MatSolverType stype,MatType mtype,MatFactorType ftype,PetscBool *foundmtype,PetscBool *foundstype,PetscErrorCode(**f)(Mat,MatFactorType,Mat*))
{ return MatSolverTypeGet(stype,mtype,ftype,foundmtype,foundstype,f); }

/* Logging support */
#define    MAT_FILE_CLASSID 1211216    /* used to indicate matrices in binary files */
PETSC_EXTERN PetscClassId MAT_CLASSID;
PETSC_EXTERN PetscClassId MAT_COLORING_CLASSID;
PETSC_EXTERN PetscClassId MAT_FDCOLORING_CLASSID;
PETSC_EXTERN PetscClassId MAT_TRANSPOSECOLORING_CLASSID;
PETSC_EXTERN PetscClassId MAT_PARTITIONING_CLASSID;
PETSC_EXTERN PetscClassId MAT_COARSEN_CLASSID;
PETSC_EXTERN PetscClassId MAT_NULLSPACE_CLASSID;
PETSC_EXTERN PetscClassId MATMFFD_CLASSID;

/*E
    MatReuse - Indicates if matrices obtained from a previous call to MatCreateSubMatrices(), MatCreateSubMatrix(), MatConvert() or several other functions
     are to be reused to store the new matrix values.

$  MAT_INITIAL_MATRIX - create a new matrix
$  MAT_REUSE_MATRIX - reuse the matrix created with a previous call that used MAT_INITIAL_MATRIX
$  MAT_INPLACE_MATRIX - replace the first input matrix with the new matrix (not applicable to all functions)
$  MAT_IGNORE_MATRIX - do not create a new matrix or reuse a give matrix, just ignore that matrix argument (not applicable to all functions)

    Level: beginner

   Any additions/changes here MUST also be made in include/petsc/finclude/petscmat.h

.seealso: MatCreateSubMatrices(), MatCreateSubMatrix(), MatDestroyMatrices(), MatConvert()
E*/
typedef enum {MAT_INITIAL_MATRIX,MAT_REUSE_MATRIX,MAT_IGNORE_MATRIX,MAT_INPLACE_MATRIX} MatReuse;

/*E
    MatCreateSubMatrixOption - Indicates if matrices obtained from a call to MatCreateSubMatrices()
     include the matrix values. Currently it is only used by MatGetSeqNonzerostructure().

    Level: beginner

.seealso: MatGetSeqNonzerostructure()
E*/
typedef enum {MAT_DO_NOT_GET_VALUES,MAT_GET_VALUES} MatCreateSubMatrixOption;

PETSC_EXTERN PetscErrorCode MatInitializePackage(void);

PETSC_EXTERN PetscErrorCode MatCreate(MPI_Comm,Mat*);
PETSC_EXTERN PetscErrorCode MatSetSizes(Mat,PetscInt,PetscInt,PetscInt,PetscInt);
PETSC_EXTERN PetscErrorCode MatSetType(Mat,MatType);
PETSC_EXTERN PetscErrorCode MatSetFromOptions(Mat);
PETSC_STATIC_INLINE PetscErrorCode MatViewFromOptions(Mat A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}
PETSC_EXTERN PetscErrorCode MatRegister(const char[],PetscErrorCode(*)(Mat));
PETSC_EXTERN PetscErrorCode MatRegisterBaseName(const char[],const char[],const char[]);
PETSC_EXTERN PetscErrorCode MatSetOptionsPrefix(Mat,const char[]);
PETSC_EXTERN PetscErrorCode MatAppendOptionsPrefix(Mat,const char[]);
PETSC_EXTERN PetscErrorCode MatGetOptionsPrefix(Mat,const char*[]);
PETSC_EXTERN PetscErrorCode MatSetErrorIfFailure(Mat,PetscBool);

PETSC_EXTERN PetscFunctionList MatList;
PETSC_EXTERN PetscFunctionList MatColoringList;
PETSC_EXTERN PetscFunctionList MatPartitioningList;

/*E
    MatStructure - Indicates if two matrices have the same nonzero structure

    Level: beginner

   Any additions/changes here MUST also be made in include/petsc/finclude/petscmat.h

.seealso: MatCopy(), MatAXPY()
E*/
typedef enum {DIFFERENT_NONZERO_PATTERN,SUBSET_NONZERO_PATTERN,SAME_NONZERO_PATTERN} MatStructure;

#if defined PETSC_HAVE_MKL_SPARSE
PETSC_EXTERN PetscErrorCode MatCreateBAIJMKL(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);
PETSC_EXTERN PetscErrorCode MatCreateSeqBAIJMKL(MPI_Comm comm,PetscInt bs,PetscInt m,PetscInt n,PetscInt nz,const PetscInt nnz[],Mat *A);
#endif

PETSC_EXTERN PetscErrorCode MatCreateSeqSELL(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);
PETSC_EXTERN PetscErrorCode MatCreateSELL(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);
PETSC_EXTERN PetscErrorCode MatSeqSELLSetPreallocation(Mat,PetscInt,const PetscInt[]);
PETSC_EXTERN PetscErrorCode MatMPISELLSetPreallocation(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[]);

PETSC_EXTERN PetscErrorCode MatCreateSeqDense(MPI_Comm,PetscInt,PetscInt,PetscScalar[],Mat*);
PETSC_EXTERN PetscErrorCode MatCreateDense(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar[],Mat*);
PETSC_EXTERN PetscErrorCode MatCreateSeqAIJ(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);
PETSC_EXTERN PetscErrorCode MatCreateAIJ(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);
PETSC_EXTERN PetscErrorCode MatCreateMPIAIJWithArrays(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[],Mat *);
PETSC_EXTERN PetscErrorCode MatCreateMPIAIJWithSplitArrays(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt[],PetscInt[],PetscScalar[],PetscInt[],PetscInt[],PetscScalar[],Mat*);
PETSC_EXTERN PetscErrorCode MatCreateMPIAIJWithSeqAIJ(MPI_Comm,Mat,Mat,const PetscInt[],Mat*);

PETSC_EXTERN PetscErrorCode MatCreateSeqBAIJ(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);
PETSC_EXTERN PetscErrorCode MatCreateBAIJ(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);
PETSC_EXTERN PetscErrorCode MatCreateMPIBAIJWithArrays(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[],Mat*);

PETSC_EXTERN PetscErrorCode MatCreateMPIAdj(MPI_Comm,PetscInt,PetscInt,PetscInt[],PetscInt[],PetscInt[],Mat*);
PETSC_EXTERN PetscErrorCode MatCreateSeqSBAIJ(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);

PETSC_EXTERN PetscErrorCode MatCreateSBAIJ(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);
PETSC_EXTERN PetscErrorCode MatCreateMPISBAIJWithArrays(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[],Mat *);
PETSC_EXTERN PetscErrorCode MatSeqSBAIJSetPreallocationCSR(Mat,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[]);
PETSC_EXTERN PetscErrorCode MatMPISBAIJSetPreallocationCSR(Mat,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[]);
PETSC_EXTERN PetscErrorCode MatXAIJSetPreallocation(Mat,PetscInt,const PetscInt[],const PetscInt[],const PetscInt[],const PetscInt[]);

PETSC_EXTERN PetscErrorCode MatCreateShell(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,void *,Mat*);
PETSC_EXTERN PetscErrorCode MatCreateNormal(Mat,Mat*);
PETSC_EXTERN PetscErrorCode MatCreateNormalHermitian(Mat,Mat*);
PETSC_EXTERN PetscErrorCode MatCreateLRC(Mat,Mat,Vec,Mat,Mat*);
PETSC_EXTERN PetscErrorCode MatLRCGetMats(Mat,Mat*,Mat*,Vec*,Mat*);
PETSC_EXTERN PetscErrorCode MatCreateIS(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,ISLocalToGlobalMapping,ISLocalToGlobalMapping,Mat*);
PETSC_EXTERN PetscErrorCode MatCreateSeqAIJCRL(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);
PETSC_EXTERN PetscErrorCode MatCreateMPIAIJCRL(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);

PETSC_EXTERN PetscErrorCode MatCreateScatter(MPI_Comm,VecScatter,Mat*);
PETSC_EXTERN PetscErrorCode MatScatterSetVecScatter(Mat,VecScatter);
PETSC_EXTERN PetscErrorCode MatScatterGetVecScatter(Mat,VecScatter*);
PETSC_EXTERN PetscErrorCode MatCreateBlockMat(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt*,Mat*);
PETSC_EXTERN PetscErrorCode MatCompositeAddMat(Mat,Mat);
PETSC_EXTERN PetscErrorCode MatCompositeMerge(Mat);
PETSC_EXTERN PetscErrorCode MatCreateComposite(MPI_Comm,PetscInt,const Mat*,Mat*);
typedef enum {MAT_COMPOSITE_ADDITIVE,MAT_COMPOSITE_MULTIPLICATIVE} MatCompositeType;
PETSC_EXTERN PetscErrorCode MatCompositeSetType(Mat,MatCompositeType);

PETSC_EXTERN PetscErrorCode MatCreateFFT(MPI_Comm,PetscInt,const PetscInt[],MatType,Mat*);
PETSC_EXTERN PetscErrorCode MatCreateSeqCUFFT(MPI_Comm,PetscInt,const PetscInt[],Mat*);

PETSC_EXTERN PetscErrorCode MatCreateTranspose(Mat,Mat*);
PETSC_EXTERN PetscErrorCode MatTransposeGetMat(Mat,Mat*);
PETSC_EXTERN PetscErrorCode MatCreateHermitianTranspose(Mat,Mat*);
PETSC_EXTERN PetscErrorCode MatCreateSubMatrixVirtual(Mat,IS,IS,Mat*);
PETSC_EXTERN PetscErrorCode MatSubMatrixVirtualUpdate(Mat,Mat,IS,IS);
PETSC_EXTERN PetscErrorCode MatCreateLocalRef(Mat,IS,IS,Mat*);

#if defined(PETSC_HAVE_HYPRE)
PETSC_EXTERN PetscErrorCode MatHYPRESetPreallocation(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[]);
#endif

PETSC_EXTERN PetscErrorCode MatPythonSetType(Mat,const char[]);

PETSC_EXTERN PetscErrorCode MatResetPreallocation(Mat);
PETSC_EXTERN PetscErrorCode MatSetUp(Mat);
PETSC_EXTERN PetscErrorCode MatDestroy(Mat*);
PETSC_EXTERN PetscErrorCode MatGetNonzeroState(Mat,PetscObjectState*);

PETSC_EXTERN PetscErrorCode MatConjugate(Mat);
PETSC_EXTERN PetscErrorCode MatRealPart(Mat);
PETSC_EXTERN PetscErrorCode MatImaginaryPart(Mat);
PETSC_EXTERN PetscErrorCode MatGetDiagonalBlock(Mat,Mat*);
PETSC_EXTERN PetscErrorCode MatGetTrace(Mat,PetscScalar*);
PETSC_EXTERN PetscErrorCode MatInvertBlockDiagonal(Mat,const PetscScalar **);
PETSC_EXTERN PetscErrorCode MatInvertVariableBlockDiagonal(Mat,PetscInt,const PetscInt*,PetscScalar*);
PETSC_EXTERN PetscErrorCode MatInvertBlockDiagonalMat(Mat,Mat);

/* ------------------------------------------------------------*/
PETSC_EXTERN PetscErrorCode MatSetValues(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
PETSC_EXTERN PetscErrorCode MatSetValuesBlocked(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
PETSC_EXTERN PetscErrorCode MatSetValuesRow(Mat,PetscInt,const PetscScalar[]);
PETSC_EXTERN PetscErrorCode MatSetValuesRowLocal(Mat,PetscInt,const PetscScalar[]);
PETSC_EXTERN PetscErrorCode MatSetValuesBatch(Mat,PetscInt,PetscInt,PetscInt[],const PetscScalar[]);
PETSC_EXTERN PetscErrorCode MatSetRandom(Mat,PetscRandom);

/*S
     MatStencil - Data structure (C struct) for storing information about a single row or
        column of a matrix as indexed on an associated grid. These are arguments to MatSetStencil() and MatSetBlockStencil()

   The i,j, and k represent the logical coordinates over the entire grid (for 2 and 1 dimensional problems the k and j entries are ignored).
   The c represents the the degrees of freedom at each grid point (the dof argument to DMDASetDOF()). If dof is 1 then this entry is ignored.

   For stencil access to vectors see DMDAVecGetArray(), DMDAVecGetArrayF90().

   Fortran usage is different, see MatSetValuesStencil() for details.

   Level: beginner

  Concepts: matrix; linear operator

.seealso:  MatSetValuesStencil(), MatSetStencil(), MatSetValuesBlockedStencil(), DMDAVecGetArray(), DMDAVecGetArrayF90()
S*/
typedef struct {
  PetscInt k,j,i,c;
} MatStencil;

PETSC_EXTERN PetscErrorCode MatSetValuesStencil(Mat,PetscInt,const MatStencil[],PetscInt,const MatStencil[],const PetscScalar[],InsertMode);
PETSC_EXTERN PetscErrorCode MatSetValuesBlockedStencil(Mat,PetscInt,const MatStencil[],PetscInt,const MatStencil[],const PetscScalar[],InsertMode);
PETSC_EXTERN PetscErrorCode MatSetStencil(Mat,PetscInt,const PetscInt[],const PetscInt[],PetscInt);

/*E
    MatAssemblyType - Indicates if the matrix is now to be used, or if you plan
     to continue to add or insert values to it

    Level: beginner

.seealso: MatAssemblyBegin(), MatAssemblyEnd()
E*/
typedef enum {MAT_FLUSH_ASSEMBLY=1,MAT_FINAL_ASSEMBLY=0} MatAssemblyType;
PETSC_EXTERN PetscErrorCode MatAssemblyBegin(Mat,MatAssemblyType);
PETSC_EXTERN PetscErrorCode MatAssemblyEnd(Mat,MatAssemblyType);
PETSC_EXTERN PetscErrorCode MatAssembled(Mat,PetscBool *);



/*E
    MatOption - Options that may be set for a matrix and its behavior or storage

    Level: beginner

   Any additions/changes here MUST also be made in include/petsc/finclude/petscmat.h
   Any additions/changes here must also be made in src/mat/interface/dlregismat.c in MatOptions[]

   Developer Notes:
    Entries that are negative need not be called collectively by all processes.

.seealso: MatSetOption()
E*/
typedef enum {MAT_OPTION_MIN = -3,
              MAT_UNUSED_NONZERO_LOCATION_ERR = -2,
              MAT_ROW_ORIENTED = -1,
              MAT_SYMMETRIC = 1,
              MAT_STRUCTURALLY_SYMMETRIC = 2,
              MAT_NEW_DIAGONALS = 3,
              MAT_IGNORE_OFF_PROC_ENTRIES = 4,
              MAT_USE_HASH_TABLE = 5,
              MAT_KEEP_NONZERO_PATTERN = 6,
              MAT_IGNORE_ZERO_ENTRIES = 7,
              MAT_USE_INODES = 8,
              MAT_HERMITIAN = 9,
              MAT_SYMMETRY_ETERNAL = 10,
              MAT_NEW_NONZERO_LOCATION_ERR = 11,
              MAT_IGNORE_LOWER_TRIANGULAR = 12,
              MAT_ERROR_LOWER_TRIANGULAR = 13,
              MAT_GETROW_UPPERTRIANGULAR = 14,
              MAT_SPD = 15,
              MAT_NO_OFF_PROC_ZERO_ROWS = 16,
              MAT_NO_OFF_PROC_ENTRIES = 17,
              MAT_NEW_NONZERO_LOCATIONS = 18,
              MAT_NEW_NONZERO_ALLOCATION_ERR = 19,
              MAT_SUBSET_OFF_PROC_ENTRIES = 20,
              MAT_SUBMAT_SINGLEIS = 21,
              MAT_STRUCTURE_ONLY = 22,
              MAT_OPTION_MAX = 23} MatOption;

PETSC_EXTERN const char *const *MatOptions;
PETSC_EXTERN PetscErrorCode MatSetOption(Mat,MatOption,PetscBool);
PETSC_EXTERN PetscErrorCode MatGetOption(Mat,MatOption,PetscBool*);
PETSC_EXTERN PetscErrorCode MatGetType(Mat,MatType*);

PETSC_EXTERN PetscErrorCode MatGetValues(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],PetscScalar[]);
PETSC_EXTERN PetscErrorCode MatGetRow(Mat,PetscInt,PetscInt *,const PetscInt *[],const PetscScalar*[]);
PETSC_EXTERN PetscErrorCode MatRestoreRow(Mat,PetscInt,PetscInt *,const PetscInt *[],const PetscScalar*[]);
PETSC_EXTERN PetscErrorCode MatGetRowUpperTriangular(Mat);
PETSC_EXTERN PetscErrorCode MatRestoreRowUpperTriangular(Mat);
PETSC_EXTERN PetscErrorCode MatGetColumnVector(Mat,Vec,PetscInt);
PETSC_EXTERN PetscErrorCode MatSeqAIJGetArray(Mat,PetscScalar *[]);
PETSC_EXTERN PetscErrorCode MatSeqAIJRestoreArray(Mat,PetscScalar *[]);
PETSC_EXTERN PetscErrorCode MatSeqAIJGetMaxRowNonzeros(Mat,PetscInt*);
PETSC_EXTERN PetscErrorCode MatSeqAIJSetValuesLocalFast(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
PETSC_EXTERN PetscErrorCode MatSeqAIJSetType(Mat,MatType);
PETSC_EXTERN PetscErrorCode MatSeqAIJRegister(const char[],PetscErrorCode (*)(Mat,MatType,MatReuse,Mat *));
PETSC_EXTERN PetscFunctionList MatSeqAIJList;
PETSC_EXTERN PetscErrorCode MatSeqSBAIJGetArray(Mat,PetscScalar *[]);
PETSC_EXTERN PetscErrorCode MatSeqSBAIJRestoreArray(Mat,PetscScalar *[]);
PETSC_EXTERN PetscErrorCode MatDenseGetArray(Mat,PetscScalar *[]);
PETSC_EXTERN PetscErrorCode MatDenseRestoreArray(Mat,PetscScalar *[]);
PETSC_EXTERN PetscErrorCode MatDensePlaceArray(Mat,const PetscScalar[]);
PETSC_EXTERN PetscErrorCode MatDenseResetArray(Mat);
PETSC_EXTERN PetscErrorCode MatDenseGetArrayRead(Mat,const PetscScalar *[]);
PETSC_EXTERN PetscErrorCode MatDenseRestoreArrayRead(Mat,const PetscScalar *[]);
PETSC_EXTERN PetscErrorCode MatGetBlockSize(Mat,PetscInt *);
PETSC_EXTERN PetscErrorCode MatSetBlockSize(Mat,PetscInt);
PETSC_EXTERN PetscErrorCode MatGetBlockSizes(Mat,PetscInt *,PetscInt *);
PETSC_EXTERN PetscErrorCode MatSetBlockSizes(Mat,PetscInt,PetscInt);
PETSC_EXTERN PetscErrorCode MatSetBlockSizesFromMats(Mat,Mat,Mat);
PETSC_EXTERN PetscErrorCode MatSetVariableBlockSizes(Mat,PetscInt,PetscInt*);
PETSC_EXTERN PetscErrorCode MatGetVariableBlockSizes(Mat,PetscInt*,const PetscInt**);

PETSC_EXTERN PetscErrorCode MatDenseGetColumn(Mat,PetscInt,PetscScalar *[]);
PETSC_EXTERN PetscErrorCode MatDenseRestoreColumn(Mat,PetscScalar *[]);

PETSC_EXTERN PetscErrorCode MatMult(Mat,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatMultDiagonalBlock(Mat,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatMultAdd(Mat,Vec,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatMultTranspose(Mat,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatMultHermitianTranspose(Mat,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatIsTranspose(Mat,Mat,PetscReal,PetscBool *);
PETSC_EXTERN PetscErrorCode MatIsHermitianTranspose(Mat,Mat,PetscReal,PetscBool *);
PETSC_EXTERN PetscErrorCode MatMultTransposeAdd(Mat,Vec,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatMultHermitianTransposeAdd(Mat,Vec,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatMultConstrained(Mat,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatMultTransposeConstrained(Mat,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatMatSolve(Mat,Mat,Mat);
PETSC_EXTERN PetscErrorCode MatMatSolveTranspose(Mat,Mat,Mat);
PETSC_EXTERN PetscErrorCode MatMatTransposeSolve(Mat,Mat,Mat);
PETSC_EXTERN PetscErrorCode MatResidual(Mat,Vec,Vec,Vec);

/*E
    MatDuplicateOption - Indicates if a duplicated sparse matrix should have
  its numerical values copied over or just its nonzero structure.

    Level: beginner

   Any additions/changes here MUST also be made in include/petsc/finclude/petscmat.h

$   MAT_DO_NOT_COPY_VALUES    - Create a matrix using the same nonzero pattern as the original matrix,
$                               with zeros for the numerical values.
$   MAT_COPY_VALUES           - Create a matrix with the same nonzero pattern as the original matrix
$                               and with the same numerical values.
$   MAT_SHARE_NONZERO_PATTERN - Create a matrix that shares the nonzero structure with the previous matrix
$                               and does not copy it, using zeros for the numerical values. The parent and
$                               child matrices will share their index (i and j) arrays, and you cannot
$                               insert new nonzero entries into either matrix.

Notes:
    Many matrix types (including SeqAIJ) do not support the MAT_SHARE_NONZERO_PATTERN optimization; in
this case the behavior is as if MAT_DO_NOT_COPY_VALUES has been specified.

.seealso: MatDuplicate()
E*/
typedef enum {MAT_DO_NOT_COPY_VALUES,MAT_COPY_VALUES,MAT_SHARE_NONZERO_PATTERN} MatDuplicateOption;

PETSC_EXTERN PetscErrorCode MatConvert(Mat,MatType,MatReuse,Mat*);
PETSC_EXTERN PetscErrorCode MatDuplicate(Mat,MatDuplicateOption,Mat*);


PETSC_EXTERN PetscErrorCode MatCopy(Mat,Mat,MatStructure);
PETSC_EXTERN PetscErrorCode MatView(Mat,PetscViewer);
PETSC_EXTERN PetscErrorCode MatIsSymmetric(Mat,PetscReal,PetscBool *);
PETSC_EXTERN PetscErrorCode MatIsStructurallySymmetric(Mat,PetscBool *);
PETSC_EXTERN PetscErrorCode MatIsHermitian(Mat,PetscReal,PetscBool *);
PETSC_EXTERN PetscErrorCode MatIsSymmetricKnown(Mat,PetscBool *,PetscBool *);
PETSC_EXTERN PetscErrorCode MatIsHermitianKnown(Mat,PetscBool *,PetscBool *);
PETSC_EXTERN PetscErrorCode MatMissingDiagonal(Mat,PetscBool  *,PetscInt *);
PETSC_EXTERN PetscErrorCode MatLoad(Mat, PetscViewer);

PETSC_EXTERN PetscErrorCode MatGetRowIJ(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt*,const PetscInt *[],const PetscInt *[],PetscBool  *);
PETSC_EXTERN PetscErrorCode MatRestoreRowIJ(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt *,const PetscInt *[],const PetscInt *[],PetscBool  *);
PETSC_EXTERN PetscErrorCode MatGetColumnIJ(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt*,const PetscInt *[],const PetscInt *[],PetscBool  *);
PETSC_EXTERN PetscErrorCode MatRestoreColumnIJ(Mat,PetscInt,PetscBool ,PetscBool ,PetscInt *,const PetscInt *[],const PetscInt *[],PetscBool  *);

/*S
     MatInfo - Context of matrix information, used with MatGetInfo()

   In Fortran this is simply a double precision array of dimension MAT_INFO_SIZE

   Level: intermediate

  Concepts: matrix^nonzero information

.seealso:  MatGetInfo(), MatInfoType
S*/
typedef struct {
  PetscLogDouble block_size;                         /* block size */
  PetscLogDouble nz_allocated,nz_used,nz_unneeded;   /* number of nonzeros */
  PetscLogDouble memory;                             /* memory allocated */
  PetscLogDouble assemblies;                         /* number of matrix assemblies called */
  PetscLogDouble mallocs;                            /* number of mallocs during MatSetValues() */
  PetscLogDouble fill_ratio_given,fill_ratio_needed; /* fill ratio for LU/ILU */
  PetscLogDouble factor_mallocs;                     /* number of mallocs during factorization */
} MatInfo;

/*E
    MatInfoType - Indicates if you want information about the local part of the matrix,
     the entire parallel matrix or the maximum over all the local parts.

    Level: beginner

   Any additions/changes here MUST also be made in include/petsc/finclude/petscmat.h

.seealso: MatGetInfo(), MatInfo
E*/
typedef enum {MAT_LOCAL=1,MAT_GLOBAL_MAX=2,MAT_GLOBAL_SUM=3} MatInfoType;
PETSC_EXTERN PetscErrorCode MatGetInfo(Mat,MatInfoType,MatInfo*);
PETSC_EXTERN PetscErrorCode MatGetDiagonal(Mat,Vec);
PETSC_EXTERN PetscErrorCode MatGetRowMax(Mat,Vec,PetscInt[]);
PETSC_EXTERN PetscErrorCode MatGetRowMin(Mat,Vec,PetscInt[]);
PETSC_EXTERN PetscErrorCode MatGetRowMaxAbs(Mat,Vec,PetscInt[]);
PETSC_EXTERN PetscErrorCode MatGetRowMinAbs(Mat,Vec,PetscInt[]);
PETSC_EXTERN PetscErrorCode MatGetRowSum(Mat,Vec);
PETSC_EXTERN PetscErrorCode MatTranspose(Mat,MatReuse,Mat*);
PETSC_EXTERN PetscErrorCode MatHermitianTranspose(Mat,MatReuse,Mat*);
PETSC_EXTERN PetscErrorCode MatPermute(Mat,IS,IS,Mat*);
PETSC_EXTERN PetscErrorCode MatDiagonalScale(Mat,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatDiagonalSet(Mat,Vec,InsertMode);

PETSC_EXTERN PetscErrorCode MatEqual(Mat,Mat,PetscBool*);
PETSC_EXTERN PetscErrorCode MatMultEqual(Mat,Mat,PetscInt,PetscBool*);
PETSC_EXTERN PetscErrorCode MatMultAddEqual(Mat,Mat,PetscInt,PetscBool*);
PETSC_EXTERN PetscErrorCode MatMultTransposeEqual(Mat,Mat,PetscInt,PetscBool*);
PETSC_EXTERN PetscErrorCode MatMultTransposeAddEqual(Mat,Mat,PetscInt,PetscBool*);
PETSC_EXTERN PetscErrorCode MatMatMultEqual(Mat,Mat,Mat,PetscInt,PetscBool*);
PETSC_EXTERN PetscErrorCode MatTransposeMatMultEqual(Mat,Mat,Mat,PetscInt,PetscBool*);

PETSC_EXTERN PetscErrorCode MatNorm(Mat,NormType,PetscReal*);
PETSC_EXTERN PetscErrorCode MatGetColumnNorms(Mat,NormType,PetscReal*);
PETSC_EXTERN PetscErrorCode MatZeroEntries(Mat);
PETSC_EXTERN PetscErrorCode MatZeroRows(Mat,PetscInt,const PetscInt [],PetscScalar,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatZeroRowsIS(Mat,IS,PetscScalar,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatZeroRowsStencil(Mat,PetscInt,const MatStencil [],PetscScalar,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatZeroRowsColumnsStencil(Mat,PetscInt,const MatStencil[],PetscScalar,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatZeroRowsColumns(Mat,PetscInt,const PetscInt [],PetscScalar,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatZeroRowsColumnsIS(Mat,IS,PetscScalar,Vec,Vec);

PETSC_EXTERN PetscErrorCode MatGetSize(Mat,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode MatGetLocalSize(Mat,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode MatGetOwnershipRange(Mat,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode MatGetOwnershipRanges(Mat,const PetscInt**);
PETSC_EXTERN PetscErrorCode MatGetOwnershipRangeColumn(Mat,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode MatGetOwnershipRangesColumn(Mat,const PetscInt**);
PETSC_EXTERN PetscErrorCode MatGetOwnershipIS(Mat,IS*,IS*);

PETSC_EXTERN PetscErrorCode MatCreateSubMatrices(Mat,PetscInt,const IS[],const IS[],MatReuse,Mat *[]);
PETSC_DEPRECATED("Use MatCreateSubMatrices()") PETSC_STATIC_INLINE PetscErrorCode MatGetSubMatrices(Mat mat,PetscInt n,const IS irow[],const IS icol[],MatReuse scall,Mat *submat[]) {return MatCreateSubMatrices(mat,n,irow,icol,scall,submat);}
PETSC_EXTERN PetscErrorCode MatCreateSubMatricesMPI(Mat,PetscInt,const IS[],const IS[],MatReuse,Mat *[]);
PETSC_DEPRECATED("Use MatCreateSubMatricesMPI()") PETSC_STATIC_INLINE PetscErrorCode MatGetSubMatricesMPI(Mat mat,PetscInt n,const IS irow[],const IS icol[],MatReuse scall,Mat *submat[]) {return MatCreateSubMatricesMPI(mat,n,irow,icol,scall,submat);}
PETSC_EXTERN PetscErrorCode MatDestroyMatrices(PetscInt,Mat *[]);
PETSC_EXTERN PetscErrorCode MatDestroySubMatrices(PetscInt,Mat *[]);
PETSC_EXTERN PetscErrorCode MatCreateSubMatrix(Mat,IS,IS,MatReuse,Mat *);
PETSC_DEPRECATED("Use MatCreateSubMatrix()") PETSC_STATIC_INLINE PetscErrorCode MatGetSubMatrix(Mat mat,IS isrow,IS iscol,MatReuse cll,Mat *newmat) {return MatCreateSubMatrix(mat,isrow,iscol,cll,newmat);}
PETSC_EXTERN PetscErrorCode MatGetLocalSubMatrix(Mat,IS,IS,Mat*);
PETSC_EXTERN PetscErrorCode MatRestoreLocalSubMatrix(Mat,IS,IS,Mat*);
PETSC_EXTERN PetscErrorCode MatGetSeqNonzeroStructure(Mat,Mat*);
PETSC_EXTERN PetscErrorCode MatDestroySeqNonzeroStructure(Mat*);

PETSC_EXTERN PetscErrorCode MatCreateMPIAIJSumSeqAIJ(MPI_Comm,Mat,PetscInt,PetscInt,MatReuse,Mat*);
PETSC_EXTERN PetscErrorCode MatCreateMPIAIJSumSeqAIJSymbolic(MPI_Comm,Mat,PetscInt,PetscInt,Mat*);
PETSC_EXTERN PetscErrorCode MatCreateMPIAIJSumSeqAIJNumeric(Mat,Mat);
PETSC_EXTERN PetscErrorCode MatMPIAIJGetLocalMat(Mat,MatReuse,Mat*);
PETSC_EXTERN PetscErrorCode MatMPIAIJGetLocalMatCondensed(Mat,MatReuse,IS*,IS*,Mat*);
PETSC_EXTERN PetscErrorCode MatGetBrowsOfAcols(Mat,Mat,MatReuse,IS*,IS*,Mat*);
PETSC_EXTERN PetscErrorCode MatGetGhosts(Mat, PetscInt *,const PetscInt *[]);

PETSC_EXTERN PetscErrorCode MatIncreaseOverlap(Mat,PetscInt,IS[],PetscInt);
PETSC_EXTERN PetscErrorCode MatIncreaseOverlapSplit(Mat mat,PetscInt n,IS is[],PetscInt ov);
PETSC_EXTERN PetscErrorCode MatMPIAIJSetUseScalableIncreaseOverlap(Mat,PetscBool);

PETSC_EXTERN PetscErrorCode MatMatMult(Mat,Mat,MatReuse,PetscReal,Mat*);
PETSC_EXTERN PetscErrorCode MatMatMultSymbolic(Mat,Mat,PetscReal,Mat*);
PETSC_EXTERN PetscErrorCode MatMatMultNumeric(Mat,Mat,Mat);

PETSC_EXTERN PetscErrorCode MatMatMatMult(Mat,Mat,Mat,MatReuse,PetscReal,Mat*);
PETSC_EXTERN PetscErrorCode MatGalerkin(Mat,Mat,Mat,MatReuse,PetscReal,Mat*);

PETSC_EXTERN PetscErrorCode MatPtAP(Mat,Mat,MatReuse,PetscReal,Mat*);
PETSC_EXTERN PetscErrorCode MatPtAPSymbolic(Mat,Mat,PetscReal,Mat*);
PETSC_EXTERN PetscErrorCode MatPtAPNumeric(Mat,Mat,Mat);
PETSC_EXTERN PetscErrorCode MatRARt(Mat,Mat,MatReuse,PetscReal,Mat*);
PETSC_EXTERN PetscErrorCode MatRARtSymbolic(Mat,Mat,PetscReal,Mat*);
PETSC_EXTERN PetscErrorCode MatRARtNumeric(Mat,Mat,Mat);

PETSC_EXTERN PetscErrorCode MatTransposeMatMult(Mat,Mat,MatReuse,PetscReal,Mat*);
PETSC_EXTERN PetscErrorCode MatTransposetMatMultSymbolic(Mat,Mat,PetscReal,Mat*);
PETSC_EXTERN PetscErrorCode MatTransposetMatMultNumeric(Mat,Mat,Mat);
PETSC_EXTERN PetscErrorCode MatMatTransposeMult(Mat,Mat,MatReuse,PetscReal,Mat*);
PETSC_EXTERN PetscErrorCode MatMatTransposeMultSymbolic(Mat,Mat,PetscReal,Mat*);
PETSC_EXTERN PetscErrorCode MatMatTransposeMultNumeric(Mat,Mat,Mat);

PETSC_EXTERN PetscErrorCode MatAXPY(Mat,PetscScalar,Mat,MatStructure);
PETSC_EXTERN PetscErrorCode MatAYPX(Mat,PetscScalar,Mat,MatStructure);

PETSC_EXTERN PetscErrorCode MatScale(Mat,PetscScalar);
PETSC_EXTERN PetscErrorCode MatShift(Mat,PetscScalar);

PETSC_EXTERN PetscErrorCode MatSetLocalToGlobalMapping(Mat,ISLocalToGlobalMapping,ISLocalToGlobalMapping);
PETSC_EXTERN PetscErrorCode MatGetLocalToGlobalMapping(Mat,ISLocalToGlobalMapping*,ISLocalToGlobalMapping*);
PETSC_EXTERN PetscErrorCode MatGetLayouts(Mat,PetscLayout*,PetscLayout*);
PETSC_EXTERN PetscErrorCode MatZeroRowsLocal(Mat,PetscInt,const PetscInt [],PetscScalar,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatZeroRowsLocalIS(Mat,IS,PetscScalar,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatZeroRowsColumnsLocal(Mat,PetscInt,const PetscInt [],PetscScalar,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatZeroRowsColumnsLocalIS(Mat,IS,PetscScalar,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatSetValuesLocal(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
PETSC_EXTERN PetscErrorCode MatSetValuesBlockedLocal(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[],const PetscScalar[],InsertMode);

PETSC_EXTERN PetscErrorCode MatStashSetInitialSize(Mat,PetscInt,PetscInt);
PETSC_EXTERN PetscErrorCode MatStashGetInfo(Mat,PetscInt*,PetscInt*,PetscInt*,PetscInt*);

PETSC_EXTERN PetscErrorCode MatInterpolate(Mat,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatInterpolateAdd(Mat,Vec,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatRestrict(Mat,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatCreateVecs(Mat,Vec*,Vec*);
PETSC_DEPRECATED("Use MatCreateVecs()") PETSC_STATIC_INLINE PetscErrorCode MatGetVecs(Mat mat,Vec *x,Vec *y) {return MatCreateVecs(mat,x,y);}
PETSC_EXTERN PetscErrorCode MatCreateRedundantMatrix(Mat,PetscInt,MPI_Comm,MatReuse,Mat*);
PETSC_EXTERN PetscErrorCode MatGetMultiProcBlock(Mat,MPI_Comm,MatReuse,Mat*);
PETSC_EXTERN PetscErrorCode MatFindZeroDiagonals(Mat,IS*);
PETSC_EXTERN PetscErrorCode MatFindOffBlockDiagonalEntries(Mat,IS*);
PETSC_EXTERN PetscErrorCode MatCreateMPIMatConcatenateSeqMat(MPI_Comm,Mat,PetscInt,MatReuse,Mat*);

/*MC
   MatSetValue - Set a single entry into a matrix.

   Not collective

   Synopsis:
     #include <petscmat.h>
     PetscErrorCode MatSetValue(Mat m,PetscInt row,PetscInt col,PetscScalar value,InsertMode mode)

   Input Parameters:
+  m - the matrix
.  row - the row location of the entry
.  col - the column location of the entry
.  value - the value to insert
-  mode - either INSERT_VALUES or ADD_VALUES

   Notes:
   For efficiency one should use MatSetValues() and set several or many
   values simultaneously if possible.

   Level: beginner

.seealso: MatSetValues(), MatSetValueLocal()
M*/
PETSC_STATIC_INLINE PetscErrorCode MatSetValue(Mat v,PetscInt i,PetscInt j,PetscScalar va,InsertMode mode) {return MatSetValues(v,1,&i,1,&j,&va,mode);}

PETSC_STATIC_INLINE PetscErrorCode MatGetValue(Mat v,PetscInt i,PetscInt j,PetscScalar *va) {return MatGetValues(v,1,&i,1,&j,va);}

PETSC_STATIC_INLINE PetscErrorCode MatSetValueLocal(Mat v,PetscInt i,PetscInt j,PetscScalar va,InsertMode mode) {return MatSetValuesLocal(v,1,&i,1,&j,&va,mode);}

/*MC
   MatPreallocateInitialize - Begins the block of code that will count the number of nonzeros per
       row in a matrix providing the data that one can use to correctly preallocate the matrix.

   Synopsis:
   #include <petscmat.h>
   PetscErrorCode MatPreallocateInitialize(MPI_Comm comm, PetscInt nrows, PetscInt ncols, PetscInt *dnz, PetscInt *onz)

   Collective on MPI_Comm

   Input Parameters:
+  comm - the communicator that will share the eventually allocated matrix
.  nrows - the number of LOCAL rows in the matrix
-  ncols - the number of LOCAL columns in the matrix

   Output Parameters:
+  dnz - the array that will be passed to the matrix preallocation routines
-  onz - the other array passed to the matrix preallocation routines

   Level: intermediate

   Notes:
    See Users-Manual: ch_performance for more details.

   Do not malloc or free dnz and onz, that is handled internally by these routines

   Use MatPreallocateInitializeSymmetric() for symmetric matrices (MPISBAIJ matrices)

   This is a MACRO not a function because it has a leading { that is closed by PetscPreallocateFinalize().

  Concepts: preallocation^Matrix

.seealso: MatPreallocateFinalize(), MatPreallocateSet(), MatPreallocateSymmetricSetBlock(), MatPreallocateSetLocal(),
          MatPreallocateInitializeSymmetric(), MatPreallocateSymmetricSetLocalBlock()
M*/
#define MatPreallocateInitialize(comm,nrows,ncols,dnz,onz) 0; \
{ \
  PetscErrorCode _4_ierr; PetscInt __nrows = (nrows),__ncols = (ncols),__rstart,__start,__end; \
  _4_ierr = PetscCalloc2((size_t)__nrows,&dnz,(size_t)__nrows,&onz);CHKERRQ(_4_ierr); \
  __start = 0; __end = __start;                                         \
  _4_ierr = MPI_Scan(&__ncols,&__end,1,MPIU_INT,MPI_SUM,comm);CHKERRQ(_4_ierr); __start = __end - __ncols;\
  _4_ierr = MPI_Scan(&__nrows,&__rstart,1,MPIU_INT,MPI_SUM,comm);CHKERRQ(_4_ierr); __rstart = __rstart - __nrows;

/*MC
   MatPreallocateSetLocal - Indicates the locations (rows and columns) in the matrix where nonzeros will be
       inserted using a local number of the rows and columns

   Synopsis:
   #include <petscmat.h>
   PetscErrorCode MatPreallocateSetLocal(ISLocalToGlobalMappping map,PetscInt nrows, PetscInt *rows,PetscInt ncols, PetscInt *cols,PetscInt *dnz, PetscInt *onz)

   Not Collective

   Input Parameters:
+  map - the row mapping from local numbering to global numbering
.  nrows - the number of rows indicated
.  rows - the indices of the rows
.  cmap - the column mapping from local to global numbering
.  ncols - the number of columns in the matrix
.  cols - the columns indicated
.  dnz - the array that will be passed to the matrix preallocation routines
-  onz - the other array passed to the matrix preallocation routines

   Level: intermediate

   Notes:
    See Users-Manual: ch_performance for more details.

   Do not malloc or free dnz and onz, that is handled internally by these routines

  Concepts: preallocation^Matrix

.seealso: MatPreallocateFinalize(), MatPreallocateSet(), MatPreallocateSymmetricSetBlock()
          MatPreallocateInitialize(), MatPreallocateSymmetricSetLocalBlock(), MatPreallocateSetLocalRemoveDups()
M*/
#define MatPreallocateSetLocal(rmap,nrows,rows,cmap,ncols,cols,dnz,onz) 0; \
{\
  PetscInt __l;\
  _4_ierr = ISLocalToGlobalMappingApply(rmap,nrows,rows,rows);CHKERRQ(_4_ierr);\
  _4_ierr = ISLocalToGlobalMappingApply(cmap,ncols,cols,cols);CHKERRQ(_4_ierr);\
  for (__l=0;__l<nrows;__l++) {\
    _4_ierr = MatPreallocateSet((rows)[__l],ncols,cols,dnz,onz);CHKERRQ(_4_ierr);\
  }\
}

/*MC
   MatPreallocateSetLocalRemoveDups - Indicates the locations (rows and columns) in the matrix where nonzeros will be
       inserted using a local number of the rows and columns. This version removes any duplicate columns in cols

   Synopsis:
   #include <petscmat.h>
   PetscErrorCode MatPreallocateSetLocalRemoveDups(ISLocalToGlobalMappping map,PetscInt nrows, PetscInt *rows,PetscInt ncols, PetscInt *cols,PetscInt *dnz, PetscInt *onz)

   Not Collective

   Input Parameters:
+  map - the row mapping from local numbering to global numbering
.  nrows - the number of rows indicated
.  rows - the indices of the rows (these values are mapped to the global values)
.  cmap - the column mapping from local to global numbering
.  ncols - the number of columns in the matrix   (this value will be changed if duplicate columns are found)
.  cols - the columns indicated (these values are mapped to the global values, they are then sorted and duplicates removed)
.  dnz - the array that will be passed to the matrix preallocation routines
-  onz - the other array passed to the matrix preallocation routines

   Level: intermediate

   Notes:
    See Users-Manual: ch_performance for more details.

   Do not malloc or free dnz and onz, that is handled internally by these routines

  Concepts: preallocation^Matrix

.seealso: MatPreallocateFinalize(), MatPreallocateSet(), MatPreallocateSymmetricSetBlock()
          MatPreallocateInitialize(), MatPreallocateSymmetricSetLocalBlock(), MatPreallocateSetLocal()
M*/
#define MatPreallocateSetLocalRemoveDups(rmap,nrows,rows,cmap,ncols,cols,dnz,onz) 0; \
{\
  PetscInt __l;\
  _4_ierr = ISLocalToGlobalMappingApply(rmap,nrows,rows,rows);CHKERRQ(_4_ierr);\
  _4_ierr = ISLocalToGlobalMappingApply(cmap,ncols,cols,cols);CHKERRQ(_4_ierr);\
  _4_ierr = PetscSortRemoveDupsInt(&ncols,cols);CHKERRQ(_4_ierr);\
  for (__l=0;__l<nrows;__l++) {\
    _4_ierr = MatPreallocateSet((rows)[__l],ncols,cols,dnz,onz);CHKERRQ(_4_ierr);\
  }\
}

/*MC
   MatPreallocateSetLocalBlock - Indicates the locations (rows and columns) in the matrix where nonzeros will be
       inserted using a local number of the rows and columns

   Synopsis:
   #include <petscmat.h>
   PetscErrorCode MatPreallocateSetLocalBlock(ISLocalToGlobalMappping map,PetscInt nrows, PetscInt *rows,PetscInt ncols, PetscInt *cols,PetscInt *dnz, PetscInt *onz)

   Not Collective

   Input Parameters:
+  map - the row mapping from local numbering to global numbering
.  nrows - the number of rows indicated
.  rows - the indices of the rows
.  cmap - the column mapping from local to global numbering
.  ncols - the number of columns in the matrix
.  cols - the columns indicated
.  dnz - the array that will be passed to the matrix preallocation routines
-  onz - the other array passed to the matrix preallocation routines

   Level: intermediate

   Notes:
    See Users-Manual: ch_performance for more details.

   Do not malloc or free dnz and onz, that is handled internally by these routines

  Concepts: preallocation^Matrix

.seealso: MatPreallocateFinalize(), MatPreallocateSet(), MatPreallocateSymmetricSetBlock()
          MatPreallocateInitialize(), MatPreallocateSymmetricSetLocalBlock()
M*/
#define MatPreallocateSetLocalBlock(rmap,nrows,rows,cmap,ncols,cols,dnz,onz) 0; \
{\
  PetscInt __l;\
  _4_ierr = ISLocalToGlobalMappingApplyBlock(rmap,nrows,rows,rows);CHKERRQ(_4_ierr);\
  _4_ierr = ISLocalToGlobalMappingApplyBlock(cmap,ncols,cols,cols);CHKERRQ(_4_ierr);\
  for (__l=0;__l<nrows;__l++) {\
    _4_ierr = MatPreallocateSet((rows)[__l],ncols,cols,dnz,onz);CHKERRQ(_4_ierr);\
  }\
}

/*MC
   MatPreallocateSymmetricSetLocalBlock - Indicates the locations (rows and columns) in the matrix where nonzeros will be
       inserted using a local number of the rows and columns

   Synopsis:
   #include <petscmat.h>
   PetscErrorCode MatPreallocateSymmetricSetLocalBlock(ISLocalToGlobalMappping map,PetscInt nrows, PetscInt *rows,PetscInt ncols, PetscInt *cols,PetscInt *dnz, PetscInt *onz)

   Not Collective

   Input Parameters:
+  map - the mapping between local numbering and global numbering
.  nrows - the number of rows indicated
.  rows - the indices of the rows
.  ncols - the number of columns in the matrix
.  cols - the columns indicated
.  dnz - the array that will be passed to the matrix preallocation routines
-  onz - the other array passed to the matrix preallocation routines

   Level: intermediate

   Notes:
    See Users-Manual: ch_performance for more details.

   Do not malloc or free dnz and onz that is handled internally by these routines

  Concepts: preallocation^Matrix

.seealso: MatPreallocateFinalize(), MatPreallocateSet()
          MatPreallocateInitialize(),  MatPreallocateSetLocal()
M*/
#define MatPreallocateSymmetricSetLocalBlock(map,nrows,rows,ncols,cols,dnz,onz) 0;\
{\
  PetscInt __l;\
  _4_ierr = ISLocalToGlobalMappingApplyBlock(map,nrows,rows,rows);CHKERRQ(_4_ierr);\
  _4_ierr = ISLocalToGlobalMappingApplyBlock(map,ncols,cols,cols);CHKERRQ(_4_ierr);\
  for (__l=0;__l<nrows;__l++) {\
    _4_ierr = MatPreallocateSymmetricSetBlock((rows)[__l],ncols,cols,dnz,onz);CHKERRQ(_4_ierr);\
  }\
}
/*MC
   MatPreallocateSet - Indicates the locations (rows and columns) in the matrix where nonzeros will be
       inserted using a local number of the rows and columns

   Synopsis:
   #include <petscmat.h>
   PetscErrorCode MatPreallocateSet(PetscInt nrows, PetscInt *rows,PetscInt ncols, PetscInt *cols,PetscInt *dnz, PetscInt *onz)

   Not Collective

   Input Parameters:
+  row - the row
.  ncols - the number of columns in the matrix
-  cols - the columns indicated

   Output Parameters:
+  dnz - the array that will be passed to the matrix preallocation routines
-  onz - the other array passed to the matrix preallocation routines

   Level: intermediate

   Notes:
    See Users-Manual: ch_performance for more details.

   Do not malloc or free dnz and onz that is handled internally by these routines

   This is a MACRO not a function because it uses variables declared in MatPreallocateInitialize().

  Concepts: preallocation^Matrix

.seealso: MatPreallocateFinalize(), MatPreallocateSet(), MatPreallocateSymmetricSetBlock()
          MatPreallocateInitialize(), MatPreallocateSetLocal()
M*/
#define MatPreallocateSet(row,nc,cols,dnz,onz) 0;\
{ PetscInt __i; \
  if (row < __rstart) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Trying to set preallocation for row %D less than first local row %D",row,__rstart);\
  if (row >= __rstart+__nrows) SETERRQ2(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Trying to set preallocation for row %D greater than last local row %D",row,__rstart+__nrows-1);\
  for (__i=0; __i<nc; __i++) {\
    if ((cols)[__i] < __start || (cols)[__i] >= __end) onz[row - __rstart]++; \
    else if (dnz[row - __rstart] < __ncols) dnz[row - __rstart]++;\
  }\
}

/*MC
   MatPreallocateSymmetricSetBlock - Indicates the locations (rows and columns) in the matrix where nonzeros will be
       inserted using a local number of the rows and columns

   Synopsis:
   #include <petscmat.h>
   PetscErrorCode MatPreallocateSymmetricSetBlock(PetscInt nrows, PetscInt *rows,PetscInt ncols, PetscInt *cols,PetscInt *dnz, PetscInt *onz)

   Not Collective

   Input Parameters:
+  nrows - the number of rows indicated
.  rows - the indices of the rows
.  ncols - the number of columns in the matrix
.  cols - the columns indicated
.  dnz - the array that will be passed to the matrix preallocation routines
-  onz - the other array passed to the matrix preallocation routines

   Level: intermediate

   Notes:
    See Users-Manual: ch_performance for more details.

   Do not malloc or free dnz and onz that is handled internally by these routines

   This is a MACRO not a function because it uses variables declared in MatPreallocateInitialize().

  Concepts: preallocation^Matrix

.seealso: MatPreallocateFinalize(), MatPreallocateSet(),  MatPreallocateInitialize(),
          MatPreallocateSymmetricSetLocalBlock(), MatPreallocateSetLocal()
M*/
#define MatPreallocateSymmetricSetBlock(row,nc,cols,dnz,onz) 0;\
{ PetscInt __i; \
  for (__i=0; __i<nc; __i++) {\
    if (cols[__i] >= __end) onz[row - __rstart]++; \
    else if (cols[__i] >= row && dnz[row - __rstart] < __ncols) dnz[row - __rstart]++;\
  }\
}

/*MC
   MatPreallocateLocation -  An alternative to MatPreallocateSet() that puts the nonzero locations into the matrix if it exists

   Synopsis:
   #include <petscmat.h>
   PetscErrorCode MatPreallocateLocations(Mat A,PetscInt row,PetscInt ncols,PetscInt *cols,PetscInt *dnz,PetscInt *onz)

   Not Collective

   Input Parameters:
.  A - matrix
.  row - row where values exist (must be local to this process)
.  ncols - number of columns
.  cols - columns with nonzeros
.  dnz - the array that will be passed to the matrix preallocation routines
-  onz - the other array passed to the matrix preallocation routines

   Level: intermediate

   Notes:
    See Users-Manual: ch_performance for more details.

   Do not malloc or free dnz and onz that is handled internally by these routines

   This is a MACRO not a function because it uses a bunch of variables private to the MatPreallocation.... routines.

  Concepts: preallocation^Matrix

.seealso: MatPreallocateInitialize(), MatPreallocateSet(), MatPreallocateSymmetricSetBlock(), MatPreallocateSetLocal(),
          MatPreallocateSymmetricSetLocalBlock()
M*/
#define MatPreallocateLocation(A,row,ncols,cols,dnz,onz) 0;if (A) {ierr = MatSetValues(A,1,&row,ncols,cols,NULL,INSERT_VALUES);CHKERRQ(ierr);} else {ierr =  MatPreallocateSet(row,ncols,cols,dnz,onz);CHKERRQ(ierr);}


/*MC
   MatPreallocateFinalize - Ends the block of code that will count the number of nonzeros per
       row in a matrix providing the data that one can use to correctly preallocate the matrix.

   Synopsis:
   #include <petscmat.h>
   PetscErrorCode MatPreallocateFinalize(PetscInt *dnz, PetscInt *onz)

   Collective on MPI_Comm

   Input Parameters:
+  dnz - the array that was be passed to the matrix preallocation routines
-  onz - the other array passed to the matrix preallocation routines

   Level: intermediate

   Notes:
    See Users-Manual: ch_performance for more details.

   Do not malloc or free dnz and onz that is handled internally by these routines

   This is a MACRO not a function because it closes the { started in MatPreallocateInitialize().

  Concepts: preallocation^Matrix

.seealso: MatPreallocateInitialize(), MatPreallocateSet(), MatPreallocateSymmetricSetBlock(), MatPreallocateSetLocal(),
          MatPreallocateSymmetricSetLocalBlock()
M*/
#define MatPreallocateFinalize(dnz,onz) 0;_4_ierr = PetscFree2(dnz,onz);CHKERRQ(_4_ierr);}

/* Routines unique to particular data structures */
PETSC_EXTERN PetscErrorCode MatShellGetContext(Mat,void *);

PETSC_EXTERN PetscErrorCode MatInodeAdjustForInodes(Mat,IS*,IS*);
PETSC_EXTERN PetscErrorCode MatInodeGetInodeSizes(Mat,PetscInt *,PetscInt *[],PetscInt *);

PETSC_EXTERN PetscErrorCode MatSeqAIJSetColumnIndices(Mat,PetscInt[]);
PETSC_EXTERN PetscErrorCode MatSeqBAIJSetColumnIndices(Mat,PetscInt[]);
PETSC_EXTERN PetscErrorCode MatCreateSeqAIJWithArrays(MPI_Comm,PetscInt,PetscInt,PetscInt[],PetscInt[],PetscScalar[],Mat*);
PETSC_EXTERN PetscErrorCode MatCreateSeqBAIJWithArrays(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt[],PetscInt[],PetscScalar[],Mat*);
PETSC_EXTERN PetscErrorCode MatCreateSeqSBAIJWithArrays(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt[],PetscInt[],PetscScalar[],Mat*);
PETSC_EXTERN PetscErrorCode MatCreateSeqAIJFromTriple(MPI_Comm,PetscInt,PetscInt,PetscInt[],PetscInt[],PetscScalar[],Mat*,PetscInt,PetscBool);

#define MAT_SKIP_ALLOCATION -4

PETSC_EXTERN PetscErrorCode MatSeqBAIJSetPreallocation(Mat,PetscInt,PetscInt,const PetscInt[]);
PETSC_EXTERN PetscErrorCode MatSeqSBAIJSetPreallocation(Mat,PetscInt,PetscInt,const PetscInt[]);
PETSC_EXTERN PetscErrorCode MatSeqAIJSetPreallocation(Mat,PetscInt,const PetscInt[]);

PETSC_EXTERN PetscErrorCode MatMPIBAIJSetPreallocation(Mat,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[]);
PETSC_EXTERN PetscErrorCode MatMPISBAIJSetPreallocation(Mat,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[]);
PETSC_EXTERN PetscErrorCode MatMPIAIJSetPreallocation(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[]);
PETSC_EXTERN PetscErrorCode MatSeqAIJSetPreallocationCSR(Mat,const PetscInt [],const PetscInt [],const PetscScalar []);
PETSC_EXTERN PetscErrorCode MatSeqBAIJSetPreallocationCSR(Mat,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[]);
PETSC_EXTERN PetscErrorCode MatMPIAIJSetPreallocationCSR(Mat,const PetscInt[],const PetscInt[],const PetscScalar[]);
PETSC_EXTERN PetscErrorCode MatMPIBAIJSetPreallocationCSR(Mat,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[]);
PETSC_EXTERN PetscErrorCode MatMPIAdjSetPreallocation(Mat,PetscInt[],PetscInt[],PetscInt[]);
PETSC_EXTERN PetscErrorCode MatMPIAdjToSeq(Mat,Mat*);
PETSC_EXTERN PetscErrorCode MatMPIDenseSetPreallocation(Mat,PetscScalar[]);
PETSC_EXTERN PetscErrorCode MatSeqDenseSetPreallocation(Mat,PetscScalar[]);
PETSC_EXTERN PetscErrorCode MatMPIAIJGetSeqAIJ(Mat,Mat*,Mat*,const PetscInt*[]);
PETSC_EXTERN PetscErrorCode MatMPIBAIJGetSeqBAIJ(Mat,Mat*,Mat*,const PetscInt*[]);
PETSC_EXTERN PetscErrorCode MatMPIAdjCreateNonemptySubcommMat(Mat,Mat*);

PETSC_EXTERN PetscErrorCode MatSeqDenseSetLDA(Mat,PetscInt);
PETSC_EXTERN PetscErrorCode MatDenseGetLocalMatrix(Mat,Mat*);

PETSC_EXTERN PetscErrorCode MatStoreValues(Mat);
PETSC_EXTERN PetscErrorCode MatRetrieveValues(Mat);

PETSC_EXTERN PetscErrorCode MatDAADSetCtx(Mat,void*);

PETSC_EXTERN PetscErrorCode MatFindNonzeroRows(Mat,IS*);
PETSC_EXTERN PetscErrorCode MatFindZeroRows(Mat,IS*);
/*
  These routines are not usually accessed directly, rather solving is
  done through the KSP and PC interfaces.
*/

/*J
    MatOrderingType - String with the name of a PETSc matrix ordering

   Level: beginner

.seealso: MatGetOrdering()
J*/
typedef const char* MatOrderingType;
#define MATORDERINGNATURAL     "natural"
#define MATORDERINGND          "nd"
#define MATORDERING1WD         "1wd"
#define MATORDERINGRCM         "rcm"
#define MATORDERINGQMD         "qmd"
#define MATORDERINGROWLENGTH   "rowlength"
#define MATORDERINGWBM         "wbm"
#define MATORDERINGSPECTRAL    "spectral"
#define MATORDERINGAMD         "amd"            /* only works if UMFPACK is installed with PETSc */

PETSC_EXTERN PetscErrorCode MatGetOrdering(Mat,MatOrderingType,IS*,IS*);
PETSC_EXTERN PetscErrorCode MatGetOrderingList(PetscFunctionList*);
PETSC_EXTERN PetscErrorCode MatOrderingRegister(const char[],PetscErrorCode(*)(Mat,MatOrderingType,IS*,IS*));
PETSC_EXTERN PetscFunctionList MatOrderingList;

PETSC_EXTERN PetscErrorCode MatReorderForNonzeroDiagonal(Mat,PetscReal,IS,IS);
PETSC_EXTERN PetscErrorCode MatCreateLaplacian(Mat,PetscReal,PetscBool,Mat*);

/*S
    MatFactorShiftType - Numeric Shift.

   Level: beginner

S*/
typedef enum {MAT_SHIFT_NONE,MAT_SHIFT_NONZERO,MAT_SHIFT_POSITIVE_DEFINITE,MAT_SHIFT_INBLOCKS} MatFactorShiftType;
PETSC_EXTERN const char *const MatFactorShiftTypes[];
PETSC_EXTERN const char *const MatFactorShiftTypesDetail[];

/*S
    MatFactorError - indicates what type of error in matrix factor

    Level: beginner

    Developer Notes:
    Any additions/changes here MUST also be made in include/petsc/finclude/petscmat.h

.seealso: MatGetFactor()
S*/
typedef enum {MAT_FACTOR_NOERROR,MAT_FACTOR_STRUCT_ZEROPIVOT,MAT_FACTOR_NUMERIC_ZEROPIVOT,MAT_FACTOR_OUTMEMORY,MAT_FACTOR_OTHER} MatFactorError;

PETSC_EXTERN PetscErrorCode MatFactorGetError(Mat,MatFactorError*);
PETSC_EXTERN PetscErrorCode MatFactorClearError(Mat);
PETSC_EXTERN PetscErrorCode MatFactorGetErrorZeroPivot(Mat,PetscReal*,PetscInt*);

/*S
   MatFactorInfo - Data passed into the matrix factorization routines, and information about the resulting factorization

   In Fortran these are simply double precision arrays of size MAT_FACTORINFO_SIZE, that is use
$     MatFactorInfo  info(MAT_FACTORINFO_SIZE)

   Notes:
    These are not usually directly used by users, instead use PC type of LU, ILU, CHOLESKY or ICC.

      You can use MatFactorInfoInitialize() to set default values.

   Level: developer

.seealso: MatLUFactorSymbolic(), MatILUFactorSymbolic(), MatCholeskyFactorSymbolic(), MatICCFactorSymbolic(), MatICCFactor(),
          MatFactorInfoInitialize()

S*/
typedef struct {
  PetscReal     diagonal_fill;  /* force diagonal to fill in if initially not filled */
  PetscReal     usedt;
  PetscReal     dt;             /* drop tolerance */
  PetscReal     dtcol;          /* tolerance for pivoting */
  PetscReal     dtcount;        /* maximum nonzeros to be allowed per row */
  PetscReal     fill;           /* expected fill, nonzeros in factored matrix/nonzeros in original matrix */
  PetscReal     levels;         /* ICC/ILU(levels) */
  PetscReal     pivotinblocks;  /* for BAIJ and SBAIJ matrices pivot in factorization on blocks, default 1.0
                                   factorization may be faster if do not pivot */
  PetscReal     zeropivot;      /* pivot is called zero if less than this */
  PetscReal     shifttype;      /* type of shift added to matrix factor to prevent zero pivots */
  PetscReal     shiftamount;     /* how large the shift is */
} MatFactorInfo;

PETSC_EXTERN PetscErrorCode MatFactorInfoInitialize(MatFactorInfo*);
PETSC_EXTERN PetscErrorCode MatCholeskyFactor(Mat,IS,const MatFactorInfo*);
PETSC_EXTERN PetscErrorCode MatCholeskyFactorSymbolic(Mat,Mat,IS,const MatFactorInfo*);
PETSC_EXTERN PetscErrorCode MatCholeskyFactorNumeric(Mat,Mat,const MatFactorInfo*);
PETSC_EXTERN PetscErrorCode MatLUFactor(Mat,IS,IS,const MatFactorInfo*);
PETSC_EXTERN PetscErrorCode MatILUFactor(Mat,IS,IS,const MatFactorInfo*);
PETSC_EXTERN PetscErrorCode MatLUFactorSymbolic(Mat,Mat,IS,IS,const MatFactorInfo*);
PETSC_EXTERN PetscErrorCode MatILUFactorSymbolic(Mat,Mat,IS,IS,const MatFactorInfo*);
PETSC_EXTERN PetscErrorCode MatICCFactorSymbolic(Mat,Mat,IS,const MatFactorInfo*);
PETSC_EXTERN PetscErrorCode MatICCFactor(Mat,IS,const MatFactorInfo*);
PETSC_EXTERN PetscErrorCode MatLUFactorNumeric(Mat,Mat,const MatFactorInfo*);
PETSC_EXTERN PetscErrorCode MatGetInertia(Mat,PetscInt*,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode MatSolve(Mat,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatForwardSolve(Mat,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatBackwardSolve(Mat,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatSolveAdd(Mat,Vec,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatSolveTranspose(Mat,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatSolveTransposeAdd(Mat,Vec,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatSolves(Mat,Vecs,Vecs);
PETSC_EXTERN PetscErrorCode MatSetUnfactored(Mat);

typedef enum {MAT_FACTOR_SCHUR_UNFACTORED, MAT_FACTOR_SCHUR_FACTORED, MAT_FACTOR_SCHUR_INVERTED} MatFactorSchurStatus;
PETSC_EXTERN PetscErrorCode MatFactorSetSchurIS(Mat,IS);
PETSC_EXTERN PetscErrorCode MatFactorGetSchurComplement(Mat,Mat*,MatFactorSchurStatus*);
PETSC_EXTERN PetscErrorCode MatFactorRestoreSchurComplement(Mat,Mat*,MatFactorSchurStatus);
PETSC_EXTERN PetscErrorCode MatFactorInvertSchurComplement(Mat);
PETSC_EXTERN PetscErrorCode MatFactorCreateSchurComplement(Mat,Mat*,MatFactorSchurStatus*);
PETSC_EXTERN PetscErrorCode MatFactorSolveSchurComplement(Mat,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatFactorSolveSchurComplementTranspose(Mat,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatFactorFactorizeSchurComplement(Mat);

/*E
    MatSORType - What type of (S)SOR to perform

    Level: beginner

   May be bitwise ORd together

   Any additions/changes here MUST also be made in include/petsc/finclude/petscmat.h

   MatSORType may be bitwise ORd together, so do not change the numbers

.seealso: MatSOR()
E*/
typedef enum {SOR_FORWARD_SWEEP=1,SOR_BACKWARD_SWEEP=2,SOR_SYMMETRIC_SWEEP=3,
              SOR_LOCAL_FORWARD_SWEEP=4,SOR_LOCAL_BACKWARD_SWEEP=8,
              SOR_LOCAL_SYMMETRIC_SWEEP=12,SOR_ZERO_INITIAL_GUESS=16,
              SOR_EISENSTAT=32,SOR_APPLY_UPPER=64,SOR_APPLY_LOWER=128} MatSORType;
PETSC_EXTERN PetscErrorCode MatSOR(Mat,Vec,PetscReal,MatSORType,PetscReal,PetscInt,PetscInt,Vec);

/*
    These routines are for efficiently computing Jacobians via finite differences.
*/

/*S
     MatColoring - Object for managing the coloring of matrices.

   Level: beginner

  Concepts: matrix, coloring

.seealso:  MatFDColoringCreate() ISColoring MatFDColoring
S*/
typedef struct _p_MatColoring* MatColoring;
/*J
    MatColoringType - String with the name of a PETSc matrix coloring

   Level: beginner

.seealso: MatColoringSetType(), MatColoring
J*/

typedef const  char*           MatColoringType;
#define MATCOLORINGJP      "jp"
#define MATCOLORINGPOWER   "power"
#define MATCOLORINGNATURAL "natural"
#define MATCOLORINGSL      "sl"
#define MATCOLORINGLF      "lf"
#define MATCOLORINGID      "id"
#define MATCOLORINGGREEDY  "greedy"

/*E
   MatColoringWeightType - Type of weight scheme

    Not Collective

+   MAT_COLORING_RANDOM  - Random weights
.   MAT_COLORING_LEXICAL - Lexical weighting based upon global numbering.
-   MAT_COLORING_LF      - Last-first weighting.

    Level: intermediate

   Any additions/changes here MUST also be made in include/petsc/finclude/petscmat.h
E*/
typedef enum {MAT_COLORING_WEIGHT_RANDOM,MAT_COLORING_WEIGHT_LEXICAL,MAT_COLORING_WEIGHT_LF,MAT_COLORING_WEIGHT_SL} MatColoringWeightType;

PETSC_EXTERN PetscErrorCode MatColoringCreate(Mat,MatColoring*);
PETSC_EXTERN PetscErrorCode MatColoringGetDegrees(Mat,PetscInt,PetscInt*);
PETSC_EXTERN PetscErrorCode MatColoringDestroy(MatColoring*);
PETSC_EXTERN PetscErrorCode MatColoringView(MatColoring,PetscViewer);
PETSC_EXTERN PetscErrorCode MatColoringSetType(MatColoring,MatColoringType);
PETSC_EXTERN PetscErrorCode MatColoringSetFromOptions(MatColoring);
PETSC_EXTERN PetscErrorCode MatColoringSetDistance(MatColoring,PetscInt);
PETSC_EXTERN PetscErrorCode MatColoringGetDistance(MatColoring,PetscInt*);
PETSC_EXTERN PetscErrorCode MatColoringSetMaxColors(MatColoring,PetscInt);
PETSC_EXTERN PetscErrorCode MatColoringGetMaxColors(MatColoring,PetscInt*);
PETSC_EXTERN PetscErrorCode MatColoringApply(MatColoring,ISColoring*);
PETSC_EXTERN PetscErrorCode MatColoringRegister(const char[],PetscErrorCode(*)(MatColoring));
PETSC_EXTERN PetscErrorCode MatColoringPatch(Mat,PetscInt,PetscInt,ISColoringValue[],ISColoring*);
PETSC_EXTERN PetscErrorCode MatColoringSetWeightType(MatColoring,MatColoringWeightType);
PETSC_EXTERN PetscErrorCode MatColoringSetWeights(MatColoring,PetscReal*,PetscInt*);
PETSC_EXTERN PetscErrorCode MatColoringCreateWeights(MatColoring,PetscReal **,PetscInt **lperm);
PETSC_EXTERN PetscErrorCode MatColoringTest(MatColoring,ISColoring);
PETSC_DEPRECATED("Use MatColoringTest()") PETSC_STATIC_INLINE PetscErrorCode MatColoringTestValid(MatColoring matcoloring,ISColoring iscoloring) {return MatColoringTest(matcoloring,iscoloring);}
PETSC_EXTERN PetscErrorCode MatISColoringTest(Mat,ISColoring);

/*S
     MatFDColoring - Object for computing a sparse Jacobian via finite differences
        and coloring

   Level: beginner

  Concepts: coloring, sparse Jacobian, finite differences

.seealso:  MatFDColoringCreate()
S*/
typedef struct _p_MatFDColoring* MatFDColoring;

PETSC_EXTERN PetscErrorCode MatFDColoringCreate(Mat,ISColoring,MatFDColoring *);
PETSC_EXTERN PetscErrorCode MatFDColoringDestroy(MatFDColoring*);
PETSC_EXTERN PetscErrorCode MatFDColoringView(MatFDColoring,PetscViewer);
PETSC_EXTERN PetscErrorCode MatFDColoringSetFunction(MatFDColoring,PetscErrorCode (*)(void),void*);
PETSC_EXTERN PetscErrorCode MatFDColoringGetFunction(MatFDColoring,PetscErrorCode (**)(void),void**);
PETSC_EXTERN PetscErrorCode MatFDColoringSetParameters(MatFDColoring,PetscReal,PetscReal);
PETSC_EXTERN PetscErrorCode MatFDColoringSetFromOptions(MatFDColoring);
PETSC_EXTERN PetscErrorCode MatFDColoringApply(Mat,MatFDColoring,Vec,void *);
PETSC_EXTERN PetscErrorCode MatFDColoringSetF(MatFDColoring,Vec);
PETSC_EXTERN PetscErrorCode MatFDColoringGetPerturbedColumns(MatFDColoring,PetscInt*,const PetscInt*[]);
PETSC_EXTERN PetscErrorCode MatFDColoringSetUp(Mat,ISColoring,MatFDColoring);
PETSC_EXTERN PetscErrorCode MatFDColoringSetBlockSize(MatFDColoring,PetscInt,PetscInt);


/*S
     MatTransposeColoring - Object for computing a sparse matrix product C=A*B^T via coloring

   Level: beginner

  Concepts: coloring, sparse matrix product

.seealso:  MatTransposeColoringCreate()
S*/
typedef struct _p_MatTransposeColoring* MatTransposeColoring;

PETSC_EXTERN PetscErrorCode MatTransposeColoringCreate(Mat,ISColoring,MatTransposeColoring *);
PETSC_EXTERN PetscErrorCode MatTransColoringApplySpToDen(MatTransposeColoring,Mat,Mat);
PETSC_EXTERN PetscErrorCode MatTransColoringApplyDenToSp(MatTransposeColoring,Mat,Mat);
PETSC_EXTERN PetscErrorCode MatTransposeColoringDestroy(MatTransposeColoring*);

/*
    These routines are for partitioning matrices: currently used only
  for adjacency matrix, MatCreateMPIAdj().
*/

/*S
     MatPartitioning - Object for managing the partitioning of a matrix or graph

   Level: beginner

  Concepts: partitioning

.seealso:  MatPartitioningCreate(), MatPartitioningType
S*/
typedef struct _p_MatPartitioning* MatPartitioning;

/*J
    MatPartitioningType - String with the name of a PETSc matrix partitioning

   Level: beginner
dm
.seealso: MatPartitioningCreate(), MatPartitioning
J*/
typedef const char* MatPartitioningType;
#define MATPARTITIONINGCURRENT  "current"
#define MATPARTITIONINGAVERAGE   "average"
#define MATPARTITIONINGSQUARE   "square"
#define MATPARTITIONINGPARMETIS "parmetis"
#define MATPARTITIONINGCHACO    "chaco"
#define MATPARTITIONINGPARTY    "party"
#define MATPARTITIONINGPTSCOTCH "ptscotch"
#define MATPARTITIONINGHIERARCH  "hierarch"


PETSC_EXTERN PetscErrorCode MatPartitioningCreate(MPI_Comm,MatPartitioning*);
PETSC_EXTERN PetscErrorCode MatPartitioningSetType(MatPartitioning,MatPartitioningType);
PETSC_EXTERN PetscErrorCode MatPartitioningSetNParts(MatPartitioning,PetscInt);
PETSC_EXTERN PetscErrorCode MatPartitioningSetAdjacency(MatPartitioning,Mat);
PETSC_EXTERN PetscErrorCode MatPartitioningSetVertexWeights(MatPartitioning,const PetscInt[]);
PETSC_EXTERN PetscErrorCode MatPartitioningSetPartitionWeights(MatPartitioning,const PetscReal []);
PETSC_EXTERN PetscErrorCode MatPartitioningApply(MatPartitioning,IS*);
PETSC_EXTERN PetscErrorCode MatPartitioningApplyND(MatPartitioning,IS*);
PETSC_EXTERN PetscErrorCode MatPartitioningDestroy(MatPartitioning*);

PETSC_EXTERN PetscErrorCode MatPartitioningRegister(const char[],PetscErrorCode (*)(MatPartitioning));



PETSC_EXTERN PetscErrorCode MatPartitioningView(MatPartitioning,PetscViewer);
PETSC_STATIC_INLINE PetscErrorCode MatPartitioningViewFromOptions(MatPartitioning A,PetscObject B,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,B,name);}

PETSC_EXTERN PetscErrorCode MatPartitioningSetFromOptions(MatPartitioning);
PETSC_EXTERN PetscErrorCode MatPartitioningGetType(MatPartitioning,MatPartitioningType*);

PETSC_EXTERN PetscErrorCode MatPartitioningParmetisSetRepartition(MatPartitioning part);
PETSC_EXTERN PetscErrorCode MatPartitioningParmetisSetCoarseSequential(MatPartitioning);
PETSC_EXTERN PetscErrorCode MatPartitioningParmetisGetEdgeCut(MatPartitioning, PetscInt *);

typedef enum { MP_CHACO_MULTILEVEL=1,MP_CHACO_SPECTRAL=2,MP_CHACO_LINEAR=4,MP_CHACO_RANDOM=5,MP_CHACO_SCATTERED=6 } MPChacoGlobalType;
PETSC_EXTERN const char *const MPChacoGlobalTypes[];
typedef enum { MP_CHACO_KERNIGHAN=1,MP_CHACO_NONE=2 } MPChacoLocalType;
PETSC_EXTERN const char *const MPChacoLocalTypes[];
typedef enum { MP_CHACO_LANCZOS=0,MP_CHACO_RQI=1 } MPChacoEigenType;
PETSC_EXTERN const char *const MPChacoEigenTypes[];

PETSC_EXTERN PetscErrorCode MatPartitioningChacoSetGlobal(MatPartitioning,MPChacoGlobalType);
PETSC_EXTERN PetscErrorCode MatPartitioningChacoGetGlobal(MatPartitioning,MPChacoGlobalType*);
PETSC_EXTERN PetscErrorCode MatPartitioningChacoSetLocal(MatPartitioning,MPChacoLocalType);
PETSC_EXTERN PetscErrorCode MatPartitioningChacoGetLocal(MatPartitioning,MPChacoLocalType*);
PETSC_EXTERN PetscErrorCode MatPartitioningChacoSetCoarseLevel(MatPartitioning,PetscReal);
PETSC_EXTERN PetscErrorCode MatPartitioningChacoSetEigenSolver(MatPartitioning,MPChacoEigenType);
PETSC_EXTERN PetscErrorCode MatPartitioningChacoGetEigenSolver(MatPartitioning,MPChacoEigenType*);
PETSC_EXTERN PetscErrorCode MatPartitioningChacoSetEigenTol(MatPartitioning,PetscReal);
PETSC_EXTERN PetscErrorCode MatPartitioningChacoGetEigenTol(MatPartitioning,PetscReal*);
PETSC_EXTERN PetscErrorCode MatPartitioningChacoSetEigenNumber(MatPartitioning,PetscInt);
PETSC_EXTERN PetscErrorCode MatPartitioningChacoGetEigenNumber(MatPartitioning,PetscInt*);

#define MP_PARTY_OPT "opt"
#define MP_PARTY_LIN "lin"
#define MP_PARTY_SCA "sca"
#define MP_PARTY_RAN "ran"
#define MP_PARTY_GBF "gbf"
#define MP_PARTY_GCF "gcf"
#define MP_PARTY_BUB "bub"
#define MP_PARTY_DEF "def"
PETSC_EXTERN PetscErrorCode MatPartitioningPartySetGlobal(MatPartitioning,const char*);
#define MP_PARTY_HELPFUL_SETS "hs"
#define MP_PARTY_KERNIGHAN_LIN "kl"
#define MP_PARTY_NONE "no"
PETSC_EXTERN PetscErrorCode MatPartitioningPartySetLocal(MatPartitioning,const char*);
PETSC_EXTERN PetscErrorCode MatPartitioningPartySetCoarseLevel(MatPartitioning,PetscReal);
PETSC_EXTERN PetscErrorCode MatPartitioningPartySetBipart(MatPartitioning,PetscBool);
PETSC_EXTERN PetscErrorCode MatPartitioningPartySetMatchOptimization(MatPartitioning,PetscBool);

typedef enum { MP_PTSCOTCH_DEFAULT,MP_PTSCOTCH_QUALITY,MP_PTSCOTCH_SPEED,MP_PTSCOTCH_BALANCE,MP_PTSCOTCH_SAFETY,MP_PTSCOTCH_SCALABILITY } MPPTScotchStrategyType;
PETSC_EXTERN const char *const MPPTScotchStrategyTypes[];

PETSC_EXTERN PetscErrorCode MatPartitioningPTScotchSetImbalance(MatPartitioning,PetscReal);
PETSC_EXTERN PetscErrorCode MatPartitioningPTScotchGetImbalance(MatPartitioning,PetscReal*);
PETSC_EXTERN PetscErrorCode MatPartitioningPTScotchSetStrategy(MatPartitioning,MPPTScotchStrategyType);
PETSC_EXTERN PetscErrorCode MatPartitioningPTScotchGetStrategy(MatPartitioning,MPPTScotchStrategyType*);

/*
 * hierarchical partitioning
 */
PETSC_EXTERN PetscErrorCode MatPartitioningHierarchicalGetFineparts(MatPartitioning,IS*);
PETSC_EXTERN PetscErrorCode MatPartitioningHierarchicalGetCoarseparts(MatPartitioning,IS*);
PETSC_EXTERN PetscErrorCode MatPartitioningHierarchicalSetNcoarseparts(MatPartitioning,PetscInt);
PETSC_EXTERN PetscErrorCode MatPartitioningHierarchicalSetNfineparts(MatPartitioning, PetscInt);

PETSC_EXTERN PetscErrorCode MatMeshToVertexGraph(Mat,PetscInt,Mat*);
PETSC_EXTERN PetscErrorCode MatMeshToCellGraph(Mat,PetscInt,Mat*);

/*
    If you add entries here you must also add them to petsc/finclude/petscmat.h
*/
typedef enum { MATOP_SET_VALUES=0,
               MATOP_GET_ROW=1,
               MATOP_RESTORE_ROW=2,
               MATOP_MULT=3,
               MATOP_MULT_ADD=4,
               MATOP_MULT_TRANSPOSE=5,
               MATOP_MULT_TRANSPOSE_ADD=6,
               MATOP_SOLVE=7,
               MATOP_SOLVE_ADD=8,
               MATOP_SOLVE_TRANSPOSE=9,
               MATOP_SOLVE_TRANSPOSE_ADD=10,
               MATOP_LUFACTOR=11,
               MATOP_CHOLESKYFACTOR=12,
               MATOP_SOR=13,
               MATOP_TRANSPOSE=14,
               MATOP_GETINFO=15,
               MATOP_EQUAL=16,
               MATOP_GET_DIAGONAL=17,
               MATOP_DIAGONAL_SCALE=18,
               MATOP_NORM=19,
               MATOP_ASSEMBLY_BEGIN=20,
               MATOP_ASSEMBLY_END=21,
               MATOP_SET_OPTION=22,
               MATOP_ZERO_ENTRIES=23,
               MATOP_ZERO_ROWS=24,
               MATOP_LUFACTOR_SYMBOLIC=25,
               MATOP_LUFACTOR_NUMERIC=26,
               MATOP_CHOLESKY_FACTOR_SYMBOLIC=27,
               MATOP_CHOLESKY_FACTOR_NUMERIC=28,
               MATOP_SETUP_PREALLOCATION=29,
               MATOP_ILUFACTOR_SYMBOLIC=30,
               MATOP_ICCFACTOR_SYMBOLIC=31,
               MATOP_GET_DIAGONAL_BLOCK=32,
               /* MATOP_PLACEHOLDER_33=33, */
               MATOP_DUPLICATE=34,
               MATOP_FORWARD_SOLVE=35,
               MATOP_BACKWARD_SOLVE=36,
               MATOP_ILUFACTOR=37,
               MATOP_ICCFACTOR=38,
               MATOP_AXPY=39,
               MATOP_CREATE_SUBMATRICES=40,
               MATOP_INCREASE_OVERLAP=41,
               MATOP_GET_VALUES=42,
               MATOP_COPY=43,
               MATOP_GET_ROW_MAX=44,
               MATOP_SCALE=45,
               MATOP_SHIFT=46,
               MATOP_DIAGONAL_SET=47,
               MATOP_ZERO_ROWS_COLUMNS=48,
               MATOP_SET_RANDOM=49,
               MATOP_GET_ROW_IJ=50,
               MATOP_RESTORE_ROW_IJ=51,
               MATOP_GET_COLUMN_IJ=52,
               MATOP_RESTORE_COLUMN_IJ=53,
               MATOP_FDCOLORING_CREATE=54,
               MATOP_COLORING_PATCH=55,
               MATOP_SET_UNFACTORED=56,
               MATOP_PERMUTE=57,
               MATOP_SET_VALUES_BLOCKED=58,
               MATOP_CREATE_SUBMATRIX=59,
               MATOP_DESTROY=60,
               MATOP_VIEW=61,
               MATOP_CONVERT_FROM=62,
               MATOP_MATMAT_MULT=63,
               MATOP_MATMAT_MULT_SYMBOLIC=64,
               MATOP_MATMAT_MULT_NUMERIC=65,
               MATOP_SET_LOCAL_TO_GLOBAL_MAP=66,
               MATOP_SET_VALUES_LOCAL=67,
               MATOP_ZERO_ROWS_LOCAL=68,
               MATOP_GET_ROW_MAX_ABS=69,
               MATOP_GET_ROW_MIN_ABS=70,
               MATOP_CONVERT=71,
               MATOP_SET_COLORING=72,
               /* MATOP_PLACEHOLDER_73=73, */
               MATOP_SET_VALUES_ADIFOR=74,
               MATOP_FD_COLORING_APPLY=75,
               MATOP_SET_FROM_OPTIONS=76,
               MATOP_MULT_CONSTRAINED=77,
               MATOP_MULT_TRANSPOSE_CONSTRAIN=78,
               MATOP_FIND_ZERO_DIAGONALS=79,
               MATOP_MULT_MULTIPLE=80,
               MATOP_SOLVE_MULTIPLE=81,
               MATOP_GET_INERTIA=82,
               MATOP_LOAD=83,
               MATOP_IS_SYMMETRIC=84,
               MATOP_IS_HERMITIAN=85,
               MATOP_IS_STRUCTURALLY_SYMMETRIC=86,
               MATOP_SET_VALUES_BLOCKEDLOCAL=87,
               MATOP_CREATE_VECS=88,
               MATOP_MAT_MULT=89,
               MATOP_MAT_MULT_SYMBOLIC=90,
               MATOP_MAT_MULT_NUMERIC=91,
               MATOP_PTAP=92,
               MATOP_PTAP_SYMBOLIC=93,
               MATOP_PTAP_NUMERIC=94,
               MATOP_MAT_TRANSPOSE_MULT=95,
               MATOP_MAT_TRANSPOSE_MULT_SYMBO=96,
               MATOP_MAT_TRANSPOSE_MULT_NUMER=97,
               /* MATOP_PLACEHOLDER_98=98, */
               /* MATOP_PLACEHOLDER_99=99, */
               /* MATOP_PLACEHOLDER_100=100, */
               /* MATOP_PLACEHOLDER_101=101, */
               MATOP_CONJUGATE=102,
               /* MATOP_PLACEHOLDER_103=103, */
               MATOP_SET_VALUES_ROW=104,
               MATOP_REAL_PART=105,
               MATOP_IMAGINARY_PART=106,
               MATOP_GET_ROW_UPPER_TRIANGULAR=107,
               MATOP_RESTORE_ROW_UPPER_TRIANG=108,
               MATOP_MAT_SOLVE=109,
               MATOP_MAT_SOLVE_TRANSPOSE=110,
               MATOP_GET_ROW_MIN=111,
               MATOP_GET_COLUMN_VECTOR=112,
               MATOP_MISSING_DIAGONAL=113,
               MATOP_GET_SEQ_NONZERO_STRUCTUR=114,
               MATOP_CREATE=115,
               MATOP_GET_GHOSTS=116,
               MATOP_GET_LOCAL_SUB_MATRIX=117,
               MATOP_RESTORE_LOCALSUB_MATRIX=118,
               MATOP_MULT_DIAGONAL_BLOCK=119,
               MATOP_HERMITIAN_TRANSPOSE=120,
               MATOP_MULT_HERMITIAN_TRANSPOSE=121,
               MATOP_MULT_HERMITIAN_TRANS_ADD=122,
               MATOP_GET_MULTI_PROC_BLOCK=123,
               MATOP_FIND_NONZERO_ROWS=124,
               MATOP_GET_COLUMN_NORMS=125,
               MATOP_INVERT_BLOCK_DIAGONAL=126,
               /* MATOP_PLACEHOLDER_127=127, */
               MATOP_CREATE_SUB_MATRICES_MPI=128,
               MATOP_SET_VALUES_BATCH=129,
               MATOP_TRANSPOSE_MAT_MULT=130,
               MATOP_TRANSPOSE_MAT_MULT_SYMBO=131,
               MATOP_TRANSPOSE_MAT_MULT_NUMER=132,
               MATOP_TRANSPOSE_COLORING_CREAT=133,
               MATOP_TRANS_COLORING_APPLY_SPT=134,
               MATOP_TRANS_COLORING_APPLY_DEN=135,
               MATOP_RART=136,
               MATOP_RART_SYMBOLIC=137,
               MATOP_RART_NUMERIC=138,
               MATOP_SET_BLOCK_SIZES=139,
               MATOP_AYPX=140,
               MATOP_RESIDUAL=141,
               MATOP_FDCOLORING_SETUP=142,
               MATOP_MPICONCATENATESEQ=144,
               MATOP_DESTROYSUBMATRICES=145
             } MatOperation;
PETSC_EXTERN PetscErrorCode MatSetOperation(Mat,MatOperation,void(*)(void));
PETSC_EXTERN PetscErrorCode MatGetOperation(Mat,MatOperation,void(**)(void));
PETSC_EXTERN PetscErrorCode MatHasOperation(Mat,MatOperation,PetscBool *);
PETSC_EXTERN PetscErrorCode MatHasCongruentLayouts(Mat,PetscBool *);

PETSC_EXTERN PetscErrorCode MatShellSetOperation(Mat,MatOperation,void(*)(void));
PETSC_EXTERN PetscErrorCode MatShellGetOperation(Mat,MatOperation,void(**)(void));
PETSC_EXTERN PetscErrorCode MatShellSetContext(Mat,void*);
PETSC_EXTERN PetscErrorCode MatShellTestMult(Mat,PetscErrorCode (*)(void*,Vec,Vec),Vec,void*,PetscBool*);
PETSC_EXTERN PetscErrorCode MatShellTestMultTranspose(Mat,PetscErrorCode (*)(void*,Vec,Vec),Vec,void*,PetscBool*);
PETSC_EXTERN PetscErrorCode MatShellSetManageScalingShifts(Mat);

/*
   Codes for matrices stored on disk. By default they are
   stored in a universal format. By changing the format with
   PetscViewerPushFormat(viewer,PETSC_VIEWER_NATIVE); the matrices will
   be stored in a way natural for the matrix, for example dense matrices
   would be stored as dense. Matrices stored this way may only be
   read into matrices of the same type.
*/
#define MATRIX_BINARY_FORMAT_DENSE -1

PETSC_EXTERN PetscErrorCode MatMPIBAIJSetHashTableFactor(Mat,PetscReal);

PETSC_EXTERN PetscErrorCode MatISSetLocalMatType(Mat,MatType);
PETSC_EXTERN PetscErrorCode MatISSetPreallocation(Mat,PetscInt,const PetscInt[],PetscInt,const PetscInt[]);
PETSC_EXTERN PetscErrorCode MatISSetUpSF(Mat);
PETSC_EXTERN PetscErrorCode MatISStoreL2L(Mat,PetscBool);
PETSC_EXTERN PetscErrorCode MatISFixLocalEmpty(Mat,PetscBool);
PETSC_EXTERN PetscErrorCode MatISGetLocalMat(Mat,Mat*);
PETSC_EXTERN PetscErrorCode MatISRestoreLocalMat(Mat,Mat*);
PETSC_EXTERN PetscErrorCode MatISSetLocalMat(Mat,Mat);
PETSC_EXTERN PETSC_DEPRECATED("Use the MatConvert() interface") PetscErrorCode MatISGetMPIXAIJ(Mat,MatReuse,Mat*);

/*S
     MatNullSpace - Object that removes a null space from a vector, i.e.
         orthogonalizes the vector to a subsapce

   Level: advanced

  Concepts: matrix; linear operator, null space

  Users manual sections:
.   sec_singular

.seealso:  MatNullSpaceCreate()
S*/
typedef struct _p_MatNullSpace* MatNullSpace;

PETSC_EXTERN PetscErrorCode MatNullSpaceCreate(MPI_Comm,PetscBool ,PetscInt,const Vec[],MatNullSpace*);
PETSC_EXTERN PetscErrorCode MatNullSpaceSetFunction(MatNullSpace,PetscErrorCode (*)(MatNullSpace,Vec,void*),void*);
PETSC_EXTERN PetscErrorCode MatNullSpaceDestroy(MatNullSpace*);
PETSC_EXTERN PetscErrorCode MatNullSpaceRemove(MatNullSpace,Vec);
PETSC_EXTERN PetscErrorCode MatGetNullSpace(Mat, MatNullSpace *);
PETSC_EXTERN PetscErrorCode MatGetTransposeNullSpace(Mat, MatNullSpace *);
PETSC_EXTERN PetscErrorCode MatSetTransposeNullSpace(Mat,MatNullSpace);
PETSC_EXTERN PetscErrorCode MatSetNullSpace(Mat,MatNullSpace);
PETSC_EXTERN PetscErrorCode MatSetNearNullSpace(Mat,MatNullSpace);
PETSC_EXTERN PetscErrorCode MatGetNearNullSpace(Mat,MatNullSpace*);
PETSC_EXTERN PetscErrorCode MatNullSpaceTest(MatNullSpace,Mat,PetscBool  *);
PETSC_EXTERN PetscErrorCode MatNullSpaceView(MatNullSpace,PetscViewer);
PETSC_EXTERN PetscErrorCode MatNullSpaceGetVecs(MatNullSpace,PetscBool*,PetscInt*,const Vec**);
PETSC_EXTERN PetscErrorCode MatNullSpaceCreateRigidBody(Vec,MatNullSpace*);

PETSC_EXTERN PetscErrorCode MatReorderingSeqSBAIJ(Mat,IS);
PETSC_EXTERN PetscErrorCode MatMPISBAIJSetHashTableFactor(Mat,PetscReal);
PETSC_EXTERN PetscErrorCode MatSeqSBAIJSetColumnIndices(Mat,PetscInt *);
PETSC_EXTERN PetscErrorCode MatSeqBAIJInvertBlockDiagonal(Mat);

PETSC_EXTERN PetscErrorCode MatCreateMAIJ(Mat,PetscInt,Mat*);
PETSC_EXTERN PetscErrorCode MatMAIJRedimension(Mat,PetscInt,Mat*);
PETSC_EXTERN PetscErrorCode MatMAIJGetAIJ(Mat,Mat*);

PETSC_EXTERN PetscErrorCode MatComputeExplicitOperator(Mat,Mat*);
PETSC_EXTERN PetscErrorCode MatComputeExplicitOperatorTranspose(Mat,Mat*);

PETSC_EXTERN PetscErrorCode MatDiagonalScaleLocal(Mat,Vec);

PETSC_EXTERN PetscErrorCode MatCreateMFFD(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,Mat*);
PETSC_EXTERN PetscErrorCode MatMFFDSetBase(Mat,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatMFFDSetFunction(Mat,PetscErrorCode(*)(void*,Vec,Vec),void*);
PETSC_EXTERN PetscErrorCode MatMFFDSetFunctioni(Mat,PetscErrorCode (*)(void*,PetscInt,Vec,PetscScalar*));
PETSC_EXTERN PetscErrorCode MatMFFDSetFunctioniBase(Mat,PetscErrorCode (*)(void*,Vec));
PETSC_EXTERN PetscErrorCode MatMFFDSetHHistory(Mat,PetscScalar[],PetscInt);
PETSC_EXTERN PetscErrorCode MatMFFDResetHHistory(Mat);
PETSC_EXTERN PetscErrorCode MatMFFDSetFunctionError(Mat,PetscReal);
PETSC_EXTERN PetscErrorCode MatMFFDSetPeriod(Mat,PetscInt);
PETSC_EXTERN PetscErrorCode MatMFFDGetH(Mat,PetscScalar *);
PETSC_EXTERN PetscErrorCode MatMFFDSetOptionsPrefix(Mat,const char[]);
PETSC_EXTERN PetscErrorCode MatMFFDCheckPositivity(void*,Vec,Vec,PetscScalar*);
PETSC_EXTERN PetscErrorCode MatMFFDSetCheckh(Mat,PetscErrorCode (*)(void*,Vec,Vec,PetscScalar*),void*);

/*S
    MatMFFD - A data structured used to manage the computation of the h differencing parameter for matrix-free
              Jacobian vector products

    Notes:
    MATMFFD is a specific MatType which uses the MatMFFD data structure

           MatMFFD*() methods actually take the Mat as their first argument. Not a MatMFFD data structure

    Level: developer

.seealso: MATMFFD, MatCreateMFFD(), MatMFFDSetFuction(), MatMFFDSetType(), MatMFFDRegister()
S*/
typedef struct _p_MatMFFD* MatMFFD;

/*J
    MatMFFDType - algorithm used to compute the h used in computing matrix-vector products via differencing of the function

   Level: beginner

.seealso: MatMFFDSetType(), MatMFFDRegister()
J*/
typedef const char* MatMFFDType;
#define MATMFFD_DS  "ds"
#define MATMFFD_WP  "wp"

PETSC_EXTERN PetscErrorCode MatMFFDSetType(Mat,MatMFFDType);
PETSC_EXTERN PetscErrorCode MatMFFDRegister(const char[],PetscErrorCode (*)(MatMFFD));

PETSC_EXTERN PetscErrorCode MatMFFDDSSetUmin(Mat,PetscReal);
PETSC_EXTERN PetscErrorCode MatMFFDWPSetComputeNormU(Mat,PetscBool );

PETSC_EXTERN PetscErrorCode MatFDColoringSetType(MatFDColoring,MatMFFDType);

PETSC_EXTERN PetscErrorCode PetscViewerMathematicaPutMatrix(PetscViewer, PetscInt, PetscInt, PetscReal *);
PETSC_EXTERN PetscErrorCode PetscViewerMathematicaPutCSRMatrix(PetscViewer, PetscInt, PetscInt, PetscInt *, PetscInt *, PetscReal *);

/*
   PETSc interface to MUMPS
*/
#ifdef PETSC_HAVE_MUMPS
PETSC_EXTERN PetscErrorCode MatMumpsSetIcntl(Mat,PetscInt,PetscInt);
PETSC_EXTERN PetscErrorCode MatMumpsGetIcntl(Mat,PetscInt,PetscInt*);
PETSC_EXTERN PetscErrorCode MatMumpsSetCntl(Mat,PetscInt,PetscReal);
PETSC_EXTERN PetscErrorCode MatMumpsGetCntl(Mat,PetscInt,PetscReal*);

PETSC_EXTERN PetscErrorCode MatMumpsGetInfo(Mat,PetscInt,PetscInt*);
PETSC_EXTERN PetscErrorCode MatMumpsGetInfog(Mat,PetscInt,PetscInt*);
PETSC_EXTERN PetscErrorCode MatMumpsGetRinfo(Mat,PetscInt,PetscReal*);
PETSC_EXTERN PetscErrorCode MatMumpsGetRinfog(Mat,PetscInt,PetscReal*);
PETSC_EXTERN PetscErrorCode MatMumpsGetInverse(Mat,Mat);
PETSC_EXTERN PetscErrorCode MatMumpsGetInverseTranspose(Mat,Mat);
#endif

/*
   PETSc interface to Mkl_Pardiso
*/
#ifdef PETSC_HAVE_MKL_PARDISO
PETSC_EXTERN PetscErrorCode MatMkl_PardisoSetCntl(Mat,PetscInt,PetscInt);
#endif

/*
   PETSc interface to Mkl_CPardiso
*/
#ifdef PETSC_HAVE_MKL_CPARDISO
PETSC_EXTERN PetscErrorCode MatMkl_CPardisoSetCntl(Mat,PetscInt,PetscInt);
#endif

/*
   PETSc interface to SUPERLU
*/
#ifdef PETSC_HAVE_SUPERLU
PETSC_EXTERN PetscErrorCode MatSuperluSetILUDropTol(Mat,PetscReal);
#endif

/*
   PETSc interface to SUPERLU_DIST
*/
#ifdef PETSC_HAVE_SUPERLU_DIST
PETSC_EXTERN PetscErrorCode MatSuperluDistGetDiagU(Mat,PetscScalar*);
#endif

/*
   PETSc interface to STRUMPACK
*/
#ifdef PETSC_HAVE_STRUMPACK
/*E
    MatSTRUMPACKReordering - sparsity reducing ordering to be used in STRUMPACK

    Level: intermediate
E*/
typedef enum {MAT_STRUMPACK_NATURAL=0,
              MAT_STRUMPACK_METIS=1,
              MAT_STRUMPACK_PARMETIS=2,
              MAT_STRUMPACK_SCOTCH=3,
              MAT_STRUMPACK_PTSCOTCH=4,
              MAT_STRUMPACK_RCM=5} MatSTRUMPACKReordering;
PETSC_EXTERN PetscErrorCode MatSTRUMPACKSetReordering(Mat,MatSTRUMPACKReordering);
PETSC_EXTERN PetscErrorCode MatSTRUMPACKSetColPerm(Mat,PetscBool);
PETSC_EXTERN PetscErrorCode MatSTRUMPACKSetHSSRelTol(Mat,PetscReal);
PETSC_DEPRECATED ("Use MatSTRUMPACKSetHSSRelTol()") PETSC_STATIC_INLINE PetscErrorCode MatSTRUMPACKSetHSSRelCompTol(Mat mat,PetscReal rtol) {return MatSTRUMPACKSetHSSRelTol(mat,rtol);}
PETSC_EXTERN PetscErrorCode MatSTRUMPACKSetHSSAbsTol(Mat,PetscReal);
PETSC_EXTERN PetscErrorCode MatSTRUMPACKSetHSSMinSepSize(Mat,PetscInt);
PETSC_DEPRECATED ("Use MatSTRUMPACKSetHSSMinSepSize()") PETSC_STATIC_INLINE PetscErrorCode MatSTRUMPACKSetHSSMinSize(Mat mat,PetscInt hssminsize) {return MatSTRUMPACKSetHSSMinSepSize(mat,hssminsize);}
PETSC_EXTERN PetscErrorCode MatSTRUMPACKSetHSSMaxRank(Mat,PetscInt);
PETSC_EXTERN PetscErrorCode MatSTRUMPACKSetHSSLeafSize(Mat,PetscInt);
#endif


#ifdef PETSC_HAVE_CUDA
/*E
    MatCUSPARSEStorageFormat - indicates the storage format for CUSPARSE (GPU)
    matrices.

    Not Collective

+   MAT_CUSPARSE_CSR - Compressed Sparse Row
.   MAT_CUSPARSE_ELL - Ellpack (requires CUDA 4.2 or later).
-   MAT_CUSPARSE_HYB - Hybrid, a combination of Ellpack and Coordinate format (requires CUDA 4.2 or later).

    Level: intermediate

   Any additions/changes here MUST also be made in include/petsc/finclude/petscmat.h

.seealso: MatCUSPARSESetFormat(), MatCUSPARSEFormatOperation
E*/

typedef enum {MAT_CUSPARSE_CSR, MAT_CUSPARSE_ELL, MAT_CUSPARSE_HYB} MatCUSPARSEStorageFormat;

/* these will be strings associated with enumerated type defined above */
PETSC_EXTERN const char *const MatCUSPARSEStorageFormats[];

/*E
    MatCUSPARSEFormatOperation - indicates the operation of CUSPARSE (GPU)
    matrices whose operation should use a particular storage format.

    Not Collective

+   MAT_CUSPARSE_MULT_DIAG - sets the storage format for the diagonal matrix in the parallel MatMult
.   MAT_CUSPARSE_MULT_OFFDIAG - sets the storage format for the offdiagonal matrix in the parallel MatMult
.   MAT_CUSPARSE_MULT - sets the storage format for the entire matrix in the serial (single GPU) MatMult
-   MAT_CUSPARSE_ALL - sets the storage format for all CUSPARSE (GPU) matrices

    Level: intermediate

.seealso: MatCUSPARSESetFormat(), MatCUSPARSEStorageFormat
E*/
typedef enum {MAT_CUSPARSE_MULT_DIAG, MAT_CUSPARSE_MULT_OFFDIAG, MAT_CUSPARSE_MULT, MAT_CUSPARSE_ALL} MatCUSPARSEFormatOperation;

PETSC_EXTERN PetscErrorCode MatCreateSeqAIJCUSPARSE(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);
PETSC_EXTERN PetscErrorCode MatCreateAIJCUSPARSE(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);
PETSC_EXTERN PetscErrorCode MatCUSPARSESetFormat(Mat,MatCUSPARSEFormatOperation,MatCUSPARSEStorageFormat);
#endif

#if defined(PETSC_HAVE_VIENNACL)
PETSC_EXTERN PetscErrorCode MatCreateSeqAIJViennaCL(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],Mat*);
PETSC_EXTERN PetscErrorCode MatCreateAIJViennaCL(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Mat*);
#endif

/*
   PETSc interface to FFTW
*/
#if defined(PETSC_HAVE_FFTW)
PETSC_EXTERN PetscErrorCode VecScatterPetscToFFTW(Mat,Vec,Vec);
PETSC_EXTERN PetscErrorCode VecScatterFFTWToPetsc(Mat,Vec,Vec);
PETSC_EXTERN PetscErrorCode MatCreateVecsFFTW(Mat,Vec*,Vec*,Vec*);
#endif

PETSC_EXTERN PetscErrorCode MatCreateNest(MPI_Comm,PetscInt,const IS[],PetscInt,const IS[],const Mat[],Mat*);
PETSC_EXTERN PetscErrorCode MatNestGetSize(Mat,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode MatNestGetISs(Mat,IS[],IS[]);
PETSC_EXTERN PetscErrorCode MatNestGetLocalISs(Mat,IS[],IS[]);
PETSC_EXTERN PetscErrorCode MatNestGetSubMats(Mat,PetscInt*,PetscInt*,Mat***);
PETSC_EXTERN PetscErrorCode MatNestGetSubMat(Mat,PetscInt,PetscInt,Mat*);
PETSC_EXTERN PetscErrorCode MatNestSetVecType(Mat,VecType);
PETSC_EXTERN PetscErrorCode MatNestSetSubMats(Mat,PetscInt,const IS[],PetscInt,const IS[],const Mat[]);
PETSC_EXTERN PetscErrorCode MatNestSetSubMat(Mat,PetscInt,PetscInt,Mat);

PETSC_EXTERN PetscErrorCode MatChop(Mat,PetscReal);
PETSC_EXTERN PetscErrorCode MatComputeBandwidth(Mat,PetscReal,PetscInt*);

PETSC_EXTERN PetscErrorCode MatSubdomainsCreateCoalesce(Mat,PetscInt,PetscInt*,IS**);

PETSC_EXTERN PetscErrorCode MatPreallocatorPreallocate(Mat,PetscBool,Mat);

PETSC_INTERN PetscErrorCode MatHeaderMerge(Mat,Mat*);
PETSC_EXTERN PetscErrorCode MatHeaderReplace(Mat,Mat*);

#endif
