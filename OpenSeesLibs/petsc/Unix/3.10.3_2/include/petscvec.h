/*
    Defines the vector component of PETSc. Vectors generally represent
  degrees of freedom for finite element/finite difference functions
  on a grid. They have more mathematical structure then simple arrays.
*/

#ifndef __PETSCVEC_H
#define __PETSCVEC_H
#include <petscis.h>
#include <petscviewer.h>

/*S
     Vec - Abstract PETSc vector object

   Level: beginner

  Concepts: field variables, unknowns, arrays

.seealso:  VecCreate(), VecType, VecSetType()
S*/
typedef struct _p_Vec*         Vec;

/*S
     VecScatter - Object used to manage communication of data
       between vectors in parallel. Manages both scatters and gathers

   Level: beginner

  Concepts: scatter

.seealso:  VecScatterCreate(), VecScatterBegin(), VecScatterEnd()
S*/
typedef struct _p_VecScatter*  VecScatter;

/*E
  ScatterMode - Determines the direction of a scatter

  Level: beginner

.seealso: VecScatter, VecScatterBegin(), VecScatterEnd()
E*/
typedef enum {SCATTER_FORWARD=0, SCATTER_REVERSE=1, SCATTER_FORWARD_LOCAL=2, SCATTER_REVERSE_LOCAL=3, SCATTER_LOCAL=2} ScatterMode;

/*MC
    SCATTER_FORWARD - Scatters the values as dictated by the VecScatterCreate() call

    Level: beginner

.seealso: VecScatter, ScatterMode, VecScatterCreate(), VecScatterBegin(), VecScatterEnd(), SCATTER_REVERSE, SCATTER_FORWARD_LOCAL,
          SCATTER_REVERSE_LOCAL

M*/

/*MC
    SCATTER_REVERSE - Moves the values in the opposite direction then the directions indicated in
         in the VecScatterCreate()

    Level: beginner

.seealso: VecScatter, ScatterMode, VecScatterCreate(), VecScatterBegin(), VecScatterEnd(), SCATTER_FORWARD, SCATTER_FORWARD_LOCAL,
          SCATTER_REVERSE_LOCAL

M*/

/*MC
    SCATTER_FORWARD_LOCAL - Scatters the values as dictated by the VecScatterCreate() call except NO parallel communication
       is done. Any variables that have be moved between processes are ignored

    Level: developer

.seealso: VecScatter, ScatterMode, VecScatterCreate(), VecScatterBegin(), VecScatterEnd(), SCATTER_REVERSE, SCATTER_FORWARD,
          SCATTER_REVERSE_LOCAL

M*/

/*MC
    SCATTER_REVERSE_LOCAL - Moves the values in the opposite direction then the directions indicated in
         in the VecScatterCreate()  except NO parallel communication
       is done. Any variables that have be moved between processes are ignored

    Level: developer

.seealso: VecScatter, ScatterMode, VecScatterCreate(), VecScatterBegin(), VecScatterEnd(), SCATTER_FORWARD, SCATTER_FORWARD_LOCAL,
          SCATTER_REVERSE

M*/

/*J
    VecType - String with the name of a PETSc vector

   Level: beginner

.seealso: VecSetType(), Vec, VecCreate(), VecDestroy()
J*/
typedef const char* VecType;
#define VECSEQ         "seq"
#define VECMPI         "mpi"
#define VECSTANDARD    "standard"   /* seq on one process and mpi on several */
#define VECSHARED      "shared"
#define VECSEQVIENNACL "seqviennacl"
#define VECMPIVIENNACL "mpiviennacl"
#define VECVIENNACL    "viennacl"   /* seqviennacl on one process and mpiviennacl on several */
#define VECSEQCUDA     "seqcuda"
#define VECMPICUDA     "mpicuda"
#define VECCUDA        "cuda"       /* seqcuda on one process and mpicuda on several */
#define VECNEST        "nest"
#define VECNODE        "node"       /* use on-node shared memory */

/*J
    VecScatterType - String with the name of a PETSc vector scatter type

   Level: beginner

.seealso: VecScatterSetType(), VecScatter, VecScatterCreate(), VecScatterDestroy()
J*/
typedef const char* VecScatterType;
#define VECSCATTERSEQ       "seq"
#define VECSCATTERMPI1      "mpi1"
#define VECSCATTERMPI3      "mpi3"     /* use MPI3 on-node shared memory */
#define VECSCATTERMPI3NODE  "mpi3node" /* use MPI3 on-node shared memory for vector type VECNODE */

/* Dynamic creation and loading functions */
PETSC_EXTERN PetscFunctionList VecScatterList;
PETSC_EXTERN PetscErrorCode VecScatterSetType(VecScatter, VecScatterType);
PETSC_EXTERN PetscErrorCode VecScatterGetType(VecScatter, VecScatterType *);
PETSC_EXTERN PetscErrorCode VecScatterSetFromOptions(VecScatter);
PETSC_EXTERN PetscErrorCode VecScatterRegister(const char[],PetscErrorCode (*)(VecScatter));
PETSC_EXTERN PetscErrorCode VecScatterCreate(Vec,IS,Vec,IS,VecScatter*);
PETSC_EXTERN PetscErrorCode VecScatterInitializePackage(void);
PETSC_EXTERN PetscErrorCode VecScatterFinalizePackage(void);

/* Logging support */
#define    REAL_FILE_CLASSID 1211213
#define    VEC_FILE_CLASSID 1211214
PETSC_EXTERN PetscClassId VEC_CLASSID;
PETSC_EXTERN PetscClassId VEC_SCATTER_CLASSID;


PETSC_EXTERN PetscErrorCode VecInitializePackage(void);
PETSC_EXTERN PetscErrorCode VecFinalizePackage(void);

PETSC_EXTERN PetscErrorCode VecCreate(MPI_Comm,Vec*);
PETSC_EXTERN PetscErrorCode VecCreateSeq(MPI_Comm,PetscInt,Vec*);
PETSC_EXTERN PetscErrorCode VecCreateMPI(MPI_Comm,PetscInt,PetscInt,Vec*);
PETSC_EXTERN PetscErrorCode VecCreateSeqWithArray(MPI_Comm,PetscInt,PetscInt,const PetscScalar[],Vec*);
PETSC_EXTERN PetscErrorCode VecCreateMPIWithArray(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscScalar[],Vec*);
PETSC_EXTERN PetscErrorCode VecCreateShared(MPI_Comm,PetscInt,PetscInt,Vec*);
PETSC_EXTERN PetscErrorCode VecCreateNode(MPI_Comm,PetscInt,PetscInt,Vec*);

PETSC_EXTERN PetscErrorCode VecSetFromOptions(Vec);
PETSC_STATIC_INLINE PetscErrorCode VecViewFromOptions(Vec A,PetscObject B,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,B,name);}

PETSC_EXTERN PetscErrorCode VecSetUp(Vec);
PETSC_EXTERN PetscErrorCode VecDestroy(Vec*);
PETSC_EXTERN PetscErrorCode VecZeroEntries(Vec);
PETSC_EXTERN PetscErrorCode VecSetOptionsPrefix(Vec,const char[]);
PETSC_EXTERN PetscErrorCode VecAppendOptionsPrefix(Vec,const char[]);
PETSC_EXTERN PetscErrorCode VecGetOptionsPrefix(Vec,const char*[]);

PETSC_EXTERN PetscErrorCode VecSetSizes(Vec,PetscInt,PetscInt);

PETSC_EXTERN PetscErrorCode VecDotNorm2(Vec,Vec,PetscScalar*,PetscReal*);
PETSC_EXTERN PetscErrorCode VecDot(Vec,Vec,PetscScalar*);
PETSC_EXTERN PetscErrorCode VecDotRealPart(Vec,Vec,PetscReal*);
PETSC_EXTERN PetscErrorCode VecTDot(Vec,Vec,PetscScalar*);
PETSC_EXTERN PetscErrorCode VecMDot(Vec,PetscInt,const Vec[],PetscScalar[]);
PETSC_EXTERN PetscErrorCode VecMTDot(Vec,PetscInt,const Vec[],PetscScalar[]);
PETSC_EXTERN PetscErrorCode VecGetSubVector(Vec,IS,Vec*);
PETSC_EXTERN PetscErrorCode VecRestoreSubVector(Vec,IS,Vec*);

/*E
    NormType - determines what type of norm to compute

    Level: beginner

.seealso: VecNorm(), VecNormBegin(), VecNormEnd(), MatNorm()
E*/
typedef enum {NORM_1=0,NORM_2=1,NORM_FROBENIUS=2,NORM_INFINITY=3,NORM_1_AND_2=4} NormType;
PETSC_EXTERN const char *const NormTypes[];
#define NORM_MAX NORM_INFINITY

/*MC
     NORM_1 - the one norm, ||v|| = sum_i | v_i |. ||A|| = max_j || v_*j ||, maximum column sum

   Level: beginner

.seealso:  NormType, MatNorm(), VecNorm(), VecNormBegin(), VecNormEnd(), NORM_2, NORM_FROBENIUS,
           NORM_INFINITY, NORM_1_AND_2

M*/

/*MC
     NORM_2 - the two norm, ||v|| = sqrt(sum_i |v_i|^2) (vectors only)

   Level: beginner

.seealso:  NormType, MatNorm(), VecNorm(), VecNormBegin(), VecNormEnd(), NORM_1, NORM_FROBENIUS,
           NORM_INFINITY, NORM_1_AND_2

M*/

/*MC
     NORM_FROBENIUS - ||A|| = sqrt(sum_ij |A_ij|^2), same as NORM_2 for vectors

   Level: beginner

.seealso:  NormType, MatNorm(), VecNorm(), VecNormBegin(), VecNormEnd(), NORM_1, NORM_2,
           NORM_INFINITY, NORM_1_AND_2

M*/

/*MC
     NORM_INFINITY - ||v|| = max_i |v_i|. ||A|| = max_i || v_i* ||, maximum row sum

   Level: beginner

.seealso:  NormType, MatNorm(), VecNorm(), VecNormBegin(), VecNormEnd(), NORM_1, NORM_2,
           NORM_FROBINIUS, NORM_1_AND_2

M*/

/*MC
     NORM_1_AND_2 - computes both the 1 and 2 norm of a vector

   Level: beginner

.seealso:  NormType, MatNorm(), VecNorm(), VecNormBegin(), VecNormEnd(), NORM_1, NORM_2,
           NORM_FROBINIUS, NORM_INFINITY

M*/

/*MC
     NORM_MAX - see NORM_INFINITY

   Level: beginner

M*/

PETSC_EXTERN PetscErrorCode VecNorm(Vec,NormType,PetscReal *);
PETSC_EXTERN PetscErrorCode VecNormAvailable(Vec,NormType,PetscBool *,PetscReal *);
PETSC_EXTERN PetscErrorCode VecNormalize(Vec,PetscReal *);
PETSC_EXTERN PetscErrorCode VecSum(Vec,PetscScalar*);
PETSC_EXTERN PetscErrorCode VecMax(Vec,PetscInt*,PetscReal *);
PETSC_EXTERN PetscErrorCode VecMin(Vec,PetscInt*,PetscReal *);
PETSC_EXTERN PetscErrorCode VecScale(Vec,PetscScalar);
PETSC_EXTERN PetscErrorCode VecCopy(Vec,Vec);
PETSC_EXTERN PetscErrorCode VecSetRandom(Vec,PetscRandom);
PETSC_EXTERN PetscErrorCode VecSet(Vec,PetscScalar);
PETSC_EXTERN PetscErrorCode VecSetInf(Vec);
PETSC_EXTERN PetscErrorCode VecSwap(Vec,Vec);
PETSC_EXTERN PetscErrorCode VecAXPY(Vec,PetscScalar,Vec);
PETSC_EXTERN PetscErrorCode VecAXPBY(Vec,PetscScalar,PetscScalar,Vec);
PETSC_EXTERN PetscErrorCode VecMAXPY(Vec,PetscInt,const PetscScalar[],Vec[]);
PETSC_EXTERN PetscErrorCode VecAYPX(Vec,PetscScalar,Vec);
PETSC_EXTERN PetscErrorCode VecWAXPY(Vec,PetscScalar,Vec,Vec);
PETSC_EXTERN PetscErrorCode VecAXPBYPCZ(Vec,PetscScalar,PetscScalar,PetscScalar,Vec,Vec);
PETSC_EXTERN PetscErrorCode VecPointwiseMax(Vec,Vec,Vec);
PETSC_EXTERN PetscErrorCode VecPointwiseMaxAbs(Vec,Vec,Vec);
PETSC_EXTERN PetscErrorCode VecPointwiseMin(Vec,Vec,Vec);
PETSC_EXTERN PetscErrorCode VecPointwiseMult(Vec,Vec,Vec);
PETSC_EXTERN PetscErrorCode VecPointwiseDivide(Vec,Vec,Vec);
PETSC_EXTERN PetscErrorCode VecMaxPointwiseDivide(Vec,Vec,PetscReal*);
PETSC_EXTERN PetscErrorCode VecShift(Vec,PetscScalar);
PETSC_EXTERN PetscErrorCode VecReciprocal(Vec);
PETSC_EXTERN PetscErrorCode VecPermute(Vec, IS, PetscBool );
PETSC_EXTERN PetscErrorCode VecSqrtAbs(Vec);
PETSC_EXTERN PetscErrorCode VecLog(Vec);
PETSC_EXTERN PetscErrorCode VecExp(Vec);
PETSC_EXTERN PetscErrorCode VecAbs(Vec);
PETSC_EXTERN PetscErrorCode VecDuplicate(Vec,Vec*);
PETSC_EXTERN PetscErrorCode VecDuplicateVecs(Vec,PetscInt,Vec*[]);
PETSC_EXTERN PetscErrorCode VecDestroyVecs(PetscInt, Vec*[]);
PETSC_EXTERN PetscErrorCode VecStrideNormAll(Vec,NormType,PetscReal[]);
PETSC_EXTERN PetscErrorCode VecStrideMaxAll(Vec,PetscInt [],PetscReal []);
PETSC_EXTERN PetscErrorCode VecStrideMinAll(Vec,PetscInt [],PetscReal []);
PETSC_EXTERN PetscErrorCode VecStrideScaleAll(Vec,const PetscScalar[]);
PETSC_EXTERN PetscErrorCode VecUniqueEntries(Vec,PetscInt*,PetscScalar**);

PETSC_EXTERN PetscErrorCode VecStrideNorm(Vec,PetscInt,NormType,PetscReal*);
PETSC_EXTERN PetscErrorCode VecStrideMax(Vec,PetscInt,PetscInt *,PetscReal *);
PETSC_EXTERN PetscErrorCode VecStrideMin(Vec,PetscInt,PetscInt *,PetscReal *);
PETSC_EXTERN PetscErrorCode VecStrideScale(Vec,PetscInt,PetscScalar);
PETSC_EXTERN PetscErrorCode VecStrideSet(Vec,PetscInt,PetscScalar);


PETSC_EXTERN PetscErrorCode VecStrideGather(Vec,PetscInt,Vec,InsertMode);
PETSC_EXTERN PetscErrorCode VecStrideScatter(Vec,PetscInt,Vec,InsertMode);
PETSC_EXTERN PetscErrorCode VecStrideGatherAll(Vec,Vec[],InsertMode);
PETSC_EXTERN PetscErrorCode VecStrideScatterAll(Vec[],Vec,InsertMode);

PETSC_EXTERN PetscErrorCode VecStrideSubSetScatter(Vec,PetscInt,const PetscInt[],const PetscInt[],Vec,InsertMode);
PETSC_EXTERN PetscErrorCode VecStrideSubSetGather(Vec,PetscInt,const PetscInt[],const PetscInt[],Vec,InsertMode);

PETSC_EXTERN PetscErrorCode VecSetValues(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
PETSC_EXTERN PetscErrorCode VecGetValues(Vec,PetscInt,const PetscInt[],PetscScalar[]);
PETSC_EXTERN PetscErrorCode VecAssemblyBegin(Vec);
PETSC_EXTERN PetscErrorCode VecAssemblyEnd(Vec);
PETSC_EXTERN PetscErrorCode VecStashSetInitialSize(Vec,PetscInt,PetscInt);
PETSC_EXTERN PetscErrorCode VecStashView(Vec,PetscViewer);
PETSC_EXTERN PetscErrorCode VecStashViewFromOptions(Vec,PetscObject,const char[]);
PETSC_EXTERN PetscErrorCode VecStashGetInfo(Vec,PetscInt*,PetscInt*,PetscInt*,PetscInt*);

/*MC
   VecSetValue - Set a single entry into a vector.

   Synopsis:
   #include <petscvec.h>
   PetscErrorCode VecSetValue(Vec v,PetscInt row,PetscScalar value, InsertMode mode);

   Not Collective

   Input Parameters:
+  v - the vector
.  row - the row location of the entry
.  value - the value to insert
-  mode - either INSERT_VALUES or ADD_VALUES

   Notes:
   For efficiency one should use VecSetValues() and set several or
   many values simultaneously if possible.

   These values may be cached, so VecAssemblyBegin() and VecAssemblyEnd()
   MUST be called after all calls to VecSetValue() have been completed.

   VecSetValue() uses 0-based indices in Fortran as well as in C.

   Level: beginner

.seealso: VecSetValues(), VecAssemblyBegin(), VecAssemblyEnd(), VecSetValuesBlockedLocal(), VecSetValueLocal()
M*/
PETSC_STATIC_INLINE PetscErrorCode VecSetValue(Vec v,PetscInt i,PetscScalar va,InsertMode mode) {return VecSetValues(v,1,&i,&va,mode);}


PETSC_EXTERN PetscErrorCode VecSetBlockSize(Vec,PetscInt);
PETSC_EXTERN PetscErrorCode VecGetBlockSize(Vec,PetscInt*);
PETSC_EXTERN PetscErrorCode VecSetValuesBlocked(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);

/* Dynamic creation and loading functions */
PETSC_EXTERN PetscFunctionList VecList;
PETSC_EXTERN PetscErrorCode VecSetType(Vec, VecType);
PETSC_EXTERN PetscErrorCode VecGetType(Vec, VecType *);
PETSC_EXTERN PetscErrorCode VecRegister(const char[],PetscErrorCode (*)(Vec));

PETSC_EXTERN PetscErrorCode VecScatterCreate(Vec,IS,Vec,IS,VecScatter *);
PETSC_EXTERN PetscErrorCode VecScatterCreateEmpty(MPI_Comm,VecScatter *);
PETSC_EXTERN PetscErrorCode VecScatterBegin(VecScatter,Vec,Vec,InsertMode,ScatterMode);
PETSC_EXTERN PetscErrorCode VecScatterEnd(VecScatter,Vec,Vec,InsertMode,ScatterMode);
PETSC_EXTERN PetscErrorCode VecScatterDestroy(VecScatter*);
PETSC_EXTERN PetscErrorCode VecScatterCopy(VecScatter,VecScatter *);
PETSC_EXTERN PetscErrorCode VecScatterView(VecScatter,PetscViewer);
PETSC_STATIC_INLINE PetscErrorCode VecScatterViewFromOptions(VecScatter A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}
PETSC_EXTERN PetscErrorCode VecScatterRemap(VecScatter,PetscInt[],PetscInt[]);
PETSC_EXTERN PetscErrorCode VecScatterGetMerged(VecScatter,PetscBool *);

PETSC_EXTERN PetscErrorCode VecGetArray4d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar****[]);
PETSC_EXTERN PetscErrorCode VecRestoreArray4d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar****[]);
PETSC_EXTERN PetscErrorCode VecGetArray3d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar***[]);
PETSC_EXTERN PetscErrorCode VecRestoreArray3d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar***[]);
PETSC_EXTERN PetscErrorCode VecGetArray2d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar**[]);
PETSC_EXTERN PetscErrorCode VecRestoreArray2d(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar**[]);
PETSC_EXTERN PetscErrorCode VecGetArray1d(Vec,PetscInt,PetscInt,PetscScalar *[]);
PETSC_EXTERN PetscErrorCode VecRestoreArray1d(Vec,PetscInt,PetscInt,PetscScalar *[]);

PETSC_EXTERN PetscErrorCode VecGetArray4dRead(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar****[]);
PETSC_EXTERN PetscErrorCode VecRestoreArray4dRead(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar****[]);
PETSC_EXTERN PetscErrorCode VecGetArray3dRead(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar***[]);
PETSC_EXTERN PetscErrorCode VecRestoreArray3dRead(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar***[]);
PETSC_EXTERN PetscErrorCode VecGetArray2dRead(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar**[]);
PETSC_EXTERN PetscErrorCode VecRestoreArray2dRead(Vec,PetscInt,PetscInt,PetscInt,PetscInt,PetscScalar**[]);
PETSC_EXTERN PetscErrorCode VecGetArray1dRead(Vec,PetscInt,PetscInt,PetscScalar *[]);
PETSC_EXTERN PetscErrorCode VecRestoreArray1dRead(Vec,PetscInt,PetscInt,PetscScalar *[]);

PETSC_EXTERN PetscErrorCode VecPlaceArray(Vec,const PetscScalar[]);
PETSC_EXTERN PetscErrorCode VecResetArray(Vec);
PETSC_EXTERN PetscErrorCode VecReplaceArray(Vec,const PetscScalar[]);
PETSC_EXTERN PetscErrorCode VecGetArrays(const Vec[],PetscInt,PetscScalar**[]);
PETSC_EXTERN PetscErrorCode VecRestoreArrays(const Vec[],PetscInt,PetscScalar**[]);

PETSC_EXTERN PetscErrorCode VecView(Vec,PetscViewer);
PETSC_EXTERN PetscErrorCode VecEqual(Vec,Vec,PetscBool *);
PETSC_EXTERN PetscErrorCode VecLoad(Vec, PetscViewer);

PETSC_EXTERN PetscErrorCode VecGetSize(Vec,PetscInt*);
PETSC_EXTERN PetscErrorCode VecGetLocalSize(Vec,PetscInt*);
PETSC_EXTERN PetscErrorCode VecGetOwnershipRange(Vec,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode VecGetOwnershipRanges(Vec,const PetscInt *[]);

PETSC_EXTERN PetscErrorCode VecSetLocalToGlobalMapping(Vec,ISLocalToGlobalMapping);
PETSC_EXTERN PetscErrorCode VecSetValuesLocal(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);

/*MC
   VecSetValueLocal - Set a single entry into a vector using the local numbering

   Synopsis:
   #include <petscvec.h>
   PetscErrorCode VecSetValueLocal(Vec v,PetscInt row,PetscScalar value, InsertMode mode);

   Not Collective

   Input Parameters:
+  v - the vector
.  row - the row location of the entry
.  value - the value to insert
-  mode - either INSERT_VALUES or ADD_VALUES

   Notes:
   For efficiency one should use VecSetValues() and set several or
   many values simultaneously if possible.

   These values may be cached, so VecAssemblyBegin() and VecAssemblyEnd()
   MUST be called after all calls to VecSetValues() have been completed.

   VecSetValues() uses 0-based indices in Fortran as well as in C.

   Level: beginner

.seealso: VecSetValues(), VecAssemblyBegin(), VecAssemblyEnd(), VecSetValuesBlockedLocal(), VecSetValue()
M*/
PETSC_STATIC_INLINE PetscErrorCode VecSetValueLocal(Vec v,PetscInt i,PetscScalar va,InsertMode mode) {return VecSetValuesLocal(v,1,&i,&va,mode);}

PETSC_EXTERN PetscErrorCode VecSetValuesBlockedLocal(Vec,PetscInt,const PetscInt[],const PetscScalar[],InsertMode);
PETSC_EXTERN PetscErrorCode VecGetLocalToGlobalMapping(Vec,ISLocalToGlobalMapping*);

PETSC_EXTERN PetscErrorCode VecDotBegin(Vec,Vec,PetscScalar *);
PETSC_EXTERN PetscErrorCode VecDotEnd(Vec,Vec,PetscScalar *);
PETSC_EXTERN PetscErrorCode VecTDotBegin(Vec,Vec,PetscScalar *);
PETSC_EXTERN PetscErrorCode VecTDotEnd(Vec,Vec,PetscScalar *);
PETSC_EXTERN PetscErrorCode VecNormBegin(Vec,NormType,PetscReal *);
PETSC_EXTERN PetscErrorCode VecNormEnd(Vec,NormType,PetscReal *);

PETSC_EXTERN PetscErrorCode VecMDotBegin(Vec,PetscInt,const Vec[],PetscScalar[]);
PETSC_EXTERN PetscErrorCode VecMDotEnd(Vec,PetscInt,const Vec[],PetscScalar[]);
PETSC_EXTERN PetscErrorCode VecMTDotBegin(Vec,PetscInt,const Vec[],PetscScalar[]);
PETSC_EXTERN PetscErrorCode VecMTDotEnd(Vec,PetscInt,const Vec[],PetscScalar[]);
PETSC_EXTERN PetscErrorCode PetscCommSplitReductionBegin(MPI_Comm);


typedef enum {VEC_IGNORE_OFF_PROC_ENTRIES,VEC_IGNORE_NEGATIVE_INDICES,VEC_SUBSET_OFF_PROC_ENTRIES} VecOption;
PETSC_EXTERN PetscErrorCode VecSetOption(Vec,VecOption,PetscBool );

PETSC_EXTERN PetscErrorCode VecGetArray(Vec,PetscScalar**);
PETSC_EXTERN PetscErrorCode VecGetArrayRead(Vec,const PetscScalar**);
PETSC_EXTERN PetscErrorCode VecRestoreArray(Vec,PetscScalar**);
PETSC_EXTERN PetscErrorCode VecRestoreArrayRead(Vec,const PetscScalar**);
PETSC_EXTERN PetscErrorCode VecGetLocalVector(Vec,Vec);
PETSC_EXTERN PetscErrorCode VecRestoreLocalVector(Vec,Vec);
PETSC_EXTERN PetscErrorCode VecGetLocalVectorRead(Vec,Vec);
PETSC_EXTERN PetscErrorCode VecRestoreLocalVectorRead(Vec,Vec);

/*@C
   VecGetArrayPair - Accesses a pair of pointers for two vectors that may be common. When not common the first is read only

   Logically Collective on Vec

   Input Parameters:
+  x - the vector
-  y - the second vector

   Output Parameters:
+  xv - location to put pointer to the first array
-  yv - location to put pointer to the second array

   Level: developer

   Not available from Fortran

.seealso: VecGetArray(), VecGetArrayRead(), VecRestoreArrayPair()

@*/
PETSC_STATIC_INLINE PetscErrorCode VecGetArrayPair(Vec x,Vec y,PetscScalar **xv,PetscScalar **yv)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecGetArray(y,yv);CHKERRQ(ierr);
  if (x != y) {
    ierr = VecGetArrayRead(x,(const PetscScalar **)xv);CHKERRQ(ierr);
  } else {
    *xv = *yv;
  }
  PetscFunctionReturn(0);
}

/*@C
   VecRestoreArrayPair - Returns a pair of pointers for two vectors that may be common. When not common the first is read only

   Logically Collective on Vec

   Input Parameters:
+  x - the vector
-  y - the second vector

   Output Parameters:
+  xv - location to put pointer to the first array
-  yv - location to put pointer to the second array

   Level: developer

   Not available from Fortran

.seealso: VecGetArray(), VecGetArrayRead(), VecGetArrayPair()

@*/
PETSC_STATIC_INLINE PetscErrorCode VecRestoreArrayPair(Vec x,Vec y,PetscScalar **xv,PetscScalar **yv)
{
  PetscErrorCode ierr;

  PetscFunctionBegin;
  ierr = VecRestoreArray(y,yv);CHKERRQ(ierr);
  if (x != y) {
    ierr = VecRestoreArrayRead(x,(const PetscScalar **)xv);CHKERRQ(ierr);
  }
  PetscFunctionReturn(0);
}

#if defined(PETSC_USE_DEBUG)
PETSC_EXTERN PetscErrorCode VecLockGet(Vec,PetscInt*);
PETSC_EXTERN PetscErrorCode VecLockPush(Vec);
PETSC_EXTERN PetscErrorCode VecLockPop(Vec);
#define VecLocked(x,arg) do {PetscInt _st; PetscErrorCode __ierr = VecLockGet(x,&_st); CHKERRQ(__ierr); if (_st > 0) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE," Vec is locked read only, argument # %d",arg);} while (0)
#else
#define VecLockGet(x,arg)     *(arg) = 0
#define VecLockPush(x)        0
#define VecLockPop(x)         0
#define VecLocked(x,arg)
#endif

PETSC_EXTERN PetscErrorCode VecValidValues(Vec,PetscInt,PetscBool);

/*
    These numbers need to match the entries in
  the function table in vecimpl.h
*/
typedef enum { VECOP_DUPLICATE = 0, VECOP_VIEW = 33, VECOP_LOAD = 41, VECOP_VIEWNATIVE = 68, VECOP_LOADNATIVE = 69 } VecOperation;
PETSC_EXTERN PetscErrorCode VecSetOperation(Vec,VecOperation,void(*)(void));

/*
     Routines for dealing with ghosted vectors:
  vectors with ghost elements at the end of the array.
*/
PETSC_EXTERN PetscErrorCode VecMPISetGhost(Vec,PetscInt,const PetscInt[]);
PETSC_EXTERN PetscErrorCode VecCreateGhost(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],Vec*);
PETSC_EXTERN PetscErrorCode VecCreateGhostWithArray(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscScalar[],Vec*);
PETSC_EXTERN PetscErrorCode VecCreateGhostBlock(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],Vec*);
PETSC_EXTERN PetscErrorCode VecCreateGhostBlockWithArray(MPI_Comm,PetscInt,PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscScalar[],Vec*);
PETSC_EXTERN PetscErrorCode VecGhostGetLocalForm(Vec,Vec*);
PETSC_EXTERN PetscErrorCode VecGhostRestoreLocalForm(Vec,Vec*);
PETSC_EXTERN PetscErrorCode VecGhostIsLocalForm(Vec,Vec,PetscBool*);
PETSC_EXTERN PetscErrorCode VecGhostUpdateBegin(Vec,InsertMode,ScatterMode);
PETSC_EXTERN PetscErrorCode VecGhostUpdateEnd(Vec,InsertMode,ScatterMode);

PETSC_EXTERN PetscErrorCode VecConjugate(Vec);

PETSC_EXTERN PetscErrorCode VecScatterCreateToAll(Vec,VecScatter*,Vec*);
PETSC_EXTERN PetscErrorCode VecScatterCreateToZero(Vec,VecScatter*,Vec*);

PETSC_EXTERN PetscErrorCode ISComplementVec(IS,Vec,IS*);
PETSC_EXTERN PetscErrorCode VecPow(Vec, PetscScalar);
PETSC_EXTERN PetscErrorCode VecMedian(Vec, Vec, Vec, Vec);
PETSC_EXTERN PetscErrorCode VecWhichInactive(Vec, Vec, Vec, Vec, PetscBool, IS *);
PETSC_EXTERN PetscErrorCode VecWhichBetween(Vec, Vec, Vec, IS *);
PETSC_EXTERN PetscErrorCode VecWhichBetweenOrEqual(Vec, Vec, Vec, IS *);
PETSC_EXTERN PetscErrorCode VecWhichGreaterThan(Vec, Vec, IS * );
PETSC_EXTERN PetscErrorCode VecWhichLessThan(Vec, Vec, IS *);
PETSC_EXTERN PetscErrorCode VecWhichEqual(Vec, Vec, IS *);
PETSC_EXTERN PetscErrorCode VecISAXPY(Vec, IS, PetscScalar,Vec);
PETSC_EXTERN PetscErrorCode VecISCopy(Vec, IS, ScatterMode, Vec);
PETSC_EXTERN PetscErrorCode VecISSet(Vec,IS, PetscScalar);
PETSC_EXTERN PetscErrorCode VecBoundGradientProjection(Vec, Vec, Vec, Vec, Vec);
PETSC_EXTERN PetscErrorCode VecStepBoundInfo(Vec,Vec,Vec,Vec,PetscReal*, PetscReal*,PetscReal*);
PETSC_EXTERN PetscErrorCode VecStepMax(Vec, Vec, PetscReal *);
PETSC_EXTERN PetscErrorCode VecStepMaxBounded(Vec,Vec,Vec,Vec,PetscReal*);

PETSC_EXTERN PetscErrorCode PetscViewerMathematicaGetVector(PetscViewer, Vec);
PETSC_EXTERN PetscErrorCode PetscViewerMathematicaPutVector(PetscViewer, Vec);

/*S
     Vecs - Collection of vectors where the data for the vectors is stored in
            one contiguous memory

   Level: advanced

   Notes:
    Temporary construct for handling multiply right hand side solves

    This is faked by storing a single vector that has enough array space for
    n vectors

  Concepts: parallel decomposition

S*/
        struct _n_Vecs  {PetscInt n; Vec v;};
typedef struct _n_Vecs* Vecs;
PETSC_EXTERN PetscErrorCode VecsDestroy(Vecs);
PETSC_EXTERN PetscErrorCode VecsCreateSeq(MPI_Comm,PetscInt,PetscInt,Vecs*);
PETSC_EXTERN PetscErrorCode VecsCreateSeqWithArray(MPI_Comm,PetscInt,PetscInt,PetscScalar*,Vecs*);
PETSC_EXTERN PetscErrorCode VecsDuplicate(Vecs,Vecs*);

#if defined(PETSC_HAVE_VIENNACL)
typedef struct _p_PetscViennaCLIndices* PetscViennaCLIndices;
PETSC_EXTERN PetscErrorCode PetscViennaCLIndicesCreate(PetscInt, PetscInt*,PetscInt, PetscInt*,PetscViennaCLIndices*);
PETSC_EXTERN PetscErrorCode PetscViennaCLIndicesDestroy(PetscViennaCLIndices*);
PETSC_EXTERN PetscErrorCode VecViennaCLCopyToGPUSome_Public(Vec,PetscViennaCLIndices);
PETSC_EXTERN PetscErrorCode VecViennaCLCopyFromGPUSome_Public(Vec,PetscViennaCLIndices);
PETSC_EXTERN PetscErrorCode VecCreateSeqViennaCL(MPI_Comm,PetscInt,Vec*);
PETSC_EXTERN PetscErrorCode VecCreateMPIViennaCL(MPI_Comm,PetscInt,PetscInt,Vec*);
#endif
#if defined(PETSC_HAVE_VECCUDA)
typedef struct _p_PetscCUDAIndices* PetscCUDAIndices;
typedef struct _p_VecScatterCUDAIndices_StoS* VecScatterCUDAIndices_StoS;
typedef struct _p_VecScatterCUDAIndices_PtoP* VecScatterCUDAIndices_PtoP;
PETSC_EXTERN PetscErrorCode VecCUDACopyToGPUSome_Public(Vec,PetscCUDAIndices);
PETSC_EXTERN PetscErrorCode VecCUDACopyFromGPUSome_Public(Vec,PetscCUDAIndices);
PETSC_EXTERN PetscErrorCode VecScatterInitializeForGPU(VecScatter,Vec,ScatterMode);
PETSC_EXTERN PetscErrorCode VecScatterFinalizeForGPU(VecScatter);
PETSC_EXTERN PetscErrorCode VecCreateSeqCUDA(MPI_Comm,PetscInt,Vec*);
PETSC_EXTERN PetscErrorCode VecCreateSeqCUDAWithArray(MPI_Comm,PetscInt,PetscInt,const PetscScalar*,Vec*);
PETSC_EXTERN PetscErrorCode VecCreateMPICUDA(MPI_Comm,PetscInt,PetscInt,Vec*);
PETSC_EXTERN PetscErrorCode VecCreateMPICUDAWithArray(MPI_Comm,PetscInt,PetscInt,PetscInt,const PetscScalar*,Vec*);
#endif

PETSC_EXTERN PetscErrorCode VecNestGetSubVecs(Vec,PetscInt*,Vec**);
PETSC_EXTERN PetscErrorCode VecNestGetSubVec(Vec,PetscInt,Vec*);
PETSC_EXTERN PetscErrorCode VecNestSetSubVecs(Vec,PetscInt,PetscInt*,Vec*);
PETSC_EXTERN PetscErrorCode VecNestSetSubVec(Vec,PetscInt,Vec);
PETSC_EXTERN PetscErrorCode VecCreateNest(MPI_Comm,PetscInt,IS*,Vec*,Vec*);
PETSC_EXTERN PetscErrorCode VecNestGetSize(Vec,PetscInt*);

PETSC_EXTERN PetscErrorCode PetscOptionsGetVec(PetscOptions,const char[],const char[],Vec,PetscBool*);
PETSC_EXTERN PetscErrorCode VecChop(Vec,PetscReal);

PETSC_EXTERN PetscErrorCode VecGetLayout(Vec,PetscLayout*);
PETSC_EXTERN PetscErrorCode VecSetLayout(Vec,PetscLayout);

PETSC_EXTERN PetscErrorCode PetscSectionVecView(PetscSection, Vec, PetscViewer);
PETSC_EXTERN PetscErrorCode VecGetValuesSection(Vec, PetscSection, PetscInt, PetscScalar **);
PETSC_EXTERN PetscErrorCode VecSetValuesSection(Vec, PetscSection, PetscInt, PetscScalar [], InsertMode);
PETSC_EXTERN PetscErrorCode PetscSectionVecNorm(PetscSection, PetscSection, Vec, NormType, PetscReal []);

PETSC_EXTERN PetscErrorCode PetscSFCreateFromZero(MPI_Comm,Vec,PetscSF*);

/*S
  VecTagger - Object used to manage the tagging of a subset of indices based on the values of a vector.  The
              motivating application is the selection of cells for refinement or coarsening based on vector containing
              the values in an error indicator metric.

  Level: advanced
S*/
typedef struct _p_VecTagger *VecTagger;

/*J
  VecTaggerType - String with the name of a VecTagger type

  Level: advanced
J*/
typedef const char* VecTaggerType;
/* tag where the vector values are in a box of explicitly defined values */
#define VECTAGGERABSOLUTE   "absolute"
/* tag where the vector values are in a box of values relative to the set of all values in the vector */
#define VECTAGGERRELATIVE   "relative"
/* tag where the vector values are in a relative range of the *cumulative distribution* of values in the vector */
#define VECTAGGERCDF        "cdf"
/* tag a vector as the union of other tags */
#define VECTAGGEROR         "or"
/* tag a vector as the intersection of other tags */
#define VECTAGGERAND        "and"

PETSC_EXTERN PetscClassId VEC_TAGGER_CLASSID;
PETSC_EXTERN PetscFunctionList VecTaggerList;
PETSC_EXTERN PetscErrorCode VecTaggerRegister(const char[],PetscErrorCode (*) (VecTagger));

PETSC_EXTERN PetscErrorCode VecTaggerCreate(MPI_Comm,VecTagger *);
PETSC_EXTERN PetscErrorCode VecTaggerSetBlockSize(VecTagger,PetscInt);
PETSC_EXTERN PetscErrorCode VecTaggerGetBlockSize(VecTagger,PetscInt*);
PETSC_EXTERN PetscErrorCode VecTaggerSetType(VecTagger,VecTaggerType);
PETSC_EXTERN PetscErrorCode VecTaggerGetType(VecTagger,VecTaggerType *);
PETSC_EXTERN PetscErrorCode VecTaggerSetInvert(VecTagger,PetscBool);
PETSC_EXTERN PetscErrorCode VecTaggerGetInvert(VecTagger,PetscBool*);
PETSC_EXTERN PetscErrorCode VecTaggerSetFromOptions(VecTagger);
PETSC_EXTERN PetscErrorCode VecTaggerSetUp(VecTagger);
PETSC_EXTERN PetscErrorCode VecTaggerView(VecTagger,PetscViewer);
PETSC_EXTERN PetscErrorCode VecTaggerComputeIS(VecTagger,Vec,IS *);
PETSC_EXTERN PetscErrorCode VecTaggerDestroy(VecTagger *);

/*S
   VecTaggerBox - A box range used to tag values.  For real scalars, this is just a closed interval; for complex scalars, the box is the closed region in the complex plane
   such that real(min) <= real(z) <= real(max) and imag(min) <= imag(z) <= imag(max).  INF is an acceptable endpoint.

   Level: beginner

.seealso: VecTaggerComputeIntervals()
S*/
typedef struct {
  PetscScalar min;
  PetscScalar max;
} VecTaggerBox;
PETSC_EXTERN PetscErrorCode VecTaggerComputeBoxes(VecTagger,Vec,PetscInt *,VecTaggerBox **);


PETSC_EXTERN PetscErrorCode VecTaggerAbsoluteSetBox(VecTagger,VecTaggerBox *);
PETSC_EXTERN PetscErrorCode VecTaggerAbsoluteGetBox(VecTagger,const VecTaggerBox **);

PETSC_EXTERN PetscErrorCode VecTaggerRelativeSetBox(VecTagger,VecTaggerBox *);
PETSC_EXTERN PetscErrorCode VecTaggerRelativeGetBox(VecTagger,const VecTaggerBox **);

PETSC_EXTERN PetscErrorCode VecTaggerCDFSetBox(VecTagger,VecTaggerBox *);
PETSC_EXTERN PetscErrorCode VecTaggerCDFGetBox(VecTagger,const VecTaggerBox **);

/*E
  VecTaggerCDFMethod - Determines what method is used to compute absolute values from cumulative distribution values (e.g., what value is the preimage of .95 in the cdf).  Relevant only in parallel: in serial it is directly computed.

  Level: advanced
.seealso: VecTaggerCDFSetMethod(), VecTaggerCDFMethods
E*/
typedef enum {VECTAGGER_CDF_GATHER,VECTAGGER_CDF_ITERATIVE,VECTAGGER_CDF_NUM_METHODS} VecTaggerCDFMethod;
PETSC_EXTERN const char *const VecTaggerCDFMethods[];

PETSC_EXTERN PetscErrorCode VecTaggerCDFSetMethod(VecTagger,VecTaggerCDFMethod);
PETSC_EXTERN PetscErrorCode VecTaggerCDFGetMethod(VecTagger,VecTaggerCDFMethod*);
PETSC_EXTERN PetscErrorCode VecTaggerCDFIterativeSetTolerances(VecTagger,PetscInt,PetscReal,PetscReal);
PETSC_EXTERN PetscErrorCode VecTaggerCDFIterativeGetTolerances(VecTagger,PetscInt*,PetscReal*,PetscReal*);

PETSC_EXTERN PetscErrorCode VecTaggerOrSetSubs(VecTagger,PetscInt,VecTagger*,PetscCopyMode);
PETSC_EXTERN PetscErrorCode VecTaggerOrGetSubs(VecTagger,PetscInt*,VecTagger**);

PETSC_EXTERN PetscErrorCode VecTaggerAndSetSubs(VecTagger,PetscInt,VecTagger*,PetscCopyMode);
PETSC_EXTERN PetscErrorCode VecTaggerAndGetSubs(VecTagger,PetscInt*,VecTagger**);

PETSC_EXTERN PetscErrorCode VecTaggerInitializePackage(void);
PETSC_EXTERN PetscErrorCode VecTaggerFinalizePackage(void);

#endif
