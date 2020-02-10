/*
   An index set is a generalization of a subset of integers.  Index sets
   are used for defining scatters and gathers.
*/
#if !defined(__PETSCIS_H)
#define __PETSCIS_H
#include <petscsys.h>
#include <petscsftypes.h>
#include <petscistypes.h>    /*I  "petscis.h" I*/

#define IS_FILE_CLASSID 1211218
PETSC_EXTERN PetscClassId IS_CLASSID;

PETSC_EXTERN PetscErrorCode ISInitializePackage(void);

/*J
    ISType - String with the name of a PETSc index set type

   Level: beginner

.seealso: ISSetType(), IS, ISCreate(), ISRegister()
J*/
typedef const char* ISType;
#define ISGENERAL      "general"
#define ISSTRIDE       "stride"
#define ISBLOCK        "block"

/* Dynamic creation and loading functions */
PETSC_EXTERN PetscFunctionList ISList;
PETSC_EXTERN PetscErrorCode ISSetType(IS, ISType);
PETSC_EXTERN PetscErrorCode ISGetType(IS, ISType *);
PETSC_EXTERN PetscErrorCode ISRegister(const char[],PetscErrorCode (*)(IS));
PETSC_EXTERN PetscErrorCode ISCreate(MPI_Comm,IS*);

/*
    Default index set data structures that PETSc provides.
*/
PETSC_EXTERN PetscErrorCode ISCreateGeneral(MPI_Comm,PetscInt,const PetscInt[],PetscCopyMode,IS *);
PETSC_EXTERN PetscErrorCode ISGeneralSetIndices(IS,PetscInt,const PetscInt[],PetscCopyMode);
PETSC_EXTERN PetscErrorCode ISCreateBlock(MPI_Comm,PetscInt,PetscInt,const PetscInt[],PetscCopyMode,IS *);
PETSC_EXTERN PetscErrorCode ISBlockSetIndices(IS,PetscInt,PetscInt,const PetscInt[],PetscCopyMode);
PETSC_EXTERN PetscErrorCode ISCreateStride(MPI_Comm,PetscInt,PetscInt,PetscInt,IS *);
PETSC_EXTERN PetscErrorCode ISStrideSetStride(IS,PetscInt,PetscInt,PetscInt);

PETSC_EXTERN PetscErrorCode ISDestroy(IS*);
PETSC_EXTERN PetscErrorCode ISSetPermutation(IS);
PETSC_EXTERN PetscErrorCode ISPermutation(IS,PetscBool *);
PETSC_EXTERN PetscErrorCode ISSetIdentity(IS);
PETSC_EXTERN PetscErrorCode ISIdentity(IS,PetscBool *);
PETSC_EXTERN PetscErrorCode ISContiguousLocal(IS,PetscInt,PetscInt,PetscInt*,PetscBool*);

PETSC_EXTERN PetscErrorCode ISGetIndices(IS,const PetscInt *[]);
PETSC_EXTERN PetscErrorCode ISRestoreIndices(IS,const PetscInt *[]);
PETSC_EXTERN PetscErrorCode ISGetTotalIndices(IS,const PetscInt *[]);
PETSC_EXTERN PetscErrorCode ISRestoreTotalIndices(IS,const PetscInt *[]);
PETSC_EXTERN PetscErrorCode ISGetNonlocalIndices(IS,const PetscInt *[]);
PETSC_EXTERN PetscErrorCode ISRestoreNonlocalIndices(IS,const PetscInt *[]);
PETSC_EXTERN PetscErrorCode ISGetNonlocalIS(IS, IS *is);
PETSC_EXTERN PetscErrorCode ISRestoreNonlocalIS(IS, IS *is);
PETSC_EXTERN PetscErrorCode ISGetSize(IS,PetscInt *);
PETSC_EXTERN PetscErrorCode ISGetLocalSize(IS,PetscInt *);
PETSC_EXTERN PetscErrorCode ISInvertPermutation(IS,PetscInt,IS*);
PETSC_EXTERN PetscErrorCode ISView(IS,PetscViewer);
PETSC_STATIC_INLINE PetscErrorCode ISViewFromOptions(IS A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}
PETSC_EXTERN PetscErrorCode ISLoad(IS,PetscViewer);
PETSC_EXTERN PetscErrorCode ISEqual(IS,IS,PetscBool  *);
PETSC_EXTERN PetscErrorCode ISEqualUnsorted(IS,IS,PetscBool *);
PETSC_EXTERN PetscErrorCode ISSort(IS);
PETSC_EXTERN PetscErrorCode ISSortRemoveDups(IS);
PETSC_EXTERN PetscErrorCode ISSorted(IS,PetscBool  *);
PETSC_EXTERN PetscErrorCode ISDifference(IS,IS,IS*);
PETSC_EXTERN PetscErrorCode ISSum(IS,IS,IS*);
PETSC_EXTERN PetscErrorCode ISExpand(IS,IS,IS*);
PETSC_EXTERN PetscErrorCode ISIntersect(IS,IS,IS*);
PETSC_EXTERN PetscErrorCode ISGetMinMax(IS,PetscInt*,PetscInt*);

PETSC_EXTERN PetscErrorCode ISLocate(IS,PetscInt,PetscInt*);
PETSC_EXTERN PetscErrorCode ISGetPointRange(IS,PetscInt*,PetscInt*,const PetscInt**);
PETSC_EXTERN PetscErrorCode ISRestorePointRange(IS,PetscInt*,PetscInt*,const PetscInt**);
PETSC_EXTERN PetscErrorCode ISGetPointSubrange(IS,PetscInt,PetscInt,const PetscInt*);

PETSC_EXTERN PetscErrorCode ISBlockGetIndices(IS,const PetscInt *[]);
PETSC_EXTERN PetscErrorCode ISBlockRestoreIndices(IS,const PetscInt *[]);
PETSC_EXTERN PetscErrorCode ISBlockGetLocalSize(IS,PetscInt *);
PETSC_EXTERN PetscErrorCode ISBlockGetSize(IS,PetscInt *);
PETSC_EXTERN PetscErrorCode ISGetBlockSize(IS,PetscInt*);
PETSC_EXTERN PetscErrorCode ISSetBlockSize(IS,PetscInt);

PETSC_EXTERN PetscErrorCode ISStrideGetInfo(IS,PetscInt *,PetscInt*);

PETSC_EXTERN PetscErrorCode ISToGeneral(IS);

PETSC_EXTERN PetscErrorCode ISDuplicate(IS,IS*);
PETSC_EXTERN PetscErrorCode ISCopy(IS,IS);
PETSC_EXTERN PetscErrorCode ISAllGather(IS,IS*);
PETSC_EXTERN PetscErrorCode ISComplement(IS,PetscInt,PetscInt,IS*);
PETSC_EXTERN PetscErrorCode ISConcatenate(MPI_Comm,PetscInt,const IS[],IS*);
PETSC_EXTERN PetscErrorCode ISListToPair(MPI_Comm,PetscInt, IS[],IS*,IS*);
PETSC_EXTERN PetscErrorCode ISPairToList(IS,IS,PetscInt*, IS *[]);
PETSC_EXTERN PetscErrorCode ISEmbed(IS,IS,PetscBool,IS*);
PETSC_EXTERN PetscErrorCode ISSortPermutation(IS,PetscBool,IS*);
PETSC_EXTERN PetscErrorCode ISOnComm(IS,MPI_Comm,PetscCopyMode,IS*);
PETSC_EXTERN PetscErrorCode ISRenumber(IS,IS,PetscInt*,IS*);

/* --------------------------------------------------------------------------*/
PETSC_EXTERN PetscClassId IS_LTOGM_CLASSID;

/*E
    ISGlobalToLocalMappingMode - Indicates if missing global indices are

   IS_GTOLM_MASK - missing global indices are replaced with -1
   IS_GTOLM_DROP - missing global indices are dropped

   Level: beginner

.seealso: ISGlobalToLocalMappingApplyBlock(), ISGlobalToLocalMappingApply()

E*/
typedef enum {IS_GTOLM_MASK,IS_GTOLM_DROP} ISGlobalToLocalMappingMode;

/*J
    ISLocalToGlobalMappingType - String with the name of a mapping method

   Level: beginner

.seealso: ISLocalToGlobalMappingSetType(), ISLocalToGlobalSetFromOptions()
J*/
typedef const char* ISLocalToGlobalMappingType;
#define ISLOCALTOGLOBALMAPPINGBASIC "basic"
#define ISLOCALTOGLOBALMAPPINGHASH  "hash"

PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingSetType(ISLocalToGlobalMapping,ISLocalToGlobalMappingType);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingRegister(const char[],PetscErrorCode (*)(ISLocalToGlobalMapping));
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingRegisterAll(void);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingCreate(MPI_Comm,PetscInt,PetscInt,const PetscInt[],PetscCopyMode,ISLocalToGlobalMapping*);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingCreateIS(IS,ISLocalToGlobalMapping *);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingCreateSF(PetscSF,PetscInt,ISLocalToGlobalMapping*);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingSetFromOptions(ISLocalToGlobalMapping);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingSetUp(ISLocalToGlobalMapping);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingView(ISLocalToGlobalMapping,PetscViewer);
PETSC_STATIC_INLINE PetscErrorCode ISLocalToGlobalMappingViewFromOptions(ISLocalToGlobalMapping A,PetscObject B,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,B,name);}

PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingDestroy(ISLocalToGlobalMapping*);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingApply(ISLocalToGlobalMapping,PetscInt,const PetscInt[],PetscInt[]);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingApplyBlock(ISLocalToGlobalMapping,PetscInt,const PetscInt[],PetscInt[]);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingApplyIS(ISLocalToGlobalMapping,IS,IS*);
PETSC_EXTERN PetscErrorCode ISGlobalToLocalMappingApply(ISLocalToGlobalMapping,ISGlobalToLocalMappingMode,PetscInt,const PetscInt[],PetscInt*,PetscInt[]);
PETSC_EXTERN PetscErrorCode ISGlobalToLocalMappingApplyBlock(ISLocalToGlobalMapping,ISGlobalToLocalMappingMode,PetscInt,const PetscInt[],PetscInt*,PetscInt[]);
PETSC_EXTERN PetscErrorCode ISGlobalToLocalMappingApplyIS(ISLocalToGlobalMapping,ISGlobalToLocalMappingMode,IS,IS*);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingGetSize(ISLocalToGlobalMapping,PetscInt*);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingGetNodeInfo(ISLocalToGlobalMapping,PetscInt*,PetscInt*[],PetscInt**[]);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingRestoreNodeInfo(ISLocalToGlobalMapping,PetscInt*,PetscInt*[],PetscInt**[]);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingGetInfo(ISLocalToGlobalMapping,PetscInt*,PetscInt*[],PetscInt*[],PetscInt**[]);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingRestoreInfo(ISLocalToGlobalMapping,PetscInt*,PetscInt*[],PetscInt*[],PetscInt**[]);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingGetBlockInfo(ISLocalToGlobalMapping,PetscInt*,PetscInt*[],PetscInt*[],PetscInt**[]);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingRestoreBlockInfo(ISLocalToGlobalMapping,PetscInt*,PetscInt*[],PetscInt*[],PetscInt**[]);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingGetIndices(ISLocalToGlobalMapping,const PetscInt**);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingRestoreIndices(ISLocalToGlobalMapping,const PetscInt**);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingGetBlockIndices(ISLocalToGlobalMapping,const PetscInt**);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingRestoreBlockIndices(ISLocalToGlobalMapping,const PetscInt**);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingConcatenate(MPI_Comm,PetscInt,const ISLocalToGlobalMapping[],ISLocalToGlobalMapping*);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingGetBlockSize(ISLocalToGlobalMapping,PetscInt*);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingSetBlockSize(ISLocalToGlobalMapping,PetscInt);
PETSC_EXTERN PetscErrorCode ISLocalToGlobalMappingDuplicate(ISLocalToGlobalMapping,ISLocalToGlobalMapping*);


/* --------------------------------------------------------------------------*/
/*E
    ISColoringType - determines if the coloring is for the entire parallel grid/graph/matrix
                     or for just the local ghosted portion

    Level: beginner

$   IS_COLORING_GLOBAL - does not include the colors for ghost points, this is used when the function
$                        is called synchronously in parallel. This requires generating a "parallel coloring".
$   IS_COLORING_LOCAL - includes colors for ghost points, this is used when the function can be called
$                         separately on individual processes with the ghost points already filled in. Does not
$                         require a "parallel coloring", rather each process colors its local + ghost part.
$                         Using this can result in much less parallel communication. Currently only works 
$                         with DMDA and if you call MatFDColoringSetFunction() with the local function.

.seealso: DMCreateColoring()
E*/
typedef enum {IS_COLORING_GLOBAL,IS_COLORING_LOCAL} ISColoringType;
PETSC_EXTERN const char *const ISColoringTypes[];
typedef unsigned PETSC_IS_COLOR_VALUE_TYPE ISColoringValue;
PETSC_EXTERN PetscErrorCode ISAllGatherColors(MPI_Comm,PetscInt,ISColoringValue*,PetscInt*,ISColoringValue*[]);

PETSC_EXTERN PetscErrorCode ISColoringCreate(MPI_Comm,PetscInt,PetscInt,const ISColoringValue[],PetscCopyMode,ISColoring*);
PETSC_EXTERN PetscErrorCode ISColoringDestroy(ISColoring*);
PETSC_EXTERN PetscErrorCode ISColoringView(ISColoring,PetscViewer);
PETSC_EXTERN PetscErrorCode ISColoringViewFromOptions(ISColoring,PetscObject,const char[]);
PETSC_EXTERN PetscErrorCode ISColoringGetIS(ISColoring,PetscInt*,IS*[]);
PETSC_EXTERN PetscErrorCode ISColoringRestoreIS(ISColoring,IS*[]);
PETSC_EXTERN PetscErrorCode ISColoringReference(ISColoring);
PETSC_EXTERN PetscErrorCode ISColoringSetType(ISColoring,ISColoringType);


/* --------------------------------------------------------------------------*/
PETSC_EXTERN PetscErrorCode ISBuildTwoSided(IS,IS,IS*);
PETSC_EXTERN PetscErrorCode ISPartitioningToNumbering(IS,IS*);
PETSC_EXTERN PetscErrorCode ISPartitioningCount(IS,PetscInt,PetscInt[]);

PETSC_EXTERN PetscErrorCode ISCompressIndicesGeneral(PetscInt,PetscInt,PetscInt,PetscInt,const IS[],IS[]);
PETSC_EXTERN PetscErrorCode ISCompressIndicesSorted(PetscInt,PetscInt,PetscInt,const IS[],IS[]);
PETSC_EXTERN PetscErrorCode ISExpandIndicesGeneral(PetscInt,PetscInt,PetscInt,PetscInt,const IS[],IS[]);


struct _n_PetscLayout{
  MPI_Comm               comm;
  PetscInt               n,N;         /* local, global vector size */
  PetscInt               rstart,rend; /* local start, local end + 1 */
  PetscInt               *range;      /* the offset of each processor */
  PetscInt               bs;          /* number of elements in each block (generally for multi-component
                                       * problems). Defaults to -1 and can be arbitrarily lazy so always use
                                       * PetscAbs(map->bs) when accessing directly and expecting result to be
                                       * positive. Do NOT multiply above numbers by bs */
  PetscInt               refcnt;      /* MPI Vecs obtained with VecDuplicate() and from MatCreateVecs() reuse map of input object */
  ISLocalToGlobalMapping mapping;     /* mapping used in Vec/MatSetValuesLocal() */
};

/*@C
     PetscLayoutFindOwner - Find the owning rank for a global index

    Not Collective

   Input Parameters:
+    map - the layout
-    idx - global index to find the owner of

   Output Parameter:
.    owner - the owning rank

   Level: developer

    Fortran Notes:
      Not available from Fortran

@*/
PETSC_STATIC_INLINE PetscErrorCode PetscLayoutFindOwner(PetscLayout map,PetscInt idx,PetscInt *owner)
{
  PetscErrorCode ierr;
  PetscMPIInt    lo = 0,hi,t;

  PetscFunctionBegin;
  *owner = -1;                  /* GCC erroneously issues warning about possibly uninitialized use when error condition */
  if (!((map->n >= 0) && (map->N >= 0) && (map->range))) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"PetscLayoutSetUp() must be called first");
  if (idx < 0 || idx > map->N) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Index %D is out of range",idx);
  ierr = MPI_Comm_size(map->comm,&hi);CHKERRQ(ierr);
  while (hi - lo > 1) {
    t = lo + (hi - lo) / 2;
    if (idx < map->range[t]) hi = t;
    else                     lo = t;
  }
  *owner = lo;
  PetscFunctionReturn(0);
}

/*@C
     PetscLayoutFindOwnerIndex - Find the owning rank and the local index for a global index

    Not Collective

   Input Parameters:
+    map   - the layout
-    idx   - global index to find the owner of

   Output Parameter:
+    owner - the owning rank
-    lidx  - local index used by the owner for idx

   Level: developer

    Fortran Notes:
      Not available from Fortran

@*/
PETSC_STATIC_INLINE PetscErrorCode PetscLayoutFindOwnerIndex(PetscLayout map,PetscInt idx,PetscInt *owner, PetscInt *lidx)
{
  PetscErrorCode ierr;
  PetscMPIInt    lo = 0,hi,t;

  PetscFunctionBegin;
  if (!((map->n >= 0) && (map->N >= 0) && (map->range))) SETERRQ(PETSC_COMM_SELF,PETSC_ERR_ARG_WRONGSTATE,"PetscLayoutSetUp() must be called first");
  if (idx < 0 || idx > map->N) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_ARG_OUTOFRANGE,"Index %D is out of range",idx);
  ierr = MPI_Comm_size(map->comm,&hi);CHKERRQ(ierr);
  while (hi - lo > 1) {
    t = lo + (hi - lo) / 2;
    if (idx < map->range[t]) hi = t;
    else                     lo = t;
  }
  if (owner) *owner = lo;
  if (lidx) *lidx  = idx-map->range[lo];
  PetscFunctionReturn(0);
}

PETSC_EXTERN PetscErrorCode PetscLayoutCreate(MPI_Comm,PetscLayout*);
PETSC_EXTERN PetscErrorCode PetscLayoutSetUp(PetscLayout);
PETSC_EXTERN PetscErrorCode PetscLayoutDestroy(PetscLayout*);
PETSC_EXTERN PetscErrorCode PetscLayoutDuplicate(PetscLayout,PetscLayout*);
PETSC_EXTERN PetscErrorCode PetscLayoutReference(PetscLayout,PetscLayout*);
PETSC_EXTERN PetscErrorCode PetscLayoutSetLocalSize(PetscLayout,PetscInt);
PETSC_EXTERN PetscErrorCode PetscLayoutGetLocalSize(PetscLayout,PetscInt *);
PETSC_EXTERN PetscErrorCode PetscLayoutSetSize(PetscLayout,PetscInt);
PETSC_EXTERN PetscErrorCode PetscLayoutGetSize(PetscLayout,PetscInt *);
PETSC_EXTERN PetscErrorCode PetscLayoutSetBlockSize(PetscLayout,PetscInt);
PETSC_EXTERN PetscErrorCode PetscLayoutGetBlockSize(PetscLayout,PetscInt*);
PETSC_EXTERN PetscErrorCode PetscLayoutGetRange(PetscLayout,PetscInt *,PetscInt *);
PETSC_EXTERN PetscErrorCode PetscLayoutGetRanges(PetscLayout,const PetscInt *[]);
PETSC_EXTERN PetscErrorCode PetscLayoutCompare(PetscLayout,PetscLayout,PetscBool*);
PETSC_EXTERN PetscErrorCode PetscLayoutSetISLocalToGlobalMapping(PetscLayout,ISLocalToGlobalMapping);
PETSC_EXTERN PetscErrorCode PetscSFSetGraphLayout(PetscSF,PetscLayout,PetscInt,const PetscInt*,PetscCopyMode,const PetscInt*);

PETSC_EXTERN PetscClassId PETSC_SECTION_CLASSID;

PETSC_EXTERN PetscErrorCode PetscSectionCreate(MPI_Comm,PetscSection*);
PETSC_EXTERN PetscErrorCode PetscSectionClone(PetscSection, PetscSection*);
PETSC_EXTERN PetscErrorCode PetscSectionCopy(PetscSection, PetscSection);
PETSC_EXTERN PetscErrorCode PetscSectionCompare(PetscSection, PetscSection, PetscBool*);
PETSC_EXTERN PetscErrorCode PetscSectionGetNumFields(PetscSection, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscSectionSetNumFields(PetscSection, PetscInt);
PETSC_EXTERN PetscErrorCode PetscSectionGetFieldName(PetscSection, PetscInt, const char *[]);
PETSC_EXTERN PetscErrorCode PetscSectionSetFieldName(PetscSection, PetscInt, const char []);
PETSC_EXTERN PetscErrorCode PetscSectionGetFieldComponents(PetscSection, PetscInt, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscSectionSetFieldComponents(PetscSection, PetscInt, PetscInt);
PETSC_EXTERN PetscErrorCode PetscSectionGetChart(PetscSection, PetscInt *, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscSectionSetChart(PetscSection, PetscInt, PetscInt);
PETSC_EXTERN PetscErrorCode PetscSectionGetPermutation(PetscSection, IS *);
PETSC_EXTERN PetscErrorCode PetscSectionSetPermutation(PetscSection, IS);
PETSC_EXTERN PetscErrorCode PetscSectionGetDof(PetscSection, PetscInt, PetscInt*);
PETSC_EXTERN PetscErrorCode PetscSectionSetDof(PetscSection, PetscInt, PetscInt);
PETSC_EXTERN PetscErrorCode PetscSectionAddDof(PetscSection, PetscInt, PetscInt);
PETSC_EXTERN PetscErrorCode PetscSectionGetFieldDof(PetscSection, PetscInt, PetscInt, PetscInt*);
PETSC_EXTERN PetscErrorCode PetscSectionSetFieldDof(PetscSection, PetscInt, PetscInt, PetscInt);
PETSC_EXTERN PetscErrorCode PetscSectionAddFieldDof(PetscSection, PetscInt, PetscInt, PetscInt);
PETSC_EXTERN PetscErrorCode PetscSectionHasConstraints(PetscSection, PetscBool *);
PETSC_EXTERN PetscErrorCode PetscSectionGetConstraintDof(PetscSection, PetscInt, PetscInt*);
PETSC_EXTERN PetscErrorCode PetscSectionSetConstraintDof(PetscSection, PetscInt, PetscInt);
PETSC_EXTERN PetscErrorCode PetscSectionAddConstraintDof(PetscSection, PetscInt, PetscInt);
PETSC_EXTERN PetscErrorCode PetscSectionGetFieldConstraintDof(PetscSection, PetscInt, PetscInt, PetscInt*);
PETSC_EXTERN PetscErrorCode PetscSectionSetFieldConstraintDof(PetscSection, PetscInt, PetscInt, PetscInt);
PETSC_EXTERN PetscErrorCode PetscSectionAddFieldConstraintDof(PetscSection, PetscInt, PetscInt, PetscInt);
PETSC_EXTERN PetscErrorCode PetscSectionGetConstraintIndices(PetscSection, PetscInt, const PetscInt**);
PETSC_EXTERN PetscErrorCode PetscSectionSetConstraintIndices(PetscSection, PetscInt, const PetscInt*);
PETSC_EXTERN PetscErrorCode PetscSectionGetFieldConstraintIndices(PetscSection, PetscInt, PetscInt, const PetscInt**);
PETSC_EXTERN PetscErrorCode PetscSectionSetFieldConstraintIndices(PetscSection, PetscInt, PetscInt, const PetscInt*);
PETSC_EXTERN PetscErrorCode PetscSectionSetUpBC(PetscSection);
PETSC_EXTERN PetscErrorCode PetscSectionSetUp(PetscSection);
PETSC_EXTERN PetscErrorCode PetscSectionGetMaxDof(PetscSection, PetscInt*);
PETSC_EXTERN PetscErrorCode PetscSectionGetStorageSize(PetscSection, PetscInt*);
PETSC_EXTERN PetscErrorCode PetscSectionGetConstrainedStorageSize(PetscSection, PetscInt*);
PETSC_EXTERN PetscErrorCode PetscSectionGetOffset(PetscSection, PetscInt, PetscInt*);
PETSC_EXTERN PetscErrorCode PetscSectionSetOffset(PetscSection, PetscInt, PetscInt);
PETSC_EXTERN PetscErrorCode PetscSectionGetFieldOffset(PetscSection, PetscInt, PetscInt, PetscInt*);
PETSC_EXTERN PetscErrorCode PetscSectionSetFieldOffset(PetscSection, PetscInt, PetscInt, PetscInt);
PETSC_EXTERN PetscErrorCode PetscSectionGetOffsetRange(PetscSection, PetscInt *, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscSectionView(PetscSection, PetscViewer);
PETSC_STATIC_INLINE PetscErrorCode PetscSectionViewFromOptions(PetscSection A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}
PETSC_EXTERN PetscErrorCode PetscSectionReset(PetscSection);
PETSC_EXTERN PetscErrorCode PetscSectionDestroy(PetscSection*);
PETSC_EXTERN PetscErrorCode PetscSectionCreateGlobalSection(PetscSection, PetscSF, PetscBool, PetscBool, PetscSection *);
PETSC_EXTERN PetscErrorCode PetscSectionCreateGlobalSectionCensored(PetscSection, PetscSF, PetscBool, PetscInt, const PetscInt [], PetscSection *);
PETSC_EXTERN PetscErrorCode PetscSectionCreateSubsection(PetscSection, PetscInt, const PetscInt [], PetscSection *);
PETSC_EXTERN PetscErrorCode PetscSectionCreateSupersection(PetscSection[], PetscInt, PetscSection *);
PETSC_EXTERN PetscErrorCode PetscSectionCreateSubmeshSection(PetscSection, IS, PetscSection *);
PETSC_EXTERN PetscErrorCode PetscSectionGetPointLayout(MPI_Comm, PetscSection, PetscLayout *);
PETSC_EXTERN PetscErrorCode PetscSectionGetValueLayout(MPI_Comm, PetscSection, PetscLayout *);
PETSC_EXTERN PetscErrorCode PetscSectionPermute(PetscSection, IS, PetscSection *);
PETSC_EXTERN PetscErrorCode PetscSectionGetField(PetscSection, PetscInt, PetscSection *);
PETSC_EXTERN PetscErrorCode PetscSectionSetUseFieldOffsets(PetscSection, PetscBool);
PETSC_EXTERN PetscErrorCode PetscSectionGetUseFieldOffsets(PetscSection, PetscBool *);

PETSC_EXTERN PetscErrorCode PetscSectionSetClosureIndex(PetscSection, PetscObject, PetscSection, IS);
PETSC_EXTERN PetscErrorCode PetscSectionGetClosureIndex(PetscSection, PetscObject, PetscSection *, IS *);
PETSC_EXTERN PetscErrorCode PetscSectionSetClosurePermutation(PetscSection, PetscObject, IS);
PETSC_EXTERN PetscErrorCode PetscSectionGetClosurePermutation(PetscSection, PetscObject, IS *);
PETSC_EXTERN PetscErrorCode PetscSectionGetClosureInversePermutation(PetscSection, PetscObject, IS *);

PETSC_EXTERN PetscClassId PETSC_SECTION_SYM_CLASSID;

/*J
  PetscSectionSymType - String with the name of a PetscSectionSym type.

  Level: developer

  Notes:
    PetscSectionSym has no default implementation, but is used by DM in PetscSectionSymCreateLabel().

.seealso: PetscSectionSymSetType(), PetscSectionSym, PetscSectionSymCreate(), PetscSectionSymRegister()
J*/
typedef const char *PetscSectionSymType;

PETSC_EXTERN PetscFunctionList PetscSectionSymList;
PETSC_EXTERN PetscErrorCode PetscSectionSymSetType(PetscSectionSym, PetscSectionSymType);
PETSC_EXTERN PetscErrorCode PetscSectionSymGetType(PetscSectionSym, PetscSectionSymType*);
PETSC_EXTERN PetscErrorCode PetscSectionSymRegister(const char[],PetscErrorCode (*)(PetscSectionSym));

PETSC_EXTERN PetscErrorCode PetscSectionSymCreate(MPI_Comm, PetscSectionSym*);
PETSC_EXTERN PetscErrorCode PetscSectionSymDestroy(PetscSectionSym*);
PETSC_EXTERN PetscErrorCode PetscSectionSymView(PetscSectionSym,PetscViewer);

PETSC_EXTERN PetscErrorCode PetscSectionSetSym(PetscSection, PetscSectionSym);
PETSC_EXTERN PetscErrorCode PetscSectionGetSym(PetscSection, PetscSectionSym*);
PETSC_EXTERN PetscErrorCode PetscSectionSetFieldSym(PetscSection, PetscInt, PetscSectionSym);
PETSC_EXTERN PetscErrorCode PetscSectionGetFieldSym(PetscSection, PetscInt, PetscSectionSym*);

PETSC_EXTERN PetscErrorCode PetscSectionGetPointSyms(PetscSection, PetscInt, const PetscInt *, const PetscInt ***, const PetscScalar ***);
PETSC_EXTERN PetscErrorCode PetscSectionRestorePointSyms(PetscSection, PetscInt, const PetscInt *, const PetscInt ***, const PetscScalar ***);
PETSC_EXTERN PetscErrorCode PetscSectionGetFieldPointSyms(PetscSection, PetscInt, PetscInt, const PetscInt *, const PetscInt ***, const PetscScalar ***);
PETSC_EXTERN PetscErrorCode PetscSectionRestoreFieldPointSyms(PetscSection, PetscInt, PetscInt, const PetscInt *, const PetscInt ***, const PetscScalar ***);


/* PetscSF support */
PETSC_EXTERN PetscErrorCode PetscSFConvertPartition(PetscSF, PetscSection, IS, ISLocalToGlobalMapping *, PetscSF *);
PETSC_EXTERN PetscErrorCode PetscSFCreateRemoteOffsets(PetscSF, PetscSection, PetscSection, PetscInt **);
PETSC_EXTERN PetscErrorCode PetscSFDistributeSection(PetscSF, PetscSection, PetscInt **, PetscSection);
PETSC_EXTERN PetscErrorCode PetscSFCreateSectionSF(PetscSF, PetscSection, PetscInt [], PetscSection, PetscSF *);

#endif
