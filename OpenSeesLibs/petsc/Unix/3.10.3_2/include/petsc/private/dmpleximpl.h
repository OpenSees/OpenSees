#if !defined(_PLEXIMPL_H)
#define _PLEXIMPL_H

#include <petscmat.h>       /*I      "petscmat.h"          I*/
#include <petscdmplex.h> /*I      "petscdmplex.h"    I*/
#include <petscbt.h>
#include <petscsf.h>
#include <petsc/private/dmimpl.h>

PETSC_EXTERN PetscLogEvent DMPLEX_Interpolate;
PETSC_EXTERN PetscLogEvent DMPLEX_Partition;
PETSC_EXTERN PetscLogEvent DMPLEX_Distribute;
PETSC_EXTERN PetscLogEvent DMPLEX_DistributeCones;
PETSC_EXTERN PetscLogEvent DMPLEX_DistributeLabels;
PETSC_EXTERN PetscLogEvent DMPLEX_DistributeSF;
PETSC_EXTERN PetscLogEvent DMPLEX_DistributeOverlap;
PETSC_EXTERN PetscLogEvent DMPLEX_DistributeField;
PETSC_EXTERN PetscLogEvent DMPLEX_DistributeData;
PETSC_EXTERN PetscLogEvent DMPLEX_Migrate;
PETSC_EXTERN PetscLogEvent DMPLEX_InterpolateSF;
PETSC_EXTERN PetscLogEvent DMPLEX_GlobalToNaturalBegin;
PETSC_EXTERN PetscLogEvent DMPLEX_GlobalToNaturalEnd;
PETSC_EXTERN PetscLogEvent DMPLEX_NaturalToGlobalBegin;
PETSC_EXTERN PetscLogEvent DMPLEX_NaturalToGlobalEnd;
PETSC_EXTERN PetscLogEvent DMPLEX_Stratify;
PETSC_EXTERN PetscLogEvent DMPLEX_Preallocate;
PETSC_EXTERN PetscLogEvent DMPLEX_ResidualFEM;
PETSC_EXTERN PetscLogEvent DMPLEX_JacobianFEM;
PETSC_EXTERN PetscLogEvent DMPLEX_InterpolatorFEM;
PETSC_EXTERN PetscLogEvent DMPLEX_InjectorFEM;
PETSC_EXTERN PetscLogEvent DMPLEX_IntegralFEM;
PETSC_EXTERN PetscLogEvent DMPLEX_CreateGmsh;

PETSC_EXTERN PetscBool      PetscPartitionerRegisterAllCalled;
PETSC_EXTERN PetscErrorCode PetscPartitionerRegisterAll(void);

typedef enum {REFINER_NOOP = 0,
              REFINER_SIMPLEX_1D,
              REFINER_SIMPLEX_2D,
              REFINER_HYBRID_SIMPLEX_2D,
              REFINER_SIMPLEX_TO_HEX_2D,
              REFINER_HYBRID_SIMPLEX_TO_HEX_2D,
              REFINER_HEX_2D,
              REFINER_HYBRID_HEX_2D,
              REFINER_SIMPLEX_3D,
              REFINER_HYBRID_SIMPLEX_3D,
              REFINER_SIMPLEX_TO_HEX_3D,
              REFINER_HYBRID_SIMPLEX_TO_HEX_3D,
              REFINER_HEX_3D,
              REFINER_HYBRID_HEX_3D} CellRefiner;

typedef struct _PetscPartitionerOps *PetscPartitionerOps;
struct _PetscPartitionerOps {
  PetscErrorCode (*setfromoptions)(PetscOptionItems*,PetscPartitioner);
  PetscErrorCode (*setup)(PetscPartitioner);
  PetscErrorCode (*view)(PetscPartitioner,PetscViewer);
  PetscErrorCode (*destroy)(PetscPartitioner);
  PetscErrorCode (*partition)(PetscPartitioner, DM, PetscInt, PetscInt, PetscInt[], PetscInt[], PetscSection, IS *);
};

struct _p_PetscPartitioner {
  PETSCHEADER(struct _PetscPartitionerOps);
  void           *data;             /* Implementation object */
  PetscInt        height;           /* Height of points to partition into non-overlapping subsets */
};

typedef struct {
  PetscInt dummy;
} PetscPartitioner_Chaco;

typedef struct {
  PetscInt ptype;
} PetscPartitioner_ParMetis;

typedef struct {
  PetscSection section;   /* Sizes for each partition */
  IS           partition; /* Points in each partition */
  PetscBool    random;    /* Flag for a random partition */
} PetscPartitioner_Shell;

typedef struct {
  PetscInt dummy;
} PetscPartitioner_Simple;

typedef struct {
  PetscInt dummy;
} PetscPartitioner_Gather;

/* Utility struct to store the contents of a Gmsh file in memory */
typedef struct {
  PetscInt dim;      /* Entity dimension */
  PetscInt id;       /* Element number */
  PetscInt numNodes; /* Size of node array */
  int nodes[8];      /* Node array */
  PetscInt numTags;  /* Size of tag array */
  int tags[4];       /* Tag array */
  PetscInt cellType; /* Cell type */
} GmshElement;

/* Utility struct to store the contents of a Fluent file in memory */
typedef struct {
  int          index;    /* Type of section */
  unsigned int zoneID;
  unsigned int first;
  unsigned int last;
  int          type;
  int          nd;       /* Either ND or element-type */
  void        *data;
} FluentSection;

struct _PetscGridHash {
  PetscInt     dim;
  PetscReal    lower[3];    /* The lower-left corner */
  PetscReal    upper[3];    /* The upper-right corner */
  PetscReal    extent[3];   /* The box size */
  PetscReal    h[3];        /* The subbox size */
  PetscInt     n[3];        /* The number of subboxes */
  PetscSection cellSection; /* Offsets for cells in each subbox*/
  IS           cells;       /* List of cells in each subbox */
  DMLabel      cellsSparse; /* Sparse storage for cell map */
};

typedef struct {
  PetscInt             refct;

  /* Sieve */
  PetscSection         coneSection;       /* Layout of cones (inedges for DAG) */
  PetscInt             maxConeSize;       /* Cached for fast lookup */
  PetscInt            *cones;             /* Cone for each point */
  PetscInt            *coneOrientations;  /* Orientation of each cone point, means cone traveral should start on point 'o', and if negative start on -(o+1) and go in reverse */
  PetscSection         supportSection;    /* Layout of cones (inedges for DAG) */
  PetscInt             maxSupportSize;    /* Cached for fast lookup */
  PetscInt            *supports;          /* Cone for each point */
  PetscBool            refinementUniform; /* Flag for uniform cell refinement */
  PetscReal            refinementLimit;   /* Maximum volume for refined cell */
  PetscErrorCode     (*refinementFunc)(const PetscReal [], PetscReal *); /* Function giving the maximum volume for refined cell */
  PetscInt             hybridPointMax[8]; /* Allow segregation of some points, each dimension has a divider (used in VTK output and refinement) */

  PetscInt            *facesTmp;          /* Work space for faces operation */

  /* Hierarchy */
  PetscBool            regularRefinement; /* This flag signals that we are a regular refinement of coarseMesh */

  /* Generation */
  char                *tetgenOpts;
  char                *triangleOpts;
  PetscPartitioner     partitioner;
  PetscBool            partitionBalance;  /* Evenly divide partition overlap when distributing */
  PetscBool            remeshBd;

  /* Submesh */
  DMLabel              subpointMap;       /* Label each original mesh point in the submesh with its depth, subpoint are the implicit numbering */

  /* Labels and numbering */
  PetscObjectState     depthState;        /* State of depth label, so that we can determine if a user changes it */
  IS                   globalVertexNumbers;
  IS                   globalCellNumbers;

  /* Constraints */
  PetscSection         anchorSection;      /* maps constrained points to anchor points */
  IS                   anchorIS;           /* anchors indexed by the above section */
  PetscErrorCode     (*createanchors)(DM); /* automatically compute anchors (probably from tree constraints) */
  PetscErrorCode     (*computeanchormatrix)(DM,PetscSection,PetscSection,Mat);

  /* Tree: automatically construct constraints for hierarchically non-conforming meshes */
  PetscSection         parentSection;     /* dof == 1 if point has parent */
  PetscInt            *parents;           /* point to parent */
  PetscInt            *childIDs;          /* point to child ID */
  PetscSection         childSection;      /* inverse of parent section */
  PetscInt            *children;          /* point to children */
  DM                   referenceTree;     /* reference tree to which child ID's refer */
  PetscErrorCode      (*getchildsymmetry)(DM,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt,PetscInt*,PetscInt*);

  /* MATIS support */
  PetscSection         subdomainSection;

  /* Adjacency */
  PetscBool            useAnchors;        /* Replace constrained points with their anchors in adjacency lists */
  PetscErrorCode      (*useradjacency)(DM,PetscInt,PetscInt*,PetscInt[],void*); /* User callback for adjacency */
  void                *useradjacencyctx;  /* User context for callback */

  /* Projection */
  PetscInt             maxProjectionHeight; /* maximum height of cells used in DMPlexProject functions */

  /* Output */
  PetscInt             vtkCellHeight;            /* The height of cells for output, default is 0 */
  PetscReal            scale[NUM_PETSC_UNITS];   /* The scale for each SI unit */

  /* Geometry */
  PetscReal            minradius;         /* Minimum distance from cell centroid to face */
  PetscBool            useHashLocation;   /* Use grid hashing for point location */
  PetscGridHash        lbox;              /* Local box for searching */

  /* Debugging */
  PetscBool            printSetValues;
  PetscInt             printFEM;
  PetscInt             printL2;
  PetscReal            printTol;
} DM_Plex;

PETSC_EXTERN PetscErrorCode DMPlexVTKWriteAll_VTU(DM,PetscViewer);
PETSC_EXTERN PetscErrorCode VecView_Plex_Local(Vec,PetscViewer);
PETSC_EXTERN PetscErrorCode VecView_Plex_Native(Vec,PetscViewer);
PETSC_EXTERN PetscErrorCode VecView_Plex(Vec,PetscViewer);
PETSC_EXTERN PetscErrorCode VecLoad_Plex_Local(Vec,PetscViewer);
PETSC_EXTERN PetscErrorCode VecLoad_Plex_Native(Vec,PetscViewer);
PETSC_EXTERN PetscErrorCode VecLoad_Plex(Vec,PetscViewer);
PETSC_INTERN PetscErrorCode DMPlexGetFieldType_Internal(DM, PetscSection, PetscInt, PetscInt *, PetscInt *, PetscViewerVTKFieldType *);
PETSC_INTERN PetscErrorCode DMPlexView_GLVis(DM,PetscViewer);
PETSC_INTERN PetscErrorCode DMSetUpGLVisViewer_Plex(PetscObject,PetscViewer);
#if defined(PETSC_HAVE_HDF5)
PETSC_EXTERN PetscErrorCode VecView_Plex_Local_HDF5(Vec, PetscViewer);
PETSC_EXTERN PetscErrorCode VecView_Plex_HDF5(Vec, PetscViewer);
PETSC_EXTERN PetscErrorCode VecLoad_Plex_HDF5(Vec, PetscViewer);
PETSC_EXTERN PetscErrorCode VecView_Plex_HDF5_Native(Vec, PetscViewer);
PETSC_EXTERN PetscErrorCode VecLoad_Plex_HDF5_Native(Vec, PetscViewer);
PETSC_EXTERN PetscErrorCode DMPlexView_HDF5(DM, PetscViewer);
PETSC_EXTERN PetscErrorCode DMPlexLoad_HDF5(DM, PetscViewer);
#endif

PETSC_INTERN PetscErrorCode DMSetFromOptions_NonRefinement_Plex(PetscOptionItems *, DM);
PETSC_INTERN PetscErrorCode DMCoarsen_Plex(DM, MPI_Comm, DM *);
PETSC_INTERN PetscErrorCode DMCoarsenHierarchy_Plex(DM, PetscInt, DM []);
PETSC_INTERN PetscErrorCode DMRefine_Plex(DM, MPI_Comm, DM *);
PETSC_INTERN PetscErrorCode DMRefineHierarchy_Plex(DM, PetscInt, DM []);
PETSC_INTERN PetscErrorCode DMAdaptLabel_Plex(DM, DMLabel, DM *);
PETSC_INTERN PetscErrorCode DMAdaptMetric_Plex(DM, Vec, DMLabel, DM *);
PETSC_INTERN PetscErrorCode DMPlexInsertBoundaryValues_Plex(DM, PetscBool, Vec, PetscReal, Vec, Vec, Vec);
PETSC_INTERN PetscErrorCode DMProjectFunctionLocal_Plex(DM,PetscReal,PetscErrorCode(**)(PetscInt,PetscReal,const PetscReal[],PetscInt,PetscScalar *,void *),void **,InsertMode,Vec);
PETSC_INTERN PetscErrorCode DMProjectFunctionLabelLocal_Plex(DM,PetscReal,DMLabel,PetscInt,const PetscInt[],PetscInt,const PetscInt[],PetscErrorCode(**)(PetscInt,PetscReal,const PetscReal[],PetscInt,PetscScalar *,void *),void **,InsertMode,Vec);
PETSC_INTERN PetscErrorCode DMProjectFieldLocal_Plex(DM,PetscReal,Vec,void (**)(PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[],const PetscScalar[],const PetscScalar[],const PetscInt[],const PetscInt[],const PetscScalar[],const PetscScalar[],const PetscScalar[],PetscReal,const PetscReal[],PetscInt,const PetscScalar[],PetscScalar[]),InsertMode,Vec);
PETSC_INTERN PetscErrorCode DMProjectFieldLabelLocal_Plex(DM,PetscReal,DMLabel,PetscInt,const PetscInt[],PetscInt,const PetscInt[],Vec,void (**)(PetscInt,PetscInt,PetscInt,const PetscInt[],const PetscInt[],const PetscScalar[],const PetscScalar[],const PetscScalar[],const PetscInt[],const PetscInt[],const PetscScalar[],const PetscScalar[],const PetscScalar[],PetscReal,const PetscReal[],PetscInt,const PetscScalar[],PetscScalar[]),InsertMode,Vec);
PETSC_INTERN PetscErrorCode DMComputeL2Diff_Plex(DM,PetscReal,PetscErrorCode(**)(PetscInt,PetscReal,const PetscReal[],PetscInt,PetscScalar *,void *),void **,Vec,PetscReal *);
PETSC_INTERN PetscErrorCode DMComputeL2GradientDiff_Plex(DM,PetscReal,PetscErrorCode(**)(PetscInt,PetscReal,const PetscReal[], const PetscReal[],PetscInt,PetscScalar *,void *),void **,Vec,const PetscReal [],PetscReal *);
PETSC_INTERN PetscErrorCode DMComputeL2FieldDiff_Plex(DM,PetscReal,PetscErrorCode(**)(PetscInt,PetscReal,const PetscReal[],PetscInt,PetscScalar *,void *),void **,Vec,PetscReal *);
PETSC_INTERN PetscErrorCode DMLocatePoints_Plex(DM, Vec, DMPointLocationType, PetscSF);

PETSC_INTERN PetscErrorCode DMPlexBuildFromCellList_Internal(DM, PetscInt, PetscInt, PetscInt, PetscInt, const int[], PetscBool);
PETSC_INTERN PetscErrorCode DMPlexBuildFromCellList_Parallel_Internal(DM, PetscInt, PetscInt, PetscInt, PetscInt, const int[], PetscBool, PetscSF *);
PETSC_INTERN PetscErrorCode DMPlexBuildCoordinates_Internal(DM, PetscInt, PetscInt, PetscInt, const double[]);
PETSC_INTERN PetscErrorCode DMPlexBuildCoordinates_Parallel_Internal(DM, PetscInt, PetscInt, PetscInt, PetscSF, const PetscReal[]);
PETSC_INTERN PetscErrorCode DMPlexLoadLabels_HDF5_Internal(DM, PetscViewer);
PETSC_INTERN PetscErrorCode DMPlexView_HDF5_Internal(DM, PetscViewer);
PETSC_INTERN PetscErrorCode DMPlexLoad_HDF5_Internal(DM, PetscViewer);
PETSC_INTERN PetscErrorCode DMPlexLoad_HDF5_Xdmf_Internal(DM, PetscViewer);
PETSC_INTERN PetscErrorCode VecView_Plex_HDF5_Internal(Vec, PetscViewer);
PETSC_INTERN PetscErrorCode VecView_Plex_HDF5_Native_Internal(Vec, PetscViewer);
PETSC_INTERN PetscErrorCode VecView_Plex_Local_HDF5_Internal(Vec, PetscViewer);
PETSC_INTERN PetscErrorCode VecLoad_Plex_HDF5_Internal(Vec, PetscViewer);
PETSC_INTERN PetscErrorCode VecLoad_Plex_HDF5_Native_Internal(Vec, PetscViewer);
/* TODO Make these INTERN */
PETSC_EXTERN PetscErrorCode DMPlexView_ExodusII_Internal(DM, int, PetscInt);
PETSC_EXTERN PetscErrorCode VecViewPlex_ExodusII_Nodal_Internal(Vec, int, int);
PETSC_EXTERN PetscErrorCode VecLoadPlex_ExodusII_Nodal_Internal(Vec, int, int);
PETSC_EXTERN PetscErrorCode VecViewPlex_ExodusII_Zonal_Internal(Vec, int, int);
PETSC_EXTERN PetscErrorCode VecLoadPlex_ExodusII_Zonal_Internal(Vec, int, int);
PETSC_INTERN PetscErrorCode DMPlexVTKGetCellType_Internal(DM,PetscInt,PetscInt,PetscInt*);
PETSC_INTERN PetscErrorCode DMPlexGetAdjacency_Internal(DM,PetscInt,PetscBool,PetscBool,PetscBool,PetscInt*,PetscInt*[]);
PETSC_INTERN PetscErrorCode DMPlexGetFaces_Internal(DM,PetscInt,PetscInt,PetscInt*,PetscInt*,const PetscInt*[]);
PETSC_INTERN PetscErrorCode DMPlexGetRawFaces_Internal(DM,PetscInt,PetscInt,const PetscInt[], PetscInt*,PetscInt*,const PetscInt*[]);
PETSC_INTERN PetscErrorCode DMPlexRestoreFaces_Internal(DM,PetscInt,PetscInt,PetscInt*,PetscInt*,const PetscInt*[]);
PETSC_INTERN PetscErrorCode DMPlexRefineUniform_Internal(DM,CellRefiner,DM*);
PETSC_INTERN PetscErrorCode DMPlexGetCellRefiner_Internal(DM,CellRefiner*);
PETSC_INTERN PetscErrorCode CellRefinerGetAffineTransforms_Internal(CellRefiner, PetscInt *, PetscReal *[], PetscReal *[], PetscReal *[]);
PETSC_INTERN PetscErrorCode CellRefinerRestoreAffineTransforms_Internal(CellRefiner, PetscInt *, PetscReal *[], PetscReal *[], PetscReal *[]);
PETSC_INTERN PetscErrorCode CellRefinerInCellTest_Internal(CellRefiner, const PetscReal[], PetscBool *);
PETSC_INTERN PetscErrorCode DMPlexInvertCell_Internal(PetscInt, PetscInt, PetscInt[]);
PETSC_INTERN PetscErrorCode DMPlexVecSetFieldClosure_Internal(DM, PetscSection, Vec, PetscBool[], PetscInt, PetscInt, const PetscInt[], const PetscScalar[], InsertMode);
PETSC_INTERN PetscErrorCode DMPlexProjectConstraints_Internal(DM, Vec, Vec);
PETSC_EXTERN PetscErrorCode DMPlexCreateReferenceTree_Union(DM,DM,const char *,DM*);
PETSC_EXTERN PetscErrorCode DMPlexComputeInterpolatorTree(DM,DM,PetscSF,PetscInt *,Mat);
PETSC_EXTERN PetscErrorCode DMPlexComputeInjectorTree(DM,DM,PetscSF,PetscInt *,Mat);
PETSC_EXTERN PetscErrorCode DMPlexAnchorsModifyMat(DM,PetscSection,PetscInt,PetscInt,const PetscInt[],const PetscInt ***,const PetscScalar[],PetscInt*,PetscInt*,PetscInt*[],PetscScalar*[],PetscInt[],PetscBool);
PETSC_EXTERN PetscErrorCode indicesPoint_private(PetscSection,PetscInt,PetscInt,PetscInt *,PetscBool,PetscInt,PetscInt []);
PETSC_EXTERN PetscErrorCode indicesPointFields_private(PetscSection,PetscInt,PetscInt,PetscInt [],PetscBool,PetscInt,PetscInt []);
PETSC_INTERN PetscErrorCode DMPlexLocatePoint_Internal(DM,PetscInt,const PetscScalar [],PetscInt,PetscInt *);

PETSC_INTERN PetscErrorCode DMPlexCreateCellNumbering_Internal(DM, PetscBool, IS *);
PETSC_INTERN PetscErrorCode DMPlexCreateVertexNumbering_Internal(DM, PetscBool, IS *);
PETSC_INTERN PetscErrorCode DMPlexRefine_Internal(DM, DMLabel, DM *);
PETSC_INTERN PetscErrorCode DMPlexCoarsen_Internal(DM, DMLabel, DM *);

/* invert dihedral symmetry: return a^-1,
 * using the representation described in
 * DMPlexGetConeOrientation() */
PETSC_STATIC_INLINE PetscInt DihedralInvert(PetscInt N, PetscInt a)
{
  return (a <= 0) ? a : (N - a);
}

/* invert dihedral symmetry: return b * a,
 * using the representation described in
 * DMPlexGetConeOrientation() */
PETSC_STATIC_INLINE PetscInt DihedralCompose(PetscInt N, PetscInt a, PetscInt b)
{
  if (!N) return 0;
  return  (a >= 0) ?
         ((b >= 0) ? ((a + b) % N) : -(((a - b - 1) % N) + 1)) :
         ((b >= 0) ? -(((N - b - a - 1) % N) + 1) : ((N + b - a) % N));
}

/* swap dihedral symmetries: return b * a^-1,
 * using the representation described in
 * DMPlexGetConeOrientation() */
PETSC_STATIC_INLINE PetscInt DihedralSwap(PetscInt N, PetscInt a, PetscInt b)
{
  return DihedralCompose(N,DihedralInvert(N,a),b);
}

PETSC_EXTERN PetscErrorCode DMPlexComputeResidual_Internal(DM, IS , PetscReal, Vec, Vec, PetscReal, Vec, void *);
PETSC_EXTERN PetscErrorCode DMPlexComputeJacobian_Internal(DM, IS, PetscReal, PetscReal, Vec, Vec, Mat, Mat, void *);
PETSC_EXTERN PetscErrorCode DMPlexReconstructGradients_Internal(DM, PetscFV, PetscInt, PetscInt, Vec, Vec, Vec, Vec);

PETSC_STATIC_INLINE void DMPlex_Invert2D_Internal(PetscReal invJ[], PetscReal J[], PetscReal detJ)
{
  const PetscReal invDet = 1.0/detJ;

  invJ[0] =  invDet*J[3];
  invJ[1] = -invDet*J[1];
  invJ[2] = -invDet*J[2];
  invJ[3] =  invDet*J[0];
  (void)PetscLogFlops(5.0);
}

PETSC_STATIC_INLINE void DMPlex_Invert3D_Internal(PetscReal invJ[], PetscReal J[], PetscReal detJ)
{
  const PetscReal invDet = 1.0/detJ;

  invJ[0*3+0] = invDet*(J[1*3+1]*J[2*3+2] - J[1*3+2]*J[2*3+1]);
  invJ[0*3+1] = invDet*(J[0*3+2]*J[2*3+1] - J[0*3+1]*J[2*3+2]);
  invJ[0*3+2] = invDet*(J[0*3+1]*J[1*3+2] - J[0*3+2]*J[1*3+1]);
  invJ[1*3+0] = invDet*(J[1*3+2]*J[2*3+0] - J[1*3+0]*J[2*3+2]);
  invJ[1*3+1] = invDet*(J[0*3+0]*J[2*3+2] - J[0*3+2]*J[2*3+0]);
  invJ[1*3+2] = invDet*(J[0*3+2]*J[1*3+0] - J[0*3+0]*J[1*3+2]);
  invJ[2*3+0] = invDet*(J[1*3+0]*J[2*3+1] - J[1*3+1]*J[2*3+0]);
  invJ[2*3+1] = invDet*(J[0*3+1]*J[2*3+0] - J[0*3+0]*J[2*3+1]);
  invJ[2*3+2] = invDet*(J[0*3+0]*J[1*3+1] - J[0*3+1]*J[1*3+0]);
  (void)PetscLogFlops(37.0);
}

PETSC_STATIC_INLINE void DMPlex_Det2D_Internal(PetscReal *detJ, const PetscReal J[])
{
  *detJ = J[0]*J[3] - J[1]*J[2];
  (void)PetscLogFlops(3.0);
}

PETSC_STATIC_INLINE void DMPlex_Det3D_Internal(PetscReal *detJ, const PetscReal J[])
{
  *detJ = (J[0*3+0]*(J[1*3+1]*J[2*3+2] - J[1*3+2]*J[2*3+1]) +
           J[0*3+1]*(J[1*3+2]*J[2*3+0] - J[1*3+0]*J[2*3+2]) +
           J[0*3+2]*(J[1*3+0]*J[2*3+1] - J[1*3+1]*J[2*3+0]));
  (void)PetscLogFlops(12.0);
}

PETSC_STATIC_INLINE void DMPlex_Det2D_Scalar_Internal(PetscReal *detJ, const PetscScalar J[])
{
  *detJ = PetscRealPart(J[0])*PetscRealPart(J[3]) - PetscRealPart(J[1])*PetscRealPart(J[2]);
  (void)PetscLogFlops(3.0);
}

PETSC_STATIC_INLINE void DMPlex_Det3D_Scalar_Internal(PetscReal *detJ, const PetscScalar J[])
{
  *detJ = (PetscRealPart(J[0*3+0])*(PetscRealPart(J[1*3+1])*PetscRealPart(J[2*3+2]) - PetscRealPart(J[1*3+2])*PetscRealPart(J[2*3+1])) +
           PetscRealPart(J[0*3+1])*(PetscRealPart(J[1*3+2])*PetscRealPart(J[2*3+0]) - PetscRealPart(J[1*3+0])*PetscRealPart(J[2*3+2])) +
           PetscRealPart(J[0*3+2])*(PetscRealPart(J[1*3+0])*PetscRealPart(J[2*3+1]) - PetscRealPart(J[1*3+1])*PetscRealPart(J[2*3+0])));
  (void)PetscLogFlops(12.0);
}

PETSC_STATIC_INLINE void DMPlex_WaxpyD_Internal(PetscInt dim, PetscReal a, const PetscReal *x, const PetscReal *y, PetscReal *w) {PetscInt d; for (d = 0; d < dim; ++d) w[d] = a*x[d] + y[d];}

PETSC_STATIC_INLINE PetscReal DMPlex_DotD_Internal(PetscInt dim, const PetscScalar *x, const PetscReal *y) {PetscReal sum = 0.0; PetscInt d; for (d = 0; d < dim; ++d) sum += PetscRealPart(x[d])*y[d]; return sum;}

PETSC_STATIC_INLINE PetscReal DMPlex_DotRealD_Internal(PetscInt dim, const PetscReal *x, const PetscReal *y) {PetscReal sum = 0.0; PetscInt d; for (d = 0; d < dim; ++d) sum += x[d]*y[d]; return sum;}

PETSC_STATIC_INLINE PetscReal DMPlex_NormD_Internal(PetscInt dim, const PetscReal *x) {PetscReal sum = 0.0; PetscInt d; for (d = 0; d < dim; ++d) sum += x[d]*x[d]; return PetscSqrtReal(sum);}

PETSC_INTERN PetscErrorCode DMPlexGetPointDualSpaceFEM(DM,PetscInt,PetscInt,PetscDualSpace *);
PETSC_INTERN PetscErrorCode DMPlexGetIndicesPoint_Internal(PetscSection,PetscInt,PetscInt,PetscInt *,PetscBool,const PetscInt[],PetscInt[]);
PETSC_INTERN PetscErrorCode DMPlexGetIndicesPointFields_Internal(PetscSection,PetscInt,PetscInt,PetscInt[],PetscBool,const PetscInt***,PetscInt,PetscInt[]);

PETSC_EXTERN PetscErrorCode DMSNESGetFEGeom(DMField, IS, PetscQuadrature, PetscBool, PetscFEGeom **);
PETSC_EXTERN PetscErrorCode DMSNESRestoreFEGeom(DMField, IS, PetscQuadrature, PetscBool, PetscFEGeom **);
PETSC_EXTERN PetscErrorCode DMPlexComputeJacobian_Patch_Internal(DM, PetscSection, PetscSection, IS, PetscReal, PetscReal, Vec, Vec, Mat, Mat, void *);
PETSC_INTERN PetscErrorCode DMCreateSubDomainDM_Plex(DM,DMLabel,PetscInt,IS*,DM*);

#endif /* _PLEXIMPL_H */
