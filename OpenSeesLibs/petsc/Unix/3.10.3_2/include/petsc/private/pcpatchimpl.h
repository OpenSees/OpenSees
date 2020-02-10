/*
      Data structure used for Patch preconditioner.
*/
#if !defined(__PATCH_IMPL)
#define __PATCH_IMPL
#include <petsc/private/pcimpl.h>
#include <petsc/private/hashseti.h>
#include <petsc/private/hashmapi.h>
#include <petscksp.h>

typedef struct {
  /* Topology */
  PCPatchConstructType ctype;              /* Algorithm for patch construction */
  PetscErrorCode     (*patchconstructop)(void*, DM, PetscInt, PetscHSetI); /* patch construction */
  PetscErrorCode     (*userpatchconstructionop)(PC, PetscInt*, IS**, IS*, void* ctx);
  void                *userpatchconstructctx;
  IS                  *userIS;
  PetscInt             npatch;             /* Number of patches */
  PetscBool            user_patches;       /* Flag for user construction of patches */
  PetscInt             dim, codim;         /* Dimension or codimension of mesh points to loop over; only one of them can be set */
  PetscSection         cellCounts;         /* Maps patch -> # cells in patch */
  IS                   cells;              /* [patch][cell in patch]: Cell number */
  PetscSection         cellNumbering;      /* Plex: NULL Firedrake: Numbering of cells in DM */
  PetscSection         pointCounts;        /* Maps patch -> # points with dofs in patch */
  IS                   points;             /* [patch][point in patch]: Point number */
  /* Dof layout */
  PetscBool            combined;           /* Use a combined space with all fields */
  PetscInt             nsubspaces;         /* Number of fields */
  PetscSF              defaultSF;          /* Combined SF mapping process local to global */
  PetscSection        *dofSection;         /* ?? For each field, patch -> # dofs in patch */
  PetscInt            *subspaceOffsets;    /* Plex: NULL Firedrake: offset of each field in concatenated process local numbering for mixed spaces */
  PetscInt           **cellNodeMap;        /* [field][cell][dof in cell]: global dofs in cell TODO Free this after its use in PCPatchCreateCellPatchDiscretisationInfo() */
  IS                   dofs;               /* [patch][cell in patch][dof in cell]: patch local dof */
  IS                   offs;               /* [patch][point in patch]: patch local offset (same layout as 'points', used for filling up patchSection) */
  PetscSection         patchSection;       /* Maps points -> patch local dofs */
  IS                   globalBcNodes;      /* Global dofs constrained by global Dirichlet conditions TODO Replace these with process local constrained dofs */
  IS                   ghostBcNodes;       /* Global dofs constrained by global Dirichlet conditions on this process and possibly others (patch overlaps boundary) */
  PetscSection         gtolCounts;         /* ?? Indices to extract from local to patch vectors */
  IS                   gtol;
  PetscInt            *bs;                 /* [field] block size per field (can come from global operators?) */
  PetscInt            *nodesPerCell;       /* [field] Dofs per cell TODO Change "node" to "dof" everywhere */
  PetscInt             totalDofsPerCell;   /* Dofs per cell counting all fields */
  PetscInt             exclude_subspace;   /* If you don't want any other dofs from a particular subspace you can exclude them with this.
                                                Used for Vanka in Stokes, for example, to eliminate all pressure dofs not on the vertex
                                                you're building the patch around */
  PetscInt             vankadim;           /* In Vanka construction, should we eliminate any entities of a certain dimension on the initial patch? */
  PetscInt             ignoredim;          /* In Vanka construction, should we eliminate any entities of a certain dimension on the boundary? */
  /* Patch system assembly */
  PetscErrorCode     (*usercomputeop)(PC, PetscInt, Mat, IS, PetscInt, const PetscInt *, void *);
  IS                   cellIS;             /* Temporary IS for each cell patch */
  void                *usercomputectx;
  PetscBool            save_operators;     /* Save all operators (or create/destroy one at a time?) */
  PetscBool            partition_of_unity; /* Weight updates by dof multiplicity? */
  /* Patch solves */
  KSP                 *ksp;                /* Solvers for each patch TODO Do we need a new KSP for each patch? */
  Mat                 *mat;                /* System matrix for each patch */
  MatType              sub_mat_type;       /* Matrix type for patch systems */
  Vec                 *patchX, *patchY;    /* RHS and solution for each patch */
  Vec                 *patch_dof_weights;  /* Weighting for dof in each patch */
  Vec                  localX, localY;     /* ??? */
  Vec                  dof_weights;        /* In how many patches does each dof lie? */
  PetscBool            symmetrise_sweep;   /* Should we sweep forwards->backwards, backwards->forwards? */
  PetscBool            optionsSet;         /* SetFromOptions was called on this PC */
  IS                   iterationSet;       /* Index set specifying how we iterate over patches */
  /* Monitoring */
  PetscBool            viewPatches;        /* View information about patch construction */
  PetscBool            viewCells;          /* View cells for each patch */
  PetscViewer          viewerCells;        /*   Viewer for patch cells */
  PetscViewerFormat    formatCells;        /*   Format for patch cells */
  PetscBool            viewPoints;         /* View points for each patch */
  PetscViewer          viewerPoints;       /*   Viewer for patch points */
  PetscViewerFormat    formatPoints;       /*   Format for patch points */
  PetscBool            viewSection;        /* View global section for each patch */
  PetscViewer          viewerSection;      /*   Viewer for patch sections */
  PetscViewerFormat    formatSection;      /*   Format for patch sections */
  PetscBool            viewMatrix;         /* View matrix for each patch */
  PetscViewer          viewerMatrix;       /*   Viewer for patch matrix */
  PetscViewerFormat    formatMatrix;       /*   Format for patch matrix */
} PC_PATCH;

PETSC_EXTERN PetscLogEvent PC_Patch_CreatePatches;
PETSC_EXTERN PetscLogEvent PC_Patch_ComputeOp;
PETSC_EXTERN PetscLogEvent PC_Patch_Solve;
PETSC_EXTERN PetscLogEvent PC_Patch_Scatter;
PETSC_EXTERN PetscLogEvent PC_Patch_Apply;
PETSC_EXTERN PetscLogEvent PC_Patch_Prealloc;

#endif
