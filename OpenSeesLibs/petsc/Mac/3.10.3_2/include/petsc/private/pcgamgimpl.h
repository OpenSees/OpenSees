#if !defined(__GAMG_IMPL)
#define __GAMG_IMPL
#include <petsc/private/pcimpl.h>
#include <petsc/private/pcmgimpl.h>                    /*I "petscksp.h" I*/
#include <petscmatcoarsen.h>                           /*I "petscmatcoarsen.h" I*/

struct _PCGAMGOps {
  PetscErrorCode (*graph)(PC, Mat, Mat*);
  PetscErrorCode (*coarsen)(PC, Mat*, PetscCoarsenData**);
  PetscErrorCode (*prolongator)(PC, Mat, Mat, PetscCoarsenData*, Mat*);
  PetscErrorCode (*optprolongator)(PC, Mat, Mat*);
  PetscErrorCode (*createlevel)(PC, Mat, PetscInt, Mat *, Mat *, PetscMPIInt *, IS *, PetscBool);
  PetscErrorCode (*createdefaultdata)(PC, Mat); /* for data methods that have a default (SA) */
  PetscErrorCode (*setfromoptions)(PetscOptionItems*,PC);
  PetscErrorCode (*destroy)(PC);
  PetscErrorCode (*view)(PC,PetscViewer);
};
#define PETSC_GAMG_MAXLEVELS 30
/* Private context for the GAMG preconditioner */
typedef struct gamg_TAG {
  PCGAMGType type;
  PetscInt  Nlevels;
  PetscInt  setup_count;
  PetscBool repart;
  PetscBool reuse_prol;
  PetscBool use_aggs_in_asm;
  PetscBool use_parallel_coarse_grid_solver;
  PetscInt  min_eq_proc;
  PetscInt  coarse_eq_limit;
  PetscReal threshold_scale;
  PetscInt  current_level; /* stash construction state */
  PetscReal threshold[PETSC_GAMG_MAXLEVELS]; /* common quatity to many AMG methods so keep it up here */

  /* these 4 are all related to the method data and should be in the subctx */
  PetscInt  data_sz;      /* nloc*data_rows*data_cols */
  PetscInt  data_cell_rows;
  PetscInt  data_cell_cols;
  PetscInt  orig_data_cell_rows;
  PetscInt  orig_data_cell_cols;
  PetscReal *data;          /* [data_sz] blocked vector of vertex data on fine grid (coordinates/nullspace) */
  PetscReal *orig_data;          /* cache data */

  struct _PCGAMGOps *ops;
  char *gamg_type_name;

  void *subctx;
} PC_GAMG;

PetscErrorCode PCReset_MG(PC);

/* hooks create derivied classes */
PetscErrorCode  PCCreateGAMG_GEO(PC);
PetscErrorCode  PCCreateGAMG_AGG(PC);
PetscErrorCode  PCCreateGAMG_Classical(PC);

PetscErrorCode PCDestroy_GAMG(PC);

/* helper methods */
PetscErrorCode PCGAMGCreateGraph(Mat, Mat*);
PetscErrorCode PCGAMGFilterGraph(Mat*, PetscReal, PetscBool);
PetscErrorCode PCGAMGGetDataWithGhosts(Mat, PetscInt, PetscReal[],PetscInt*, PetscReal **);

#if defined PETSC_USE_LOG
#define PETSC_GAMG_USE_LOG
enum tag {SET1,SET2,GRAPH,GRAPH_MAT,GRAPH_FILTER,GRAPH_SQR,SET4,SET5,SET6,FIND_V,SET7,SET8,SET9,SET10,SET11,SET12,SET13,SET14,SET15,SET16,NUM_SET};
#if defined PETSC_GAMG_USE_LOG
PETSC_EXTERN PetscLogEvent petsc_gamg_setup_events[NUM_SET];
#endif
PETSC_EXTERN PetscLogEvent PC_GAMGGraph_AGG;
PETSC_EXTERN PetscLogEvent PC_GAMGGraph_GEO;
PETSC_EXTERN PetscLogEvent PC_GAMGCoarsen_AGG;
PETSC_EXTERN PetscLogEvent PC_GAMGCoarsen_GEO;
PETSC_EXTERN PetscLogEvent PC_GAMGProlongator_AGG;
PETSC_EXTERN PetscLogEvent PC_GAMGProlongator_GEO;
PETSC_EXTERN PetscLogEvent PC_GAMGOptProlongator_AGG;
#endif

typedef struct _PCGAMGHashTable {
  PetscInt *table;
  PetscInt *data;
  PetscInt size;
} PCGAMGHashTable;


PETSC_INTERN PetscErrorCode PCGAMGHashTableCreate(PetscInt, PCGAMGHashTable*);
PETSC_INTERN PetscErrorCode PCGAMGHashTableDestroy(PCGAMGHashTable*);
PETSC_INTERN PetscErrorCode PCGAMGHashTableAdd(PCGAMGHashTable*,PetscInt,PetscInt);

#define GAMG_HASH(key) (PetscInt)((((PetscInt64)7)*(PetscInt64)key)%(PetscInt64)a_tab->size)
PETSC_STATIC_INLINE PetscErrorCode PCGAMGHashTableFind(PCGAMGHashTable *a_tab, PetscInt a_key, PetscInt *a_data)
{
  PetscInt kk,idx;

  PetscFunctionBegin;
  if (a_key<0) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"Negative key %D.",a_key);
  for (kk = 0, idx = GAMG_HASH(a_key); kk < a_tab->size; kk++, idx = (idx==(a_tab->size-1)) ? 0 : idx + 1) {
    if (a_tab->table[idx] == a_key) {
      *a_data = a_tab->data[idx];
      break;
    } else if (a_tab->table[idx] == -1) {
      /* not here */
      *a_data = -1;
      break;
    }
  }
  if (kk==a_tab->size) SETERRQ1(PETSC_COMM_SELF,PETSC_ERR_USER,"key %D not found in table",a_key);
  PetscFunctionReturn(0);
}

#endif

