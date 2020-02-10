
#if !defined(_SWARMIMPL_H)
#define _SWARMIMPL_H

#include <petscvec.h> /*I "petscvec.h" I*/
#include <petscmat.h>       /*I      "petscmat.h"          I*/
#include <petscdmswarm.h> /*I      "petscdmswarm.h"    I*/
#include <petsc/private/dmimpl.h>

PETSC_EXTERN PetscLogEvent DMSWARM_Migrate;
PETSC_EXTERN PetscLogEvent DMSWARM_SetSizes;
PETSC_EXTERN PetscLogEvent DMSWARM_AddPoints;
PETSC_EXTERN PetscLogEvent DMSWARM_RemovePoints;
PETSC_EXTERN PetscLogEvent DMSWARM_Sort;
PETSC_EXTERN PetscLogEvent DMSWARM_DataExchangerTopologySetup;
PETSC_EXTERN PetscLogEvent DMSWARM_DataExchangerBegin;
PETSC_EXTERN PetscLogEvent DMSWARM_DataExchangerEnd;
PETSC_EXTERN PetscLogEvent DMSWARM_DataExchangerSendCount;
PETSC_EXTERN PetscLogEvent DMSWARM_DataExchangerPack;

typedef struct _p_DMSwarmDataField* DMSwarmDataField;
typedef struct _p_DMSwarmDataBucket* DMSwarmDataBucket;
typedef struct _p_DMSwarmSort* DMSwarmSort;

typedef struct {
  DMSwarmDataBucket db;

  PetscBool field_registration_initialized;
  PetscBool field_registration_finalized;
  /* DMSwarmProjectMethod *swarm_project;*/ /* swarm, geometry, result */

  /* PetscInt overlap; */
  /* PetscErrorCode (*update_overlap)(void); */

  char      vec_field_name[PETSC_MAX_PATH_LEN];
  PetscBool vec_field_set;
  PetscInt  vec_field_bs,vec_field_nlocal;

  PetscBool          issetup;
  DMSwarmType        swarm_type;
  DMSwarmMigrateType migrate_type;
  DMSwarmCollectType collect_type;

  DM        dmcell;

  PetscBool migrate_error_on_missing_point;

  PetscBool collect_view_active;
  PetscInt  collect_view_reset_nlocal;
  DMSwarmSort sort_context;
} DM_Swarm;

typedef struct {
  PetscInt point_index;
  PetscInt cell_index;
} SwarmPoint;

struct _p_DMSwarmSort {
  PetscBool isvalid;
  PetscInt ncells,npoints;
  PetscInt *pcell_offsets;
  SwarmPoint *list;
};


PETSC_INTERN PetscErrorCode DMSwarmMigrate_Push_Basic(DM, PetscBool);
PETSC_INTERN PetscErrorCode DMSwarmMigrate_CellDMScatter(DM,PetscBool);
PETSC_INTERN PetscErrorCode DMSwarmMigrate_CellDMExact(DM,PetscBool);

#endif /* _SWARMIMPL_H */
