#if !defined(__PETSCDMSWARM_H)
#define __PETSCDMSWARM_H

#include <petscdm.h>

/*E
   DMSwarmType - Defines the type of swarm
 
   DMSWARM_BASIC defines N entries of varied data-types which the user may register.

   DMSWARM_PIC is suitable for particle-in-cell methods. Configured as DMSWARM_PIC, the swarm will be aware of, another DM which serves as the background mesh. Fields specific to particle-in-cell methods are registered by default. These include spatial coordinates, a unique identifier, a cell index and an index for the owning rank. The background mesh will (by default) define the spatial decomposition of the points defined in the swarm. DMSWARM_PIC provides support for particle-in-cell operations such as defining initial point coordinates, communicating particles between sub-domains, projecting particle data fields on to the mesh.

   Level: beginner

.seealso: DMSwarmSetType()
E*/
typedef enum {
  DMSWARM_BASIC=0,
  DMSWARM_PIC
} DMSwarmType;

typedef enum {
  DMSWARM_MIGRATE_BASIC=0,
  DMSWARM_MIGRATE_DMCELLNSCATTER,
  DMSWARM_MIGRATE_DMCELLEXACT,
  DMSWARM_MIGRATE_USER
} DMSwarmMigrateType;

typedef enum {
  DMSWARM_COLLECT_BASIC=0,
  DMSWARM_COLLECT_DMDABOUNDINGBOX,
  DMSWARM_COLLECT_GENERAL,
  DMSWARM_COLLECT_USER
} DMSwarmCollectType;

/*E
   DMSwarmPICLayoutType - Defines the method used to define particle coordinates within each cell. The layouts are constructured using the reference cell geometry
 
   DMSWARMPIC_LAYOUT_REGULAR defines points on a regular ijk mesh.
   When using DMSWARMPIC_LAYOUT_REGULAR, the fill_param defines the number of points in each spatial direction.
 
   DMSWARMPIC_LAYOUT_GAUSS defines points using an npoint Gauss-Legendre tensor product quadrature rule.
   When using DMSWARMPIC_LAYOUT_GAUSS, the fill_param defines the number of quadrature points in each spatial direction.
 
   DMSWARMPIC_LAYOUT_SUBDIVISION defines points on the centroid of a sub-divided reference cell.
   When using DMSWARMPIC_LAYOUT_SUBDIVISION, the fill_param defines the number times the reference cell is sub-divided.
 
   Level: beginner
 
.seealso DMSwarmInsertPointsUsingCellDM()
E*/
typedef enum {
  DMSWARMPIC_LAYOUT_REGULAR=0,
  DMSWARMPIC_LAYOUT_GAUSS,
  DMSWARMPIC_LAYOUT_SUBDIVISION
} DMSwarmPICLayoutType;

PETSC_EXTERN const char* DMSwarmTypeNames[];
PETSC_EXTERN const char* DMSwarmMigrateTypeNames[];
PETSC_EXTERN const char* DMSwarmCollectTypeNames[];

PETSC_EXTERN const char DMSwarmField_pid[];
PETSC_EXTERN const char DMSwarmField_rank[];
PETSC_EXTERN const char DMSwarmPICField_coor[];
PETSC_EXTERN const char DMSwarmPICField_cellid[];

PETSC_EXTERN PetscErrorCode DMSwarmCreateGlobalVectorFromField(DM,const char[],Vec*);
PETSC_EXTERN PetscErrorCode DMSwarmDestroyGlobalVectorFromField(DM,const char[],Vec*);
PETSC_EXTERN PetscErrorCode DMSwarmCreateLocalVectorFromField(DM,const char[],Vec*);
PETSC_EXTERN PetscErrorCode DMSwarmDestroyLocalVectorFromField(DM,const char[],Vec*);

PETSC_EXTERN PetscErrorCode DMSwarmInitializeFieldRegister(DM);
PETSC_EXTERN PetscErrorCode DMSwarmFinalizeFieldRegister(DM);
PETSC_EXTERN PetscErrorCode DMSwarmSetLocalSizes(DM,PetscInt,PetscInt);
PETSC_EXTERN PetscErrorCode DMSwarmRegisterPetscDatatypeField(DM,const char[],PetscInt,PetscDataType);
PETSC_EXTERN PetscErrorCode DMSwarmRegisterUserStructField(DM,const char[],size_t);
PETSC_EXTERN PetscErrorCode DMSwarmRegisterUserDatatypeField(DM,const char[],size_t,PetscInt);
PETSC_EXTERN PetscErrorCode DMSwarmGetField(DM,const char[],PetscInt*,PetscDataType*,void**);
PETSC_EXTERN PetscErrorCode DMSwarmRestoreField(DM,const char[],PetscInt*,PetscDataType*,void**);

PETSC_EXTERN PetscErrorCode DMSwarmVectorDefineField(DM,const char[]);

PETSC_EXTERN PetscErrorCode DMSwarmAddPoint(DM);
PETSC_EXTERN PetscErrorCode DMSwarmAddNPoints(DM,PetscInt);
PETSC_EXTERN PetscErrorCode DMSwarmRemovePoint(DM);
PETSC_EXTERN PetscErrorCode DMSwarmRemovePointAtIndex(DM,PetscInt);
PETSC_EXTERN PetscErrorCode DMSwarmCopyPoint(DM dm,PetscInt,PetscInt);

PETSC_EXTERN PetscErrorCode DMSwarmGetLocalSize(DM,PetscInt*);
PETSC_EXTERN PetscErrorCode DMSwarmGetSize(DM,PetscInt*);
PETSC_EXTERN PetscErrorCode DMSwarmMigrate(DM,PetscBool);

PETSC_EXTERN PetscErrorCode DMSwarmCollectViewCreate(DM);
PETSC_EXTERN PetscErrorCode DMSwarmCollectViewDestroy(DM);
PETSC_EXTERN PetscErrorCode DMSwarmSetCellDM(DM,DM);
PETSC_EXTERN PetscErrorCode DMSwarmGetCellDM(DM,DM*);

PETSC_EXTERN PetscErrorCode DMSwarmSetType(DM,DMSwarmType);

PETSC_EXTERN PetscErrorCode DMSwarmSetPointsUniformCoordinates(DM,PetscReal*,PetscReal*,PetscInt*,InsertMode);
PETSC_EXTERN PetscErrorCode DMSwarmSetPointCoordinates(DM,PetscInt,PetscReal*,PetscBool,InsertMode);
PETSC_EXTERN PetscErrorCode DMSwarmInsertPointsUsingCellDM(DM,DMSwarmPICLayoutType,PetscInt);
PETSC_EXTERN PetscErrorCode DMSwarmSetPointCoordinatesCellwise(DM,PetscInt,PetscReal*);
PETSC_EXTERN PetscErrorCode DMSwarmViewFieldsXDMF(DM,const char*,PetscInt,const char**);
PETSC_EXTERN PetscErrorCode DMSwarmViewXDMF(DM,const char*);

PETSC_EXTERN PetscErrorCode DMSwarmSortGetAccess(DM);
PETSC_EXTERN PetscErrorCode DMSwarmSortRestoreAccess(DM);
PETSC_EXTERN PetscErrorCode DMSwarmSortGetPointsPerCell(DM,PetscInt,PetscInt*,PetscInt**);
PETSC_EXTERN PetscErrorCode DMSwarmSortGetNumberOfPointsPerCell(DM,PetscInt,PetscInt*);
PETSC_EXTERN PetscErrorCode DMSwarmSortGetIsValid(DM,PetscBool*);
PETSC_EXTERN PetscErrorCode DMSwarmSortGetSizes(DM,PetscInt*,PetscInt*);

PETSC_EXTERN PetscErrorCode DMSwarmProjectFields(DM,PetscInt,const char**,Vec**,PetscBool);

#endif

