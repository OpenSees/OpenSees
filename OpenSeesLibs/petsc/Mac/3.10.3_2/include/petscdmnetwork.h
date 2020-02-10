/*
  DMNetwork, for parallel unstructured network problems.
*/
#if !defined(__PETSCDMNETWORK_H)
#define __PETSCDMNETWORK_H

#include <petscdm.h>
#include <petscviewer.h>

/*
  DMNetworkComponentGenericDataType - This is the data type that PETSc uses for storing the component data.
            For compatibility with PetscSF, which is used for data distribution, its declared as PetscInt.
	    To get the user-specific data type, one needs to cast it to the appropriate type.
*/
typedef PetscInt DMNetworkComponentGenericDataType;

PETSC_EXTERN PetscErrorCode DMNetworkCreate(MPI_Comm,DM*);
PETSC_EXTERN PetscErrorCode DMNetworkSetSizes(DM,PetscInt,PetscInt,PetscInt[],PetscInt[],PetscInt[],PetscInt[]);
PETSC_EXTERN PetscErrorCode DMNetworkSetEdgeList(DM,PetscInt*[],PetscInt*[]);
PETSC_EXTERN PetscErrorCode DMNetworkLayoutSetUp(DM);
PETSC_EXTERN PetscErrorCode DMNetworkRegisterComponent(DM,const char*,PetscInt,PetscInt*);
PETSC_EXTERN PetscErrorCode DMNetworkGetVertexRange(DM,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode DMNetworkGetEdgeRange(DM,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode DMNetworkAddComponent(DM,PetscInt,PetscInt,void*);
PETSC_EXTERN PetscErrorCode DMNetworkGetComponent(DM,PetscInt,PetscInt,PetscInt*,void**);
PETSC_EXTERN PetscErrorCode DMNetworkGetNumComponents(DM,PetscInt,PetscInt*);
PETSC_EXTERN PetscErrorCode DMNetworkGetVariableOffset(DM,PetscInt,PetscInt*);
PETSC_EXTERN PetscErrorCode DMNetworkGetVariableGlobalOffset(DM,PetscInt,PetscInt*);
PETSC_EXTERN PetscErrorCode DMNetworkGetEdgeOffset(DM,PetscInt,PetscInt*);
PETSC_EXTERN PetscErrorCode DMNetworkGetVertexOffset(DM,PetscInt,PetscInt*);
PETSC_EXTERN PetscErrorCode DMNetworkAddNumVariables(DM,PetscInt,PetscInt);
PETSC_EXTERN PetscErrorCode DMNetworkGetNumVariables(DM,PetscInt,PetscInt*);
PETSC_EXTERN PetscErrorCode DMNetworkSetNumVariables(DM,PetscInt,PetscInt);
PETSC_EXTERN PetscErrorCode DMNetworkAssembleGraphStructures(DM);
PETSC_EXTERN PetscErrorCode PetscSFGetSubSF(PetscSF,ISLocalToGlobalMapping,PetscSF*);
PETSC_EXTERN PetscErrorCode DMNetworkDistribute(DM*,PetscInt);
PETSC_EXTERN PetscErrorCode DMNetworkGetSupportingEdges(DM,PetscInt,PetscInt*,const PetscInt*[]);
PETSC_EXTERN PetscErrorCode DMNetworkGetConnectedVertices(DM,PetscInt,const PetscInt*[]);
PETSC_EXTERN PetscErrorCode DMNetworkIsGhostVertex(DM,PetscInt,PetscBool*);
PETSC_EXTERN PetscErrorCode DMNetworkEdgeSetMatrix(DM,PetscInt,Mat[]);
PETSC_EXTERN PetscErrorCode DMNetworkVertexSetMatrix(DM,PetscInt,Mat[]);
PETSC_EXTERN PetscErrorCode DMNetworkHasJacobian(DM,PetscBool,PetscBool);
PETSC_EXTERN PetscErrorCode DMNetworkGetPlex(DM,DM*);
PETSC_EXTERN PetscErrorCode DMNetworkGetGlobalEdgeIndex(DM,PetscInt,PetscInt*);
PETSC_EXTERN PetscErrorCode DMNetworkGetGlobalVertexIndex(DM,PetscInt,PetscInt*);

PETSC_EXTERN PetscErrorCode DMNetworkGetSubnetworkInfo(DM,PetscInt,PetscInt*,PetscInt*,const PetscInt**,const PetscInt**);

typedef struct _p_DMNetworkMonitorList *DMNetworkMonitorList;
struct _p_DMNetworkMonitorList
{
  PetscViewer viewer;
  Vec         v;
  PetscInt    element;
  PetscInt    nodes;
  PetscInt    start;
  PetscInt    blocksize;
  DMNetworkMonitorList next;
};

typedef struct _p_DMNetworkMonitor *DMNetworkMonitor;
struct _p_DMNetworkMonitor
{
  MPI_Comm             comm;
  DM                   network;
  DMNetworkMonitorList firstnode;
};

PETSC_EXTERN PetscErrorCode DMNetworkMonitorCreate(DM,DMNetworkMonitor*);
PETSC_EXTERN PetscErrorCode DMNetworkMonitorDestroy(DMNetworkMonitor*);
PETSC_EXTERN PetscErrorCode DMNetworkMonitorPop(DMNetworkMonitor);
PETSC_EXTERN PetscErrorCode DMNetworkMonitorAdd(DMNetworkMonitor,const char*,PetscInt,PetscInt,PetscInt,PetscInt,PetscReal,PetscReal,PetscBool);
PETSC_EXTERN PetscErrorCode DMNetworkMonitorView(DMNetworkMonitor,Vec);

#endif
