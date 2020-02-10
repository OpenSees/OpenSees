
#ifndef __CHARACTERISTICIMPL_H
#define __CHARACTERISTICIMPL_H

#include <petsccharacteristic.h>
#include <petsc/private/petscimpl.h>

/* Logging support */
PETSC_EXTERN PetscClassId CHARACTERISTIC_CLASSID;
PETSC_EXTERN PetscBool        CharacteristicRegisterAllCalled;
PETSC_EXTERN PetscErrorCode   CharacteristicRegisterAll(void);
PETSC_EXTERN PetscLogEvent CHARACTERISTIC_SetUp;
PETSC_EXTERN PetscLogEvent CHARACTERISTIC_Solve;
PETSC_EXTERN PetscLogEvent CHARACTERISTIC_QueueSetup;
PETSC_EXTERN PetscLogEvent CHARACTERISTIC_DAUpdate;
PETSC_EXTERN PetscLogEvent CHARACTERISTIC_HalfTimeLocal;
PETSC_EXTERN PetscLogEvent CHARACTERISTIC_HalfTimeRemote;
PETSC_EXTERN PetscLogEvent CHARACTERISTIC_HalfTimeExchange;
PETSC_EXTERN PetscLogEvent CHARACTERISTIC_FullTimeLocal;
PETSC_EXTERN PetscLogEvent CHARACTERISTIC_FullTimeRemote;
PETSC_EXTERN PetscLogEvent CHARACTERISTIC_FullTimeExchange;

#define MAX_COMPONENTS 10

typedef struct _p_Item {
  int           proc; /* Relative processor from which data is required (mapped to absolute by neighbors) */
  int           i, j; /* The vertex for which we need field values */
  PetscScalar   x, y; /* Coordinates of a point on the characteristic */
  PetscScalar   u, v; /* Velocity of a point on the characteristic */
  PetscScalar   field[MAX_COMPONENTS]; /* Field being advected */
} CharacteristicPointDA2D;

typedef CharacteristicPointDA2D *Queue;

struct _CharacteristicOps {
  PetscErrorCode (*view)(Characteristic, PetscViewer);
  PetscErrorCode (*destroy)(Characteristic);
  PetscErrorCode (*setup)(Characteristic);
};

struct _p_Characteristic {
  PETSCHEADER(struct _CharacteristicOps);
  PetscInt     setupcalled;
  PetscBool    structured;      /* Flag for mesh type */
  PetscInt     numIds;          /* Number of integers necessary to identify a mesh element */
  /* Velocity interpolation structures */
  DM           velocityDA;      /* DM for the velocity field */
  Vec          velocity;        /* Velocity field at t_n */
  Vec          velocityOld;     /* Velocity field at t_n-1 */
  PetscInt     numVelocityComp; /* Number of velocity components (should be the mesh dimension) */
  PetscInt    *velocityComp;    /* Components of the velocity in the DM */
  PetscErrorCode (*velocityInterp)(Vec, PetscReal [], PetscInt, PetscInt [], PetscScalar [], void *);
  PetscErrorCode (*velocityInterpLocal)(void *, PetscReal [], PetscInt, PetscInt [], PetscScalar [], void *);
  void        *velocityCtx;     /* User context for velocity inteprolation */
  /* Field interpolation structures */
  DM           fieldDA;         /* DM for the field field */
  Vec          field;           /* Field field at t_n */
  Vec          fieldOld;        /* Field field at t_n-1 */
  PetscInt     numFieldComp;    /* Number of field components (should be the mesh dimension) */
  PetscInt    *fieldComp;       /* Components of the field in the DM */
  PetscErrorCode (*fieldInterp)(Vec, PetscReal [], PetscInt, PetscInt [], PetscScalar [], void *);
  PetscErrorCode (*fieldInterpLocal)(void *, PetscReal [], PetscInt, PetscInt [], PetscScalar [], void *);
  void        *fieldCtx;        /* User context for field inteprolation */
  /* Communication structures*/
  MPI_Datatype itemType;        /* Type corresponding to the item struct */
  Queue        queue;
  PetscInt     queueSize;
  PetscInt     queueMax;
  Queue        queueLocal;      /* Queue of Items to receive from other processes */
  PetscInt     queueLocalSize;
  PetscInt     queueLocalMax;
  Queue        queueRemote;     /* Queue of Items to send to other processes */
  PetscInt     queueRemoteSize;
  PetscInt     queueRemoteMax;
  PetscInt     numNeighbors;    /* Number of neighboring processes */
  PetscMPIInt *neighbors;       /* Ranks of neighbors */
  PetscInt    *needCount;       /* Number of Items requested from other processes */
  PetscInt    *localOffsets;    /* Offset into queue for each process (Prefix sums of need_count) */
  PetscInt    *fillCount;       /* Number of Items requested by other processes */
  PetscInt    *remoteOffsets;   /* Offset into remote queue for each process (Prefix sums of fill_count) */
  MPI_Request *request;         /* Requests for sizes/velocities/fields from other processes */
  MPI_Status  *status;          /* Status structues for the persistent requests */
  void        *data;            /* Holder for implementation class */
};

PETSC_EXTERN PetscErrorCode CharacteristicSetNeighbors(Characteristic, PetscInt, PetscMPIInt []);
PETSC_EXTERN PetscErrorCode CharacteristicAddPoint(Characteristic, CharacteristicPointDA2D *);
PETSC_EXTERN PetscErrorCode CharacteristicSendCoordinatesBegin(Characteristic);
PETSC_EXTERN PetscErrorCode CharacteristicSendCoordinatesEnd(Characteristic);
PETSC_EXTERN PetscErrorCode CharacteristicGetValuesBegin(Characteristic);
PETSC_EXTERN PetscErrorCode CharacteristicGetValuesEnd(Characteristic);

#endif /*__CHARACTERISTICIMPL_H*/
