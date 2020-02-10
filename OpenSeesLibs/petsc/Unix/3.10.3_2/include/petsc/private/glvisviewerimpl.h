#ifndef _GLVISIMPL_H
#define _GLVISIMPL_H

#include <petscviewer.h>
#include <petscsys.h>

struct _n_PetscViewerGLVisVecInfo {
  char* fec_type; /* the output of FiniteElementCollection::Name() */
};
typedef struct _n_PetscViewerGLVisVecInfo *PetscViewerGLVisVecInfo;

struct _n_PetscViewerGLVisInfo {
  PetscBool enabled;  /* whether or not to visualize data from the process (it works, but it currently misses a public API) */
  PetscBool init;     /* whether or not the popup window has been initialized (must be done after having sent the data the first time) */
  PetscReal pause;    /* pause argument */
  char*     fmt;      /* format */
};
typedef struct _n_PetscViewerGLVisInfo *PetscViewerGLVisInfo;

typedef enum {PETSCVIEWERGLVIS_DISCONNECTED, PETSCVIEWERGLVIS_CONNECTED, PETSCVIEWERGLVIS_DISABLED} PetscViewerGLVisStatus;

PETSC_EXTERN PetscErrorCode PetscViewerGLVisPause_Private(PetscViewer);
PETSC_EXTERN PetscErrorCode PetscViewerGLVisSetDM_Private(PetscViewer,PetscObject);
PETSC_EXTERN PetscErrorCode PetscViewerGLVisGetDM_Private(PetscViewer,PetscObject*);
PETSC_EXTERN PetscErrorCode PetscViewerGLVisInitWindow_Private(PetscViewer,PetscBool,PetscInt,const char*);
PETSC_EXTERN PetscErrorCode PetscViewerGLVisGetStatus_Private(PetscViewer,PetscViewerGLVisStatus*);
PETSC_EXTERN PetscErrorCode PetscViewerGLVisGetType_Private(PetscViewer,PetscViewerGLVisType*);
PETSC_EXTERN PetscErrorCode PetscViewerGLVisGetWindow_Private(PetscViewer,PetscInt,PetscViewer*);
PETSC_EXTERN PetscErrorCode PetscViewerGLVisRestoreWindow_Private(PetscViewer,PetscInt,PetscViewer*);
PETSC_EXTERN PetscErrorCode PetscViewerGLVisGetFields_Private(PetscViewer,PetscInt*,const char**[],PetscInt*[],PetscErrorCode(**)(PetscObject,PetscInt,PetscObject[],void*),PetscObject*[],void**);
PETSC_EXTERN PetscErrorCode PetscViewerGLVisGetDMWindow_Private(PetscViewer,PetscViewer*);
#endif
