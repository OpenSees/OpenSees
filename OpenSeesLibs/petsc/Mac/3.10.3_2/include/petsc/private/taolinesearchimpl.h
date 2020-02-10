#ifndef __TAOLINESEARCH_IMPL_H
#define __TAOLINESEARCH_IMPL_H
#include <petscvec.h>
#include <petsc/private/petscimpl.h>
#include <petsctaolinesearch.h>

typedef struct _TaoLineSearchOps *TaoLineSearchOps;
struct _TaoLineSearchOps {
    PetscErrorCode (*computeobjective)(TaoLineSearch, Vec, PetscReal*, void*);
    PetscErrorCode (*computegradient)(TaoLineSearch, Vec, Vec, void*);
    PetscErrorCode (*computeobjectiveandgradient)(TaoLineSearch, Vec, PetscReal *, Vec, void*);
    PetscErrorCode (*computeobjectiveandgts)(TaoLineSearch, Vec, Vec, PetscReal*, PetscReal*,void*);
    PetscErrorCode (*setup)(TaoLineSearch);
    PetscErrorCode (*apply)(TaoLineSearch,Vec,PetscReal*,Vec,Vec);
    PetscErrorCode (*view)(TaoLineSearch,PetscViewer);
    PetscErrorCode (*setfromoptions)(PetscOptionItems*,TaoLineSearch);
    PetscErrorCode (*reset)(TaoLineSearch);
    PetscErrorCode (*destroy)(TaoLineSearch);
    PetscErrorCode (*monitor)(TaoLineSearch);
};

struct _p_TaoLineSearch {
    PETSCHEADER(struct _TaoLineSearchOps);
    void *userctx_func;
    void *userctx_grad;
    void *userctx_funcgrad;
    void *userctx_funcgts;
    PetscBool usemonitor;
    PetscViewer viewer;
    
    PetscBool setupcalled;
    PetscBool usegts;
    PetscBool usetaoroutines;
    PetscBool hasobjective;
    PetscBool hasgradient;
    PetscBool hasobjectiveandgradient;
    void *data;

    /* bounds used for some line searches */
    Vec lower;
    Vec upper;
    PetscInt bounded;

    Vec start_x;
    Vec stepdirection;
    PetscReal f_fullstep;
    PetscReal new_f;
    Vec new_x;
    Vec new_g;

    PetscReal step;
    PetscReal initstep;

    PetscInt max_funcs;
    PetscInt nfeval;
    PetscInt ngeval;
    PetscInt nfgeval;
    TaoLineSearchConvergedReason reason;

    PetscReal rtol;      /* relative tol for acceptable step (rtol>0) */
    PetscReal ftol;      /* tol for sufficient decr. condition (ftol>0) */
    PetscReal gtol;      /* tol for curvature condition (gtol>0)*/
    PetscReal stepmin;   /* lower bound for step */
    PetscReal stepmax;   /* upper bound for step */

    Tao tao;
};

PETSC_EXTERN PetscLogEvent TAOLINESEARCH_Apply;
PETSC_EXTERN PetscLogEvent TAOLINESEARCH_Eval;
#endif
