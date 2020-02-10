#if !defined(_PETSCDSIMPL_H)
#define _PETSCDSIMPL_H

#include <petscds.h>
#include <petsc/private/petscimpl.h>

PETSC_EXTERN PetscBool      PetscDSRegisterAllCalled;
PETSC_EXTERN PetscErrorCode PetscDSRegisterAll(void);

typedef struct _n_DSBoundary *DSBoundary;

struct _n_DSBoundary {
  const char *name;
  const char *labelname;
  DMBoundaryConditionType type;
  PetscInt    field;
  PetscInt    numcomps;
  PetscInt   *comps;
  void      (*func)(void);
  PetscInt    numids;
  PetscInt   *ids;
  void       *ctx;
  DSBoundary  next;
};

typedef struct _PetscDSOps *PetscDSOps;
struct _PetscDSOps {
  PetscErrorCode (*setfromoptions)(PetscDS);
  PetscErrorCode (*setup)(PetscDS);
  PetscErrorCode (*view)(PetscDS,PetscViewer);
  PetscErrorCode (*destroy)(PetscDS);
};

struct _p_PetscDS {
  PETSCHEADER(struct _PetscDSOps);
  void        *data;              /* Implementation object */
  PetscDS     *subprobs;          /* The subspaces for each dimension */
  PetscBool    setup;             /* Flag for setup */
  PetscInt     Nf;                /* The number of solution fields */
  PetscBool   *implicit;          /* Flag for implicit or explicit solve for each field */
  PetscBool    defaultAdj[2];     /* [use cone() or support() first, use the transitive closure] for the case of no fields */
  PetscBool   *adjacency;         /* Flags for defining variable influence (adjacency) for each field [use cone() or support() first, use the transitive closure] */
  PetscBool    useJacPre;         /* Flag for using the Jacobian preconditioner */
  PetscObject *disc;              /* The discretization for each solution field (PetscFE, PetscFV, etc.) */
  PetscPointFunc   *obj;          /* Scalar integral (like an objective function) */
  PetscPointFunc   *f;            /* Weak form integrands for F, f_0, f_1 */
  PetscPointJac    *g;            /* Weak form integrands for J = dF/du, g_0, g_1, g_2, g_3 */
  PetscPointJac    *gp;           /* Weak form integrands for preconditioner for J, g_0, g_1, g_2, g_3 */
  PetscPointJac    *gt;           /* Weak form integrands for dF/du_t, g_0, g_1, g_2, g_3 */
  PetscBdPointFunc *fBd;          /* Weak form boundary integrands F_bd, f_0, f_1 */
  PetscBdPointJac  *gBd;          /* Weak form boundary integrands J_bd = dF_bd/du, g_0, g_1, g_2, g_3 */
  PetscRiemannFunc *r;            /* Riemann solvers */
  PetscPointFunc   *update;       /* Direct update of field coefficients */
  PetscSimplePointFunc *exactSol; /* Exact solutions for each field */
  PetscInt          numConstants; /* Number of constants passed to point functions */
  PetscScalar      *constants;    /* Array of constants passed to point functions */
  void       **ctx;               /* User contexts for each field */
  PetscInt     dimEmbed;          /* The real space coordinate dimension */
  /* Computed sizes */
  PetscInt     totDim;            /* Total system dimension */
  PetscInt     totComp;           /* Total field components */
  PetscInt    *Nc;                /* Number of components for each field */
  PetscInt    *Nb;                /* Number of basis functions for each field */
  PetscInt    *off;               /* Offsets for each field */
  PetscInt    *offDer;            /* Derivative offsets for each field */
  PetscReal  **basis;             /* Default basis tabulation for each field */
  PetscReal  **basisDer;          /* Default basis derivative tabulation for each field */
  PetscReal  **basisFace;         /* Basis tabulation for each local face and field */
  PetscReal  **basisDerFace;      /* Basis derivative tabulation for each local face and field */
  /* Work space */
  PetscScalar *u;                 /* Field evaluation */
  PetscScalar *u_t;               /* Field time derivative evaluation */
  PetscScalar *u_x;               /* Field gradient evaluation */
  PetscScalar *refSpaceDer;       /* Workspace for computing derivative in the reference coordinates */
  PetscReal   *x;                 /* Workspace for computing real coordinates */
  PetscScalar *f0, *f1;           /* Point evaluations of weak form residual integrands */
  PetscScalar *g0, *g1, *g2, *g3; /* Point evaluations of weak form Jacobian integrands */
  DSBoundary   boundary;          /* Linked list of boundary conditions */
};

typedef struct {
  PetscInt dummy; /* */
} PetscDS_Basic;

#endif
