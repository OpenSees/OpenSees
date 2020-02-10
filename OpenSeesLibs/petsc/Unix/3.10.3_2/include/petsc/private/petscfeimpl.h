#if !defined(_PETSCFEIMPL_H)
#define _PETSCFEIMPL_H

#include <petscfe.h>
#include <petscds.h>
#include <petsc/private/petscimpl.h>
#include <petsc/private/dmpleximpl.h>

PETSC_EXTERN PetscBool PetscSpaceRegisterAllCalled;
PETSC_EXTERN PetscBool PetscDualSpaceRegisterAllCalled;
PETSC_EXTERN PetscBool PetscFERegisterAllCalled;
PETSC_EXTERN PetscErrorCode PetscSpaceRegisterAll(void);
PETSC_EXTERN PetscErrorCode PetscDualSpaceRegisterAll(void);
PETSC_EXTERN PetscErrorCode PetscFERegisterAll(void);

PETSC_EXTERN PetscBool FEcite;
PETSC_EXTERN const char FECitation[];

typedef struct _PetscSpaceOps *PetscSpaceOps;
struct _PetscSpaceOps {
  PetscErrorCode (*setfromoptions)(PetscOptionItems*,PetscSpace);
  PetscErrorCode (*setup)(PetscSpace);
  PetscErrorCode (*view)(PetscSpace,PetscViewer);
  PetscErrorCode (*destroy)(PetscSpace);

  PetscErrorCode (*getdimension)(PetscSpace,PetscInt*);
  PetscErrorCode (*evaluate)(PetscSpace,PetscInt,const PetscReal*,PetscReal*,PetscReal*,PetscReal*);
  PetscErrorCode (*getheightsubspace)(PetscSpace,PetscInt,PetscSpace *);
};

struct _p_PetscSpace {
  PETSCHEADER(struct _PetscSpaceOps);
  void                   *data;          /* Implementation object */
  PetscInt                degree;        /* The approximation order of the space */
  PetscInt                maxDegree;     /* The containing approximation order of the space */
  PetscInt                Nc;            /* The number of components */
  PetscInt                Nv;            /* The number of variables in the space, e.g. x and y */
  PetscInt                dim;           /* The dimension of the space */
  DM                      dm;            /* Shell to use for temp allocation */
};

typedef struct {
  PetscBool   symmetric;    /* Use only symmetric polynomials */
  PetscBool   tensor;       /* Flag for tensor product */
  PetscInt   *degrees;      /* Degrees of single variable which we need to compute */
  PetscBool   setupCalled;
  PetscSpace *subspaces;    /* Subspaces for each dimension */
} PetscSpace_Poly;

typedef struct {
  PetscSpace *tensspaces;
  PetscInt    numTensSpaces;
  PetscInt    dim;
  PetscBool   uniform;
  PetscBool   setupCalled;
  PetscSpace *heightsubspaces;    /* Height subspaces */
} PetscSpace_Tensor;

typedef struct {
  PetscQuadrature quad;         /* The points defining the space */
} PetscSpace_Point;

typedef struct _PetscDualSpaceOps *PetscDualSpaceOps;
struct _PetscDualSpaceOps {
  PetscErrorCode (*setfromoptions)(PetscOptionItems*,PetscDualSpace);
  PetscErrorCode (*setup)(PetscDualSpace);
  PetscErrorCode (*view)(PetscDualSpace,PetscViewer);
  PetscErrorCode (*destroy)(PetscDualSpace);

  PetscErrorCode (*duplicate)(PetscDualSpace,PetscDualSpace*);
  PetscErrorCode (*getdimension)(PetscDualSpace,PetscInt*);
  PetscErrorCode (*getnumdof)(PetscDualSpace,const PetscInt**);
  PetscErrorCode (*getheightsubspace)(PetscDualSpace,PetscInt,PetscDualSpace *);
  PetscErrorCode (*getpointsubspace)(PetscDualSpace,PetscInt,PetscDualSpace *);
  PetscErrorCode (*getsymmetries)(PetscDualSpace,const PetscInt****,const PetscScalar****);
  PetscErrorCode (*apply)(PetscDualSpace, PetscInt, PetscReal, PetscFEGeom *, PetscInt, PetscErrorCode (*)(PetscInt, PetscReal, const PetscReal [], PetscInt, PetscScalar *, void *), void *, PetscScalar *);
  PetscErrorCode (*applyall)(PetscDualSpace, const PetscScalar *, PetscScalar *);
  PetscErrorCode (*createallpoints)(PetscDualSpace, PetscQuadrature *);
};

struct _p_PetscDualSpace {
  PETSCHEADER(struct _PetscDualSpaceOps);
  void            *data;       /* Implementation object */
  DM               dm;         /* The integration region K */
  PetscInt         order;      /* The approximation order of the space */
  PetscInt         Nc;         /* The number of components */
  PetscQuadrature *functional; /* The basis of functionals for this space */
  PetscQuadrature  allPoints;
  PetscBool        setupcalled;
};

typedef struct {
  PetscInt       *numDof;
  PetscBool       simplexCell;
  PetscBool       tensorSpace;
  PetscBool       continuous;
  PetscInt        height;
  PetscDualSpace *subspaces;
  PetscInt     ***symmetries;
  PetscInt        numSelfSym;
  PetscInt        selfSymOff;
} PetscDualSpace_Lag;

typedef struct {
  PetscInt  dim;
  PetscInt *numDof;
} PetscDualSpace_Simple;

typedef struct _PetscFEOps *PetscFEOps;
struct _PetscFEOps {
  PetscErrorCode (*setfromoptions)(PetscOptionItems*,PetscFE);
  PetscErrorCode (*setup)(PetscFE);
  PetscErrorCode (*view)(PetscFE,PetscViewer);
  PetscErrorCode (*destroy)(PetscFE);
  PetscErrorCode (*getdimension)(PetscFE,PetscInt*);
  PetscErrorCode (*gettabulation)(PetscFE,PetscInt,const PetscReal*,PetscReal*,PetscReal*,PetscReal*);
  /* Element integration */
  PetscErrorCode (*integrate)(PetscFE, PetscDS, PetscInt, PetscInt, PetscFEGeom *, const PetscScalar[], PetscDS, const PetscScalar[], PetscScalar[]);
  PetscErrorCode (*integratebd)(PetscFE, PetscDS, PetscInt, PetscBdPointFunc, PetscInt, PetscFEGeom *, const PetscScalar[], PetscDS, const PetscScalar[], PetscScalar[]);
  PetscErrorCode (*integrateresidual)(PetscFE, PetscDS, PetscInt, PetscInt, PetscFEGeom *, const PetscScalar[], const PetscScalar[], PetscDS, const PetscScalar[], PetscReal, PetscScalar[]);
  PetscErrorCode (*integratebdresidual)(PetscFE, PetscDS, PetscInt, PetscInt, PetscFEGeom *, const PetscScalar[], const PetscScalar[], PetscDS, const PetscScalar[], PetscReal, PetscScalar[]);
  PetscErrorCode (*integratejacobianaction)(PetscFE, PetscDS, PetscInt, PetscInt, PetscFEGeom *, const PetscScalar[], const PetscScalar[], PetscDS, const PetscScalar[], PetscReal, PetscReal, PetscScalar[]);
  PetscErrorCode (*integratejacobian)(PetscFE, PetscDS, PetscFEJacobianType, PetscInt, PetscInt, PetscInt, PetscFEGeom *, const PetscScalar[], const PetscScalar[], PetscDS, const PetscScalar[], PetscReal, PetscReal, PetscScalar[]);
  PetscErrorCode (*integratebdjacobian)(PetscFE, PetscDS, PetscInt, PetscInt, PetscInt, PetscFEGeom *, const PetscScalar[], const PetscScalar[], PetscDS, const PetscScalar[], PetscReal, PetscReal, PetscScalar[]);
};

struct _p_PetscFE {
  PETSCHEADER(struct _PetscFEOps);
  void           *data;                  /* Implementation object */
  PetscSpace      basisSpace;            /* The basis space P */
  PetscDualSpace  dualSpace;             /* The dual space P' */
  PetscInt        numComponents;         /* The number of field components */
  PetscQuadrature quadrature;            /* Suitable quadrature on K */
  PetscQuadrature faceQuadrature;        /* Suitable face quadrature on \partial K */
  PetscFE        *subspaces;             /* Subspaces for each dimension */
  PetscReal      *invV;                  /* Change of basis matrix, from prime to nodal basis set */
  PetscReal      *B,  *D,  *H;           /* Tabulation of basis and derivatives at quadrature points */
  PetscReal      *Bf, *Df, *Hf;          /* Tabulation of basis and derivatives at quadrature points on each face */
  PetscReal      *F;                     /* Tabulation of basis at face centroids */
  PetscInt        blockSize, numBlocks;  /* Blocks are processed concurrently */
  PetscInt        batchSize, numBatches; /* A batch is made up of blocks, Batches are processed in serial */
  PetscBool       setupcalled;
};

typedef struct {
  PetscInt cellType;
} PetscFE_Basic;

#ifdef PETSC_HAVE_OPENCL

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

typedef struct {
  cl_platform_id   pf_id;
  cl_device_id     dev_id;
  cl_context       ctx_id;
  cl_command_queue queue_id;
  PetscDataType    realType;
  PetscLogEvent    residualEvent;
  PetscInt         op; /* ANDY: Stand-in for real equation code generation */
} PetscFE_OpenCL;
#endif

typedef struct {
  CellRefiner   cellRefiner;    /* The cell refiner defining the cell division */
  PetscInt      numSubelements; /* The number of subelements */
  PetscReal    *v0;             /* The affine transformation for each subelement */
  PetscReal    *jac, *invjac;
  PetscInt     *embedding;      /* Map from subelements dofs to element dofs */
} PetscFE_Composite;

/* Utility functions */
PETSC_STATIC_INLINE void CoordinatesRefToReal(PetscInt dimReal, PetscInt dimRef, const PetscReal xi0[], const PetscReal v0[], const PetscReal J[], const PetscReal xi[], PetscReal x[])
{
  PetscInt d, e;

  for (d = 0; d < dimReal; ++d) {
    x[d] = v0[d];
    for (e = 0; e < dimRef; ++e) {
      x[d] += J[d*dimReal+e]*(xi[e] - xi0[e]);
    }
  }
}

PETSC_STATIC_INLINE void CoordinatesRealToRef(PetscInt dimReal, PetscInt dimRef, const PetscReal xi0[], const PetscReal v0[], const PetscReal invJ[], const PetscReal x[], PetscReal xi[])
{
  PetscInt d, e;

  for (d = 0; d < dimRef; ++d) {
    xi[d] = xi0[d];
    for (e = 0; e < dimReal; ++e) {
      xi[d] += invJ[d*dimReal+e]*(x[e] - v0[e]);
    }
  }
}

PETSC_STATIC_INLINE void EvaluateFieldJets(PetscInt dim, PetscInt Nf, const PetscInt Nb[], const PetscInt Nc[], PetscInt q, PetscReal *basisField[], PetscReal *basisFieldDer[], PetscScalar refSpaceDer[], const PetscReal invJ[], const PetscScalar coefficients[], const PetscScalar coefficients_t[], PetscScalar u[], PetscScalar u_x[], PetscScalar u_t[])
{
  PetscInt dOffset = 0, fOffset = 0, f;

  for (f = 0; f < Nf; ++f) {
    const PetscInt   Nbf = Nb[f], Ncf = Nc[f];
    const PetscReal *Bq = &basisField[f][q*Nbf*Ncf];
    const PetscReal *Dq = &basisFieldDer[f][q*Nbf*Ncf*dim];
    PetscInt         b, c, d, e;

    for (c = 0; c < Ncf; ++c)     u[fOffset+c] = 0.0;
    for (d = 0; d < dim*Ncf; ++d) refSpaceDer[d] = 0.0;
    for (b = 0; b < Nbf; ++b) {
      for (c = 0; c < Ncf; ++c) {
        const PetscInt cidx = b*Ncf+c;

        u[fOffset+c] += Bq[cidx]*coefficients[dOffset+b];
        for (d = 0; d < dim; ++d) refSpaceDer[c*dim+d] += Dq[cidx*dim+d]*coefficients[dOffset+b];
      }
    }
    for (c = 0; c < Ncf; ++c) for (d = 0; d < dim; ++d) for (e = 0, u_x[(fOffset+c)*dim+d] = 0.0; e < dim; ++e) u_x[(fOffset+c)*dim+d] += invJ[e*dim+d]*refSpaceDer[c*dim+e];
    if (u_t) {
      for (c = 0; c < Ncf; ++c) u_t[fOffset+c] = 0.0;
      for (b = 0; b < Nbf; ++b) {
        for (c = 0; c < Ncf; ++c) {
          const PetscInt cidx = b*Ncf+c;

          u_t[fOffset+c] += Bq[cidx]*coefficients_t[dOffset+b];
        }
      }
    }
#if 0
    for (c = 0; c < Ncf; ++c) {
      ierr = PetscPrintf(PETSC_COMM_SELF, "    u[%d,%d]: %g\n", f, c, PetscRealPart(u[fOffset+c]));CHKERRQ(ierr);
      if (u_t) {ierr = PetscPrintf(PETSC_COMM_SELF, "    u_t[%d,%d]: %g\n", f, c, PetscRealPart(u_t[fOffset+c]));CHKERRQ(ierr);}
      for (d = 0; d < dim; ++d) {
        ierr = PetscPrintf(PETSC_COMM_SELF, "    gradU[%d,%d]_%c: %g\n", f, c, 'x'+d, PetscRealPart(u_x[(fOffset+c)*dim+d]));CHKERRQ(ierr);
      }
    }
#endif
    fOffset += Ncf;
    dOffset += Nbf;
  }
}

PETSC_STATIC_INLINE PetscErrorCode EvaluateFaceFields(PetscDS prob, PetscInt field, PetscInt faceLoc, const PetscScalar coefficients[], PetscScalar u[])
{
  PetscFE        fe;
  PetscReal     *faceBasis;
  PetscInt       Nb, Nc, b, c;
  PetscErrorCode ierr;

  if (!prob) return 0;
  ierr = PetscDSGetDiscretization(prob, field, (PetscObject *) &fe);CHKERRQ(ierr);
  ierr = PetscFEGetDimension(fe, &Nb);CHKERRQ(ierr);
  ierr = PetscFEGetNumComponents(fe, &Nc);CHKERRQ(ierr);
  ierr = PetscFEGetFaceCentroidTabulation(fe, &faceBasis);CHKERRQ(ierr);
  for (c = 0; c < Nc; ++c) {u[c] = 0.0;}
  for (b = 0; b < Nb; ++b) {
    for (c = 0; c < Nc; ++c) {
      const PetscInt cidx = b*Nc+c;

      u[c] += coefficients[cidx]*faceBasis[faceLoc*Nb*Nc+cidx];
    }
  }
  return 0;
}

PETSC_STATIC_INLINE void TransformF(PetscInt dim, PetscInt Nc, PetscInt q, const PetscReal invJ[], PetscReal detJ, const PetscReal quadWeights[], PetscScalar refSpaceDer[], PetscScalar f0[], PetscScalar f1[])
{
  const PetscReal w = detJ*quadWeights[q];
  PetscInt        c, d, e;

  if (f0) for (c = 0; c < Nc; ++c) f0[q*Nc+c] *= w;
  if (f1) {
    for (c = 0; c < Nc; ++c) {
      for (d = 0; d < dim; ++d) {
        f1[(q*Nc + c)*dim+d] = 0.0;
        for (e = 0; e < dim; ++e) f1[(q*Nc + c)*dim+d] += invJ[d*dim+e]*refSpaceDer[c*dim+e];
        f1[(q*Nc + c)*dim+d] *= w;
      }
    }
  }
#if 0
  if (debug > 1) {
    for (c = 0; c < Nc; ++c) {
      ierr = PetscPrintf(PETSC_COMM_SELF, "    f0[%d]: %g\n", c, PetscRealPart(f0[q*Nc+c]));CHKERRQ(ierr);
      for (d = 0; d < dim; ++d) {
        ierr = PetscPrintf(PETSC_COMM_SELF, "    f1[%d]_%c: %g\n", c, 'x'+d, PetscRealPart(f1[(q*Nc + c)*dim+d]));CHKERRQ(ierr);
      }
    }
  }
#endif
}

PETSC_STATIC_INLINE void UpdateElementVec(PetscInt dim, PetscInt Nq, PetscInt Nb, PetscInt Nc, PetscReal basis[], PetscReal basisDer[], PetscScalar f0[], PetscScalar f1[], PetscScalar elemVec[])
{
  PetscInt b, c;

  for (b = 0; b < Nb; ++b) {
    elemVec[b] = 0.0;

    for (c = 0; c < Nc; ++c) {
      const PetscInt cidx = b*Nc+c;
      PetscInt       q;

      for (q = 0; q < Nq; ++q) {
        PetscInt d;

        elemVec[b] += basis[q*Nb*Nc+cidx]*f0[q*Nc+c];
        for (d = 0; d < dim; ++d) elemVec[b] += basisDer[(q*Nb*Nc+cidx)*dim+d]*f1[(q*Nc+c)*dim+d];
      }
    }
  }
#if 0
  if (debug > 1) {
    for (b = 0; b < Nb/Nc; ++b) {
      for (c = 0; c < Nc; ++c) {
        ierr = PetscPrintf(PETSC_COMM_SELF, "    elemVec[%d,%d]: %g\n", b, c, PetscRealPart(elemVec[b*Nc+c]));CHKERRQ(ierr);
      }
    }
  }
#endif
}

PETSC_STATIC_INLINE PetscErrorCode PetscFEInterpolate_Static(PetscFE fe, const PetscScalar x[], PetscInt q, PetscScalar interpolant[])
{
  PetscReal     *basis;
  PetscInt       Nb, Nc, fc, f;
  PetscErrorCode ierr;

  PetscFunctionBeginHot;
  ierr = PetscFEGetDimension(fe, &Nb);CHKERRQ(ierr);
  ierr = PetscFEGetNumComponents(fe, &Nc);CHKERRQ(ierr);
  ierr = PetscFEGetDefaultTabulation(fe, &basis, NULL, NULL);CHKERRQ(ierr);
  for (fc = 0; fc < Nc; ++fc) {
    interpolant[fc] = 0.0;
    for (f = 0; f < Nb; ++f) {
      interpolant[fc] += x[f]*basis[(q*Nb + f)*Nc + fc];
    }
  }
  PetscFunctionReturn(0);
}

PETSC_STATIC_INLINE PetscErrorCode PetscFEInterpolateGradient_Static(PetscFE fe, const PetscScalar x[], PetscInt dim, const PetscReal invJ[], const PetscReal n[], PetscInt q, PetscScalar interpolant[])
{
  PetscReal     *basisDer;
  PetscReal      realSpaceDer[3];
  PetscScalar    compGradient[3];
  PetscInt       Nb, Nc, fc, f, d, g;
  PetscErrorCode ierr;

  PetscFunctionBeginHot;
  ierr = PetscFEGetDimension(fe, &Nb);CHKERRQ(ierr);
  ierr = PetscFEGetNumComponents(fe, &Nc);CHKERRQ(ierr);
  ierr = PetscFEGetDefaultTabulation(fe, NULL, &basisDer, NULL);CHKERRQ(ierr);
  for (fc = 0; fc < Nc; ++fc) {
    interpolant[fc] = 0.0;
    for (d = 0; d < dim; ++d) compGradient[d] = 0.0;
    for (f = 0; f < Nb; ++f) {

      for (d = 0; d < dim; ++d) {
        realSpaceDer[d] = 0.0;
        for (g = 0; g < dim; ++g) {
          realSpaceDer[d] += invJ[dim*dim*q + g*dim+d]*basisDer[((q*Nb + f)*Nc + fc)*dim + g];
        }
        compGradient[d] += x[f]*realSpaceDer[d];
      }
    }
    if (n) {
      for (d = 0; d < dim; ++d) interpolant[fc] += compGradient[d]*n[d];
    } else {
      for (d = 0; d < dim; ++d) interpolant[fc*dim+d] = compGradient[d];
    }
  }
  PetscFunctionReturn(0);
}

PETSC_STATIC_INLINE PetscErrorCode PetscFEInterpolateFieldAndGradient_Static(PetscFE fe, const PetscScalar x[], PetscInt dim, const PetscReal invJ[], PetscInt q, PetscScalar interpolant[], PetscScalar interpolantGrad[])
{
  PetscReal     *basis, *basisDer;
  PetscReal      realSpaceDer[3];
  PetscInt       Nb, Nc, fc, f, d, g;
  PetscErrorCode ierr;

  PetscFunctionBeginHot;
  ierr = PetscFEGetDimension(fe, &Nb);CHKERRQ(ierr);
  ierr = PetscFEGetNumComponents(fe, &Nc);CHKERRQ(ierr);
  ierr = PetscFEGetDefaultTabulation(fe, &basis, &basisDer, NULL);CHKERRQ(ierr);
  for (fc = 0; fc < Nc; ++fc) {
    interpolant[fc] = 0.0;
    for (d = 0; d < dim; ++d) interpolantGrad[fc*dim+d] = 0.0;
    for (f = 0; f < Nb; ++f) {
      interpolant[fc] += x[f]*basis[(q*Nb + f)*Nc + fc];
      for (d = 0; d < dim; ++d) {
        realSpaceDer[d] = 0.0;
        for (g = 0; g < dim; ++g) {
          realSpaceDer[d] += invJ[dim*dim*q + g*dim+d]*basisDer[((q*Nb + f)*Nc + fc)*dim + g];
        }
        interpolantGrad[fc*dim+d] += x[f]*realSpaceDer[d];
      }
    }
  }
  PetscFunctionReturn(0);
}

PETSC_EXTERN PetscErrorCode PetscFEGetDimension_Basic(PetscFE, PetscInt *);
PETSC_EXTERN PetscErrorCode PetscFEIntegrateResidual_Basic(PetscFE, PetscDS, PetscInt, PetscInt, PetscFEGeom *, const PetscScalar [], const PetscScalar [], PetscDS, const PetscScalar [], PetscReal, PetscScalar []);
PETSC_EXTERN PetscErrorCode PetscFEIntegrateBdResidual_Basic(PetscFE, PetscDS, PetscInt, PetscInt, PetscFEGeom *, const PetscScalar [], const PetscScalar [], PetscDS, const PetscScalar [], PetscReal, PetscScalar[]);
PETSC_EXTERN PetscErrorCode PetscFEIntegrateJacobian_Basic(PetscFE, PetscDS, PetscFEJacobianType, PetscInt, PetscInt, PetscInt, PetscFEGeom *, const PetscScalar [], const PetscScalar [], PetscDS, const PetscScalar [], PetscReal, PetscReal, PetscScalar []);
#endif
