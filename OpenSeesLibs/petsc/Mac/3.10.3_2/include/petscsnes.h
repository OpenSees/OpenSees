/*
    User interface for the nonlinear solvers package.
*/
#if !defined(__PETSCSNES_H)
#define __PETSCSNES_H
#include <petscksp.h>
#include <petscdmtypes.h>
#include <petscfvtypes.h>
#include <petscdmdatypes.h>

/*S
     SNES - Abstract PETSc object that manages all nonlinear solves

   Level: beginner

  Concepts: nonlinear solvers

.seealso:  SNESCreate(), SNESSetType(), SNESType, TS, KSP, KSP, PC, SNESDestroy()
S*/
typedef struct _p_SNES* SNES;

/*J
    SNESType - String with the name of a PETSc SNES method.

   Level: beginner

.seealso: SNESSetType(), SNES, SNESCreate(), SNESDestroy(), SNESSetFromOptions()
J*/
typedef const char* SNESType;
#define SNESNEWTONLS     "newtonls"
#define SNESNEWTONTR     "newtontr"
#define SNESPYTHON       "python"
#define SNESNRICHARDSON  "nrichardson"
#define SNESKSPONLY      "ksponly"
#define SNESVINEWTONRSLS "vinewtonrsls"
#define SNESVINEWTONSSLS "vinewtonssls"
#define SNESNGMRES       "ngmres"
#define SNESQN           "qn"
#define SNESSHELL        "shell"
#define SNESNGS          "ngs"
#define SNESNCG          "ncg"
#define SNESFAS          "fas"
#define SNESMS           "ms"
#define SNESNASM         "nasm"
#define SNESANDERSON     "anderson"
#define SNESASPIN        "aspin"
#define SNESCOMPOSITE    "composite"

/* Logging support */
PETSC_EXTERN PetscClassId SNES_CLASSID;
PETSC_EXTERN PetscClassId DMSNES_CLASSID;

PETSC_EXTERN PetscErrorCode SNESInitializePackage(void);

PETSC_EXTERN PetscErrorCode SNESCreate(MPI_Comm,SNES*);
PETSC_EXTERN PetscErrorCode SNESReset(SNES);
PETSC_EXTERN PetscErrorCode SNESDestroy(SNES*);
PETSC_EXTERN PetscErrorCode SNESSetType(SNES,SNESType);
PETSC_EXTERN PetscErrorCode SNESMonitor(SNES,PetscInt,PetscReal);
PETSC_EXTERN PetscErrorCode SNESMonitorSet(SNES,PetscErrorCode(*)(SNES,PetscInt,PetscReal,void*),void *,PetscErrorCode (*)(void**));
PETSC_EXTERN PetscErrorCode SNESMonitorSetFromOptions(SNES,const char[],const char[],const char [],PetscErrorCode (*)(SNES,PetscInt,PetscReal,PetscViewerAndFormat*),PetscErrorCode (*)(SNES,PetscViewerAndFormat*));
PETSC_EXTERN PetscErrorCode SNESMonitorCancel(SNES);
PETSC_EXTERN PetscErrorCode SNESMonitorSAWs(SNES,PetscInt,PetscReal,void*);
PETSC_EXTERN PetscErrorCode SNESMonitorSAWsCreate(SNES,void**);
PETSC_EXTERN PetscErrorCode SNESMonitorSAWsDestroy(void**);
PETSC_EXTERN PetscErrorCode SNESSetConvergenceHistory(SNES,PetscReal[],PetscInt[],PetscInt,PetscBool );
PETSC_EXTERN PetscErrorCode SNESGetConvergenceHistory(SNES,PetscReal*[],PetscInt *[],PetscInt *);
PETSC_EXTERN PetscErrorCode SNESSetUp(SNES);
PETSC_EXTERN PetscErrorCode SNESSolve(SNES,Vec,Vec);
PETSC_EXTERN PetscErrorCode SNESSetErrorIfNotConverged(SNES,PetscBool );
PETSC_EXTERN PetscErrorCode SNESGetErrorIfNotConverged(SNES,PetscBool  *);

PETSC_EXTERN PetscErrorCode SNESSetWorkVecs(SNES,PetscInt);

PETSC_EXTERN PetscErrorCode SNESAddOptionsChecker(PetscErrorCode (*)(SNES));

PETSC_EXTERN PetscErrorCode SNESSetUpdate(SNES, PetscErrorCode (*)(SNES, PetscInt));


PETSC_EXTERN PetscErrorCode SNESRegister(const char[],PetscErrorCode (*)(SNES));

PETSC_EXTERN PetscErrorCode SNESGetKSP(SNES,KSP*);
PETSC_EXTERN PetscErrorCode SNESSetKSP(SNES,KSP);
PETSC_EXTERN PetscErrorCode SNESSetSolution(SNES,Vec);
PETSC_EXTERN PetscErrorCode SNESGetSolution(SNES,Vec*);
PETSC_EXTERN PetscErrorCode SNESGetSolutionUpdate(SNES,Vec*);
PETSC_EXTERN PetscErrorCode SNESGetRhs(SNES,Vec*);
PETSC_EXTERN PetscErrorCode SNESView(SNES,PetscViewer);
PETSC_EXTERN PetscErrorCode SNESLoad(SNES,PetscViewer);
PETSC_STATIC_INLINE PetscErrorCode SNESViewFromOptions(SNES A,PetscObject obj,const char name[]) {return PetscObjectViewFromOptions((PetscObject)A,obj,name);}
PETSC_EXTERN PetscErrorCode SNESReasonView(SNES,PetscViewer);
PETSC_EXTERN PetscErrorCode SNESReasonViewFromOptions(SNES);

#define SNES_FILE_CLASSID 1211224

PETSC_EXTERN PetscErrorCode SNESSetOptionsPrefix(SNES,const char[]);
PETSC_EXTERN PetscErrorCode SNESAppendOptionsPrefix(SNES,const char[]);
PETSC_EXTERN PetscErrorCode SNESGetOptionsPrefix(SNES,const char*[]);
PETSC_EXTERN PetscErrorCode SNESSetFromOptions(SNES);
PETSC_EXTERN PetscErrorCode SNESResetFromOptions(SNES);

PETSC_EXTERN PetscErrorCode SNESSetUseMatrixFree(SNES,PetscBool,PetscBool);
PETSC_EXTERN PetscErrorCode SNESGetUseMatrixFree(SNES,PetscBool*,PetscBool*);
PETSC_EXTERN PetscErrorCode MatCreateSNESMF(SNES,Mat*);
PETSC_EXTERN PetscErrorCode MatSNESMFGetSNES(Mat,SNES*);
PETSC_EXTERN PetscErrorCode MatSNESMFSetReuseBase(Mat,PetscBool);
PETSC_EXTERN PetscErrorCode MatSNESMFGetReuseBase(Mat,PetscBool*);
PETSC_EXTERN PetscErrorCode MatMFFDComputeJacobian(SNES,Vec,Mat,Mat,void*);

PETSC_EXTERN PetscErrorCode MatDAADSetSNES(Mat,SNES);

PETSC_EXTERN PetscErrorCode SNESGetType(SNES,SNESType*);
PETSC_EXTERN PetscErrorCode SNESMonitorDefault(SNES,PetscInt,PetscReal,PetscViewerAndFormat *);
PETSC_EXTERN PetscErrorCode SNESMonitorScaling(SNES,PetscInt,PetscReal,PetscViewerAndFormat *);
PETSC_EXTERN PetscErrorCode SNESMonitorRange(SNES,PetscInt,PetscReal,PetscViewerAndFormat *);
PETSC_EXTERN PetscErrorCode SNESMonitorRatio(SNES,PetscInt,PetscReal,PetscViewerAndFormat *);
PETSC_EXTERN PetscErrorCode SNESMonitorRatioSetUp(SNES,PetscViewerAndFormat*);
PETSC_EXTERN PetscErrorCode SNESMonitorSolution(SNES,PetscInt,PetscReal,PetscViewerAndFormat *);
PETSC_EXTERN PetscErrorCode SNESMonitorResidual(SNES,PetscInt,PetscReal,PetscViewerAndFormat *);
PETSC_EXTERN PetscErrorCode SNESMonitorSolutionUpdate(SNES,PetscInt,PetscReal,PetscViewerAndFormat *);
PETSC_EXTERN PetscErrorCode SNESMonitorDefaultShort(SNES,PetscInt,PetscReal,PetscViewerAndFormat *);
PETSC_EXTERN PetscErrorCode SNESMonitorDefaultField(SNES,PetscInt,PetscReal,PetscViewerAndFormat *);
PETSC_EXTERN PetscErrorCode SNESMonitorJacUpdateSpectrum(SNES,PetscInt,PetscReal,PetscViewerAndFormat *);
PETSC_EXTERN PetscErrorCode SNESMonitorFields(SNES,PetscInt,PetscReal,PetscViewerAndFormat *);
PETSC_EXTERN PetscErrorCode KSPMonitorSNES(KSP,PetscInt,PetscReal,void*);
PETSC_EXTERN PetscErrorCode KSPMonitorSNESLGResidualNormCreate(MPI_Comm,const char[],const char[],int,int,int,int,PetscObject**);
PETSC_EXTERN PetscErrorCode KSPMonitorSNESLGResidualNorm(KSP,PetscInt,PetscReal,PetscObject*);
PETSC_EXTERN PetscErrorCode KSPMonitorSNESLGResidualNormDestroy(PetscObject**);

PETSC_EXTERN PetscErrorCode SNESSetTolerances(SNES,PetscReal,PetscReal,PetscReal,PetscInt,PetscInt);
PETSC_EXTERN PetscErrorCode SNESSetDivergenceTolerance(SNES,PetscReal);
PETSC_EXTERN PetscErrorCode SNESGetTolerances(SNES,PetscReal*,PetscReal*,PetscReal*,PetscInt*,PetscInt*);
PETSC_EXTERN PetscErrorCode SNESGetDivergenceTolerance(SNES,PetscReal*);
PETSC_EXTERN PetscErrorCode SNESSetTrustRegionTolerance(SNES,PetscReal);
PETSC_EXTERN PetscErrorCode SNESGetForceIteration(SNES,PetscBool*);
PETSC_EXTERN PetscErrorCode SNESSetForceIteration(SNES,PetscBool);
PETSC_EXTERN PetscErrorCode SNESGetIterationNumber(SNES,PetscInt*);
PETSC_EXTERN PetscErrorCode SNESSetIterationNumber(SNES,PetscInt);

PETSC_EXTERN PetscErrorCode SNESGetNonlinearStepFailures(SNES,PetscInt*);
PETSC_EXTERN PetscErrorCode SNESSetMaxNonlinearStepFailures(SNES,PetscInt);
PETSC_EXTERN PetscErrorCode SNESGetMaxNonlinearStepFailures(SNES,PetscInt*);
PETSC_EXTERN PetscErrorCode SNESGetNumberFunctionEvals(SNES,PetscInt*);

PETSC_EXTERN PetscErrorCode SNESSetLagPreconditioner(SNES,PetscInt);
PETSC_EXTERN PetscErrorCode SNESGetLagPreconditioner(SNES,PetscInt*);
PETSC_EXTERN PetscErrorCode SNESSetLagJacobian(SNES,PetscInt);
PETSC_EXTERN PetscErrorCode SNESGetLagJacobian(SNES,PetscInt*);
PETSC_EXTERN PetscErrorCode SNESSetLagPreconditionerPersists(SNES,PetscBool);
PETSC_EXTERN PetscErrorCode SNESSetLagJacobianPersists(SNES,PetscBool);
PETSC_EXTERN PetscErrorCode SNESSetGridSequence(SNES,PetscInt);
PETSC_EXTERN PetscErrorCode SNESGetGridSequence(SNES,PetscInt*);

PETSC_EXTERN PetscErrorCode SNESGetLinearSolveIterations(SNES,PetscInt*);
PETSC_EXTERN PetscErrorCode SNESGetLinearSolveFailures(SNES,PetscInt*);
PETSC_EXTERN PetscErrorCode SNESSetMaxLinearSolveFailures(SNES,PetscInt);
PETSC_EXTERN PetscErrorCode SNESGetMaxLinearSolveFailures(SNES,PetscInt*);
PETSC_EXTERN PetscErrorCode SNESSetCountersReset(SNES,PetscBool);

PETSC_EXTERN PetscErrorCode SNESKSPSetUseEW(SNES,PetscBool );
PETSC_EXTERN PetscErrorCode SNESKSPGetUseEW(SNES,PetscBool *);
PETSC_EXTERN PetscErrorCode SNESKSPSetParametersEW(SNES,PetscInt,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal);
PETSC_EXTERN PetscErrorCode SNESKSPGetParametersEW(SNES,PetscInt*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscReal*);

#include <petscdrawtypes.h>
PETSC_EXTERN PetscErrorCode SNESMonitorLGCreate(MPI_Comm,const char[],const char[],int,int,int,int,PetscDrawLG*);
PETSC_EXTERN PetscErrorCode SNESMonitorLGResidualNorm(SNES,PetscInt,PetscReal,void*);
PETSC_EXTERN PetscErrorCode SNESMonitorLGRange(SNES,PetscInt,PetscReal,void*);

PETSC_EXTERN PetscErrorCode SNESSetApplicationContext(SNES,void *);
PETSC_EXTERN PetscErrorCode SNESGetApplicationContext(SNES,void *);
PETSC_EXTERN PetscErrorCode SNESSetComputeApplicationContext(SNES,PetscErrorCode (*)(SNES,void**),PetscErrorCode (*)(void**));

PETSC_EXTERN PetscErrorCode SNESPythonSetType(SNES,const char[]);

PETSC_EXTERN PetscErrorCode SNESSetFunctionDomainError(SNES);
PETSC_EXTERN PetscErrorCode SNESGetFunctionDomainError(SNES, PetscBool *);

/*E
    SNESConvergedReason - reason a SNES method was said to
         have converged or diverged

   Level: beginner

   The two most common reasons for divergence are
$   1) an incorrectly coded or computed Jacobian or
$   2) failure or lack of convergence in the linear system (in this case we recommend
$      testing with -pc_type lu to eliminate the linear solver as the cause of the problem).

   Diverged Reasons:
.    SNES_DIVERGED_LOCAL_MIN - this can only occur when using the line-search variant of SNES.
       The line search wants to minimize Q(alpha) = 1/2 || F(x + alpha s) ||^2_2  this occurs
       at Q'(alpha) = s^T F'(x+alpha s)^T F(x+alpha s) = 0. If s is the Newton direction - F'(x)^(-1)F(x) then
       you get Q'(alpha) = -F(x)^T F'(x)^(-1)^T F'(x+alpha s)F(x+alpha s); when alpha = 0
       Q'(0) = - ||F(x)||^2_2 which is always NEGATIVE if F'(x) is invertible. This means the Newton
       direction is a descent direction and the line search should succeed if alpha is small enough.

       If F'(x) is NOT invertible AND F'(x)^T F(x) = 0 then Q'(0) = 0 and the Newton direction
       is NOT a descent direction so the line search will fail. All one can do at this point
       is change the initial guess and try again.

       An alternative explanation: Newton's method can be regarded as replacing the function with
       its linear approximation and minimizing the 2-norm of that. That is F(x+s) approx F(x) + F'(x)s
       so we minimize || F(x) + F'(x) s ||^2_2; do this using Least Squares. If F'(x) is invertible then
       s = - F'(x)^(-1)F(x) otherwise F'(x)^T F'(x) s = -F'(x)^T F(x). If F'(x)^T F(x) is NOT zero then there
       exists a nontrival (that is F'(x)s != 0) solution to the equation and this direction is
       s = - [F'(x)^T F'(x)]^(-1) F'(x)^T F(x) so Q'(0) = - F(x)^T F'(x) [F'(x)^T F'(x)]^(-T) F'(x)^T F(x)
       = - (F'(x)^T F(x)) [F'(x)^T F'(x)]^(-T) (F'(x)^T F(x)). Since we are assuming (F'(x)^T F(x)) != 0
       and F'(x)^T F'(x) has no negative eigenvalues Q'(0) < 0 so s is a descent direction and the line
       search should succeed for small enough alpha.

       Note that this RARELY happens in practice. Far more likely the linear system is not being solved
       (well enough?) or the Jacobian is wrong.

   SNES_DIVERGED_MAX_IT means that the solver reached the maximum number of iterations without satisfying any
   convergence criteria. SNES_CONVERGED_ITS means that SNESConvergedSkip() was chosen as the convergence test;
   thus the usual convergence criteria have not been checked and may or may not be satisfied.

   Developer Notes:
    this must match petsc/finclude/petscsnes.h

       The string versions of these are in SNESConvergedReasons, if you change any value here you must
     also adjust that array.

   Each reason has its own manual page.

.seealso: SNESSolve(), SNESGetConvergedReason(), KSPConvergedReason, SNESSetConvergenceTest()
E*/
typedef enum {/* converged */
              SNES_CONVERGED_FNORM_ABS         =  2, /* ||F|| < atol */
              SNES_CONVERGED_FNORM_RELATIVE    =  3, /* ||F|| < rtol*||F_initial|| */
              SNES_CONVERGED_SNORM_RELATIVE    =  4, /* Newton computed step size small; || delta x || < stol || x ||*/
              SNES_CONVERGED_ITS               =  5, /* maximum iterations reached */
              SNES_CONVERGED_TR_DELTA          =  7,
              /* diverged */
              SNES_DIVERGED_FUNCTION_DOMAIN     = -1, /* the new x location passed the function is not in the domain of F */
              SNES_DIVERGED_FUNCTION_COUNT      = -2,
              SNES_DIVERGED_LINEAR_SOLVE        = -3, /* the linear solve failed */
              SNES_DIVERGED_FNORM_NAN           = -4,
              SNES_DIVERGED_MAX_IT              = -5,
              SNES_DIVERGED_LINE_SEARCH         = -6, /* the line search failed */
              SNES_DIVERGED_INNER               = -7, /* inner solve failed */
              SNES_DIVERGED_LOCAL_MIN           = -8, /* || J^T b || is small, implies converged to local minimum of F() */
              SNES_DIVERGED_DTOL                = -9, /* || F || > divtol*||F_initial|| */

              SNES_CONVERGED_ITERATING          =  0} SNESConvergedReason;
PETSC_EXTERN const char *const*SNESConvergedReasons;

/*MC
     SNES_CONVERGED_FNORM_ABS - 2-norm(F) <= abstol

   Level: beginner

.seealso:  SNESSolve(), SNESGetConvergedReason(), SNESConvergedReason, SNESSetTolerances()

M*/

/*MC
     SNES_CONVERGED_FNORM_RELATIVE - 2-norm(F) <= rtol*2-norm(F(x_0)) where x_0 is the initial guess

   Level: beginner

.seealso:  SNESSolve(), SNESGetConvergedReason(), SNESConvergedReason, SNESSetTolerances()

M*/

/*MC
     SNES_CONVERGED_SNORM_RELATIVE - The 2-norm of the last step <= stol * 2-norm(x) where x is the current
          solution and stol is the 4th argument to SNESSetTolerances()

     Options Database Keys:
      -snes_stol <stol> - the step tolerance

   Level: beginner

.seealso:  SNESSolve(), SNESGetConvergedReason(), SNESConvergedReason, SNESSetTolerances()

M*/

/*MC
     SNES_DIVERGED_FUNCTION_COUNT - The user provided function has been called more times then the final
         argument to SNESSetTolerances()

   Level: beginner

.seealso:  SNESSolve(), SNESGetConvergedReason(), SNESConvergedReason, SNESSetTolerances()

M*/

/*MC
     SNES_DIVERGED_DTOL - The norm of the function has increased by a factor of divtol set with SNESSetDivergenceTolerance()

   Level: beginner

.seealso:  SNESSolve(), SNESGetConvergedReason(), SNESConvergedReason, SNESSetTolerances(), SNESSetDivergenceTolerance()

M*/

/*MC
     SNES_DIVERGED_FNORM_NAN - the 2-norm of the current function evaluation is not-a-number (NaN), this
      is usually caused by a division of 0 by 0.

   Level: beginner

.seealso:  SNESSolve(), SNESGetConvergedReason(), SNESConvergedReason, SNESSetTolerances()

M*/

/*MC
     SNES_DIVERGED_MAX_IT - SNESSolve() has reached the maximum number of iterations requested

   Level: beginner

.seealso:  SNESSolve(), SNESGetConvergedReason(), SNESConvergedReason, SNESSetTolerances()

M*/

/*MC
     SNES_DIVERGED_LINE_SEARCH - The line search has failed. This only occurs for a SNES solvers that use a line search

   Level: beginner

.seealso:  SNESSolve(), SNESGetConvergedReason(), SNESConvergedReason, SNESSetTolerances(), SNESLineSearch

M*/

/*MC
     SNES_DIVERGED_LOCAL_MIN - the algorithm seems to have stagnated at a local minimum that is not zero.
        See the manual page for SNESConvergedReason for more details

   Level: beginner

.seealso:  SNESSolve(), SNESGetConvergedReason(), SNESConvergedReason, SNESSetTolerances()

M*/

/*MC
     SNES_CONERGED_ITERATING - this only occurs if SNESGetConvergedReason() is called during the SNESSolve()

   Level: beginner

.seealso:  SNESSolve(), SNESGetConvergedReason(), SNESConvergedReason, SNESSetTolerances()

M*/

PETSC_EXTERN PetscErrorCode SNESSetConvergenceTest(SNES,PetscErrorCode (*)(SNES,PetscInt,PetscReal,PetscReal,PetscReal,SNESConvergedReason*,void*),void*,PetscErrorCode (*)(void*));
PETSC_EXTERN PetscErrorCode SNESConvergedDefault(SNES,PetscInt,PetscReal,PetscReal,PetscReal,SNESConvergedReason*,void*);
PETSC_EXTERN PetscErrorCode SNESConvergedSkip(SNES,PetscInt,PetscReal,PetscReal,PetscReal,SNESConvergedReason*,void*);
PETSC_EXTERN PetscErrorCode SNESGetConvergedReason(SNES,SNESConvergedReason*);
PETSC_EXTERN PetscErrorCode SNESSetConvergedReason(SNES,SNESConvergedReason);

PETSC_DEPRECATED("Use SNESConvergedSkip()") PETSC_STATIC_INLINE void SNESSkipConverged(void) { /* never called */ }
#define SNESSkipConverged (SNESSkipConverged, SNESConvergedSkip)

/* --------- Solving systems of nonlinear equations --------------- */
PETSC_EXTERN PetscErrorCode SNESSetFunction(SNES,Vec,PetscErrorCode (*)(SNES,Vec,Vec,void*),void*);
PETSC_EXTERN PetscErrorCode SNESGetFunction(SNES,Vec*,PetscErrorCode (**)(SNES,Vec,Vec,void*),void**);
PETSC_EXTERN PetscErrorCode SNESComputeFunction(SNES,Vec,Vec);
PETSC_EXTERN PetscErrorCode SNESSetInitialFunction(SNES,Vec);

PETSC_EXTERN PetscErrorCode SNESSetJacobian(SNES,Mat,Mat,PetscErrorCode (*)(SNES,Vec,Mat,Mat,void*),void*);
PETSC_EXTERN PetscErrorCode SNESGetJacobian(SNES,Mat*,Mat*,PetscErrorCode (**)(SNES,Vec,Mat,Mat,void*),void**);
PETSC_EXTERN PetscErrorCode SNESObjectiveComputeFunctionDefaultFD(SNES,Vec,Vec,void*);
PETSC_EXTERN PetscErrorCode SNESComputeJacobianDefault(SNES,Vec,Mat,Mat,void*);
PETSC_EXTERN PetscErrorCode SNESComputeJacobianDefaultColor(SNES,Vec,Mat,Mat,void*);
PETSC_EXTERN PetscErrorCode SNESSetComputeInitialGuess(SNES,PetscErrorCode (*)(SNES,Vec,void*),void*);
PETSC_EXTERN PetscErrorCode SNESSetPicard(SNES,Vec,PetscErrorCode (*)(SNES,Vec,Vec,void*),Mat,Mat,PetscErrorCode (*)(SNES,Vec,Mat,Mat,void*),void*);
PETSC_EXTERN PetscErrorCode SNESGetPicard(SNES,Vec*,PetscErrorCode (**)(SNES,Vec,Vec,void*),Mat*,Mat*,PetscErrorCode (**)(SNES,Vec,Mat,Mat,void*),void**);
PETSC_EXTERN PetscErrorCode SNESPicardComputeFunction(SNES,Vec,Vec,void *);
PETSC_EXTERN PetscErrorCode SNESPicardComputeJacobian(SNES,Vec,Mat,Mat,void*);

PETSC_EXTERN PetscErrorCode SNESSetObjective(SNES,PetscErrorCode (*)(SNES,Vec,PetscReal *,void*),void*);
PETSC_EXTERN PetscErrorCode SNESGetObjective(SNES,PetscErrorCode (**)(SNES,Vec,PetscReal *,void*),void**);
PETSC_EXTERN PetscErrorCode SNESComputeObjective(SNES,Vec,PetscReal *);

/*E
    SNESNormSchedule - Frequency with which the norm is computed

   Level: advanced

   Support for these is highly dependent on the solver.

   Notes:
   This is primarily used to turn off extra norm and function computation
   when the solvers are composed.

.seealso: SNESSolve(), SNESGetConvergedReason(), KSPSetNormType(),
          KSPSetConvergenceTest(), KSPSetPCSide()
E*/

typedef enum {SNES_NORM_DEFAULT            = -1,
              SNES_NORM_NONE               =  0,
              SNES_NORM_ALWAYS             =  1,
              SNES_NORM_INITIAL_ONLY       =  2,
              SNES_NORM_FINAL_ONLY         =  3,
              SNES_NORM_INITIAL_FINAL_ONLY =  4} SNESNormSchedule;
PETSC_EXTERN const char *const*const SNESNormSchedules;

/*MC
    SNES_NORM_NONE - Don't compute function and its L2 norm.

   Level: advanced

    Notes:
    This is most useful for stationary solvers with a fixed number of iterations used as smoothers.

.seealso: SNESNormSchedule, SNESSetNormSchedule(), SNES_NORM_DEFAULT
M*/

/*MC
    SNES_NORM_ALWAYS - Compute the function and its L2 norm at each iteration.

   Level: advanced

    Notes:
    Most solvers will use this no matter what norm type is passed to them.

.seealso: SNESNormSchedule, SNESSetNormSchedule(), SNES_NORM_NONE
M*/

/*MC
    SNES_NORM_INITIAL_ONLY - Compute the function and its L2 at iteration 0, but do not update it.

   Level: advanced

   Notes:
   This method is useful in composed methods, when a true solution might actually be found before SNESSolve() is called.
   This option enables the solve to abort on the zeroth iteration if this is the case.

   For solvers that require the computation of the L2 norm of the function as part of the method, this merely cancels
   the norm computation at the last iteration (if possible).

.seealso: SNESNormSchedule, SNESSetNormSchedule(), SNES_NORM_FINAL_ONLY, SNES_NORM_INITIAL_FINAL_ONLY
M*/

/*MC
    SNES_NORM_FINAL_ONLY - Compute the function and its L2 norm on only the final iteration.

   Level: advanced

   Notes:
   For solvers that require the computation of the L2 norm of the function as part of the method, behaves
   exactly as SNES_NORM_DEFAULT.  This method is useful when the function is gotten after SNESSolve and
   used in subsequent computation for methods that do not need the norm computed during the rest of the
   solution procedure.

.seealso: SNESNormSchedule, SNESSetNormSchedule(), SNES_NORM_INITIAL_ONLY, SNES_NORM_INITIAL_FINAL_ONLY
M*/

/*MC
    SNES_NORM_INITIAL_FINAL_ONLY - Compute the function and its L2 norm on only the initial and final iterations.

   Level: advanced

   Notes:
   This method combines the benefits of SNES_NORM_INITIAL_ONLY and SNES_NORM_FINAL_ONLY.

.seealso: SNESNormSchedule, SNESSetNormSchedule(), SNES_NORM_SNES_NORM_INITIAL_ONLY, SNES_NORM_FINAL_ONLY
M*/

PETSC_EXTERN PetscErrorCode SNESSetNormSchedule(SNES,SNESNormSchedule);
PETSC_EXTERN PetscErrorCode SNESGetNormSchedule(SNES,SNESNormSchedule*);
PETSC_EXTERN PetscErrorCode SNESSetFunctionNorm(SNES,PetscReal);
PETSC_EXTERN PetscErrorCode SNESGetFunctionNorm(SNES,PetscReal*);

/*E
    SNESFunctionType - Type of function computed

   Level: advanced

   Support for these is highly dependent on the solver.

.seealso: SNESSolve(), SNESGetConvergedReason(), KSPSetNormType(),
          KSPSetConvergenceTest(), KSPSetPCSide()
E*/
typedef enum {SNES_FUNCTION_DEFAULT          = -1,
              SNES_FUNCTION_UNPRECONDITIONED =  0,
              SNES_FUNCTION_PRECONDITIONED   =  1} SNESFunctionType;
PETSC_EXTERN const char *const*const SNESFunctionTypes;

PETSC_EXTERN PetscErrorCode SNESSetFunctionType(SNES,SNESFunctionType);
PETSC_EXTERN PetscErrorCode SNESGetFunctionType(SNES,SNESFunctionType*);

PETSC_EXTERN PetscErrorCode SNESSetNGS(SNES,PetscErrorCode (*)(SNES,Vec,Vec,void*),void*);
PETSC_EXTERN PetscErrorCode SNESGetNGS(SNES,PetscErrorCode (**)(SNES,Vec,Vec,void*),void**);
PETSC_EXTERN PetscErrorCode SNESSetUseNGS(SNES,PetscBool);
PETSC_EXTERN PetscErrorCode SNESGetUseNGS(SNES,PetscBool *);
PETSC_EXTERN PetscErrorCode SNESComputeNGS(SNES,Vec,Vec);

PETSC_EXTERN PetscErrorCode SNESNGSSetSweeps(SNES,PetscInt);
PETSC_EXTERN PetscErrorCode SNESNGSGetSweeps(SNES,PetscInt *);
PETSC_EXTERN PetscErrorCode SNESNGSSetTolerances(SNES,PetscReal,PetscReal,PetscReal,PetscInt);
PETSC_EXTERN PetscErrorCode SNESNGSGetTolerances(SNES,PetscReal*,PetscReal*,PetscReal*,PetscInt*);

PETSC_EXTERN PetscErrorCode SNESUpdateCheckJacobian(SNES,PetscInt);

PETSC_EXTERN PetscErrorCode SNESSetAlwaysComputesFinalResidual(SNES,PetscBool);
PETSC_EXTERN PetscErrorCode SNESGetAlwaysComputesFinalResidual(SNES,PetscBool*);

PETSC_EXTERN PetscErrorCode SNESShellGetContext(SNES,void**);
PETSC_EXTERN PetscErrorCode SNESShellSetContext(SNES,void*);
PETSC_EXTERN PetscErrorCode SNESShellSetSolve(SNES,PetscErrorCode (*)(SNES,Vec));

/* --------- Routines specifically for line search methods --------------- */


/*S
     SNESLineSearch - Abstract PETSc object that manages line-search operations

   Level: beginner

  Concepts: nonlinear solvers, line search

.seealso:  SNESLineSearchCreate(), SNESLineSearchSetType(), SNES
S*/
typedef struct _p_LineSearch* SNESLineSearch;

/*J
    SNESLineSearchType - String with the name of a PETSc line search method

   Level: beginner

.seealso: SNESLineSearchSetType(), SNES
J*/
typedef const char* SNESLineSearchType;
#define SNESLINESEARCHBT                 "bt"
#define SNESLINESEARCHNLEQERR            "nleqerr"
#define SNESLINESEARCHBASIC              "basic"
#define SNESLINESEARCHL2                 "l2"
#define SNESLINESEARCHCP                 "cp"
#define SNESLINESEARCHSHELL              "shell"

PETSC_EXTERN PetscFunctionList SNESList;
PETSC_EXTERN PetscClassId      SNESLINESEARCH_CLASSID;
PETSC_EXTERN PetscFunctionList SNESLineSearchList;

#define SNES_LINESEARCH_ORDER_LINEAR    1
#define SNES_LINESEARCH_ORDER_QUADRATIC 2
#define SNES_LINESEARCH_ORDER_CUBIC     3

PETSC_EXTERN_TYPEDEF typedef PetscErrorCode (*SNESLineSearchVIProjectFunc)(SNES,Vec);
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode (*SNESLineSearchVINormFunc)(SNES,Vec,Vec,PetscReal *);
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode (*SNESLineSearchApplyFunc)(SNESLineSearch);
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode (*SNESLineSearchUserFunc)(SNESLineSearch, void *);

PETSC_EXTERN PetscErrorCode SNESLineSearchCreate(MPI_Comm, SNESLineSearch*);
PETSC_EXTERN PetscErrorCode SNESLineSearchReset(SNESLineSearch);
PETSC_EXTERN PetscErrorCode SNESLineSearchView(SNESLineSearch,PetscViewer);
PETSC_EXTERN PetscErrorCode SNESLineSearchDestroy(SNESLineSearch *);
PETSC_EXTERN PetscErrorCode SNESLineSearchSetType(SNESLineSearch, SNESLineSearchType);
PETSC_EXTERN PetscErrorCode SNESLineSearchSetFromOptions(SNESLineSearch);
PETSC_EXTERN PetscErrorCode SNESLineSearchSetFunction(SNESLineSearch,PetscErrorCode (*)(SNES,Vec,Vec));
PETSC_EXTERN PetscErrorCode SNESLineSearchSetUp(SNESLineSearch);
PETSC_EXTERN PetscErrorCode SNESLineSearchApply(SNESLineSearch, Vec, Vec, PetscReal *, Vec);
PETSC_EXTERN PetscErrorCode SNESLineSearchPreCheck(SNESLineSearch,Vec,Vec,PetscBool *);
PETSC_EXTERN PetscErrorCode SNESLineSearchPostCheck(SNESLineSearch,Vec,Vec,Vec,PetscBool *,PetscBool *);
PETSC_EXTERN PetscErrorCode SNESLineSearchSetWorkVecs(SNESLineSearch, PetscInt);

/* set the functions for precheck and postcheck */

PETSC_EXTERN PetscErrorCode SNESLineSearchSetPreCheck(SNESLineSearch, PetscErrorCode (*)(SNESLineSearch,Vec,Vec,PetscBool*,void*), void *ctx);
PETSC_EXTERN PetscErrorCode SNESLineSearchSetPostCheck(SNESLineSearch, PetscErrorCode (*)(SNESLineSearch,Vec,Vec,Vec,PetscBool *,PetscBool *,void*), void *ctx);

PETSC_EXTERN PetscErrorCode SNESLineSearchGetPreCheck(SNESLineSearch, PetscErrorCode (**)(SNESLineSearch,Vec,Vec,PetscBool*,void*), void **ctx);
PETSC_EXTERN PetscErrorCode SNESLineSearchGetPostCheck(SNESLineSearch, PetscErrorCode (**)(SNESLineSearch,Vec,Vec,Vec,PetscBool *,PetscBool *,void*), void **ctx);

/* set the functions for VI-specific line search operations */

PETSC_EXTERN PetscErrorCode SNESLineSearchSetVIFunctions(SNESLineSearch, SNESLineSearchVIProjectFunc, SNESLineSearchVINormFunc);
PETSC_EXTERN PetscErrorCode SNESLineSearchGetVIFunctions(SNESLineSearch, SNESLineSearchVIProjectFunc*, SNESLineSearchVINormFunc*);

/* pointers to the associated SNES in order to be able to get the function evaluation out */
PETSC_EXTERN PetscErrorCode SNESLineSearchSetSNES(SNESLineSearch,SNES);
PETSC_EXTERN PetscErrorCode SNESLineSearchGetSNES(SNESLineSearch,SNES*);

/* set and get the parameters and vectors */
PETSC_EXTERN PetscErrorCode SNESLineSearchGetTolerances(SNESLineSearch,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscInt*);
PETSC_EXTERN PetscErrorCode SNESLineSearchSetTolerances(SNESLineSearch,PetscReal,PetscReal,PetscReal,PetscReal,PetscReal,PetscInt);

PETSC_EXTERN PetscErrorCode SNESLineSearchPreCheckPicard(SNESLineSearch,Vec,Vec,PetscBool*,void*);

PETSC_EXTERN PetscErrorCode SNESLineSearchGetLambda(SNESLineSearch,PetscReal*);
PETSC_EXTERN PetscErrorCode SNESLineSearchSetLambda(SNESLineSearch,PetscReal);

PETSC_EXTERN PetscErrorCode SNESLineSearchGetDamping(SNESLineSearch,PetscReal*);
PETSC_EXTERN PetscErrorCode SNESLineSearchSetDamping(SNESLineSearch,PetscReal);

PETSC_EXTERN PetscErrorCode SNESLineSearchGetOrder(SNESLineSearch,PetscInt *order);
PETSC_EXTERN PetscErrorCode SNESLineSearchSetOrder(SNESLineSearch,PetscInt order);

/*E
    SNESLineSearchReason - if line search has succeeded or failed and why

   Level: intermediate

   Developer Notes:
    this must match petsc/finclude/petscsnes.h

   Developer Note: The string versions of these are in SNESLineSearchReasons, if you change any value here you must
     also adjust that array.

.seealso: SNESSolve(), SNESGetConvergedReason(), KSPConvergedReason, SNESSetConvergenceTest()
E*/
typedef enum {SNES_LINESEARCH_SUCCEEDED,
              SNES_LINESEARCH_FAILED_NANORINF,
              SNES_LINESEARCH_FAILED_DOMAIN,
              SNES_LINESEARCH_FAILED_REDUCT,       /* INSUFFICENT REDUCTION */
              SNES_LINESEARCH_FAILED_USER,
              SNES_LINESEARCH_FAILED_FUNCTION} SNESLineSearchReason;

PETSC_EXTERN PetscErrorCode SNESLineSearchGetReason(SNESLineSearch, SNESLineSearchReason*);
PETSC_EXTERN PetscErrorCode SNESLineSearchSetReason(SNESLineSearch, SNESLineSearchReason);

PETSC_EXTERN PetscErrorCode SNESLineSearchGetVecs(SNESLineSearch,Vec*,Vec*,Vec*,Vec*,Vec*);
PETSC_EXTERN PetscErrorCode SNESLineSearchSetVecs(SNESLineSearch,Vec,Vec,Vec,Vec,Vec);

PETSC_EXTERN PetscErrorCode SNESLineSearchGetNorms(SNESLineSearch, PetscReal *, PetscReal *, PetscReal *);
PETSC_EXTERN PetscErrorCode SNESLineSearchSetNorms(SNESLineSearch, PetscReal, PetscReal, PetscReal);
PETSC_EXTERN PetscErrorCode SNESLineSearchComputeNorms(SNESLineSearch);
PETSC_EXTERN PetscErrorCode SNESLineSearchSetComputeNorms(SNESLineSearch, PetscBool);

PETSC_EXTERN PetscErrorCode SNESLineSearchMonitor(SNESLineSearch);
PETSC_EXTERN PetscErrorCode SNESLineSearchMonitorSet(SNESLineSearch,PetscErrorCode(*)(SNESLineSearch,void*),void *,PetscErrorCode (*)(void**));
PETSC_EXTERN PetscErrorCode SNESLineSearchMonitorSetFromOptions(SNESLineSearch,const char[],const char[],const char[],PetscErrorCode (*)(SNESLineSearch,PetscViewerAndFormat*),PetscErrorCode (*)(SNESLineSearch,PetscViewerAndFormat*));
PETSC_EXTERN PetscErrorCode SNESLineSearchMonitorCancel(SNESLineSearch);
PETSC_EXTERN PetscErrorCode SNESLineSearchMonitorUpdate(SNESLineSearch,PetscViewerAndFormat*);
PETSC_EXTERN PetscErrorCode SNESLineSearchSetDefaultMonitor(SNESLineSearch,PetscViewer);
PETSC_EXTERN PetscErrorCode SNESLineSearchGetDefaultMonitor(SNESLineSearch,PetscViewer*);
PETSC_EXTERN PetscErrorCode SNESLineSearchMonitorSolutionUpdate(SNESLineSearch,PetscViewerAndFormat*);

PETSC_EXTERN PetscErrorCode SNESLineSearchAppendOptionsPrefix(SNESLineSearch, const char prefix[]);
PETSC_EXTERN PetscErrorCode SNESLineSearchGetOptionsPrefix(SNESLineSearch, const char *prefix[]);


/* Shell interface functions */
PETSC_EXTERN PetscErrorCode SNESLineSearchShellSetUserFunc(SNESLineSearch,SNESLineSearchUserFunc,void*);
PETSC_EXTERN PetscErrorCode SNESLineSearchShellGetUserFunc(SNESLineSearch,SNESLineSearchUserFunc*,void**);

/* BT interface functions */
PETSC_EXTERN PetscErrorCode SNESLineSearchBTSetAlpha(SNESLineSearch, PetscReal);
PETSC_EXTERN PetscErrorCode SNESLineSearchBTGetAlpha(SNESLineSearch, PetscReal*);

/*register line search types */
PETSC_EXTERN PetscErrorCode SNESLineSearchRegister(const char[],PetscErrorCode(*)(SNESLineSearch));

/* Routines for VI solver */
PETSC_EXTERN PetscErrorCode SNESVISetVariableBounds(SNES,Vec,Vec);
PETSC_EXTERN PetscErrorCode SNESVISetComputeVariableBounds(SNES, PetscErrorCode (*)(SNES,Vec,Vec));
PETSC_EXTERN PetscErrorCode SNESVIGetInactiveSet(SNES,IS*);
PETSC_EXTERN PetscErrorCode SNESVIGetActiveSetIS(SNES,Vec,Vec,IS*);
PETSC_EXTERN PetscErrorCode SNESVIComputeInactiveSetFnorm(SNES,Vec,Vec,PetscReal*);
PETSC_EXTERN PetscErrorCode SNESVISetRedundancyCheck(SNES,PetscErrorCode(*)(SNES,IS,IS*,void*),void*);

PETSC_EXTERN PetscErrorCode SNESTestLocalMin(SNES);

/* Should this routine be private? */
PETSC_EXTERN PetscErrorCode SNESComputeJacobian(SNES,Vec,Mat,Mat);
PETSC_EXTERN PetscErrorCode SNESTestJacobian(SNES);

PETSC_EXTERN PetscErrorCode SNESSetDM(SNES,DM);
PETSC_EXTERN PetscErrorCode SNESGetDM(SNES,DM*);
PETSC_EXTERN PetscErrorCode SNESSetNPC(SNES,SNES);
PETSC_EXTERN PetscErrorCode SNESGetNPC(SNES,SNES*);
PETSC_EXTERN PetscErrorCode SNESHasNPC(SNES,PetscBool*);
PETSC_EXTERN PetscErrorCode SNESApplyNPC(SNES,Vec,Vec,Vec);
PETSC_EXTERN PetscErrorCode SNESGetNPCFunction(SNES,Vec,PetscReal*);
PETSC_EXTERN PetscErrorCode SNESComputeFunctionDefaultNPC(SNES,Vec,Vec);
PETSC_EXTERN PetscErrorCode SNESSetNPCSide(SNES,PCSide);
PETSC_EXTERN PetscErrorCode SNESGetNPCSide(SNES,PCSide*);
PETSC_EXTERN PetscErrorCode SNESSetLineSearch(SNES,SNESLineSearch);
PETSC_EXTERN PetscErrorCode SNESGetLineSearch(SNES,SNESLineSearch*);
PETSC_EXTERN PetscErrorCode SNESRestrictHookAdd(SNES,PetscErrorCode (*)(SNES,SNES,void*),void*);
PETSC_EXTERN PetscErrorCode SNESRestrictHooksRun(SNES,SNES);

PETSC_DEPRECATED("Use SNESGetLineSearch()") PETSC_STATIC_INLINE PetscErrorCode SNESGetSNESLineSearch(SNES snes,SNESLineSearch *ls) {return SNESGetLineSearch(snes,ls);}
PETSC_DEPRECATED("Use SNESSetLineSearch()") PETSC_STATIC_INLINE PetscErrorCode SNESSetSNESLineSearch(SNES snes,SNESLineSearch ls) {return SNESSetLineSearch(snes,ls);}

PETSC_EXTERN PetscErrorCode SNESSetUpMatrices(SNES);
PETSC_EXTERN PetscErrorCode DMSNESSetFunction(DM,PetscErrorCode(*)(SNES,Vec,Vec,void*),void*);
PETSC_EXTERN PetscErrorCode DMSNESGetFunction(DM,PetscErrorCode(**)(SNES,Vec,Vec,void*),void**);
PETSC_EXTERN PetscErrorCode DMSNESSetNGS(DM,PetscErrorCode(*)(SNES,Vec,Vec,void*),void*);
PETSC_EXTERN PetscErrorCode DMSNESGetNGS(DM,PetscErrorCode(**)(SNES,Vec,Vec,void*),void**);
PETSC_EXTERN PetscErrorCode DMSNESSetJacobian(DM,PetscErrorCode(*)(SNES,Vec,Mat,Mat,void*),void*);
PETSC_EXTERN PetscErrorCode DMSNESGetJacobian(DM,PetscErrorCode(**)(SNES,Vec,Mat,Mat,void*),void**);
PETSC_EXTERN PetscErrorCode DMSNESSetPicard(DM,PetscErrorCode(*)(SNES,Vec,Vec,void*),PetscErrorCode(*)(SNES,Vec,Mat,Mat,void*),void*);
PETSC_EXTERN PetscErrorCode DMSNESGetPicard(DM,PetscErrorCode(**)(SNES,Vec,Vec,void*),PetscErrorCode(**)(SNES,Vec,Mat,Mat,void*),void**);
PETSC_EXTERN PetscErrorCode DMSNESSetObjective(DM,PetscErrorCode (*)(SNES,Vec,PetscReal *,void*),void*);
PETSC_EXTERN PetscErrorCode DMSNESGetObjective(DM,PetscErrorCode (**)(SNES,Vec,PetscReal *,void*),void**);
PETSC_EXTERN PetscErrorCode DMCopyDMSNES(DM,DM);

PETSC_EXTERN_TYPEDEF typedef PetscErrorCode (*DMDASNESFunction)(DMDALocalInfo*,void*,void*,void*);
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode (*DMDASNESJacobian)(DMDALocalInfo*,void*,Mat,Mat,void*);
PETSC_EXTERN_TYPEDEF typedef PetscErrorCode (*DMDASNESObjective)(DMDALocalInfo*,void*,PetscReal*,void*);

PETSC_EXTERN PetscErrorCode DMDASNESSetFunctionLocal(DM,InsertMode,DMDASNESFunction,void*);
PETSC_EXTERN PetscErrorCode DMDASNESSetJacobianLocal(DM,DMDASNESJacobian,void*);
PETSC_EXTERN PetscErrorCode DMDASNESSetObjectiveLocal(DM,DMDASNESObjective,void*);
PETSC_EXTERN PetscErrorCode DMDASNESSetPicardLocal(DM,InsertMode,PetscErrorCode (*)(DMDALocalInfo*,void*,void*,void*),PetscErrorCode (*)(DMDALocalInfo*,void*,Mat,Mat,void*),void*);

PETSC_EXTERN PetscErrorCode DMPlexSNESGetGeometryFVM(DM,Vec*,Vec*,PetscReal*);
PETSC_EXTERN PetscErrorCode DMPlexSNESGetGradientDM(DM,PetscFV,DM*);
PETSC_EXTERN PetscErrorCode DMPlexGetCellFields(DM, IS, Vec, Vec, Vec, PetscScalar **, PetscScalar **, PetscScalar **);
PETSC_EXTERN PetscErrorCode DMPlexRestoreCellFields(DM, IS, Vec, Vec, Vec, PetscScalar **, PetscScalar **, PetscScalar **);
PETSC_EXTERN PetscErrorCode DMPlexGetFaceFields(DM, PetscInt, PetscInt, Vec, Vec, Vec, Vec, Vec, PetscInt *, PetscScalar **, PetscScalar **);
PETSC_EXTERN PetscErrorCode DMPlexRestoreFaceFields(DM, PetscInt, PetscInt, Vec, Vec, Vec, Vec, Vec, PetscInt *, PetscScalar **, PetscScalar **);
PETSC_EXTERN PetscErrorCode DMPlexGetFaceGeometry(DM, PetscInt, PetscInt, Vec, Vec, PetscInt *, PetscFVFaceGeom **, PetscReal **);
PETSC_EXTERN PetscErrorCode DMPlexRestoreFaceGeometry(DM, PetscInt, PetscInt, Vec, Vec, PetscInt *, PetscFVFaceGeom **, PetscReal **);

PETSC_EXTERN PetscErrorCode DMSNESSetBoundaryLocal(DM,PetscErrorCode (*)(DM,Vec,void*),void*);
PETSC_EXTERN PetscErrorCode DMSNESSetFunctionLocal(DM,PetscErrorCode (*)(DM,Vec,Vec,void*),void*);
PETSC_EXTERN PetscErrorCode DMSNESSetJacobianLocal(DM,PetscErrorCode (*)(DM,Vec,Mat,Mat,void*),void*);

/* Routines for Multiblock solver */
PETSC_EXTERN PetscErrorCode SNESMultiblockSetFields(SNES, const char [], PetscInt, const PetscInt *);
PETSC_EXTERN PetscErrorCode SNESMultiblockSetIS(SNES, const char [], IS);
PETSC_EXTERN PetscErrorCode SNESMultiblockSetBlockSize(SNES, PetscInt);
PETSC_EXTERN PetscErrorCode SNESMultiblockSetType(SNES, PCCompositeType);

/*J
    SNESMSType - String with the name of a PETSc SNESMS method.

   Level: intermediate

.seealso: SNESMSSetType(), SNES
J*/
typedef const char* SNESMSType;
#define SNESMSM62       "m62"
#define SNESMSEULER     "euler"
#define SNESMSJAMESON83 "jameson83"
#define SNESMSVLTP21    "vltp21"
#define SNESMSVLTP31    "vltp31"
#define SNESMSVLTP41    "vltp41"
#define SNESMSVLTP51    "vltp51"
#define SNESMSVLTP61    "vltp61"

PETSC_EXTERN PetscErrorCode SNESMSRegister(SNESMSType,PetscInt,PetscInt,PetscReal,const PetscReal[],const PetscReal[],const PetscReal[]);
PETSC_EXTERN PetscErrorCode SNESMSSetType(SNES,SNESMSType);
PETSC_EXTERN PetscErrorCode SNESMSFinalizePackage(void);
PETSC_EXTERN PetscErrorCode SNESMSInitializePackage(void);
PETSC_EXTERN PetscErrorCode SNESMSRegisterDestroy(void);

/* routines for NGMRES solver */

typedef enum {
  SNES_NGMRES_RESTART_NONE       = 0,
  SNES_NGMRES_RESTART_PERIODIC   = 1,
  SNES_NGMRES_RESTART_DIFFERENCE = 2} SNESNGMRESRestartType;
PETSC_EXTERN const char *const SNESNGMRESRestartTypes[];

typedef enum {
  SNES_NGMRES_SELECT_NONE       = 0,
  SNES_NGMRES_SELECT_DIFFERENCE = 1,
  SNES_NGMRES_SELECT_LINESEARCH = 2} SNESNGMRESSelectType;
PETSC_EXTERN const char *const SNESNGMRESSelectTypes[];

PETSC_EXTERN PetscErrorCode SNESNGMRESSetRestartType(SNES, SNESNGMRESRestartType);
PETSC_EXTERN PetscErrorCode SNESNGMRESSetSelectType(SNES, SNESNGMRESSelectType);
PETSC_EXTERN PetscErrorCode SNESNGMRESSetRestartFmRise(SNES, PetscBool);
PETSC_EXTERN PetscErrorCode SNESNGMRESGetRestartFmRise(SNES, PetscBool*);

/* routines for NCG solver */

typedef enum {
  SNES_NCG_FR    = 0,
  SNES_NCG_PRP   = 1,
  SNES_NCG_HS    = 2,
  SNES_NCG_DY    = 3,
  SNES_NCG_CD    = 4} SNESNCGType;
PETSC_EXTERN const char *const SNESNCGTypes[];

PETSC_EXTERN PetscErrorCode SNESNCGSetType(SNES, SNESNCGType);

typedef enum {SNES_QN_SCALE_DEFAULT    = 0,
              SNES_QN_SCALE_NONE       = 1,
              SNES_QN_SCALE_SHANNO     = 2,
              SNES_QN_SCALE_LINESEARCH = 3,
              SNES_QN_SCALE_JACOBIAN   = 4} SNESQNScaleType;
PETSC_EXTERN const char *const SNESQNScaleTypes[];
typedef enum {SNES_QN_RESTART_DEFAULT  = 0,
              SNES_QN_RESTART_NONE     = 1,
              SNES_QN_RESTART_POWELL   = 2,
              SNES_QN_RESTART_PERIODIC = 3} SNESQNRestartType;
PETSC_EXTERN const char *const SNESQNRestartTypes[];
typedef enum {SNES_QN_LBFGS      = 0,
              SNES_QN_BROYDEN    = 1,
              SNES_QN_BADBROYDEN = 2
             } SNESQNType;
PETSC_EXTERN const char *const SNESQNTypes[];

PETSC_EXTERN PetscErrorCode SNESQNSetType(SNES, SNESQNType);
PETSC_EXTERN PetscErrorCode SNESQNSetScaleType(SNES, SNESQNScaleType);
PETSC_EXTERN PetscErrorCode SNESQNSetRestartType(SNES, SNESQNRestartType);

PETSC_EXTERN PetscErrorCode SNESNASMGetType(SNES,PCASMType*);
PETSC_EXTERN PetscErrorCode SNESNASMSetType(SNES,PCASMType);
PETSC_EXTERN PetscErrorCode SNESNASMGetSubdomains(SNES,PetscInt*,SNES**,VecScatter**,VecScatter**,VecScatter**);
PETSC_EXTERN PetscErrorCode SNESNASMSetSubdomains(SNES,PetscInt,SNES*,VecScatter*,VecScatter*,VecScatter*);
PETSC_EXTERN PetscErrorCode SNESNASMSetDamping(SNES,PetscReal);
PETSC_EXTERN PetscErrorCode SNESNASMGetDamping(SNES,PetscReal*);
PETSC_EXTERN PetscErrorCode SNESNASMGetSubdomainVecs(SNES,PetscInt*,Vec**,Vec**,Vec**,Vec**);
PETSC_EXTERN PetscErrorCode SNESNASMSetComputeFinalJacobian(SNES,PetscBool);
PETSC_EXTERN PetscErrorCode SNESNASMGetSNES(SNES,PetscInt,SNES *);
PETSC_EXTERN PetscErrorCode SNESNASMGetNumber(SNES,PetscInt*);
PETSC_EXTERN PetscErrorCode SNESNASMSetWeight(SNES,Vec);

typedef enum {SNES_COMPOSITE_ADDITIVE,SNES_COMPOSITE_MULTIPLICATIVE,SNES_COMPOSITE_ADDITIVEOPTIMAL} SNESCompositeType;
PETSC_EXTERN const char *const SNESCompositeTypes[];

PETSC_EXTERN PetscErrorCode SNESCompositeSetType(SNES,SNESCompositeType);
PETSC_EXTERN PetscErrorCode SNESCompositeAddSNES(SNES,SNESType);
PETSC_EXTERN PetscErrorCode SNESCompositeGetSNES(SNES,PetscInt,SNES *);
PETSC_EXTERN PetscErrorCode SNESCompositeGetNumber(SNES,PetscInt*);
PETSC_EXTERN PetscErrorCode SNESCompositeSetDamping(SNES,PetscInt,PetscReal);

/*E
    SNESFASType - Determines the type of nonlinear multigrid method that is run.

   Level: beginner

   Values:
+  SNES_FAS_MULTIPLICATIVE (default) - traditional V or W cycle as determined by SNESFASSetCycles()
.  SNES_FAS_ADDITIVE                 - additive FAS cycle
.  SNES_FAS_FULL                     - full FAS cycle
-  SNES_FAS_KASKADE                  - Kaskade FAS cycle
.seealso: PCMGSetType(), PCMGType

E*/
typedef enum { SNES_FAS_MULTIPLICATIVE, SNES_FAS_ADDITIVE, SNES_FAS_FULL, SNES_FAS_KASKADE } SNESFASType;
PETSC_EXTERN const char *const  SNESFASTypes[];

/* called on the finest level FAS instance*/
PETSC_EXTERN PetscErrorCode SNESFASSetType(SNES, SNESFASType);
PETSC_EXTERN PetscErrorCode SNESFASGetType(SNES, SNESFASType*);
PETSC_EXTERN PetscErrorCode SNESFASSetLevels(SNES, PetscInt, MPI_Comm *);
PETSC_EXTERN PetscErrorCode SNESFASGetLevels(SNES, PetscInt *);
PETSC_EXTERN PetscErrorCode SNESFASGetCycleSNES(SNES, PetscInt, SNES*);
PETSC_EXTERN PetscErrorCode SNESFASSetNumberSmoothUp(SNES, PetscInt);
PETSC_EXTERN PetscErrorCode SNESFASSetNumberSmoothDown(SNES, PetscInt);
PETSC_EXTERN PetscErrorCode SNESFASSetCycles(SNES, PetscInt);
PETSC_EXTERN PetscErrorCode SNESFASSetMonitor(SNES, PetscViewerAndFormat *, PetscBool);
PETSC_EXTERN PetscErrorCode SNESFASSetLog(SNES, PetscBool);

PETSC_EXTERN PetscErrorCode SNESFASSetGalerkin(SNES, PetscBool);
PETSC_EXTERN PetscErrorCode SNESFASGetGalerkin(SNES, PetscBool*);
PETSC_EXTERN PetscErrorCode SNESFASGalerkinFunctionDefault(SNES,Vec,Vec,void*);

/* called on any level -- "Cycle" FAS instance */
PETSC_EXTERN PetscErrorCode SNESFASCycleGetSmoother(SNES, SNES*);
PETSC_EXTERN PetscErrorCode SNESFASCycleGetSmootherUp(SNES, SNES*);
PETSC_EXTERN PetscErrorCode SNESFASCycleGetSmootherDown(SNES, SNES*);
PETSC_EXTERN PetscErrorCode SNESFASCycleGetCorrection(SNES, SNES*);
PETSC_EXTERN PetscErrorCode SNESFASCycleGetInterpolation(SNES, Mat*);
PETSC_EXTERN PetscErrorCode SNESFASCycleGetRestriction(SNES, Mat*);
PETSC_EXTERN PetscErrorCode SNESFASCycleGetInjection(SNES, Mat*);
PETSC_EXTERN PetscErrorCode SNESFASCycleGetRScale(SNES, Vec*);
PETSC_EXTERN PetscErrorCode SNESFASCycleSetCycles(SNES, PetscInt);
PETSC_EXTERN PetscErrorCode SNESFASCycleIsFine(SNES, PetscBool*);

/* called on the (outer) finest level FAS to set/get parameters on any level instance */
PETSC_EXTERN PetscErrorCode SNESFASSetInterpolation(SNES, PetscInt, Mat);
PETSC_EXTERN PetscErrorCode SNESFASGetInterpolation(SNES, PetscInt, Mat*);
PETSC_EXTERN PetscErrorCode SNESFASSetRestriction(SNES, PetscInt, Mat);
PETSC_EXTERN PetscErrorCode SNESFASGetRestriction(SNES, PetscInt, Mat*);
PETSC_EXTERN PetscErrorCode SNESFASSetInjection(SNES, PetscInt, Mat);
PETSC_EXTERN PetscErrorCode SNESFASGetInjection(SNES, PetscInt, Mat*);
PETSC_EXTERN PetscErrorCode SNESFASSetRScale(SNES, PetscInt, Vec);
PETSC_EXTERN PetscErrorCode SNESFASGetRScale(SNES, PetscInt, Vec*);
PETSC_EXTERN PetscErrorCode SNESFASSetContinuation(SNES,PetscBool);

PETSC_EXTERN PetscErrorCode SNESFASGetSmoother(SNES,     PetscInt, SNES*);
PETSC_EXTERN PetscErrorCode SNESFASGetSmootherUp(SNES,   PetscInt, SNES*);
PETSC_EXTERN PetscErrorCode SNESFASGetSmootherDown(SNES, PetscInt, SNES*);
PETSC_EXTERN PetscErrorCode SNESFASGetCoarseSolve(SNES, SNES*);

/* parameters for full FAS */
PETSC_EXTERN PetscErrorCode SNESFASFullSetDownSweep(SNES,PetscBool);
PETSC_EXTERN PetscErrorCode SNESFASCreateCoarseVec(SNES,Vec*);
PETSC_EXTERN PetscErrorCode SNESFASRestrict(SNES,Vec,Vec);

PETSC_EXTERN PetscErrorCode DMSNESCheckFromOptions(SNES,Vec,PetscErrorCode (**)(PetscInt,PetscReal,const PetscReal[],PetscInt,PetscScalar*,void*),void**);

#endif
