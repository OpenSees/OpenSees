/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Written: Fernando Sarabia, Ph.D., P.E.
// Created: 08/2026
//
// Description: This file contains the class definition for UmatMaterial.
// UmatMaterial is a generic wrapper NDMaterial that drives an external
// Abaqus-UMAT Fortran subroutine.  The wrapper itself contains NO
// model-specific code: the UMAT entry symbol, the material constants
// (PROPS) and the state variables (STATEV) are all supplied on the
// command line, so arbitrary (including GPL-licensed) user materials
// can be compiled separately and bound at run time.
//
// Scope: a 3D small-strain UMAT adapter, not an Abaqus emulation.  UMAT
// quantities with a meaningful OpenSees equivalent (STRAN, DSTRAN, STRESS,
// STATEV, PROPS, DDSDDE, TIME/DTIME, PNEWDT, CELENT) are translated
// directly; Abaqus analysis metadata that the NDMaterial interface cannot
// supply (COORDS, DROT, DFGRD0/1, NOEL, NPT, LAYER, KSPT, JSTEP) is passed
// as documented, deterministic placeholders (see invokeUmat).  UMATs that
// base their response on those quantities are outside the supported scope.
//
//   nDMaterial UMAT $tag -symbol <fortran_symbol> -nprops N -props p1 .. pN
//                        -nstatv M <-statev v1 .. vM>
//                        <-initial_stress s11 s22 s33 s12 s13 s23>
//                        <-library /path/to/lib.so> <-rho r>
//                        <-check_tangent h>
//
//   M may be zero for stateless UMATs; -statev is then omitted.
//   -check_tangent enables a development-only central-difference check of
//   DDSDDE with strain perturbation h (diagnostic report only; the FD
//   tangent is never substituted and never fails the analysis).
//
// Symbol resolution:
//   - with -library: dlopen(library, RTLD_NOW | RTLD_GLOBAL) then dlsym
//   - without: dlsym(RTLD_DEFAULT, symbol), falling back to a dlopen
//     (RTLD_NOLOAD) of the shared object this class itself lives in
//     (needed because Python loads extension modules RTLD_LOCAL).
//   - on Windows: LoadLibrary/GetProcAddress mirror the above, with
//     GetModuleHandle(NULL) and the module containing this class as the
//     no-library fallbacks.  UNTESTED on Windows (compile-guarded only).
//   Both "symbol" and "symbol_" (gfortran trailing underscore) are tried.
//
// Call contract (3D only): NTENS=6, NDI=3, NSHR=3, tension-positive with
// engineering shear strains.  The tangent is NOT symmetrized:
// hypoplastic/barodetic UMATs return unsymmetric tangents, so use a
// non-symmetric solver (e.g. system UmfPack or FullGeneral) with these
// materials.
//
// sendSelf/recvSelf serialize the construction recipe plus the committed
// state and re-resolve the UMAT symbol on the receiving process, so the
// material works across processes; the shared library (if any) must be
// present at the same path on every node.
//
// IMPLEMENTATION NOTES
//   - Voigt map: OpenSees order [11,22,33,12,23,31] vs Abaqus order
//     [11,22,33,12,13,23] differ only in components 5<->6; the swap (its
//     own inverse) is applied to strain in, stress out, and both indices
//     of DDSDDE.  DDSDDE arrives column-major (Fortran); no sign flips.
//   - State (FeapMaterial pattern): committed STATEV/stress/strain
//     snapshots plus trial copies.  Every setTrialStrain re-runs the UMAT
//     from the COMMITTED snapshot with the total strain increment since
//     the last commit; commitState copies trial->committed;
//     revertToLastCommit discards the trial state; revertToStart restores
//     the -statev / -initial_stress values.  getStress/getTangent never
//     run the UMAT (reference-returning accessors cannot report a
//     constitutive failure): they serve the last successful trial
//     evaluation, else the committed stress / the cached initial tangent.
//     Constitutive evaluation and its status live in setTrialStrain.
//   - PNEWDT: a returned PNEWDT < 1 is translated into a failed material
//     update (setTrialStrain returns -1).  The multiplier itself is NOT
//     propagated: step reduction and retry, where supported, remain the
//     responsibility of the analysis procedure driving the model.
//   - CELENT: the parent element's characteristic length, obtained via
//     the ops_TheActiveElement backdoor (ASDConcrete3D precedent) on the
//     first setTrialStrain and cached per instance (1.0 fallback).
//   - TIME/DTIME: DTIME is the domain's current time step (the ops_Dt
//     global): real dt under transient integrators, the pseudo-time step
//     (load-factor increment) under static ones, and 1.0 while unset
//     (before the first step).  TIME is the sum of the committed DTIMEs,
//     advanced on commitState so re-runs from the committed state replay
//     the same window.  Rate-independent models are unaffected;
//     rate-dependent UMATs receive the domain's (transient) time.  Both
//     are OpenSees translations, not reconstructed Abaqus step
//     bookkeeping: DTIME follows the domain increment and TIME advances
//     only on commit.
//   - CMNAME: OpenSees has no separate Abaqus material-name field, so
//     CMNAME carries the UMAT symbol as a deterministic identifier.
//   - Abaqus host metadata (COORDS, DROT, DFGRD0/1, NOEL, NPT, LAYER,
//     KSPT, JSTEP): deterministic placeholders only -- see invokeUmat.
//   - XIT: the Abaqus fatal-error utility is provided by UmatMaterial.cpp
//     as `xit_`; while a UMAT call is in flight a setjmp/longjmp trap
//     turns it into a -1 failure instead of terminating the process.
//     The trap is process-global state: UMAT evaluations are assumed not
//     to run concurrently from multiple threads.
//   - Serialization: one ID (tag, nprops, nstatv, string lengths) + one
//     Vector (rho, kinc, timeC, PROPS, initial STATEV/stress, committed
//     STATEV/stress/strain) + a second ID carrying the symbol and library
//     strings char-by-char (FileDatastore implements sendID/recvID but
//     not sendMsg/recvMsg, so the strings must not travel as a Message);
//     the dlsym pointer itself never travels.
//   - Windows: LoadLibrary/GetProcAddress branch mirrors the POSIX path;
//     compile-guarded only, UNTESTED on Windows.

#ifndef UmatMaterial_h
#define UmatMaterial_h

#include <NDMaterial.h>
#include <Vector.h>
#include <Matrix.h>
#include <stddef.h>

// Full Abaqus UMAT signature, gfortran calling convention: everything by
// reference, hidden CMNAME character-length argument appended at the end.
// KSTEP is passed as a 4-int array so that both scalar-KSTEP (2009-vintage)
// and JSTEP(4) (Abaqus >= 2017) UMATs read valid memory.
typedef void (*umat_fp)(
    double *stress, double *statev, double *ddsdde,
    double *sse, double *spd, double *scd,
    double *rpl, double *ddsddt, double *drplde, double *drpldt,
    double *stran, double *dstran, double *time, double *dtime,
    double *temp, double *dtemp, double *predef, double *dpred,
    char *cmname, int *ndi, int *nshr, int *ntens, int *nstatv,
    double *props, int *nprops, double *coords, double *drot,
    double *pnewdt, double *celent, double *dfgrd0, double *dfgrd1,
    int *noel, int *npt, int *layer, int *kspt, int *kstep, int *kinc,
    size_t cmname_len);

class UmatMaterial : public NDMaterial {
 public:
  UmatMaterial(int tag, const char *symbol, const char *library,
               int nprops, const double *props,
               int nstatv, const double *statev0,
               const double *initialStress,  // Abaqus order, may be 0
               double rho);
  UmatMaterial();  // for FEM_ObjectBroker (unusable until recvSelf)
  ~UmatMaterial();

  const char *getClassType(void) const { return "UmatMaterial"; }

  int setTrialStrain(const Vector &strain);
  const Vector &getStrain(void);
  const Vector &getStress(void);
  const Matrix &getTangent(void);
  const Matrix &getInitialTangent(void);
  double getRho(void);

  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);

  NDMaterial *getCopy(void);
  NDMaterial *getCopy(const char *type);
  const char *getType(void) const;
  int getOrder(void) const;

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel,
               FEM_ObjectBroker &theBroker);

  // runtime PROPS updates via the parameter API: "prop<i>" (1-based index
  // into PROPS).  Elements forward setParameter to their materials, so
  // per-element updates work: setParameter('-val', v, '-ele', tag, 'prop3').
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int responseID, Information &info);

  void Print(OPS_Stream &s, int flag = 0);

  // resolves the UMAT entry point; returns 0 on success
  int resolveSymbol(void);

  // enables the -check_tangent development diagnostic (h <= 0 disables)
  void setTangentCheck(double h) { checkTangentH = h; }

 private:
  // core UMAT call from the committed snapshot to strainTrial, into
  // caller-owned buffers (side-effect-free w.r.t. member state); returns
  // 0 on success, -1 on XIT / non-finite output / invalid or <1 PNEWDT
  int callUmat(const double *strainTrial, double *stressOut,
               double *statevOut, double *tangentOut);

  // runs the UMAT to the current trial strain; returns 0 on success, -1
  // on failure (trial stress/STATEV then restored to committed)
  int invokeUmat(void);

  // evaluates the UMAT once at the initial state and caches tangent0
  void computeInitialTangent(void);

  // -check_tangent diagnostic: central-difference DDSDDE verification
  void checkTangentFD(void);

  char *symbol;    // Fortran entry symbol (as typed by the user)
  char *library;   // optional shared-library path (0 if none)
  umat_fp umat;    // resolved entry point
  void *libHandle; // dlopen handle (0 when using the global/self scope)

  int nprops;
  int nstatv;
  double *props;    // material constants, PROPS (updatable via setParameter)
  double *statev0;  // initial STATEV (revertToStart target)
  double stress0[6]; // initial stress, Abaqus order (revertToStart target)

  // committed and trial state (all stress/strain in ABAQUS Voigt order)
  double *statevC;   // committed STATEV
  double *statevT;   // trial STATEV
  double stressC[6]; // committed stress
  double stressT[6]; // trial stress
  double strainC[6]; // committed strain
  double strainT[6]; // trial strain
  double tangentT[36]; // trial tangent (row-major, OpenSees order)
  double tangent0[36]; // cached initial-state tangent (getInitialTangent)

  double rho;
  int kinc;           // commit counter (Abaqus KINC, 1-based)
  bool tangentValid;  // tangentT/stressT reflect strainT
  bool tangent0Valid; // tangent0 has been computed

  // per-element characteristic length passed as CELENT, cached from
  // ops_TheActiveElement on the first setTrialStrain (1.0 fallback)
  double celent;
  bool celentCached;

  // domain time: dtimeT is the DTIME of the in-flight increment (the
  // ops_Dt global -- the domain's current time step -- or 1.0 while
  // unset); timeC is the committed TIME, the sum of the DTIMEs committed
  // via commitState
  double timeC;
  double dtimeT;

  // -check_tangent strain perturbation; 0 = diagnostic disabled
  double checkTangentH;

  static Vector sigma;    // 6, OpenSees order
  static Vector epsilon;  // 6, OpenSees order
  static Matrix tangent;  // 6x6, OpenSees order
};

#endif
