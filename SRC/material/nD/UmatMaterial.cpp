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
// Description: Generic Abaqus-UMAT wrapper NDMaterial.  See
// UmatMaterial.h for the command syntax, the Voigt mapping and the
// state semantics.

#include "UmatMaterial.h"

#include <Channel.h>
#include <Element.h>
#include <ID.h>
#include <Information.h>
#include <OPS_Globals.h>
#include <Parameter.h>
#include <classTags.h>
#include <elementAPI.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>

#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#include <csetjmp>

// ---------------------------------------------------------------------------
// Abaqus utility routine XIT.  User materials call XIT to terminate the
// analysis on fatal errors (bad parameters etc.); the host is expected to
// provide it.  gfortran mangles `call XIT` / `call xit()` to `xit_`.
// While a UMAT invocation is in flight we longjmp back into invokeUmat so
// the fatal error surfaces as a setTrialStrain failure (-1) instead of
// killing the hosting (e.g. Python) process; a call outside that window
// (should not happen) exits.
// The jump buffer and armed flag are process-global: UMAT evaluations are
// assumed not to run concurrently from multiple threads.
// ---------------------------------------------------------------------------
static jmp_buf umatXitJmpBuf;
static bool umatXitArmed = false;

extern "C" void xit_(void) {
  opserr << "UMAT called XIT (fatal error in user material)\n";
  if (umatXitArmed)
    longjmp(umatXitJmpBuf, 1);
  exit(-1);
}

// OpenSees Voigt order [11,22,33,12,23,31] <-> Abaqus [11,22,33,12,13,23]:
// components 5 and 6 swap (the map is its own inverse).
static const int VMAP[6] = {0, 1, 2, 3, 5, 4};

Vector UmatMaterial::sigma(6);
Vector UmatMaterial::epsilon(6);
Matrix UmatMaterial::tangent(6, 6);

static char *copyString(const char *s) {
  if (s == 0)
    return 0;
  char *c = new char[strlen(s) + 1];
  strcpy(c, s);
  return c;
}

// anchor whose address locates the shared object this code lives in
static void umatMaterialAnchor(void) {}

UmatMaterial::UmatMaterial(int tag, const char *sym, const char *lib,
                           int np, const double *p, int ns,
                           const double *sv0, const double *sig0, double r)
    : NDMaterial(tag, ND_TAG_UmatMaterial), symbol(copyString(sym)),
      library(copyString(lib)), umat(0), libHandle(0), nprops(np),
      nstatv(ns), props(0), statev0(0), statevC(0), statevT(0), rho(r),
      kinc(1), tangentValid(false), tangent0Valid(false), celent(1.0),
      celentCached(false), timeC(0.0), dtimeT(1.0), checkTangentH(0.0) {
  props = new double[nprops];
  for (int i = 0; i < nprops; i++)
    props[i] = p[i];

  // at least one slot so a stateless (NSTATV=0) UMAT still receives a
  // valid STATEV pointer
  int nsv = (nstatv > 0) ? nstatv : 1;
  statev0 = new double[nsv];
  statevC = new double[nsv];
  statevT = new double[nsv];
  statev0[0] = statevC[0] = statevT[0] = 0.0;
  for (int i = 0; i < nstatv; i++) {
    statev0[i] = (sv0 != 0) ? sv0[i] : 0.0;
    statevC[i] = statev0[i];
    statevT[i] = statev0[i];
  }

  for (int i = 0; i < 6; i++) {
    stress0[i] = (sig0 != 0) ? sig0[i] : 0.0;
    stressC[i] = stress0[i];
    stressT[i] = stress0[i];
    strainC[i] = 0.0;
    strainT[i] = 0.0;
  }
  for (int i = 0; i < 36; i++)
    tangentT[i] = tangent0[i] = 0.0;
}

UmatMaterial::UmatMaterial()
    : NDMaterial(0, ND_TAG_UmatMaterial), symbol(0), library(0), umat(0),
      libHandle(0), nprops(0), nstatv(0), props(0), statev0(0), statevC(0),
      statevT(0), rho(0.0), kinc(1), tangentValid(false),
      tangent0Valid(false), celent(1.0), celentCached(false), timeC(0.0),
      dtimeT(1.0), checkTangentH(0.0) {
  for (int i = 0; i < 6; i++) {
    stress0[i] = stressC[i] = stressT[i] = 0.0;
    strainC[i] = strainT[i] = 0.0;
  }
  for (int i = 0; i < 36; i++)
    tangentT[i] = tangent0[i] = 0.0;
}

UmatMaterial::~UmatMaterial() {
  if (symbol != 0)
    delete[] symbol;
  if (library != 0)
    delete[] library;
  if (props != 0)
    delete[] props;
  if (statev0 != 0)
    delete[] statev0;
  if (statevC != 0)
    delete[] statevC;
  if (statevT != 0)
    delete[] statevT;
  // libHandle is deliberately never dlclose'd / FreeLibrary'd: other
  // material instances (copies) may share the resolved function pointer.
}

#ifdef _WIN32
// LoadLibrary/GetProcAddress mirror of the POSIX path below.
// UNTESTED on Windows: compile-guarded only, no Windows build available.
static umat_fp findUmatSymbol(HMODULE handle, const char *name) {
  if (handle == 0)
    return 0;
  void *p = (void *)GetProcAddress(handle, name);
  if (p == 0) {
    // gfortran appends a trailing underscore
    char decorated[256];
    snprintf(decorated, sizeof(decorated), "%s_", name);
    p = (void *)GetProcAddress(handle, decorated);
  }
  return (umat_fp)p;
}
#else
static umat_fp findUmatSymbol(void *handle, const char *name) {
  void *p = dlsym(handle, name);
  if (p == 0) {
    // gfortran appends a trailing underscore
    char decorated[256];
    snprintf(decorated, sizeof(decorated), "%s_", name);
    p = dlsym(handle, decorated);
  }
  return (umat_fp)p;
}
#endif

int UmatMaterial::resolveSymbol(void) {
  if (symbol == 0) {
    opserr << "UmatMaterial::resolveSymbol - no -symbol given\n";
    return -1;
  }

#ifdef _WIN32
  if (library != 0) {
    HMODULE h = LoadLibraryA(library);
    if (h == 0) {
      opserr << "UmatMaterial::resolveSymbol - LoadLibrary(\"" << library
             << "\") failed (error " << (int)GetLastError() << ")\n";
      return -1;
    }
    libHandle = (void *)h;
    umat = findUmatSymbol(h, symbol);
    if (umat == 0) {
      opserr << "UmatMaterial::resolveSymbol - symbol \"" << symbol
             << "\" (also tried \"" << symbol << "_\") not found in library "
             << library << "\n";
      return -1;
    }
    return 0;
  }

  // no library: Windows has no global symbol scope; try the executable,
  // then the module this class itself was loaded from (compiled-in
  // Fortran via OPS_UMAT_SOURCES)
  umat = findUmatSymbol(GetModuleHandleA(0), symbol);
  if (umat == 0) {
    HMODULE self = 0;
    if (GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                               GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                           (LPCSTR)&umatMaterialAnchor, &self) != 0)
      umat = findUmatSymbol(self, symbol);
  }

  if (umat == 0) {
    opserr << "UmatMaterial::resolveSymbol - symbol \"" << symbol
           << "\" (also tried \"" << symbol << "_\") not found in the "
           << "executable or in the OpenSees module itself; compile the "
           << "Fortran source in via -DOPS_UMAT_SOURCES or pass "
           << "-library /path/to/lib.dll\n";
    return -1;
  }
  return 0;
#else
  if (library != 0) {
    libHandle = dlopen(library, RTLD_NOW | RTLD_GLOBAL);
    if (libHandle == 0) {
      opserr << "UmatMaterial::resolveSymbol - dlopen(\"" << library
             << "\") failed: " << dlerror() << "\n";
      return -1;
    }
    umat = findUmatSymbol(libHandle, symbol);
    if (umat == 0) {
      opserr << "UmatMaterial::resolveSymbol - symbol \"" << symbol
             << "\" (also tried \"" << symbol << "_\") not found in library "
             << library << "\n";
      return -1;
    }
    return 0;
  }

  // no library: search the global scope first ...
  umat = findUmatSymbol(RTLD_DEFAULT, symbol);

  // ... then the shared object this class itself was loaded from (Python
  // loads extension modules RTLD_LOCAL, so compiled-in Fortran symbols are
  // in our .so's dynamic table but not in the global scope).
  if (umat == 0) {
    Dl_info info;
    if (dladdr((void *)&umatMaterialAnchor, &info) != 0 &&
        info.dli_fname != 0) {
      void *self = dlopen(info.dli_fname, RTLD_LAZY | RTLD_NOLOAD);
      if (self != 0)
        umat = findUmatSymbol(self, symbol);
    }
  }

  if (umat == 0) {
    opserr << "UmatMaterial::resolveSymbol - symbol \"" << symbol
           << "\" (also tried \"" << symbol << "_\") not found in the "
           << "global scope or in the OpenSees module itself; compile the "
           << "Fortran source in via -DOPS_UMAT_SOURCES or pass "
           << "-library /path/to/lib.so\n";
    return -1;
  }
  return 0;
#endif
}

// Core UMAT call: runs the UMAT from the COMMITTED snapshot to the given
// trial strain (Abaqus order), writing into the caller-supplied buffers.
// stressOut[6]/statevOut[nstatv] are preloaded from the committed state
// here; tangentOut[36] receives DDSDDE mapped to OpenSees order,
// row-major.  No member state is touched beyond the buffers passed in,
// so scratch evaluations (initial tangent, finite-difference check) are
// side-effect-free.  Returns 0 on success, -1 on XIT / non-finite
// output / invalid or <1 PNEWDT.
int UmatMaterial::callUmat(const double *strainTrial, double *stressOut,
                           double *statevOut, double *tangentOut) {
  if (umat == 0) {
    opserr << "UmatMaterial::callUmat - UMAT symbol not resolved\n";
    return -1;
  }

  for (int i = 0; i < 6; i++)
    stressOut[i] = stressC[i];
  for (int i = 0; i < nstatv; i++)
    statevOut[i] = statevC[i];

  double stran[6], dstran[6];
  for (int i = 0; i < 6; i++) {
    stran[i] = strainC[i];
    dstran[i] = strainTrial[i] - strainC[i];
  }

  double ddsdde[36];
  for (int i = 0; i < 36; i++)
    ddsdde[i] = 0.0;

  // DTIME is the domain's current time step (dtimeT, set from the ops_Dt
  // global by invokeUmat): the real dt under a transient integrator, the
  // pseudo-time step (load-factor increment) under a static one --
  // OpenSees domain time IS pseudo-time there -- and 1.0 while
  // unset/non-positive (before the first step).  TIME is the sum of the
  // DTIMEs committed so far, so re-runs from the committed snapshot
  // replay the same time window.  Rate-independent models ignore both;
  // rate-dependent UMATs receive the domain's (transient) time.
  double dtime = dtimeT;
  double time[2];
  time[0] = time[1] = timeC;

  double sse = 0.0, spd = 0.0, scd = 0.0;
  double rpl = 0.0, drpldt = 0.0;
  double ddsddt[6] = {0, 0, 0, 0, 0, 0};
  double drplde[6] = {0, 0, 0, 0, 0, 0};
  double temp = 0.0, dtemp = 0.0;
  double predef[1] = {0.0}, dpred[1] = {0.0};

  // OpenSees has no separate Abaqus material-name field: CMNAME carries
  // the UMAT symbol as a deterministic identifier.  UMATs that dispatch
  // internally on a distinct Abaqus CMNAME need an adapter extension.
  char cmname[80];
  memset(cmname, ' ', sizeof(cmname));
  size_t symLen = strlen(symbol);
  if (symLen > sizeof(cmname))
    symLen = sizeof(cmname);
  memcpy(cmname, symbol, symLen);

  int ndi = 3, nshr = 3, ntens = 6;
  int nstatvL = nstatv, npropsL = nprops;

  // Abaqus host metadata with no NDMaterial equivalent: deterministic
  // placeholders, not a reproduction of the Abaqus element/integration-
  // point environment.  Supported UMATs must not depend on COORDS, DROT,
  // DFGRD0/1, NOEL, NPT, LAYER, KSPT or JSTEP.  Rotation and deformation-
  // gradient quantities are identity because the adapter targets the
  // small-strain NDMaterial interface.
  double coords[3] = {0.0, 0.0, 0.0};
  double drot[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1};   // identity
  double dfgrd0[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // identity
  double dfgrd1[9] = {1, 0, 0, 0, 1, 0, 0, 0, 1}; // identity
  double celentL = celent; // per-element characteristic length (see setTrialStrain)

  // never pass a null PNEWDT: the target models write it
  double pnewdt = 1.0;

  int noel = this->getTag();
  int npt = 1;
  int layer = 0, kspt = 0;
  int kstep[4] = {1, 0, 0, 0}; // covers scalar KSTEP and JSTEP(4) readers
  int kincL = kinc;

  // arm the XIT trap: a UMAT calling the Abaqus XIT utility longjmps back
  // here and the call is reported as a failure
  umatXitArmed = true;
  if (setjmp(umatXitJmpBuf) != 0) {
    umatXitArmed = false;
    opserr << "UmatMaterial::callUmat - UMAT \"" << symbol
           << "\" (tag " << this->getTag() << ") terminated via XIT\n";
    return -1;
  }

  umat(stressOut, statevOut, ddsdde, &sse, &spd, &scd, &rpl, ddsddt, drplde,
       &drpldt, stran, dstran, time, &dtime, &temp, &dtemp, predef, dpred,
       cmname, &ndi, &nshr, &ntens, &nstatvL, props, &npropsL, coords, drot,
       &pnewdt, &celentL, dfgrd0, dfgrd1, &noel, &npt, &layer, &kspt, kstep,
       &kincL, sizeof(cmname));

  umatXitArmed = false;

  // Fortran DDSDDE(I,J) = dsig_I/deps_J, column-major: element (I,J) is at
  // ddsdde[(I-1) + 6*(J-1)].  Map both indices through the 5<->6 swap.
  // NOT symmetrized (the target models have unsymmetric tangents).
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      tangentOut[6 * i + j] = ddsdde[VMAP[i] + 6 * VMAP[j]];

  // Validate all constitutive outputs before accepting the result: a
  // failed evaluation must never be served as a usable trial state.
  for (int i = 0; i < 6; i++)
    if (!std::isfinite(stressOut[i])) {
      opserr << "UmatMaterial::callUmat - non-finite stress returned by "
             << "UMAT \"" << symbol << "\" (tag " << this->getTag() << ")\n";
      return -1;
    }
  for (int i = 0; i < 36; i++)
    if (!std::isfinite(tangentOut[i])) {
      opserr << "UmatMaterial::callUmat - non-finite tangent returned by "
             << "UMAT \"" << symbol << "\" (tag " << this->getTag() << ")\n";
      return -1;
    }

  // a non-finite or non-positive PNEWDT is an invalid return, not a
  // cutback request (NaN < 1.0 is false, so it must be caught explicitly)
  if (!std::isfinite(pnewdt) || pnewdt <= 0.0) {
    opserr << "UmatMaterial::callUmat - invalid PNEWDT returned by UMAT \""
           << symbol << "\" (tag " << this->getTag() << ")\n";
    return -1;
  }

  // PNEWDT < 1 is the UMAT asking the host for a cutback
  if (pnewdt < 1.0)
    return -1;

  return 0;
}

// Runs the UMAT from the committed snapshot to the current trial strain.
// OpenSees calls setTrialStrain repeatedly while iterating toward
// equilibrium, so every evaluation restarts from the last committed
// constitutive state (Kratos/FeapMaterial semantics) -- rejected trial
// iterations must not accumulate irreversible STATEV or stress updates.
// On failure the trial stress/STATEV are restored to the committed
// values, so the accessors always serve a defined state.
int UmatMaterial::invokeUmat(void) {
  // the previous trial result is stale from here on; only a fully
  // validated evaluation may mark the tangent usable again
  tangentValid = false;

  dtimeT = (ops_Dt > 0.0) ? ops_Dt : 1.0; // see the TIME note in callUmat

  if (this->callUmat(strainT, stressT, statevT, tangentT) != 0) {
    for (int i = 0; i < 6; i++)
      stressT[i] = stressC[i];
    for (int i = 0; i < nstatv; i++)
      statevT[i] = statevC[i];
    return -1;
  }

  if (checkTangentH > 0.0)
    this->checkTangentFD();

  tangentValid = true;
  return 0;
}

// Optional development diagnostic (-check_tangent h): central-difference
// verification of the UMAT Jacobian,
//   dsig_i/deps_j ~ (sig_i(eps + h e_j) - sig_i(eps - h e_j)) / (2 h),
// where h is the user-supplied strain perturbation, applied verbatim (no
// silent rescaling).  Every +h/-h probe runs from the exact same
// committed stress/STATEV snapshot as the trial evaluation itself
// (callUmat re-preloads them; scratch buffers keep the trial state
// untouched).  DIAGNOSTIC ONLY: the FD tangent is never substituted for
// DDSDDE and a mismatch never fails the analysis -- it reports the
// worst absolute discrepancy and its relative size
//   e_rel = max|C_UMAT - C_FD| / max(max|C_UMAT|, max|C_FD|, eps)
// when e_rel exceeds 10% (FD across a path-dependent model is itself
// approximate: treat the report as a smell, not a proof).
void UmatMaterial::checkTangentFD(void) {
  int nsv = (nstatv > 0) ? nstatv : 1;
  double *sv = new double[nsv];
  double sp[6], sm[6], tanScratch[36], strainP[6], fd[36];
  double h = checkTangentH;

  bool complete = true;
  for (int j = 0; j < 6 && complete; j++) {
    for (int i = 0; i < 6; i++)
      strainP[i] = strainT[i];

    strainP[j] = strainT[j] + h;
    if (this->callUmat(strainP, sp, sv, tanScratch) != 0)
      complete = false;
    strainP[j] = strainT[j] - h;
    if (complete && this->callUmat(strainP, sm, sv, tanScratch) != 0)
      complete = false;

    for (int i = 0; complete && i < 6; i++)
      fd[6 * i + j] = (sp[i] - sm[i]) / (2.0 * h); // Abaqus order (i,j)
  }
  delete[] sv;

  if (!complete) {
    opserr << "UmatMaterial::checkTangentFD - tag " << this->getTag()
           << ": a perturbed evaluation failed, check skipped\n";
    return;
  }

  // compare in OpenSees order against the accepted DDSDDE
  double ref = 1.0e-30, worst = 0.0;
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++) {
      double t = tangentT[6 * i + j];
      double f = fd[6 * VMAP[i] + VMAP[j]];
      double d = t - f;
      if (t < 0.0)
        t = -t;
      if (f < 0.0)
        f = -f;
      if (d < 0.0)
        d = -d;
      if (t > ref)
        ref = t;
      if (f > ref)
        ref = f;
      if (d > worst)
        worst = d;
    }
  if (worst > 0.10 * ref)
    opserr << "UmatMaterial::checkTangentFD - tag " << this->getTag()
           << ": max |DDSDDE - FD| = " << worst << ", relative "
           << 100.0 * worst / ref << "% (h = " << h << ")\n";
}

int UmatMaterial::setTrialStrain(const Vector &strain) {
  // Cache the parent element's characteristic length as CELENT on the
  // first call, via the ops_TheActiveElement backdoor (same mechanism as
  // ASDConcrete3DMaterial::setTrialStrain).  Materials are per-Gauss-point
  // copies, so caching per instance is per-element-correct and keeps
  // CELENT stable across subsequent constitutive iterations.
  if (!celentCached) {
    if (ops_TheActiveElement != 0)
      celent = ops_TheActiveElement->getCharacteristicLength();
    celentCached = true;
  }

  // OpenSees -> Abaqus component order
  for (int i = 0; i < 6; i++)
    strainT[VMAP[i]] = strain(i);
  return this->invokeUmat();
}

const Vector &UmatMaterial::getStrain(void) {
  for (int i = 0; i < 6; i++)
    epsilon(i) = strainT[VMAP[i]];
  return epsilon;
}

// The accessors never run the UMAT: reference-returning signatures cannot
// propagate a constitutive failure (PNEWDT, XIT, non-finite output), so
// evaluation and its status reporting live in setTrialStrain/invokeUmat.
// stressT always holds a defined state -- the last successful trial
// evaluation, else the committed/initial stress (construction, revert,
// recvSelf and failed evaluations all leave it there).
const Vector &UmatMaterial::getStress(void) {
  for (int i = 0; i < 6; i++)
    sigma(i) = stressT[VMAP[i]];
  return sigma;
}

const Matrix &UmatMaterial::getTangent(void) {
  // OpenSees forms the first tangent BEFORE any setTrialStrain call (and
  // revert/recvSelf/updateParameter also clear the trial tangent), so
  // when no valid trial evaluation exists serve the cached initial
  // tangent instead of an all-zero matrix; its dedicated evaluation
  // checks its own result
  if (!tangentValid)
    return this->getInitialTangent();

  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      tangent(i, j) = tangentT[6 * i + j];
  return tangent;
}

// Evaluate the UMAT once at the INITIAL state (initial STATEV/stress, zero
// strain, zero increment) and cache DDSDDE as the initial tangent.  The
// live committed/trial state is stashed and restored around the call, so
// the query leaves the constitutive history untouched.
void UmatMaterial::computeInitialTangent(void) {
  if (umat == 0)
    return; // symbol not resolved yet (broker-constructed): keep zeros

  // stash the live state
  int nsv = (nstatv > 0) ? nstatv : 1;
  double *sv = new double[2 * nsv];
  memcpy(sv, statevC, nsv * sizeof(double));
  memcpy(sv + nsv, statevT, nsv * sizeof(double));
  double sC[6], sT[6], eC[6], eT[6], tT[36];
  memcpy(sC, stressC, sizeof(sC));
  memcpy(sT, stressT, sizeof(sT));
  memcpy(eC, strainC, sizeof(eC));
  memcpy(eT, strainT, sizeof(eT));
  memcpy(tT, tangentT, sizeof(tT));
  bool tv = tangentValid;
  int k = kinc;
  double tC = timeC, dtT = dtimeT;

  // point the material at its initial state and run the regular UMAT
  // path (XIT trap, Voigt mapping, output validation)
  for (int i = 0; i < nstatv; i++)
    statevC[i] = statev0[i];
  for (int i = 0; i < 6; i++) {
    stressC[i] = stress0[i];
    strainC[i] = 0.0;
    strainT[i] = 0.0;
  }
  kinc = 1;
  timeC = 0.0;

  if (this->invokeUmat() == 0) {
    memcpy(tangent0, tangentT, sizeof(tangent0));
    tangent0Valid = true;
  } else {
    opserr << "UmatMaterial::getInitialTangent - UMAT \"" << symbol
           << "\" (tag " << this->getTag() << ") failed at the initial "
           << "state; falling back to the current tangent\n";
  }

  // restore the live state
  memcpy(statevC, sv, nsv * sizeof(double));
  memcpy(statevT, sv + nsv, nsv * sizeof(double));
  delete[] sv;
  memcpy(stressC, sC, sizeof(sC));
  memcpy(stressT, sT, sizeof(sT));
  memcpy(strainC, eC, sizeof(eC));
  memcpy(strainT, eT, sizeof(eT));
  memcpy(tangentT, tT, sizeof(tT));
  tangentValid = tv;
  kinc = k;
  timeC = tC;
  dtimeT = dtT;
}

// Returns the tangent of the INITIAL material state, evaluated once and
// cached: independent of the subsequent trial/committed history, for
// algorithms that explicitly request the initial stiffness.
const Matrix &UmatMaterial::getInitialTangent(void) {
  if (!tangent0Valid)
    this->computeInitialTangent();
  if (!tangent0Valid) {
    // intentional fallback (warned in computeInitialTangent): the last
    // trial tangent stands in because no true initial stiffness is
    // available.  Served directly from tangentT -- NOT via getTangent(),
    // which falls back here when no trial tangent exists.
    for (int i = 0; i < 6; i++)
      for (int j = 0; j < 6; j++)
        tangent(i, j) = tangentT[6 * i + j];
    return tangent;
  }

  for (int i = 0; i < 6; i++)
    for (int j = 0; j < 6; j++)
      tangent(i, j) = tangent0[6 * i + j];
  return tangent;
}

double UmatMaterial::getRho(void) { return rho; }

int UmatMaterial::commitState(void) {
  for (int i = 0; i < 6; i++) {
    stressC[i] = stressT[i];
    strainC[i] = strainT[i];
  }
  for (int i = 0; i < nstatv; i++)
    statevC[i] = statevT[i];
  timeC += dtimeT; // commit the time window actually used (see invokeUmat)
  kinc++;
  return 0;
}

int UmatMaterial::revertToLastCommit(void) {
  for (int i = 0; i < 6; i++) {
    stressT[i] = stressC[i];
    strainT[i] = strainC[i];
  }
  for (int i = 0; i < nstatv; i++)
    statevT[i] = statevC[i];
  tangentValid = false;
  return 0;
}

int UmatMaterial::revertToStart(void) {
  for (int i = 0; i < 6; i++) {
    stressC[i] = stressT[i] = stress0[i];
    strainC[i] = strainT[i] = 0.0;
  }
  for (int i = 0; i < nstatv; i++)
    statevC[i] = statevT[i] = statev0[i];
  kinc = 1;
  timeC = 0.0;
  dtimeT = 1.0;
  tangentValid = false;
  return 0;
}

// A true material clone, including the current committed and trial
// constitutive state.  OpenSees may request copies after the material has
// evolved, so getCopy() must not reset the copy to its initial STATEV,
// stress, strain, time or tangent state.
NDMaterial *UmatMaterial::getCopy(void) {
  UmatMaterial *copy =
      new UmatMaterial(this->getTag(), symbol, library, nprops, props,
                       nstatv, statev0, stress0, rho);
  // share the resolved entry point
  copy->umat = umat;
  copy->libHandle = 0; // owner semantics stay with the original

  for (int i = 0; i < nstatv; i++) {
    copy->statevC[i] = statevC[i];
    copy->statevT[i] = statevT[i];
  }
  for (int i = 0; i < 6; i++) {
    copy->stressC[i] = stressC[i];
    copy->stressT[i] = stressT[i];
    copy->strainC[i] = strainC[i];
    copy->strainT[i] = strainT[i];
  }
  for (int i = 0; i < 36; i++) {
    copy->tangentT[i] = tangentT[i];
    copy->tangent0[i] = tangent0[i];
  }
  copy->tangentValid = tangentValid;
  copy->tangent0Valid = tangent0Valid;
  copy->kinc = kinc;
  copy->timeC = timeC;
  copy->dtimeT = dtimeT;
  copy->checkTangentH = checkTangentH;
  // CELENT stays uncached: each copy belongs to its own element's
  // integration point and must pick up that element's length
  return copy;
}

NDMaterial *UmatMaterial::getCopy(const char *type) {
  if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0)
    return this->getCopy();

  opserr << "UmatMaterial::getCopy - type \"" << type
         << "\" not supported (3D only: the wrapped UMATs require NDI=3)\n";
  return 0;
}

const char *UmatMaterial::getType(void) const { return "ThreeDimensional"; }

int UmatMaterial::getOrder(void) const { return 6; }

// runtime PROPS updates: "prop<i>" addresses PROPS(i), 1-based.  Elements
// forward setParameter to their materials, so per-element updates work:
// setParameter('-val', v, '-ele', tag, 'prop3').
int UmatMaterial::setParameter(const char **argv, int argc,
                               Parameter &param) {
  if (argc < 1 || argv[0] == 0)
    return -1;
  if (strncmp(argv[0], "prop", 4) == 0) {
    int idx = atoi(argv[0] + 4);
    if (idx >= 1 && idx <= nprops)
      return param.addObject(1000 + idx, this);
  }
  return -1;
}

int UmatMaterial::updateParameter(int responseID, Information &info) {
  if (responseID > 1000 && responseID <= 1000 + nprops) {
    props[responseID - 1001] = info.theDouble;
    // both cached tangents were computed with the old PROPS: the next
    // setTrialStrain reruns the UMAT from the committed snapshot, and
    // getInitialTangent re-evaluates the initial state on demand
    tangentValid = false;
    tangent0Valid = false;
    return 0;
  }
  return -1;
}

// What travels across the Channel is the construction RECIPE (tag, sizes,
// PROPS, initial STATEV/stress, rho, the symbol and library strings) plus
// the committed state -- never the dlsym function pointer.  Three sends:
// an ID (tag, nprops, nstatv, string lengths), a Vector (numerics), and a
// second ID carrying the strings char-by-char (FileDatastore implements
// sendID/recvID but not sendMsg/recvMsg, so a Message leg would break
// op.database()/op.save() for every UMAT-bearing domain).  recvSelf
// re-resolves the symbol through the same dlopen/dlsym path the
// constructor uses, so the shared library (if any) must be present at the
// same path on every process.
int UmatMaterial::sendSelf(int commitTag, Channel &theChannel) {
  int dbTag = this->getDbTag();

  int symLen = (symbol != 0) ? (int)strlen(symbol) : 0;
  int libLen = (library != 0) ? (int)strlen(library) : 0;

  ID idata(5);
  idata(0) = this->getTag();
  idata(1) = nprops;
  idata(2) = nstatv;
  idata(3) = symLen;
  idata(4) = libLen;
  if (theChannel.sendID(dbTag, commitTag, idata) < 0) {
    opserr << "UmatMaterial::sendSelf - failed to send ID data\n";
    return -1;
  }

  // rho, kinc, timeC, checkTangentH, PROPS, initial STATEV, initial
  // stress, committed STATEV, committed stress, committed strain
  Vector vdata(4 + nprops + 2 * nstatv + 18);
  int pos = 0;
  vdata(pos++) = rho;
  vdata(pos++) = (double)kinc;
  vdata(pos++) = timeC;
  vdata(pos++) = checkTangentH;
  for (int i = 0; i < nprops; i++)
    vdata(pos++) = props[i];
  for (int i = 0; i < nstatv; i++)
    vdata(pos++) = statev0[i];
  for (int i = 0; i < 6; i++)
    vdata(pos++) = stress0[i];
  for (int i = 0; i < nstatv; i++)
    vdata(pos++) = statevC[i];
  for (int i = 0; i < 6; i++)
    vdata(pos++) = stressC[i];
  for (int i = 0; i < 6; i++)
    vdata(pos++) = strainC[i];
  if (theChannel.sendVector(dbTag, commitTag, vdata) < 0) {
    opserr << "UmatMaterial::sendSelf - failed to send Vector data\n";
    return -1;
  }

  // symbol and library strings packed into an ID, char -> int (lengths
  // travel in the first ID; see the note above on FileDatastore)
  if (symLen + libLen > 0) {
    ID sid(symLen + libLen);
    for (int i = 0; i < symLen; i++)
      sid(i) = (int)(unsigned char)symbol[i];
    for (int i = 0; i < libLen; i++)
      sid(symLen + i) = (int)(unsigned char)library[i];
    if (theChannel.sendID(dbTag, commitTag, sid) < 0) {
      opserr << "UmatMaterial::sendSelf - failed to send string data\n";
      return -1;
    }
  }

  return 0;
}

int UmatMaterial::recvSelf(int commitTag, Channel &theChannel,
                           FEM_ObjectBroker &theBroker) {
  int dbTag = this->getDbTag();

  ID idata(5);
  if (theChannel.recvID(dbTag, commitTag, idata) < 0) {
    opserr << "UmatMaterial::recvSelf - failed to recv ID data\n";
    return -1;
  }
  this->setTag(idata(0));
  nprops = idata(1);
  nstatv = idata(2);
  int symLen = idata(3);
  int libLen = idata(4);

  // (re)allocate the owned arrays to the received sizes
  if (props != 0)
    delete[] props;
  if (statev0 != 0)
    delete[] statev0;
  if (statevC != 0)
    delete[] statevC;
  if (statevT != 0)
    delete[] statevT;
  int nsv = (nstatv > 0) ? nstatv : 1; // valid pointer even when NSTATV=0
  props = new double[nprops];
  statev0 = new double[nsv];
  statevC = new double[nsv];
  statevT = new double[nsv];
  statev0[0] = statevC[0] = statevT[0] = 0.0;

  Vector vdata(4 + nprops + 2 * nstatv + 18);
  if (theChannel.recvVector(dbTag, commitTag, vdata) < 0) {
    opserr << "UmatMaterial::recvSelf - failed to recv Vector data\n";
    return -1;
  }
  int pos = 0;
  rho = vdata(pos++);
  kinc = (int)vdata(pos++);
  timeC = vdata(pos++);
  checkTangentH = vdata(pos++);
  for (int i = 0; i < nprops; i++)
    props[i] = vdata(pos++);
  for (int i = 0; i < nstatv; i++)
    statev0[i] = vdata(pos++);
  for (int i = 0; i < 6; i++)
    stress0[i] = vdata(pos++);
  for (int i = 0; i < nstatv; i++)
    statevC[i] = vdata(pos++);
  for (int i = 0; i < 6; i++)
    stressC[i] = vdata(pos++);
  for (int i = 0; i < 6; i++)
    strainC[i] = vdata(pos++);

  // trial state starts at the received committed state
  for (int i = 0; i < nstatv; i++)
    statevT[i] = statevC[i];
  for (int i = 0; i < 6; i++) {
    stressT[i] = stressC[i];
    strainT[i] = strainC[i];
  }
  dtimeT = 1.0;
  tangentValid = false;
  tangent0Valid = false; // re-evaluate on the receiving process
  celent = 1.0;
  celentCached = false; // re-cache from the receiving process's elements

  char *sdata = new char[symLen + libLen + 1];
  if (symLen + libLen > 0) {
    // strings travel as an ID (see sendSelf: FileDatastore lacks recvMsg)
    ID sid(symLen + libLen);
    if (theChannel.recvID(dbTag, commitTag, sid) < 0) {
      opserr << "UmatMaterial::recvSelf - failed to recv string data\n";
      delete[] sdata;
      return -1;
    }
    for (int i = 0; i < symLen + libLen; i++)
      sdata[i] = (char)sid(i);
  }
  sdata[symLen + libLen] = '\0';
  if (symbol != 0)
    delete[] symbol;
  symbol = new char[symLen + 1];
  memcpy(symbol, sdata, symLen);
  symbol[symLen] = '\0';
  if (library != 0)
    delete[] library;
  library = 0;
  if (libLen > 0) {
    library = new char[libLen + 1];
    memcpy(library, sdata + symLen, libLen);
    library[libLen] = '\0';
  }
  delete[] sdata;

  // re-resolve on this process (the pointer itself never travels)
  umat = 0;
  libHandle = 0;
  return this->resolveSymbol();
}

void UmatMaterial::Print(OPS_Stream &s, int flag) {
  s << "UmatMaterial, tag: " << this->getTag() << "\n";
  s << "  symbol: " << (symbol != 0 ? symbol : "(none)") << "\n";
  s << "  library: " << (library != 0 ? library : "(compiled in)") << "\n";
  s << "  nprops: " << nprops << ", nstatv: " << nstatv << "\n";
  s << "  rho: " << rho << "\n";
}

//
// nDMaterial UMAT $tag -symbol <sym> -nprops N -props p1..pN -nstatv M
//                      <-statev v1..vM> <-initial_stress s11 s22 s33 s12 s13 s23>
//                      <-library /path/lib.so> <-rho r> <-check_tangent h>
//
void *OPS_UmatMaterial(void) {
  if (OPS_GetNumRemainingInputArgs() < 7) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: nDMaterial UMAT tag? -symbol sym? "
           << "-nprops N? -props p1..pN? -nstatv M? <-statev v1..vM> "
           << "<-initial_stress s11 s22 s33 s12 s13 s23> "
           << "<-library path?> <-rho r?> <-check_tangent h?>\n";
    return 0;
  }

  int tag = 0;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) {
    opserr << "WARNING invalid integer tag: nDMaterial UMAT\n";
    return 0;
  }

  char *symbol = 0;
  char *library = 0;
  int nprops = 0, nstatv = 0;
  double *props = 0;
  double *statev = 0;
  double initialStress[6] = {0, 0, 0, 0, 0, 0};
  bool haveInitialStress = false;
  bool haveNstatv = false;
  double rho = 0.0;
  double checkTangentH = 0.0;
  bool ok = true;

  while (ok && OPS_GetNumRemainingInputArgs() > 0) {
    const char *opt = OPS_GetString();

    if (strcmp(opt, "-symbol") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING -symbol requires a value: nDMaterial UMAT\n";
        ok = false;
        break;
      }
      symbol = copyString(OPS_GetString());

    } else if (strcmp(opt, "-library") == 0) {
      if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "WARNING -library requires a value: nDMaterial UMAT\n";
        ok = false;
        break;
      }
      library = copyString(OPS_GetString());

    } else if (strcmp(opt, "-nprops") == 0) {
      numData = 1;
      if (OPS_GetIntInput(&numData, &nprops) < 0 || nprops <= 0) {
        opserr << "WARNING invalid -nprops: nDMaterial UMAT\n";
        ok = false;
        break;
      }

    } else if (strcmp(opt, "-props") == 0) {
      if (nprops <= 0) {
        opserr << "WARNING -nprops N must come before -props: nDMaterial UMAT\n";
        ok = false;
        break;
      }
      props = new double[nprops];
      numData = nprops;
      if (OPS_GetDoubleInput(&numData, props) < 0) {
        opserr << "WARNING expected " << nprops
               << " values after -props: nDMaterial UMAT\n";
        ok = false;
        break;
      }

    } else if (strcmp(opt, "-nstatv") == 0) {
      numData = 1;
      // zero is allowed: stateless UMATs need no solution-dependent
      // state variables (-statev is then omitted)
      if (OPS_GetIntInput(&numData, &nstatv) < 0 || nstatv < 0) {
        opserr << "WARNING invalid -nstatv: nDMaterial UMAT\n";
        ok = false;
        break;
      }
      haveNstatv = true;

    } else if (strcmp(opt, "-statev") == 0) {
      if (nstatv <= 0) {
        opserr << "WARNING -nstatv M must come before -statev: nDMaterial UMAT\n";
        ok = false;
        break;
      }
      statev = new double[nstatv];
      numData = nstatv;
      if (OPS_GetDoubleInput(&numData, statev) < 0) {
        opserr << "WARNING expected " << nstatv
               << " values after -statev: nDMaterial UMAT\n";
        ok = false;
        break;
      }

    } else if (strcmp(opt, "-initial_stress") == 0) {
      numData = 6;
      if (OPS_GetDoubleInput(&numData, initialStress) < 0) {
        opserr << "WARNING expected 6 values (Abaqus order "
               << "s11 s22 s33 s12 s13 s23) after -initial_stress: "
               << "nDMaterial UMAT\n";
        ok = false;
        break;
      }
      haveInitialStress = true;

    } else if (strcmp(opt, "-rho") == 0) {
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &rho) < 0) {
        opserr << "WARNING invalid -rho: nDMaterial UMAT\n";
        ok = false;
        break;
      }

    } else if (strcmp(opt, "-check_tangent") == 0) {
      // development diagnostic: FD-verify DDSDDE with strain
      // perturbation h (see UmatMaterial::checkTangentFD)
      numData = 1;
      if (OPS_GetDoubleInput(&numData, &checkTangentH) < 0 ||
          checkTangentH <= 0.0) {
        opserr << "WARNING invalid -check_tangent (want h > 0): "
               << "nDMaterial UMAT\n";
        ok = false;
        break;
      }

    } else {
      opserr << "WARNING unknown option " << opt << ": nDMaterial UMAT\n";
      ok = false;
      break;
    }
  }

  if (ok && symbol == 0) {
    opserr << "WARNING -symbol is required: nDMaterial UMAT\n";
    ok = false;
  }
  if (ok && (nprops <= 0 || props == 0)) {
    opserr << "WARNING -nprops and -props are required: nDMaterial UMAT\n";
    ok = false;
  }
  if (ok && !haveNstatv) {
    opserr << "WARNING -nstatv is required: nDMaterial UMAT\n";
    ok = false;
  }

  UmatMaterial *theMaterial = 0;
  if (ok) {
    theMaterial = new UmatMaterial(
        tag, symbol, library, nprops, props, nstatv, statev,
        haveInitialStress ? initialStress : 0, rho);
    theMaterial->setTangentCheck(checkTangentH);
    if (theMaterial->resolveSymbol() != 0) {
      delete theMaterial;
      theMaterial = 0;
    }
  }

  if (symbol != 0)
    delete[] symbol;
  if (library != 0)
    delete[] library;
  if (props != 0)
    delete[] props;
  if (statev != 0)
    delete[] statev;

  return theMaterial;
}
