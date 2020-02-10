#ifndef __ARPACK_H__
#define __ARPACK_H__

#include "arpackdef.h"

#ifdef __cplusplus
extern "C" {
#endif

void cnaupd_c(a_int* ido, char const* bmat, a_int n, char const* which, a_int nev, float tol, float _Complex* resid, a_int ncv, float _Complex* v, a_int ldv, a_int* iparam, a_int* ipntr, float _Complex* workd, float _Complex* workl, a_int lworkl, float* rwork, a_int* info);
void cneupd_c(bool rvec, char const* howmny, a_int const* select, float _Complex* d, float _Complex* z, a_int ldz, float _Complex sigma, float _Complex* workev, char const* bmat, a_int n, char const* which, a_int nev, float tol, float _Complex* resid, a_int ncv, float _Complex* v, a_int ldv, a_int* iparam, a_int* ipntr, float _Complex* workd, float _Complex* workl, a_int lworkl, float* rwork, a_int* info);
void dnaupd_c(a_int* ido, char const* bmat, a_int n, char const* which, a_int nev, double tol, double* resid, a_int ncv, double* v, a_int ldv, a_int* iparam, a_int* ipntr, double* workd, double* workl, a_int lworkl, a_int* info);
void dneupd_c(bool rvec, char const* howmny, a_int const* select, double* dr, double* di, double* z, a_int ldz, double sigmar, double sigmai, double * workev, char const* bmat, a_int n, char const* which, a_int nev, double tol, double* resid, a_int ncv, double* v, a_int ldv, a_int* iparam, a_int* ipntr, double* workd, double* workl, a_int lworkl, a_int* info);
void dsaupd_c(a_int* ido, char const* bmat, a_int n, char const* which, a_int nev, double tol, double* resid, a_int ncv, double* v, a_int ldv, a_int* iparam, a_int* ipntr, double* workd, double* workl, a_int lworkl, a_int* info);
void dseupd_c(bool rvec, char const* howmny, a_int const* select, double* d, double* z, a_int ldz, double sigma, char const* bmat, a_int n, char const* which, a_int nev, double tol, double* resid, a_int ncv, double* v, a_int ldv, a_int* iparam, a_int* ipntr, double* workd, double* workl, a_int lworkl, a_int* info);
void snaupd_c(a_int* ido, char const* bmat, a_int n, char const* which, a_int nev, float tol, float* resid, a_int ncv, float* v, a_int ldv, a_int* iparam, a_int* ipntr, float* workd, float* workl, a_int lworkl, a_int* info);
void sneupd_c(bool rvec, char const* howmny, a_int const* select, float* dr, float* di, float* z, a_int ldz, float sigmar, float sigmai, float * workev, char const* bmat, a_int n, char const* which, a_int nev, float tol, float* resid, a_int ncv, float* v, a_int ldv, a_int* iparam, a_int* ipntr, float* workd, float* workl, a_int lworkl, a_int* info);
void ssaupd_c(a_int* ido, char const* bmat, a_int n, char const* which, a_int nev, float tol, float* resid, a_int ncv, float* v, a_int ldv, a_int* iparam, a_int* ipntr, float* workd, float* workl, a_int lworkl, a_int* info);
void sseupd_c(bool rvec, char const* howmny, a_int const* select, float* d, float* z, a_int ldz, float sigma, char const* bmat, a_int n, char const* which, a_int nev, float tol, float* resid, a_int ncv, float* v, a_int ldv, a_int* iparam, a_int* ipntr, float* workd, float* workl, a_int lworkl, a_int* info);
void znaupd_c(a_int* ido, char const* bmat, a_int n, char const* which, a_int nev, double tol, double _Complex* resid, a_int ncv, double _Complex* v, a_int ldv, a_int* iparam, a_int* ipntr, double _Complex* workd, double _Complex* workl, a_int lworkl, double* rwork, a_int* info);
void zneupd_c(bool rvec, char const* howmny, a_int const* select, double _Complex* d, double _Complex* z, a_int ldz, double _Complex sigma, double _Complex* workev, char const* bmat, a_int n, char const* which, a_int nev, double tol, double _Complex* resid, a_int ncv, double _Complex* v, a_int ldv, a_int* iparam, a_int* ipntr, double _Complex* workd, double _Complex* workl, a_int lworkl, double* rwork, a_int* info);

#ifdef  __cplusplus
}
#endif

#endif
