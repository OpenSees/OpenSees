#ifndef __PARPACK_H__
#define __PARPACK_H__

#include "arpackdef.h"

/*
 * IMPORTANT: MPI communicators MUST be passed from C to Fortran using MPI_Comm_c2f.
 *            MPI_Fint MCW = MPI_Comm_c2f(MPI_COMM_WORLD);
 */

#include "mpi.h"

#ifdef __cplusplus
extern "C" {
#endif

void pcnaupd_c(MPI_Fint comm, a_int* ido, char const* bmat, a_int n, char const* which, a_int nev, float tol, float _Complex* resid, a_int ncv, float _Complex* v, a_int ldv, a_int* iparam, a_int* ipntr, float _Complex* workd, float _Complex* workl, a_int lworkl, float _Complex* rwork, a_int* info);
void pcneupd_c(MPI_Fint comm, bool rvec, char const* howmny, a_int const* select, float _Complex* d, float _Complex* z, a_int ldz, float _Complex sigma, float _Complex* workev, char const* bmat, a_int n, char const* which, a_int nev, float tol, float _Complex* resid, a_int ncv, float _Complex* v, a_int ldv, a_int* iparam, a_int* ipntr, float _Complex* workd, float _Complex* workl, a_int lworkl, float _Complex* rwork, a_int* info);
void pdnaupd_c(MPI_Fint comm, a_int* ido, char const* bmat, a_int n, char const* which, a_int nev, double tol, double* resid, a_int ncv, double* v, a_int ldv, a_int* iparam, a_int* ipntr, double* workd, double* workl, a_int lworkl, a_int* info);
void pdneupd_c(MPI_Fint comm, bool rvec, char const* howmny, a_int const* select, double* dr, double* di, double* z, a_int ldz, double sigmar, double sigmai, double * workev, char const* bmat, a_int n, char const* which, a_int nev, double tol, double* resid, a_int ncv, double* v, a_int ldv, a_int* iparam, a_int* ipntr, double* workd, double* workl, a_int lworkl, a_int* info);
void pdsaupd_c(MPI_Fint comm, a_int* ido, char const* bmat, a_int n, char const* which, a_int nev, double tol, double* resid, a_int ncv, double* v, a_int ldv, a_int* iparam, a_int* ipntr, double* workd, double* workl, a_int lworkl, a_int* info);
void pdseupd_c(MPI_Fint comm, bool rvec, char const* howmny, a_int const* select, double* d, double* z, a_int ldz, double sigma, char const* bmat, a_int n, char const* which, a_int nev, double tol, double* resid, a_int ncv, double* v, a_int ldv, a_int* iparam, a_int* ipntr, double* workd, double* workl, a_int lworkl, a_int* info);
void psnaupd_c(MPI_Fint comm, a_int* ido, char const* bmat, a_int n, char const* which, a_int nev, float tol, float* resid, a_int ncv, float* v, a_int ldv, a_int* iparam, a_int* ipntr, float* workd, float* workl, a_int lworkl, a_int* info);
void psneupd_c(MPI_Fint comm, bool rvec, char const* howmny, a_int const* select, float* dr, float* di, float* z, a_int ldz, float sigmar, float sigmai, float * workev, char const* bmat, a_int n, char const* which, a_int nev, float tol, float* resid, a_int ncv, float* v, a_int ldv, a_int* iparam, a_int* ipntr, float* workd, float* workl, a_int lworkl, a_int* info);
void pssaupd_c(MPI_Fint comm, a_int* ido, char const* bmat, a_int n, char const* which, a_int nev, float tol, float* resid, a_int ncv, float* v, a_int ldv, a_int* iparam, a_int* ipntr, float* workd, float* workl, a_int lworkl, a_int* info);
void psseupd_c(MPI_Fint comm, bool rvec, char const* howmny, a_int const* select, float* d, float* z, a_int ldz, float sigma, char const* bmat, a_int n, char const* which, a_int nev, float tol, float* resid, a_int ncv, float* v, a_int ldv, a_int* iparam, a_int* ipntr, float* workd, float* workl, a_int lworkl, a_int* info);
void pznaupd_c(MPI_Fint comm, a_int* ido, char const* bmat, a_int n, char const* which, a_int nev, double tol, double _Complex* resid, a_int ncv, double _Complex* v, a_int ldv, a_int* iparam, a_int* ipntr, double _Complex* workd, double _Complex* workl, a_int lworkl, double _Complex* rwork, a_int* info);
void pzneupd_c(MPI_Fint comm, bool rvec, char const* howmny, a_int const* select, double _Complex* d, double _Complex* z, a_int ldz, double _Complex sigma, double _Complex* workev, char const* bmat, a_int n, char const* which, a_int nev, double tol, double _Complex* resid, a_int ncv, double _Complex* v, a_int ldv, a_int* iparam, a_int* ipntr, double _Complex* workd, double _Complex* workl, a_int lworkl, double _Complex* rwork, a_int* info);

#ifdef  __cplusplus
}
#endif

#endif
