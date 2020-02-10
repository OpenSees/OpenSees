#ifndef __PARPACK_HPP__
#define __PARPACK_HPP__

#include "arpackdef.h"

#include "arpack.hpp"

#include <mpi.h>

#include <complex.h>
#include <complex>

namespace arpack {

namespace internal {
#include "parpack.h"
}  // namespace internal

inline void saupd(MPI_Fint comm, a_int& ido, bmat const bmat_option, a_int n,
                  which const which_option, a_int nev, float tol, float* resid,
                  a_int ncv, float* v, a_int ldv, a_int* iparam, a_int* ipntr,
                  float* workd, float* workl, a_int lworkl, a_int& info) {
  internal::pssaupd_c(comm, &ido, internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol, resid,
                      ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void seupd(MPI_Fint comm, bool rvec, howmny const howmny_option,
                  a_int* select, float* d, float* z, a_int ldz, float sigma,
                  bmat const bmat_option, a_int n, which const which_option,
                  a_int nev, float tol, float* resid, a_int ncv, float* v, a_int ldv,
                  a_int* iparam, a_int* ipntr, float* workd, float* workl,
                  a_int lworkl, a_int& info) {
  internal::psseupd_c(comm, rvec, internal::convert_to_char(howmny_option),
                      select, d, z, ldz, sigma,
                      internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol, resid,
                      ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void saupd(MPI_Fint comm, a_int& ido, bmat const bmat_option, a_int n,
                  which const which_option, a_int nev, double tol, double* resid,
                  a_int ncv, double* v, a_int ldv, a_int* iparam, a_int* ipntr,
                  double* workd, double* workl, a_int lworkl, a_int& info) {
  internal::pdsaupd_c(comm, &ido, internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol, resid,
                      ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void seupd(MPI_Fint comm, bool rvec, howmny const howmny_option,
                  a_int* select, double* d, double* z, a_int ldz, double sigma,
                  bmat const bmat_option, a_int n, which const which_option,
                  a_int nev, double tol, double* resid, a_int ncv, double* v,
                  a_int ldv, a_int* iparam, a_int* ipntr, double* workd,
                  double* workl, a_int lworkl, a_int& info) {
  internal::pdseupd_c(comm, rvec, internal::convert_to_char(howmny_option),
                      select, d, z, ldz, sigma,
                      internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol, resid,
                      ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void naupd(MPI_Fint comm, a_int& ido, bmat const bmat_option, a_int n,
                  which const which_option, a_int nev, float tol, float* resid,
                  a_int ncv, float* v, a_int ldv, a_int* iparam, a_int* ipntr,
                  float* workd, float* workl, a_int lworkl, a_int& info) {
  internal::psnaupd_c(comm, &ido, internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol, resid,
                      ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void neupd(MPI_Fint comm, bool rvec, howmny const howmny_option,
                  a_int* select, float* dr, float* di, float* z, a_int ldz,
                  float sigmar, float sigmai, float * workev,
                  bmat const bmat_option, a_int n,
                  which const which_option, a_int nev, float tol, float* resid,
                  a_int ncv, float* v, a_int ldv, a_int* iparam, a_int* ipntr,
                  float* workd, float* workl, a_int lworkl, a_int& info) {
  internal::psneupd_c(comm, rvec, internal::convert_to_char(howmny_option),
                      select, dr, di, z, ldz, sigmar, sigmai, workev,
                      internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol, resid,
                      ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void naupd(MPI_Fint comm, a_int& ido, bmat const bmat_option, a_int n,
                  which const which_option, a_int nev, double tol, double* resid,
                  a_int ncv, double* v, a_int ldv, a_int* iparam, a_int* ipntr,
                  double* workd, double* workl, a_int lworkl, a_int& info) {
  internal::pdnaupd_c(comm, &ido, internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol, resid,
                      ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void neupd(MPI_Fint comm, bool rvec, howmny const howmny_option,
                  a_int* select, double* dr, double* di, double* z, a_int ldz,
                  double sigmar, double sigmai, double * workev,
                  bmat const bmat_option, a_int n,
                  which const which_option, a_int nev, double tol, double* resid,
                  a_int ncv, double* v, a_int ldv, a_int* iparam, a_int* ipntr,
                  double* workd, double* workl, a_int lworkl, a_int& info) {
  internal::pdneupd_c(comm, rvec, internal::convert_to_char(howmny_option),
                      select, dr, di, z, ldz, sigmar, sigmai, workev,
                      internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol, resid,
                      ncv, v, ldv, iparam, ipntr, workd, workl, lworkl, &info);
}

inline void naupd(MPI_Fint comm, a_int& ido, bmat const bmat_option, a_int n,
                  which const which_option, a_int nev, float tol,
                  std::complex<float>* resid, a_int ncv, std::complex<float>* v,
                  a_int ldv, a_int* iparam, a_int* ipntr, std::complex<float>* workd,
                  std::complex<float>* workl, a_int lworkl,
                  std::complex<float>* rwork, a_int& info) {
  internal::pcnaupd_c(comm, &ido, internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol,
                      reinterpret_cast<_Complex float*>(resid), ncv,
                      reinterpret_cast<_Complex float*>(v), ldv, iparam, ipntr,
                      reinterpret_cast<_Complex float*>(workd),
                      reinterpret_cast<_Complex float*>(workl), lworkl,
                      reinterpret_cast<_Complex float*>(rwork), &info);
}

inline void neupd(MPI_Fint comm, bool rvec, howmny const howmny_option,
                  a_int* select, std::complex<float>* d, std::complex<float>* z,
                  a_int ldz, std::complex<float> sigma,
                  std::complex<float>* workev, bmat const bmat_option, a_int n,
                  which const which_option, a_int nev, float tol,
                  std::complex<float>* resid, a_int ncv, std::complex<float>* v,
                  a_int ldv, a_int* iparam, a_int* ipntr, std::complex<float>* workd,
                  std::complex<float>* workl, a_int lworkl,
                  std::complex<float>* rwork, a_int& info)

{
  internal::pcneupd_c(comm, rvec, internal::convert_to_char(howmny_option),
                      select, reinterpret_cast<_Complex float*>(d),
                      reinterpret_cast<_Complex float*>(z), ldz,
                      std::real(sigma) + _Complex_I * std::imag(sigma),
                      reinterpret_cast<_Complex float*>(workev),
                      internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol,
                      reinterpret_cast<_Complex float*>(resid), ncv,
                      reinterpret_cast<_Complex float*>(v), ldv, iparam, ipntr,
                      reinterpret_cast<_Complex float*>(workd),
                      reinterpret_cast<_Complex float*>(workl), lworkl,
                      reinterpret_cast<_Complex float*>(rwork), &info);
}

inline void naupd(MPI_Fint comm, a_int& ido, bmat const bmat_option, a_int n,
                  which const which_option, a_int nev, double tol,
                  std::complex<double>* resid, a_int ncv, std::complex<double>* v,
                  a_int ldv, a_int* iparam, a_int* ipntr, std::complex<double>* workd,
                  std::complex<double>* workl, a_int lworkl,
                  std::complex<double>* rwork, a_int& info) {
  internal::pznaupd_c(comm, &ido, internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol,
                      reinterpret_cast<_Complex double*>(resid), ncv,
                      reinterpret_cast<_Complex double*>(v), ldv, iparam, ipntr,
                      reinterpret_cast<_Complex double*>(workd),
                      reinterpret_cast<_Complex double*>(workl), lworkl,
                      reinterpret_cast<_Complex double*>(rwork), &info);
}

inline void neupd(MPI_Fint comm, bool rvec, howmny const howmny_option,
                  a_int* select, std::complex<double>* d, std::complex<double>* z,
                  a_int ldz, std::complex<double> sigma,
                  std::complex<double>* workev, bmat const bmat_option, a_int n,
                  which const which_option, a_int nev, double tol,
                  std::complex<double>* resid, a_int ncv, std::complex<double>* v,
                  a_int ldv, a_int* iparam, a_int* ipntr, std::complex<double>* workd,
                  std::complex<double>* workl, a_int lworkl,
                  std::complex<double>* rwork, a_int& info) {
  internal::pzneupd_c(comm, rvec, internal::convert_to_char(howmny_option),
                      select, reinterpret_cast<_Complex double*>(d),
                      reinterpret_cast<_Complex double*>(z), ldz,
                      std::real(sigma) + _Complex_I * std::imag(sigma),
                      reinterpret_cast<_Complex double*>(workev),
                      internal::convert_to_char(bmat_option), n,
                      internal::convert_to_char(which_option), nev, tol,
                      reinterpret_cast<_Complex double*>(resid), ncv,
                      reinterpret_cast<_Complex double*>(v), ldv, iparam, ipntr,
                      reinterpret_cast<_Complex double*>(workd),
                      reinterpret_cast<_Complex double*>(workl), lworkl,
                      reinterpret_cast<_Complex double*>(rwork), &info);
}
}  // namespace arpack
#endif
