//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, Claudio M. Perez
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
#ifndef blasdecl_H
#define blasdecl_H

#ifndef _WIN32
# define  DAXPY  daxpy_
# define  DSCAL  dscal_
# define  DGEMV  dgemv_
// Level 3
# define  DGETRF dgetrf_  
# define  DGETRI dgetri_  
# define  DGEMM  dgemm_   
// Lapack
# define  DGESV  dgesv_
# define  DGETRS dgetrs_
# define  DGBSV  dgbsv_
# define  DGBTRS dgbtrs_
# define  DPBSV  dpbsv_
# define  DPBTRS dpbtrs_
#endif

extern "C" {
  // Level 1
  // x = alpha * x + beta * y
  void DAXPY (int*, double*, double*, const int*, double*, const int*);
  void DSCAL (int*, double*, double*, const int*);
  void DGEMV (const char* trans, int* M, int* N,
              double* alpha,
              const double* A, int* lda,
              double* X, int* incX,
              double* beta,
              double* Y, int* incY);
// Level 3
#if 0
  void DGESV(int *N, int *NRHS, double *A, int *LDA, 
            int *iPiv, double *B, int *LDB, int *INFO);
  void DGETRS(char *TRANS, unsigned int sizeT,
            int *N, int *NRHS, double *A, int *LDA, 
            int *iPiv, double *B, int *LDB, int *INFO);
#endif

  void DGETRF(int *M, int *N, double *A, int *LDA, 
              int *iPiv, int *INFO);

  void DGETRI(int *N, double *A, int *LDA, 
                     int *iPiv, double *Work, int *WORKL, int *INFO);

  void DGEMM(const char* transA, const char* transB, int* M, int* N, int* K,
             double* alpha,
             double* A, const int* lda,
             double* B, const int* ldb,
             double* beta,
             double* C, const int* ldc);

//
// Lapack
//
// FullGen
int  DGESV(int *N, int *NRHS, double *A, int *LDA, 
          int *iPiv, double *B, int *LDB, int *INFO);
         
int  DGETRS(char *TRANS,
           int *N, int *NRHS, double *A, int *LDA, 
           int *iPiv, double *B, int *LDB, int *INFO);

// BandGen
int DGBSV(int *N, int *KL, int *KU, int *NRHS, double *A, 
          int *LDA, int *iPiv, double *B, int *LDB, 
          int *INFO);

int DGBTRS(char *TRANS, 
           int *N, int *KL, int *KU, int *NRHS,
           double *A, int *LDA, int *iPiv, 
           double *B, int *LDB, int *INFO);

// BandSPD
int  DPBSV(char *UPLO,
          int *N, int *KD, int *NRHS, 
          double *A, int *LDA, double *B, int *LDB, 
          int *INFO);

int  DPBTRS(char *UPLO,
           int *N, int *KD, int *NRHS, 
           double *A, int *LDA, double *B, int *LDB, 
           int *INFO);
}

#endif // blasdecl_H
