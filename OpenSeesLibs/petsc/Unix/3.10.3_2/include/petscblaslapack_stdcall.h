/*
  This file deals with
     1) STDCALL BLAS/LAPACK calling conventions
     2) capitalization of BLAS/LAPACK names
     3) character strings passing conventions where the string length is passed in the stack immediately after the string
        as opposed to at the end of the function arguements.

  This issue only comes up on some Microsoft Windows compilers.
*/
#if !defined(_BLASLAPACK_STDCALL_H)
#define _BLASLAPACK_STDCALL_H

/* PETSC_STDCALL is defined on some Microsoft Windows systems and is used for functions compiled by the Fortran compiler */
#if !defined(PETSC_STDCALL)
#define PETSC_STDCALL
#endif

#if !defined(PETSC_USE_COMPLEX)
# if defined(PETSC_BLASLAPACK_SINGLEISDOUBLE) || defined(PETSC_USE_REAL_SINGLE)
/* Real single precision without character string arguments. */
#  define LAPACKgeqrf_ SGEQRF
#  define LAPACKorgqr_ SORGQR
#  define LAPACKgetrf_ SGETRF
#  define LAPACKgetri_ SGETRI
#  define BLASdot_     SDOT
#  define BLASdotu_    SDOT
#  define BLASnrm2_    SNRM2
#  define BLASscal_    SSCAL
#  define BLAScopy_    SCOPY
#  define BLASswap_    SSWAP
#  define BLASaxpy_    SAXPY
#  define BLASasum_    SASUM
#  define BLAStrmv_    STRMV
#  define LAPACKpttrf_ SPTTRF
#  define LAPACKpttrs_ SPTTRS
#  define LAPACKstein_ SSTEIN
#  define LAPACKgesv_  SGESV
#  define LAPACKgelss_ SGELSS
#  define LAPACKtgsen_ STGSEN
/* Real single precision with character string arguments. */
#  define LAPACKpotrf_(a,b,c,d,e)                   SPOTRF((a),1,(b),(c),(d),(e))
PETSC_EXTERN void PETSC_STDCALL SPOTRF(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKpotrs_(a,b,c,d,e,f,g,h)             SPOTRS((a),1,(b),(c),(d),(e),(f),(g),(h))
PETSC_EXTERN void PETSC_STDCALL SPOTRS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKpotri_(a,b,c,d,e)                   SPOTRI((a),1,(b),(c),(d),(e))
PETSC_EXTERN void PETSC_STDCALL SPOTRI(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsytrf_(a,b,c,d,e,f,g,h)             SSYTRF((a),1,(b),(c),(d),(e),(f),(g),(h))
PETSC_EXTERN void PETSC_STDCALL SSYTRF(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsytrs_(a,b,c,d,e,f,g,h,i)           SSYTRS((a),1,(b),(c),(d),(e),(f),(g),(h),(i))
PETSC_EXTERN void PETSC_STDCALL SSYTRS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsytri_(a,b,c,d,e,f,g)               SSYTRI((a),1,(b),(c),(d),(e),(f),(g))
PETSC_EXTERN void PETSC_STDCALL SSYTRI(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*);
#  define LAPACKgels_(a,b,c,d,e,f,g,h,i,j,k)        SGELS((a),1,(b),(c),(d),(e),(f),(g),(h),(i),(j),(k))
PETSC_EXTERN void PETSC_STDCALL SGELS(const char*,int,const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*);
#  define BLASgemm_(a,b,c,d,e,f,g,h,i,j,k,l,m)      SGEMM((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m))
PETSC_EXTERN void PETSC_STDCALL SGEMM(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*);
#  define BLASsymm_(a,b,c,d,e,f,g,h,i,j,k,l)        SSYMM((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l))
PETSC_EXTERN void PETSC_STDCALL SSYMM(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*);
#  define BLASsyrk_(a,b,c,d,e,f,g,h,i,j)            SSYRK((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j))
PETSC_EXTERN void PETSC_STDCALL SSYRK(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*);
#  define BLASsyr2k_(a,b,c,d,e,f,g,h,i,j,k,l)       SSYR2K((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l))
PETSC_EXTERN void PETSC_STDCALL SSYR2K(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*);
#  define BLAStrsm_(a,b,c,d,e,f,g,h,i,j,k,l,m)      STRSM((a),1,(b),1,(c),1,(d),1,(e),(f),(g),(h),(i),(j),(k))
PETSC_EXTERN void PETSC_STDCALL STRSM(const char*,int,const char*,int,const char*,int,const char*,int,const PetscBLASInt*,const PetscBLASInt*,const PetscScalar*,const PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*);
#  define LAPACKgetrs_(a,b,c,d,e,f,g,h,i)           SGETRS((a),1,(b),(c),(d),(e),(f),(g),(h),(i))
PETSC_EXTERN void PETSC_STDCALL SGETRS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define BLASgemv_(a,b,c,d,e,f,g,h,i,j,k)          SGEMV((a),1,(b),(c),(d),(e),(f),(g),(h),(i),(j),(k))
PETSC_EXTERN void PETSC_STDCALL SGEMV(const char*,int,const PetscBLASInt*,const PetscBLASInt*,const PetscScalar*,const PetscScalar*,const PetscBLASInt*,const PetscScalar *,const PetscBLASInt*,const PetscScalar*,PetscScalar*,const PetscBLASInt*);
#  define LAPACKgeev_(a,b,c,d,e,f,g,h,i,j,k,l,m,n)  SGEEV((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n))
PETSC_EXTERN void PETSC_STDCALL SGEEV(const char*,int,const char*,int,PetscBLASInt*,PetscScalar *,PetscBLASInt*,PetscReal*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKgesvd_(a,b,c,d,e,f,g,h,i,j,k,l,m,n) SGESVD((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n))
PETSC_EXTERN void PETSC_STDCALL SGESVD(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar *,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);

#  define LAPACKsyev_(a,b,c,d,e,f,g,h,i)            SSYEV((a),(b),1,(c),1,(d),(e),(f),(g),(h),(i))
PETSC_EXTERN void PETSC_STDCALL SSYEV(const char*,const char*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsyevx_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t) SSYEVX((a),(b),1,(c),1,(d),1,(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r),(s),(t))
PETSC_EXTERN void PETSC_STDCALL SSYEVX(const char*,const char*,const char*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);

#  define LAPACKsygv_(a,b,c,d,e,f,g,h,i,j,k,l)      SSYGV((a),(b),1,(c),1,(d),(e),(f),(g),(h),(i),(j),(k),(l))
PETSC_EXTERN void PETSC_STDCALL SSYGV(PetscBLASInt*,const char*,int,const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsygvx_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w) SSYGVX((a),(b),1,(c),1,(d),1,(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r),(s),(t),(u),(v),(w))
PETSC_EXTERN void PETSC_STDCALL SSYGVX(PetscBLASInt*,const char*,int,const char*,int,const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);

#  define LAPACKstebz_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) SSTEBZ((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r))
PETSC_EXTERN void PETSC_STDCALL SSTEBZ(const char*,int,const char*,int,PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsteqr_(a,b,c,d,e,f,g,h)                     SSTEQR((a),1,(b),(c),(d),(e),(f),(g),(h))
PETSC_EXTERN void PETSC_STDCALL SSTEQR(const char*,int,PetscBLASInt*,PetscReal*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*);
PETSC_EXTERN void PETSC_STDCALL SPTTRS(PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKgerfs_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q) SGERFS((a),1,(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q))
PETSC_EXTERN void PETSC_STDCALL SGERFS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKtrsen_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) STRSEN((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r))
PETSC_EXTERN void PETSC_STDCALL STRSEN(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt*,PetscReal*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);

#  define LAPACKgges_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u) SGEES((a),1,(b),1,(c),1,(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r),(s),(t),(u))
PETSC_EXTERN void PETSC_STDCALL SGGES(const char*,int,const char*,int,const char*,int,PetscBLASInt(*)(),PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKhseqr_(a,b,c,d,e,f,g,h,i,j,k,l,m,n) SHSEQR((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n))
PETSC_EXTERN void PETSC_STDCALL SHSEQR(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
# else
/* Real double precision without character string arguments. */
#  define LAPACKgeqrf_ DGEQRF
#  define LAPACKorgqr_ DORGQR
#  define LAPACKgetrf_ DGETRF
#  define LAPACKgetri_ DGETRI
#  define BLASdot_     DDOT
#  define BLASdotu_    DDOT
#  define BLASnrm2_    DNRM2
#  define BLASscal_    DSCAL
#  define BLAScopy_    DCOPY
#  define BLASswap_    DSWAP
#  define BLASaxpy_    DAXPY
#  define BLASasum_    DASUM
#  define BLAStrmv_    DTRMV
#  define LAPACKpttrf_ DPTTRF
#  define LAPACKpttrs_ DPTTRS
#  define LAPACKstein_ DSTEIN
#  define LAPACKgesv_  DGESV
#  define LAPACKgelss_ DGELSS
#  define LAPACKtgsen_ DTGSEN
/* Real double precision with character string arguments. */
#  define LAPACKpotrf_(a,b,c,d,e)                   DPOTRF((a),1,(b),(c),(d),(e))
PETSC_EXTERN void PETSC_STDCALL DPOTRF(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKpotrs_(a,b,c,d,e,f,g,h)             DPOTRS((a),1,(b),(c),(d),(e),(f),(g),(h))
PETSC_EXTERN void PETSC_STDCALL DPOTRS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKpotri_(a,b,c,d,e)                   DPOTRI((a),1,(b),(c),(d),(e))
PETSC_EXTERN void PETSC_STDCALL DPOTRI(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsytrf_(a,b,c,d,e,f,g,h)             DSYTRF((a),1,(b),(c),(d),(e),(f),(g),(h))
PETSC_EXTERN void PETSC_STDCALL DSYTRF(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsytrs_(a,b,c,d,e,f,g,h,i)           DSYTRS((a),1,(b),(c),(d),(e),(f),(g),(h),(i))
PETSC_EXTERN void PETSC_STDCALL DSYTRS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsytri_(a,b,c,d,e,f,g)               DSYTRI((a),1,(b),(c),(d),(e),(f),(g))
PETSC_EXTERN void PETSC_STDCALL DSYTRI(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*);
#  define BLASgemv_(a,b,c,d,e,f,g,h,i,j,k)          DGEMV((a),1,(b),(c),(d),(e),(f),(g),(h),(i),(j),(k))
PETSC_EXTERN void PETSC_STDCALL DGEMV(const char*,int,const PetscBLASInt*,const PetscBLASInt*,const PetscScalar*,const PetscScalar*,const PetscBLASInt*,const PetscScalar *,const PetscBLASInt*,const PetscScalar*,PetscScalar*,const PetscBLASInt*);
#  define LAPACKgetrs_(a,b,c,d,e,f,g,h,i)           DGETRS((a),1,(b),(c),(d),(e),(f),(g),(h),(i))
PETSC_EXTERN void PETSC_STDCALL DGETRS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKgels_(a,b,c,d,e,f,g,h,i,j,k)        DGELS((a),1,(b),(c),(d),(e),(f),(g),(h),(i),(j),(k))
PETSC_EXTERN void PETSC_STDCALL DGELS(const char*,int,const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*);
#  define BLASgemm_(a,b,c,d,e,f,g,h,i,j,k,l,m)      DGEMM((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m))
PETSC_EXTERN void PETSC_STDCALL DGEMM(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*);
#  define BLASsymm_(a,b,c,d,e,f,g,h,i,j,k,l)        DSYMM((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l))
PETSC_EXTERN void PETSC_STDCALL DSYMM(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*);
#  define BLASsyrk_(a,b,c,d,e,f,g,h,i,j)            DSYRK((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j))
PETSC_EXTERN void PETSC_STDCALL DSYRK(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*);
#  define BLASsyr2k_(a,b,c,d,e,f,g,h,i,j,k,l)       DSYR2K((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l))
PETSC_EXTERN void PETSC_STDCALL DSYR2K(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*);
#  define BLAStrsm_(a,b,c,d,e,f,g,h,i,j,k,l,m)      DTRSM((a),1,(b),1,(c),1,(d),1,(e),(f),(g),(h),(i),(j),(k))
PETSC_EXTERN void PETSC_STDCALL DTRSM(const char*,int,const char*,int,const char*,int,const char*,int,const PetscBLASInt*,const PetscBLASInt*,const PetscScalar*,const PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*);
#  define LAPACKgesvd_(a,b,c,d,e,f,g,h,i,j,k,l,m,n) DGESVD((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n))
PETSC_EXTERN void PETSC_STDCALL DGESVD(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar *,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKgeev_(a,b,c,d,e,f,g,h,i,j,k,l,m,n)  DGEEV((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n))
PETSC_EXTERN void PETSC_STDCALL DGEEV(const char*,int,const char*,int,PetscBLASInt*,PetscScalar *,PetscBLASInt*,PetscReal*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);

#  define LAPACKsyev_(a,b,c,d,e,f,g,h,i)            DSYEV((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i))
PETSC_EXTERN void PETSC_STDCALL DSYEV(const char*,int,const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsyevx_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t) DSYEVX((a),1,(b),1,(c),1,(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r),(s),(t))
PETSC_EXTERN void PETSC_STDCALL DSYEVX(const char*,int,const char*,int,const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);

#  define LAPACKsygv_(a,b,c,d,e,f,g,h,i,j,k,l)      DSYGV((a),(b),1,(c),1,(d),(e),(f),(g),(h),(i),(j),(k),(l))
PETSC_EXTERN void PETSC_STDCALL DSYGV(PetscBLASInt*,const char*,int,const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsygvx_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w) DSYGVX((a),(b),1,(c),1,(d),1,(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r),(s),(t),(u),(v),(w))
PETSC_EXTERN void PETSC_STDCALL DSYGVX(PetscBLASInt*,const char*,int,const char*,int,const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKstebz_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) DSTEBZ((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r))
PETSC_EXTERN void PETSC_STDCALL DSTEBZ(const char*,int,const char*,int,PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsteqr_(a,b,c,d,e,f,g,h)                     DSTEQR((a),1,(b),(c),(d),(e),(f),(g),(h))
PETSC_EXTERN void PETSC_STDCALL DSTEQR(const char*,int,PetscBLASInt*,PetscReal*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*);
PETSC_EXTERN void PETSC_STDCALL DPTTRS(PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKgerfs_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q) DGERFS((a),1,(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q))
PETSC_EXTERN void PETSC_STDCALL DGERFS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKtrsen_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) DTRSEN((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r))
PETSC_EXTERN void PETSC_STDCALL DTRSEN(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt*,PetscReal*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKgges_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u) DGGES((a),1,(b),1,(c),1,(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r),(s),(t),(u))
PETSC_EXTERN void PETSC_STDCALL DGGES(const char*,int,const char*,int,const char*,int,PetscBLASInt(*)(),PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKhseqr_(a,b,c,d,e,f,g,h,i,j,k,l,m,n) DHSEQR((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n))
PETSC_EXTERN void PETSC_STDCALL DHSEQR(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
# endif
#else
# if defined(PETSC_BLASLAPACK_SINGLEISDOUBLE) || defined(PETSC_USE_REAL_SINGLE)
/* Complex single precision without character string arguments. */
#  define LAPACKgeqrf_ CGEQRF
#  define LAPACKorgqr_ CUNGQR
#  define LAPACKgetrf_ CGETRF
#  define LAPACKgetri_ CGETRI
/* #  define BLASdot_     CDOTC */
/* #  define BLASdotu_    CDOTU */
#  define BLASnrm2_    SCNRM2
#  define BLASscal_    CSCAL
#  define BLAScopy_    CCOPY
#  define BLASswap_    CSWAP
#  define BLASaxpy_    CAXPY
#  define BLASasum_    SCASUM
#  define BLAStrmv_    CTRMV
#  define LAPACKpttrf_ CPTTRF
#  define LAPACKstein_ CSTEIN
#  define LAPACKgesv_  CGESV
#  define LAPACKgelss_ CGELSS
#  define LAPACKtgsen_ CTGSEN
/* Complex single precision with character string arguments. */
#  define LAPACKpotrf_(a,b,c,d,e)                   CPOTRF((a),1,(b),(c),(d),(e))
PETSC_EXTERN void PETSC_STDCALL CPOTRF(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKpotrs_(a,b,c,d,e,f,g,h)             CPOTRS((a),1,(b),(c),(d),(e),(f),(g),(h))
PETSC_EXTERN void PETSC_STDCALL CPOTRS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKpotri_(a,b,c,d,e)                   CPOTRI((a),1,(b),(c),(d),(e))
PETSC_EXTERN void PETSC_STDCALL CPOTRI(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsytrf_(a,b,c,d,e,f,g,h)             CSYTRF((a),1,(b),(c),(d),(e),(f),(g),(h))
PETSC_EXTERN void PETSC_STDCALL CSYTRF(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsytrs_(a,b,c,d,e,f,g,h,i)           CSYTRS((a),1,(b),(c),(d),(e),(f),(g),(h),(i))
PETSC_EXTERN void PETSC_STDCALL CSYTRS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsytri_(a,b,c,d,e,f,g)               CSYTRI((a),1,(b),(c),(d),(e),(f),(g))
PETSC_EXTERN void PETSC_STDCALL CSYTRI(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*);
#  define LAPACKhetrf_(a,b,c,d,e,f,g,h)             CHETRF((a),1,(b),(c),(d),(e),(f),(g),(h))
PETSC_EXTERN void PETSC_STDCALL CHETRF(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKhetrs_(a,b,c,d,e,f,g,h,i)           CHETRS((a),1,(b),(c),(d),(e),(f),(g),(h),(i))
PETSC_EXTERN void PETSC_STDCALL CHETRS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKhetri_(a,b,c,d,e,f,g)               CHETRI((a),1,(b),(c),(d),(e),(f),(g))
PETSC_EXTERN void PETSC_STDCALL CHETRI(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*);
#  define BLASgemv_(a,b,c,d,e,f,g,h,i,j,k)          CGEMV((a),1,(b),(c),(d),(e),(f),(g),(h),(i),(j),(k))
PETSC_EXTERN void PETSC_STDCALL CGEMV(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar *,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*);
#  define LAPACKgetrs_(a,b,c,d,e,f,g,h,i)           CGETRS((a),1,(b),(c),(d),(e),(f),(g),(h),(i))
PETSC_EXTERN void PETSC_STDCALL CGETRS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKgels_(a,b,c,d,e,f,g,h,i,j,k)        CGELS((a),1,(b),(c),(d),(e),(f),(g),(h),(i),(j),(k))
PETSC_EXTERN void PETSC_STDCALL CGELS(const char*,int,const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*);
#  define BLASgemm_(a,b,c,d,e,f,g,h,i,j,k,l,m)      CGEMM((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m))
PETSC_EXTERN void PETSC_STDCALL CGEMM(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*);
#  define BLASsymm_(a,b,c,d,e,f,g,h,i,j,k,l)        CSYMM((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l))
PETSC_EXTERN void PETSC_STDCALL CSYMM(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*);
#  define BLASsyrk_(a,b,c,d,e,f,g,h,i,j)            CSYRK((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j))
PETSC_EXTERN void PETSC_STDCALL CSYRK(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*);
#  define BLASsyr2k_(a,b,c,d,e,f,g,h,i,j,k,l)       CSYR2K((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l))
PETSC_EXTERN void PETSC_STDCALL CSYR2K(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*);
#  define BLAStrsm_(a,b,c,d,e,f,g,h,i,j,k,l,m)      CTRSM((a),1,(b),1,(c),1,(d),1,(e),(f),(g),(h),(i),(j),(k))
PETSC_EXTERN void PETSC_STDCALL CTRSM(const char*,int,const char*,int,const char*,int,const char*,int,const PetscBLASInt*,const PetscBLASInt*,const PetscScalar*,const PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*);
#  define LAPACKgeev_(a,b,c,d,e,f,g,h,i,j,k,l,m,n)  CGEEV((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n))
PETSC_EXTERN void PETSC_STDCALL CGEEV(const char*,int,const char*,int,PetscBLASInt*,PetscScalar *,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*);
#  define LAPACKgesvd_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) CGESVD((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o))
PETSC_EXTERN void PETSC_STDCALL CGESVD(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar *,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*);

#  define LAPACKsyev_(a,b,c,d,e,f,g,h,i,j)          CHEEV((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j))
PETSC_EXTERN void PETSC_STDCALL CHEEV(const char*,int,const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*);
#  define LAPACKsyevx_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u) CHEEVX((a),1,(b),1,(c),1,(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r),(s),(t),(u))
PETSC_EXTERN void PETSC_STDCALL CHEEVX(const char*,int,const char*,int,const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*);

#  define LAPACKsygv_(a,b,c,d,e,f,g,h,i,j,k,l,m)      CHEGV((a),(b),1,(c),1,(d),(e),(f),(g),(h),(i),(j),(k),(l),(m))
PETSC_EXTERN void PETSC_STDCALL CHEGV(PetscBLASInt*,const char*,int,const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*);
#  define LAPACKsygvx_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x) CHEGVX((a),(b),1,(c),1,(d),1,(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r),(s),(t),(u),(v),(w),(x))
PETSC_EXTERN void PETSC_STDCALL CHEGVX(PetscBLASInt*,const char*,int,const char*,int,const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsteqr_(a,b,c,d,e,f,g,h)                     CSTEQR((a),1,(b),(c),(d),(e),(f),(g),(h))
PETSC_EXTERN void PETSC_STDCALL CSTEQR(const char*,int,PetscBLASInt*,PetscReal*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*);
#  define LAPACKpttrs_(a,b,c,d,e,f,g,h) CPTTRS((a),1,(b),(c),(d),(e),(f),(g),(h))
PETSC_EXTERN void PETSC_STDCALL CPTTRS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKgerfs_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) CGERFS((a),1,(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r))
PETSC_EXTERN void PETSC_STDCALL CGERFS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscReal*,PetscBLASInt*);
#  define LAPACKtrsen_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) CTRSEN((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o))
PETSC_EXTERN void PETSC_STDCALL CTRSEN(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKgges_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u) CGGES((a),1,(b),1,(c),1,(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r),(s),(t),(u))
PETSC_EXTERN void PETSC_STDCALL CGGES(const char*,int,const char*,int,const char*,int,PetscBLASInt(*)(),PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKhseqr_(a,b,c,d,e,f,g,h,i,j,k,l,m) CHSEQR((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m))
PETSC_EXTERN void PETSC_STDCALL CHSEQR(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
# else
/* Complex double precision without character string arguments */
#  define LAPACKgeqrf_ ZGEQRF
#  define LAPACKorgqr_ ZUNGQR
#  define LAPACKgetrf_ ZGETRF
#  define LAPACKgetri_ ZGETRI
/* #  define BLASdot_     ZDOTC */
/* #  define BLASdotu_    ZDOTU */
#  define BLASnrm2_    DZNRM2
#  define BLASscal_    ZSCAL
#  define BLAScopy_    ZCOPY
#  define BLASswap_    ZSWAP
#  define BLASaxpy_    ZAXPY
#  define BLASasum_    DZASUM
#  define BLAStrmv_    ZTRMV
#  define LAPACKpttrf_ ZPTTRF
#  define LAPACKstein_ ZSTEIN
#  define LAPACKgesv_  ZGESV
#  define LAPACKgelss_ ZGELSS
#  define LAPACKtgsen_ ZTGSEN
/* Complex double precision with character string arguments */
#  define LAPACKpotrf_(a,b,c,d,e)                   ZPOTRF((a),1,(b),(c),(d),(e))
PETSC_EXTERN void PETSC_STDCALL ZPOTRF(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKpotrs_(a,b,c,d,e,f,g,h)             ZPOTRS((a),1,(b),(c),(d),(e),(f),(g),(h))
PETSC_EXTERN void PETSC_STDCALL ZPOTRS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKpotri_(a,b,c,d,e)                   ZPOTRI((a),1,(b),(c),(d),(e))
PETSC_EXTERN void PETSC_STDCALL ZPOTRI(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsytrf_(a,b,c,d,e,f,g,h)             ZSYTRF((a),1,(b),(c),(d),(e),(f),(g),(h))
PETSC_EXTERN void PETSC_STDCALL ZSYTRF(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsytrs_(a,b,c,d,e,f,g,h,i)           ZSYTRS((a),1,(b),(c),(d),(e),(f),(g),(h),(i))
PETSC_EXTERN void PETSC_STDCALL ZSYTRS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsytri_(a,b,c,d,e,f,g)               ZSYTRI((a),1,(b),(c),(d),(e),(f),(g))
PETSC_EXTERN void PETSC_STDCALL ZSYTRI(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*);
#  define LAPACKhetrf_(a,b,c,d,e,f,g,h)             ZHETRF((a),1,(b),(c),(d),(e),(f),(g),(h))
PETSC_EXTERN void PETSC_STDCALL ZHETRF(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKhetrs_(a,b,c,d,e,f,g,h,i)           ZHETRS((a),1,(b),(c),(d),(e),(f),(g),(h),(i))
PETSC_EXTERN void PETSC_STDCALL ZHETRS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKhetri_(a,b,c,d,e,f,g)               ZHETRI((a),1,(b),(c),(d),(e),(f),(g))
PETSC_EXTERN void PETSC_STDCALL ZHETRI(const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*);
#  define BLASgemv_(a,b,c,d,e,f,g,h,i,j,k)          ZGEMV((a),1,(b),(c),(d),(e),(f),(g),(h),(i),(j),(k))
PETSC_EXTERN void PETSC_STDCALL ZGEMV(const char*,const int,const PetscBLASInt*,const PetscBLASInt*,const PetscScalar*,const PetscScalar*,const PetscBLASInt*,const PetscScalar *,const PetscBLASInt*,const PetscScalar*,PetscScalar*,const PetscBLASInt*);
#  define LAPACKgetrs_(a,b,c,d,e,f,g,h,i)           ZGETRS((a),1,(b),(c),(d),(e),(f),(g),(h),(i))
PETSC_EXTERN void PETSC_STDCALL ZGETRS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKgels_(a,b,c,d,e,f,g,h,i,j,k)        ZGELS((a),1,(b),(c),(d),(e),(f),(g),(h),(i),(j),(k))
PETSC_EXTERN void PETSC_STDCALL ZGELS(const char*,int,const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*);
#  define BLASgemm_(a,b,c,d,e,f,g,h,i,j,k,l,m)      ZGEMM((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m))
PETSC_EXTERN void PETSC_STDCALL ZGEMM(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*);
#  define BLASsymm_(a,b,c,d,e,f,g,h,i,j,k,l)        ZSYMM((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l))
PETSC_EXTERN void PETSC_STDCALL ZSYMM(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*);
#  define BLASsyrk_(a,b,c,d,e,f,g,h,i,j)            ZSYRK((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j))
PETSC_EXTERN void PETSC_STDCALL ZSYRK(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*);
#  define BLASsyr2k_(a,b,c,d,e,f,g,h,i,j,k,l)       ZSYR2K((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l))
PETSC_EXTERN void PETSC_STDCALL ZSYR2K(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*);
#  define BLAStrsm_(a,b,c,d,e,f,g,h,i,j,k,l,m)      ZTRSM((a),1,(b),1,(c),1,(d),1,(e),(f),(g),(h),(i),(j),(k))
PETSC_EXTERN void PETSC_STDCALL ZTRSM(const char*,int,const char*,int,const char*,int,const char*,int,const PetscBLASInt*,const PetscBLASInt*,const PetscScalar*,const PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*);
#  define LAPACKgeev_(a,b,c,d,e,f,g,h,i,j,k,l,m,n)  ZGEEV((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n))
PETSC_EXTERN void PETSC_STDCALL ZGEEV(const char*,int,const char*,int,PetscBLASInt*,PetscScalar *,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*);
#  define LAPACKgesvd_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) ZGESVD((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o))
PETSC_EXTERN void PETSC_STDCALL ZGESVD(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar *,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*);

#  define LAPACKsyev_(a,b,c,d,e,f,g,h,i,j)           ZHEEV((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j))
PETSC_EXTERN void PETSC_STDCALL ZHEEV(const char*,int,const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*);
#  define LAPACKsyevx_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u)    ZHEEVX((a),1,(b),1,(c),1,(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r),(s),(t),(u))
PETSC_EXTERN void PETSC_STDCALL ZHEEVX(const char*,int,const char*,int,const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);

#  define LAPACKsygv_(a,b,c,d,e,f,g,h,i,j,k,l,m)      ZHEGV((a),(b),1,(c),1,(d),(e),(f),(g),(h),(i),(j),(k),(l),(m))
PETSC_EXTERN void PETSC_STDCALL ZHEGV(PetscBLASInt*,const char*,int,const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*);
#  define LAPACKsygvx_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x)    ZHEGVX((a),(b),1,(c),1,(d),1,(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r),(s),(t),(u),(v),(w),(x))
PETSC_EXTERN void PETSC_STDCALL ZHEGVX(PetscBLASInt*,const char*,int,const char*,int,const char*,int,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKsteqr_(a,b,c,d,e,f,g,h)                     ZSTEQR((a),1,(b),(c),(d),(e),(f),(g),(h))
PETSC_EXTERN void PETSC_STDCALL ZSTEQR(const char*,int,PetscBLASInt*,PetscReal*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*);
#  define LAPACKpttrs_(a,b,c,d,e,f,g,h)             ZPTTRS((a),1,(b),(c),(d),(e),(f),(g),(h))
PETSC_EXTERN void PETSC_STDCALL ZPTTRS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKgerfs_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r) ZGERFS((a),1,(b),(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r))
PETSC_EXTERN void PETSC_STDCALL ZGERFS(const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscReal*,PetscBLASInt*);
#  define LAPACKtrsen_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o) ZTRSEN((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o))
PETSC_EXTERN void PETSC_STDCALL ZTRSEN(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKgges_(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u) ZGGES((a),1,(b),1,(c),1,(d),(e),(f),(g),(h),(i),(j),(k),(l),(m),(n),(o),(p),(q),(r),(s),(t),(u))
PETSC_EXTERN void PETSC_STDCALL ZGEES(const char*,int,const char*,int,const char*,int,PetscBLASInt(*)(),PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
#  define LAPACKhseqr_(a,b,c,d,e,f,g,h,i,j,k,l,m) ZHSEQR((a),1,(b),1,(c),(d),(e),(f),(g),(h),(i),(j),(k),(l),(m))
PETSC_EXTERN void PETSC_STDCALL ZHSEQR(const char*,int,const char*,int,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
# endif
#endif

PETSC_EXTERN void PETSC_STDCALL LAPACKgetrf_(PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
PETSC_EXTERN void PETSC_STDCALL LAPACKgetri_(PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
PETSC_EXTERN void PETSC_STDCALL LAPACKgeqrf_(PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
PETSC_EXTERN void PETSC_STDCALL LAPACKorgqr_(PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscBLASInt*);
PETSC_EXTERN void PETSC_STDCALL LAPACKpttrf_(PetscBLASInt*,PetscReal*,PetscScalar*,PetscBLASInt*);
PETSC_EXTERN void PETSC_STDCALL LAPACKstein_(PetscBLASInt*,PetscReal*,PetscReal*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
PETSC_EXTERN void PETSC_STDCALL LAPACKgesv_(const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*);

#if defined(PETSC_USE_COMPLEX)
PETSC_EXTERN void PETSC_STDCALL LAPACKgelss_(const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,const PetscReal*,PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,PetscBLASInt*);
PETSC_EXTERN void PETSC_STDCALL LAPACKtgsen_(PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscScalar*,PetscScalar*,PetscBLASInt*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscReal*,PetscReal*,PetscScalar*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
#else
PETSC_EXTERN void PETSC_STDCALL LAPACKgelss_(const PetscBLASInt*,const PetscBLASInt*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscReal*,const PetscReal*,PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscBLASInt*);
PETSC_EXTERN void PETSC_STDCALL LAPACKtgsen_(PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscBLASInt*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscReal*,PetscReal*,PetscReal*,PetscReal*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*,PetscBLASInt*);
#endif


PETSC_EXTERN PetscScalar PETSC_STDCALL BLASdot_(const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*);
PETSC_EXTERN PetscReal PETSC_STDCALL BLASnrm2_(const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*);
PETSC_EXTERN void PETSC_STDCALL BLASscal_(const PetscBLASInt*,const PetscScalar*,PetscScalar*,const PetscBLASInt*);
PETSC_EXTERN void PETSC_STDCALL BLAScopy_(const PetscBLASInt*,const PetscScalar*,PetscBLASInt*,const PetscScalar*,const PetscBLASInt*);
PETSC_EXTERN void PETSC_STDCALL BLASswap_(const PetscBLASInt*,PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*);
PETSC_EXTERN void PETSC_STDCALL BLASaxpy_(const PetscBLASInt*,const PetscScalar*,const PetscScalar*,const PetscBLASInt*,PetscScalar*,const PetscBLASInt*);
PETSC_EXTERN PetscReal PETSC_STDCALL BLASasum_(const PetscBLASInt*,const PetscScalar*,const PetscBLASInt*);

#endif
