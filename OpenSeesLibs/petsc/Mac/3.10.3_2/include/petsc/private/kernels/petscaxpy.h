
/*
    PetscKernelAXPY -  X = X + alpha * Y

   Input Parameters:
+    X, Y - arrays
.    alpha - scalar
-    n - length of arrays

   Also PetscKernelAXPY2(), PetscKernelAXPY3(), PetscKernelAXPY4()

*/

#ifndef PetscKernelAXPY

#if defined(PETSC_USE_FORTRAN_KERNEL_MAXPY)
#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define fortranmaxpy4_ FORTRANMAXPY4
#define fortranmaxpy3_ FORTRANMAXPY3
#define fortranmaxpy2_ FORTRANMAXPY2
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define fortranmaxpy4_ fortranmaxpy4
#define fortranmaxpy3_ fortranmaxpy3
#define fortranmaxpy2_ fortranmaxpy2
#endif
PETSC_EXTERN void fortranmaxpy4_(void*,void*,void*,void*,void*,const void*,const void*,const void*,const void*,PetscInt*);
PETSC_EXTERN void fortranmaxpy3_(void*,void*,void*,void*,const void*,const void*,const void*,PetscInt*);
PETSC_EXTERN void fortranmaxpy2_(void*,void*,void*,const void*,const void*,PetscInt*);
#endif
#include <petscblaslapack.h>

#if defined(PETSC_USE_FORTRAN_KERNEL_MAXPY)
#define PetscKernelAXPY(U,a1,p1,n)                   {PetscBLASInt one=1; PetscBLASInt nn = (PetscBLASInt) n;  PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&nn,&a1,p1,&one,U,&one));}
#define PetscKernelAXPY2(U,a1,a2,p1,p2,n)            {fortranmaxpy2_(U,&a1,&a2,p1,p2,&n);}
#define PetscKernelAXPY3(U,a1,a2,a3,p1,p2,p3,n)      {fortranmaxpy3_(U,&a1,&a2,&a3,p1,p2,p3,&n);}
#define PetscKernelAXPY4(U,a1,a2,a3,a4,p1,p2,p3,p4,n){fortranmaxpy4_(U,&a1,&a2,&a3,&a4,p1,p2,p3,p4,&n);}

#elif defined(PETSC_USE_UNROLL_KERNELS)

#define PetscKernelAXPY(U,Alpha,P,n) {\
  switch (n & 0x3) {\
  case 3: *U++    += Alpha * *P++;\
  case 2: *U++    += Alpha * *P++;\
  case 1: *U++    += Alpha * *P++;\
  n -= 4;case 0: break;}while (n>0) {U[0] += Alpha * P[0];U[1] += Alpha * P[1]; U[2] += Alpha * P[2]; U[3] += Alpha * P[3];U += 4; P += 4; n -= 4;}}
#define PetscKernelAXPY2(U,a1,a2,p1,p2,n) {\
  switch (n & 0x3) {\
  case 3: *U++    += a1 * *p1++ + a2 * *p2++;\
  case 2: *U++    += a1 * *p1++ + a2 * *p2++;\
  case 1: *U++    += a1 * *p1++ + a2 * *p2++;\
  n -= 4;case 0: break;}\
  while (n>0) {U[0]+=a1*p1[0]+a2*p2[0];U[1]+=a1*p1[1]+a2*p2[1]; U[2]+=a1*p1[2]+a2*p2[2];U[3]+=a1*p1[3]+a2*p2[3];U+=4;p1+=4;p2+=4;n -= 4;}}
#define PetscKernelAXPY3(U,a1,a2,a3,p1,p2,p3,n) {\
  switch (n & 0x3) {\
  case 3: *U++    += a1 * *p1++ + a2 * *p2++ + a3 * *p3++;\
  case 2: *U++    += a1 * *p1++ + a2 * *p2++ + a3 * *p3++;\
  case 1: *U++    += a1 * *p1++ + a2 * *p2++ + a3 * *p3++;\
  n -= 4;case 0:break;}\
  while (n>0) {U[0]+=a1*p1[0]+a2*p2[0]+a3*p3[0];U[1]+=a1*p1[1]+a2*p2[1]+a3*p3[1];U[2]+=a1*p1[2]+a2*p2[2]+a3*p3[2];U[3]+=a1*p1[3]+a2*p2[3]+a3*p3[3];U+=4;p1+=4;p2+=4;p3+=4;n-=4;}}
#define PetscKernelAXPY4(U,a1,a2,a3,a4,p1,p2,p3,p4,n) {\
  switch (n & 0x3) {\
  case 3: *U++    += a1 * *p1++ + a2 * *p2++ + a3 * *p3++ + a4 * *p4++;\
  case 2: *U++    += a1 * *p1++ + a2 * *p2++ + a3 * *p3++ + a4 * *p4++;\
  case 1: *U++    += a1 * *p1++ + a2 * *p2++ + a3 * *p3++ + a4 * *p4++;\
  n -= 4;case 0:break;}\
  while (n>0) {U[0]+=a1*p1[0]+a2*p2[0]+a3*p3[0]+a4*p4[0];U[1]+=a1*p1[1]+a2*p2[1]+a3*p3[1]+a4*p4[1];U[2]+=a1*p1[2]+a2*p2[2]+a3*p3[2]+a4*p4[2];U[3]+=a1*p1[3]+a2*p2[3]+a3*p3[3]+a4*p4[3];U+=4;p1+=4;p2+=4;p3+=4;p4+=4;n-=4;}}

#elif defined(PETSC_USE_WHILE_KERNELS)

#define PetscKernelAXPY(U,a1,p1,n)  {while (n--) *U++ += a1 * *p1++;}
#define PetscKernelAXPY2(U,a1,a2,p1,p2,n)  {while (n--) *U++ += a1 * *p1++ + a2 * *p2++;}
#define PetscKernelAXPY3(U,a1,a2,a3,p1,p2,p3,n) {while (n--) *U++ += a1 * *p1++ + a2 * *p2++ + a3 * *p3++;}
#define PetscKernelAXPY4(U,a1,a2,a3,a4,p1,p2,p3,p4,n) {while (n--) *U++ += a1 * *p1++ + a2 * *p2++ + a3 * *p3++ + a4 * *p4++;}

#elif defined(PETSC_USE_BLAS_KERNELS)

#define PetscKernelAXPY(U,a1,p1,n)  {PetscBLASInt one=1; PetscBLASInt nn = (PetscBLASInt) n; PetscStackCallBLAS("BLASaxpy",BLASaxpy_(&nn,&a1,p1,&one,U,&one));}
#define PetscKernelAXPY2(U,a1,a2,p1,p2,n){PetscKernelAXPY(U,a1,p1,n); PetscKernelAXPY(U,a2,p2,n);}
#define PetscKernelAXPY3(U,a1,a2,a3,p1,p2,p3,n){PetscKernelAXPY2(U,a1,a2,p1,p2,n); PetscKernelAXPY(U,a3,p3,n);}
#define PetscKernelAXPY4(U,a1,a2,a3,a4,p1,p2,p3,p4,n){PetscKernelAXPY2(U,a1,a2,p1,p2,n); PetscKernelAXPY2(U,a3,a4,p3,p4,n);}

#elif defined(PETSC_USE_FOR_KERNELS)

#define PetscKernelAXPY(U,a1,p1,n)  {PetscInt __i;PetscScalar __s1,__s2; \
  for (__i=0;__i<n-1;__i+=2){\
   __s1=a1*p1[__i];__s2=a1*p1[__i+1]; __s1+=U[__i];__s2+=U[__i+1];U[__i]=__s1;U[__i+1]=__s2;}\
  if (n & 0x1) U[__i] += a1 * p1[__i];}
#define PetscKernelAXPY2(U,a1,a2,p1,p2,n) {PetscInt __i; for (__i=0;__i<n;__i++) U[__i] += a1 * p1[__i] + a2 * p2[__i];}
#define PetscKernelAXPY3(U,a1,a2,a3,p1,p2,p3,n){PetscInt __i; for (__i=0;__i<n;__i++)U[__i]+=a1*p1[__i]+a2*p2[__i]+a3*p3[__i];}
#define PetscKernelAXPY4(U,a1,a2,a3,a4,p1,p2,p3,p4,n){PetscInt __i;for (__i=0;__i<n;__i++)U[__i]+=a1*p1[__i]+a2*p2[__i]+a3*p3[__i]+a4*p4[__i];}

#else

#define PetscKernelAXPY(U,a1,p1,n)  {PetscInt __i;PetscScalar _a1=a1; for (__i=0;__i<n;__i++)U[__i]+=_a1 * p1[__i];}
#define PetscKernelAXPY2(U,a1,a2,p1,p2,n) {PetscInt __i; for (__i=0;__i<n;__i++)U[__i] += a1 * p1[__i] + a2 * p2[__i];}
#define PetscKernelAXPY3(U,a1,a2,a3,p1,p2,p3,n){PetscInt __i; for (__i=0;__i<n;__i++)U[__i]+=a1*p1[__i]+a2*p2[__i]+a3*p3[__i];}
#define PetscKernelAXPY4(U,a1,a2,a3,a4,p1,p2,p3,p4,n){PetscInt __i;for (__i=0;__i<n;__i++)U[__i]+=a1*p1[__i]+a2*p2[__i]+a3*p3[__i]+a4*p4[__i];}

#endif

#endif
