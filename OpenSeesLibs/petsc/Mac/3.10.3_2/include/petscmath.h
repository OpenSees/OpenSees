/*

      PETSc mathematics include file. Defines certain basic mathematical
    constants and functions for working with single, double, and quad precision
    floating point numbers as well as complex single and double.

    This file is included by petscsys.h and should not be used directly.

*/

#if !defined(__PETSCMATH_H)
#define __PETSCMATH_H
#include <math.h>

/*

     Defines operations that are different for complex and real numbers;
   note that one cannot mix the use of complex and real in the same
   PETSc program. All PETSc objects in one program are built around the object
   PetscScalar which is either always a real or a complex.

*/

#if defined(PETSC_USE_REAL_SINGLE)
#define MPIU_REAL   MPI_FLOAT
typedef float PetscReal;
#define PetscRoundReal(a)   roundf(a)
#define PetscSqrtReal(a)    sqrtf(a)
#define PetscExpReal(a)     expf(a)
#define PetscLogReal(a)     logf(a)
#define PetscLog10Real(a)   log10f(a)
#ifdef PETSC_HAVE_LOG2
#define PetscLog2Real(a)    log2f(a)
#endif
#define PetscSinReal(a)     sinf(a)
#define PetscCosReal(a)     cosf(a)
#define PetscTanReal(a)     tanf(a)
#define PetscAsinReal(a)    asinf(a)
#define PetscAcosReal(a)    acosf(a)
#define PetscAtanReal(a)    atanf(a)
#define PetscAtan2Real(a,b) atan2f(a,b)
#define PetscSinhReal(a)    sinhf(a)
#define PetscCoshReal(a)    coshf(a)
#define PetscTanhReal(a)    tanhf(a)
#define PetscPowReal(a,b)   powf(a,b)
#define PetscCeilReal(a)    ceilf(a)
#define PetscFloorReal(a)   floorf(a)
#define PetscFmodReal(a,b)  fmodf(a,b)
#define PetscTGamma(a)      tgammaf(a)
#elif defined(PETSC_USE_REAL_DOUBLE)
#define MPIU_REAL   MPI_DOUBLE
typedef double PetscReal;
#define PetscRoundReal(a)   round(a)
#define PetscSqrtReal(a)    sqrt(a)
#define PetscExpReal(a)     exp(a)
#define PetscLogReal(a)     log(a)
#define PetscLog10Real(a)   log10(a)
#ifdef PETSC_HAVE_LOG2
#define PetscLog2Real(a)    log2(a)
#endif
#define PetscSinReal(a)     sin(a)
#define PetscCosReal(a)     cos(a)
#define PetscTanReal(a)     tan(a)
#define PetscAsinReal(a)    asin(a)
#define PetscAcosReal(a)    acos(a)
#define PetscAtanReal(a)    atan(a)
#define PetscAtan2Real(a,b) atan2(a,b)
#define PetscSinhReal(a)    sinh(a)
#define PetscCoshReal(a)    cosh(a)
#define PetscTanhReal(a)    tanh(a)
#define PetscPowReal(a,b)   pow(a,b)
#define PetscCeilReal(a)    ceil(a)
#define PetscFloorReal(a)   floor(a)
#define PetscFmodReal(a,b)  fmod(a,b)
#define PetscTGamma(a)      tgamma(a)
#elif defined(PETSC_USE_REAL___FLOAT128)
#if defined(__cplusplus)
extern "C" {
#endif
#include <quadmath.h>
#if defined(__cplusplus)
}
#endif
PETSC_EXTERN MPI_Datatype MPIU___FLOAT128 PetscAttrMPITypeTag(__float128);
#define MPIU_REAL MPIU___FLOAT128
typedef __float128 PetscReal;
#define PetscRoundReal(a)   roundq(a)
#define PetscSqrtReal(a)    sqrtq(a)
#define PetscExpReal(a)     expq(a)
#define PetscLogReal(a)     logq(a)
#define PetscLog10Real(a)   log10q(a)
#ifdef PETSC_HAVE_LOG2
#define PetscLog2Real(a)    log2q(a)
#endif
#define PetscSinReal(a)     sinq(a)
#define PetscCosReal(a)     cosq(a)
#define PetscTanReal(a)     tanq(a)
#define PetscAsinReal(a)    asinq(a)
#define PetscAcosReal(a)    acosq(a)
#define PetscAtanReal(a)    atanq(a)
#define PetscAtan2Real(a,b) atan2q(a,b)
#define PetscSinhReal(a)    sinhq(a)
#define PetscCoshReal(a)    coshq(a)
#define PetscTanhReal(a)    tanhq(a)
#define PetscPowReal(a,b)   powq(a,b)
#define PetscCeilReal(a)    ceilq(a)
#define PetscFloorReal(a)   floorq(a)
#define PetscFmodReal(a,b)  fmodq(a,b)
#define PetscTGamma(a)      tgammaq(a)
#elif defined(PETSC_USE_REAL___FP16)
PETSC_EXTERN MPI_Datatype MPIU___FP16 PetscAttrMPITypeTag(__fp16);
#define MPIU_REAL MPIU___FP16
typedef __fp16 PetscReal;
#define PetscRoundReal(a)   roundf(a)
#define PetscSqrtReal(a)    sqrtf(a)
#define PetscExpReal(a)     expf(a)
#define PetscLogReal(a)     logf(a)
#define PetscLog10Real(a)   log10f(a)
#ifdef PETSC_HAVE_LOG2
#define PetscLog2Real(a)    log2f(a)
#endif
#define PetscSinReal(a)     sinf(a)
#define PetscCosReal(a)     cosf(a)
#define PetscTanReal(a)     tanf(a)
#define PetscAsinReal(a)    asinf(a)
#define PetscAcosReal(a)    acosf(a)
#define PetscAtanReal(a)    atanf(a)
#define PetscAtan2Real(a,b) atan2f(a,b)
#define PetscSinhReal(a)    sinhf(a)
#define PetscCoshReal(a)    coshf(a)
#define PetscTanhReal(a)    tanhf(a)
#define PetscPowReal(a,b)   powf(a,b)
#define PetscCeilReal(a)    ceilf(a)
#define PetscFloorReal(a)   floorf(a)
#define PetscFmodReal(a,b)  fmodf(a,b)
#define PetscTGamma(a)      tgammaf(a)
#endif /* PETSC_USE_REAL_* */

/*
    Complex number definitions
 */
#if defined(__cplusplus) && defined(PETSC_HAVE_CXX_COMPLEX) && !defined(PETSC_USE_REAL___FLOAT128)
#if !defined(PETSC_SKIP_COMPLEX)
#define PETSC_HAVE_COMPLEX 1
/* C++ support of complex number */
#if defined(PETSC_HAVE_VECCUDA) && __CUDACC_VER_MAJOR__ > 6
/* complex headers in thrust only available in CUDA 7.0 and above */
#define complexlib thrust
#include <thrust/complex.h>
#else
#define complexlib std
#include <complex>
#endif

#define PetscRealPartComplex(a)      (a).real()
#define PetscImaginaryPartComplex(a) (a).imag()
#define PetscAbsComplex(a)           complexlib::abs(a)
#define PetscConjComplex(a)          complexlib::conj(a)
#define PetscSqrtComplex(a)          complexlib::sqrt(a)
#define PetscPowComplex(a,b)         complexlib::pow(a,b)
#define PetscExpComplex(a)           complexlib::exp(a)
#define PetscLogComplex(a)           complexlib::log(a)
#define PetscSinComplex(a)           complexlib::sin(a)
#define PetscCosComplex(a)           complexlib::cos(a)
#define PetscAsinComplex(a)          complexlib::asin(a)
#define PetscAcosComplex(a)          complexlib::acos(a)
#if defined(PETSC_HAVE_TANCOMPLEX)
#define PetscTanComplex(a)           complexlib::tan(a)
#else
#define PetscTanComplex(a)           PetscSinComplex(a)/PetscCosComplex(a)
#endif
#define PetscSinhComplex(a)          complexlib::sinh(a)
#define PetscCoshComplex(a)          complexlib::cosh(a)
#if defined(PETSC_HAVE_TANHCOMPLEX)
#define PetscTanhComplex(a)          complexlib::tanh(a)
#else
#define PetscTanhComplex(a)          PetscSinhComplex(a)/PetscCoshComplex(a)
#endif

#if defined(PETSC_USE_REAL_SINGLE)
typedef complexlib::complex<float> PetscComplex;
#if defined(PETSC_USE_CXX_COMPLEX_FLOAT_WORKAROUND)
static inline PetscComplex operator+(const PetscComplex& lhs, const double& rhs) { return lhs + float(rhs); }
static inline PetscComplex operator+(const double& lhs, const PetscComplex& rhs) { return float(lhs) + rhs; }
static inline PetscComplex operator-(const PetscComplex& lhs, const double& rhs) { return lhs - float(rhs); }
static inline PetscComplex operator-(const double& lhs, const PetscComplex& rhs) { return float(lhs) - rhs; }
static inline PetscComplex operator*(const PetscComplex& lhs, const double& rhs) { return lhs * float(rhs); }
static inline PetscComplex operator*(const double& lhs, const PetscComplex& rhs) { return float(lhs) * rhs; }
static inline PetscComplex operator/(const PetscComplex& lhs, const double& rhs) { return lhs / float(rhs); }
static inline PetscComplex operator/(const double& lhs, const PetscComplex& rhs) { return float(lhs) / rhs; }
static inline bool operator==(const PetscComplex& lhs, const double& rhs) { return lhs.imag() == float(0) && lhs.real() == float(rhs); }
static inline bool operator==(const double& lhs, const PetscComplex& rhs) { return rhs.imag() == float(0) && rhs.real() == float(lhs); }
static inline bool operator!=(const PetscComplex& lhs, const double& rhs) { return lhs.imag() != float(0) || lhs.real() != float(rhs); }
static inline bool operator!=(const double& lhs, const PetscComplex& rhs) { return rhs.imag() != float(0) || rhs.real() != float(lhs); }
#endif  /* PETSC_USE_CXX_COMPLEX_FLOAT_WORKAROUND */
#elif defined(PETSC_USE_REAL_DOUBLE)
typedef complexlib::complex<double> PetscComplex;
#if defined(PETSC_USE_CXX_COMPLEX_FLOAT_WORKAROUND)
static inline PetscComplex operator+(const PetscComplex& lhs, const PetscInt& rhs) { return lhs + double(rhs); }
static inline PetscComplex operator+(const PetscInt& lhs, const PetscComplex& rhs) { return double(lhs) + rhs; }
static inline PetscComplex operator-(const PetscComplex& lhs, const PetscInt& rhs) { return lhs - double(rhs); }
static inline PetscComplex operator-(const PetscInt& lhs, const PetscComplex& rhs) { return double(lhs) - rhs; }
static inline PetscComplex operator*(const PetscComplex& lhs, const PetscInt& rhs) { return lhs * double(rhs); }
static inline PetscComplex operator*(const PetscInt& lhs, const PetscComplex& rhs) { return double(lhs) * rhs; }
static inline PetscComplex operator/(const PetscComplex& lhs, const PetscInt& rhs) { return lhs / double(rhs); }
static inline PetscComplex operator/(const PetscInt& lhs, const PetscComplex& rhs) { return double(lhs) / rhs; }
static inline bool operator==(const PetscComplex& lhs, const PetscInt& rhs) { return lhs.imag() == double(0) && lhs.real() == double(rhs); }
static inline bool operator==(const PetscInt& lhs, const PetscComplex& rhs) { return rhs.imag() == double(0) && rhs.real() == double(lhs); }
static inline bool operator!=(const PetscComplex& lhs, const PetscInt& rhs) { return lhs.imag() != double(0) || lhs.real() != double(rhs); }
static inline bool operator!=(const PetscInt& lhs, const PetscComplex& rhs) { return rhs.imag() != double(0) || rhs.real() != double(lhs); }
#endif  /* PETSC_USE_CXX_COMPLEX_FLOAT_WORKAROUND */
#elif defined(PETSC_USE_REAL___FLOAT128)
typedef complexlib::complex<__float128> PetscComplex; /* Notstandard and not expected to work, use __complex128 */
PETSC_EXTERN MPI_Datatype MPIU___COMPLEX128;
#endif  /* PETSC_USE_REAL_ */
#endif  /* ! PETSC_SKIP_COMPLEX */

#elif defined(PETSC_HAVE_C99_COMPLEX) && !defined(PETSC_USE_REAL___FP16)
#if !defined(PETSC_SKIP_COMPLEX)
#define PETSC_HAVE_COMPLEX 1
#include <complex.h>

#if defined(PETSC_USE_REAL_SINGLE) || defined(PETSC_USE_REAL___FP16)
typedef float _Complex PetscComplex;

#define PetscRealPartComplex(a)      crealf(a)
#define PetscImaginaryPartComplex(a) cimagf(a)
#define PetscAbsComplex(a)           cabsf(a)
#define PetscConjComplex(a)          conjf(a)
#define PetscSqrtComplex(a)          csqrtf(a)
#define PetscPowComplex(a,b)         cpowf(a,b)
#define PetscExpComplex(a)           cexpf(a)
#define PetscLogComplex(a)           clogf(a)
#define PetscSinComplex(a)           csinf(a)
#define PetscCosComplex(a)           ccosf(a)
#define PetscAsinComplex(a)          casinf(a)
#define PetscAcosComplex(a)          cacosf(a)
#define PetscTanComplex(a)           ctanf(a)
#define PetscSinhComplex(a)          csinhf(a)
#define PetscCoshComplex(a)          ccoshf(a)
#define PetscTanhComplex(a)          ctanhf(a)

#elif defined(PETSC_USE_REAL_DOUBLE)
typedef double _Complex PetscComplex;

#define PetscRealPartComplex(a)      creal(a)
#define PetscImaginaryPartComplex(a) cimag(a)
#define PetscAbsComplex(a)           cabs(a)
#define PetscConjComplex(a)          conj(a)
#define PetscSqrtComplex(a)          csqrt(a)
#define PetscPowComplex(a,b)         cpow(a,b)
#define PetscExpComplex(a)           cexp(a)
#define PetscLogComplex(a)           clog(a)
#define PetscSinComplex(a)           csin(a)
#define PetscCosComplex(a)           ccos(a)
#define PetscAsinComplex(a)          casin(a)
#define PetscAcosComplex(a)          cacos(a)
#define PetscTanComplex(a)           ctan(a)
#define PetscSinhComplex(a)          csinh(a)
#define PetscCoshComplex(a)          ccosh(a)
#define PetscTanhComplex(a)          ctanh(a)

#elif defined(PETSC_USE_REAL___FLOAT128)
typedef __complex128 PetscComplex;
PETSC_EXTERN MPI_Datatype MPIU___COMPLEX128 PetscAttrMPITypeTag(__complex128);

#define PetscRealPartComplex(a)      crealq(a)
#define PetscImaginaryPartComplex(a) cimagq(a)
#define PetscAbsComplex(a)           cabsq(a)
#define PetscConjComplex(a)          conjq(a)
#define PetscSqrtComplex(a)          csqrtq(a)
#define PetscPowComplex(a,b)         cpowq(a,b)
#define PetscExpComplex(a)           cexpq(a)
#define PetscLogComplex(a)           clogq(a)
#define PetscSinComplex(a)           csinq(a)
#define PetscCosComplex(a)           ccosq(a)
#define PetscAsinComplex(a)          casinq(a)
#define PetscAcosComplex(a)          cacosq(a)
#define PetscTanComplex(a)           ctanq(a)
#define PetscSinhComplex(a)          csinhq(a)
#define PetscCoshComplex(a)          ccoshq(a)
#define PetscTanhComplex(a)          ctanhq(a)

#endif /* PETSC_USE_REAL_* */
#elif (defined(PETSC_USE_COMPLEX) && !defined(PETSC_SKIP_COMPLEX))
#error "PETSc was configured --with-scalar-type=complex, but a language-appropriate complex library is not available"
#endif /* !PETSC_SKIP_COMPLEX */
#endif /* (__cplusplus && PETSC_HAVE_CXX_COMPLEX) else-if (!__cplusplus && PETSC_HAVE_C99_COMPLEX) */

#if defined(PETSC_HAVE_COMPLEX)
#if defined(PETSC_HAVE_MPI_C_DOUBLE_COMPLEX)
#define MPIU_C_DOUBLE_COMPLEX MPI_C_DOUBLE_COMPLEX
#define MPIU_C_COMPLEX MPI_C_COMPLEX
#else
# if defined(__cplusplus) && defined(PETSC_HAVE_CXX_COMPLEX) && !defined(PETSC_USE_REAL___FLOAT128)
  typedef complexlib::complex<double> petsc_mpiu_c_double_complex;
  typedef complexlib::complex<float> petsc_mpiu_c_complex;
# elif !defined(__cplusplus) && defined(PETSC_HAVE_C99_COMPLEX)
  typedef double _Complex petsc_mpiu_c_double_complex;
  typedef float _Complex petsc_mpiu_c_complex;
# else
  typedef struct {double real,imag;} petsc_mpiu_c_double_complex;
  typedef struct {float real,imag;} petsc_mpiu_c_complex;
# endif
PETSC_EXTERN MPI_Datatype MPIU_C_DOUBLE_COMPLEX PetscAttrMPITypeTagLayoutCompatible(petsc_mpiu_c_double_complex);
PETSC_EXTERN MPI_Datatype MPIU_C_COMPLEX PetscAttrMPITypeTagLayoutCompatible(petsc_mpiu_c_complex);
#endif /* PETSC_HAVE_MPI_C_DOUBLE_COMPLEX */
#endif /* PETSC_HAVE_COMPLEX */

#if defined(PETSC_HAVE_COMPLEX)
#  if defined(PETSC_USE_REAL_SINGLE) || defined(PETSC_USE_REAL___FP16)
#    define MPIU_COMPLEX MPIU_C_COMPLEX
#  elif defined(PETSC_USE_REAL_DOUBLE)
#    define MPIU_COMPLEX MPIU_C_DOUBLE_COMPLEX
#  elif defined(PETSC_USE_REAL___FLOAT128)
#    define MPIU_COMPLEX MPIU___COMPLEX128
#  endif /* PETSC_USE_REAL_* */
#endif

#if (defined(PETSC_USE_COMPLEX) && !defined(PETSC_SKIP_COMPLEX))
typedef PetscComplex PetscScalar;

/*MC
   PetscRealPart - Returns the real part of a PetscScalar

   Synopsis:
   #include <petscmath.h>
   PetscScalar PetscRealPart(PetscScalar v)

   Not Collective

   Input Parameter:
.  v - value to find the real part of

   Level: beginner

.seealso: PetscScalar, PetscImaginaryPart(), PetscMax(), PetscClipInterval(), PetscAbsInt(), PetscAbsReal(), PetscSqr()

M*/
#define PetscRealPart(a)      PetscRealPartComplex(a)

/*MC
   PetscImaginaryPart - Returns the imaginary part of a PetscScalar

   Synopsis:
   #include <petscmath.h>
   PetscScalar PetscImaginaryPart(PetscScalar v)

   Not Collective

   Input Parameter:
.  v - value to find the imaginary part of

   Level: beginner

   Notes:
       If PETSc was configured for real numbers then this always returns the value 0

.seealso: PetscScalar, PetscRealPart(), PetscMax(), PetscClipInterval(), PetscAbsInt(), PetscAbsReal(), PetscSqr()

M*/
#define PetscImaginaryPart(a) PetscImaginaryPartComplex(a)

#define PetscAbsScalar(a)     PetscAbsComplex(a)
#define PetscConj(a)          PetscConjComplex(a)
#define PetscSqrtScalar(a)    PetscSqrtComplex(a)
#define PetscPowScalar(a,b)   PetscPowComplex(a,b)
#define PetscExpScalar(a)     PetscExpComplex(a)
#define PetscLogScalar(a)     PetscLogComplex(a)
#define PetscSinScalar(a)     PetscSinComplex(a)
#define PetscCosScalar(a)     PetscCosComplex(a)
#define PetscAsinScalar(a)    PetscAsinComplex(a)
#define PetscAcosScalar(a)    PetscAcosComplex(a)
#define PetscTanScalar(a)     PetscTanComplex(a)
#define PetscSinhScalar(a)    PetscSinhComplex(a)
#define PetscCoshScalar(a)    PetscCoshComplex(a)
#define PetscTanhScalar(a)    PetscTanhComplex(a)
#define MPIU_SCALAR MPIU_COMPLEX

/*
    real number definitions
 */
#else /* PETSC_USE_COMPLEX */
typedef PetscReal PetscScalar;
#define MPIU_SCALAR MPIU_REAL

#define PetscRealPart(a)      (a)
#define PetscImaginaryPart(a) ((PetscReal)0.)
PETSC_STATIC_INLINE PetscReal PetscAbsScalar(PetscScalar a) {return a < 0.0 ? -a : a;}
#define PetscConj(a)          (a)
#if !defined(PETSC_USE_REAL___FLOAT128)  && !defined(PETSC_USE_REAL___FP16)
#define PetscSqrtScalar(a)    sqrt(a)
#define PetscPowScalar(a,b)   pow(a,b)
#define PetscExpScalar(a)     exp(a)
#define PetscLogScalar(a)     log(a)
#define PetscSinScalar(a)     sin(a)
#define PetscCosScalar(a)     cos(a)
#define PetscAsinScalar(a)    asin(a)
#define PetscAcosScalar(a)    acos(a)
#define PetscTanScalar(a)     tan(a)
#define PetscSinhScalar(a)    sinh(a)
#define PetscCoshScalar(a)    cosh(a)
#define PetscTanhScalar(a)    tanh(a)
#elif defined(PETSC_USE_REAL___FP16)
#define PetscSqrtScalar(a)    sqrtf(a)
#define PetscPowScalar(a,b)   powf(a,b)
#define PetscExpScalar(a)     expf(a)
#define PetscLogScalar(a)     logf(a)
#define PetscSinScalar(a)     sinf(a)
#define PetscCosScalar(a)     cosf(a)
#define PetscAsinScalar(a)    asinf(a)
#define PetscAcosScalar(a)    acosf(a)
#define PetscTanScalar(a)     tanf(a)
#define PetscSinhScalar(a)    sinhf(a)
#define PetscCoshScalar(a)    coshf(a)
#define PetscTanhScalar(a)    tanhf(a)
#else /* PETSC_USE_REAL___FLOAT128 */
#define PetscSqrtScalar(a)    sqrtq(a)
#define PetscPowScalar(a,b)   powq(a,b)
#define PetscExpScalar(a)     expq(a)
#define PetscLogScalar(a)     logq(a)
#define PetscSinScalar(a)     sinq(a)
#define PetscCosScalar(a)     cosq(a)
#define PetscAsinScalar(a)    asinq(a)
#define PetscAcosScalar(a)    acosq(a)
#define PetscTanScalar(a)     tanq(a)
#define PetscSinhScalar(a)    sinhq(a)
#define PetscCoshScalar(a)    coshq(a)
#define PetscTanhScalar(a)    tanhq(a)
#endif /* PETSC_USE_REAL___FLOAT128 */

#endif /* PETSC_USE_COMPLEX */

#define PetscSign(a) (((a) >= 0) ? ((a) == 0 ? 0 : 1) : -1)
#define PetscSignReal(a) (((a) >= 0.0) ? ((a) == 0.0 ? 0.0 : 1.0) : -1.0)
#define PetscAbs(a)  (((a) >= 0) ? (a) : (-(a)))

/* --------------------------------------------------------------------------*/

/*
   Certain objects may be created using either single or double precision.
   This is currently not used.
*/
typedef enum { PETSC_SCALAR_DOUBLE,PETSC_SCALAR_SINGLE, PETSC_SCALAR_LONG_DOUBLE, PETSC_SCALAR_HALF } PetscScalarPrecision;

#if defined(PETSC_HAVE_COMPLEX)
/* PETSC_i is the imaginary number, i */
PETSC_EXTERN PetscComplex PETSC_i;

/* Try to do the right thing for complex number construction: see

  http://www.open-std.org/jtc1/sc22/wg14/www/docs/n1464.htm

  for details
*/
PETSC_STATIC_INLINE PetscComplex PetscCMPLX(PetscReal x, PetscReal y)
{
#if   defined(__cplusplus) && defined(PETSC_HAVE_CXX_COMPLEX) && !defined(PETSC_USE_REAL___FLOAT128)
  return PetscComplex(x,y);
#elif defined(_Imaginary_I)
  return x + y * _Imaginary_I;
#else
  { /* In both C99 and C11 (ISO/IEC 9899, Section 6.2.5),

       "For each floating type there is a corresponding real type, which is always a real floating
       type. For real floating types, it is the same type. For complex types, it is the type given
       by deleting the keyword _Complex from the type name."

       So type punning should be portable. */
    union { PetscComplex z; PetscReal f[2]; } uz;

    uz.f[0] = x;
    uz.f[1] = y;
    return uz.z;
  }
#endif
}
#endif


/*MC
   PetscMin - Returns minimum of two numbers

   Synopsis:
   #include <petscmath.h>
   type PetscMin(type v1,type v2)

   Not Collective

   Input Parameter:
+  v1 - first value to find minimum of
-  v2 - second value to find minimum of

   Notes:
    type can be integer or floating point value

   Level: beginner

.seealso: PetscMax(), PetscClipInterval(), PetscAbsInt(), PetscAbsReal(), PetscSqr()

M*/
#define PetscMin(a,b)   (((a)<(b)) ?  (a) : (b))

/*MC
   PetscMax - Returns maxium of two numbers

   Synopsis:
   #include <petscmath.h>
   type max PetscMax(type v1,type v2)

   Not Collective

   Input Parameter:
+  v1 - first value to find maximum of
-  v2 - second value to find maximum of

   Notes:
    type can be integer or floating point value

   Level: beginner

.seealso: PetscMin(), PetscClipInterval(), PetscAbsInt(), PetscAbsReal(), PetscSqr()

M*/
#define PetscMax(a,b)   (((a)<(b)) ?  (b) : (a))

/*MC
   PetscClipInterval - Returns a number clipped to be within an interval

   Synopsis:
   #include <petscmath.h>
   type clip PetscClipInterval(type x,type a,type b)

   Not Collective

   Input Parameter:
+  x - value to use if within interval (a,b)
.  a - lower end of interval
-  b - upper end of interval

   Notes:
    type can be integer or floating point value

   Level: beginner

.seealso: PetscMin(), PetscMax(), PetscAbsInt(), PetscAbsReal(), PetscSqr()

M*/
#define PetscClipInterval(x,a,b)   (PetscMax((a),PetscMin((x),(b))))

/*MC
   PetscAbsInt - Returns the absolute value of an integer

   Synopsis:
   #include <petscmath.h>
   int abs PetscAbsInt(int v1)

   Not Collective

   Input Parameter:
.   v1 - the integer

   Level: beginner

.seealso: PetscMax(), PetscMin(), PetscAbsReal(), PetscSqr()

M*/
#define PetscAbsInt(a)  (((a)<0)   ? (-(a)) : (a))

/*MC
   PetscAbsReal - Returns the absolute value of an real number

   Synopsis:
   #include <petscmath.h>
   Real abs PetscAbsReal(PetscReal v1)

   Not Collective

   Input Parameter:
.   v1 - the double


   Level: beginner

.seealso: PetscMax(), PetscMin(), PetscAbsInt(), PetscSqr()

M*/
#if defined(PETSC_USE_REAL_SINGLE)
#define PetscAbsReal(a) fabsf(a)
#elif defined(PETSC_USE_REAL_DOUBLE)
#define PetscAbsReal(a) fabs(a)
#elif defined(PETSC_USE_REAL___FLOAT128)
#define PetscAbsReal(a) fabsq(a)
#elif defined(PETSC_USE_REAL___FP16)
#define PetscAbsReal(a) fabsf(a)
#endif

/*MC
   PetscSqr - Returns the square of a number

   Synopsis:
   #include <petscmath.h>
   type sqr PetscSqr(type v1)

   Not Collective

   Input Parameter:
.   v1 - the value

   Notes:
    type can be integer or floating point value

   Level: beginner

.seealso: PetscMax(), PetscMin(), PetscAbsInt(), PetscAbsReal()

M*/
#define PetscSqr(a)     ((a)*(a))

/* ----------------------------------------------------------------------------*/

#if defined(PETSC_USE_REAL_SINGLE)
#define PetscRealConstant(constant) constant##F
#elif defined(PETSC_USE_REAL___FLOAT128)
#define PetscRealConstant(constant) constant##Q
#else
#define PetscRealConstant(constant) constant
#endif

/*
     Basic constants
*/
#define PETSC_PI   PetscRealConstant(3.1415926535897932384626433832795029)
#define PETSC_PHI  PetscRealConstant(1.6180339887498948482045868343656381)

#if !defined(PETSC_USE_64BIT_INDICES)
#define PETSC_MAX_INT            2147483647
#define PETSC_MIN_INT            (-PETSC_MAX_INT - 1)
#else
#define PETSC_MAX_INT            9223372036854775807L
#define PETSC_MIN_INT            (-PETSC_MAX_INT - 1)
#endif

#if defined(PETSC_USE_REAL_SINGLE)
#  define PETSC_MAX_REAL                3.40282346638528860e+38F
#  define PETSC_MIN_REAL                (-PETSC_MAX_REAL)
#  define PETSC_MACHINE_EPSILON         1.19209290e-07F
#  define PETSC_SQRT_MACHINE_EPSILON    3.45266983e-04F
#  define PETSC_SMALL                   1.e-5F
#elif defined(PETSC_USE_REAL_DOUBLE)
#  define PETSC_MAX_REAL                1.7976931348623157e+308
#  define PETSC_MIN_REAL                (-PETSC_MAX_REAL)
#  define PETSC_MACHINE_EPSILON         2.2204460492503131e-16
#  define PETSC_SQRT_MACHINE_EPSILON    1.490116119384766e-08
#  define PETSC_SMALL                   1.e-10
#elif defined(PETSC_USE_REAL___FLOAT128)
#  define PETSC_MAX_REAL                FLT128_MAX
#  define PETSC_MIN_REAL                (-FLT128_MAX)
#  define PETSC_MACHINE_EPSILON         FLT128_EPSILON
#  define PETSC_SQRT_MACHINE_EPSILON    1.38777878078144567552953958511352539e-17Q
#  define PETSC_SMALL                   1.e-20Q
#elif defined(PETSC_USE_REAL___FP16)  /* maybe should use single precision values for these? */
#  define PETSC_MAX_REAL                65504.
#  define PETSC_MIN_REAL                (-PETSC_MAX_REAL)
#  define PETSC_MACHINE_EPSILON         .00097656
#  define PETSC_SQRT_MACHINE_EPSILON    .0312
#  define PETSC_SMALL                   5.e-3
#endif

#define PETSC_INFINITY               (PETSC_MAX_REAL/4)
#define PETSC_NINFINITY              (-PETSC_INFINITY)

PETSC_EXTERN PetscBool PetscIsInfReal(PetscReal);
PETSC_EXTERN PetscBool PetscIsNanReal(PetscReal);
PETSC_EXTERN PetscBool PetscIsNormalReal(PetscReal);
PETSC_STATIC_INLINE PetscBool PetscIsInfOrNanReal(PetscReal v) {return PetscIsInfReal(v) || PetscIsNanReal(v) ? PETSC_TRUE : PETSC_FALSE;}
PETSC_STATIC_INLINE PetscBool PetscIsInfScalar(PetscScalar v) {return PetscIsInfReal(PetscAbsScalar(v));}
PETSC_STATIC_INLINE PetscBool PetscIsNanScalar(PetscScalar v) {return PetscIsNanReal(PetscAbsScalar(v));}
PETSC_STATIC_INLINE PetscBool PetscIsInfOrNanScalar(PetscScalar v) {return PetscIsInfOrNanReal(PetscAbsScalar(v));}
PETSC_STATIC_INLINE PetscBool PetscIsNormalScalar(PetscScalar v) {return PetscIsNormalReal(PetscAbsScalar(v));}

PETSC_EXTERN PetscBool PetscEqualReal(PetscReal,PetscReal);
PETSC_EXTERN PetscBool PetscEqualScalar(PetscScalar,PetscScalar);

/*
    These macros are currently hardwired to match the regular data types, so there is no support for a different
    MatScalar from PetscScalar. We left the MatScalar in the source just in case we use it again.
 */
#define MPIU_MATSCALAR MPIU_SCALAR
typedef PetscScalar MatScalar;
typedef PetscReal MatReal;

struct petsc_mpiu_2scalar {PetscScalar a,b;};
PETSC_EXTERN MPI_Datatype MPIU_2SCALAR PetscAttrMPITypeTagLayoutCompatible(struct petsc_mpiu_2scalar);
#if defined(PETSC_USE_64BIT_INDICES) || !defined(MPI_2INT)
struct petsc_mpiu_2int {PetscInt a,b;};
PETSC_EXTERN MPI_Datatype MPIU_2INT PetscAttrMPITypeTagLayoutCompatible(struct petsc_mpiu_2int);
#else
#define MPIU_2INT MPI_2INT
#endif

PETSC_STATIC_INLINE PetscInt PetscPowInt(PetscInt base,PetscInt power)
{
  PetscInt result = 1;
  while (power) {
    if (power & 1) result *= base;
    power >>= 1;
    base *= base;
  }
  return result;
}

PETSC_STATIC_INLINE PetscReal PetscPowRealInt(PetscReal base,PetscInt power)
{
  PetscReal result = 1;
  if (power < 0) {
    power = -power;
    base  = ((PetscReal)1)/base;
  }
  while (power) {
    if (power & 1) result *= base;
    power >>= 1;
    base *= base;
  }
  return result;
}

PETSC_STATIC_INLINE PetscScalar PetscPowScalarInt(PetscScalar base,PetscInt power)
{
  PetscScalar result = 1;
  if (power < 0) {
    power = -power;
    base  = ((PetscReal)1)/base;
  }
  while (power) {
    if (power & 1) result *= base;
    power >>= 1;
    base *= base;
  }
  return result;
}

PETSC_STATIC_INLINE PetscScalar PetscPowScalarReal(PetscScalar base,PetscReal power)
{
  PetscScalar cpower = power;
  return PetscPowScalar(base,cpower);
}

#ifndef PETSC_HAVE_LOG2
PETSC_STATIC_INLINE PetscReal PetscLog2Real(PetscReal n)
{
  return PetscLogReal(n)/PetscLogReal(2.0);
}
#endif
#endif
