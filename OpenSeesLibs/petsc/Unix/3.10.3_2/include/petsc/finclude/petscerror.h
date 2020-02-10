
!
!  Include file for Fortran error codes
!    These are also in include/petscerror.h
!
#if !defined (__PETSCERRORDEF_H)
#define __PETSCERRORDEF_H

#define PETSC_ERR_MEM              55
#define PETSC_ERR_SUP              56
#define PETSC_ERR_SUP_SYS          57
#define PETSC_ERR_ORDER            58
#define PETSC_ERR_SIG              59
#define PETSC_ERR_FP               72
#define PETSC_ERR_COR              74
#define PETSC_ERR_LIB              76
#define PETSC_ERR_PLIB             77
#define PETSC_ERR_MEMC             78
#define PETSC_ERR_CONV_FAILED      82
#define PETSC_ERR_USER             83
#define PETSC_ERR_SYS              88
#define PETSC_ERR_POINTER          70
#define PETSC_ERR_MPI_LIB_INCOMP   87

#define PETSC_ERR_ARG_SIZ          60
#define PETSC_ERR_ARG_IDN          61
#define PETSC_ERR_ARG_WRONG        62
#define PETSC_ERR_ARG_CORRUPT      64
#define PETSC_ERR_ARG_OUTOFRANGE   63
#define PETSC_ERR_ARG_BADPTR       68
#define PETSC_ERR_ARG_NOTSAMETYPE  69
#define PETSC_ERR_ARG_NOTSAMECOMM  80
#define PETSC_ERR_ARG_WRONGSTATE   73
#define PETSC_ERR_ARG_TYPENOTSET   89
#define PETSC_ERR_ARG_INCOMP       75
#define PETSC_ERR_ARG_NULL         85
#define PETSC_ERR_ARG_UNKNOWN_TYPE 86

#define PETSC_ERR_FILE_OPEN        65
#define PETSC_ERR_FILE_READ        66
#define PETSC_ERR_FILE_WRITE       67
#define PETSC_ERR_FILE_UNEXPECTED  79

#define PETSC_ERR_MAT_LU_ZRPVT     71
#define PETSC_ERR_MAT_CH_ZRPVT     81

#define PETSC_ERR_INT_OVERFLOW     84

#define PETSC_ERR_FLOP_COUNT       90
#define PETSC_ERR_NOT_CONVERGED    91
#define PETSC_ERR_MISSING_FACTOR   92
#define PETSC_ERR_OPT_OVERWRITE    93

#endif
