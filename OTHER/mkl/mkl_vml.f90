! file: mkl_vml.fi
!===============================================================================
! Copyright 2006-2018 Intel Corporation All Rights Reserved.
!
! The source code,  information  and material  ("Material") contained  herein is
! owned by Intel Corporation or its  suppliers or licensors,  and  title to such
! Material remains with Intel  Corporation or its  suppliers or  licensors.  The
! Material  contains  proprietary  information  of  Intel or  its suppliers  and
! licensors.  The Material is protected by  worldwide copyright  laws and treaty
! provisions.  No part  of  the  Material   may  be  used,  copied,  reproduced,
! modified, published,  uploaded, posted, transmitted,  distributed or disclosed
! in any way without Intel's prior express written permission.  No license under
! any patent,  copyright or other  intellectual property rights  in the Material
! is granted to  or  conferred  upon  you,  either   expressly,  by implication,
! inducement,  estoppel  or  otherwise.  Any  license   under such  intellectual
! property rights must be express and approved by Intel in writing.
!
! Unless otherwise agreed by Intel in writing,  you may not remove or alter this
! notice or  any  other  notice   embedded  in  Materials  by  Intel  or Intel's
! suppliers or licensors in any way.
!===============================================================================

!++
!  Fortran 90 VML interface.
!--

!++
!  PARAMETER DEFINITIONS
!  Parameter definitions for VML mode and VML error status.
!
!  VML mode controls VML function accuracy, floating-point settings (rounding
!  mode and precision) and VML error handling options. Default VML mode is
!  VML_HA | VML_ERRMODE_DEFAULT, i.e. VML high accuracy functions are
!  called, and current floating-point precision and the rounding mode is used.
!
!  Error status macros are used for error classification.
!
!  NOTE: A | B means bitwise OR operation with A and B
!--

!  VML FUNCTION ACCURACY CONTROL
!  VML_HA - when VML_HA is set, high accuracy VML functions are called
!  VML_LA - when VML_LA is set, low accuracy VML functions are called
!  VML_EP - when VML_EP is set, enhanced performance VML functions are called
!
!  NOTE: VML_HA, VML_LA and VML_EP must not be used in combination
      INTEGER(KIND=4) VML_LA
      INTEGER(KIND=4) VML_HA
      INTEGER(KIND=4) VML_EP
      PARAMETER (VML_LA = Z"00000001")
      PARAMETER (VML_HA = Z"00000002")
      PARAMETER (VML_EP = Z"00000003")

!  SETTING OPTIMAL FLOATING-POINT PRECISION AND ROUNDING MODE
!  Definitions below are to set optimal floating-point control word
!  (precision and rounding mode).
!
!  For their correct work, VML functions change floating-point precision and
!  rounding mode (if necessary). Since control word changing is typically
!  expensive operation, it is recommended to set precision and rounding mode
!  to optimal values before VML function calls.
!
!  VML_FLOAT_CONSISTENT  - use this value if the calls are typically to single
!                          precision VML functions
!  VML_DOUBLE_CONSISTENT - use this value if the calls are typically to double
!                          precision VML functions
!  VML_RESTORE           - restore original floating-point precision and
!                          rounding mode
!  VML_DEFAULT_PRECISION - use default (current) floating-point precision and
!                          rounding mode
!  NOTE: VML_FLOAT_CONSISTENT, VML_DOUBLE_CONSISTENT, VML_RESTORE and
!        VML_DEFAULT_PRECISION must not be used in combination
      INTEGER(KIND=4) VML_DEFAULT_PRECISION
      INTEGER(KIND=4) VML_FLOAT_CONSISTENT
      INTEGER(KIND=4) VML_DOUBLE_CONSISTENT
      INTEGER(KIND=4) VML_RESTORE
      PARAMETER (VML_DEFAULT_PRECISION = Z"00000000")
      PARAMETER (VML_FLOAT_CONSISTENT  = Z"00000010")
      PARAMETER (VML_DOUBLE_CONSISTENT = Z"00000020")
      PARAMETER (VML_RESTORE           = Z"00000030")

!  VML ERROR HANDLING CONTROL
!  Macros below are used to control VML error handler.
!
!  VML_ERRMODE_IGNORE   - ignore errors
!  VML_ERRMODE_ERRNO    - errno variable is set on error
!  VML_ERRMODE_STDERR   - error description text is written to stderr on error
!  VML_ERRMODE_EXCEPT   - exception is raised on error
!  VML_ERRMODE_CALLBACK - user's error handler function is called on error
!  VML_ERRMODE_NOERR    - ignore errors and do not update status
!  VML_ERRMODE_DEFAULT  - errno variable is set, exceptions are raised and
!                         user's error handler is called on error
!  NOTE: VML_ERRMODE_IGNORE must not be used in combination with
!        VML_ERRMODE_ERRNO, VML_ERRMODE_STDERR, VML_ERRMODE_EXCEPT,
!        VML_ERRMODE_CALLBACK and VML_ERRMODE_DEFAULT.
!  NOTE: VML_ERRMODE_NOERR must not be used in combination with any
!        other VML_ERRMODE setting.
      INTEGER(KIND=4) VML_ERRMODE_IGNORE
      INTEGER(KIND=4) VML_ERRMODE_ERRNO
      INTEGER(KIND=4) VML_ERRMODE_STDERR
      INTEGER(KIND=4) VML_ERRMODE_EXCEPT
      INTEGER(KIND=4) VML_ERRMODE_CALLBACK
      INTEGER(KIND=4) VML_ERRMODE_NOERR
      INTEGER(KIND=4) VML_ERRMODE_DEFAULT
      PARAMETER (VML_ERRMODE_IGNORE   = Z"00000100")
      PARAMETER (VML_ERRMODE_ERRNO    = Z"00000200")
      PARAMETER (VML_ERRMODE_STDERR   = Z"00000400")
      PARAMETER (VML_ERRMODE_EXCEPT   = Z"00000800")
      PARAMETER (VML_ERRMODE_CALLBACK = Z"00001000")
      PARAMETER (VML_ERRMODE_NOERR    = Z"00002000")
      PARAMETER (VML_ERRMODE_DEFAULT  = IOR(VML_ERRMODE_ERRNO,          &
     &           IOR(VML_ERRMODE_CALLBACK,VML_ERRMODE_EXCEPT)))

!  ACCURACY, FLOATING-POINT CONTROL AND ERROR HANDLING MASKS
!  Accuracy, floating-point and error handling control are packed in
!  the VML mode variable. Macros below are useful to extract accuracy and/or
!  floating-point control and/or error handling control settings.
!
!  VML_ACCURACY_MASK           - extract accuracy bits
!  VML_FPUMODE_MASK            - extract floating-point control bits
!  VML_ERRMODE_MASK            - extract error handling control bits
!                                (including error callback bits)
!  VML_ERRMODE_STDHANDLER_MASK - extract error handling control bits
!                                (not including error callback bits)
!  VML_ERRMODE_CALLBACK_MASK   - extract error callback bits
      INTEGER(KIND=4) VML_ACCURACY_MASK
      INTEGER(KIND=4) VML_FPUMODE_MASK
      INTEGER(KIND=4) VML_ERRMODE_MASK
      INTEGER(KIND=4) VML_ERRMODE_STDHANDLER_MASK
      INTEGER(KIND=4) VML_ERRMODE_CALLBACK_MASK
      PARAMETER (VML_ACCURACY_MASK = Z"0000000f")
      PARAMETER (VML_FPUMODE_MASK  = Z"000000f0")
      PARAMETER (VML_ERRMODE_MASK  = Z"0000ff00")
      PARAMETER (VML_ERRMODE_STDHANDLER_MASK = Z"00000f00")
      PARAMETER (VML_ERRMODE_CALLBACK_MASK = Z"0000f000")

!  ERROR STATUS PARAMETER DEFINITIONS
!  VML_STATUS_OK        - no errors
!  VML_STATUS_BADSIZE   - array dimension is not positive
!  VML_STATUS_BADMEM    - invalid pointer passed
!  VML_STATUS_ERRDOM    - at least one of arguments is out of function domain
!  VML_STATUS_SING      - at least one of arguments caused singularity
!  VML_STATUS_OVERFLOW  - at least one of arguments caused overflow
!  VML_STATUS_UNDERFLOW - at least one of arguments caused underflow
!  VML_STATUS_ACCURACYWARNING - function doesn't support set accuracy mode,
!                               lower accuracy mode was used instead
      INTEGER(KIND=4) VML_STATUS_OK
      INTEGER(KIND=4) VML_STATUS_BADSIZE
      INTEGER(KIND=4) VML_STATUS_BADMEM
      INTEGER(KIND=4) VML_STATUS_ERRDOM
      INTEGER(KIND=4) VML_STATUS_SING
      INTEGER(KIND=4) VML_STATUS_OVERFLOW
      INTEGER(KIND=4) VML_STATUS_UNDERFLOW
      INTEGER(KIND=4) VML_STATUS_ACCURACYWARNING
      PARAMETER (VML_STATUS_OK        = 0)
      PARAMETER (VML_STATUS_BADSIZE   = -1)
      PARAMETER (VML_STATUS_BADMEM    = -2)
      PARAMETER (VML_STATUS_ERRDOM    = 1)
      PARAMETER (VML_STATUS_SING      = 2)
      PARAMETER (VML_STATUS_OVERFLOW  = 3)
      PARAMETER (VML_STATUS_UNDERFLOW = 4)
      PARAMETER (VML_STATUS_ACCURACYWARNING = 1000)

!++
!  TYPE DEFINITIONS
!--

!  ERROR CALLBACK CONTEXT.
!  Error callback context structure is used in a user's error callback
!  function with the following interface:
!
!  Error callback context fields:
!  ICODE        - error status
!  IINDEX       - index of bad argument
!  DBA1         - 1-st argument value, at which error occured
!  DBA2         - 2-nd argument value, at which error occured
!                 (2-argument functions only)
!  DBR1         - 1-st resulting value
!  DBR2         - 2-nd resulting value (2-result functions only)
!  CFUNCNAME    - function name, for which error occured
!  IFUNCNAMELEN - length of function name
      !dec$ options /warn=noalignment
      TYPE ERROR_STRUCTURE
            SEQUENCE
            INTEGER(KIND=4) ICODE
            INTEGER(KIND=4) IINDEX
            REAL(KIND=8)    DBA1
            REAL(KIND=8)    DBA2
            REAL(KIND=8)    DBR1
            REAL(KIND=8)    DBR2
            CHARACTER(64)   CFUNCNAME
            INTEGER(KIND=4) IFUNCNAMELEN
            REAL(KIND=8)    DBA1IM
            REAL(KIND=8)    DBA2IM
            REAL(KIND=8)    DBR1IM
            REAL(KIND=8)    DBR2IM
      END TYPE ERROR_STRUCTURE
      !dec$ end options

!++
!  VML ELEMENTARY FUNCTION INTERFACES.
!--

!  Absolute value: r[i] = |a[i]|

      INTERFACE
        SUBROUTINE vsabs(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdabs(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsabs(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdabs(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex absolute value: r[i] = |a[i]|

      INTERFACE
        SUBROUTINE vcabs(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          REAL(KIND=4),INTENT(OUT)    :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzabs(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          REAL(KIND=8),INTENT(OUT)    :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcabs(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          REAL(KIND=4),INTENT(OUT)    :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzabs(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          REAL(KIND=8),INTENT(OUT)    :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Argument of complex value: r[i] = carg(a[i])

      INTERFACE
        SUBROUTINE vcarg(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          REAL(KIND=4),INTENT(OUT)    :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzarg(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          REAL(KIND=8),INTENT(OUT)    :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcarg(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          REAL(KIND=4),INTENT(OUT)    :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzarg(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          REAL(KIND=8),INTENT(OUT)    :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Addition: r[i] = a[i] + b[i]

      INTERFACE
        SUBROUTINE vsadd(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdadd(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsadd(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdadd(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex addition: r[i] = a[i] + b[i]

      INTERFACE
        SUBROUTINE vcadd(n,a,b,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n),b(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzadd(n,a,b,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n),b(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcadd(n,a,b,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n),b(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzadd(n,a,b,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n),b(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Subtraction: r[i] = a[i] - b[i]

      INTERFACE
        SUBROUTINE vssub(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdsub(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmssub(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdsub(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex subtraction: r[i] = a[i] - b[i]

      INTERFACE
        SUBROUTINE vcsub(n,a,b,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n),b(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzsub(n,a,b,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n),b(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcsub(n,a,b,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n),b(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzsub(n,a,b,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n),b(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Reciprocal: r[i] = 1.0 / a[i]

      INTERFACE
        SUBROUTINE vsinv(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdinv(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsinv(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdinv(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Square Root: r[i] = a[i]^0.5

      INTERFACE
        SUBROUTINE vssqrt(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdsqrt(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmssqrt(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdsqrt(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex Square Root: r[i] = a[i]^0.5

      INTERFACE
        SUBROUTINE vcsqrt(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzsqrt(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcsqrt(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzsqrt(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Reciprocal Square Root: r[i] = 1/a[i]^0.5

      INTERFACE
        SUBROUTINE vsinvsqrt(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdinvsqrt(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsinvsqrt(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdinvsqrt(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Cube Root: r[i] = a[i]^(1/3)

      INTERFACE
        SUBROUTINE vscbrt(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdcbrt(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmscbrt(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdcbrt(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Reciprocal Cube Root: r[i] = 1/a[i]^(1/3)

      INTERFACE
        SUBROUTINE vsinvcbrt(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdinvcbrt(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsinvcbrt(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdinvcbrt(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Squaring: r[i] = a[i]^2

      INTERFACE
        SUBROUTINE vssqr(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdsqr(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmssqr(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdsqr(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Exponential Function: r[i] = e^a[i]

      INTERFACE
        SUBROUTINE vsexp(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdexp(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsexp(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdexp(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex Exponential Function: r[i] = e^a[i]

      INTERFACE
        SUBROUTINE vcexp(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzexp(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcexp(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzexp(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Exponential of arguments decreased by 1: r[i] = e^(a[i]-1)

      INTERFACE
        SUBROUTINE vsexpm1(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdexpm1(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsexpm1(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdexpm1(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Logarithm (base e): r[i] = ln(a[i])

      INTERFACE
        SUBROUTINE vsln(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdln(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsln(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdln(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex Logarithm (base e): r[i] = ln(a[i])

      INTERFACE
        SUBROUTINE vcln(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzln(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcln(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzln(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Logarithm (base 10): r[i] = lg(a[i])

      INTERFACE
        SUBROUTINE vslog10(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdlog10(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmslog10(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdlog10(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Logarithm (base e) of arguments increased by 1: r[i] = log(1+a[i])

      INTERFACE
        SUBROUTINE vslog1p(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdlog1p(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmslog1p(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdlog1p(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex Logarithm (base 10): r[i] = lg(a[i])

      INTERFACE
        SUBROUTINE vclog10(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzlog10(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmclog10(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzlog10(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Sine: r[i] = sin(a[i])

      INTERFACE
        SUBROUTINE vssin(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdsin(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmssin(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdsin(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex Sine: r[i] = sin(a[i])

      INTERFACE
        SUBROUTINE vcsin(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzsin(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcsin(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzsin(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Cosine: r[i] = cos(a[i])

      INTERFACE
        SUBROUTINE vscos(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdcos(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmscos(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdcos(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex Cosine: r[i] = cos(a[i])

      INTERFACE
        SUBROUTINE vccos(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzcos(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmccos(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzcos(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Tangent: r[i] = tan(a[i])

      INTERFACE
        SUBROUTINE vstan(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdtan(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmstan(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdtan(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex Tangent: r[i] = tan(a[i])

      INTERFACE
        SUBROUTINE vctan(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vztan(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmctan(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmztan(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Hyperbolic Sine: r[i] = sh(a[i])

      INTERFACE
        SUBROUTINE vssinh(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdsinh(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmssinh(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdsinh(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex Hyperbolic Sine: r[i] = sh(a[i])

      INTERFACE
        SUBROUTINE vcsinh(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzsinh(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcsinh(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzsinh(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Hyperbolic Cosine: r[i] = ch(a[i])

      INTERFACE
        SUBROUTINE vscosh(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdcosh(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmscosh(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdcosh(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex Hyperbolic Cosine: r[i] = ch(a[i])

      INTERFACE
        SUBROUTINE vccosh(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzcosh(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmccosh(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzcosh(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Hyperbolic Tangent: r[i] = th(a[i])

      INTERFACE
        SUBROUTINE vstanh(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdtanh(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmstanh(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdtanh(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex Hyperbolic Tangent: r[i] = th(a[i])

      INTERFACE
        SUBROUTINE vctanh(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vztanh(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmctanh(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmztanh(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Arc Cosine: r[i] = arccos(a[i])

      INTERFACE
        SUBROUTINE vsacos(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdacos(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsacos(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdacos(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex Arc Cosine: r[i] = arccos(a[i])

      INTERFACE
        SUBROUTINE vcacos(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzacos(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcacos(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzacos(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Arc Sine: r[i] = arcsin(a[i])

      INTERFACE
        SUBROUTINE vsasin(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdasin(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsasin(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdasin(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex Arc Sine: r[i] = arcsin(a[i])

      INTERFACE
        SUBROUTINE vcasin(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzasin(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcasin(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzasin(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Arc Tangent: r[i] = arctan(a[i])

      INTERFACE
        SUBROUTINE vsatan(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdatan(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsatan(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdatan(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex Arc Tangent: r[i] = arctan(a[i])

      INTERFACE
        SUBROUTINE vcatan(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzatan(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcatan(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzatan(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Hyperbolic Arc Cosine: r[i] = arcch(a[i])

      INTERFACE
        SUBROUTINE vsacosh(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdacosh(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsacosh(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdacosh(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex Hyperbolic Arc Cosine: r[i] = arcch(a[i])

      INTERFACE
        SUBROUTINE vcacosh(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzacosh(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcacosh(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzacosh(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Hyperbolic Arc Sine: r[i] = arcsh(a[i])

      INTERFACE
        SUBROUTINE vsasinh(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdasinh(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsasinh(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdasinh(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex Hyperbolic Arc Sine: r[i] = arcsh(a[i])

      INTERFACE
        SUBROUTINE vcasinh(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzasinh(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcasinh(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzasinh(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Hyperbolic Arc Tangent: r[i] = arcth(a[i])

      INTERFACE
        SUBROUTINE vsatanh(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdatanh(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsatanh(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdatanh(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex Hyperbolic Arc Tangent: r[i] = arcth(a[i])

      INTERFACE
        SUBROUTINE vcatanh(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzatanh(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcatanh(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzatanh(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Error Function: r[i] = erf(a[i])

      INTERFACE
        SUBROUTINE vserf(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vderf(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmserf(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmderf(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Inverse error function: r[i] = erfinv(a[i])

      INTERFACE
        SUBROUTINE vserfinv(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vderfinv(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmserfinv(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmderfinv(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Square root of the sum of the squares: r[i] = hypot(a[i],b[i])

      INTERFACE
        SUBROUTINE vshypot(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(IN)    :: b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdhypot(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(IN)    :: b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmshypot(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(IN)    :: b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdhypot(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(IN)    :: b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complementary Error Function: r[i] = 1 - erf(a[i])

      INTERFACE
        SUBROUTINE vserfc(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vderfc(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmserfc(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmderfc(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Inverse complementary error function: r[i] = erfcinv(a[i])

      INTERFACE
        SUBROUTINE vserfcinv(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vderfcinv(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmserfcinv(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmderfcinv(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Cumulative normal distribution function: r[i] = cdfnorm(a[i])

      INTERFACE
        SUBROUTINE vscdfnorm(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdcdfnorm(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmscdfnorm(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdcdfnorm(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Inverse cumulative normal distribution function:
!  r[i] = cdfnorminv(a[i])

      INTERFACE
        SUBROUTINE vscdfnorminv(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdcdfnorminv(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmscdfnorminv(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdcdfnorminv(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Logarithm (base e) of the absolute value of gamma function:
!  r[i] = lgamma(a[i])

      INTERFACE
        SUBROUTINE vslgamma(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdlgamma(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmslgamma(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdlgamma(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Gamma function: r[i] = tgamma(a[i])

      INTERFACE
        SUBROUTINE vstgamma(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdtgamma(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmstgamma(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdtgamma(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Arc Tangent of a/b: r[i] = arctan(a[i]/b[i])

      INTERFACE
        SUBROUTINE vsatan2(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdatan2(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsatan2(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdatan2(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Multiplicaton: r[i] = a[i] * b[i]

      INTERFACE
        SUBROUTINE vsmul(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdmul(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsmul(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdmul(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex multiplicaton: r[i] = a[i] * b[i]

      INTERFACE
        SUBROUTINE vcmul(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          COMPLEX(KIND=4),INTENT(IN)    :: a(n),b(n)
          COMPLEX(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzmul(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          COMPLEX(KIND=8),INTENT(IN)    :: a(n),b(n)
          COMPLEX(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcmul(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          COMPLEX(KIND=4),INTENT(IN)    :: a(n),b(n)
          COMPLEX(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN)    :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzmul(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          COMPLEX(KIND=8),INTENT(IN)    :: a(n),b(n)
          COMPLEX(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN)    :: mode
        END SUBROUTINE
      END INTERFACE

!  Division: r[i] = a[i] / b[i]

      INTERFACE
        SUBROUTINE vsdiv(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vddiv(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsdiv(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmddiv(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex division: r[i] = a[i] / b[i]

      INTERFACE
        SUBROUTINE vcdiv(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          COMPLEX(KIND=4),INTENT(IN)    :: a(n),b(n)
          COMPLEX(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzdiv(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          COMPLEX(KIND=8),INTENT(IN)    :: a(n),b(n)
          COMPLEX(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcdiv(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          COMPLEX(KIND=4),INTENT(IN)    :: a(n),b(n)
          COMPLEX(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN)    :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzdiv(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          COMPLEX(KIND=8),INTENT(IN)    :: a(n),b(n)
          COMPLEX(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN)    :: mode
        END SUBROUTINE
      END INTERFACE

!  Power Function: r[i] = a[i]^b[i]

      INTERFACE
        SUBROUTINE vspow(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdpow(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmspow(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdpow(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex Power Function: r[i] = a[i]^b[i]

      INTERFACE
        SUBROUTINE vcpow(n,a,b,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n),b(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzpow(n,a,b,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n),b(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcpow(n,a,b,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n),b(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzpow(n,a,b,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n),b(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Power Function: r[i] = a[i]^(3/2)

      INTERFACE
        SUBROUTINE vspow3o2(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdpow3o2(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmspow3o2(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdpow3o2(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Power Function: r[i] = a[i]^(2/3)

      INTERFACE
        SUBROUTINE vspow2o3(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdpow2o3(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmspow2o3(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdpow2o3(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Power Function with Fixed Degree: r[i] = a[i]^b

      INTERFACE
        SUBROUTINE vspowx(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: b
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdpowx(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: b
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmspowx(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: b
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdpowx(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: b
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex Power Function with Fixed Degree: r[i] = a[i]^b

      INTERFACE
        SUBROUTINE vcpowx(n,a,b,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: b
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzpowx(n,a,b,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: b
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcpowx(n,a,b,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: b
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzpowx(n,a,b,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: b
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Sine & Cosine: r1[i] = sin(a[i]), r2[i]=cos(a[i])

      INTERFACE
        SUBROUTINE vssincos(n,a,r1,r2)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r1(n),r2(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdsincos(n,a,r1,r2)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r1(n),r2(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmssincos(n,a,r1,r2,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r1(n),r2(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdsincos(n,a,r1,r2,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r1(n),r2(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Linear fraction: r[i] = (a[i]*scalea + shifta)/(b[i]*scaleb + shiftb)

      INTERFACE
        SUBROUTINE vslinearfrac(n,a,b,scalea,shifta,scaleb,shiftb,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(IN)    :: scalea,shifta,scaleb,shiftb
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdlinearfrac(n,a,b,scalea,shifta,scaleb,shiftb,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(IN)    :: scalea,shifta,scaleb,shiftb
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmslinearfrac(n,a,b,scalea,shifta,scaleb,shiftb,r,   &
     &                           mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(IN)    :: scalea,shifta,scaleb,shiftb
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdlinearfrac(n,a,b,scalea,shifta,scaleb,shiftb,r,   &
     &                           mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(IN)    :: scalea,shifta,scaleb,shiftb
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Integer value rounded towards plus infinity: r[i] = ceil(a[i])

      INTERFACE
        SUBROUTINE vsceil(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdceil(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsceil(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdceil(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Integer value rounded towards minus infinity: r[i] = floor(a[i])

      INTERFACE
        SUBROUTINE vsfloor(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdfloor(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsfloor(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdfloor(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Signed fraction part

      INTERFACE
        SUBROUTINE vsfrac(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdfrac(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsfrac(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdfrac(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE


!  Truncated integer value and the remaining fraction part

      INTERFACE
        SUBROUTINE vsmodf(n,a,r1,r2)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r1(n),r2(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdmodf(n,a,r1,r2)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r1(n),r2(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsmodf(n,a,r1,r2,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r1(n),r2(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdmodf(n,a,r1,r2,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r1(n),r2(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Rounded integer value in the current rounding mode:
!  r[i] = nearbyint(a[i])

      INTERFACE
        SUBROUTINE vsnearbyint(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdnearbyint(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsnearbyint(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdnearbyint(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Rounded integer value in the current rounding mode with inexact
!  result exception raised for rach changed value: r[i] = rint(a[i])

      INTERFACE
        SUBROUTINE vsrint(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdrint(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsrint(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdrint(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Value rounded to the nearest integer: r[i] = round(a[i])

      INTERFACE
        SUBROUTINE vsround(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdround(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsround(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdround(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Truncated integer value and the temaining fraction part:
!  r[i] = trunc(a[i])

      INTERFACE
        SUBROUTINE vstrunc(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdtrunc(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmstrunc(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdtrunc(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Element by element conjugation: r[i] = conj(a[i])

      INTERFACE
        SUBROUTINE vcconj(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzconj(n,a,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcconj(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzconj(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Element by element multiplication of vector A element and
!  conjugated vector B element: r[i] = mulbyconj(a[i],b[i])

      INTERFACE
        SUBROUTINE vcmulbyconj(n,a,b,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n),b(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzmulbyconj(n,a,b,r)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n),b(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmcmulbyconj(n,a,b,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=4),INTENT(IN)  :: a(n),b(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzmulbyconj(n,a,b,r,mode)
          INTEGER,INTENT(IN)  :: n
          COMPLEX(KIND=8),INTENT(IN)  :: a(n),b(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Complex exponent of real vector elements: r[i] = CIS(a[i])

      INTERFACE
        SUBROUTINE vccis(n,a,r)
          INTEGER,INTENT(IN)  :: n
          REAL(KIND=4),INTENT(IN)     :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzcis(n,a,r)
          INTEGER,INTENT(IN)  :: n
          REAL(KIND=8),INTENT(IN)     :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmccis(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          REAL(KIND=4),INTENT(IN)     :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmzcis(n,a,r,mode)
          INTEGER,INTENT(IN)  :: n
          REAL(KIND=8),INTENT(IN)     :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: r(n)
          INTEGER(KIND=8),INTENT(IN)  :: mode
        END SUBROUTINE
      END INTERFACE

!  Exponential integral: r[i] = e1(a[i])

      INTERFACE
        SUBROUTINE vsexpint1(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdexpint1(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsexpint1(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdexpint1(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Cosine: r[i] = cos(a[i]*PI)

      INTERFACE
        SUBROUTINE vscospi(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdcospi(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmscospi(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdcospi(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Sine: r[i] = sin(a[i]*PI)

      INTERFACE
        SUBROUTINE vssinpi(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdsinpi(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmssinpi(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdsinpi(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Tangent: r[i] = tan(a[i]*PI)

      INTERFACE
        SUBROUTINE vstanpi(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdtanpi(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmstanpi(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdtanpi(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Arc cosine: r[i] = acos(a[i])/PI

      INTERFACE
        SUBROUTINE vsacospi(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdacospi(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsacospi(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdacospi(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Arc sine: r[i] = asin(a[i])/PI

      INTERFACE
        SUBROUTINE vsasinpi(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdasinpi(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsasinpi(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdasinpi(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Arc tangent: r[i] = atan(a[i])/PI

      INTERFACE
        SUBROUTINE vsatanpi(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdatanpi(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsatanpi(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdatanpi(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Cosine: r[i] = cos(a[i]*PI/180)

      INTERFACE
        SUBROUTINE vscosd(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdcosd(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmscosd(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdcosd(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Sine: r[i] = sin(a[i]*PI/180)

      INTERFACE
        SUBROUTINE vssind(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdsind(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmssind(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdsind(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Tangent: r[i] = tan(a[i]*PI/180)

      INTERFACE
        SUBROUTINE vstand(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdtand(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmstand(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdtand(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Exponential function (base 2): r[i] = 2^a[i]

      INTERFACE
        SUBROUTINE vsexp2(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdexp2(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsexp2(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdexp2(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Exponential function (base 10): r[i] = 10^a[i]

      INTERFACE
        SUBROUTINE vsexp10(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdexp10(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsexp10(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdexp10(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Logarithm (base 2): r[i] = lb(a[i])

      INTERFACE
        SUBROUTINE vslog2(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdlog2(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmslog2(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdlog2(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Computes the exponent: r[i] = logb(a[i])

      INTERFACE
        SUBROUTINE vslogb(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdlogb(n,a,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmslogb(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdlogb(n,a,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

! Arc Tangent of a/b Divided by PI: r[i] = arctan(a[i]/b[i])/PI

      INTERFACE
        SUBROUTINE vsatan2pi(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdatan2pi(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsatan2pi(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdatan2pi(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Power Function with a[i]>=0: r[i] = a[i]^b[i]

      INTERFACE
        SUBROUTINE vspowr(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdpowr(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmspowr(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdpowr(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Remainder Function: r[i] = remainder(a[i], b[i])

      INTERFACE
        SUBROUTINE vsremainder(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdremainder(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsremainder(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdremainder(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Modulus Function: r[i] = fmod(a[i], b[i])

      INTERFACE
        SUBROUTINE vsfmod(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdfmod(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsfmod(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdfmod(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Copy Sign Function: r[i] = copysign(a[i], b[i])

      INTERFACE
        SUBROUTINE vscopysign(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdcopysign(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmscopysign(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdcopysign(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Next After Function: r[i] = nextafter(a[i], b[i])

      INTERFACE
        SUBROUTINE vsnextafter(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdnextafter(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsnextafter(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdnextafter(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Positive Difference Function: r[i] = fdim(a[i], b[i])

      INTERFACE
        SUBROUTINE vsfdim(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdfdim(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsfdim(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdfdim(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Maximum Function: r[i] = fmax(a[i], b[i])

      INTERFACE
        SUBROUTINE vsfmax(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdfmax(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsfmax(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdfmax(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Minimum Function: r[i] = fmin(a[i], b[i])

      INTERFACE
        SUBROUTINE vsfmin(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdfmin(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsfmin(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdfmin(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Maximum Magnitude Function: r[i] = maxmag(a[i], b[i])

      INTERFACE
        SUBROUTINE vsmaxmag(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdmaxmag(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsmaxmag(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdmaxmag(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!  Minimum magnitude function: r[i] = minmag(a[i], b[i])

      INTERFACE
        SUBROUTINE vsminmag(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdminmag(n,a,b,r)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmsminmag(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=4),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=4),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vmdminmag(n,a,b,r,mode)
          INTEGER,INTENT(IN) :: n
          REAL(KIND=8),INTENT(IN)    :: a(n),b(n)
          REAL(KIND=8),INTENT(OUT)   :: r(n)
          INTEGER(KIND=8),INTENT(IN) :: mode
        END SUBROUTINE
      END INTERFACE

!++
!  VML PACK FUNCTION INTERFACES.
!--

!  Positive Increment Indexing
      INTERFACE
        SUBROUTINE vspacki(n,a,incra,y)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: incra
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdpacki(n,a,incra,y)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: incra
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vcpacki(n,a,incra,y)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: incra
          COMPLEX(KIND=4),INTENT(IN)    :: a(n)
          COMPLEX(KIND=4),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzpacki(n,a,incra,y)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: incra
          COMPLEX(KIND=8),INTENT(IN)    :: a(n)
          COMPLEX(KIND=8),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

!  Index Vector Indexing
      INTERFACE
        SUBROUTINE vspackv(n,a,ia,y)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: ia(n)
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdpackv(n,a,ia,y)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: ia(n)
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vcpackv(n,a,ia,y)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: ia(n)
          COMPLEX(KIND=4),INTENT(IN)    :: a(n)
          COMPLEX(KIND=4),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzpackv(n,a,ia,y)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: ia(n)
          COMPLEX(KIND=8),INTENT(IN)    :: a(n)
          COMPLEX(KIND=8),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

!  Mask Vector Indexing
      INTERFACE
        SUBROUTINE vspackm(n,a,ma,y)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: ma(n)
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdpackm(n,a,ma,y)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: ma(n)
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vcpackm(n,a,ma,y)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: ma(n)
          COMPLEX(KIND=4),INTENT(IN)    :: a(n)
          COMPLEX(KIND=4),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzpackm(n,a,ma,y)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: ma(n)
          COMPLEX(KIND=8),INTENT(IN)    :: a(n)
          COMPLEX(KIND=8),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

!++
!  VML UNPACK FUNCTION DECLARATIONS.
!--

!  Positive Increment Indexing
      INTERFACE
        SUBROUTINE vsunpacki(n,a,y,incry)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: incry
          REAL(KIND=4),INTENT(IN)  :: a(n)
          REAL(KIND=4),INTENT(OUT) :: y(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdunpacki(n,a,y,incry)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: incry
          REAL(KIND=8),INTENT(IN)  :: a(n)
          REAL(KIND=8),INTENT(OUT) :: y(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vcunpacki(n,a,y,incry)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: incry
          COMPLEX(KIND=4),INTENT(IN)  :: a(n)
          COMPLEX(KIND=4),INTENT(OUT) :: y(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzunpacki(n,a,y,incry)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: incry
          COMPLEX(KIND=8),INTENT(IN)  :: a(n)
          COMPLEX(KIND=8),INTENT(OUT) :: y(n)
        END SUBROUTINE
      END INTERFACE

!  Index Vector Indexing
      INTERFACE
        SUBROUTINE vsunpackv(n,a,y,iy)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: iy(n)
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdunpackv(n,a,y,iy)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: iy(n)
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vcunpackv(n,a,y,iy)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: iy(n)
          COMPLEX(KIND=4),INTENT(IN)    :: a(n)
          COMPLEX(KIND=4),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzunpackv(n,a,y,iy)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: iy(n)
          COMPLEX(KIND=8),INTENT(IN)    :: a(n)
          COMPLEX(KIND=8),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

!  Mask Vector Indexing
      INTERFACE
        SUBROUTINE vsunpackm(n,a,y,my)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: my(n)
          REAL(KIND=4),INTENT(IN)    :: a(n)
          REAL(KIND=4),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vdunpackm(n,a,y,my)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: my(n)
          REAL(KIND=8),INTENT(IN)    :: a(n)
          REAL(KIND=8),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vcunpackm(n,a,y,my)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: my(n)
          COMPLEX(KIND=4),INTENT(IN)    :: a(n)
          COMPLEX(KIND=4),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

      INTERFACE
        SUBROUTINE vzunpackm(n,a,y,my)
          INTEGER,INTENT(IN) :: n
          INTEGER,INTENT(IN) :: my(n)
          COMPLEX(KIND=8),INTENT(IN)    :: a(n)
          COMPLEX(KIND=8),INTENT(OUT)   :: y(n)
        END SUBROUTINE
      END INTERFACE

!++
!  VML ERROR HANDLING FUNCTION DECLARATIONS.
!--

!  Set VML Error Status
      INTERFACE
        INTEGER(KIND=4) FUNCTION vmlseterrstatus(err)
          INTEGER,INTENT(IN) :: err
        END FUNCTION
      END INTERFACE

!  Get VML Error Status
      INTERFACE
        INTEGER(KIND=4) FUNCTION vmlgeterrstatus()
        END FUNCTION
      END INTERFACE

!  Clear VML Error Status
      INTERFACE
        INTEGER(KIND=4) FUNCTION vmlclearerrstatus()
        END FUNCTION
      END INTERFACE

!  Set VML Error Callback Function
      INTERFACE
        INTEGER(KIND=4) FUNCTION vmlseterrorcallback(cb)
          INTEGER,EXTERNAL :: cb
        END FUNCTION vmlseterrorcallback
      END INTERFACE

!  Get VML Error Callback Function
      INTERFACE
        INTEGER(KIND=4) FUNCTION vmlgeterrorcallback()
        END FUNCTION vmlgeterrorcallback
      END INTERFACE

!  Reset VML Error Callback Function
      INTERFACE
        INTEGER(KIND=4) FUNCTION vmlclearerrorcallback()
        END FUNCTION vmlclearerrorcallback
      END INTERFACE

!++
!  VML MODE FUNCTION DECLARATIONS.
!--

!  Set VML Mode
      INTERFACE
        INTEGER(KIND=4) FUNCTION vmlsetmode(n)
          INTEGER,INTENT(IN) :: n
        END FUNCTION
      END INTERFACE

!  Get VML Mode
      INTERFACE
        INTEGER(KIND=4) FUNCTION vmlgetmode()
        END FUNCTION
      END INTERFACE
