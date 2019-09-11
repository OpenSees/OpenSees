
!  --------------------------------------------------
!  THIS SUBROUTINE IS CALLED FROM C++
!  --------------------------------------------------

SUBROUTINE SDM3D (STRESS_CURRENT, &
                  STRAIN_CURRENT, &
                  STRAIN_NEXT,    &
                  MODEL_PARAMETER,&
                  SSL_VOID_RATIO, &
                  SSL_PRESSURE,   &
                  HSL_VOID_RATIO, &
                  HSL_PRESSURE,   &
                  HARDENING_PARAMETER_REAL, &
                  HARDENING_PARAMETER_INT, &
                  ANISOTROPY, &
                  TANGENT)
                  
    USE ALL_INTERFACES_3D
    IMPLICIT NONE

    REAL(KIND=DBL),DIMENSION(6),INTENT(INOUT)::STRESS_CURRENT
    REAL(KIND=DBL),DIMENSION(6),INTENT(IN)::STRAIN_NEXT
    REAL(KIND=DBL),DIMENSION(6),INTENT(INOUT)::STRAIN_CURRENT
    REAL(KIND=DBL),DIMENSION(16),INTENT(INOUT)::MODEL_PARAMETER
    REAL(KIND=DBL),DIMENSION(10),INTENT(IN)::SSL_VOID_RATIO
    REAL(KIND=DBL),DIMENSION(10),INTENT(IN)::SSL_PRESSURE
    REAL(KIND=DBL),DIMENSION(10),INTENT(IN)::HSL_VOID_RATIO
    REAL(KIND=DBL),DIMENSION(10),INTENT(IN)::HSL_PRESSURE    
    REAL(KIND=DBL),DIMENSION(13*NSURFACE+5),INTENT(INOUT)::HARDENING_PARAMETER_REAL
    INTEGER,DIMENSION(2),INTENT(INOUT)::HARDENING_PARAMETER_INT
    REAL(KIND=DBL),DIMENSION(6),INTENT(IN)::ANISOTROPY
    REAL(KIND=DBL),DIMENSION(6,6),INTENT(INOUT)::TANGENT
           
    CALL MODEL_3D(STRESS_CURRENT, &
                  STRAIN_CURRENT, &
                  STRAIN_NEXT,    &
                  MODEL_PARAMETER,&
                  SSL_VOID_RATIO, &
                  SSL_PRESSURE,   &
                  HSL_VOID_RATIO, &
                  HSL_PRESSURE,   &
                  HARDENING_PARAMETER_REAL, &
                  HARDENING_PARAMETER_INT, &
                  ANISOTROPY,   &
                  TANGENT)
                 
END SUBROUTINE SDM3D                 
