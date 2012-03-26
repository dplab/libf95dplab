!
MODULE OptimModule
  !
  !  Optimization routines.
  !
  !  ==================== Other Modules (USE statements)
  !
  USE IntrinsicTypesModule, RK => REAL_KIND
  USE ConstantsModule, GCONJ => RK_GOLD_CONJ

  IMPLICIT NONE
  !
  !  ==================== Public Entities
  !
  !  variables, procedures, constants, derived types and namelist groups
  !
  PRIVATE   ! all objects are private unless declared otherwise
  !
CONTAINS ! ============================================= MODULE PROCEDURES
  !
  ! ==================================================== BEGIN:  minSearchGolden
  !
  SUBROUTINE minSearchGolden(f, interval, fmin, xmin, tolx)
    !
    ! Use golden section search to find a minimum.
    !
    ! In the future, more sophisticated convergence criteria could be specified.
    ! The current implementation only accepts an absolute tolerance on x, but it
    ! would be nice to have relative tolerances both on x and on the value.  For 
    ! that reason, the tolerance arguemnt is made optional so that other tolerance
    ! settings can be easily incorporated.
    !
    ! * ARGS *
    !
    ! f is the function to be minimizedi
    ! interval is a 2-vector
    !          giving endpoints of interval containing min; this is updated
    !          in this routine
    ! fmin is a scalar, the minimum value found
    ! xmin is a scalar, the location of the minimizing point
    ! tolx is a scalar, the absolute tolerance on x
    !
    INTERFACE
      FUNCTION f(x) RESULT(res)
      USE IntrinsicTypesModule, RK => REAL_KIND
        REAL(RK), INTENT(IN)  :: x
        REAL(RK) :: res
      END FUNCTION f
    END INTERFACE
    !
    REAL(RK), INTENT(IN OUT) :: interval(2)
    REAL(RK), INTENT(OUT) :: fmin, xmin
    REAL(RK), INTENT(IN) :: tolx
    !
    OPTIONAL :: tolx
    !
    !  ========== Locals
    !
    REAL(RK) :: lk, rk, ak, bk, dk, fa, fb
    !
    !  ================================================== Executable Code
    !
    lk = interval(1)
    rk = interval(2)
    dk = rk - lk
    ak = rk - GCONJ * dk
    bk = lk + GCONJ * dk
    fa = f(ak)
    fb = f(bk)
    fmin = MIN(fa, fb)
    !
    ! Loop until tolerance is reached.
    !
    ! * Interval is reduced by a constant factor each time.
    !
    DO WHILE (dk < tolx)
      IF (fa <= fb) THEN
        ! lk = lk
        rk = bk
        dk = rk - lk
        ak = rk - GCONJ * dk
        bk = ak
        fa = f(ak)
        fb = fa
      ELSE
        lk = ak
        ! rk = rk
        dk = rk - lk
        ak = bk
        bk = lk + GCONJ * dk
        fa = fb
        fb = f(bk)
      END IF
    END DO
    !
    ! Update output variables on return.
    !
    fmin = MIN(fa, fb)
    IF (fa <= fb) THEN
      interval(1) = lk
      xmin = ak
      interval(2) = bk
    ELSE
      interval(1) = ak
      xmin = bk
      interval(2) = rk
    END IF
    !
  END SUBROUTINE minSearchGolden
  !
  ! ====================================================   END:  minSearchGolden

  !
END MODULE OptimModule
