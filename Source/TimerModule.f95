MODULE TimerModule 
  !
  !  This module handles the Timer object.
  !
  !  * Uses CPU_TIME instead of SYSTEM_CLOCK
  !
  !  USAGE:
  !
  !  SUBROUTINE TimerInit(self, name, ACTIVATE)
  !  SUBROUTINE TimerStart(timer)
  !  SUBROUTINE TimerStop(timer)
  !  TimerWrite(timer, unit)
  !
  !--------------Other Modules
  !
  IMPLICIT NONE
  !
  !--------------Types
  !
  !  Basic type for timer.
  !
  TYPE TimerType
    CHARACTER(LEN=32) :: name = ''
    LOGICAL :: active = .TRUE., failed = .FALSE.
    INTEGER :: starts = 0, stops = 0
    REAL :: time_beg = 0.0, time_end = 0.0, time_total = 0.0
  END TYPE TimerType
  !
  !-------------- Public Entities
  !
  PRIVATE  ! all objects are private unless declared otherwise
  !
  PUBLIC :: TimerType
  !
  PUBLIC TimerInit, TimerActivate, TimerStart, TimerStop, TimerWrite
  !
CONTAINS ! ===============================================================
  !
  SUBROUTINE TimerInit(self, name, ACTIVATE)
    !
    !  Initialization for the Timer object.
    !
    TYPE(TimerType), INTENT(IN OUT) :: self
    CHARACTER(LEN=*), INTENT(IN)    :: name
    LOGICAL, INTENT(IN), OPTIONAL   :: ACTIVATE
    !
    !  ACTIVATE -- optional argument to activate the timer immediately
    !
    !--------------*-------------------------------------------------------
    !
    self % name   = name
    !
    if (PRESENT(ACTIVATE)) then
      self % active = ACTIVATE
    end if
    !
  END SUBROUTINE TimerInit
  !
  !
  SUBROUTINE TimerActivate(timer, ACTIVATE)
    !
    !  Activate timer.
    !
    !--------------*-------------------------------------------------------
    !
    !  Arguments.
    !
    TYPE(TimerType), INTENT(IN OUT) :: timer
    LOGICAL, INTENT(IN), OPTIONAL :: ACTIVATE
    !
    !--------------*-------------------------------------------------------
    !
    IF (PRESENT(ACTIVATE)) THEN
      timer % active = ACTIVATE
    ELSE
      timer % active = .TRUE.
    END IF
    !
  END SUBROUTINE TimerActivate
  !
  !
  SUBROUTINE TimerStart(timer)
    !
    !  Start the timer.  
    !
    TYPE(TimerType), INTENT(IN OUT) :: timer
    !
    !  timer -- the timer to start
    !
    !--------------*-------------------------------------------------------
    !
    if (timer % active) then
      CALL CPU_TIME(timer % time_beg)
      IF (timer % time_beg < 0.0) THEN
        timer % active = .FALSE.
        timer % failed = .TRUE.
      ELSE
        timer % starts = timer % starts + 1
      END IF
    end if
    !
  END SUBROUTINE TimerStart
  !
  !
  SUBROUTINE TimerStop(timer)
    !
    !  Stop timer and compute elapsed time.
    !
    TYPE(TimerType), INTENT(IN OUT) :: timer
    !
    !  timer -- the timer to stop
    !
    !--------------*-------------------------------------------------------
    !
    if (timer % active) then
      CALL CPU_TIME(timer % time_end)
      timer % time_total = timer % time_total + timer % time_end - timer % time_beg
      timer % stops = timer % stops + 1
    end if
    !
  END SUBROUTINE TimerStop
  !
  !
  SUBROUTINE TimerWrite(timer, unit, level)
    !
    !  Print statistics for each timer.
    !
    TYPE(TimerType), INTENT(IN) :: timer
    INTEGER, INTENT(IN) :: unit
    INTEGER, INTENT(IN), OPTIONAL :: level
    !
    !  Locals.
    !
    CHARACTER(LEN=4) :: INDENT='    '
    INTEGER :: mylevel
    CHARACTER(LEN=32) :: myfmt
    !
    !--------------*-------------------------------------------------------
    !
    write(unit, '(a)') 'Statistics for timer:  '//timer % name
    if (timer % active) then
      WRITE(unit, '(a)')  INDENT // 'timer is currently active'
    else
      WRITE(unit, '(a)')  INDENT // 'timer is currently inactive'
    end if
    !
    if (timer % failed) then
      WRITE(unit, '(a)')  INDENT // '*** timer failed'
    end if
    !
    WRITE(unit, '(a,i0)') INDENT // 'timer starts:  ', timer % starts
    WRITE(unit, '(a,i0)') INDENT // 'timer stops:  ', timer % stops
    WRITE(unit, '(a,f0.4)') INDENT // 'total time in seconds:  ', timer % time_total
    !
    if (timer % stops > 0) then
      WRITE(unit, '(a,f0.4)') INDENT // 'average time per call:  ', timer % time_total/REAL(timer % stops)
    ELSE
      WRITE(unit, '(a)') INDENT // 'timer has not been used'
    end if
    !
    !  Summary line.
    !
    IF (PRESENT(level)) THEN
      mylevel = 2*level
    ELSE
      mylevel = 0
    END IF
    !
    IF (mylevel > 0) THEN
      WRITE(myfmt, '(a,i0,a)') '(a,a35,', mylevel, 'x,f0.3,a,i0,a,f0.3)'
    ELSE
      myfmt = '(a,a35,f0.3,a,i0,a,f0.3)'
    END IF

    WRITE(unit, myfmt) 'timer summary for ', ADJUSTR(timer % name) // ':  ', &
         &   timer % time_total, '/', timer % stops, '=',timer % time_total/REAL(timer % stops)
    !
  END SUBROUTINE TimerWrite
  !
END MODULE TimerModule
