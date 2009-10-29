!  $Id$
!
MODULE CommandModule
  !
  !  Tools for processing command strings.
  !
  !--------------Other Modules-------------------------------------------
  !
  USE StringsModule, ONLY: getWord

  IMPLICIT NONE
  !
  !--------------Data----------------------------------------------------
  !
  !
  !--------------Types---------------------------------------------------
  !

  !
  !--------------Public Data---------------------------------------------
  !
  PRIVATE   ! all objects are private unless declared otherwise
  !
  PUBLIC :: ExecCommand
  !
CONTAINS
  !
  FUNCTION ExecCommand(command, cmdLine, callback, status) RESULT(Match)
    !
    !  Call the callback, if command name matches.
    !  Return true on a match and false otherwise.
    !
    ! ========== Arguments
    !
    ! cmdLine  -- string with command and arguments
    ! callback -- subroutine to call if command matches
    !
    CHARACTER(LEN=*), INTENT(IN) :: command
    CHARACTER(LEN=*), INTENT(IN) :: cmdLine
    !
    INTERFACE
      SUBROUTINE callback(a, s)
        CHARACTER(LEN=*), INTENT(IN)  :: a
        INTEGER,          INTENT(OUT) :: s
      END SUBROUTINE callback
    END INTERFACE
    !
    INTEGER, INTENT(OUT) :: status
    !  
    !  Result.
    !  
    LOGICAL :: Match
    !
    ! ========== Locals
    !
    CHARACTER(LEN=LEN(cmdLine)) :: word, args
    INTEGER :: iStat
    !
    !--------------*-------------------------------------------------------
    !
    status = 0
    args = cmdLine
    !
    CALL getWord(args, word, iStat)
    IF (word == command) THEN
      Match = .TRUE.
      CALL callback(args, status)
    ELSE
      Match = .FALSE.
    END IF
    !
  END FUNCTION ExecCommand
  !
END MODULE CommandModule
