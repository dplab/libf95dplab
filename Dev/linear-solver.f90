!  $Id$
!
PROGRAM LinearSolver
  !
  !  Test PETSc linear solver interface.
  !
  USE LinearSolverPETScModule

  IMPLICIT NONE
  !
  !  ========== Locals
  !
  !  Parameters
  !
  INTEGER, PARAMETER :: RK = KIND(1.0d0)
  CHARACTER(LEN=*), PARAMETER :: NAME = 'Test Solver'
  !
  !  The solver
  !
  TYPE(LinearSolverType) :: solver  
  !
  !  Array dimensions and data.
  !
  INTEGER :: nDOFpe, nElem, nDOFgl, numBC, status
  INTEGER :: myEmin, myEmax ! elements on this process
  INTEGER, ALLOCATABLE :: conn(:, :), bcNodes(:), locSizes(:)
  !
  REAL(RK), ALLOCATABLE :: eMats(:, :, :), eRhs(:, :)
  REAL(RK), ALLOCATABLE :: sol(:)
  REAL(RK), ALLOCATABLE  :: ARHS(:)
  REAL(RK), ALLOCATABLE  :: BCVALS(:)
  REAL(RK), ALLOCATABLE  :: SOL_GLOBAL(:)
  !
  ! ================================================== Executable Code
  !
  !  Read Matrices and RHS.
  !
  CALL LinearSolverPETScInit()
  !
  call ReadInput()
  !
  CALL LinearSolverCreate(solver, &
       &   NAME, conn, bcnodes, locsizes)

  CALL LinearSolverSolve(solver, sol, eMats,&
       &   ERHS=erhs, STATUS=status)
  !
  WRITE(6, *) 'The solution is:  '
  WRITE(6, *) sol

  CALL LinearSolverPETScFinalize()

CONTAINS ! ========================= Internal Procedures

  ! ================================ BEGIN:  ReadInput
  !
  SUBROUTINE ReadInput()
    !
    !  Read input.
    !
    !  Comment lines are expected before each block.
    !
    ! ========== Arguments
    !

    !
    ! ========== Locals
    !
    INTEGER :: unit = 11
    CHARACTER(LEN=128) :: line
    !
    ! ============================== Executable Code
    !
    OPEN(UNIT=unit, FILE='input.txt', ACTION='READ')
    !
    READ(unit, '(a)') line;  WRITE(6, *) myrank, ':', TRIM(line)
    READ(unit, '(a)') line;  WRITE(6, *) myrank, ':', TRIM(line)
    READ(unit, *) nDOFpe, nElem, nDOFgl, numBC

    call SetLocalSizes()

    ALLOCATE(conn(nDOFpe, nElem), bcNodes(numBC))
    ALLOCATE(eMats(nDOFpe, nDOFpe, nElem), &
         &   eRhs(nDOFpe, nElem),&
         &   sol(locsizes(myrank)))
    !
    READ(unit, '(a)') line;  WRITE(6, *) myrank, ':', TRIM(line)
    READ(unit, *) conn

    READ(unit, '(a)') line;  WRITE(6, *) myrank, ':', TRIM(line)
    READ(unit, *) bcnodes

    READ(unit, '(a)') line;  WRITE(6, *) myrank, ':', TRIM(line)
    READ(unit, *) eMats
    WRITE(6, *) eMats(:, :, myEmin:myEmax)

    READ(unit, '(a)') line;  WRITE(6, *) myrank, ':', TRIM(line)
    READ(unit, *) eRhs
    WRITE(6, *) eRhs(:, myEmin:myEmax)
    !
    !
  END SUBROUTINE ReadInput
  !
  ! =================================   END:  ReadInput
  ! ================================ BEGIN:  SetLocalSizes
  !
  SUBROUTINE SetLocalSizes()
    !
    !  Partition data among processes.
    !
    ! ========== Locals
    !
    INTEGER :: perProc, numLeft, p, EperProc
    !
    ! ============================== Executable Code
    !
    ALLOCATE(locsizes(0:numProcs-1))
    perProc = nDOFgl / numProcs
    numLeft = nDOFgl - perProc*numProcs
    locsizes = perProc
    locsizes(1:numLeft) = perProc + 1
    PRINT *, 'local sizes:  ', locsizes
    !
    !  Set element ranges.
    !
    EperProc = nElem / numProcs
    numLeft  = nElem - EperProc*numProcs
    myEmin = 1
    DO p=0, myrank
      IF (myRank == p) THEN
        myEmax = myEmin + EperProc - 1
        IF (p < numLeft) myEmax = myEmax + 1
      ELSE
        myEmin = myEmin + EperProc
        IF (p < numLeft) myEmin = myEmin + 1
      END IF
    END DO
    PRINT *, 'my rank and elements:  ', myrank, ':  ', myEmin, myEmax
    !
  END SUBROUTINE SetLocalSizes
  !
  ! =================================   END:  SetLocalSizes
  !
END PROGRAM LinearSolver
