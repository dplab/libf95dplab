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
    OPEN(UNIT=unit, FILE='input.txt')
    !
    READ(unit, '(a)') line;  WRITE(6, *) TRIM(line)
    READ(unit, '(a)') line;  WRITE(6, *) TRIM(line)
    READ(unit, *) nDOFpe, nElem, nDOFgl, numBC

    call SetLocalSizes()

    ALLOCATE(conn(nDOFpe, nElem), bcNodes(numBC))
    ALLOCATE(eMats(nDOFpe, nDOFpe, nElem), &
         &   eRhs(nDOFpe, nElem),&
         &   sol(locsizes(myrank)))
    !
    READ(unit, '(a)') line;  WRITE(6, *) TRIM(line)
    READ(unit, *) conn

    READ(unit, '(a)') line;  WRITE(6, *) TRIM(line)
    READ(unit, *) bcnodes

    READ(unit, '(a)') line;  WRITE(6, *) TRIM(line)
    READ(unit, *) eMats
    WRITE(6, *) eMats

    READ(unit, '(a)') line;  WRITE(6, *) TRIM(line)
    READ(unit, *) eRhs
    WRITE(6, *) eRhs
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
    INTEGER :: perProc, numLeft
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
  END SUBROUTINE SetLocalSizes
  !
  ! =================================   END:  SetLocalSizes
  !
END PROGRAM LinearSolver
