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
  call ReadInput()
  !
  CALL LinearSolverCreate(solver, &
       &   NAME, conn, bcnodes, locsizes)

  CALL LinearSolverSolve(solver, sol, eMats,&
       &   ERHS=erhs, STATUS=status)
  !
  WRITE(6, *) 'The solution is:  '
  WRITE(6, *) sol

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
    CHARACTER(LEN=128) :: line
    !
    ! ============================== Executable Code
    !
    READ(5, '(a)') line;  WRITE(6, *) TRIM(line)
    READ(5, '(a)') line;  WRITE(6, *) TRIM(line)
    READ(5, *) nDOFpe, nElem, nDOFgl, numBC
    ALLOCATE(locsizes(1))
    locsizes(1) = nDOFgl
    ALLOCATE(conn(nDOFpe, nElem), bcNodes(numBC))
    ALLOCATE(eMats(nDOFpe, nDOFpe, nElem), &
         &   eRhs(nDOFpe, nElem),&
         &   sol(locsizes(1)))
    !
    READ(5, '(a)') line;  WRITE(6, *) TRIM(line)
    READ(5, *) conn

    READ(5, '(a)') line;  WRITE(6, *) TRIM(line)
    READ(5, *) bcnodes

    READ(5, '(a)') line;  WRITE(6, *) TRIM(line)
    READ(5, *) eMats
    WRITE(6, *) eMats

    READ(5, '(a)') line;  WRITE(6, *) TRIM(line)
    READ(5, *) eRhs
    WRITE(6, *) eRhs
    !
    !
  END SUBROUTINE ReadInput
  !
  ! =================================   END:  ReadInput
  !
END PROGRAM LinearSolver
