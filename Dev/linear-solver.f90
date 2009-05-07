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
  INTEGER :: nDOFpe, nElem, nDOFgl, numBC
  INTEGER. ALLOCATABLE :: conn(:, :), bcNodes(:), locSize
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
       &   NAME, conn, bcnodes, locsize)

  CALL LinearSolverSolver(solver)
  SUBROUTINE LinearSolverSolve(solver, sol, eMats,&
       &   ERHS=erhs, STATUS=status)
    !

  CONTAINS ! ========================= Internal Procedures

    ! ================================ BEGIN:  ReadInput
    !
    SUBROUTINE ReadInput()
      !
      !  Read input.
      !
      ! ========== Arguments
      !

      !
      ! ========== Locals
      !

      !
      ! ============================== Executable Code
      !
      READ(5, *) nDOFpe, nElem, nDOFgl, numBC
      locsize = nDOFgl
      ALLOCATE(conn(nDOFpe, nElem), bcNodes(numBC))
      ALLOCATE(eMats(nDOFpe, nDOFpe, nElem), &
           &   eRhs(nDOFpe, nElem),&
           &   sol(locsize))
      READ(5, *) eMats
      READ(5, *) eRhs
      !
      !
    END SUBROUTINE ReadInput
    !
    ! =================================   END:  ReadInput
  !
END PROGRAM LinearSolver
