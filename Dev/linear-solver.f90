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
  !  Options
  !
  LOGICAL :: homoBC, useARHS
  !
  !  The solver
  !
  TYPE(LinearSolverType) :: solver  
  !
  !  Array dimensions and data.
  !
  INTEGER :: nDOFpe, nElem, nDOFgl, numBC, status
  INTEGER :: myEmin, myEmax ! elements on this process
  INTEGER :: myDmin, myDmax
  !
  INTEGER, ALLOCATABLE :: conn(:, :), bcNodes(:), locSizes(:)
  !
  REAL(RK), ALLOCATABLE :: eMats(:, :, :), eRhs(:, :)
  REAL(RK), ALLOCATABLE :: sol(:)
  REAL(RK), ALLOCATABLE :: arhs(:)
  REAL(RK), ALLOCATABLE :: bcVals(:)
  !
  REAL(RK), ALLOCATABLE :: solGlobal(:)
  !
  ! ================================================== Executable Code
  !
  !  Read Matrices and RHS.
  !
  CALL LinearSolverPETScInit()
  !
  call ReadInput()
  !
  CALL LinearSolverCreate(solver, NAME,&
       &   conn(:, myEmin:myEmax), bcnodes, locsizes)
  !
  !  Call solver according to options.
  !
  IF (homoBC) THEN
    IF (useARHS) THEN
      CALL LinearSolverSolve(solver, sol, eMats(:, :, myEmin:myEmax), &
           &   ERHS=erhs(:, myEmin:myEmax), ARHS=ARHS(myDmin:myDmax), &
           &   SOL_GLOBAL=solGlobal, STATUS=status)
    ELSE
      CALL LinearSolverSolve(solver, sol, eMats(:, :, myEmin:myEmax), &
           &   ERHS=erhs(:, myEmin:myEmax), &
           &   SOL_GLOBAL=solGlobal, STATUS=status)
    END IF
  ELSE
    IF (useARHS) THEN
      CALL LinearSolverSolve(solver, sol, eMats(:, :, myEmin:myEmax), &
           &   ERHS=erhs(:, myEmin:myEmax), ARHS=ARHS(myDmin:myDmax), BCVALS=bcVals,&
           &   SOL_GLOBAL=solGlobal, STATUS=status)
    ELSE
      CALL LinearSolverSolve(solver, sol, eMats(:, :, myEmin:myEmax), &
           &   ERHS=erhs(:, myEmin:myEmax), BCVALS=bcVals,&
           &   SOL_GLOBAL=solGlobal, STATUS=status)
    END IF
  END IF
  !
  IF (myRank == 0) THEN
    WRITE(6, *) 'The solution is:  '
    WRITE(6, *) solGlobal
  END IF

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
    !
    ! ============================== Executable Code
    !
    OPEN(UNIT=unit, FILE='input.txt', ACTION='READ')
    !
    CALL ReadComment(unit)  ! problem title/description
    CALL ReadComment(unit)
    READ(unit, *) homoBC, useARHS

    CALL ReadComment(unit)
    READ(unit, *) nDOFpe, nElem, nDOFgl, numBC

    call SetLocalSizes()

    ALLOCATE(conn(nDOFpe, nElem), bcNodes(numBC))
    ALLOCATE(eMats(nDOFpe, nDOFpe, nElem), &
         &   eRhs(nDOFpe, nElem),&
         &   sol(locsizes(myrank)),&
         &   solGlobal(nDOFgl))
    !
    CALL ReadComment(unit)
    READ(unit, *) conn

    CALL ReadComment(unit)
    READ(unit, *) bcnodes

    CALL ReadComment(unit)
    READ(unit, *) eMats
    WRITE(6, *) eMats(:, :, myEmin:myEmax)

    CALL ReadComment(unit)
    READ(unit, *) eRhs
    WRITE(6, *) eRhs(:, myEmin:myEmax)

    IF (.NOT. homoBC) THEN
      ALLOCATE(bcVals(numBC))
      CALL ReadComment(unit)
      READ(unit, *) bcVals
      IF (myRank == 0) THEN
        WRITE(6, '(i0,":  arhs=",5(f12.2,1x))') myRank, bcVals
      END IF
    END IF
    !
    IF (useARHS) THEN
      ALLOCATE(ARHS(nDOFgl))
      CALL ReadComment(unit)
      READ(unit, *) ARHS
      WRITE(6, '(i0,":  arhs=",5f12.2)') myRank, ARHS(myDmin:myDmax)
    END IF
    !
  END SUBROUTINE ReadInput
  !
  ! =================================   END:  ReadInput
  ! ==================================================== BEGIN:  ReadComment
  !
  SUBROUTINE ReadComment(unit)
    !
    !  Read and echo a comment line.
    !
    ! ========== Arguments
    !
    !  .     REQUIRED ARGS
    !
    INTEGER, INTENT(IN) :: unit
    !
    ! ========== Locals
    !
    CHARACTER(LEN=*), PARAMETER :: myNAME = 'ReadComment'
    !
    CHARACTER(LEN=128) :: line
    !
    ! ================================================== Executable Code
    !
    READ(unit, '(a)') line;  
    WRITE(6, '(i0,a)') myrank, ':  '//TRIM(line)
    !
  END SUBROUTINE ReadComment
  !
  ! ====================================================   END:  ReadComment
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
    locsizes(0:(numLeft-1)) = perProc + 1
    IF (myRank == 0) THEN
      WRITE(6, '(i0,a)') myrank, ':  '//'local sizes:  '
      PRINT *, locsizes
    END IF
    !
    !  Set element ranges.
    !
    EperProc = nElem / numProcs
    numLeft  = nElem - EperProc*numProcs
    myEmin = 1; myDmin = 1
    DO p=0, myrank
      IF (myRank == p) THEN
        myEmax = myEmin + EperProc - 1
        IF (p < numLeft) myEmax = myEmax + 1
        !
        myDmax = myDmin + locsizes(p) - 1
      ELSE
        myEmin = myEmin + EperProc
        IF (p < numLeft) myEmin = myEmin + 1
        myDmin = myDmin + locsizes(p)
      END IF
    END DO
    WRITE(6, '(i0,a,1x,i0,1x,i0)') myrank, ':  '//'my elements', myEmin, myEmax
    WRITE(6, '(i0,a,1x,i0,1x,i0)') myrank, ':  '//'my degrees of freedom', myDmin, myDmax
    !
  END SUBROUTINE SetLocalSizes
  !
  ! =================================   END:  SetLocalSizes
  !
END PROGRAM LinearSolver
