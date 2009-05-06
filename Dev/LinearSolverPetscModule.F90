!  $Id$
!
MODULE LinearSolverPETScModule 
  !
  !  This module interfaces the PETSc solver with our
  !  parallel finite element codes.
  !
  !  NOTES:
  !
  ! * Max iterations for a system is set to twice the total degrees of freedom.
  !
  USE IntrinsicTypesModule,  RK=>REAL_KIND, IK=>INTEGER_KIND, LK=>LOGICAL_KIND
  !
  IMPLICIT NONE
  !
  ! ==================== PETSc Includes
  !
#include "include/finclude/petsc.h"
#include "include/finclude/petscvec.h"
#include "include/finclude/petscvec.h90"
#include "include/finclude/petscmat.h"
#include "include/finclude/petscpc.h"
#include "include/finclude/petscksp.h"

#include "include/finclude/petscsys.h"
#include "include/finclude/petscviewer.h"
  !
  ! ==================== Config
  !
  LOGICAL,          PARAMETER :: DFLT_SYMMETRY = .TRUE.
  INTEGER,          PARAMETER :: PREALLOCATION_FACTOR = 4 ! for estimating DOF per row
  DOUBLE PRECISION, PARAMETER :: DEFAULT_ERROR_TOL = 1.0d-12
  !
  ! ========== Return Values
  !
  INTEGER, PARAMETER :: RETURN_SUCCESS = 0, RETURN_FAILURE = 1  ! generic errors
  !
  ! ==================== Types
  !
  !  Helper type for connectivity and bcs.
  !  ------------------------------------
  !
  TYPE SystemConn
    !  Signed Connectivity
    INTEGER, ALLOCATABLE :: sconn(:, :), bcs(:)
  END TYPE SystemConn
  !
  !  Main type for linear solver.
  !  ---------------------------
  !
  TYPE LinearSolverType
    !
    CHARACTER(LEN=64) :: name
    !
    !  PETSc types
    !
    Vec :: x, b, x0, ax0
    Mat :: A
    !
    !  l_numDOF -- number of DOF on this (local) process
    !  g_num0   -- starting (global) node number [1-based]
    !
    INTEGER :: l_numDOF, g_num0
    !
    Type(SystemConn) :: sysConn
    LOGICAL          :: symmetry = DFLT_SYMMETRY
    !
    !  Various settings for PETSc solvers.
    !
    DOUBLE PRECISION :: error_tol = DEFAULT_ERROR_TOL
    INTEGER          :: maxits
    !
  END TYPE LinearSolverType
  !
  ! ==================== Module Data
  !
  LOGICAL :: PETScInitialized = .FALSE.
  !
  !  Empty array for passing as unused argument.
  !
  DOUBLE PRECISION :: NULL2D(0, 0)
  !
  !  PETSc/Parallel items.
  !
  INTEGER :: numprocs = 0, myrank = -1
  !
  ! ==================== Public Entities
  !
  PRIVATE
  !
  ! ========== Data
  !
  PUBLIC :: numprocs, myrank
  PUBLIC :: RETURN_SUCCESS, RETURN_FAILURE
  !
  ! numprocs -- number of processes
  ! myrank   -- rank of current process (myid in parallel_mod)
  !
  ! ========== Derived Types
  !
  PUBLIC :: LinearSolverType
  !
  ! ========== Procedures
  !
  PUBLIC :: LinearSolverCreate, LinearSolverSolve
  !
CONTAINS ! ============================================= MODULE PROCEDURES
  !
  ! ==================================================== BEGIN:  LinearSolverPETScInit
  SUBROUTINE LinearSolverPETScInit(status)
    !
    !  Module initialization.  (PRIVATE)
    !
    !  Call PETScInitialize and determine number of processes and rank.
    !
    !  Don't forget to call PetscFinalize at the end of the program.
    !
    ! ========== Arguments
    !
    !  status -- return value:  0 for success
    !
    INTEGER, INTENT(OUT) :: status
    !
    ! ========== Locals
    !
    INTEGER :: ierr
    !
    ! ================================================== Executable Code
    !
    call PetscInitialize(PETSC_NULL_CHARACTER, status)
    !
    call MPI_Comm_size(PETSC_COMM_WORLD, numprocs,ierr)
    call MPI_Comm_rank(PETSC_COMM_WORLD, myrank,  ierr)
    !
    PETScInitialized = .TRUE.
    !
  END SUBROUTINE LinearSolverPETScInit
  ! ====================================================   END:  LinearSolverPETScInit
  ! ==================================================== BEGIN:  LinearSolverCreate
  !
  SUBROUTINE LinearSolverCreate(self, name, conn, bcnodes, locsizes, SYMM)
    !
    !  Create an instance of a linear solver.
    !
    ! ========== Arguments
    !
    !  name     -- the name of the linear system (for messages)
    !  locsizes -- number of degrees of freedom on each processor
    !  conn     -- 1-based connectivity section (for DOF's of the system)
    !           -- each process passes some section of the connectivity
    !  bcnodes  -- 1-based global list of DOF numbers with essential BC's
    !           -- each process needs to pass the complete list of DOF numbers.
    !  locsizes -- number of degrees of freedom on each processor
    !  SYMM     -- symmetric system indicator (optional) 
    !  status   -- return value (optional):  0 for success (not used at the moment)
    !
    TYPE(LinearSolverType), INTENT(IN OUT) :: self
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    !
    INTEGER, INTENT(IN) :: conn(:, :)
    INTEGER, INTENT(IN) :: bcnodes(:)
    INTEGER, INTENT(IN) :: locsizes(0:)
    LOGICAL, INTENT(IN), OPTIONAL :: SYMM
    !
    INTEGER, INTENT(OUT), OPTIONAL :: STATUS
    !
    ! ========== Locals
    !
    INTEGER :: iproc, ndof_loc, ndof_glo, aStat
    INTEGER :: rowAlloc, IERR
    !
    INTEGER, ALLOCATABLE :: bcsign(:)
    !
    ! ================================================== Executable Code
    !
    !  Initialize module, if not already.
    !
    IF (.NOT. PETScInitialized) THEN
      CALL LinearSolverPETScInit()
    END IF
    !
    !  Set status.
    !
    IF (PRESENT(STATUS)) THEN
      STATUS = 0
    END IF
    !
    !  Name and symmetry.
    !
    self % name = name
    !
    IF (PRESENT(SYMMETRY)) THEN
      self % symmetry = SYMMETRY
    END IF
    !
    ! ========== Local array sizes and ranges
    !
    ALLOCATE(self % localsizes (0:(numprocs-1)), &
         &   self % section_beg(0:(numprocs-1))  )    
    !
    self % localsizes  = locsizes
    self % section_beg(0) = 1
    do iproc=1, numprocs - 1
      self % section_beg(iproc) = self % section_beg(iproc - 1) + locsizes(iproc - 1)
    end do
    !
    ! ========== Connectivity and BCs
    !
    ndof_glo = SUM(locsizes)
    ndof_loc = locsizes(myrank)
    !
    dofpe  = SIZE(conn, 1)
    numel  = SIZE(conn, 2)
    numbc  = SIZE(bcnodes)
    !
    ALLOCATE(self % sysConn % bcs(numbc), &
         &   self % sysConn % sconn(dofpe, numel),&
         &   bcsign(ndof_glo),&
         &   STAT=aStat)
    ! TODO:  check aStat and deallocate if allocation failed
    !
    self % sysconn % bcs = bcnodes - 1 ! make zero-based
    !
    !  This makes the connectivity negative for nodes with BCs.
    !  These are then ignored by PETSc. [checked 2009-04-30, see MatSetValues docs]
    !
    !  NOTE:  I multiply the 1-based node numbers by -1 for nodes WITH BCs and then
    !         subtract one from all node numbers to make the resulting list  0-based.  
    !         Subtracting one before changing sign does not work properly for node 1(0).
    !
    bcsign = 1
    bcsign(bcnodes) = -1
    sysconn % sconn = conn*&
         &   RESHAPE(bcsign(RESHAPE(conn, (/dofpe*numel/))), &
         &   (/dofpe, numel/)) - 1 ! changes sign for BC dof
    !
    DEALLOCATE(bcsign)
    !
    ! ========== Make PETSc objects
    !
    ! Matrix A
    !
    rowAlloc = PREALLOCATION_FACTOR*dofpe
    call MatCreateMPIAIJ(PETSC_COMM_WORLD, ndof_loc, ndof_loc,&
         &   ndof_glo, ndof_glo, &
         &   rowAlloc, PETSC_NULL_INTEGER, rowAlloc, PETSC_NULL_INTEGER, &
         &   self % A, IERR)
    !
    if (self % symmetry) then
      CALL MatSetOption(self % A, MAT_SYMMETRIC, IERR)
    end if
    !
    ! Vectors x, b, x0, ax0
    !
    call VecCreateMPI(PETSC_COMM_WORLD, ndof_loc, ndof_glo, self % x, IERR)
    call VecDuplicate(self % x, self % b, IERR)
    call VecDuplicate(self % x, self % x0, IERR)
    call VecDuplicate(self % x, self % ax0, IERR)
    !
    ! ========== Other parameters
    !
    self % maxits = 2 * ndof_glo
    !
  END SUBROUTINE LinearSolverCreate
  !
  ! ==================================================== LinearSolverCreate:  END
  ! ==================================================== BEGIN:  LinearSolverSolve
  !
  SUBROUTINE LinearSolverSolve(sol, eMats,&
       &   ERHS, ARHS, &
       &   BCVALS, SOL_GLOBAL, STATUS)
    !
    !  Solve a linear system using PETSc routines.
    !
    !  This can be called with both RHS args, e.g.  ERHS and ARHS, 
    !  as in isaiah/vplas().  BCVALS are optional in the sense that
    !  they are set to zero if not explicitly passed.
    !
    ! ========== Arguments
    !
    !  .     REQUIRED ARGS
    !  sol   -- the local solution vector (section on this process)
    !  eMats -- the local elemental stiffness matrices
    !
    !  .     OPTIONAL ARGS
    !  ERHS       -- the local elemental right-hand sides
    !  ARHS       -- assembled right-hand side (1-based array);
    !             -- mutually exclusive with ERHS
    !  BCVALS     -- boundary values to enforce (usually nonzero)
    !  SOL_GLOBAL -- global solution vector
    !  STATUS     -- return value; convergence status from PetSc
    !
    DOUBLE PRECISION, INTENT(IN OUT) :: sol(:)
    DOUBLE PRECISION, INTENT(IN)     :: emats(:, :, :)
    !
    DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: ERHS(:, :)
    DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: ARHS(:)
    DOUBLE PRECISION, INTENT(IN),  OPTIONAL :: BCVALS(:)
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: SOL_GLOBAL(:)
    !
    INTEGER, INTENT(OUT), OPTIONAL :: STATUS
    !
    ! ========== Locals
    !
    INTEGER :: IERR, numbc, mySTATUS
    INTEGER, POINTER :: bcn(:)
    !
    ! ================================================== Executable Code
    !
    CALL SetStatus(RETURN_SUCCESS)
    !
    IF (.NOT. PRESENT(ERHS) .AND. .NOT. PRESENT(ARHS)) THEN 
      ! Neither are present.
      CALL SetStatus(RETURN_FAILURE);  RETURN
    END IF
    !
    ! ========== Form RHS
    !
    !  . If portion is assembled already, set to that; otherwise zero.
    !
    if (PRESENT(arhs)) then
      !
      call SetBFromARHS(self, arhs)
    else 
      !
      call VecSet   (self % b, 0.0d0, IERR)
      call VecAssemblyBegin(self % b, IERR)
      call VecAssemblyEnd  (self % b, IERR)
      !
    end if
    ! 
    !  . Add in unassembled contribution
    !  -- if nonzero BCs, do not ignore contributions of those nodes, since
    !     they will be backsolved
    !
    IF (PRESENT(ERHS)) THEN
      CALL BFromERHS(self, ERHS, USE_ALL_DOF=PRESENT(BCVALS))
    END IF
    !
    ! ========== Handle BCs
    !
    !  Nonhomogeneous BCs are handled by creating a RHS that satisfies
    !  the BCs, but not the system.  Then, subtracting, you get a homeogeneous
    !  system.
    !
    if (PRESENT(BCVALS)) then
      !
      bcn   => self % sysconn % bcs
      numbc =  SIZE(bcn, 1)
      !
      call FormA(sysid, EMATS, APPLY_BCS=.FALSE.)
      !
      !  Form x0 satisfying bcs.
      !
      call VecSet(self % x0, 0.0d0, IERR)
      !
      IF (myrank == 0) THEN
        call VecSetValues(self % x0, numbc, bcn, BCVALS, &
             &  INSERT_VALUES, IERR)
      end if
      !
      call VecAssemblyBegin(self % x0, IERR)
      call VecAssemblyEnd  (self % x0, IERR)
      !
      call MatMult(self %   A, self % x0,   self % ax0, IERR)
      call VecAXPY(self % ax0,    -1.0d0,   self % b,   IERR)
      !
      !  Prepare to solve homogeneous problem:  A(x-x0) = b-b0
      !
      if (myrank == 0) then
        call VecSetValues(self % b, numbc, bcn, &
             &  SPREAD(0.0d0, DIM=1, NCOPIES=numbc), &
             &  INSERT_VALUES, IERR)
      end if
      !
      call VecAssemblyBegin(self % b, IERR)
      call VecAssemblyEnd  (self % b, IERR)
      !
    end if
    !
    ! ========== Solve System
    !
    !  A is now formed for a homogeneous problem.
    !
    call FormA(self, eMats, APPLY_BCS=.TRUE.)
    call SolveSystem(self, mySTATUS)
    !
    IF (mySTATUS /= RETURN_SUCCESS) THEN
      call SetStatus(RETURN_FAILURE); RETURN 
    end if
    !
    !  At this point, the solution succeeded
    !
    if (PRESENT(BCVALS)) then ! add back in x0
      CALL VecAXPY(self % x, 1.0d0, self % x0, IERR)
    endif
    !
    IF (PRESENT(SOL_GLOBAL)) THEN
      call SetGlobal()
    END IF

  CONTAINS ! ========================= Internal Procedures

    ! ================================ BEGIN:  SetStatus
    !
    SUBROUTINE SetStatus(stat)
      !
      !  Set status, if present
      !
      ! ========== Arguments
      !
      INTEGER, INTENT(IN) :: stat
      !
      ! ============================== Executable Code
      !
      IF (PRESENT(STATUS)) THEN
        STATUS = stat
      END IF
      !
    END SUBROUTINE SetStatus
    !
    ! ================================   END:  SetStatus
    ! ================================ BEGIN:  SetGlobal
    !
    SUBROUTINE SetGlobal()
      !
      !  Set global solution.
      !
      !
      ! ========== Locals
      !

      !
      ! ============================== Executable Code
      !
    INTEGER :: ierr, dofpe, elmin, elmax, 
    INTEGER :: elem, node, its, ibc, vecsize, i, iproc
    !
    PetscScalar  :: xptr(1)
    PetscOffset  :: offs
    PetscScalar, POINTER :: xx_v(:) => NULL()
    !DOUBLE PRECISION, ALLOCATABLE :: xx_v(:)
    !
    INTEGER :: mpistat(MPI_STATUS_SIZE), mystatus
    !
    !  *** End:
    INTEGER :: iTmp
    DOUBLE PRECISION :: vMin, vMax

      vecsize = localsizes(myrank, sysid)
      !allocate(xx_v(vecsize))
      call VecGetArrayF90(self % x, xx_v, ierr)
      PRINT *, '*** GetArrayF90:  ', ierr

      CALL VecMin(self % x, iTmp, vMin, ierr)
      CALL VecMax(self % x, iTmp, vMax, ierr)
      PRINT *, '*** max(x):  ', vMin, vMax
      PRINT *, '*** max(xx_v):  ', MAXVAL(xx_v), MINVAL(xx_v)

      sol = xx_v
      call VecRestoreArrayF90(self % x, xx_v, ierr)
      PRINT *, '*** RestoreArrayF90:  ', ierr
      !deallocate(xx_v)
      WRITE(99, *) '---sol---'
      WRITE(99, *) sol
      !
      !using VecGetArray!!
      !using VecGetArray!vecsize = localsizes(myrank, sysid)
      !using VecGetArray!call VecGetArray(self % x, xptr, offs, ierr)
      !using VecGetArray!do i=1, vecsize
      !using VecGetArray!   sol(i) = xptr(offs + i)
      !using VecGetArray!end do
      !using VecGetArray!call VecRestoreArray(self % x, xptr, offs, ierr)
      !
      !  Now construct global solution if requested.
      !
      SOL_GLOBAL = sol
      !
      !
    END SUBROUTINE SetGlobal
    !
    ! =================================   END:  SetGlobal

    !
  END SUBROUTINE LinearSolverSolve
  !
  ! ====================================================   END:  LinearSolverSolve
  ! ==================================================== BEGIN:  SetBfromARHS
  !
  SUBROUTINE SetBfromARHS(self, arhs)
    !
    !  (PRIVATE) Set the PETSc right-hand side, b, from assembled RHS, arhs.
    !
    ! ========== Arguments
    !
    !  .     REQUIRED ARGS
    !  self -- the LinearSolver instance
    !  arhs -- assembled right-hand side (1-based array, local to process)
    !
    TYPE(LinearSolverType), INTENT(IN OUT) :: self
    DOUBLE PRECISION,       INTENT(IN)     :: arhs(:)
    !
    ! ========== Locals
    !
    INTEGER :: i, n, ierr
    !
    ! ================================================== Executable Code
    !
    do i=1, self % l_numDOF
      n = self % g_num0 + i - 2  ! make 0-based
      CALL VecSetValues(self % b, 1, n, arhs(i), INSERT_VALUES, IERR)
    end do
    !
    call VecAssemblyBegin(self % b, IERR)
    call VecAssemblyEnd  (self % b, IERR)
    !
  END SUBROUTINE SetBfromARHS
  !
  ! ====================================================   END:  SetBfromARHS
  ! ==================================================== BEGIN:  SetBfromERHS
  !
  SUBROUTINE SetBfromERHS(self Erhs, USE_ALL_DOF)
    !
    !  Set the PETSc right-hand side from elemental right-hand sides.
    !
    !  NOTES:  
    !
    !  *  This appends to existing values.
    !  *  CHECK THIS:  assembles all nodes in both cases, but only applies
    !     BCs if USE_ALL_DOF is false.
    !
    ! ========== Arguments
    !
    !  .     REQUIRED ARGS
    !  self -- the LinearSolverType
    !  Erhs -- elemental right-hand side (references a 1-based array)
    !  USE_ALL_DOF -- flag to assemble all DOF regardless of BC
    !
    !  .     OPTIONAL ARGS
    !
    TYPE(LinearSolverType), INTENT(IN OUT) :: self
    !
    DOUBLE PRECISION, INTENT(IN) :: Erhs(:, :)
    LOGICAL, INTENT(IN) :: USE_ALL_DOF
    !
    ! ========== Locals
    !
    INTEGER :: dofpe, elem, elmin, elmax, numbc, IERR
    INTEGER, POINTER :: sconn(:, :), bcn(:)
    !
    ! ================================================== Executable Code
    !
    dofpe = SIZE  (erhs, 1)
    elmin = LBOUND(erhs, 2)
    elmax = UBOUND(erhs, 2)
    !
    sconn => self % sysconn % sconn
    bcn   => self % sysconn % bcs
    numbc = SIZE(bcn, 1)
    !
    if (USE_ALL_DOF) then
      !
      !  Assemble all nodes, regardless of applied BCs.
      !
      do elem=elmin, elmax
        call VecSetValues(self % b, dofpe, abs(1 + sconn(:, elem)) - 1, &
             &  erhs(:, elem), &
             &  ADD_VALUES, IERR)
      end do
      !
      call VecAssemblyBegin(self % b, IERR)
      call VecAssemblyEnd  (self % b, IERR)
      !
    else
      !
      !  Only assemble nodes without applied BCs.
      !
      do elem=elmin, elmax
        call VecSetValues(self % b, dofpe, abs(1 + sconn(:, elem)) - 1, &
             &  erhs(:, elem), &
             &  ADD_VALUES, IERR)
        !
      end do
      !
      call VecAssemblyBegin(self % b, IERR)
      call VecAssemblyEnd  (self % b, IERR)
      !
      if (myrank == 0) then
        CALL VecSetValues(self % b, numbc, bcn, &
             &  SPREAD(0.0d0, DIM=1, NCOPIES=numbc), &
             &  INSERT_VALUES, IERR)
      end if
      !
      call VecAssemblyBegin(self % b, IERR)
      call VecAssemblyEnd  (self % b, IERR)
      !
    end if
    !
  END SUBROUTINE SetBfromERHS
  !
  ! ====================================================   END:  SetBfromERHS
  ! ==================================================== BEGIN:  FormA
  !
  SUBROUTINE FormA(self, emats, APPLY_BCS)
    !
    !  Form PETSc matrix.
    !
    ! ========== Arguments
    !
    !  .     REQUIRED ARGS
    !  self  -- the LinearSolver instance
    !  eMats -- elemental matrices
    !  APPLY_BCS -- flag for applying boundary conditions
    !            -- if true, set up the matrix for applied boundary conditions;
    !            -- otherwise, assemble the matrix for all degrees of freedom.
    !
    !  .     OPTIONAL ARGS
    !
    TYEP(LinearSolverType), INTENT(IN OUT) :: self
    !
    DOUBLE PRECISION, INTENT(IN) :: emats(:, :, :)
    LOGICAL, INTENT(IN) :: APPLY_BCS
    !
    ! ========== Locals
    !
    INTEGER :: dofpe, elem, elmin, elmax, numbc, ibc, IERR
    INTEGER, POINTER :: sconn(:, :), bcn(:) ! pointers to strucutre components
    !
    ! ================================================== Executable Code
    !
    dofpe = SIZE(emats, 1)
    elmin = LBOUND(emats, 3)
    elmax = UBOUND(emats, 3)
    !
    sconn => self % sysconn % sconn
    bcn   => self % sysconn % bcs
    numbc = SIZE(bcn, 1)

    call MatZeroEntries(A(sysid), ierr)

    if (APPLY_BCS) then
      !
      !  Negative DOF number is ignored by PETSc.
      !
      do elem=elmin, elmax
        call MatSetValues(self % A, &
             &   dofpe, sconn(:, elem),&
             &   dofpe, sconn(:, elem),&
             &   TRANSPOSE(emats(:, :, elem)), ADD_VALUES, IERR)
      end do
      !
      !  Only one process need apply the BCs.
      !
      if (myrank == 0) then
        do ibc=1, numbc
          call MatSetValues(self % A, &
               &  1, bcn(ibc), 1, bcn(ibc), &
               &  1.0d0, ADD_VALUES, ierr)
        end do
      end if
      !
    else
      !
      !  Use absolute value of connectivity to ensure
      !  all DOF are assembled.
      !
      do elem=elmin, elmax
        call MatSetValues(A(sysid), &
             &   dofpe, ABS(1 + sconn(:, elem)) - 1,&
             &   dofpe, ABS(1 + sconn(:, elem)) - 1,&
             &   TRANSPOSE(emats(:, :, elem)), ADD_VALUES, IERR)
      end do
    endif
    !
    call MatAssemblyBegin(self % A, MAT_FINAL_ASSEMBLY, IERR)
    call MatAssemblyEnd  (self % A, MAT_FINAL_ASSEMBLY, IERR)
    !
  END SUBROUTINE FormA
  !
  ! ====================================================   END:  FormA
  ! ==================================================== BEGIN:  SolveSystem
  !
  SUBROUTINE SolveSystem(self, STATUS)
    !
    !  Solve the PETSc system.
    !
    !  A and b are already constructed.  Currently this
    !  only used Jacobi preconditioner.
    !
    ! ========== Arguments
    !
    !  .     REQUIRED ARGS
    !  self -- the LinearSolverType 
    !
    !  .     OPTIONAL ARGS
    !  STATUS -- return value 
    !
    TYPE(LinearSolverType), INTENT(IN OUT) :: self
    !
    INTEGER, INTENT(OUT), OPTIONAL :: STATUS
    !
    ! ========== Locals
    !
    !  PETSc objects.
    !
    KSP         ksp
    KSPType     ksptype
    PC          pc
    PCType      ptype
    !
    KSPConvergedReason reason
    !PetscViewer viewer 
    !
    INTEGER :: IERR, its
    !
    ! ================================================== Executable Code
    !
    !
    !  Form PETSc solver context.
    !
    CALL KSPCreate(PETSC_COMM_WORLD, ksp, IERR)
    call KSPSetOperators(ksp, self % A, self % A, DIFFERENT_NONZERO_PATTERN, IERR)
    !
    !
    !  Set PETSc options.
    !
    !  . ksptype (Krylov space method)
    !
    if (self % symmetry) then
      ksptype = KSPCG
    else
      ksptype = KSPBICG
    end if
    call KSPSetType(ksp, ksptype, IERR)
    !
    !  . pctype (Preconditioner)
    ! 
    ptype = PCJACOBI
    CALL KSPGetPC(ksp, pc, IERR)
    call PCSetType(pc, ptype, IERR)
    !
    ! . error tolerance
    !
    call KSPSetTolerances(ksp, self % error_tol, &
         &   PETSC_DEFAULT_DOUBLE_PRECISION,&
         &   PETSC_DEFAULT_DOUBLE_PRECISION, &
         &   self % maxIts, IERR)
    !
    !  . command line options
    !
    CALL KSPSetFromOptions(ksp, IERR)
    !
    !  Solve system and check result.
    !
    CALL KSPSetUp(ksp, IERR)
    call KSPSolve(ksp, self % b, self % x, IERR) ! solve!

    call KSPGetConvergedReason(ksp, reason, IERR)
    CALL KSPGetIterationNumber(ksp, its, IERR)
    !
    PRINT *, 'Linear system '//TRIM(ADJUSTL(self % name))// ':  iterations ', its, , ', code ', reason
    IF (reason < 0) THEN
      CALL DescribeDivergedReason(0, reason)
    END IF
    !
    if (PRESENT(STATUS)) then
      STATUS = reason
    end if
    !
    call KSPDestroy(ksp, IERR)
    !
  END SUBROUTINE SolveSystem
  !
  ! ====================================================   END:  SolveSystem
  ! ==================================================== BEGIN:  DescribeDivergedReason
  SUBROUTINE DescribeDivergedReason(u, reason)
    !
    !  Print text for error code.
    !
    !  u      -- unit to write to
    !  reason -- PETSc failure code
    !
    ! ========== Arguments
    !
    INTEGER, INTENT(IN) :: u, reason    
    !
    ! ========== Locals
    !
    CHARACTER(LEN=*), PARAMETER :: eMessage = '**** PETSc solver failed:  '
    CHARACTER(LEN=100)          :: divReason
    !
    ! ================================================== Executable Code
    !
    SELECT CASE(reason)   
      !
    CASE (KSP_DIVERGED_NULL)   
      divReason = 'KSP_DIVERGED_NULL'
    CASE (KSP_DIVERGED_ITS)
      divReason = 'Ran out of iterations before any convergence criteria was reached '
    CASE (KSP_DIVERGED_DTOL)   
      divReason = 'norm(r) >= dtol*norm(b)'
    CASE (KSP_DIVERGED_BREAKDOWN)   
      divReason = 'A breakdown in the Krylov method was detected so the method &
           &could not continue to enlarge the Krylov space. '
    CASE (KSP_DIVERGED_BREAKDOWN_BICG)   
      divReason = 'A breakdown in the KSPBICG method was detected so the method &
           &could not continue to enlarge the Krylov space.'
    CASE (KSP_DIVERGED_NONSYMMETRIC)   
      divReason = 'It appears the operator or preconditioner is not symmetric and &
           &this Krylov method (KSPCG, KSPMINRES, KSPCR) requires symmetry '
    CASE (KSP_DIVERGED_INDEFINITE_PC)   
      divReason = 'It appears the preconditioner is indefinite'
    CASE (KSP_DIVERGED_NAN)   
      divReason = 'residual norm became Not-a-number likely do to 0/0)'
    CASE (KSP_DIVERGED_INDEFINITE_MAT)   
      divReason = 'KSP_DIVERGED_INDEFINITE_MAT'
      !
    CASE DEFAULT    
      divReason = 'unknown reason'
      !
    END SELECT

    WRITE(u, '(a)') eMessage // TRIM(divReason)
    !
    ! ==================== Error Codes
    !
    !  CONVERGENCE REASONs
    !typedef enum {/* converged */
    !              KSP_CONVERGED_RTOL               =  2,
    !              KSP_CONVERGED_ATOL               =  3,
    !              KSP_CONVERGED_ITS                =  4,
    !              KSP_CONVERGED_CG_NEG_CURVE       =  5,
    !              KSP_CONVERGED_CG_CONSTRAINED     =  6,
    !              KSP_CONVERGED_STEP_LENGTH        =  7,
    !              KSP_CONVERGED_HAPPY_BREAKDOWN    =  8,
    !              /* diverged */
    !              KSP_DIVERGED_NULL                = -2,
    !              KSP_DIVERGED_ITS                 = -3,
    !              KSP_DIVERGED_DTOL                = -4,
    !              KSP_DIVERGED_BREAKDOWN           = -5,
    !              KSP_DIVERGED_BREAKDOWN_BICG      = -6,
    !              KSP_DIVERGED_NONSYMMETRIC        = -7,
    !              KSP_DIVERGED_INDEFINITE_PC       = -8,
    !              KSP_DIVERGED_NAN                 = -9,
    !              KSP_DIVERGED_INDEFINITE_MAT      = -10,
    ! 
    !              KSP_CONVERGED_ITERATING          =  0} KSPConvergedReason;
    !
  END SUBROUTINE DescribeDivergedReason
  ! ====================================================   END:  DescribeDivergedReason

END MODULE LinearSolverPETScModule
