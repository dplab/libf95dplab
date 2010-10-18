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
  IMPLICIT NONE
  !
  !  NOTE:  Here is how to write the PETSc vectors to an ascii file.
  !
  !  PetscViewer :: v
  !  CALL PETScViewerASCIIOpen(Petsc_comm_world, "B.dat", v, IERR)
  !  CALL VecView(self % b, v, IERR)
  !  CALL VecView(self % b, v, IERR)
  !  CALL MatView(self % A, v, IERR)
  ! 
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
  DOUBLE PRECISION, PARAMETER :: DEFAULT_DIV_TOL = 1.0d8
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
    !  l_numDOF  -- number of DOF on this (local) process
    !  g_num1    -- starting (global) node number [1-based]
    !
    !  These two are needed to construct global solution.
    !
    !  l_sizes   -- array of local DOF for all processes
    !  g_offSets -- array of global DOF offsets
    !
    INTEGER :: l_numDOF, g_num1
    INTEGER, ALLOCATABLE :: l_sizes(:), g_offSets(:)
    !
    Type(SystemConn) :: sysConn
    LOGICAL          :: symmetry = DFLT_SYMMETRY
    !
    !  Various settings for PETSc solvers.
    !
    DOUBLE PRECISION :: error_tol = DEFAULT_ERROR_TOL, div_tol = DEFAULT_DIV_TOL
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
  INTEGER :: numProcs = 0, myRank = -1
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
  PUBLIC :: LinearSolverPETScInit, LinearSolverPETScFinalize
  PUBLIC :: LinearSolverCreate, LinearSolverSolve
  !
CONTAINS ! ============================================= MODULE PROCEDURES
  !
  ! ==================================================== BEGIN:  LinearSolverPETScInit
  SUBROUTINE LinearSolverPETScInit()
    !
    !  Module initialization.  (PRIVATE)
    !
    !  Call PETScInitialize and determine number of processes and rank.
    !
    !  Don't forget to call PetscFinalize at the end of the program.
    !
    ! ========== Arguments
    !
    ! ***NONE*** 
    !
    CHARACTER(LEN=*), PARAMETER :: myNAME = 'LinearSolverPETScInit'
    !
    ! ========== Locals
    !
    INTEGER :: ierr
    !
    ! ================================================== Executable Code
    !
    call PetscInitialize(PETSC_NULL_CHARACTER, IERR)
    !
    call MPI_Comm_size(PETSC_COMM_WORLD, numprocs,IERR)
    call MPI_Comm_rank(PETSC_COMM_WORLD, myrank,  IERR)
    !
    PETScInitialized = .TRUE.
    !
  END SUBROUTINE LinearSolverPETScInit
  ! ====================================================   END:  LinearSolverPETScInit
  ! ==================================================== BEGIN:  LinearSolverPETScFinalize
  !
  SUBROUTINE LinearSolverPETScFinalize()
    !
    !  Call PETScFinalize
    !
    ! ========== Arguments
    !
    ! ***NONE***
    !
    CHARACTER(LEN=*), PARAMETER :: myNAME = 'LinearSolverPETScFinalize'
    !
    ! ========== Locals
    !
    !
    INTEGER :: IERR
    !
    ! ================================================== Executable Code
    !
    CALL PetscFinalize(IERR)
    !
  END SUBROUTINE LinearSolverPETScFinalize
  !
  ! ====================================================   END:  LinearSolverPETScFinalize
  ! ==================================================== BEGIN:  LinearSolverCreate
  !
  SUBROUTINE LinearSolverCreate(self, name, conn, bcNodes, locSizes, SYMMETRY, STATUS)
    !
    !  Create an instance of a linear solver.
    !
    ! ========== Arguments
    !
    !  .     REQUIRED ARGS
    !
    !  self     -- the linear system type
    !  name     -- the name of the linear system (for messages)
    !  conn     -- 1-based connectivity section (for DOF's of the system)
    !           -- each process passes some section of the connectivity
    !  bcNodes  -- 1-based global list of DOF numbers with essential BC's
    !           -- each process needs to pass the complete list of DOF numbers.
    !  locSizes -- number of degrees of freedom on all processor
    !              (need to know both number of DOF on this processor and starting 
    !               global node number)
    !
    !  .     OPTIONAL ARGS
    !
    !  SYMM     -- symmetric system indicator (optional) 
    !  STATUS   -- return value (optional):  0 for success (not used at the moment)
    !
    TYPE(LinearSolverType), INTENT(IN OUT) :: self
    !
    CHARACTER(LEN=*), INTENT(IN) :: name
    !
    INTEGER, INTENT(IN) :: conn(:, :)
    INTEGER, INTENT(IN) :: bcnodes(:)
    INTEGER, INTENT(IN) :: locsizes(0:)
    !
    LOGICAL, INTENT(IN)  :: SYMMETRY
    INTEGER, INTENT(OUT) :: STATUS
    !
    OPTIONAL :: SYMMETRY, STATUS
    !
    ! ========== Locals
    !
    INTEGER :: iproc, ndof_loc, ndof_glo, aStat, dofpe, numEl, numBC
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
    ndof_loc = locsizes(myrank)
    ndof_glo = SUM(locsizes)
    !
    self % l_numDOF = ndof_loc
    !
    ALLOCATE(self % l_sizes(0:numProcs-1), self % g_offSets(0:numProcs-1))
    self % l_sizes   = locSizes
    self % g_offSets = 0
    !
    DO iproc=1, numProcs - 1
      self % g_offSets(iproc) = self % g_offSets(iproc - 1) + locSizes(iproc - 1)
    end do
    self % g_num1 = self % g_offSets(myRank) + 1
    !
    ! ========== Connectivity and BCs
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
    self % sysConn % bcs = bcnodes - 1 ! make zero-based
    !
    !  This makes the connectivity negative for nodes with BCs.
    !  These are then ignored by PETSc. [checked 2009-04-30, see MatSetValues docs]
    !
    !  NOTE:  I multiply the 1-based node numbers by -1 for nodes WITH BCs and then
    !         subtract one from all node numbers to make the resulting list  0-based.  
    !         Subtracting one before changing sign does not work properly for node 1(0).
    !
    bcSign = 1
    bcSign(bcnodes) = -1
    self % sysConn % sconn = conn*&
         &   RESHAPE(bcsign(RESHAPE(conn, (/dofpe*numel/))), &
         &   (/dofpe, numel/)) - 1 ! changes sign for BC dof
    !
    DEALLOCATE(bcSign)
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
  SUBROUTINE LinearSolverSolve(self, sol, eMats,&
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
    !  self  -- the linear solver instance
    !  sol   -- the local solution vector (section on this process)
    !  eMats -- the local elemental stiffness matrices
    !
    !  .     OPTIONAL ARGS
    !  ERHS       -- the elemental right-hand sides for each process
    !  ARHS       -- already-assembled right-hand side (1-based array);
    !             -- each process sends its own DOFs; both ERHS and
    !             -- ARHS can be passed together, but at least one
    !             -- is required
    !  BCVALS     -- boundary values to apply (usually nonzero);
    !                this is a global array;  only values from 
    !                process 0 are used, but this argument must
    !                be present on all other processes
    !  SOL_GLOBAL -- global solution vector (if present for one process,
    !             -- it must be present for all)
    !  STATUS     -- return value; convergence status from PetSc
    !
    TYPE(LinearSolverType), INTENT(IN OUT) :: self
    !
    DOUBLE PRECISION, INTENT(IN OUT) :: sol(:)
    DOUBLE PRECISION, INTENT(IN)     :: emats(:, :, :)
    !
    DOUBLE PRECISION, INTENT(IN)  :: ERHS(:, :)
    DOUBLE PRECISION, INTENT(IN)  :: ARHS(:)
    DOUBLE PRECISION, INTENT(IN)  :: BCVALS(:)
    DOUBLE PRECISION, INTENT(OUT) :: SOL_GLOBAL(:)
    !
    INTEGER, INTENT(OUT) :: STATUS
    !
    OPTIONAL :: ERHS, ARHS, BCVALS, SOL_GLOBAL, STATUS
    TARGET   :: self  ! so can use pointers on components
    !
    ! ========== Locals
    !
    INTEGER :: IERR, numbc, mySTATUS
    INTEGER, POINTER :: bcn(:)
    LOGICAL :: statPresent
    !
    ! ================================================== Executable Code
    !
    statPresent = PRESENT(STATUS)
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
    if (PRESENT(ARHS)) then
      !
      call SetBFromARHS(self, ARHS)
    else 
      !
      call VecSet   (self % b, 0.0d0, IERR)
      call VecAssemblyBegin(self % b, IERR)
      call VecAssemblyEnd  (self % b, IERR)
      !
    end if
    ! 
    !  . Add in unassembled contribution
    !
    IF (PRESENT(ERHS)) THEN
      CALL SetBFromERHS(self, ERHS)
    END IF
    !
    ! ========== Handle BCs
    !
    !  Nonhomogeneous BCs are handled by creating a RHS that satisfies
    !  the BCs, but not the system.  Then, subtracting, you get a homeogeneous
    !  system.
    !
    bcn   => self % sysconn % bcs
    numbc =  SIZE(bcn, 1)
    !
    if (PRESENT(BCVALS)) then
      !
      !  Set up the homogeneous problem:  A(x-x0) = b-b0
      !
      call FormA(self, EMATS, HOMOGENEOUS_BCS=.FALSE.)
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
      call MatMult(self % A, self % x0,   self % ax0, IERR)
      call VecAXPY(self % b,    -1.0d0,   self % ax0, IERR)
      !
    end if
    !
    !  Zero nodes with applied BCs.
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
    ! ========== Solve System
    !
    !  A is now formed for a homogeneous problem.
    !
    !  * Note:  mySTATUS from SolveSystem is PETSc reason (good if >= 0)
    !
    call FormA(self, eMats, HOMOGENEOUS_BCS=.TRUE.)
    call SolveSystem(self, mySTATUS)
    !
    IF (mySTATUS < 0) THEN
      call SetStatus(RETURN_FAILURE); RETURN 
    ELSE
      call SetStatus(RETURN_SUCCESS)
    end if
    !
    !  At this point, the solution succeeded
    !
    if (PRESENT(BCVALS)) then ! add back in x0
      CALL VecAXPY(self % x, 1.0d0, self % x0, IERR) ! y = y + a * x
    endif
    !
    !  Copy solution into local copy for f90, and get
    !  global solution if requested.
    !
    call GetLocal()
    !
    IF (PRESENT(SOL_GLOBAL)) THEN
      call GetGlobal()
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
      IF (statPresent) THEN
        STATUS = stat
      END IF
      !
    END SUBROUTINE SetStatus
    !
    ! ================================   END:  SetStatus
    ! ================================ BEGIN:  GetLocal
    !
    SUBROUTINE GetLocal()
      !
      !  Get local f90 solution vector.
      !
      ! ========== Arguments
      !

      !
      ! ========== Locals
      !
      INTEGER :: IERR

      PetscScalar, POINTER :: xx_v(:) => NULL()

      INTEGER :: iTmp
      DOUBLE PRECISION :: vMin, vMax
      !
      ! ============================== Executable Code
      !
      !ALLOCATE(xx_v(self % l_numDOF))
      call VecGetArrayF90(self % x, xx_v, IERR)
      sol = xx_v
      call VecRestoreArrayF90(self % x, xx_v, IERR)
      !DEALLOCATE(xx_v)
      !
    END SUBROUTINE GetLocal
    !
    ! =================================   END:  GetLocal
    ! ================================= BEGIN:  GetGlobal
    !
    SUBROUTINE GetGlobal()
      !
      !  Set global solution.
      !
      !
      ! ========== Locals
      !
      INTEGER :: IERR
      !
      ! ============================== Executable Code
      !
      CALL MPI_AllGatherV(&
           &  sol, self % l_numDOF, MPI_DOUBLE_PRECISION, &
           &  SOL_GLOBAL, self % l_sizes, self % g_offSets, MPI_DOUBLE_PRECISION,&
           &  PETSC_COMM_WORLD, IERR)
      !
    END SUBROUTINE GetGlobal
    !
    ! =================================   END:  GetGlobal

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
      n = self % g_num1 + i - 2  ! make 0-based
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
  SUBROUTINE SetBfromERHS(self, Erhs)
    !
    !  Set the PETSc right-hand side from elemental right-hand sides.
    !
    !  NOTES:  
    !
    !  *  This appends to existing values.
    !
    ! ========== Arguments
    !
    !  .     REQUIRED ARGS
    !  self -- the LinearSolverType
    !  Erhs -- elemental right-hand side (references a 1-based array)
    !
    TYPE(LinearSolverType), INTENT(IN OUT) :: self
    !
    DOUBLE PRECISION, INTENT(IN) :: Erhs(:, :)
    !
    TARGET :: self 
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
  END SUBROUTINE SetBfromERHS
  !
  ! ====================================================   END:  SetBfromERHS
  ! ==================================================== BEGIN:  FormA
  !
  SUBROUTINE FormA(self, emats, HOMOGENEOUS_BCS)
    !
    !  Form PETSc matrix.
    !
    ! ========== Arguments
    !
    !  .     REQUIRED ARGS
    !  self  -- the LinearSolver instance
    !  eMats -- elemental matrices
    !
    !  HOMOGENEOUS_BCS 
    !        -- flag indicating homogeneous (zero) boundary conditions
    !        -- if true, ignore matrix contributions from these nodes;
    !        -- otherwise, include them, as they will need to be backsolved
    !
    !
    TYPE(LinearSolverType), INTENT(IN OUT) :: self
    !
    DOUBLE PRECISION, INTENT(IN) :: emats(:, :, :)
    LOGICAL, INTENT(IN) :: HOMOGENEOUS_BCS
    !
    TARGET :: self
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

    call MatZeroEntries(self % A, ierr)

    if (HOMOGENEOUS_BCS) then
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
      !  Set the diagonal to one.  Only one process need do it.
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
      !  NONHOMOGENEOUS BCS
      !
      !  Use absolute value of connectivity to ensure
      !  all DOF are assembled.
      !
      do elem=elmin, elmax
        call MatSetValues(self % A, &
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
         &   self % div_tol, &
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
    IF (myrank == 0) THEN
      PRINT *, 'Linear system "'//TRIM(ADJUSTL(self % name))//'"'
      PRINT *, '       iterations ', its, ', code ', reason
    END IF

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
