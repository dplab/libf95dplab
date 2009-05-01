!  $Id$
!
MODULE LinearSolverPetscModule 
  !
  !  This module interfaces the Petsc solver with our
  !  parallel finite element codes.
  !
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
    !  PETSc types
    !
    Vec :: x, b, x0, ax0
    Mat :: A
    !
    !  Both arrays 'localsizes' and 'section_beg' are 0-based on 
    !  the number of processes.  The array 'localsizes' gives number 
    !  of DOF for each process; the array 'section_beg' gives the 1-based 
    !  index of the first DOF on each process.
    !
    INTEGER, ALLOCATABLE :: localsizes(:)  ! (numprocs, MAXSYS) 
    INTEGER, ALLOCATABLE :: section_beg(:) ! (numprocs, MAXSYS) 
    !
    Type(SystemConn) :: sysConn
    LOGICAL          :: symmetry = DFLT_SYMMETRY
    !
  END TYPE LinearSolverType
  !
  ! ==================== Module Data
  !
  LOGICAL :: PetscInitialized = .FALSE.
  !
  !
  !  Module variables.
  !
  LOGICAL :: sys_used(MAXSYS) = .FALSE.
  DOUBLE PRECISION :: error_tol = DEFAULT_ERROR_TOL

  !
  KSP         ksp
  KSPType     ksptype
  PC          pc
  PCType      ptype
  KSPConvergedReason reason
  PetscViewer viewer
  !
  !  Petsc items.
  !
  INTEGER :: numprocs = 0, myrank = -1
  !
  !  Empty array for passing as unused argument.
  !
  DOUBLE PRECISION :: NULL2D(0, 0)
  !
  ! ==================== Public Entities
  !
  PRIVATE
  !
  ! ========== Data
  !
  PUBLIC :: numprocs, myrank
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
  PUBLIC :: LinearSolverCreate
  !
CONTAINS ! ================================================== Module Procedures
  !
  ! ==================================================== BEGIN:  LinearSolverPetscInit
  SUBROUTINE LinearSolverPetscInit(status)
    !
    !  Module initialization.  
    !
    !  Call PetscInitialize and determine number of processes and rank.
    !
    !  Don't forget to call PetscFinalize at the end of the program.
    !
    ! ========== Arguments
    !
    INTEGER, INTENT(OUT) :: status
    !
    !  status -- return value:  0 for success
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
    PetscInitialized = .TRUE.
    !
  END SUBROUTINE LinearSolverPetscInit
  ! ====================================================   END:  LinearSolverPetscInit
  ! ==================================================== BEGIN:  LinearSolverCreate
  SUBROUTINE LinearSolverCreate(self, conn, bcnodes, locsizes, SYMM)
    !
    !  Create an instance of a linear solver.
    !
    ! ========== Arguments
    !
    TYPE(LinearSolverType), INTENT(IN OUT) :: self
    !
    INTEGER, INTENT(IN) :: conn(:, :)
    INTEGER, INTENT(IN) :: bcnodes(:)
    INTEGER, INTENT(IN) :: locsizes(0:)
    LOGICAL, INTENT(IN), OPTIONAL :: SYMM
    !
    INTEGER, INTENT(OUT), OPTIONAL :: STATUS
    !
    !  locsizes -- number of degrees of freedom on each processor
    !
    !  conn     -- 1-based connectivity section (for DOF's of the system)
    !           -- each process passes some section of the connectivity
    !
    !  bcnodes  -- 1-based global list of DOF numbers with essential BC's
    !           -- each process needs to pass the complete list of DOF numbers.
    !
    !  locsizes -- number of degrees of freedom on each processor
    !
    !  SYMM     -- symmetric system indicator (optional) 
    !
    !  status   -- return value (optional):  0 for success (not used at the moment)
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
    IF (.NOT. PetscInitialized) THEN
      CALL LinearSolverPetscInit()
    END IF
    !
    !  Set status.
    !
    IF (PRESENT(STATUS)) THEN
      STATUS = 0
    END IF
    !
    !  Set symmetry.
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
  END SUBROUTINE LinearSolverCreate
  ! ==================================================== LinearSolverCreate:  END
  !
  SUBROUTINE LinearSolverPetscSolve(sysid, sol, emats, &
       &   ERHS, ARHS, &
       &   BCVALS, SOL_GLOBAL, STATUS)
    !
    !  ***  Description:  
    !
    !  Solve linear system using Petsc routines.
    !
    !  *** Argument Declarations:
    !
    !  sysid -- system identifier
    !
    INTEGER, INTENT(IN) :: sysid
    !
    !  sol -- the local solution vector (section on this process)
    !
    DOUBLE PRECISION, INTENT(INOUT) :: sol(:)
    !
    !  emats -- the local elemental stiffness matrices
    !
    DOUBLE PRECISION, INTENT(IN) :: emats(:, :, :)
    !
    !  erhs -- the local elemental right-hand sides
    !
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: ERHS(:, :)
    !
    !  arhs -- assembled right-hand side (1-based array)
    !
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: ARHS(:)
    !
    !  NOTE:  One of 'erhs' or 'arhs' must be present, but not both.
    !
    !  bcvals -- boundary values to enforce (usually nonzero)
    !
    DOUBLE PRECISION, INTENT(IN), OPTIONAL :: BCVALS(:)
    !
    !  SOL_GLOBAL -- global solution vector
    !
    DOUBLE PRECISION, INTENT(OUT), OPTIONAL :: SOL_GLOBAL(:)
    !
    !  STATUS -- return value; convergence status from PetSc
    !
    INTEGER, INTENT(OUT), OPTIONAL :: STATUS
    !
    !
    !--------------Locals------------------------------------------------------
    !
    INTEGER :: ierr, dofpe, elmin, elmax, numbc
    INTEGER :: elem, node, its, ibc, vecsize, i, iproc
    INTEGER, POINTER :: bcn(:)
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
    !
    !--------------*-------------------------------------------------------
    !
    if (.NOT. PRESENT(ERHS) .AND. .NOT. PRESENT(ARHS)) then 
      ! Neither are present.
      RETURN ! Error status
    end if
    !
    !  Form right hand side according to input options.
    !
    if (PRESENT(arhs)) then
      !
      call BFromARHS(sysid, arhs)
    else 
      !
      call VecSet   (b(sysid), 0.0d0, ierr)
      call VecAssemblyBegin(b(sysid), ierr)
      call VecAssemblyEnd  (b(sysid), ierr)
      !
    end if
    !
    if (PRESENT(BCVALS)) then
      !
      bcn   => sysconn(sysid)%bcs
      numbc = SIZE(bcn, 1)
      !
      if (PRESENT(erhs)) then
        call BFromERHS(sysid, erhs, USE_ALL_DOF=.TRUE.)
      end if
      call FormA(sysid, emats, APPLY_BCS=.FALSE.)
      !
      !  Form x0 satisfying bcs.
      !
      call VecSet(x0(sysid), 0.0d0, ierr)
      !
      IF (myrank == 0) THEN
        call VecSetValues(x0(sysid), numbc, bcn, bcvals, &
             &  INSERT_VALUES, ierr)
      end if
      !
      call VecAssemblyBegin(x0(sysid), ierr)
      call VecAssemblyEnd  (x0(sysid), ierr)
      !
      call MatMult(A(sysid), x0(sysid), ax0(sysid), ierr)
      call VecAXPY(ax0(sysid),  -1.0d0,   b(sysid), ierr)
      !
      !  Prepare to solve homogeneous problem:  A(x-x0) = b-b0
      !
      if (myrank == 0) then
        call VecSetValues(b(sysid), numbc, bcn, &
             &  SPREAD(0.0d0, DIM=1, NCOPIES=numbc), &
             &  INSERT_VALUES, ierr)
      end if
      !
      call VecAssemblyBegin(b(sysid), ierr)
      call VecAssemblyEnd  (b(sysid), ierr)
      !
    else ! homogeneous boundary values
      if (PRESENT(ERHS)) then
        call BFromERHS(sysid, ERHS, USE_ALL_DOF=.FALSE.)
      end if
    end if
    !
    !  Now form matrix and solve system.
    !
    call FormA(sysid, emats, APPLY_BCS=.TRUE.)
    call SolveSystem(sysid, mystatus)
    !
    !  Quit if system fails.
    !
    if (mystatus < 0) then
      call closfl()
    end if
    !
    if (PRESENT(STATUS)) then
      STATUS = mystatus
    end if
    !
    if (PRESENT(BCVALS)) then ! add back in x0
      ! testing
      !was ...! CALL VecAXPY(x0(sysid), 1.0d0, x(sysid), ierr)
      CALL VecAXPY(x(sysid), 1.0d0, x0(sysid), ierr)
    endif
    !
    !  Copy solution to 'sol'
    !
    vecsize = localsizes(myrank, sysid)
    !allocate(xx_v(vecsize))
    call VecGetArrayF90(x(sysid), xx_v, ierr)
    PRINT *, '*** GetArrayF90:  ', ierr

    CALL VecMin(x(sysid), iTmp, vMin, ierr)
    CALL VecMax(x(sysid), iTmp, vMax, ierr)
    PRINT *, '*** max(x):  ', vMin, vMax
    PRINT *, '*** max(xx_v):  ', MAXVAL(xx_v), MINVAL(xx_v)

    sol = xx_v
    call VecRestoreArrayF90(x(sysid), xx_v, ierr)
    PRINT *, '*** RestoreArrayF90:  ', ierr
    !deallocate(xx_v)
    WRITE(99, *) '---sol---'
    WRITE(99, *) sol
    !
    !using VecGetArray!!
    !using VecGetArray!vecsize = localsizes(myrank, sysid)
    !using VecGetArray!call VecGetArray(x(sysid), xptr, offs, ierr)
    !using VecGetArray!do i=1, vecsize
    !using VecGetArray!   sol(i) = xptr(offs + i)
    !using VecGetArray!end do
    !using VecGetArray!call VecRestoreArray(x(sysid), xptr, offs, ierr)
    !
    !  Now construct global solution if requested.
    !
    if (PRESENT(SOL_GLOBAL)) then
      !debug!! 
      !debug!!  Get all local vectors on node 0 and return.
      !debug!!
      !debug!call MPI_AllGatherV(&
      !debug!     &  sol, vecsize, MPI_DOUBLE_PRECISION, &
      !debug!     &  SOL_GLOBAL, localsizes(:, sysid), &
      !debug!     &  section_beg(:, sysid) - 1, MPI_DOUBLE_PRECISION,&
      !debug!     &  PETSC_COMM_WORLD, ierr)
      !debug!PRINT *, '*** AllGather status:  ', ierr
      !debug!WRITE(99, *) '---sol global---'
      !debug!WRITE(99, *) SOL_GLOBAL
      SOL_GLOBAL = sol
      !
    endif
    !
  END SUBROUTINE LinearSolverPetscSolve
  !
  !  *** Program Unit:  private subroutine
  !  ***    Unit Name:  BFromARHS
  !
  !  *** Unit Declaration: 
  !
  SUBROUTINE BFromARHS(sysid, rhs)
    !
    !  ***  Description:  
    !
    !  Compute `b(sysid)' from assembled right-hand side.
    !  
    !
    !  *** Argument Declarations:
    !
    INTEGER, INTENT(IN) :: sysid
    !
    !  sysid -- system identifier
    !
    DOUBLE PRECISION, INTENT(IN) :: rhs(:)
    !
    !  rhs -- assembled right-hand side (1-based array)
    !
    !  *** End:
    !
    !  *** Locals:
    !
    INTEGER :: i, ierr
    !
    !--------------*-------------------------------------------------------
    !
    do i=1, localsizes(myrank, sysid)
      call VecSetValues(b(sysid), &
           &  1, section_beg(myrank, sysid) + i - 2, rhs(i), &
           &  INSERT_VALUES, ierr)
    end do
    !
    call VecAssemblyBegin(b(sysid), ierr)
    call VecAssemblyEnd  (b(sysid), ierr)
    !
  END SUBROUTINE BFromARHS
  !
  !  *** Program Unit:  private subroutine
  !  ***    Unit Name:  BFromERHS
  !
  !  *** Unit Declaration: 
  !
  SUBROUTINE BFromERHS(sysid, erhs, USE_ALL_DOF)
    !
    !  ***  Description:  
    !
    !  Compute `b(sysid)' from elemental right-hand sides.
    !  This appends to existing values.
    !
    !  *** Argument Declarations:
    !
    INTEGER, INTENT(IN) :: sysid
    !
    !  sysid -- system identifier
    !
    DOUBLE PRECISION, INTENT(IN) :: erhs(:, :)
    !
    !  erhs -- elemental right-hand side (references a 1-based array)
    !
    LOGICAL, INTENT(IN) :: USE_ALL_DOF
    !
    !  USE_ALL_DOF -- flag to assemble all DOF regardless of BC
    !
    !  *** End:
    !
    !  *** Locals:
    !
    INTEGER :: i, ierr, dofpe, elem, elmin, elmax, numbc
    INTEGER, POINTER :: sconn(:, :), bcn(:)
    !
    !--------------*-------------------------------------------------------
    !
    dofpe = SIZE  (erhs, 1)
    elmin = LBOUND(erhs, 2)
    elmax = UBOUND(erhs, 2)
    !
    sconn => sysconn(sysid)%sconn
    bcn   => sysconn(sysid)%bcs
    numbc = SIZE(bcn, 1)
    !
    !  NOTE:  2009-02-20
    !
    !  * Fixed second VecSetValues call, which had negative entries in sconn
    !  * Use "abs(1 + sconn) - 1", since 1 is subtracted from connectivity,
    !    even from negative values (BCs)
    !
    if (USE_ALL_DOF) then
      do elem=elmin, elmax
        call VecSetValues(b(sysid), dofpe, abs(1 + sconn(:, elem)) - 1, &
             &  erhs(:, elem), &
             &  ADD_VALUES, ierr)
      end do
      !
      call VecAssemblyBegin(b(sysid), ierr)
      call VecAssemblyEnd  (b(sysid), ierr)
      !
    else
      do elem=elmin, elmax
        call VecSetValues(b(sysid), dofpe, abs(1 + sconn(:, elem)) - 1, &
             &  erhs(:, elem), &
             &  ADD_VALUES, ierr)
        !
      end do
      !
      call VecAssemblyBegin(b(sysid), ierr)
      call VecAssemblyEnd  (b(sysid), ierr)
      !
      if (myrank == 0) then
        CALL VecSetValues(b(sysid), numbc, bcn, &
             &  SPREAD(0.0d0, DIM=1, NCOPIES=numbc), &
             &  INSERT_VALUES, ierr)
      end if
      !
      call VecAssemblyBegin(b(sysid), ierr)
      call VecAssemblyEnd  (b(sysid), ierr)
      !
    end if
    !
  END SUBROUTINE BFromERHS
  !
  !  *** Program Unit:  private subroutine
  !  ***    Unit Name:  FormA
  !
  !  *** Unit Declaration: 
  !
  SUBROUTINE FormA(sysid, emats, APPLY_BCS)
    !
    !  ***  Description:  
    !
    !  Form A(sysid) from the elemental matrices.
    !
    !  *** Argument Declarations:
    !
    INTEGER, INTENT(IN) :: sysid
    !
    !  sysid -- system identifier
    !
    DOUBLE PRECISION, INTENT(IN) :: emats(:, :, :)
    !
    !  emats -- the local elemental stiffness matrices
    !
    LOGICAL, INTENT(IN) :: APPLY_BCS
    !
    !  APPLY_BCS -- flag for applying boundar conditions
    !
    !  If true, set up the matrix for applied boundary conditions;
    !  otherwise, assemble the matrix for all degrees of freedom.
    !
    !  *** End:
    !
    !  *** Locals:
    !
    INTEGER :: dofpe, elem, elmin, elmax, numbc, ibc, ierr
    INTEGER, POINTER :: sconn(:, :), bcn(:)
    !
    !--------------*-------------------------------------------------------
    !
    dofpe = SIZE(emats, 1)
    elmin = LBOUND(emats, 3)
    elmax = UBOUND(emats, 3)
    !
    sconn => sysconn(sysid)%sconn
    bcn   => sysconn(sysid)%bcs
    numbc = SIZE(bcn, 1)

    call MatZeroEntries(A(sysid), ierr)

    if (APPLY_BCS) then
      !
      !  Negative DOF number is ignored by Petsc.
      !
      do elem=elmin, elmax
        call MatSetValues(A(sysid), &
             &   dofpe, sconn(:, elem),&
             &   dofpe, sconn(:, elem),&
             &   TRANSPOSE(emats(:, :, elem)), ADD_VALUES, ierr)
      end do
      !
      if (myrank == 0) then
        do ibc=1, numbc
          call MatSetValues(A(sysid), &
               &  1, bcn(ibc), 1, bcn(ibc), &
               &  1.0d0, ADD_VALUES, ierr)
        end do
      end if
    else
      !
      !  Use absolute value of connectivity to ensure
      !  all DOF are assembled.
      !
      do elem=elmin, elmax
        call MatSetValues(A(sysid), &
             &   dofpe, ABS(1 + sconn(:, elem)) - 1,&
             &   dofpe, ABS(1 + sconn(:, elem)) - 1,&
             &   TRANSPOSE(emats(:, :, elem)), ADD_VALUES, ierr)
      end do
    endif
    !
    call MatAssemblyBegin(A(sysid),MAT_FINAL_ASSEMBLY,ierr)
    call MatAssemblyEnd  (A(sysid),MAT_FINAL_ASSEMBLY,ierr)
    !
  END SUBROUTINE FormA
  !
  !  *** Program Unit:  private subroutine
  !  ***    Unit Name:  SolveSystem
  !
  !  *** Unit Declaration: 
  !
  SUBROUTINE SolveSystem(sysid, STATUS)
    !
    !  ***  Description:  
    !
    !  Solve the system A(sysid)*x(sysid) = b(sysid)
    !
    !  A and b are already constructed.  Currently this
    !  only used Jacobi preconditioner.
    !
    !  *** Argument Declarations:
    !
    !  sysid -- system identifier
    !
    INTEGER, INTENT(IN) :: sysid
    !
    !  STATUS -- return value 
    !
    INTEGER, INTENT(OUT), OPTIONAL :: STATUS
    !
    !  *** End:
    !
    !  *** Locals:
    !
    INTEGER :: ierr, its, maxits
    !
    !--------------*-------------------------------------------------------
    !
    !--------------* Form Petsc solver context.
    !
    CALL KSPCreate(PETSC_COMM_WORLD, ksp, ierr)
    !
    call KSPSetOperators(ksp, A(sysid), A(sysid), &
         &  DIFFERENT_NONZERO_PATTERN, ierr)
    !
    !  Set solver and preconditioner options.
    !
    if (SYMMETRY(sysid)) then
      ksptype = KSPCG
    else
      ksptype = KSPBICG
    end if
    call KSPSetType(ksp, ksptype, ierr)
    !
    ptype = PCJACOBI
    CALL KSPGetPC(ksp, pc, ierr)
    call PCSetType(pc, ptype, ierr)
    !
    maxits = 2*(section_beg(numprocs - 1, sysid) + localsizes(numprocs - 1, sysid))
    call KSPSetTolerances(ksp, error_tol, &
         &   PETSC_DEFAULT_DOUBLE_PRECISION,&
         &   PETSC_DEFAULT_DOUBLE_PRECISION, &
         &   maxits, ierr)

    CALL KSPSetFromOptions(ksp, ierr)
    CALL KSPSetUp(ksp, ierr)
    call KSPSolve(ksp, b(sysid), x(sysid), ierr) ! solve!

    call KSPGetConvergedReason(ksp, reason, ierr)
    CALL KSPGetIterationNumber(ksp, its, ierr)
    PRINT *, '   PETSc results:  sysid/its/reason: ', sysid, its, reason
    !
    if (PRESENT(STATUS)) then
      STATUS = reason
      IF (reason < 0) THEN
        CALL DescribeDivergedReason(6, reason)
      END IF
    end if
    !
    call KSPDestroy(ksp,ierr)

    !
  END SUBROUTINE SolveSystem
  !
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
