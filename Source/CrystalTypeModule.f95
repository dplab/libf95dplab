MODULE CrystalTypeModule 
!
!  *** Program Unit:  Module
!  ***    Unit Name:  CrystalTypeModule
!  ***  Description:
!
!  This module handles the basic CrystalType object.
!
!  TODO:  
!
!  *  Add vertices to object structure
!  *  Add P*P^T to "get list"
!
!  *** Use Statements:
!
USE IntrinsicTypesModule, &
     &  RK=>REAL_KIND, IK=>INTEGER_KIND
USE ConstantsModule
USE FilesModule
USE Tensor3DModule 
!
!  *** End:
!
IMPLICIT NONE
!
PRIVATE
!
!  Crystal class indicators.
!
INTEGER(IK), PARAMETER :: CLASS_FCC=1, CLASS_BCC=2, CLASS_HCP=3
INTEGER(IK) :: DECOMP_DFLT = DECOMP_MPSIM
!
!  *** Derived Types:
!
TYPE CrystalTypeType
   !
   !PRIVATE  ! allowing access to components
   !
   !CHARACTER, POINTER :: name(:)
   !
   INTEGER  :: class   ! FCC, BCC, HCP
   INTEGER  :: decomp  ! Decomposition convention
   INTEGER  :: numslip, numvertices
   !
   REAL(RK), POINTER :: schmid_3x3(:, :, :), dev(:, :), skw(:, :), pptrans(:, :, :)
   REAL(RK), POINTER :: vertices(:, :), vertices3x3(:, :, :)
   !
END TYPE CrystalTypeType
!
!--------------Public Entities
!
!
!  *** Public Types:
!
PUBLIC :: CrystalTypeType 
!
!  *** Public Data:
!
!  Below are flags indicating crystal type and decomposition convention.
!
PUBLIC :: CLASS_FCC, CLASS_BCC, CLASS_HCP
PUBLIC :: DECOMP_MPSIM, DECOMP_FEMEVPS
!
!  *** Public Procedures:
!
PUBLIC :: CrystalTypeCreate, CrystalTypeDescribe,&
     &  CrystalTypeGet
!
!  *** End:
!
!--------------*-------------------------------------------------------
!
CONTAINS 
!
!  *** Program Unit:  function
!  ***    Unit Name:  CrystalTypeCreate
!
!  *** Unit Declaration: 
!
FUNCTION CrystalTypeCreate(CTYPE, &
     &   C_OVER_A, HRATIO_HCP, VERTEX_FILE, DECOMP) &
     &   RESULT(self)
  !
  !  ***  Description:  
  !
  !  Create a CrystalTypeType object.
  !
  !  *** Argument Declarations:
  !
  !  CTYPE -- crystal type
  !
  INTEGER, INTENT(IN) :: CTYPE
  !
  !  C_OVER_A -- C over A ratio
  !
  REAL(RK), INTENT(IN), OPTIONAL :: C_OVER_A
  !
  !  HRATIO_HCP -- Ratio of pyramidal strengths to basal/prismatic
  !
  REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_HCP
  !
  !  VERTEX_FILE -- file containing list of vertices
  !
  CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: VERTEX_FILE
  !
  !  DECOMP -- decomposition convention indicator
  !
  INTEGER, INTENT(IN), OPTIONAL :: DECOMP 
  !
  !  *** Result:
  !  
  TYPE(CrystalTypeType), POINTER :: self
  !
  !  *** Locals:
  !
  INTEGER(IK) :: mystat, unit, numvert
  !
  !  *** End:
  !
  !--------------*-------------------------------------------------------
  !
  ALLOCATE(self, STAT=mystat)
  !
  self%class = CTYPE
  !
  !
  !  Expression is a scalar of type integer, character, or logical.
  !
  SELECT CASE(CTYPE)
     !
  CASE (CLASS_FCC, CLASS_BCC)
     !
     call SchmidTensors(CTYPE, self%schmid_3x3)
     !
  CASE (CLASS_HCP)
     !
     !  check whether arguments are present ...
     !
     call SchmidTensors(CTYPE, self%schmid_3x3,&
          &  C_OVER_A, HRATIO_HCP)
     !
  CASE DEFAULT 
     !
  END SELECT
  !
  self%numslip = SIZE(self%schmid_3x3, DIM=3)
  !
  if (PRESENT(DECOMP)) then
     self%decomp = DECOMP
  else
     self%decomp = DECOMP_DFLT
  end if
  !
  call CrystalTypeGet(self, &
       &   DEV=self%dev, SKW=self%skw, PPTRANS=self%pptrans)
  !
  if (PRESENT(VERTEX_FILE)) then
     unit = GetUnitNumber(VERTEX_FILE)
     read(unit, *) numvert
     ALLOCATE(self%vertices3x3(3, 3, numvert))
     read(unit, *) self%vertices3x3
     close(unit)
     !
     self%numvertices = numvert
     !
     ALLOCATE(self%vertices(5, self%numvertices))
     call Tensor3DDecompose(self%vertices3x3, DEV=self%vertices, DECOMP=self%decomp)
     !
  else
     self%numvertices = 0
  end if
  !
END FUNCTION CrystalTypeCreate
!
!  *** Program Unit:  subroutine
!  ***    Unit Name:  SchmidTensors
!
!  *** Unit Declaration: 
!
SUBROUTINE SchmidTensors(CTYPE, SCHMID, &
     &   C_OVER_A, HRATIO_HCP)
  !
  !  ***  Description:  
  !
  !  *** Arguments:
  !
  !  CTYPE -- crystal type
  !
  INTEGER, INTENT(IN) :: CTYPE
  !
  !  SCHMID -- Schmid tensors
  !
  REAL(RK), POINTER :: SCHMID(:, :, :)
  !
  !  C_OVER_A -- C over A ratio
  !
  REAL(RK), INTENT(IN), OPTIONAL :: C_OVER_A
  !
  !  HRATIO_HCP -- Ratio of pyramidal strengths to basal/prismatic
  !
  REAL(RK), INTENT(IN), OPTIONAL :: HRATIO_HCP
  !
  !  *** Locals:
  !
  !  Various factors making unit vectors.
  !
  REAL(RK), PARAMETER :: Z  =  RK_ZERO
  REAL(RK), PARAMETER :: P2 =  RK_ONE/RK_ROOT_2, P3 = RK_ONE/RK_ROOT_3
  REAL(RK), PARAMETER :: M2 = -P2, M3 = -P3
  !
  !  For BCC {211} plane normals.
  !
  REAL(RK), PARAMETER :: P6_1 = RK_ONE/RK_ROOT_6
  REAL(RK), PARAMETER :: P6_2 = RK_TWO * P6_1
  REAL(RK), PARAMETER :: M6_1 = -P6_1, M6_2 = -P6_2
  !
  !  For BCC {123} plane normals.
  !
  REAL(RK), PARAMETER :: P14_1 =   RK_ONE/RK_ROOT_14
  REAL(RK), PARAMETER :: P14_2 =   RK_TWO/RK_ROOT_14
  REAL(RK), PARAMETER :: P14_3 = RK_THREE/RK_ROOT_14
  REAL(RK), PARAMETER :: M14_1=-P14_1, M14_2=-P14_2, M14_3=-P14_3
  !     
  !  SLIP NORMAL AND DIRECTIONS. 
  ! 
  REAL(RK), PARAMETER, DIMENSION(36) :: cub_111_dat = (/&
       &   P3, P3, P3,     P3, P3, P3,    P3, P3, P3,&
       &   P3, P3, M3,     P3, P3, M3,    P3, P3, M3,&
       &   P3, M3, P3,     P3, M3, P3,    P3, M3, P3,&
       &   P3, M3, M3,     P3, M3, M3,    P3, M3, M3 &
       &  /)
  REAL(RK), PARAMETER, DIMENSION(3, 12) :: &
       &   cub_111 = RESHAPE(SOURCE=cub_111_dat, SHAPE=(/3, 12/))
  !     
  REAL(RK), PARAMETER, DIMENSION(36) :: cub_110_dat = (/&
       &  Z, P2, M2,    P2, Z, M2,    P2, M2, Z,&
       &  Z, P2, P2,    P2, Z, P2,    P2, M2, Z,&
       &  Z, P2, P2,    P2, Z, M2,    P2, P2, Z,&
       &  Z, P2, M2,    P2, Z, P2,    P2, P2, Z &
       &  /)
  REAL(RK), PARAMETER, DIMENSION(3, 12) :: &
       &   cub_110 = RESHAPE(SOURCE=cub_110_dat, SHAPE=(/3, 12/))
  !     
  REAL(RK), PARAMETER, DIMENSION(36) :: cub_211_dat = (/&
       &  M6_1, M6_1, P6_2,    M6_1, P6_2, M6_1,    P6_2, M6_1, M6_1,&
       &  M6_1, M6_1, M6_2,    M6_1, P6_2, P6_1,    P6_2, M6_1, P6_1,&
       &  M6_1, P6_1, P6_2,    M6_1, M6_2, M6_1,    P6_2, P6_1, M6_1,&
       &  M6_1, P6_1, M6_2,    M6_1, M6_2, P6_1,    P6_2, P6_1, P6_1 &
       &/)
  REAL(RK), PARAMETER, DIMENSION(3, 12) :: &
       &   cub_211 = RESHAPE(SOURCE=cub_211_dat, SHAPE=(/3, 12/))
  !     
  REAL(RK), PARAMETER, DIMENSION(36) :: cub_123_a_dat = (/&
       &  M14_1, M14_2, P14_3, M14_1, P14_3, M14_2,  P14_3, M14_1, M14_2,&
       &  M14_1, M14_2, M14_3, M14_1, P14_3, P14_2,  P14_3, M14_1, P14_2,&
       &  M14_1, P14_2, P14_3, M14_1, M14_3, M14_2,  P14_3, P14_1, M14_2,&
       &  M14_1, P14_2, M14_3, M14_1, M14_3, P14_2,  P14_3, P14_1, P14_2 &
       &  /)
  REAL(RK), PARAMETER, DIMENSION(3, 12) :: &
       &   cub_123_a = RESHAPE(SOURCE=cub_123_a_dat, SHAPE=(/3, 12/))
  !
  REAL(RK), PARAMETER, DIMENSION(36) :: cub_123_b_dat = (/&
       &   M14_2, M14_1, P14_3, M14_2, P14_3, M14_1,  P14_3, M14_2, M14_1,& 
       &   M14_2, M14_1, M14_3, M14_2, P14_3, P14_1,  P14_3, M14_2, P14_1,& 
       &   M14_2, P14_1, P14_3, M14_2, M14_3, M14_1,  P14_3, P14_2, M14_1,& 
       &   M14_2, P14_1, M14_3, M14_2, M14_3, P14_1,  P14_3, P14_2, P14_1 &
       &  /)
  REAL(RK), PARAMETER, DIMENSION(3, 12) :: &
       &   cub_123_b = RESHAPE(SOURCE=cub_123_b_dat, SHAPE=(/3, 12/))
  !
  !  HCP slip data.
  !
  !  ... parameters
  !
  REAL(RK), PARAMETER :: ONE = 1.0_RK, HALF = 0.5_RK
  REAL(RK), PARAMETER :: COS_30 = HALF*RK_ROOT_3
  !
  REAL(RK), PARAMETER, DIMENSION(3, 3) :: &
       &   hex_basal_sn = RESHAPE(&
       &      SOURCE=(/&
       &         Z, Z, ONE,     Z, Z, ONE,     Z, Z, ONE &
       &      /), &
       &      SHAPE=(/3, 3/))
  !
  REAL(RK), PARAMETER, DIMENSION(3, 3) :: &
       &   hex_basal_sd = RESHAPE(&
       &      SOURCE=(/&
       &         ONE, Z, Z,    -HALF, COS_30, Z,   -HALF,  -COS_30, Z &
       &      /), &
       &      SHAPE=(/3, 3/))
  !
  REAL(RK), PARAMETER, DIMENSION(3, 3) :: &
       &   hex_pris_sn = RESHAPE(&
       &      SOURCE=(/&
       &         Z, ONE, Z,     -COS_30, -HALF,  Z,     COS_30, -HALF,  Z &
       &      /), &
       &      SHAPE=(/3, 3/))
  !
  REAL(RK), PARAMETER, DIMENSION(3, 3) :: &
       &   hex_pris_sd = RESHAPE(&
       &      SOURCE=(/&
       &         ONE, Z, Z,    -HALF, COS_30, Z,       -HALF, -COS_30, Z &
       &      /), &
       &      SHAPE=(/3, 3/))
  !
  REAL(RK), PARAMETER, DIMENSION(4, 12) :: &
       &   hex_pyr1_sn = RESHAPE(&
       &      SOURCE=(/&
       &          1.0,  0.0, -1.0,  1.0,&
       &          1.0,  0.0, -1.0,  1.0,&
       &          0.0,  1.0, -1.0,  1.0,&
       &          0.0,  1.0, -1.0,  1.0,&
       &         -1.0,  1.0,  0.0,  1.0,&
       &         -1.0,  1.0,  0.0,  1.0,&
       &         -1.0,  0.0,  1.0,  1.0,&
       &         -1.0,  0.0,  1.0,  1.0,&
       &          0.0, -1.0,  1.0,  1.0,&
       &          0.0, -1.0,  1.0,  1.0,&
       &          1.0, -1.0,  0.0,  1.0,&
       &          1.0, -1.0,  0.0,  1.0 &
       &      /), &
       &      SHAPE=(/4,  12/))
  !
  !  NOTE:  does 3.0 need to be 3.0_RK ?
  !
  REAL(RK), PARAMETER, DIMENSION(4, 12) :: &
       &   hex_pyr1_sd = RESHAPE(&
       &      SOURCE=(/&
       &         -2.0,  1.0,  1.0,  3.0,&
       &         -1.0, -1.0,  2.0,  3.0,&
       &         -1.0, -1.0,  2.0,  3.0,&
       &          1.0, -1.0,  2.0,  3.0,& 
       &          1.0, -2.0,  1.0,  3.0,& 
       &          2.0,  1.0, -1.0,  3.0,& 
       &          2.0, -1.0, -1.0,  3.0,& 
       &          1.0,  1.0, -2.0,  3.0,& 
       &          1.0,  1.0, -2.0,  3.0,& 
       &         -1.0,  2.0, -1.0,  3.0,& 
       &         -1.0,  2.0, -1.0,  3.0,& 
       &         -2.0,  1.0,  1.0,  3.0 &
       &       /), &
       &      SHAPE=(/4,  12/))
  !
  !  Schmid tensors.
  !
  REAL(RK) :: snrm(3, 12), sdir(3, 12), rescale
  !
  !  *** End:
  !
  !--------------*-------------------------------------------------------
  !
  SELECT CASE(CTYPE)
     !
  CASE (CLASS_FCC)
     !
     ALLOCATE(SCHMID(3, 3, 12))
     !
     schmid(1, 1, :) = cub_110(1, :)*cub_111(1, :)
     schmid(2, 1, :) = cub_110(2, :)*cub_111(1, :)
     schmid(3, 1, :) = cub_110(3, :)*cub_111(1, :)

     schmid(1, 2, :) = cub_110(1, :)*cub_111(2, :)
     schmid(2, 2, :) = cub_110(2, :)*cub_111(2, :)
     schmid(3, 2, :) = cub_110(3, :)*cub_111(2, :)

     schmid(1, 3, :) = cub_110(1, :)*cub_111(3, :)
     schmid(2, 3, :) = cub_110(2, :)*cub_111(3, :)
     schmid(3, 3, :) = cub_110(3, :)*cub_111(3, :)
     !
  CASE (CLASS_BCC)
     !
     ALLOCATE(SCHMID(3, 3, 12))
     !
     schmid(1, 1, :) = cub_111(1, :)*cub_110(1, :)
     schmid(2, 1, :) = cub_111(2, :)*cub_110(1, :)
     schmid(3, 1, :) = cub_111(3, :)*cub_110(1, :)

     schmid(1, 2, :) = cub_111(1, :)*cub_110(2, :)
     schmid(2, 2, :) = cub_111(2, :)*cub_110(2, :)
     schmid(3, 2, :) = cub_111(3, :)*cub_110(2, :)

     schmid(1, 3, :) = cub_111(1, :)*cub_110(3, :)
     schmid(2, 3, :) = cub_111(2, :)*cub_110(3, :)
     schmid(3, 3, :) = cub_111(3, :)*cub_110(3, :)
     !
  CASE (CLASS_HCP)
     !
     ALLOCATE(SCHMID(3, 3, 18))
     !
     !  Basal.
     !
     schmid(1, 1, 1:3) = hex_basal_sd(1, :)*hex_basal_sn(1, :)
     schmid(2, 1, 1:3) = hex_basal_sd(2, :)*hex_basal_sn(1, :)
     schmid(3, 1, 1:3) = hex_basal_sd(3, :)*hex_basal_sn(1, :)

     schmid(1, 2, 1:3) = hex_basal_sd(1, :)*hex_basal_sn(2, :)
     schmid(2, 2, 1:3) = hex_basal_sd(2, :)*hex_basal_sn(2, :)
     schmid(3, 2, 1:3) = hex_basal_sd(3, :)*hex_basal_sn(2, :)

     schmid(1, 3, 1:3) = hex_basal_sd(1, :)*hex_basal_sn(3, :)
     schmid(2, 3, 1:3) = hex_basal_sd(2, :)*hex_basal_sn(3, :)
     schmid(3, 3, 1:3) = hex_basal_sd(3, :)*hex_basal_sn(3, :)
     !
     !  Prismatic.
     !
     schmid(1, 1, 4:6) = hex_pris_sd(1, :)*hex_pris_sn(1, :)
     schmid(2, 1, 4:6) = hex_pris_sd(2, :)*hex_pris_sn(1, :)
     schmid(3, 1, 4:6) = hex_pris_sd(3, :)*hex_pris_sn(1, :)

     schmid(1, 2, 4:6) = hex_pris_sd(1, :)*hex_pris_sn(2, :)
     schmid(2, 2, 4:6) = hex_pris_sd(2, :)*hex_pris_sn(2, :)
     schmid(3, 2, 4:6) = hex_pris_sd(3, :)*hex_pris_sn(2, :)

     schmid(1, 3, 4:6) = hex_pris_sd(1, :)*hex_pris_sn(3, :)
     schmid(2, 3, 4:6) = hex_pris_sd(2, :)*hex_pris_sn(3, :)
     schmid(3, 3, 4:6) = hex_pris_sd(3, :)*hex_pris_sn(3, :)
     !       
     !       Pyramidal 1.
     !
     rescale  = 1.0d0/HRATIO_HCP
     !
     !         Convert Miller indices to spatial directions.
     !
     snrm(1, :) = hex_pyr1_sn(1, :)
     snrm(2, :) = (2.0_RK*hex_pyr1_sn(2, :) + snrm(1, :))/RK_ROOT_3
     snrm(3, :) = hex_pyr1_sn(4, :)/C_OVER_A

     snrm = snrm / SPREAD(SOURCE=SQRT(SUM(snrm*snrm, DIM=1)),&
          &   DIM=1, NCOPIES=3)

     sdir(1, :) = 1.5_RK*hex_pyr1_sd(1, :)
     sdir(2, :) =(hex_pyr1_sd(2, :) + 0.5_RK*hex_pyr1_sd(1, :))*RK_ROOT_3
     sdir(3, :) = hex_pyr1_sd(4, :)*C_OVER_A

     sdir = sdir / SPREAD(SOURCE=SQRT(SUM(sdir*sdir, DIM=1)),&
          &   DIM=1, NCOPIES=3)


     !         
     schmid(1, 1, 7:18) = rescale*sdir(1, :)*snrm(1, :)
     schmid(2, 1, 7:18) = rescale*sdir(2, :)*snrm(1, :)
     schmid(3, 1, 7:18) = rescale*sdir(3, :)*snrm(1, :)

     schmid(1, 2, 7:18) = rescale*sdir(1, :)*snrm(2, :)
     schmid(2, 2, 7:18) = rescale*sdir(2, :)*snrm(2, :)
     schmid(3, 2, 7:18) = rescale*sdir(3, :)*snrm(2, :)

     schmid(1, 3, 7:18) = rescale*sdir(1, :)*snrm(3, :)
     schmid(2, 3, 7:18) = rescale*sdir(2, :)*snrm(3, :)
     schmid(3, 3, 7:18) = rescale*sdir(3, :)*snrm(3, :)
     !
  CASE DEFAULT 
     !
  END SELECT
  !
END SUBROUTINE SchmidTensors
!
!  *** Program Unit:  subroutine
!  ***    Unit Name:  CrystalTypeDescribe
!
!  *** Unit Declaration: 
!
SUBROUTINE CrystalTypeDescribe(self)
  !
  !  ***  Description:  
  !
  !  Print information about crystal type.
  !
  !  *** Arguments:
  !
  !  self -- this crystal type
  !
  TYPE(CrystalTypeType) :: self
  !
  !  *** Locals:
  !

  !  *** End:
  !
  !--------------*-------------------------------------------------------
  !
  write(6, *) 'class:  ', self%class
  write(6, *) 'number of slip systems:  ', self%numslip
  !
END SUBROUTINE CrystalTypeDescribe
!
!  *** Program Unit:  subroutine
!  ***    Unit Name:  CrystalTypeGet
!
!  *** Unit Declaration: 
!
SUBROUTINE CrystalTypeGet(self, DEV, SKW, PPTRANS, VERTICES,&
     &   NUMSLIP, NUMVERTICES)
  !
  !  ***  Description:  
  !  
  !  Return deviatoric parts of Schmid tensors.
  !
  !  *** Arguments:
  !
  !  self -- the CrystalType object
  !
  TYPE(CrystalTypeType) :: self
  !
  !  DEV -- deviatoric part of Schmid tensors
  !
  REAL(RK), POINTER, OPTIONAL  :: DEV(:, :)
  !
  !  SKW  -- skew part of Schmid tensors
  !
  REAL(RK), POINTER, OPTIONAL  :: SKW(:, :)
  !
  !  PPTRANS -- matrices of Schmid diads
  !
  REAL(RK), POINTER, OPTIONAL  :: PPTRANS(:, :, :)
  !
  !  VERTICES -- vertices in 5-vector form
  !
  REAL(RK), POINTER, OPTIONAL  :: VERTICES(:, :)
  !
  !  NUMSLIP -- number of slip systems
  !
  INTEGER, OPTIONAL :: NUMSLIP 
  !
  !  NUMVERTICES -- number of vertices
  !
  INTEGER, OPTIONAL :: NUMVERTICES
  !
  !  *** Locals:
  !
  INTEGER :: DSHAPE(2), SSHAPE(2), ARGPPTSHAPE(3), MYPPTSHAPE(3), i
  !
  REAL(RK), POINTER :: mydev(:, :)
  !
  !  *** End:
  !
  !--------------*-------------------------------------------------------
  !
  !
  if (PRESENT(DEV)) then
     !
     SSHAPE = (/5, self%numslip/)
     !SHAPE(self%schmid_sym)
     !
     if (ASSOCIATED(DEV)) then
        !
        !  Could check dimensions here ...
        !
        DSHAPE = SHAPE(DEV)
        if ( (DSHAPE(1) /= SSHAPE(1)) .OR. (DSHAPE(2) /= SSHAPE(2))) then
           DEALLOCATE(DEV)
           ALLOCATE(DEV(SSHAPE(1), SSHAPE(2)))
        end if
     else
        ALLOCATE(DEV(SSHAPE(1), SSHAPE(2)))
     end if
     !
     call Tensor3DDecompose(self%schmid_3x3, &
          &  DEV=DEV, DECOMP=self%decomp)
     !
  end if
  !
  if (PRESENT(SKW)) then
     !
     SSHAPE = (/3, self%numslip/)
     !
     if (ASSOCIATED(SKW)) then
        !
        !  Could check dimensions here ...
        !
        DSHAPE = SHAPE(SKW)
        if ( (DSHAPE(1) /= SSHAPE(1)) .OR. (DSHAPE(2) /= SSHAPE(2))) then
           DEALLOCATE(SKW)
           ALLOCATE(SKW(SSHAPE(1), SSHAPE(2)))
        end if
     else
        ALLOCATE(SKW(SSHAPE(1), SSHAPE(2)))
     end if
     !
     call Tensor3DDecompose(self%schmid_3x3, &
          &  SKW=SKW, DECOMP=self%decomp)
     !
  end if
  !
  if (PRESENT(PPTRANS)) then
     !
     MYPPTSHAPE = (/5, 5, self%numslip/)
     !
     if (ASSOCIATED(PPTRANS)) then
        !
        !  Could check dimensions here ...
        !
        ARGPPTSHAPE = SHAPE(PPTRANS)
        if (   (ARGPPTSHAPE(1) /= MYPPTSHAPE(1)) .OR. &
             & (ARGPPTSHAPE(2) /= MYPPTSHAPE(2)) .OR. &
             & (ARGPPTSHAPE(3) /= MYPPTSHAPE(3))  ) then
           DEALLOCATE(PPTRANS)
           ALLOCATE(PPTRANS(MYPPTSHAPE(1), MYPPTSHAPE(2), MYPPTSHAPE(3)))
        end if
     else
        ALLOCATE(PPTRANS(MYPPTSHAPE(1), MYPPTSHAPE(2), MYPPTSHAPE(3)))
     end if
     !
     ALLOCATE(mydev(5, self%numslip))
     !
     call Tensor3DDecompose(self%schmid_3x3, DEV=mydev, DECOMP=self%decomp)
     !
     do i=1, self%numslip
        PPTRANS(:, :, i) = MATMUL(&
             &  RESHAPE(mydev(:, i), SHAPE=(/5, 1/)), &
             &  RESHAPE(mydev(:, i), SHAPE=(/1, 5/)) )
     end do
     !
     DEALLOCATE(mydev)
     !
  end if
  !
  if (PRESENT(VERTICES)) then
     !
     SSHAPE = (/5, self%numvertices/)
     !
     if (ASSOCIATED(VERTICES)) then
        !
        !  Could check dimensions here ...
        !
        DSHAPE = SHAPE(VERTICES)
        if ( (DSHAPE(1) /= SSHAPE(1)) .OR. (DSHAPE(2) /= SSHAPE(2))) then
           DEALLOCATE(VERTICES)
           ALLOCATE(VERTICES(SSHAPE(1), SSHAPE(2)))
        end if
     else
        ALLOCATE(VERTICES(SSHAPE(1), SSHAPE(2)))
     end if
     !
     VERTICES = self%vertices
     !
  end if
  !
  if (PRESENT(NUMSLIP)) then
     NUMSLIP = self%numslip
  end if
  !
  if (PRESENT(NUMVERTICES)) then
     NUMVERTICES = self%numvertices
  end if
  !
  !  Need call to get vertices as 3x3.
  !
END SUBROUTINE CrystalTypeGet

END MODULE CrystalTypeModule 
