!
! Copyright (C) 2017 National Institute of Advanced Industrial Science and Technology (AIST)
! [ This code is written by Satomichi Nishihara. ]
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE eqn_lauedipole(rismt, expand, prepare, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... solve dipole part of Laue-RISM equation, which is defined as
  ! ...
  ! ...             /
  ! ...   hd1(z1) = | dz2 cd2(z2) * x21(gxy=0,z2-z1) .
  ! ...             /
  ! ...
  ! ... if the right-hand side solvent, cd2 is
  ! ...
  ! ...              D2             pi   z2
  ! ...   cd2(z2) = ---- * (1 - sin(-- * --)) .
  ! ...              2              2    z0
  ! ...
  ! ... if the left-hand side solvent, cd2 is
  ! ...
  ! ...              D2             pi   z2
  ! ...   cd2(z2) = ---- * (1 + sin(-- * --)) .
  ! ...              2              2    z0
  ! ...
  ! ... when prepare = .true., only intermediate integrations are performed as
  ! ...
  ! ...    1    /               pi   z2
  ! ...   --- * | dz2 (1 -+ sin(-- * --)) * x21(gxy=0,z2-z1)
  ! ...    2    /               2    z0
  ! ...
  ! ... , which are stored to `hdz'.
  ! ...
  ! ... when prepare = .false., integrations are concluded as
  ! ...
  ! ...              D2    /               pi   z2
  ! ...   hd1(z1) = ---- * | dz2 (1 -+ sin(-- * --)) * x21(gxy=0,z2-z1)
  ! ...              2     /               2    z0
  ! ...
  ! ... , which are added to `hgz' or `hsgz'.
  ! ...
  !
  USE cell_base, ONLY : alat, at
  USE constants, ONLY : pi, K_BOLTZMANN_RY
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE rism,      ONLY : rism_type, ITYPE_LAUERISM
  USE solvmol,   ONLY : solVs, get_nuniq_in_solVs, &
                      & iuniq_to_isite, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  LOGICAL,         INTENT(IN)    :: expand
  LOGICAL,         INTENT(IN)    :: prepare
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER               :: nq
  INTEGER               :: iq1, iq2
  INTEGER               :: iiq1, iiq2
  INTEGER               :: iv2
  INTEGER               :: isolV2
  INTEGER               :: iatom2
  INTEGER               :: iz1, iz2
  INTEGER               :: iiz1, iiz2
  INTEGER               :: izdelt
  INTEGER               :: izsta1, izsta2
  INTEGER               :: izend1, izend2
  INTEGER               :: nzint1, nzint2
  INTEGER               :: izint1, izint2
  INTEGER               :: izsolv
  REAL(DP)              :: beta
  REAL(DP)              :: qv2
  REAL(DP)              :: z
  REAL(DP)              :: z0
  REAL(DP)              :: zstep
  REAL(DP)              :: zoffs
  REAL(DP)              :: zedge
  REAL(DP)              :: ssign
  REAL(DP)              :: vline0
  REAL(DP)              :: vline1
  REAL(DP)              :: sedge
  REAL(DP)              :: cedge
  REAL(DP), ALLOCATABLE :: x21(:,:)
  REAL(DP), ALLOCATABLE :: h1(:)
  !
  REAL(DP),   PARAMETER :: STEP_FUNC_THR = 1.0E-6_DP
  !
  EXTERNAL :: dgemv
  !
  ! ... number of sites in solvents
  nq = get_nuniq_in_solVs()
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzs < rismt%cfft%dfftt%nr3) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzl < rismt%lfft%nrz) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... is one-hand ?
  IF (rismt%lfft%xright .AND. rismt%lfft%xleft) THEN
    IF (prepare) THEN
      rism%cdzs = 0.0_DP
      rism%hdz  = 0.0_DP
    ELSE
      rism%cdza = 0.0_DP
    END IF
    !
    ierr = IERR_RISM_NULL
    RETURN
  END IF
  !
  ! ... beta = 1 / (kB * T)
  beta = 1.0_DP / K_BOLTZMANN_RY / rismt%temp
  !
  ! ... set domain of z1 as index of long Z-stick (i.e. expanded cell)
  IF (rismt%lfft%xright) THEN
    IF (expand .OR. prepare) THEN
      izsta1 = rismt%lfft%izright_start
      izend1 = rismt%lfft%nrz
    ELSE
      izsta1 = rismt%lfft%izright_start
      izend1 = rismt%lfft%izcell_end
    END IF
    !
    izsolv = izsta1
    !
    IF (rismt%lfft%gxystart > 1) THEN
      ssign  = -1.0_DP
      vline0 = AIMAG(rismt%vleft(1))
      vline1 = DBLE( rismt%vleft(1)) / alat
    END IF
    !
  ELSE !IF (rismt%lfft%xleft) THEN
    IF (expand .OR. prepare) THEN
      izsta1 = 1
      izend1 = rismt%lfft%izleft_end
    ELSE
      izsta1 = rismt%lfft%izcell_start
      izend1 = rismt%lfft%izleft_end
    END IF
    !
    izsolv = izend1
    !
    IF (rismt%lfft%gxystart > 1) THEN
      ssign  = +1.0_DP
      vline0 = AIMAG(rismt%vright(1))
      vline1 = DBLE( rismt%vright(1)) / alat
    END IF
  END IF
  !
  ! ... count integral points along z1
  nzint1 = izend1 - izsta1 + 1
  !
  ! ... properties about length (in a.u.)
  z0    = alat * 0.5_DP * at(3, 3)
  zstep = alat * rismt%lfft%zstep
  zoffs = alat * (rismt%lfft%zleft + rismt%lfft%zoffset)
  zedge = zoffs + zstep * DBLE(izsolv - 1)
  !
  IF (prepare) THEN
    ! ...
    ! ... Prepare intermediate integrations
    ! ..............................................................................
    !
    ! ... set domain of z2 as index of long Z-stick (i.e. expanded cell)
    IF (rismt%lfft%xright) THEN
      izsta2 = rismt%lfft%izright_start
      izend2 = rismt%lfft%izright_end
      !
    ELSE !IF (rismt%lfft%xleft) THEN
      izsta2 = rismt%lfft%izleft_start
      izend2 = rismt%lfft%izleft_end
    END IF
    !
    ! ... count integral points along z2
    nzint2 = izend2 - izsta2 + 1
    !
    ! ... allocate working memory
    IF (nzint2 * nzint1 > 0) THEN
      ALLOCATE(x21(nzint2, nzint1))
    END IF
    !
    ! ... calculate cdzs: step function of dipole part
    IF (rismt%nrzs > 0) THEN
      rismt%cdzs = 0.0_DP
    END IF
    !
    IF (rismt%lfft%gxystart > 1) THEN
!$omp parallel do default(shared) private(iz2, iiz2, z)
      DO iz2 = izsta2, izend2
        iiz2 = iz2 - rismt%lfft%izcell_start + 1
        z = zoffs + zstep * DBLE(iz2 - 1)
        rismt%cdzs(iiz2) = 0.5_DP * (1.0_DP + ssign * SIN(0.5_DP * pi * z / z0))
      END DO
!$omp end parallel do
    END IF
    !
    ! ... calculate hdz, for all solvent-pairs
    IF (rismt%nrzl * rismt%nsite * rismt%mp_site%nsite > 0) THEN
      rismt%hdz = 0.0_DP
    END IF
    !
    DO iq1 = 1, nq
      DO iq2 = rismt%mp_site%isite_start, rismt%mp_site%isite_end
        iiq2 = iq2 - rismt%mp_site%isite_start + 1
        !
        IF (nzint2 * nzint1 > 0) THEN
          x21 = 0.0_DP
        END IF
        !
        IF (rismt%lfft%gxystart > 1) THEN
          ! ... set solvent susceptibility x21
!$omp parallel do default(shared) private(iz1, iz2, izint1, izint2, izdelt)
          DO iz2 = izsta2, izend2
            izint2 = iz2 - izsta2 + 1
            DO iz1 = izsta1, izend1
              izint1 = iz1 - izsta1 + 1
              izdelt = ABS(iz1 - iz2) + 1
              x21(izint2, izint1) = rismt%xgs(izdelt, iiq2, iq1)
            END DO
          END DO
!$omp end parallel do
          !
          ! ... convolute cdzs and x21 -> hdz
          IF (nzint2 * nzint1 > 0) THEN
            CALL dgemv('T', nzint2, nzint1, zstep, x21, nzint2, &
                     & rismt%cdzs(izsta2 - rismt%lfft%izcell_start + 1), 1, 0.0_DP, &
                     & rismt%hdz(izsta1, iiq2, iq1), 1)
          END IF
        END IF
        !
      END DO
    END DO
    !
    IF (rismt%nrzs > 0) THEN
      CALL mp_sum(rismt%cdzs, rismt%mp_site%intra_sitg_comm)
    END IF
    !
    IF (rismt%nrzl * rismt%nsite * rismt%mp_site%nsite > 0) THEN
      CALL mp_sum(rismt%hdz, rismt%mp_site%intra_sitg_comm)
    END IF
    !
    ! ... deallocate working memory
    IF (nzint1 * nzint2 > 0) THEN
      DEALLOCATE(x21)
    END IF
    !
  ELSE
    ! ...
    ! ... Calculate dipole part of total correlations
    ! ..............................................................................
    !
    ! ... allocate working memory
    IF (nzint1 > 0) THEN
      ALLOCATE(h1(nzint1))
    END IF
    !
    ! ... calculate cdza: amplitude of dipole part
    IF (rismt%nsite > 0) THEN
      rismt%cdza = 0.0_DP
    END IF
    !
    DO iq2 = rismt%mp_site%isite_start, rismt%mp_site%isite_end
      iiq2   = iq2 - rismt%mp_site%isite_start + 1
      iv2    = iuniq_to_isite(1, iq2)
      isolV2 = isite_to_isolV(iv2)
      iatom2 = isite_to_iatom(iv2)
      qv2    = solVs(isolV2)%charge(iatom2)
      !
      IF (rismt%lfft%gxystart > 1) THEN
        iiz2 = izsolv - rismt%lfft%izcell_start + 1
        sedge = 0.5_DP * (1.0_DP + ssign * SIN(0.5_DP * pi * zedge / z0))
        cedge = DBLE(rismt%csgz(iiz2, iiq2)) &
            & - beta * qv2 * DBLE(rismt%vlgz(izsolv)) &
            & + beta * qv2 * (vline1 * zedge + vline0)
        !
        IF (ABS(1.0_DP - sedge) > STEP_FUNC_THR) THEN
          rismt%cdza(iiq2) = cedge / (1.0_DP - sedge)
        ELSE
          rismt%cdza(iiq2) = cedge + rismt%cdza(iiq2) * sedge
        END IF
      END IF
    END DO
    !
    IF (rismt%nsite > 0) THEN
      CALL mp_sum(rismt%cdza, rismt%mp_site%intra_sitg_comm)
    END IF
    !
    ! ... Laue-RISM equation of dipole part
    DO iq1 = 1, nq
      IF (rismt%mp_site%isite_start <= iq1 .AND. iq1 <= rismt%mp_site%isite_end) THEN
        iiq1 = iq1 - rismt%mp_site%isite_start + 1
      ELSE
        iiq1 = 0
      END IF
      !
      IF (nzint1 > 0) THEN
        h1 = 0.0_DP
      END IF
      !
      DO iq2 = rismt%mp_site%isite_start, rismt%mp_site%isite_end
        iiq2 = iq2 - rismt%mp_site%isite_start + 1
        !
        ! ... calculate h1(z1)
        IF (rismt%lfft%gxystart > 1) THEN
!$omp parallel do default(shared) private(iz1, izint1)
          DO iz1 = izsta1, izend1
            izint1 = iz1 - izsta1 + 1
            h1(izint1) = h1(izint1) + rismt%cdza(iiq2) * rismt%hdz(iz1, iiq2, iq1)
          END DO
!$omp end parallel do
        END IF
        !
      END DO
      !
      IF (nzint1 > 0) THEN
        CALL mp_sum(h1, rismt%mp_site%inter_sitg_comm)
      END IF
      !
      IF (iiq1 > 0 .AND. rismt%lfft%gxystart > 1) THEN
        IF (expand) THEN
          ! ... add h1 -> hsgz
!$omp parallel do default(shared) private(iz1, izint1)
          DO iz1 = izsta1, izend1
            izint1 = iz1 - izsta1 + 1
            rismt%hsgz(iz1, iiq1) = rismt%hsgz(iz1, iiq1) + CMPLX(h1(izint1), 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
        ELSE
          ! ... add h1 -> hgz
!$omp parallel do default(shared) private(iz1, iiz1, izint1)
          DO iz1 = izsta1, izend1
            iiz1 = iz1 - rismt%lfft%izcell_start + 1
            izint1 = iz1 - izsta1 + 1
            rismt%hgz(iiz1, iiq1) = rismt%hgz(iiz1, iiq1) + CMPLX(h1(izint1), 0.0_DP, kind=DP)
          END DO
!$omp end parallel do
        END IF
      END IF
      !
    END DO
    !
    ! ... deallocate working memory
    IF (nzint1 > 0) THEN
      DEALLOCATE(h1)
    END IF
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE eqn_lauedipole
