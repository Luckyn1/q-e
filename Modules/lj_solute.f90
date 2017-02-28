!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE lj_setup_solU_tau(rismt, rsmax, count_only, ierr)
  !--------------------------------------------------------------------------
  !
  ! ... setup coordinate of solute's atoms,
  ! ... which can contribute to Lennard-Jones potentials
  !
  USE cell_base, ONLY : at, bg, alat
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE ions_base, ONLY : nat, tau
  USE kinds,     ONLY : DP
  USE rism,      ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  USE solute,    ONLY : solU_nat, solU_tau, solU_ljsig, isup_to_iuni
  USE solvmol,   ONLY : nsolV, solVs
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  REAL(DP),        INTENT(IN)  :: rsmax
  LOGICAL,         INTENT(IN)  :: count_only
  INTEGER,         INTENT(OUT) :: ierr
  !
  LOGICAL               :: laue
  INTEGER               :: isolV
  INTEGER               :: iatom
  INTEGER               :: ia
  INTEGER               :: im1, im2, im3
  INTEGER               :: nm1, nm2, nm3
  REAL(DP)              :: rmax
  REAL(DP)              :: sv, su, suv
  REAL(DP)              :: rm1, rm2, rm3
  REAL(DP)              :: bgnrm(3)
  REAL(DP)              :: tau_tmp(3)
  REAL(DP), ALLOCATABLE :: tau_uni(:,:)
  !
  REAL(DP), EXTERNAL :: dnrm2
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM .AND. rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... alloc memory
  ALLOCATE(tau_uni(3, nat))
  !
  ! ... set variables
  laue = .FALSE.
  IF (rismt%itype == ITYPE_LAUERISM) THEN
    laue = .TRUE.
  END IF
  !
  bgnrm(1) = dnrm2(3, bg (1, 1), 1)
  bgnrm(2) = dnrm2(3, bg (1, 2), 1)
  bgnrm(3) = dnrm2(3, bg (1, 3), 1)
  !
  sv = 0.0_DP
  DO isolV = 1, nsolV
    DO iatom = 1, solVs(isolV)%natom
      sv = MAX(sv, solVs(isolV)%ljsig(iatom))
    END DO
  END DO
  !
  ! ... count maximum of cell size
  su = 0.0_DP
  DO ia = 1, nat
    su = MAX(su, solU_ljsig(ia))
  END DO
  !
  suv  = 0.5_DP * (sv + su)
  rmax = rsmax * suv / alat
  !
  nm1 = CEILING(bgnrm(1) * rmax)
  nm2 = CEILING(bgnrm(2) * rmax)
  IF (.NOT. laue) THEN
    nm3 = CEILING(bgnrm(3) * rmax)
  ELSE
    nm3 = 0
  END IF
  !
  ! ... wrap coordinates back into cell
  tau_uni = tau
  CALL cryst_to_cart(nat, tau_uni, bg, -1)
  IF (.NOT. laue) THEN
    tau_uni = tau_uni - FLOOR(tau_uni)
  ELSE
    tau_uni(1:2, :) = tau_uni(1:2, :) - FLOOR(tau_uni(1:2, :))
  END IF
  !
  ! ... set unit cell
  solU_nat = nat
  IF (.NOT. count_only) THEN
    DO ia = 1, nat
      solU_tau(:,  ia) = tau_uni(:, ia)
      isup_to_iuni(ia) = ia
    END DO
  END IF
  !
  ! ... set neighbor cells
  DO im1 = -nm1, nm1
    DO im2 = -nm2, nm2
      DO im3 = -nm3, nm3
        !
        IF (im1 == 0 .AND. im2 == 0 .AND. im3 == 0) THEN
          CYCLE
        END IF
        !
        DO ia = 1, nat
          su   = solU_ljsig(ia)
          suv  = 0.5_DP * (sv + su)
          rmax = rsmax * suv / alat
          rm1  = bgnrm(1) * rmax
          rm2  = bgnrm(2) * rmax
          rm3  = bgnrm(3) * rmax
          !
          tau_tmp(1) = tau_uni(1, ia) + DBLE(im1)
          tau_tmp(2) = tau_uni(2, ia) + DBLE(im2)
          tau_tmp(3) = tau_uni(3, ia) + DBLE(im3)
          !
          IF (tau_tmp(1) < -rm1 .OR. (rm1 + 1.0_DP) < tau_tmp(1)) THEN
            CYCLE
          END IF
          !
          IF (tau_tmp(2) < -rm2 .OR. (rm2 + 1.0_DP) < tau_tmp(2)) THEN
            CYCLE
          END IF
          !
          IF (.NOT. laue) THEN
            IF (tau_tmp(3) < -rm3 .OR. (rm3 + 1.0_DP) < tau_tmp(3)) THEN
              CYCLE
            END IF
          END IF
          !
          solU_nat = solU_nat + 1
          IF (.NOT. count_only) THEN
            solU_tau(:,  solU_nat) = tau_tmp(:)
            isup_to_iuni(solU_nat) = ia
          END IF
        END DO
        !
      END DO
    END DO
  END DO
  !
  IF (.NOT. count_only) THEN
    CALL cryst_to_cart(solU_nat, solU_tau, at, +1)
  END IF
  !
  ! ... dealloc memory
  DEALLOCATE(tau_uni)
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE lj_setup_solU_tau
!
!--------------------------------------------------------------------------
SUBROUTINE lj_setup_solU_vlj(rismt, rsmax, ierr)
  !--------------------------------------------------------------------------
  !
  ! ... calculate solute-solvent's Lennard-Jones potential
  ! ...
  ! ...   ----           [ [  sig  ]12    [  sig  ]6 ]
  ! ...   >    4 * esp * [ [-------]   -  [-------]  ]
  ! ...   ----           [ [|r - R|]      [|r - R|]  ]
  ! ...    R
  !
  USE err_rism, ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,    ONLY : DP
  USE rism,     ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  USE solvmol,  ONLY : get_nuniq_in_solVs
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  REAL(DP),        INTENT(IN)    :: rsmax
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER :: nq
  INTEGER :: iq
  LOGICAL :: laue
  !
  ! ... number of sites in solvents
  nq = get_nuniq_in_solVs()
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM .AND. rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nr < rismt%cfft%dfftt%nnr) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... Laue-RISM or not
  laue = .FALSE.
  IF (rismt%itype == ITYPE_LAUERISM) THEN
    laue = .TRUE.
  END IF
  !
  ! ... calculate Lennard-Jones
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    CALL lj_setup_solU_vlj_x(iq, rismt, rsmax, laue)
  END DO
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE lj_setup_solU_vlj
!
!--------------------------------------------------------------------------
SUBROUTINE lj_setup_solU_vlj_x(iq, rismt, rsmax, laue)
  !--------------------------------------------------------------------------
  !
  ! ... calculate solute-solvent's Lennard-Jones potential
  ! ... for a solvent's site
  !
  USE constants, ONLY : eps6
  USE cell_base, ONLY : alat, at
  USE kinds,     ONLY : DP
  USE rism,      ONLY : rism_type
  USE solute,    ONLY : solU_nat, solU_tau, solU_ljeps, solU_ljsig, isup_to_iuni
  USE solvmol,   ONLY : iuniq_to_isite, isite_to_isolV, isite_to_iatom, solVs
  !
  IMPLICIT NONE
  !
  INTEGER,         INTENT(IN)    :: iq
  TYPE(rism_type), INTENT(INOUT) :: rismt
  REAL(DP),        INTENT(IN)    :: rsmax
  LOGICAL,         INTENT(IN)    :: laue
  !
  REAL(DP), PARAMETER :: RSMIN = eps6
  !
  INTEGER  :: iiq
  INTEGER  :: iv
  INTEGER  :: isolV
  INTEGER  :: iatom
  INTEGER  :: ir
  INTEGER  :: idx
  INTEGER  :: idx0
  INTEGER  :: i3min
  INTEGER  :: i3max
  INTEGER  :: i1, i2, i3
  INTEGER  :: n1, n2, n3
  INTEGER  :: nx1, nx2, nx3
  INTEGER  :: ia, iia
  REAL(DP) :: tau_r(3)
  REAL(DP) :: r1, r2, r3
  REAL(DP) :: ev, eu, euv
  REAL(DP) :: sv, su, suv
  REAL(DP) :: rmax
  REAL(DP) :: rmin
  REAL(DP) :: ruv2
  REAL(DP) :: xuv, yuv, zuv
  REAL(DP) :: sr2, sr6, sr12
  REAL(DP) :: vlj
  !
  ! ... FFT box
  n1  = rismt%cfft%dfftt%nr1
  n2  = rismt%cfft%dfftt%nr2
  n3  = rismt%cfft%dfftt%nr3
  nx1 = rismt%cfft%dfftt%nr1x
  nx2 = rismt%cfft%dfftt%nr2x
  nx3 = rismt%cfft%dfftt%nr3x
  !
  ! ... solvent properties
  iiq   = iq - rismt%mp_site%isite_start + 1
  iv    = iuniq_to_isite(1, iq)
  isolV = isite_to_isolV(iv)
  iatom = isite_to_iatom(iv)
  sv    = solVs(isolV)%ljsig(iatom)
  ev    = solVs(isolV)%ljeps(iatom)
  !
  ! ... calculate potential on each FFT grid
  idx0  = nx1 * nx2 * rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
  i3min = rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
  i3max = rismt%cfft%dfftt%npp(rismt%cfft%dfftt%mype + 1) + i3min
  !
!$omp parallel do default(shared) private(ir, idx, i1, i2, i3, r1, r2, r3, tau_r, vlj, &
!$omp             ia, iia, su, suv, rmax, rmin, xuv, yuv, zuv, ruv2, eu, euv, sr2, sr6, sr12)
  DO ir = 1, rismt%cfft%dfftt%nnr
    !
    ! ... create coordinate of a FFT grid
    idx = idx0 + ir - 1
    i3  = idx / (nx1 * nx2)
    IF (i3 < i3min .OR. i3 >= i3max .OR. i3 >= n3) THEN
      rismt%uljr(ir, iiq) = 0.0_DP
      CYCLE
    END IF
    !
    idx = idx - (nx1 * nx2) * i3
    i2  = idx / nx1
    IF (i2 >= n2) THEN
      rismt%uljr(ir, iiq) = 0.0_DP
      CYCLE
    END IF
    !
    idx = idx - nx1 * i2
    i1  = idx
    IF (i1 >= n1) THEN
      rismt%uljr(ir, iiq) = 0.0_DP
      CYCLE
    END IF
    !
    r1 = DBLE(i1) / DBLE(n1)
    r2 = DBLE(i2) / DBLE(n2)
    r3 = DBLE(i3) / DBLE(n3)
    IF (laue) THEN
      IF (i3 >= (n3 - (n3 / 2))) THEN
        r3 =  r3 - 1.0_DP
      END IF
    END IF
    !
    tau_r(:) = r1 * at(:, 1) + r2 * at(:, 2) + r3 * at(:, 3)
    !
    ! ... contribution from each solute's atom
    vlj = 0.0_DP
    !
    DO ia = 1, solU_nat
      iia  = isup_to_iuni(ia)
      su   = solU_ljsig(iia)
      suv  = 0.5_DP * (sv + su)
      rmax = rsmax * suv / alat
      rmin = RSMIN * suv / alat
      xuv  = tau_r(1) - solU_tau(1, ia)
      yuv  = tau_r(2) - solU_tau(2, ia)
      zuv  = tau_r(3) - solU_tau(3, ia)
      ruv2 = xuv * xuv + yuv * yuv + zuv * zuv
      IF (ruv2 > (rmax * rmax)) THEN
        CYCLE
      END IF
      IF (ruv2 < (rmin * rmin)) THEN
        ruv2 = rmin * rmin
      END IF
      !
      eu   = solU_ljeps(iia)
      euv  = SQRT(ev * eu)
      sr2  = suv * suv / ruv2 / alat / alat
      sr6  = sr2 * sr2 * sr2
      sr12 = sr6 * sr6
      vlj  = vlj + 4.0_DP * euv * (sr12 - sr6)
    END DO
    !
    rismt%uljr(ir, iiq) = vlj
    !
  END DO
!$omp end parallel do
  !
END SUBROUTINE lj_setup_solU_vlj_x
!
!--------------------------------------------------------------------------
SUBROUTINE lj_get_force(rismt, force, rsmax, ierr)
  !--------------------------------------------------------------------------
  !
  ! ... calculate solute-solvent's Lennard-Jones force
  ! ...
  ! ...           /                     [ 12*(x - X) [  sig  ]12    6*(x - X) [  sig  ]6 ]
  ! ...   - rho * | dr g(r) * 4 * esp * [ ---------- [-------]   -  --------- [-------]  ]
  ! ...           /                     [ |r - R|^2  [|r - R|]      |r - R|^2 [|r - R|]  ]
  ! ...
  !
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE ions_base, ONLY : nat
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE rism,      ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  USE solvmol,   ONLY : get_nuniq_in_solVs
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  REAL(DP),        INTENT(OUT) :: force(3, nat)
  REAL(DP),        INTENT(IN)  :: rsmax
  INTEGER,         INTENT(OUT) :: ierr
  !
  INTEGER :: nq
  INTEGER :: iq
  LOGICAL :: laue
  !
  ! ... number of sites in solvents
  nq = get_nuniq_in_solVs()
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM .AND. rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nr < rismt%cfft%dfftt%nnr) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... Laue-RISM or not
  laue = .FALSE.
  IF (rismt%itype == ITYPE_LAUERISM) THEN
    laue = .TRUE.
  END IF
  !
  ! ... calculate Lennard-Jones force
  force = 0.0_DP
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    CALL lj_get_force_x(iq, rismt, force, rsmax, laue)
  END DO
  !
  CALL mp_sum(force, rismt%mp_site%inter_sitg_comm)
  CALL mp_sum(force, rismt%mp_site%intra_sitg_comm)
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE lj_get_force
!
!--------------------------------------------------------------------------
SUBROUTINE lj_get_force_x(iq, rismt, force, rsmax, laue)
  !--------------------------------------------------------------------------
  !
  ! ... calculate solute-solvent's Lennard-Jones force
  ! ... for a solvent's site
  !
  USE constants, ONLY : eps6
  USE cell_base, ONLY : alat, at, omega
#if defined(__OPENMP)
  USE ions_base, ONLY : nat
#endif
  USE kinds,     ONLY : DP
  USE rism,      ONLY : rism_type
  USE solute,    ONLY : solU_nat, solU_tau, solU_ljeps, solU_ljsig, isup_to_iuni
  USE solvmol,   ONLY : solVs, iuniq_to_isite, iuniq_to_nsite, &
                      & isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  INTEGER,         INTENT(IN)    :: iq
  TYPE(rism_type), INTENT(IN)    :: rismt
  REAL(DP),        INTENT(INOUT) :: force(1:3, 1:*)
  REAL(DP),        INTENT(IN)    :: rsmax
  LOGICAL,         INTENT(IN)    :: laue
  !
  REAL(DP), PARAMETER :: RSMIN = eps6
  !
  INTEGER  :: iiq
  INTEGER  :: iv
  INTEGER  :: nv
  INTEGER  :: isolV
  INTEGER  :: iatom
  INTEGER  :: ir
  INTEGER  :: idx
  INTEGER  :: idx0
  INTEGER  :: i3min
  INTEGER  :: i3max
  INTEGER  :: i1, i2, i3
  INTEGER  :: n1, n2, n3
  INTEGER  :: nx1, nx2, nx3
  INTEGER  :: iz
  INTEGER  :: ia, iia
  REAL(DP) :: fac
  REAL(DP) :: weight
  REAL(DP) :: rho_right
  REAL(DP) :: rho_left
  REAL(DP) :: rhog
  REAL(DP) :: tau_r(3)
  REAL(DP) :: r1, r2, r3
  REAL(DP) :: ev, eu, euv
  REAL(DP) :: sv, su, suv
  REAL(DP) :: rmax
  REAL(DP) :: rmin
  REAL(DP) :: ruv2
  REAL(DP) :: xuv, yuv, zuv
  REAL(DP) :: sr2, sr6, sr12
#if defined(__OPENMP)
  REAL(DP), ALLOCATABLE :: fromp(:,:)
#endif
  !
  ! ... FFT box
  n1  = rismt%cfft%dfftt%nr1
  n2  = rismt%cfft%dfftt%nr2
  n3  = rismt%cfft%dfftt%nr3
  nx1 = rismt%cfft%dfftt%nr1x
  nx2 = rismt%cfft%dfftt%nr2x
  nx3 = rismt%cfft%dfftt%nr3x
  !
  ! ... solvent properties
  iiq       = iq - rismt%mp_site%isite_start + 1
  iv        = iuniq_to_isite(1, iq)
  nv        = iuniq_to_nsite(iq)
  isolV     = isite_to_isolV(iv)
  iatom     = isite_to_iatom(iv)
  sv        = solVs(isolV)%ljsig(iatom)
  ev        = solVs(isolV)%ljeps(iatom)
  rho_right = DBLE(nv) * solVs(isolV)%density
  rho_left  = DBLE(nv) * solVs(isolV)%subdensity
  !
  ! ... calculate potential on each FFT grid
  idx0  = nx1 * nx2 * rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
  i3min = rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
  i3max = rismt%cfft%dfftt%npp(rismt%cfft%dfftt%mype + 1) + i3min
  !
  weight = omega / DBLE(n1 * n2 * n3)
  !
!$omp parallel default(shared) private(ir, idx, i1, i2, i3, r1, r2, r3, tau_r, rhog, &
!$omp          ia, iia, su, suv, rmax, rmin, xuv, yuv, zuv, ruv2, eu, euv, sr2, sr6, sr12, &
!$omp          fac, fromp)
#if defined(__OPENMP)
  ALLOCATE(fromp(3, nat))
  fromp = 0.0_DP
#endif
!$omp do
  DO ir = 1, rismt%cfft%dfftt%nnr
    !
    ! ... create coordinate of a FFT grid
    idx = idx0 + ir - 1
    i3  = idx / (nx1 * nx2)
    IF (i3 < i3min .OR. i3 >= i3max .OR. i3 >= n3) THEN
      CYCLE
    END IF
    !
    idx = idx - (nx1 * nx2) * i3
    i2  = idx / nx1
    IF (i2 >= n2) THEN
      CYCLE
    END IF
    !
    idx = idx - nx1 * i2
    i1  = idx
    IF (i1 >= n1) THEN
      CYCLE
    END IF
    !
    r1 = DBLE(i1) / DBLE(n1)
    r2 = DBLE(i2) / DBLE(n2)
    r3 = DBLE(i3) / DBLE(n3)
    IF (laue) THEN
      IF (i3 >= (n3 - (n3 / 2))) THEN
        r3 =  r3 - 1.0_DP
      END IF
    END IF
    !
    tau_r(:) = r1 * at(:, 1) + r2 * at(:, 2) + r3 * at(:, 3)
    !
    IF (.NOT. laue) THEN
      rhog = rho_right * rismt%gr(ir, iiq)
    ELSE
      IF (i3 < (n3 - (n3 / 2))) THEN
        iz = i3 + (n3 / 2)
      ELSE
        iz = i3 - n3 + (n3 / 2)
      END IF
      iz = iz + rismt%lfft%izcell_start
      IF (iz > rismt%lfft%izleft_end) THEN
        rhog = rho_right * rismt%gr(ir, iiq)
      ELSE
        rhog = rho_left  * rismt%gr(ir, iiq)
      END IF
    END IF
    !
    ! ... contribution from each solute's atom
    DO ia = 1, solU_nat
      iia  = isup_to_iuni(ia)
      su   = solU_ljsig(iia)
      suv  = 0.5_DP * (sv + su)
      rmax = rsmax * suv / alat
      rmin = RSMIN * suv / alat
      xuv  = tau_r(1) - solU_tau(1, ia)
      yuv  = tau_r(2) - solU_tau(2, ia)
      zuv  = tau_r(3) - solU_tau(3, ia)
      ruv2 = xuv * xuv + yuv * yuv + zuv * zuv
      IF (ruv2 > (rmax * rmax)) THEN
        CYCLE
      END IF
      IF (ruv2 < (rmin * rmin)) THEN
        ruv2 = rmin * rmin
      END IF
      !
      eu   = solU_ljeps(iia)
      euv  = SQRT(ev * eu)
      sr2  = suv * suv / ruv2 / alat / alat
      sr6  = sr2 * sr2 * sr2
      sr12 = sr6 * sr6
      fac  = 4.0_DP * euv * (12.0_DP * sr12 - 6.0_DP * sr6) / ruv2 / alat
#if defined(__OPENMP)
      fromp(1, iia) = fromp(1, iia) - weight * rhog * fac * xuv
      fromp(2, iia) = fromp(2, iia) - weight * rhog * fac * yuv
      fromp(3, iia) = fromp(3, iia) - weight * rhog * fac * zuv
#else
      force(1, iia) = force(1, iia) - weight * rhog * fac * xuv
      force(2, iia) = force(2, iia) - weight * rhog * fac * yuv
      force(3, iia) = force(3, iia) - weight * rhog * fac * zuv
#endif
    END DO
    !
  END DO
!$omp end do
#if defined(__OPENMP)
!$omp critical
  force(1:3, 1:nat) = force(1:3, 1:nat) + fromp(1:3, 1:nat)
!$omp end critical
  DEALLOCATE(fromp)
#endif
!$omp end parallel
  !
END SUBROUTINE lj_get_force_x
!
!--------------------------------------------------------------------------
SUBROUTINE lj_get_stress(rismt, sigma, rsmax, ierr)
  !--------------------------------------------------------------------------
  !
  ! ... calculate solute-solvent's Lennard-Jones stress
  !
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE rism,      ONLY : rism_type, ITYPE_3DRISM, ITYPE_LAUERISM
  USE solvmol,   ONLY : get_nuniq_in_solVs
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)  :: rismt
  REAL(DP),        INTENT(OUT) :: sigma(3, 3)
  REAL(DP),        INTENT(IN)  :: rsmax
  INTEGER,         INTENT(OUT) :: ierr
  !
  INTEGER :: nq
  INTEGER :: iq
  LOGICAL :: laue
  !
  ! ... number of sites in solvents
  nq = get_nuniq_in_solVs()
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_3DRISM .AND. rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%mp_site%nsite < nq) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nr < rismt%cfft%dfftt%nnr) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... Laue-RISM or not
  laue = .FALSE.
  IF (rismt%itype == ITYPE_LAUERISM) THEN
    laue = .TRUE.
  END IF
  !
  ! ... calculate Lennard-Jones stress
  sigma = 0.0_DP
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    CALL lj_get_stress_x(iq, rismt, sigma, rsmax, laue)
  END DO
  !
  CALL mp_sum(sigma, rismt%mp_site%inter_sitg_comm)
  CALL mp_sum(sigma, rismt%mp_site%intra_sitg_comm)
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
END SUBROUTINE lj_get_stress
!
!--------------------------------------------------------------------------
SUBROUTINE lj_get_stress_x(iq, rismt, sigma, rsmax, laue)
  !--------------------------------------------------------------------------
  !
  ! ... calculate solute-solvent's Lennard-Jones force
  ! ... for a solvent's site
  !
  USE constants, ONLY : eps6
  USE cell_base, ONLY : alat, at, omega
  USE kinds,     ONLY : DP
  USE rism,      ONLY : rism_type
  USE solute,    ONLY : solU_nat, solU_tau, solU_ljeps, solU_ljsig, isup_to_iuni
  USE solvmol,   ONLY : solVs, iuniq_to_isite, iuniq_to_nsite, &
                      & isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  INTEGER,         INTENT(IN)    :: iq
  TYPE(rism_type), INTENT(IN)    :: rismt
  REAL(DP),        INTENT(INOUT) :: sigma(3, 3)
  REAL(DP),        INTENT(IN)    :: rsmax
  LOGICAL,         INTENT(IN)    :: laue
  !
  REAL(DP), PARAMETER :: RSMIN = eps6
  !
  INTEGER  :: iiq
  INTEGER  :: iv
  INTEGER  :: nv
  INTEGER  :: isolV
  INTEGER  :: iatom
  INTEGER  :: ir
  INTEGER  :: idx
  INTEGER  :: idx0
  INTEGER  :: i3min
  INTEGER  :: i3max
  INTEGER  :: i1, i2, i3
  INTEGER  :: n1, n2, n3
  INTEGER  :: nx1, nx2, nx3
  INTEGER  :: iz
  INTEGER  :: ia, iia
  REAL(DP) :: fac
  REAL(DP) :: weight
  REAL(DP) :: rho_right
  REAL(DP) :: rho_left
  REAL(DP) :: rhog
  REAL(DP) :: tau_r(3)
  REAL(DP) :: r1, r2, r3
  REAL(DP) :: ev, eu, euv
  REAL(DP) :: sv, su, suv
  REAL(DP) :: rmax
  REAL(DP) :: rmin
  REAL(DP) :: ruv2
  REAL(DP) :: xuv, yuv, zuv
  REAL(DP) :: sr2, sr6, sr12
#if defined(__OPENMP)
  REAL(DP) :: sgomp(3, 3)
#endif
  !
  ! ... FFT box
  n1  = rismt%cfft%dfftt%nr1
  n2  = rismt%cfft%dfftt%nr2
  n3  = rismt%cfft%dfftt%nr3
  nx1 = rismt%cfft%dfftt%nr1x
  nx2 = rismt%cfft%dfftt%nr2x
  nx3 = rismt%cfft%dfftt%nr3x
  !
  ! ... solvent properties
  iiq       = iq - rismt%mp_site%isite_start + 1
  iv        = iuniq_to_isite(1, iq)
  nv        = iuniq_to_nsite(iq)
  isolV     = isite_to_isolV(iv)
  iatom     = isite_to_iatom(iv)
  sv        = solVs(isolV)%ljsig(iatom)
  ev        = solVs(isolV)%ljeps(iatom)
  rho_right = DBLE(nv) * solVs(isolV)%density
  rho_left  = DBLE(nv) * solVs(isolV)%subdensity
  !
  ! ... calculate potential on each FFT grid
  idx0  = nx1 * nx2 * rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
  i3min = rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
  i3max = rismt%cfft%dfftt%npp(rismt%cfft%dfftt%mype + 1) + i3min
  !
  weight = omega / DBLE(n1 * n2 * n3)
  !
!$omp parallel default(shared) private(ir, idx, i1, i2, i3, r1, r2, r3, tau_r, rhog, &
!$omp          ia, iia, su, suv, rmax, rmin, xuv, yuv, zuv, ruv2, eu, euv, sr2, sr6, sr12, &
!$omp          fac, sgomp)
#if defined(__OPENMP)
  sgomp = 0.0_DP
#endif
!$omp do
  DO ir = 1, rismt%cfft%dfftt%nnr
    !
    ! ... create coordinate of a FFT grid
    idx = idx0 + ir - 1
    i3  = idx / (nx1 * nx2)
    IF (i3 < i3min .OR. i3 >= i3max .OR. i3 >= n3) THEN
      CYCLE
    END IF
    !
    idx = idx - (nx1 * nx2) * i3
    i2  = idx / nx1
    IF (i2 >= n2) THEN
      CYCLE
    END IF
    !
    idx = idx - nx1 * i2
    i1  = idx
    IF (i1 >= n1) THEN
      CYCLE
    END IF
    !
    r1 = DBLE(i1) / DBLE(n1)
    r2 = DBLE(i2) / DBLE(n2)
    r3 = DBLE(i3) / DBLE(n3)
    IF (laue) THEN
      IF (i3 >= (n3 - (n3 / 2))) THEN
        r3 =  r3 - 1.0_DP
      END IF
    END IF
    !
    tau_r(:) = r1 * at(:, 1) + r2 * at(:, 2) + r3 * at(:, 3)
    !
    IF (.NOT. laue) THEN
      rhog = rho_right * rismt%gr(ir, iiq)
    ELSE
      IF (i3 < (n3 - (n3 / 2))) THEN
        iz = i3 + (n3 / 2)
      ELSE
        iz = i3 - n3 + (n3 / 2)
      END IF
      iz = iz + rismt%lfft%izcell_start
      IF (iz > rismt%lfft%izleft_end) THEN
        rhog = rho_right * rismt%gr(ir, iiq)
      ELSE
        rhog = rho_left  * rismt%gr(ir, iiq)
      END IF
    END IF
    !
    ! ... contribution from each solute's atom
    DO ia = 1, solU_nat
      iia  = isup_to_iuni(ia)
      su   = solU_ljsig(iia)
      suv  = 0.5_DP * (sv + su)
      rmax = rsmax * suv / alat
      rmin = RSMIN * suv / alat
      xuv  = tau_r(1) - solU_tau(1, ia)
      yuv  = tau_r(2) - solU_tau(2, ia)
      zuv  = tau_r(3) - solU_tau(3, ia)
      ruv2 = xuv * xuv + yuv * yuv + zuv * zuv
      IF (ruv2 > (rmax * rmax)) THEN
        CYCLE
      END IF
      IF (ruv2 < (rmin * rmin)) THEN
        ruv2 = rmin * rmin
      END IF
      !
      eu   = solU_ljeps(iia)
      euv  = SQRT(ev * eu)
      sr2  = suv * suv / ruv2 / alat / alat
      sr6  = sr2 * sr2 * sr2
      sr12 = sr6 * sr6

#if defined(__OPENMP)
      ! TODO
      ! TODO set sgomp
      ! TODO
#else
      ! TODO
      ! TODO set sigma
      ! TODO
#endif
    END DO
    !
  END DO
!$omp end do
#if defined(__OPENMP)
!$omp critical
  sigma = sigma + sgomp
!$omp end critical
#endif
!$omp end parallel
  !
END SUBROUTINE lj_get_stress_x