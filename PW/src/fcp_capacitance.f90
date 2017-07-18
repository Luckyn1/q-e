!
! Copyright (C) 2001-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE fcp_capacitance(capacitance)
  !----------------------------------------------------------------------------
  !
  ! ... evaluate capacitance for FCP
  !
  USE cell_base,     ONLY : alat, at
  USE constants,     ONLY : fpi, e2, K_BOLTZMANN_RY
  USE esm,           ONLY : esm_bc, esm_w
  USE kinds,         ONLY : DP
  USE rism3d_facade, ONLY : rism3t, rism3d_is_laue, rism3d_is_both_hands
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: capacitance
  !
  REAL(DP) :: fac
  REAL(DP) :: epsr
  REAL(DP) :: z0
  REAL(DP) :: area_xy
  REAL(DP) :: zsol
  REAL(DP) :: rho0
  REAL(DP) :: beta
  !
  ! ... set permittivity and length of z-axis
  !
  IF (TRIM(esm_bc) == 'bc2' .OR. TRIM(esm_bc) == 'bc3' .OR. TRIM(esm_bc) == 'bc4') THEN
     !
     ! ... Parallel plate capacitor
     !
     IF (TRIM(esm_bc) == 'bc2') THEN
        fac  = 2.0_DP
     ELSE
        fac  = 1.0_DP
     END IF
     !
     epsr = 1.0_DP
     z0   = 0.5_DP * alat * at(3, 3) + esm_w
     !
  ELSE IF (TRIM(esm_bc) == 'bc1' .AND. rism3d_is_laue()) THEN
     !
     ! ... Debye-Huckel's theory
     !
     beta = 1.0_DP / K_BOLTZMANN_RY / rism3t%temp
     !
     CALL get_solvent_data(epsr, zsol, rho0)
     !
     IF (rism3d_is_both_hands()) THEN
        fac = 2.0_DP
     ELSE
        fac = 1.0_DP
     END IF
     !
     z0 = SQRT(0.5_DP * (epsr / fpi / e2) / (beta * rho0 * zsol * zsol))
     !
  ELSE
     !
     CALL errore('fcp_capacitance', 'cannot evaluate capacitance', 1)
     !
  END IF
  !
  ! ... calculate capacitance
  !
  area_xy = alat * alat * ABS(at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1))
  !
  capacitance = fac * (epsr / fpi / e2) * area_xy / z0
  !
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE get_solvent_data(epsr, zsol, rho0)
    !----------------------------------------------------------------------------
    !
    USE constants, ONLY : eps8
    USE solvmol,   ONLY : nsolV, solVs, get_nuniq_in_solVs, &
                        & iuniq_to_nsite, iuniq_to_isite, isite_to_isolV, isite_to_iatom
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT) :: epsr
    REAL(DP), INTENT(OUT) :: zsol
    REAL(DP), INTENT(OUT) :: rho0
    !
    INTEGER               :: iq, nq
    INTEGER               :: iv, nv
    INTEGER               :: isolV
    INTEGER               :: iatom
    REAL(DP)              :: qv
    REAL(DP)              :: epsv
    REAL(DP)              :: rho1
    REAL(DP)              :: rho2
    REAL(DP)              :: rhov
    REAL(DP)              :: rhot
    REAL(DP), ALLOCATABLE :: qsol(:)
    !
    ALLOCATE(qsol(nsolV))
    !
    nq = get_nuniq_in_solVs()
    !
    qsol = 0.0_DP
    !
    DO iq = 1, nq
       !
       iv    = iuniq_to_isite(1, iq)
       nv    = iuniq_to_nsite(iq)
       isolV = isite_to_isolV(iv)
       iatom = isite_to_iatom(iv)
       qv    = solVs(isolV)%charge(iatom)
       !
       qsol(isolV) = qsol(isolV) + DBLE(nv) * qv
       !
    END DO
    !
    epsr = 0.0_DP
    zsol = 0.0_DP
    rho0 = 0.0_DP
    rhot = 0.0_DP
    !
    DO isolV = 1, nsolV
       !
       rho1 = solVs(isolV)%density
       rho2 = solVs(isolV)%subdensity
       rhov = 0.5_DP * (rho1 + rho2)
       !
       IF (ABS(qsol(isolV)) > eps8) THEN
          !
          ! ... this is an ion
          !
          IF (zsol < ABS(qsol(isolV))) THEN
             !
             zsol = ABS(qsol(isolV))
             rho0 = rhov
             !
          END IF
          !
       ELSE
          !
          ! ... this is a neutral molecule
          !
          epsv = solVs(isolV)%permittivity
          !
          IF (epsv < eps8) THEN
             !
             epsv = 1.0_DP
             !
          END IF
          !
          epsr = epsr + rhov * epsv
          rhot = rhot + rhov
          !
       END IF
       !
    END DO
    !
    IF (rhot > eps8) THEN
       !
       epsr = epsr / rhot
       !
    END IF
    !
    IF (epsr < eps8) epsr = 1.0_DP
    IF (zsol < eps8) zsol = 1.0_DP
    IF (rho0 < eps8) rho0 = 1.0_DP
    !
    DEALLOCATE(qsol)
    !
  END SUBROUTINE get_solvent_data
  !
END SUBROUTINE fcp_capacitance
