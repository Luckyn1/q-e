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
  USE cell_base,   ONLY : alat, at
  USE constants,   ONLY : fpi, e2
  USE esm,         ONLY : esm_bc, esm_w
  USE fcp_module,  ONLY : solvation_radius
  USE kinds,       ONLY : DP
  USE rism_module, ONLY : lrism
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(OUT) :: capacitance
  !
  REAL(DP) :: z0
  REAL(DP) :: area_xy
  !
  IF (TRIM(esm_bc) == 'bc2' .OR. TRIM(esm_bc) == 'bc3' .OR. TRIM(esm_bc) == 'bc4') THEN
     !
     z0 = 0.5_DP * alat * at(3, 3) + esm_w
     !
  ELSE IF (lrism .AND. (solvation_radius > 0.0_DP)) THEN
     !
     z0 = solvation_radius
     !
  ELSE
     !
     z0 = 0.5_DP * alat * at(3, 3)
     !
  END IF
  !
  area_xy = alat * alat * ABS(at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1))
  !
  capacitance = (1.0_DP / fpi / e2) * area_xy / z0
  !
END SUBROUTINE fcp_capacitance
