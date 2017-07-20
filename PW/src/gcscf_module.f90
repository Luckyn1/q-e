!
! Copyright (C) 2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE gcscf_module
  !--------------------------------------------------------------------------
  !
  ! ... This module controls Grand-Canonical SCF.
  ! ... JCP 146, 114104 (2017), R.Sundararaman, W.A.Goddard-III, T.A.Arias
  !
  USE esm,             ONLY : do_comp_esm, esm_bc
  USE exx,             ONLY : x_gamma_extrapolation
  USE fcp_module,      ONLY : lfcp
  USE fixed_occ,       ONLY : tfixed_occ
  USE funct,           ONLY : exx_is_active
  USE kinds,           ONLY : DP
  USE klist,           ONLY : tot_charge, lgauss, degauss, ltetra, two_fermi_energies
  USE rism_module,     ONLY : lrism
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... variables for FCP method
  LOGICAL           :: lgcscf   = .FALSE.  ! to calculate GC-SCF method, or not
  REAL(DP)          :: gcscf_mu = 0.0_DP   ! target Fermi energy (in Ry)
  REAL(DP)          :: gcscf_g0 = 0.0_DP   ! wavelength shift for mixing (in 1/bohr)
  !
  ! ... public components
  PUBLIC :: lgcscf
  PUBLIC :: gcscf_mu
  PUBLIC :: gcscf_g0
  !
  PUBLIC :: gcscf_check
  PUBLIC :: gcscf_iosys
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE gcscf_check()
    !----------------------------------------------------------------------------
    !
    ! ... check conditions,
    ! ... and stop program if incorrect conditions.
    !
    IMPLICIT NONE
    !
    ! ... only ESM
    !
    IF (.NOT. do_comp_esm) THEN
       CALL errore('gcscf_check', 'please set assume_isolated = "esm", for GC-SCF', 1)
    END IF
    !
    ! ... cannot use PBC
    !
    IF (TRIM(esm_bc) == 'pbc') THEN
       CALL errore('gcscf_check', 'please do not set esm_bc = "pbc", for GC-SCF', 1)
    END IF
    !
    ! ... cannot use Vacuum/Slab/Vacuum
    !
    IF (TRIM(esm_bc) == 'bc1' .AND. (.NOT. lrism)) THEN
       CALL errore('gcscf_check', 'cannot use ESM-BC1 without RISM, for GC-SCF', 1)
    END IF
    !
    ! ... correct Vexx(G=0) ?
    !
    IF (exx_is_active() .AND. (.NOT. x_gamma_extrapolation)) THEN
       CALL errore('gcscf_check', 'GC-SCF calculation requires Vexx(G=0)', 1)
    END IF
    !
    ! ... cannot use FCP
    !
    IF (lfcp) THEN
       CALL errore('gcscf_check', 'cannot use FCP with GC-SCF', 1)
    END IF
    !
    ! ... only for metallic system
    !
    IF (tfixed_occ .OR. ltetra .OR. (.NOT. lgauss) .OR. (degauss <= 0.0_DP)) THEN
       CALL errore('gcscf_check', 'please set occupations = "smearing", for GC-SCF', 1)
    END IF
    !
    ! ... only one Fermi energy
    !
    IF (two_fermi_energies) THEN
       CALL errore('gcscf_check', 'please do not set tot_magnetization, for GC-SCF', 1)
    END IF
    !
  END SUBROUTINE gcscf_check
  !
  !----------------------------------------------------------------------------
  SUBROUTINE gcscf_iosys()
    !----------------------------------------------------------------------------
    !
    ! ... set variables from input file
    !
    IMPLICIT NONE
    !
    CALL iosys_gcscf()
    !
  END SUBROUTINE gcscf_iosys
  !
END MODULE gcscf_module
