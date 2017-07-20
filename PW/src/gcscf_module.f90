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
  USE constants,       ONLY : RYTOEV
  USE control_flags,   ONLY : imix
  USE esm,             ONLY : do_comp_esm, esm_bc
  USE exx,             ONLY : x_gamma_extrapolation
  USE fcp_module,      ONLY : lfcp
  USE fixed_occ,       ONLY : tfixed_occ
  USE funct,           ONLY : exx_is_active
  USE io_global,       ONLY : stdout
  USE ions_base,       ONLY : nat, ityp, zv
  USE kinds,           ONLY : DP
  USE klist,           ONLY : nks, wk, nelec, tot_charge, &
                            & lgauss, degauss, ltetra, two_fermi_energies
  USE mp,              ONLY : mp_sum
  USE mp_pools,        ONLY : inter_pool_comm
  USE rism_module,     ONLY : lrism
  USE wvfct,           ONLY : nbnd, wg
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... variables for FCP method
  LOGICAL  :: lgcscf     = .FALSE.  ! to calculate GC-SCF method, or not
  REAL(DP) :: gcscf_mu   = 0.0_DP   ! target Fermi energy (in Ry)
  REAL(DP) :: gcscf_g0   = 0.0_DP   ! wavelength shift for mixing (in 1/bohr)
  REAL(DP) :: gcscf_beta = 0.0_DP   ! mixing rate of Fermi energy
  !
  ! ... public components
  PUBLIC :: lgcscf
  PUBLIC :: gcscf_mu
  PUBLIC :: gcscf_g0
  PUBLIC :: gcscf_beta
  !
  PUBLIC :: gcscf_check
  PUBLIC :: gcscf_iosys
  PUBLIC :: gcscf_summary
  PUBLIC :: gcscf_calc_nelec
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
    ! ... only TF or local-TF mixing
    !
    IF (imix /= 1 .AND. imix /= 2) THEN
       CALL errore('gcscf_check', 'please set mixing_mode = "TF" or "local-TF", for GC-SCF', 1)
    END IF
    !
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
  !----------------------------------------------------------------------------
  SUBROUTINE gcscf_summary()
    !----------------------------------------------------------------------------
    !
    ! ... This routine writes on output the GC-SCF's information.
    !
    IMPLICIT NONE
    !
    IF (.NOT. lgcscf) RETURN
    !
    WRITE(stdout, '(/,5X,">>>> Grand-Canonical SCF is activated <<<<")' )
    WRITE(stdout, '(5X,"Initial Total Charge = ",F12.6," e"   )') tot_charge
    WRITE(stdout, '(5X,"Target Fermi Energy  = ",F12.6," Ry"  )') gcscf_mu
    WRITE(stdout, '(5X,"                     = ",F12.6," eV"  )') gcscf_mu * RYTOEV
    !
    FLUSH(stdout)
    !
  END SUBROUTINE gcscf_summary
  !
  !----------------------------------------------------------------------------
  SUBROUTINE gcscf_calc_nelec()
    !----------------------------------------------------------------------------
    !
    ! ... calculate number of electrons from weights
    !
    IMPLICIT NONE
    !
    INTEGER :: ik
    INTEGER :: ibnd
    !
    IF (.NOT. lgcscf) RETURN
    !
    nelec = 0.0_DP
    !
    DO ik = 1, nks
       !
       DO ibnd = 1, nbnd
          !
          nelec = nelec + wg(ibnd, ik)
          !
       END DO
       !
    END DO
    !
    CALL mp_sum(nelec, inter_pool_comm)
    !
    tot_charge = SUM(zv(ityp(1:nat))) - nelec
    !
  END SUBROUTINE gcscf_calc_nelec
  !
END MODULE gcscf_module
