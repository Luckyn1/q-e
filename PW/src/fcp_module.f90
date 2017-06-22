!
! Copyright (C) 2001-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE fcp_module
  !--------------------------------------------------------------------------
  !
  ! ... Original version by Minoru Otani (AIST) and Nicephore Bonnet (AIST).
  ! ...
  ! ... This module controls the Fictitious Charge Particle (FCP) for constant-mu
  ! ... method developed by N. Bonnet, T. Morishita, O. Sugino, and M. Otani
  ! ... (see PRL 109, 266101 [2012]).
  ! ...
  ! ... Constant-mu scheme with the boundary condition 'bc2' and 'bc3' enables
  ! ... description of the system connected to a potentiostat which preserves
  ! ... the Fermi energy of the system as the target Fermi energy (mu).
  ! ...
  ! ... MDIIS and BFGS algorithms are implemented by S. Nishihara (2016-2017)
  ! ...
  ! ... . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  ! ...   This module is the facade of FCP calculations.
  ! ... . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  !
  USE constants,       ONLY : fpi, e2, RYTOEV
  USE control_flags,   ONLY : lbfgs, lmd, tr2
  USE cell_base,       ONLY : alat, at
  USE dynamics_module, ONLY : dt
  USE esm,             ONLY : do_comp_esm, esm_bc, esm_w
  USE ener,            ONLY : ef
  USE exx,             ONLY : x_gamma_extrapolation
  USE fcp_dynamics,    ONLY : fcpdyn_final, fcpdyn_update, &
                            & fcpdyn_set_verlet, fcpdyn_set_proj_verlet
  USE fcp_relaxation,  ONLY : fcprlx_final, fcprlx_update, &
                            & fcprlx_set_line_min, fcprlx_set_mdiis
  USE fixed_occ,       ONLY : tfixed_occ
  USE funct,           ONLY : exx_is_active
  USE io_global,       ONLY : stdout
  USE kinds,           ONLY : DP
  USE klist,           ONLY : tot_charge, lgauss, degauss, &
                            & ltetra, two_fermi_energies
  USE relax,           ONLY : starting_scf_threshold
  USE rism_module,     ONLY : lrism
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... variables for FCP method
  LOGICAL          :: lfcp             = .FALSE.  ! to calculate FCP method, or not
  REAL(DP)         :: fcp_mu           = 0.0_DP   ! target Fermi energy (in Ry)
  REAL(DP)         :: fcp_eps          = 0.0_DP   ! convergence threshold (in Ry)
  REAL(DP)         :: fcp_eps0         = 0.0_DP   ! initial convergence threshold (in Ry)
  CHARACTER(LEN=8) :: fcp_calc         = ''       ! type of calculation ( lm | mdiis | damp | verlet )
  REAL(DP)         :: solvation_radius = 0.0_DP   ! solvation radius to estimate capacity (in bohr)
  !
  ! ... public components
  PUBLIC :: lfcp
  PUBLIC :: fcp_mu
  PUBLIC :: fcp_eps
  PUBLIC :: fcp_eps0
  PUBLIC :: fcp_calc
  PUBLIC :: solvation_radius
  !
  PUBLIC :: fcp_check
  PUBLIC :: fcp_iosys
  PUBLIC :: fcp_summary
  PUBLIC :: fcp_relax
  PUBLIC :: fcp_verlet
  PUBLIC :: fcp_terminate
  PUBLIC :: fcp_new_conv_thr
  PUBLIC :: fcp_is_dynamics
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcp_check()
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
       CALL errore('fcp_check', 'please set assume_isolated = "esm", for FCP', 1)
    END IF
    !
    ! ... cannot use PBC
    !
    IF (TRIM(esm_bc) == 'pbc') THEN
       CALL errore('fcp_check', 'please do not set esm_bc = "pbc", for FCP', 1)
    END IF
    !
    ! ... cannot use Vacuum/Slab/Vacuum
    !
    IF (TRIM(esm_bc) == 'bc1' .AND. (.NOT. lrism)) THEN
       CALL errore('fcp_check', 'cannot use ESM-BC1 without RISM, for FCP', 1)
    END IF
    !
    ! ... correct Vexx(G=0) ?
    !
    IF (exx_is_active() .AND. (.NOT. x_gamma_extrapolation)) THEN
       CALL errore('fcp_check', 'FCP calculation requires Vexx(G=0)', 1)
    END IF
    !
    ! ... only for metallic system
    !
    IF (tfixed_occ .OR. ltetra .OR. (.NOT. lgauss) .OR. (degauss <= 0.0_DP)) THEN
       CALL errore('fcp_check', 'please set occupations = "smearing", for FCP', 1)
    END IF
    !
    ! ... only one Fermi energy
    !
    IF (two_fermi_energies) THEN
       CALL errore('fcp_check', 'please do not set tot_magnetization, for FCP', 1)
    END IF
    !
    ! ... must be relax or md
    !
    IF (.NOT. (lbfgs .OR. lmd)) THEN
       CALL errore('fcp_check', 'calculation has to be relax or md, for FCP', 1)
    END IF
    !
  END SUBROUTINE fcp_check
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcp_iosys(tfcp)
    !----------------------------------------------------------------------------
    !
    ! ... set variables from input file
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: tfcp
    !
    lfcp = tfcp
    !
    IF (lfcp) THEN
       !
       CALL iosys_fcp()
       !
    END IF
    !
  END SUBROUTINE fcp_iosys
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcp_summary()
    !----------------------------------------------------------------------------
    !
    ! ... This routine writes on output the FCP's information.
    !
    IMPLICIT NONE
    !
    IF (.NOT. lfcp) RETURN
    !
    IF (fcp_is_dynamics()) THEN
       !
       WRITE(stdout, '(/,5X,">>>>> FCP Dynamics is activated <<<<<<")' )
       !
    ELSE
       !
       WRITE(stdout, '(/,5X,">>>> FCP Relaxation is activated <<<<<")' )
       !
    END IF
    !
    WRITE(stdout, '(5X,"Initial Total Charge = ",F12.6," e"   )') tot_charge
    WRITE(stdout, '(5X,"Target Fermi Energy  = ",F12.6," Ry"  )') fcp_mu
    WRITE(stdout, '(5X,"                     = ",F12.6," eV"  )') fcp_mu * RYTOEV
    !
    FLUSH(stdout)
    !
  END SUBROUTINE fcp_summary
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcp_relax(conv)
    !----------------------------------------------------------------------------
    !
    ! ... relaxation of FCP. in the current version,
    ! ... Projected-Verlet, Line-Minimization and MDIIS are implemented.
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(INOUT) :: conv
    !
    REAL(DP) :: force
    REAL(DP) :: z0
    REAL(DP) :: area_xy
    REAL(DP) :: capacity
    REAL(DP) :: step_max
    !
    IF (.NOT. lfcp) RETURN
    !
    CALL fcp_check()
    !
    ! ... evaluate maximum step
    !
    force = fcp_mu - ef
    !
    IF (TRIM(esm_bc) == 'bc2' .OR. TRIM(esm_bc) == 'bc3' .OR. TRIM(esm_bc) == 'bc4') THEN
       z0 = 0.5_DP * alat * at(3, 3) + esm_w
    ELSE IF (lrism .AND. (solvation_radius > 0.0_DP)) THEN
       z0 = solvation_radius
    ELSE
       z0 = 0.5_DP * alat * at(3, 3)
    END IF
    !
    area_xy  = alat * alat * ABS(at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1))
    !
    capacity = (1.0_DP / fpi / e2) * area_xy / z0
    !
    step_max = ABS(capacity * force)
    !
    ! ... perform each algorithm
    !
    IF (TRIM(fcp_calc) == 'lm') THEN
       !
       ! ... update nelec by Line-Minimization
       !
       CALL fcprlx_set_line_min(fcp_eps, step_max)
       !
       CALL fcprlx_update(fcp_mu, conv)
       !
       !
    ELSE IF (TRIM(fcp_calc) == 'mdiis') THEN
       !
       ! ... update nelec by MDIIS
       !
       CALL fcprlx_set_mdiis(fcp_eps, step_max)
       !
       CALL fcprlx_update(fcp_mu, conv)
       !
    ELSE IF (TRIM(fcp_calc) == 'damp') THEN
       !
       ! ... update nelec by Projected-Verlet
       !
       CALL fcpdyn_set_proj_verlet(fcp_eps, step_max)
       !
       CALL fcpdyn_update(fcp_mu, dt, conv)
       !
    ELSE
       !
       CALL errore('fcp_relax', 'incorrect calculation: ' // TRIM(fcp_calc), 1)
       !
    END IF
    !
  END SUBROUTINE fcp_relax
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcp_verlet()
    !----------------------------------------------------------------------------
    !
    ! ... dynamics of FCP, using Verlet algorithm.
    !
    IMPLICIT NONE
    !
    IF (.NOT. lfcp) RETURN
    !
    CALL fcp_check()
    !
    IF (TRIM(fcp_calc) /= 'verlet') THEN
       !
       CALL errore('fcp_verlet', 'incorrect calculation: ' // TRIM(fcp_calc), 1)
       !
    END IF
    !
    ! ... update nelec by Verlet
    !
    CALL fcpdyn_set_verlet()
    !
    CALL fcpdyn_update(fcp_mu, dt)
    !
  END SUBROUTINE fcp_verlet
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcp_terminate()
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    IF (.NOT. lfcp) RETURN
    !
    IF (fcp_is_dynamics()) THEN
       !
       ! ... finalize fcp_dynamics
       !
       CALL fcpdyn_final()
       !
    ELSE
       !
       ! ... finalize fcp_relaxation
       !
       CALL fcprlx_final()
       !
    END IF
    !
  END SUBROUTINE fcp_terminate
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcp_new_conv_thr()
    !----------------------------------------------------------------------------
    !
    ! ... update convergence threshold.
    !
    IMPLICIT NONE
    !
    REAL(DP), PARAMETER :: TR2_EXPON = 0.50_DP
    !
    IF (.NOT. lfcp) RETURN
    !
    IF (fcp_eps > 0.0_DP .AND. fcp_eps0               > 0.0_DP .AND. &
        tr2     > 0.0_DP .AND. starting_scf_threshold > 0.0_DP) THEN
       !
       fcp_eps = fcp_eps0 * ((tr2 / starting_scf_threshold) ** TR2_EXPON)
       !
    ELSE
       !
       fcp_eps = fcp_eps0
       !
    END IF
    !
  END SUBROUTINE fcp_new_conv_thr
  !
  !----------------------------------------------------------------------------
  LOGICAL FUNCTION fcp_is_dynamics()
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    IF (TRIM(fcp_calc) == 'damp' .OR. TRIM(fcp_calc) == 'verlet') THEN
       !
       fcp_is_dynamics = .TRUE.
       !
    ELSE
       !
       fcp_is_dynamics = .FALSE.
       !
    END IF
    !
  END FUNCTION fcp_is_dynamics
  !
END MODULE fcp_module
