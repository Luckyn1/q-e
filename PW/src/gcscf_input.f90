!
! Copyright (C) 2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
SUBROUTINE iosys_gcscf()
  !--------------------------------------------------------------------------
  !
  ! ...  Copy data read from input file (in subroutine "read_input_file") and
  ! ...  stored in modules input_parameters into internal modules of GC-SCF
  !
  USE constants,     ONLY : RYTOEV
  USE control_flags, ONLY : imix, diago_full_acc
  USE gcscf_module,  ONLY : gcscf_ignore_mun_ => gcscf_ignore_mun, &
                          & gcscf_mu_         => gcscf_mu,         &
                          & gcscf_gk_         => gcscf_gk,         &
                          & gcscf_gh_         => gcscf_gh,         &
                          & gcscf_beta_       => gcscf_beta,       &
                          & gcscf_eps_        => gcscf_eps,        &
                          & gcscf_check
  USE kinds,         ONLY : DP
  USE rism3d_facade, ONLY : lrism3d, conv_level
  !
  ! ... SYSTEM namelist
  !
  USE input_parameters, ONLY : gcscf_ignore_mun, &
                             & gcscf_mu, gcscf_gk, gcscf_gh, gcscf_beta
  !
  ! ... ELECTRONS namelist
  !
  USE input_parameters, ONLY : mixing_mode
  !
  IMPLICIT NONE
  !
  ! ... modify imix
  !
  IF (imix /= 1 .AND. imix /= 2) THEN
     !
     imix = 1  ! mixing_mode = 'TF'
     !imix = 2 ! mixing_mode = 'local-TF'
     !
     CALL infomsg('iosys', &
     & 'mixing_mode=' // TRIM(mixing_mode) // " is ignored, 'TF' is adopted")
     !
  END IF
  !
  ! ... modify diago_full_acc
  !
  IF (.NOT. diago_full_acc) THEN
     !
     diago_full_acc = .TRUE.
     !
     CALL infomsg('iosys', &
     & 'accurate eigenvalues are required for all states: diago_full_acc=.TRUE.')
     !
  END IF
  !
  ! ... modify conv_level
  !
  IF (lrism3d .AND. conv_level < 0.0_DP) THEN
     !
     conv_level = 0.2_DP
     !
     CALL infomsg('iosys', &
     & 'convergence of 3D-RISM is set: rism3d_conv_level=0.2')
     !
  END IF
  !
  ! ... set variables from namelist
  !
  gcscf_ignore_mun_ = gcscf_ignore_mun
  gcscf_mu_         = gcscf_mu / RYTOEV
  gcscf_gk_         = gcscf_gk
  gcscf_gh_         = gcscf_gh
  gcscf_beta_       = gcscf_beta
  gcscf_eps_        = gcscf_eps
  !
  ! ... check condition
  !
  CALL gcscf_check()
  !
END SUBROUTINE iosys_gcscf
