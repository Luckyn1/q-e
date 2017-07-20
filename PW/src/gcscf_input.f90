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
  ! ...  stored in modules input_parameters into internal modules of FCP
  !
  USE constants,    ONLY : RYTOEV
  USE gcscf_module, ONLY : lgcscf_     => lgcscf,     &
                         & gcscf_mu_   => gcscf_mu,   &
                         & gcscf_g0_   => gcscf_g0,   &
                         & gcscf_beta_ => gcscf_beta, &
                         & gcscf_check
  USE kinds,        ONLY : DP
  !
  ! ... SYSTEM namelist
  !
  USE input_parameters, ONLY : lgcscf, gcscf_mu, gcscf_g0
  !
  IMPLICIT NONE
  !
  ! ... set variables from namelist
  !
  lgcscf_ = lgcscf
  !
  IF (.NOT. lgcscf_) THEN
     !
     RETURN
     !
  END IF
  !
  gcscf_mu_   = gcscf_mu / RYTOEV
  gcscf_g0_   = gcscf_g0
  gcscf_beta_ = gcscf_beta
  !
  ! ... check condition
  !
  CALL gcscf_check()
  !
END SUBROUTINE iosys_gcscf
