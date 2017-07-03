!
! Copyright (C) 2002-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE fcp_variables
  !--------------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  LOGICAL :: &
       lfcp             = .FALSE.   ! if .TRUE. FCP is optimized too.
  !
  REAL(DP) :: &
       fcp_mu           = 0.0_DP    ! target Fermi energy,
                                    ! in Hartree
  !
  LOGICAL :: &
       lfcp_linmin      = .FALSE., &! .TRUE. if fcp_scheme = "lm"
       lfcp_newton      = .FALSE.   ! .TRUE. if fcp_scheme = "newton"
  !
  REAL(DP) :: &
       fcp_thr          = 0.001_DP  ! convergence threshold for FCP relaxation,
                                    ! in Hartree
  !
  INTEGER :: &
       fcp_ndiis        = 4         ! size of DIIS for Newton algorithm
  !
  REAL(DP) :: &
       tot_charge_first = 0.0_DP,  &! total charge of the first image
       tot_charge_last  = 0.0_DP    ! total charge of the last image
  !
  REAL(DP) :: &
       solvation_radius = 6.0_DP    ! solvation radius to estimate capacity,
                                    ! in Bohr
  !
END MODULE fcp_variables
