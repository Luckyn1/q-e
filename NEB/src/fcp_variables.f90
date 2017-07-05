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
  USE kinds, ONLY : DP
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
       solvation_radius = 6.0_DP    ! solvation radius to estimate capacity,
                                    ! in Bohr
  !
  REAL(DP), ALLOCATABLE :: &
       fcp_nelec(:),               &! the numbers of electrons
       fcp_ef(:),                   ! the Fermi energies
       fcp_dos(:)                   ! the DOSs on Fermi surfaces
  !
  CONTAINS
     !
     !--------------------------------------------------------------------------
     SUBROUTINE fcp_allocation()
       !--------------------------------------------------------------------------
       !
       USE path_variables, ONLY : num_of_images
       !
       IMPLICIT NONE
       !
       ALLOCATE( fcp_nelec( num_of_images ) )
       ALLOCATE( fcp_ef(    num_of_images ) )
       ALLOCATE( fcp_dos(   num_of_images ) )
       !
     END SUBROUTINE fcp_allocation
     !
     !--------------------------------------------------------------------------
     SUBROUTINE fcp_deallocation()
       !--------------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       IF ( ALLOCATED( fcp_nelec ) ) DEALLOCATE( fcp_nelec )
       IF ( ALLOCATED( fcp_ef ) )    DEALLOCATE( fcp_ef )
       IF ( ALLOCATED( fcp_dos ) )   DEALLOCATE( fcp_dos )
       !
     END SUBROUTINE fcp_deallocation
     !
END MODULE fcp_variables
