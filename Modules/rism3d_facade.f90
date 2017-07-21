!
! Copyright (C) 2015-2016 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE rism3d_facade
  !--------------------------------------------------------------------------
  !
  ! ... Facade (or Interface) of 3D-RISM's library.
  ! ... External codes, which utilize 3D-RISM, can access data and subroutines
  ! ... throught this module. Also, Laue-RISM is avairable.
  !
  USE constants,     ONLY : eps8, eps12
  USE control_flags, ONLY : gamma_only
  USE err_rism,      ONLY : stop_by_err_rism, IERR_RISM_NULL, &
                          & IERR_RISM_NOT_CONVERGED, IERR_RISM_NONZERO_CHARGE, &
                          & IERR_RISM_NOT_ANY_IONS
  USE gvect,         ONLY : ngl
  USE io_global,     ONLY : stdout
  USE io_rism_xml,   ONLY : write_3drism, read_3drism
  USE ions_base,     ONLY : nsp, nat
  USE kinds,         ONLY : DP
  USE mp,            ONLY : mp_sum
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE mp_images,     ONLY : intra_image_comm
  USE rism,          ONLY : rism_type, clean_rism_data, allocate_3drism, allocate_lauerism, &
                          & refresh_suscept_3drism, refresh_suscept_lauerism, &
                          & deallocate_rism, ITYPE_3DRISM, ITYPE_LAUERISM
  USE solute,        ONLY : update_solU
  USE solvmol,       ONLY : get_nuniq_in_solVs, iuniq_to_isite, iuniq_to_nsite, &
                          & isite_to_isolV, isite_to_iatom, solVs, nsolV
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... define constants
  INTEGER, PARAMETER     :: LEN_STR = 30
  !
  ! ... define variables
  LOGICAL                :: lrism3d        = .FALSE.  ! to calculate 3D-RISM, or not
  LOGICAL                :: has_any_corr   = .FALSE.  ! has nonzero correlations
  CHARACTER(LEN=LEN_STR) :: starting_corr  = ''       ! initial correlations: 'zero', 'file'
  INTEGER                :: niter          = 0        ! maximum number of iteration
  REAL(DP)               :: epsv           = 0.0_DP   ! convergence threshold
  REAL(DP)               :: starting_epsv  = 0.0_DP   ! convergence threshold (initial value)
  INTEGER                :: mdiis_size     = 0        ! size of MDIIS
  REAL(DP)               :: mdiis_step     = 0.0_DP   ! step of MDIIS
  REAL(DP)               :: ecutsolv       = 0.0_DP   ! energy cutoff for solvents on G-space (in Ry)
  INTEGER                :: conv_level     = 0        ! convergence level (0:low, 1:medium, 2:high)
  LOGICAL                :: planar_average = .FALSE.  ! calculate planar average, or not
  INTEGER                :: laue_nfit      = 0        ! #points to fit potential (for Laue-RISM)
  REAL(DP)               :: qsol           = 0.0_DP   ! total charge of solvent system (for Laue-RISM)
  REAL(DP)               :: expand_r       = -1.0_DP  ! expanding length of right (in alat, for Laue-RISM)
  REAL(DP)               :: expand_l       = -1.0_DP  ! expanding length of left (in alat, for Laue-RISM)
  REAL(DP)               :: starting_r     = 0.0_DP   ! starting position of right (in alat, for Laue-RISM)
  REAL(DP)               :: starting_l     = 0.0_DP   ! starting position of left (in alat, for Laue-RISM)
  REAL(DP)               :: buffer_r       = -1.0_DP  ! buffering length of right (in alat, for Laue-RISM)
  REAL(DP)               :: buffer_l       = -1.0_DP  ! buffering length of left (in alat, for Laue-RISM)
  LOGICAL                :: both_hands     = .FALSE.  ! apply both-hands calculation, or not (for Laue-RISM)
  INTEGER                :: ireference     = 0        ! reference type of potential (for Laue-RISM)
  !
  ! ..... reference types of potential
  INTEGER, PARAMETER :: IREFERENCE_NULL    = 0
  INTEGER, PARAMETER :: IREFERENCE_AVERAGE = 1
  INTEGER, PARAMETER :: IREFERENCE_RIGHT   = 2
  INTEGER, PARAMETER :: IREFERENCE_LEFT    = 3
  !
  ! ... define 3D-RISM's main data
  TYPE(rism_type) :: rism3t
  !
  ! ... public components
  PUBLIC :: lrism3d
  PUBLIC :: starting_corr
  PUBLIC :: niter
  PUBLIC :: epsv
  PUBLIC :: starting_epsv
  PUBLIC :: mdiis_size
  PUBLIC :: mdiis_step
  PUBLIC :: ecutsolv
  PUBLIC :: conv_level
  PUBLIC :: planar_average
  PUBLIC :: laue_nfit
  PUBLIC :: qsol
  PUBLIC :: expand_r
  PUBLIC :: expand_l
  PUBLIC :: starting_r
  PUBLIC :: starting_l
  PUBLIC :: buffer_r
  PUBLIC :: buffer_l
  PUBLIC :: both_hands
  PUBLIC :: ireference
  !
  PUBLIC :: IREFERENCE_NULL
  PUBLIC :: IREFERENCE_AVERAGE
  PUBLIC :: IREFERENCE_RIGHT
  PUBLIC :: IREFERENCE_LEFT
  !
  PUBLIC :: rism3t
  PUBLIC :: rism3d_iosys
  PUBLIC :: rism3d_summary
  PUBLIC :: rism3d_initialize
  PUBLIC :: rism3d_finalize
  PUBLIC :: rism3d_update_solute
  PUBLIC :: rism3d_prepare
  PUBLIC :: rism3d_reprepare
  PUBLIC :: rism3d_run
  PUBLIC :: rism3d_potential
  PUBLIC :: rism3d_force
  PUBLIC :: rism3d_stress
  PUBLIC :: rism3d_write_to_restart
  PUBLIC :: rism3d_read_to_restart
  PUBLIC :: rism3d_is_laue
  PUBLIC :: rism3d_set_laue
  PUBLIC :: rism3d_is_both_hands
  PUBLIC :: rism3d_printpot
  PUBLIC :: rism3d_print_clock
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism3d_iosys(trism, laue, init)
    !----------------------------------------------------------------------------
    !
    ! ... set variables
    ! ... if linit=.TRUE., initialize this module
    !
    IMPLICIT NONE
    !
    LOGICAL,           INTENT(IN) :: trism
    LOGICAL, OPTIONAL, INTENT(IN) :: laue
    LOGICAL, OPTIONAL, INTENT(IN) :: init
    !
    LOGICAL :: laue_
    LOGICAL :: init_
    !
    lrism3d = trism
    !
    IF (.NOT. lrism3d) THEN
      RETURN
    END IF
    !
    IF (PRESENT(laue)) THEN
      laue_ = laue
    ELSE
      laue_ = (rism3t%itype == ITYPE_LAUERISM)
    END IF
    !
    IF (PRESENT(init)) THEN
      init_ = init
    ELSE
      init_ = .TRUE.
    END IF
    !
    CALL iosys_3drism(laue_, init_)
    !
  END SUBROUTINE rism3d_iosys
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism3d_summary()
    !----------------------------------------------------------------------------
    !
    ! ... print conditions
    !
    IMPLICIT NONE
    !
    IF (.NOT. lrism3d) THEN
      RETURN
    END IF
    !
    CALL summary_3drism()
    !
  END SUBROUTINE rism3d_summary
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism3d_initialize(laue)
    !----------------------------------------------------------------------------
    !
    ! ... initialize this module
    !
    IMPLICIT NONE
    !
    LOGICAL, OPTIONAL, INTENT(IN) :: laue
    !
    LOGICAL  :: laue_
    INTEGER  :: nq
    REAL(DP) :: starting1_r
    REAL(DP) :: starting1_l
    REAL(DP) :: starting2_r
    REAL(DP) :: starting2_l
    !
    IF (.NOT. lrism3d) THEN
      RETURN
    END IF
    !
    IF (PRESENT(laue)) THEN
      laue_ = laue
    ELSE
      laue_ = (rism3t%itype == ITYPE_LAUERISM)
    END IF
    !
    nq = get_nuniq_in_solVs()
    !
    IF (.NOT. laue_) THEN
      CALL allocate_3drism(rism3t, nq, ecutsolv, intra_bgrp_comm, intra_image_comm)
      !
    ELSE
      starting1_r = starting_r - MAX(0.0_DP, buffer_r)
      starting1_l = starting_l + MAX(0.0_DP, buffer_l)
      starting2_r = starting_r
      starting2_l = starting_l
      CALL allocate_lauerism(rism3t, nq, ecutsolv, laue_nfit, expand_r, expand_l, &
                           & starting1_r, starting1_l, starting2_r, starting2_l, &
                           & both_hands, intra_bgrp_comm, intra_image_comm)
    END IF
    !
    IF (rism3t%itype == ITYPE_LAUERISM) THEN
      CALL check_solvent_is_neutral()
    END IF
    !
  END SUBROUTINE rism3d_initialize
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism3d_finalize(lall)
    !----------------------------------------------------------------------------
    !
    ! ... finalize this module
    !
    ! ... if lall=.TRUE., deallocate all data.
    ! ... if lall=.FALSE., deallocate data, which depend on FFT box.
    !
    IMPLICIT NONE
    !
    LOGICAL, INTENT(IN) :: lall
    !
    IF (.NOT. lrism3d) THEN
      RETURN
    END IF
    !
    CALL deallocate_rism(rism3t, lall)
    !
  END SUBROUTINE rism3d_finalize
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism3d_update_solute()
    !----------------------------------------------------------------------------
    !
    ! ... notify 3D-RISM's data that position of solute has been updated
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !
    IF (.NOT. lrism3d) THEN
      RETURN
    END IF
    !
    CALL start_clock('3DRISM_ions')
    !
    CALL update_solU(rism3t, ierr)
    !
    IF (ierr /= IERR_RISM_NULL) THEN
      CALL stop_by_err_rism('rism3d_update_solute', ierr)
    END IF
    !
    CALL stop_clock('3DRISM_ions')
    !
  END SUBROUTINE rism3d_update_solute
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism3d_prepare()
    !----------------------------------------------------------------------------
    !
    ! ... prepare 3D-RISM's iterative calculation
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !
    IF (.NOT. lrism3d) THEN
      RETURN
    END IF
    !
    CALL start_clock('3DRISM_pre')
    !
    ! ... initial calculation (L.J.-potential, solvent's susceptibility)
    CALL init_3drism(rism3t, .TRUE., both_hands, ierr)
    !
    IF (ierr /= IERR_RISM_NULL) THEN
      CALL stop_by_err_rism('rism3d_prepare', ierr)
    END IF
    !
    ! ... create initial correlation
    IF (TRIM(starting_corr) == 'file') THEN
      WRITE(stdout, '()')
      WRITE(stdout, '(5X,"Correlation function is read from file")')
      WRITE(stdout, '()')
      CALL clean_rism_data(rism3t)
      CALL rism3d_read_to_restart()
      has_any_corr = .TRUE.
      !
    ELSE !IF (TRIM(starting_corr) == 'zero') THEN
      CALL clean_rism_data(rism3t)
      has_any_corr = .FALSE.
    END IF
    !
    CALL stop_clock('3DRISM_pre')
    !
  END SUBROUTINE rism3d_prepare
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism3d_reprepare(at_old)
    !----------------------------------------------------------------------------
    !
    ! ... re-prepare 3D-RISM's iterative calculation, for Variable Cell
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: at_old(3, 3)
    !
    INTEGER :: ierr
    LOGICAL :: laue
    !
    IF (.NOT. lrism3d) THEN
      RETURN
    END IF
    !
    CALL start_clock('3DRISM_pre2')
    !
    IF (rism3t%itype == ITYPE_3DRISM) THEN
      laue = .FALSE.
    ELSE !IF (rism3t%itype == ITYPE_LAUERISM) THEN
      laue = .TRUE.
    END IF
    !
    ! ... scale G-vectors
    CALL scale_fft_3drism(rism3t%cfft, at_old, laue)
    IF (laue) THEN
      CALL scale_fft_lauerism(rism3t%lfft, at_old)
    END IF
    !
    ! ... allocate solvent's susceptibility
    IF (.NOT. laue) THEN
      CALL refresh_suscept_3drism(rism3t)
    ELSE
      CALL refresh_suscept_lauerism(rism3t, both_hands)
    END IF
    !
    ! ... calculate new solvent's susceptibility
    CALL init_3drism(rism3t, .FALSE., both_hands, ierr)
    !
    CALL stop_clock('3DRISM_pre2')
    !
  END SUBROUTINE rism3d_reprepare
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism3d_run(vrs, rhog, lconv, epsv_curr)
    !----------------------------------------------------------------------------
    !
    ! ... perform 3D-RISM's iterative calculation
    ! ...
    ! ...   vrs:       coulomb potential in R-space.
    ! ...   rhog:      electronic density in G-space, which is used only if Laue-RISM.
    ! ...   epsv_curr: convergence threshold
    !
    IMPLICIT NONE
    !
    REAL(DP),           INTENT(IN)  :: vrs(:)
    COMPLEX(DP),        INTENT(IN)  :: rhog(:)
    LOGICAL,            INTENT(OUT) :: lconv
    REAL(DP), OPTIONAL, INTENT(IN)  :: epsv_curr
    !
    INTEGER  :: ierr
    REAL(DP) :: epsv_
    REAL(DP) :: charge
    !
    IF (.NOT. lrism3d) THEN
      lconv = .FALSE.
      RETURN
    END IF
    !
    CALL start_clock('3DRISM_run')
    !
    ! ... check epsv_curr
    IF (conv_level >= 2) THEN
      ! high level
      epsv_ = epsv
    ELSE IF (PRESENT(epsv_curr)) THEN
      IF (conv_level == 1 .AND. epsv_curr > 0.0_DP) THEN
        ! medium level
        epsv_ = MAX(epsv, epsv_curr * SQRT(epsv_curr))
      ELSE
        ! low level
        epsv_ = MAX(epsv, epsv_curr)
      END IF
    ELSE
      epsv_ = epsv
    END IF
    !
    ! ... set DFT's potential
    CALL potential_3drism(rism3t, vrs, rhog, ierr)
    !
    IF (ierr /= IERR_RISM_NULL) THEN
      lconv = .FALSE.
      CALL stop_by_err_rism('rism3d_run', ierr)
    END IF
    !
    ! ... set initial guess
    IF (.NOT. has_any_corr) THEN
      CALL guess_3drism(rism3t, ierr)
      !
      IF (ierr /= IERR_RISM_NULL) THEN
        lconv = .FALSE.
        CALL stop_by_err_rism('rism3d_run', ierr)
      END IF
    END IF
    !
    ! ... calculate 3D-RISM
    IF (rism3t%itype == ITYPE_3DRISM) THEN
      CALL do_3drism(rism3t, niter, epsv_, mdiis_size, mdiis_step, '', ierr)
      !
    ELSE !IF (rism3t%itype == ITYPE_LAUERISM) THEN
      CALL charge_esm(rhog, charge)
      qsol = -charge
      !
      IF (ABS(qsol) > eps8) THEN
        CALL check_solvent_has_ions()
      END IF
      !
      CALL do_lauerism(rism3t, niter, epsv_, mdiis_size, mdiis_step, &
                     & qsol, both_hands, ireference, '', ierr)
    END IF
    !
    ! ... 3D-RISM has been converged ?
    IF (ierr == IERR_RISM_NOT_CONVERGED) THEN
      lconv = .FALSE.
      !
    ELSE IF (ierr == IERR_RISM_NULL) THEN
      lconv = .TRUE.
      !
    ELSE ! an error has been occurred
      lconv = .FALSE.
      CALL stop_by_err_rism('rism3d_run', ierr)
    END IF
    !
    ! ... here, correlation is nonzero.
    has_any_corr = .TRUE.
    !
    CALL stop_clock('3DRISM_run')
    !
  END SUBROUTINE rism3d_run
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism3d_potential(vrs, rhog)
    !----------------------------------------------------------------------------
    !
    ! ... calculate 3D-RISM's potential (both solute and solvent)
    ! ...
    ! ...   vrs:  electronic coulomb potential in R-space.
    ! ...   rhog: electronic density in G-space, which is used only if Laue-RISM.
    !
    IMPLICIT NONE
    !
    REAL(DP),    INTENT(IN)  :: vrs(:)
    COMPLEX(DP), INTENT(IN)  :: rhog(:)
    !
    INTEGER :: ierr
    !
    IF (.NOT. lrism3d) THEN
      RETURN
    END IF
    !
    CALL start_clock('3DRISM_pot')
    !
    ! ... potential from solute
    CALL potential_3drism(rism3t, vrs, rhog, ierr)
    !
    IF (ierr /= IERR_RISM_NULL) THEN
      CALL stop_by_err_rism('rism3d_potential', ierr)
    END IF
    !
    ! ... potential from solvent
    IF (rism3t%itype == ITYPE_3DRISM) THEN
      CALL solvation_3drism(rism3t, ierr)
    ELSE !IF (rism3t%itype == ITYPE_LAUERISM) THEN
      CALL solvation_lauerism(rism3t, qsol, ireference, ierr)
    END IF
    !
    IF (ierr /= IERR_RISM_NULL) THEN
      CALL stop_by_err_rism('rism3d_potential', ierr)
    END IF
    !
    CALL stop_clock('3DRISM_pot')
    !
  END SUBROUTINE rism3d_potential
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism3d_force(force, vloc)
    !----------------------------------------------------------------------------
    !
    ! ... calculate 3D-RISM's force
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT) :: force(3, nat)
    REAL(DP), INTENT(IN)  :: vloc(ngl, nsp)
    !
    INTEGER :: ierr
    !
    IF (.NOT. lrism3d) THEN
      RETURN
    END IF
    !
    CALL start_clock('3DRISM_for')
    !
    CALL solvation_force(rism3t, force, vloc, ierr)
    !
    IF (ierr /= IERR_RISM_NULL) THEN
      CALL stop_by_err_rism('rism3d_force', ierr)
    END IF
    !
    CALL stop_clock('3DRISM_for')
    !
  END SUBROUTINE rism3d_force
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism3d_stress(sigma)
    !----------------------------------------------------------------------------
    !
    ! ... calculate 3D-RISM's stress
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT) :: sigma(3, 3)
    !
    INTEGER :: ierr
    !
    IF (.NOT. lrism3d) THEN
      RETURN
    END IF
    !
    CALL start_clock('3DRISM_str')
    !
    CALL solvation_stress(rism3t, sigma, ierr)
    !
    IF (ierr /= IERR_RISM_NULL) THEN
      CALL stop_by_err_rism('rism3d_stress', ierr)
    END IF
    !
    CALL stop_clock('3DRISM_str')
    !
  END SUBROUTINE rism3d_stress
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism3d_write_to_restart(ext)
    !----------------------------------------------------------------------------
    !
    ! ... write 3D-RISM's data (for restart calculation)
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: ext
    !
    IF (.NOT. lrism3d) THEN
      RETURN
    END IF
    !
    IF (PRESENT(ext)) THEN
      CALL write_3drism(rism3t, ecutsolv, gamma_only, ext)
    ELSE
      CALL write_3drism(rism3t, ecutsolv, gamma_only)
    END IF
    !
  END SUBROUTINE rism3d_write_to_restart
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism3d_read_to_restart(ext)
    !----------------------------------------------------------------------------
    !
    ! ... read 3D-RISM's data (for restart calculation)
    !
    IMPLICIT NONE
    !
    CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: ext
    !
    IF (.NOT. lrism3d) THEN
      RETURN
    END IF
    !
    IF (PRESENT(ext)) THEN
      CALL read_3drism(rism3t, ecutsolv, ext)
    ELSE
      CALL read_3drism(rism3t, ecutsolv)
    END IF
    !
  END SUBROUTINE rism3d_read_to_restart
  !
  !----------------------------------------------------------------------------
  LOGICAL FUNCTION rism3d_is_laue()
    !----------------------------------------------------------------------------
    !
    ! ... Laue-RISM or not.
    !
    IMPLICIT NONE
    !
    IF (.NOT. lrism3d) THEN
      rism3d_is_laue = .FALSE.
      RETURN
    END IF
    !
    IF (rism3t%itype == ITYPE_LAUERISM) THEN
      rism3d_is_laue = .TRUE.
    ELSE
      rism3d_is_laue = .FALSE.
    END IF
    !
  END FUNCTION rism3d_is_laue
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism3d_set_laue()
    !----------------------------------------------------------------------------
    !
    ! ... to be Laue-RISM.
    !
    IMPLICIT NONE
    !
    rism3t%itype = ITYPE_LAUERISM
    !
  END SUBROUTINE rism3d_set_laue
  !
  !----------------------------------------------------------------------------
  LOGICAL FUNCTION rism3d_is_both_hands()
    !----------------------------------------------------------------------------
    !
    ! ... Laue-RISM with both-hands or not.
    !
    IMPLICIT NONE
    !
    IF (.NOT. lrism3d) THEN
      rism3d_is_both_hands = .FALSE.
      RETURN
    END IF
    !
    IF (rism3t%itype /= ITYPE_LAUERISM) THEN
      rism3d_is_both_hands = .FALSE.
      RETURN
    END IF
    !
    IF (rism3t%lfft%xright .AND. rism3t%lfft%xleft) THEN
      rism3d_is_both_hands = .TRUE.
    ELSE
      rism3d_is_both_hands = .FALSE.
    END IF
    !
  END FUNCTION rism3d_is_both_hands
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism3d_printpot()
    !----------------------------------------------------------------------------
    !
    ! ... print planar averaged data.
    !
    IMPLICIT NONE
    !
    INTEGER :: ierr
    !
    CHARACTER(LEN=5), PARAMETER :: ext = 'rism1'
    !
    IF (.NOT. lrism3d) THEN
      RETURN
    END IF
    !
    ierr = IERR_RISM_NULL
    IF (rism3t%itype == ITYPE_LAUERISM .OR. planar_average) THEN
      CALL print_solvavg(rism3t, ext, ierr)
    END IF
    !
    IF (ierr /= IERR_RISM_NULL) THEN
      CALL stop_by_err_rism('rism3d_printpot', ierr)
    END IF
    !
  END SUBROUTINE rism3d_printpot
  !
  !----------------------------------------------------------------------------
  SUBROUTINE rism3d_print_clock()
    !----------------------------------------------------------------------------
    !
    ! ... print clock for 3D-RISM
    !
    IMPLICIT NONE
    !
    IF (.NOT. lrism3d) THEN
      RETURN
    END IF
    !
    CALL print_clock('3DRISM_pre')
    CALL print_clock('3DRISM_pre2')
    CALL print_clock('3DRISM_run')
    CALL print_clock('3DRISM_pot')
    CALL print_clock('3DRISM_ions')
    CALL print_clock('3DRISM_for')
    CALL print_clock('3DRISM_str')
    !
  END SUBROUTINE rism3d_print_clock
  !
  !----------------------------------------------------------------------------
  SUBROUTINE check_solvent_is_neutral()
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER  :: iq
    INTEGER  :: iiq
    INTEGER  :: iv
    INTEGER  :: nv
    INTEGER  :: isolV
    INTEGER  :: iatom
    REAL(DP) :: rhotot1
    REAL(DP) :: rhotot2
    REAL(DP) :: rhov1
    REAL(DP) :: rhov2
    REAL(DP) :: qv
    !
    rhotot1 = 0.0_DP
    rhotot2 = 0.0_DP
    !
    DO iq = rism3t%mp_site%isite_start, rism3t%mp_site%isite_end
      iiq     = iq - rism3t%mp_site%isite_start + 1
      iv      = iuniq_to_isite(1, iq)
      nv      = iuniq_to_nsite(iq)
      isolV   = isite_to_isolV(iv)
      iatom   = isite_to_iatom(iv)
      rhov1   = DBLE(nv) * solVs(isolV)%density
      rhov2   = DBLE(nv) * solVs(isolV)%subdensity
      qv      = solVs(isolV)%charge(iatom)
      rhotot1 = rhotot1 + qv * rhov1
      rhotot2 = rhotot2 + qv * rhov2
    END DO
    !
    CALL mp_sum(rhotot1, rism3t%mp_site%inter_sitg_comm)
    CALL mp_sum(rhotot2, rism3t%mp_site%inter_sitg_comm)
    !
    IF (ABS(rhotot1) > eps12 .OR. ABS(rhotot2) > eps12) THEN
      CALL stop_by_err_rism('rism3d_initialize', IERR_RISM_NONZERO_CHARGE)
    END IF
    !
  END SUBROUTINE check_solvent_is_neutral
  !
  !----------------------------------------------------------------------------
  SUBROUTINE check_solvent_has_ions()
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER               :: iq
    INTEGER               :: iiq
    INTEGER               :: iv
    INTEGER               :: nv
    INTEGER               :: isolV
    INTEGER               :: iatom
    REAL(DP)              :: qv
    REAL(DP), ALLOCATABLE :: qmol(:)
    !
    ALLOCATE(qmol(nsolV))
    qmol = 0.0_DP
    !
    DO iq = rism3t%mp_site%isite_start, rism3t%mp_site%isite_end
      iiq     = iq - rism3t%mp_site%isite_start + 1
      iv      = iuniq_to_isite(1, iq)
      nv      = iuniq_to_nsite(iq)
      isolV   = isite_to_isolV(iv)
      iatom   = isite_to_iatom(iv)
      qv      = solVs(isolV)%charge(iatom)
      qmol(isolV) = qmol(isolV) + DBLE(nv) * qv
    END DO
    !
    CALL mp_sum(qmol, rism3t%mp_site%inter_sitg_comm)
    !
    IF (.NOT. ANY(ABS(qmol(:)) > eps12)) THEN
      CALL stop_by_err_rism('rism3d_initialize', IERR_RISM_NOT_ANY_IONS)
    END IF
    !
    DEALLOCATE(qmol)
    !
  END SUBROUTINE check_solvent_has_ions
  !
END MODULE rism3d_facade

