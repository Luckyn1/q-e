!
! Copyright (C) 2001-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE fcp_relaxation
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
  ! ...   This module performes relaxation of FCP.
  ! ... . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
  !
  USE constants,     ONLY : eps16, RYTOEV
  USE control_flags, ONLY : iverbosity
  USE ener,          ONLY : ef
  USE io_global,     ONLY : stdout
  USE ions_base,     ONLY : nat, ityp, zv
  USE kinds,         ONLY : DP
  USE klist,         ONLY : nelec, tot_charge, nkstot, wk, degauss, ngauss
  USE mdiis,         ONLY : mdiis_type, allocate_mdiis, deallocate_mdiis, update_by_mdiis
  USE wvfct,         ONLY : nbnd, et
  !
  IMPLICIT NONE
  SAVE
  PRIVATE
  !
  ! ... define parameters
  INTEGER, PARAMETER :: IRELAX_NULL     = 0
  INTEGER, PARAMETER :: IRELAX_LINE_MIN = 1
  INTEGER, PARAMETER :: IRELAX_MDIIS    = 2
  INTEGER, PARAMETER :: IRELAX_NEWTON   = 3
  !
  ! ... define variables
  INTEGER          :: irelax         ! type of relaxation
  INTEGER          :: iter           ! number of iteration
  REAL(DP)         :: epsf           ! convergence threshold (in Ry)
  REAL(DP)         :: step_max       ! maximum step (in number of charge)
  REAL(DP)         :: nelec_old      ! old number of electrons
  REAL(DP)         :: force_old      ! old force acting on FCP
  REAL(DP)         :: line_min_step  ! step of Line-Minimization
  LOGICAL          :: line_min_init  ! init Line-Minimization, or not ?
  INTEGER          :: mdiis_size     ! size of MDIIS
  REAL(DP)         :: mdiis_step     ! step of MDIIS
  LOGICAL          :: mdiis_init     ! init MDIIS, or not ?
  REAL(DP)         :: newton_step    ! step of Newton-Raphson
  LOGICAL          :: newton_init    ! init Newton-Raphson, or not ?
  TYPE(mdiis_type) :: mdiist         ! data of MDIIS
  !
  ! ... public components
  PUBLIC :: fcprlx_init
  PUBLIC :: fcprlx_final
  PUBLIC :: fcprlx_prm_line_min
  PUBLIC :: fcprlx_prm_mdiis
  PUBLIC :: fcprlx_prm_newton
  PUBLIC :: fcprlx_set_line_min
  PUBLIC :: fcprlx_set_mdiis
  PUBLIC :: fcprlx_set_newton
  PUBLIC :: fcprlx_update
  !
CONTAINS
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcprlx_init()
    !----------------------------------------------------------------------------
    !
    ! ... initialize this module
    !
    IMPLICIT NONE
    !
    irelax        = IRELAX_NULL
    iter          = 0
    epsf          = 0.0_DP
    step_max      = 0.0_DP
    nelec_old     = 0.0_DP
    force_old     = 0.0_DP
    line_min_step = 0.5_DP
    mdiis_size    = 4
    mdiis_step    = 0.5_DP
    mdiis_init    = .FALSE.
    newton_step   = 0.5_DP
    newton_init   = .FALSE.
    !
  END SUBROUTINE fcprlx_init
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcprlx_final()
    !----------------------------------------------------------------------------
    !
    ! ... finalize this module
    !
    IMPLICIT NONE
    !
    IF (line_min_init) THEN
       !
       line_min_init = .FALSE.
       nelec_old     = 0.0_DP
       force_old     = 0.0_DP
       !
    END IF
    !
    IF (mdiis_init) THEN
       !
       mdiis_init = .FALSE.
       nelec_old  = 0.0_DP
       CALL deallocate_mdiis(mdiist)
       !
    END IF
    !
    IF (newton_init) THEN
       !
       newton_init = .FALSE.
       nelec_old   = 0.0_DP
       !
    END IF
    !
  END SUBROUTINE fcprlx_final
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcprlx_prm_line_min(step)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: step
    !
    IF (step > 0.0_DP) THEN
       line_min_step = step
    END IF
    !
  END SUBROUTINE fcprlx_prm_line_min
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcprlx_prm_mdiis(nsize, step)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    INTEGER,  INTENT(IN) :: nsize
    REAL(DP), INTENT(IN) :: step
    !
    IF (nsize > 0) THEN
       mdiis_size = nsize
    END IF
    !
    IF (step > 0.0_DP) THEN
       mdiis_step = step
    END IF
    !
  END SUBROUTINE fcprlx_prm_mdiis
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcprlx_prm_newton(step)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: step
    !
    IF (step > 0.0_DP) THEN
       newton_step = step
    END IF
    !
  END SUBROUTINE fcprlx_prm_newton
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcprlx_set_line_min(eps, smax)
    !----------------------------------------------------------------------------
    !
    ! ... set the type of relaxation to Line-Minimization
    ! ...
    ! ... Variables:
    ! ...   eps:  convergende threshold (in Ry/e)
    ! ...   smax: maximum step (in e)
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: eps
    REAL(DP), INTENT(IN) :: smax
    !
    irelax = IRELAX_LINE_MIN
    !
    CALL set_eps_smax(eps, smax, 'fcprlx_set_line_min')
    !
  END SUBROUTINE fcprlx_set_line_min
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcprlx_set_mdiis(eps, smax)
    !----------------------------------------------------------------------------
    !
    ! ... set the type of relaxation to MDIIS
    ! ...
    ! ... Variables:
    ! ...   eps:  convergende threshold (in Ry/e)
    ! ...   smax: maximum step (in e)
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: eps
    REAL(DP), INTENT(IN) :: smax
    !
    irelax = IRELAX_MDIIS
    !
    CALL set_eps_smax(eps, smax, 'fcprlx_set_mdiis')
    !
  END SUBROUTINE fcprlx_set_mdiis
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcprlx_set_newton(eps, smax)
    !----------------------------------------------------------------------------
    !
    ! ... set the type of relaxation to Newton-Raphson
    ! ...
    ! ... Variables:
    ! ...   eps:  convergende threshold (in Ry/e)
    ! ...   smax: maximum step (in e)
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: eps
    REAL(DP), INTENT(IN) :: smax
    !
    irelax = IRELAX_NEWTON
    !
    CALL set_eps_smax(eps, smax, 'fcprlx_set_newton')
    !
  END SUBROUTINE fcprlx_set_newton
  !
  !----------------------------------------------------------------------------
  SUBROUTINE set_eps_smax(eps, smax, sname)
    !----------------------------------------------------------------------------
    !
    IMPLICIT NONE
    !
    REAL(DP),         INTENT(IN) :: eps
    REAL(DP),         INTENT(IN) :: smax
    CHARACTER(LEN=*), INTENT(IN) :: sname
    !
    IF (eps < 0.0_DP) THEN
       !
       CALL errore(TRIM(sname), 'eps is negative', 1)
       !
    END IF
    !
    IF (smax <= 0.0_DP) THEN
       !
       CALL errore(TRIM(sname), 'smax is not positive', 1)
       !
    END IF
    !
    epsf     = eps
    step_max = smax
    !
  END SUBROUTINE set_eps_smax
  !
  !----------------------------------------------------------------------------
  SUBROUTINE fcprlx_update(mu, conv)
    !----------------------------------------------------------------------------
    !
    ! ... update number of electrons
    ! ...
    ! ... Variables:
    ! ...   mu:   target Fermi energy (in Ry)
    ! ...   conv: converged, or not ?
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN)    :: mu
    LOGICAL,  INTENT(INOUT) :: conv
    !
    REAL(DP) :: force
    REAL(DP) :: tot_charge_
    !
    ! ... count up
    !
    iter = iter + 1
    !
    ! ... set variablds
    !
    force = mu - ef
    !
    tot_charge_ = tot_charge
    !
    ! ... check if convergence for FCP minimization is achieved
    !
    conv = conv .AND. (ABS(force) < epsf)
    !
    IF (conv) THEN
       !
       ! ... converged
       !
       WRITE(stdout, '(/,5X,"FCP Relaxation: convergence achieved in " &
                      & ,I3," steps")' ) iter
       WRITE(stdout, '(/,5X,"End of FCP relaxation calculation")' )
       !
    ELSE
       !
       ! ... update nelec
       !
       IF (irelax == IRELAX_LINE_MIN) THEN
          !
          CALL line_minimization(force)
          !
       ELSE IF (irelax == IRELAX_MDIIS) THEN
          !
          CALL do_mdiis(force)
          !
       ELSE IF (irelax == IRELAX_NEWTON) THEN
          !
          CALL do_newton(force)
          !
       ELSE
          !
          CALL errore('fcprlx_update', 'irelax is incorrect', 1)
          !
       END IF
       !
       ! ... update tot_charge
       !
       tot_charge = SUM(zv(ityp(1:nat))) - nelec
       !
    END IF
    !
    ! ... write information
    !
    IF (.NOT. conv) THEN
       WRITE(stdout, '(/,5X,"FCP: iteration #",I5)') iter
       WRITE(stdout, '(  5X,"FCP: Total Charge = ",F12.6,"  -> ",F12.6)') tot_charge_, tot_charge
    ELSE
       WRITE(stdout, '(/,5X,"FCP: Total Charge = ",F12.6)') tot_charge_
    END IF
    !
    WRITE(stdout, '(5X,"FCP: Fermi Energy = ",F12.6," Ry (",F12.6," eV)")') ef,    ef    * RYTOEV
    WRITE(stdout, '(5X,"FCP: Target Level = ",F12.6," Ry (",F12.6," eV)")') mu,    mu    * RYTOEV
    WRITE(stdout, '(5X,"FCP: Force on FCP = ",F12.6," Ry (",F12.6," eV)")') force, force * RYTOEV
    WRITE(stdout, '(5X,"FCP: Force Thr.   = ",F12.6," Ry (",F12.6," eV)")') epsf,  epsf  * RYTOEV
    WRITE(stdout, '()')
    !
  END SUBROUTINE fcprlx_update
  !
  !----------------------------------------------------------------------------
  SUBROUTINE line_minimization(force)
    !----------------------------------------------------------------------------
    !
    ! ... perform one step of relaxation, using Line-Minimization algorithm.
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: force
    !
    REAL(DP) :: nelec0
    !
    ! ... initialize
    !
    IF (.NOT. line_min_init) THEN
       !
       line_min_init = .TRUE.
       !
       WRITE(stdout, '(/,5X,"FCP Relaxation Calculation")')
       WRITE(stdout, '(/,5X,"FCP: Line-Minimization Algorithm is used.")')
       WRITE(stdout, '(  5X,"FCP: Initial Step = ",F7.3," Ry^-1")') line_min_step
       !
       nelec_old = nelec
       force_old = force
       !
    END IF
    !
    ! ... perform Line-Minimization
    !
    IF (ABS(force_old - force) < eps16) THEN
       !
       nelec0 = nelec + line_min_step * force
       !
    ELSE
       !
       nelec0 = (nelec * force_old - nelec_old * force) / (force_old - force)
       !
    END IF
    !
    nelec_old = nelec
    force_old = force
    !
    ! ... update number of electrons
    !
    CALL update_nelec(nelec0)
    !
  END SUBROUTINE line_minimization
  !
  !----------------------------------------------------------------------------
  SUBROUTINE do_mdiis(force)
    !----------------------------------------------------------------------------
    !
    ! ... perform one step of relaxation, using MDIIS algorithm.
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: force
    !
    REAL(DP) :: nelec0
    REAL(DP) :: nelec1(1)
    REAL(DP) :: force1(1)
    !
    ! ... initialize
    !
    IF (.NOT. mdiis_init) THEN
       !
       mdiis_init = .TRUE.
       !
       WRITE(stdout, '(/,5X,"FCP Relaxation Calculation")')
       WRITE(stdout, '(/,5X,"FCP: MDIIS Algorithm is used.")')
       WRITE(stdout, '(  5X,"FCP: MDIIS Size   = ",I3           )') mdiis_size
       WRITE(stdout, '(  5X,"FCP: MDIIS Step   = ",F7.3," Ry^-1")') mdiis_step
       !
       CALL allocate_mdiis(mdiist, mdiis_size, 1, mdiis_step, 1)
       !
    END IF
    !
    ! ... save number of electrons
    !
    nelec_old = nelec
    !
    ! ... perform MDIIS
    !
    nelec1(1) = nelec
    force1(1) = force
    CALL update_by_mdiis(mdiist, nelec1, force1)
    nelec0 = nelec1(1)
    !
    ! ... update number of electrons
    !
    CALL update_nelec(nelec0)
    !
  END SUBROUTINE do_mdiis
  !
  !----------------------------------------------------------------------------
  SUBROUTINE do_newton(force)
    !----------------------------------------------------------------------------
    !
    ! ... perform one step of relaxation, using Newton-Raphson algorithm.
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: force
    !
    REAL(DP) :: nelec0
    REAL(DP) :: hess
    !
    ! ... initialize
    !
    IF (.NOT. newton_init) THEN
       !
       newton_init = .TRUE.
       !
       WRITE(stdout, '(/,5X,"FCP Relaxation Calculation")')
       WRITE(stdout, '(/,5X,"FCP: Newton-Raphson Algorithm is used.")')
       !
    END IF
    !
    ! ... save number of electrons
    !
    nelec_old = nelec
    !
    ! ... perform Newton-Raphson
    !
    CALL calc_hess(hess)
    !
    IF () THEN
       !
       nelec0 = nelec + force * hess
       !
    ELSE
       !
       nelec0 = nelec + force * hess
       !
    END IF
    !
    ! ... update number of electrons
    !
    CALL update_nelec(nelec0)
    !
  END SUBROUTINE do_newton
  !
  !----------------------------------------------------------------------------
  SUBROUTINE update_nelec(nelec0)
    !----------------------------------------------------------------------------
    !
    ! ... update nelec within restriction
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(IN) :: nelec0
    !
    REAL(DP) :: step
    REAL(DP) :: ionic_charge
    !
    step  = nelec0 - nelec_old
    step  = MIN(step, +step_max)
    step  = MAX(step, -step_max)
    nelec = nelec_old + step
    !
    IF (iverbosity > 0) THEN
       !
       ionic_charge = SUM(zv(ityp(1:nat)))
       !
       WRITE(stdout,'(5X,"FCP: Original charge = ",F12.6)') ionic_charge - nelec_old
       WRITE(stdout,'(5X,"FCP: Expected charge = ",F12.6)') ionic_charge - nelec0
       WRITE(stdout,'(5X,"FCP: Next charge     = ",F12.6)') ionic_charge - nelec
       !
    END IF
    !
  END SUBROUTINE update_nelec
  !
  !----------------------------------------------------------------------------
  SUBROUTINE calc_hess(hess)
    !----------------------------------------------------------------------------
    !
    ! ... calculate Hessian:
    ! ...
    ! ...     d^2E/dN^2 = d(ef)/dN = 1/DOS(ef)
    ! ...
    !
    IMPLICIT NONE
    !
    REAL(DP), INTENT(OUT) :: hess
    !
    INTEGER :: ik
    INTEGER :: ibnd
    !
    REAL(DP), EXTERNAL :: w0gauss
    !
    hess = 0.0_DP
    !
    DO ik = 1, nkstot
       !
       DO ibnd = 1, nbnd
          !
          hess = hess + wk (ik) * &
               & degauss * w0gauss((ef - et(ibnd, ik)) / degauss, ngauss)
          !
       END DO
       !
    END DO
    !
  END SUBROUTINE calc_hess
  !
END MODULE fcp_relaxation
