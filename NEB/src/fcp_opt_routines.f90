!
! Copyright (C) 2003-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE fcp_opt_routines
   !---------------------------------------------------------------------------
   !
   ! ... This module contains all subroutines and functions needed for
   ! ... the optimisation of the FCPs.
   !
   ! ... Written by Carlo Sbraccia ( 2003-2006 )
   !
   ! ... Newton algorithm is implemented by S. Nishihara ( 2016-2017 )
   !
   USE kinds,          ONLY : DP
   USE io_global,      ONLY : meta_ionode, meta_ionode_id
   USE mp,             ONLY : mp_bcast
   USE mp_world,       ONLY : world_comm
   USE fcp_module,     ONLY : fcp_check
   USE fcp_variables,  ONLY : fcp_mu, lfcp_linmin, lfcp_newton, &
                              fcp_ndiis, tot_charge_first, tot_charge_last
   USE mdiis,          ONLY : mdiis_type, allocate_mdiis, deallocate_mdiis, update_by_mdiis
   USE path_variables, ONLY : num_of_images
   !
   IMPLICIT NONE
   !
   PRIVATE
   !
   REAL(DP), ALLOCATABLE :: fcp_neb_nelec(:)
   REAL(DP), ALLOCATABLE :: fcp_neb_ef(:)
   REAL(DP), ALLOCATABLE :: fcp_neb_dos(:)
   !
   ! ... variables for Line-Minimization
   REAL(DP), ALLOCATABLE :: force0(:)
   REAL(DP), ALLOCATABLE :: nelec0(:)
   LOGICAL,  ALLOCATABLE :: firstcall(:)
   !
   ! ... variables for DIIS (coupled with Newton-Raphson)
   LOGICAL          :: init_mdiis
   TYPE(mdiis_type) :: mdiist
   !
   PUBLIC :: fcp_neb_nelec, fcp_neb_ef, fcp_neb_dos
   PUBLIC :: fcp_opt_allocation, fcp_opt_deallocation, fcp_opt_perform
   !
CONTAINS
   !
   !----------------------------------------------------------------------
   SUBROUTINE fcp_opt_allocation()
      !----------------------------------------------------------------------
      !
      USE ions_base,      ONLY : nat, ityp, zv
      USE klist,          ONLY : tot_charge, nelec
#if ! defined (__XSD)
      USE pw_restart,     ONLY : pw_readfile
#else
#endif
      USE io_files,       ONLY : prefix, tmp_dir
      USE path_variables, ONLY : restart
      USE ener,           ONLY : ef
      !
      IMPLICIT NONE
      !
      REAL(DP) :: ionic_charge, nelec_, first, last
      INTEGER  :: n, i, ierr
      CHARACTER (LEN=256)   :: tmp_dir_saved
      !
      CHARACTER(LEN=6), EXTERNAL :: int_to_char
      !
      CALL fcp_check( .TRUE. )
      !
      ALLOCATE( fcp_neb_nelec( num_of_images ) )
      ALLOCATE( fcp_neb_ef   ( num_of_images ) )
      ALLOCATE( fcp_neb_dos  ( num_of_images ) )
      !
      IF ( lfcp_linmin ) THEN
         !
         ALLOCATE( force0    ( num_of_images ) )
         ALLOCATE( nelec0    ( num_of_images ) )
         ALLOCATE( firstcall ( num_of_images ) )
         !
         force0    (:) = 0.0_DP
         nelec0    (:) = 0.0_DP
         firstcall (:) = .TRUE.
         !
      ELSE IF ( lfcp_newton ) THEN
         !
         init_mdiis = .TRUE.
         CALL allocate_mdiis( mdiist, fcp_ndiis, num_of_images, 1.0_DP, 1 )
         !
      END IF
      !
      IF ( restart ) THEN
         !
         tmp_dir_saved = tmp_dir
         !
         DO i = 1, num_of_images
            !
            tmp_dir = TRIM( tmp_dir_saved ) // TRIM( prefix ) // "_" // &
                 TRIM( int_to_char( i ) ) // "/"
            !
#if ! defined (__XSD)
            CALL pw_readfile( 'reset', ierr )
            IF ( ierr /= 0 ) &
               CALL errore('fcp_opt_routines','error in pw_readfile(reset)',ABS(ierr))
            !
            CALL pw_readfile( 'nowave', ierr )
            IF ( ierr /= 0 ) &
               CALL errore('fcp_opt_routines','error in pw_readfile(nowave)',ABS(ierr))
            !
#else
            CALL errore('fcp_opt_routines','XSD implementation pending',1)
            !
#endif
            fcp_neb_nelec(i) = nelec
            fcp_neb_ef   (i) = ef
            CALL fcp_hessian( fcp_neb_dos(i) )
            !
         END DO
         !
         tmp_dir = tmp_dir_saved
         !
      ELSE
         !
         ionic_charge = SUM( zv(ityp(1:nat)) )
         nelec_ = ionic_charge - tot_charge
         !
         n = num_of_images
         first = tot_charge_first
         last  = tot_charge_last
         DO i = 1, n
            fcp_neb_nelec(i) = nelec_ - (first * (n - i) + last * (i - 1) ) / (n - 1)
         END DO
         !
         fcp_neb_ef (:) = 0.0_DP
         fcp_neb_dos(:) = 0.0_DP
         !
      END IF
      !
   END SUBROUTINE fcp_opt_allocation
   !
   !----------------------------------------------------------------------
   SUBROUTINE fcp_opt_deallocation()
      !----------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      IF ( ALLOCATED( fcp_neb_nelec ) ) DEALLOCATE( fcp_neb_nelec )
      IF ( ALLOCATED( fcp_neb_ef    ) ) DEALLOCATE( fcp_neb_ef    )
      IF ( ALLOCATED( fcp_neb_dos   ) ) DEALLOCATE( fcp_neb_dos   )
      IF ( ALLOCATED( force0        ) ) DEALLOCATE( force0        )
      IF ( ALLOCATED( nelec0        ) ) DEALLOCATE( nelec0        )
      IF ( ALLOCATED( firstcall     ) ) DEALLOCATE( firstcall     )
      !
      IF ( init_mdiis ) THEN
         !
         CALL deallocate_mdiis( mdiist )
         !
      END IF
      !
   END SUBROUTINE fcp_opt_deallocation
   !
   !----------------------------------------------------------------------
   SUBROUTINE fcp_opt_perform()
      !----------------------------------------------------------------------
      !
      USE fcp_module,    ONLY : solvation_radius
      USE fcp_variables, ONLY : solvation_radius_ => solvation_radius
      !
      IMPLICIT NONE
      !
      REAL(DP) :: step_max
      REAL(DP) :: capacitance
      !
      CALL fcp_check( .TRUE. )
      !
      ! ... evaluate maximum step
      !
      solvation_radius = solvation_radius_
      CALL fcp_capacitance( capacitance )
      !
      step_max = ABS( capacitance * 0.1_DP )
      !
      IF ( lfcp_linmin ) THEN
         !
         ! ... perform Line-Minimization
         !
         CALL fcp_line_minimization( step_max )
         !
      ELSE IF ( lfcp_newton ) THEN
         !
         ! ... perform Newton-Raphson
         !
         CALL fcp_newton( step_max )
         !
      END IF
      !
   END SUBROUTINE fcp_opt_perform
   !
   !----------------------------------------------------------------------
   SUBROUTINE fcp_line_minimization( step_max )
      !----------------------------------------------------------------------
      !
      USE constants,      ONLY : eps16
      USE path_variables, ONLY : frozen
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN) :: step_max
      !
      INTEGER  :: image
      REAL(DP) :: ef, dos, force, step
      REAL(DP) :: nelec, nelec_new
      !
      IF ( meta_ionode ) THEN
         !
         DO image = 1, num_of_images
            !
            IF ( frozen(image) ) CYCLE
            !
            nelec = fcp_neb_nelec(image)
            ef    = fcp_neb_ef   (image)
            dos   = fcp_neb_dos  (image)
            force = fcp_mu - ef
            !
            IF ( firstcall(image) ) THEN
               !
               ! ... initialization
               !
               firstcall(image) = .FALSE.
               !
               nelec0(image) = nelec
               force0(image) = force
               !
            END IF
            !
            IF ( ABS( force0(image) - force ) < eps16 ) THEN
               !
               ! ... Steepest-Descent
               !
               step = 0.0_DP
               CALL step_newton( dos, force, step )
               !
               nelec_new = nelec + step
               !
            ELSE
               !
               ! ... Line-Minimization
               !
               nelec_new = (nelec * force0(image) - nelec0(image) * force) &
                         & / (force0(image) - force)
               !
            END IF
            !
            ! ... save #electrons and force
            !
            nelec0(image) = nelec
            force0(image) = force
            !
            ! ... update #electrons
            !
            step = nelec_new - nelec
            step = MIN( step, +step_max )
            step = MAX( step, -step_max )
            !
            fcp_neb_nelec(image) = nelec + step
            !
         END DO
         !
      END IF
      !
      CALL mp_bcast( fcp_neb_nelec, meta_ionode_id, world_comm )
      !
      RETURN
      !
   END SUBROUTINE fcp_line_minimization
   !
   !----------------------------------------------------------------------
   SUBROUTINE fcp_newton( step_max )
      !----------------------------------------------------------------------
      !
      USE path_variables, ONLY : frozen
      !
      IMPLICIT NONE
      !
      REAL(DP), INTENT(IN) :: step_max
      !
      INTEGER               :: image
      REAL(DP)              :: ef, dos, force, step
      REAL(DP)              :: nelec, nelec_new
      REAL(DP), ALLOCATABLE :: nelec1(:)
      REAL(DP), ALLOCATABLE :: step1(:)
      !
      ALLOCATE(nelec1(num_of_images))
      ALLOCATE(step1(num_of_images))
      !
      IF ( meta_ionode ) THEN
         !
         DO image = 1, num_of_images
            !
            ! ... current #electrons and Newton's steps
            !
            nelec = fcp_neb_nelec(image)
            ef    = fcp_neb_ef   (image)
            dos   = fcp_neb_dos  (image)
            force = fcp_mu - ef
            !
            nelec1(image) = nelec
            CALL step_newton( dos, force, step1(image) )
            !
         END DO
         !
         ! ... apply DIIS
         !
         CALL update_by_mdiis( mdiist, nelec1, step1 )
         !
         DO image = 1, num_of_images
            !
            IF ( frozen(image) ) CYCLE
            !
            ! ... update #electrons
            !
            nelec     = fcp_neb_nelec(image)
            nelec_new = nelec1(image)
            !
            step = nelec_new - nelec
            step = MIN( step, +step_max )
            step = MAX( step, -step_max )
            !
            fcp_neb_nelec(image) = nelec + step
            !
         END DO
         !
      END IF
      !
      CALL mp_bcast( fcp_neb_nelec, meta_ionode_id, world_comm )
      !
      DEALLOCATE(nelec1)
      DEALLOCATE(step1)
      !
      RETURN
      !
   END SUBROUTINE fcp_newton
   !
   !----------------------------------------------------------------------
   SUBROUTINE step_newton( dos, force, step )
     !----------------------------------------------------------------------
     !
     USE constants,   ONLY : eps4
     USE cell_base,   ONLY : alat, at
     USE rism_module, ONLY : lrism
     !
     IMPLICIT NONE
     !
     REAL(DP), INTENT(IN)  :: dos
     REAL(DP), INTENT(IN)  :: force
     REAL(DP), INTENT(OUT) :: step
     !
     REAL(DP) :: area_xy, slope
     !
     IF ( ABS( dos ) > eps4 ) THEN
        !
        step = dos * force
        !
     ELSE
        !
        area_xy = alat * alat * ABS(at(1, 1) * at(2, 2) - at(1, 2) * at(2, 1))
        !
        IF (lrism) THEN
           !
           slope = 0.100_DP * area_xy
           !
        ELSE
           !
           slope = 0.001_DP * area_xy
           !
        END IF
        !
        step = slope * force
        !
     END IF
     !
   END SUBROUTINE step_newton
   !
END MODULE fcp_opt_routines
