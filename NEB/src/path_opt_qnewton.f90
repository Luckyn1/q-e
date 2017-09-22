!
! Copyright (C) 2017 Satomichi Nishihara
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------------
MODULE path_opt_qnewton
  !---------------------------------------------------------------------------
  !
  ! ... This module provides quasi-newton methods to optimize the reaction path.
  ! ... Two algorithms are implemented:
  ! ...   Limited memory BFGS (L-BFGS)
  ! ...   Limited memory SR1  (L-SR1)
  ! ...
  ! ... (R.H.Byrd et al., Math. Program., 63:129-156, 1994.)
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : eps8, eps16
  USE path_io_units_module, ONLY : qnew_file, iunqnew, iunpath
  USE path_variables,       ONLY : ds, pos, grad, dim1, frozen, nim => num_of_images, &
                                   qnewton_ndim, qnewton_step
  USE fcp_variables,        ONLY : fcp_mu, fcp_nelec, fcp_ef
  USE fcp_opt_routines,     ONLY : fcp_opt_scale
  USE io_global,            ONLY : meta_ionode, meta_ionode_id
  USE mp,                   ONLY : mp_bcast
  USE mp_world,             ONLY : world_comm
  !
  USE basic_algebra_routines
  !
  IMPLICIT NONE
  !
  PRIVATE
  !
  INTEGER, PARAMETER ::     &
            HESS_LBFGS = 1, &
            HESS_LSR1  = 2
  !
  PUBLIC :: qnewton_lbfgs, &
            qnewton_lsr1
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE qnewton_lbfgs( lfcp )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       LOGICAL, INTENT(IN) :: lfcp
       !
       CALL qnewton_x( HESS_LBFGS, lfcp )
       !
     END SUBROUTINE qnewton_lbfgs
     !
     !-----------------------------------------------------------------------
     SUBROUTINE qnewton_lsr1( lfcp )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       LOGICAL, INTENT(IN) :: lfcp
       !
       CALL qnewton_x( HESS_LSR1, lfcp )
       !
     END SUBROUTINE qnewton_lsr1
     !
     !-----------------------------------------------------------------------
     SUBROUTINE qnewton_x( ihess, lfcp )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       INTEGER,  INTENT(IN)  :: ihess
       LOGICAL,  INTENT(IN)  :: lfcp
       !
       REAL(DP), ALLOCATABLE :: mask(:)
       REAL(DP), ALLOCATABLE :: dx(:)
       REAL(DP), ALLOCATABLE :: x1(:), g1(:)
       REAL(DP), ALLOCATABLE :: x2(:), g2(:)
       INTEGER,  ALLOCATABLE :: map(:)
       REAL(DP), ALLOCATABLE :: s(:,:), y(:,:)
       !
       REAL(DP)              :: norm_g, norm_dx
       REAL(DP)              :: xscale
       LOGICAL               :: exists, opened
       INTEGER               :: dim2
       INTEGER               :: map1, nsave
       INTEGER               :: i, I_in, I_fin
       !
       ! ... set dimension an image
       !
       IF ( lfcp ) THEN
          dim2 = dim1 + 1
       ELSE
          dim2 = dim1
       END IF
       !
       ! ... scaling factor: #electron -> Bohr
       !
       IF ( lfcp ) THEN
          xscale = fcp_opt_scale()
          IF ( xscale < eps8 ) THEN
             xscale = 1.0_DP
          END IF
       END IF
       !
       ALLOCATE( mask( dim2*nim ) )
       ALLOCATE( dx( dim2*nim ) )
       ALLOCATE( x1( dim2*nim ), g1( dim2*nim ) )
       ALLOCATE( x2( dim2*nim ), g2( dim2*nim ) )
       ALLOCATE( map( qnewton_ndim ) )
       ALLOCATE( s( dim2*nim, qnewton_ndim ) )
       ALLOCATE( y( dim2*nim, qnewton_ndim ) )
       !
       ! ... mask at frozen images
       !
       mask(:) = 0.0_DP
       DO i = 1, nim
          I_in  = ( i - 1 )*dim2 + 1
          I_fin = i * dim2
          IF ( frozen(i) ) CYCLE
          mask(I_in:I_fin) = 1.0_DP
       END DO
       !
       ! ... current position and gradient
       !
       DO i = 1, nim
          I_in  = ( i - 1 )*dim2 + 1
          I_fin = ( i - 1 )*dim2 + dim1
          !
          x1(I_in:I_fin) = pos (:,i)
          g1(I_in:I_fin) = grad(:,i)
          !
          IF ( lfcp ) THEN
             x1(I_fin+1) = fcp_nelec(i) * xscale
             g1(I_fin+1) = (fcp_ef(i) - fcp_mu) / xscale
          END IF
       END DO
       !
       ! ... check norm of gradient
       !
       norm_g = MAXVAL( ABS( mask(:) * g1(:) ) )
       !
       IF ( meta_ionode .AND. norm_g > 0.0_DP ) THEN
          !
          ! ... read the history from the file
          !
          INQUIRE( FILE = qnew_file, EXIST = exists )
          !
          IF ( exists ) THEN
             !
             OPEN( UNIT = iunqnew, FILE = qnew_file, STATUS = "OLD" )
             !
             READ( UNIT = iunqnew, FMT = * ) i
             !
             exists = ( i == nim )
             opened = .TRUE.
             !
          END IF
          !
          IF ( exists ) THEN
             !
             READ( UNIT = iunqnew, FMT = * ) nsave
             READ( UNIT = iunqnew, FMT = * ) x2(:)
             READ( UNIT = iunqnew, FMT = * ) g2(:)
             IF ( 1 <= nsave .AND. nsave <= qnewton_ndim ) THEN
                READ( UNIT = iunqnew, FMT = * ) map(1:nsave)
                READ( UNIT = iunqnew, FMT = * ) s(:,1:nsave)
                READ( UNIT = iunqnew, FMT = * ) y(:,1:nsave)
             ELSE
                nsave  = 0
                map(:) = 0
                s(:,:) = 0.0_DP
                y(:,:) = 0.0_DP
             END IF
             !
          ELSE
             !
             nsave  = -1
             x2(:)  = 0.0_DP
             g2(:)  = 0.0_DP
             map(:) = 0
             s(:,:) = 0.0_DP
             y(:,:) = 0.0_DP
             !
          END IF
          !
          IF ( opened ) CLOSE( UNIT = iunqnew )
          !
          ! ... add current data
          !
          IF ( nsave < qnewton_ndim ) THEN
             nsave = MAX( nsave + 1, 0 )
             map1  = nsave
          ELSE
             nsave = qnewton_ndim
             map1  = map(nsave)
          END IF
          !
          DO i = nsave, 2, -1
             map(i) = map(i-1)
          END DO
          IF ( nsave > 0 ) THEN
             map(1) = map1
          END IF
          !
          DO i = 1, nsave
             IF ( map(i) < 1 .OR. qnewton_ndim < map(i) ) THEN
                map(:) = 0
                map(1) = 1
                EXIT
             END IF
          END DO
          !
          IF ( nsave > 0 ) THEN
             s(:,map(1)) = x1(:) - x2(:)
             y(:,map(1)) = g1(:) - g2(:)
          END IF
          !
          ! ... act hessian
          !
          IF ( ihess == HESS_LBFGS ) THEN
             !
             CALL lbfgs_hess( dim2*nim, nsave, dx, g1, map, s, y )
             !
          ELSE IF ( ihess == HESS_LSR1 ) THEN
             !
             CALL lsr1_hess( dim2*nim, nsave, dx, g1, map, s, y )
             !
          END IF
          !
          IF ( nsave < 0 ) THEN
             !
             nsave = 0
             dx(:) = -ds * ds * g1(:)
             !
             WRITE( UNIT = iunpath, &
                    FMT = '(/,5X,"hessian is not well-defined : history is reset",/)' )
             !
          END IF
          !
          ! ... write the history to the file
          !
          OPEN( UNIT = iunqnew, FILE = qnew_file )
          !
          WRITE( UNIT = iunqnew, FMT = * ) nim
          WRITE( UNIT = iunqnew, FMT = * ) nsave
          WRITE( UNIT = iunqnew, FMT = * ) x1(:)
          WRITE( UNIT = iunqnew, FMT = * ) g1(:)
          IF ( nsave > 0 ) THEN
             WRITE( UNIT = iunqnew, FMT = * ) map(1:nsave)
             WRITE( UNIT = iunqnew, FMT = * ) s(:,1:nsave)
             WRITE( UNIT = iunqnew, FMT = * ) y(:,1:nsave)
          END IF
          !
          CLOSE( UNIT = iunqnew )
          !
          ! ... update position
          !
          dx(:) = mask(:) * dx(:)
          !
          norm_dx = norm( dx(:) )
          !
          dx(:) = dx(:) / norm_dx * MIN( norm_dx, qnewton_step )
          !
          DO i = 1, nim
             I_in  = ( i - 1 )*dim2 + 1
             I_fin = ( i - 1 )*dim2 + dim1
             !
             pos (:,i) = pos (:,i) + dx(I_in:I_fin)
             !
             IF ( lfcp ) THEN
                fcp_nelec(i) = fcp_nelec(i) + dx(I_fin+1) / xscale
             END IF
          END DO
          !
       END IF
       !
       CALL mp_bcast( pos, meta_ionode_id, world_comm )
       !
       IF ( lfcp ) THEN
          CALL mp_bcast( fcp_nelec, meta_ionode_id, world_comm )
       END IF
       !
       DEALLOCATE( mask )
       DEALLOCATE( dx )
       DEALLOCATE( x1, g1 )
       DEALLOCATE( x2, g2 )
       DEALLOCATE( map )
       DEALLOCATE( s )
       DEALLOCATE( y )
       !
       RETURN
       !
     END SUBROUTINE qnewton_x
     !
     !-----------------------------------------------------------------------
     SUBROUTINE lbfgs_hess( ndim, nvec, dx, g, map, s, y )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       INTEGER,  INTENT(IN)    :: ndim
       INTEGER,  INTENT(INOUT) :: nvec
       REAL(DP), INTENT(OUT)   :: dx(ndim)
       REAL(DP), INTENT(IN)    :: g(ndim)
       INTEGER,  INTENT(IN)    :: map(*)
       REAL(DP), INTENT(IN)    :: s(ndim,*)
       REAL(DP), INTENT(IN)    :: y(ndim,*)
       !
       REAL(DP), ALLOCATABLE   :: rho(:), alpha(:)
       REAL(DP)                :: H0, ys, beta
       INTEGER                 :: i, ii
       !
       H0 = ds * ds
       !
       IF ( nvec < 1 ) THEN
          !
          dx(:) = -H0 * g(:)
          !
          RETURN
          !
       END IF
       !
       ALLOCATE( rho( nvec ), alpha( nvec ) )
       !
       dx(:) = g(:)
       !
       DO i = 1, nvec
          !
          ii = map(i)
          !
          ys = ( y(:,ii) .dot. s(:,ii) )
          !
          IF ( ABS( ys ) < eps16 ) THEN
             nvec = -1
             GOTO 100
          END IF
          !
          rho(i)   = 1.0_DP / ys
          alpha(i) = rho(i) * ( s(:,ii) .dot. dx(:) )
          dx(:)    = dx(:) - alpha(i) * y(:,ii)
          !
       END DO
       !
       dx(:) = H0 * dx(:)
       !
       DO i = nvec, 1, -1
          !
          ii = map(i)
          !
          beta  = rho(i) * ( y(:,ii) .dot. dx(:) )
          dx(:) = dx(:) + (alpha(i) - beta) * s(:,ii)
          !
       END DO
       !
       dx(:) = -dx(:)
       !
100    CONTINUE
       DEALLOCATE( rho, alpha )
       !
     END SUBROUTINE lbfgs_hess
     !
     !-----------------------------------------------------------------------
     SUBROUTINE lsr1_hess( ndim, nvec, dx, g, map, s, y )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       INTEGER,  INTENT(IN)    :: ndim
       INTEGER,  INTENT(INOUT) :: nvec
       REAL(DP), INTENT(OUT)   :: dx(ndim)
       REAL(DP), INTENT(IN)    :: g(ndim)
       INTEGER,  INTENT(IN)    :: map(*)
       REAL(DP), INTENT(IN)    :: s(ndim,*)
       REAL(DP), INTENT(IN)    :: y(ndim,*)
       !
       REAL(DP), ALLOCATABLE   :: M(:,:)
       REAL(DP), ALLOCATABLE   :: psi(:,:)
       REAL(DP), ALLOCATABLE   :: psig(:)
       REAL(DP)                :: H0
       INTEGER                 :: i, j, ii
       INTEGER                 :: info
       INTEGER,  ALLOCATABLE   :: ipiv(:)
       INTEGER                 :: nwork
       REAL(DP), ALLOCATABLE   :: work(:)
       !
       H0 = ds * ds
       !
       dx(:) = -H0 * g(:)
       !
       IF ( nvec < 1 ) RETURN
       !
       nwork = 64 * nvec
       !
       ALLOCATE( M( nvec, nvec ) )
       ALLOCATE( psi( ndim, nvec ) )
       ALLOCATE( psig( nvec ) )
       ALLOCATE( ipiv( nvec ) )
       ALLOCATE( work( nwork ) )
       !
       DO i = 1, nvec
          !
          ii = map(nvec - i + 1)
          !
          psi(:,i) = s(:,ii) - H0 * y(:,ii)
          !
       END DO
       !
       DO i = 1, nvec
          !
          ii = map(nvec - i + 1)
          !
          DO j = 1, i
             !
             M(i,j) = ( y(:,ii) .dot. psi(:,j) )
             M(j,i) = M(i,j)
             !
          END DO
          !
       END DO
       !
       CALL DGEMV( 'T', ndim, nvec, 1.0_DP, psi, ndim, g, 1, 0.0_DP, psig, 1 )
       !
       CALL DSYSV( 'U', nvec, 1, M, nvec, ipiv, psig, nvec, work, nwork, info )
       !
       IF ( info == 0 ) THEN
          CALL DGEMV( 'N', ndim, nvec, -1.0_DP, psi, ndim, psig, 1, 1.0_DP, dx, 1 )
       ELSE
          nvec = -1
       END IF
       !
       DEALLOCATE( M )
       DEALLOCATE( psi )
       DEALLOCATE( psig )
       DEALLOCATE( ipiv )
       DEALLOCATE( work )
       !
     END SUBROUTINE lsr1_hess
     !
END MODULE path_opt_qnewton
