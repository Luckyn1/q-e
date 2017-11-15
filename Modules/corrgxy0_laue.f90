!
! Copyright (C) 2017 National Institute of Advanced Industrial Science and Technology (AIST)
! [ This code is written by Satomichi Nishihara. ]
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE corrgxy0_laue(rismt, lextract, ar, ag0, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... sum or extract Gxy=0 terms of correlations, in Laue-RISM calculation.
  ! ...
  ! ... Variables:
  ! ...   lextract: if .TRUE.  extract Gxy=0 terms of correlations
  ! ...             if .FALSE. sum A(r) + A(gxy=0,z)
  !
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE rism,      ONLY : rism_type, ITYPE_LAUERISM
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(IN)    :: rismt
  LOGICAL,         INTENT(IN)    :: lextract
  REAL(DP),        INTENT(INOUT) :: ar(:, :)
  REAL(DP),        INTENT(INOUT) :: ag0(:, :)
  INTEGER,         INTENT(OUT)   :: ierr
  !
  ! ... check data type
  IF (rismt%itype /= ITYPE_LAUERISM) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nr < rismt%cfft%dfftt%nnr) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (rismt%nrzl < rismt%lfft%nrz) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  IF (.NOT. lextract) THEN
    !
    ! ... extract Gxy=0 terms
    CALL extract_gxy0()
    !
  ELSE
    !
    ! ... sum A(r) + A(gxy=0,z)
    CALL sum_gxy0()
    !
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
CONTAINS
  !
  SUBROUTINE extract_gxy0()
    IMPLICIT NONE
    INTEGER               :: ir
    INTEGER               :: idx
    INTEGER               :: idx0
    INTEGER               :: i3min
    INTEGER               :: i3max
    INTEGER               :: i1, i2, i3
    INTEGER               :: iz
    INTEGER               :: iiz
    INTEGER               :: isite
    REAL(DP), ALLOCATABLE :: bg0(:,:)
#if defined(__OPENMP)
    REAL(DP), ALLOCATABLE :: bg1(:,:)
#endif
    !
    IF (rismt%nsite < 1) THEN
      RETURN
    END IF
    !
    idx0 = rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x &
       & * rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
    !
    i3min = rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
    i3max = rismt%cfft%dfftt%npp(rismt%cfft%dfftt%mype + 1) + i3min
    !
    ALLOCATE(bg0(rismt%cfft%dfftt%nr3, rismt%nsite))
    bg0 = 0.0_DP
    !
!$omp parallel default(shared) private(ir, idx, i1, i2, i3, iz, isite, bg1)
#if defined(__OPENMP)
    ALLOCATE(bg1(rismt%cfft%dfftt%nr3, rismt%nsite))
    bg1 = 0.0_DP
#endif
!$omp do
    DO ir = 1, rismt%cfft%dfftt%nnr
      !
      idx = idx0 + ir - 1
      i3  = idx / (rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x)
      IF (i3 < i3min .OR. i3 >= i3max .OR. i3 >= rismt%cfft%dfftt%nr3) THEN
        CYCLE
      END IF
      !
      idx = idx - (rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x) * i3
      i2  = idx / rismt%cfft%dfftt%nr1x
      IF (i2 >= rismt%cfft%dfftt%nr2) THEN
        CYCLE
      END IF
      !
      idx = idx - rismt%cfft%dfftt%nr1x * i2
      i1  = idx
      IF (i1 >= rismt%cfft%dfftt%nr1) THEN
        CYCLE
      END IF
      !
      IF (i3 < (rismt%cfft%dfftt%nr3 - (rismt%cfft%dfftt%nr3 / 2))) THEN
        iz = i3 + (rismt%cfft%dfftt%nr3 / 2)
      ELSE
        iz = i3 - rismt%cfft%dfftt%nr3 + (rismt%cfft%dfftt%nr3 / 2)
      END IF
      iz = iz + 1
      !
      DO isite = 1, rismt%nsite
#if defined(__OPENMP)
        bg1(iz, isite) = bg1(iz, isite) + ar(ir, isite)
#else
        bg0(iz, isite) = bg0(iz, isite) + ar(ir, isite)
#endif
      END DO
      !
    END DO
!$omp end do
#if defined(__OPENMP)
!$omp critical
    bg0 = bg0 + bg1
!$omp end critical
    DEALLOCATE(bg1)
#endif
!$omp end parallel
    !
    CALL mp_sum(bg0, rismt%mp_site%intra_sitg_comm)
    !
    bg0 = bg0 / DBLE(rismt%cfft%dfftt%nr1 * rismt%cfft%dfftt%nr2)
    !
    DO isite = 1, rismt%nsite
      DO iz = rismt%lfft%izcell_start, rismt%lfft%izcell_end
        iiz = iz - rismt%lfft%izcell_start + 1
        ag0(iz, isite) = bg0(iiz, isite)
      END DO
    END DO
    !
    DEALLOCATE(bg0)
    !
  END SUBROUTINE extract_gxy0
  !
  SUBROUTINE sum_gxy0()
    IMPLICIT NONE
    INTEGER :: ir
    INTEGER :: idx
    INTEGER :: idx0
    INTEGER :: i3min
    INTEGER :: i3max
    INTEGER :: i1, i2, i3
    INTEGER :: iz
    INTEGER :: isite
    !
    IF (rismt%nsite < 1) THEN
      RETURN
    END IF
    !
    idx0 = rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x &
       & * rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
    !
    i3min = rismt%cfft%dfftt%ipp(rismt%cfft%dfftt%mype + 1)
    i3max = rismt%cfft%dfftt%npp(rismt%cfft%dfftt%mype + 1) + i3min
    !
!$omp parallel do default(shared) private(ir, idx, i1, i2, i3, iz, isite)
    DO ir = 1, rismt%cfft%dfftt%nnr
      !
      idx = idx0 + ir - 1
      i3  = idx / (rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x)
      IF (i3 < i3min .OR. i3 >= i3max .OR. i3 >= rismt%cfft%dfftt%nr3) THEN
        CYCLE
      END IF
      !
      idx = idx - (rismt%cfft%dfftt%nr1x * rismt%cfft%dfftt%nr2x) * i3
      i2  = idx / rismt%cfft%dfftt%nr1x
      IF (i2 >= rismt%cfft%dfftt%nr2) THEN
        CYCLE
      END IF
      !
      idx = idx - rismt%cfft%dfftt%nr1x * i2
      i1  = idx
      IF (i1 >= rismt%cfft%dfftt%nr1) THEN
        CYCLE
      END IF
      !
      IF (i3 < (rismt%cfft%dfftt%nr3 - (rismt%cfft%dfftt%nr3 / 2))) THEN
        iz = i3 + (rismt%cfft%dfftt%nr3 / 2)
      ELSE
        iz = i3 - rismt%cfft%dfftt%nr3 + (rismt%cfft%dfftt%nr3 / 2)
      END IF
      iz = iz + rismt%lfft%izcell_start
      !
      DO isite = 1, rismt%nsite
        ar(ir, isite) = ar(ir, isite) + ag0(iz, isite)
      END DO
      !
    END DO
!$omp end parallel do
    !
  END SUBROUTINE sum_gxy0
  !
END SUBROUTINE corrgxy0_laue
