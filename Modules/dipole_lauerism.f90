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
SUBROUTINE dipole_lauerism(rismt, ierr)
  !---------------------------------------------------------------------------
  !
  ! ... extract dipole part of direct correlations, in Laue-RISM calculation.
  !
  USE cell_base, ONLY : alat
  USE constants, ONLY : K_BOLTZMANN_RY
  USE err_rism,  ONLY : IERR_RISM_NULL, IERR_RISM_INCORRECT_DATA_TYPE
  USE kinds,     ONLY : DP
  USE mp,        ONLY : mp_sum
  USE rism,      ONLY : rism_type, ITYPE_LAUERISM
  USE solvmol,   ONLY : solVs, iuniq_to_isite, isite_to_isolV, isite_to_iatom
  !
  IMPLICIT NONE
  !
  TYPE(rism_type), INTENT(INOUT) :: rismt
  INTEGER,         INTENT(OUT)   :: ierr
  !
  INTEGER               :: iq
  INTEGER               :: iiq
  INTEGER               :: iv
  INTEGER               :: isolV
  INTEGER               :: iatom
  INTEGER               :: izsolv
  INTEGER               :: jzsolv
  REAL(DP)              :: qv
  REAL(DP)              :: beta
  REAL(DP)              :: zstep
  REAL(DP)              :: zoffs
  REAL(DP)              :: zedge
  REAL(DP)              :: vline0
  REAL(DP)              :: vline1
  REAL(DP), ALLOCATABLE :: cd0(:)
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
  IF (rismt%nrzs < rismt%cfft%dfftt%nr3) THEN
    ierr = IERR_RISM_INCORRECT_DATA_TYPE
    RETURN
  END IF
  !
  ! ... is one-hand ?
  IF (rismt%lfft%xright .AND. rismt%lfft%xleft) THEN
    !
    rismt%cdza = 0.0_DP
    rismt%cdsr = rismt%csr
    !
    ierr = IERR_RISM_NULL
    RETURN
  END IF
  !
  ! ... beta = 1 / (kB * T)
  beta = 1.0_DP / K_BOLTZMANN_RY / rismt%temp
  !
  ! ... set some variables
  IF (rismt%lfft%xright) THEN
    izsolv = rismt%lfft%izright_start
    !
    IF (rismt%lfft%gxystart > 1) THEN
      vline0 = AIMAG(rismt%vleft(1))
      vline1 = DBLE( rismt%vleft(1)) / alat
    END IF
    !
  ELSE !IF (rismt%lfft%xleft) THEN
    izsolv = rismt%lfft%izleft_end
    !
    IF (rismt%lfft%gxystart > 1) THEN
      vline0 = AIMAG(rismt%vright(1))
      vline1 = DBLE( rismt%vright(1)) / alat
    END IF
  END IF
  !
  ! ... properties about length (in a.u.)
  zstep = alat * rismt%lfft%zstep
  zoffs = alat * (rismt%lfft%zleft + rismt%lfft%zoffset)
  zedge = zoffs + zstep * DBLE(izsolv - 1)
  !
  ! ... allocate memory
  IF (rismt%nsite > 0) THEN
    ALLOCATE(cd0(rismt%nsite))
    cd0 = 0.0_DP
  END IF
  !
  ! ... calculate amplitudes of dipole part,
  ! ... which is included in short-range direct correlations
  DO iq = rismt%mp_site%isite_start, rismt%mp_site%isite_end
    iiq   = iq - rismt%mp_site%isite_start + 1
    iv    = iuniq_to_isite(1, iq)
    isolV = isite_to_isolV(iv)
    iatom = isite_to_iatom(iv)
    qv    = solVs(isolV)%charge(iatom)
    !
    IF (rismt%lfft%gxystart > 1) THEN
      jzsolv = izsolv - rismt%lfft%izcell_start + 1
      cd0(iiq) = DBLE(rismt%csgz(jzsolv, iiq)) &         ! short-range
             & - beta * qv * DBLE(rismt%vlgz(izsolv)) &  ! long-range
             & + beta * qv * (vline1 * zedge + vline0)   ! from void-region
    END IF
  END DO
  !
  IF (rismt%nsite > 0) THEN
    CALL mp_sum(cd0, rismt%mp_site%intra_sitg_comm)
  END IF
  !
  ! ... update amplitudes of dipole part of direct correlations
  IF (rismt%nsite > 0) THEN
    rismt%cdza = rismt%cdza + cd0
  END IF
  !
  ! ... update short-range direct correlations
  CALL update_short_range()
  !
  ! ... deallocate memoery
  IF (rismt%nsite > 0) THEN
    DEALLOCATE(cd0)
  END IF
  !
  ! ... normally done
  ierr = IERR_RISM_NULL
  !
CONTAINS
  !
  SUBROUTINE update_short_range()
    IMPLICIT NONE
    !
    INTEGER  :: ir
    INTEGER  :: idx
    INTEGER  :: idx0
    INTEGER  :: i3min
    INTEGER  :: i3max
    INTEGER  :: i1, i2, i3
    INTEGER  :: iz, iiz
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
!$omp parallel do default(shared) private(ir, idx, i1, i2, i3, iz, iiz)
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
      iiz = iz + 1
      iz  = iz + rismt%lfft%izcell_start
      !
      IF (iz > rismt%lfft%izright_end .OR. iz < rismt%lfft%izleft_start) THEN
        CYCLE
      END IF
      IF (iz < rismt%lfft%izright_start .AND. iz > rismt%lfft%izleft_end) THEN
        CYCLE
      END IF
      !
      rismt%csr (ir, :) = rismt%csr(ir, :) - cd0(:) * rismt%cdzs(iiz)
      rismt%csdr(ir, :) = rismt%csr(ir, :) + rismt%cdza(:) * rismt%cdzs(iiz)
      !
    END DO
!$omp end parallel do
    !
  END SUBROUTINE update_short_range
  !
END SUBROUTINE dipole_lauerism
